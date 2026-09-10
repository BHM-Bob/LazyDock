#!/usr/bin/env python
"""
Standalone ABFE run engine for LazyDock.
Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)

Migrated from BindFlow ``engine_fep.py`` (which itself is the Snakemake-free
version of the BindFlow FEP workflow).  The system building step is NOT
included here: it is done by ``prepare-abfe`` (lazydock.scripts.prepare_abfe),
which produces the ``input/{complex,ligand}`` directories consumed here.

The complete flow per ligand:
    ligand equilibrations -> complex equilibrations -> Boresch restraints
    -> FEP (ligand + complex) -> analyse -> gather.

Differences with the BindFlow engine (by design):
    * ``tools.gmx_runner`` is replaced by ``lazydock.gmx.run.Gromacs``
      (expect-based command layer, no subprocess piping).
    * No Snakemake; a step is finished if its gmx command returned 0.
    * No system building (``input/`` comes from prepare-abfe).
"""

import copy
import logging
import os
import shutil
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from lazydock.gmx.abfe import mdp
from lazydock.gmx.abfe import templates
from lazydock.gmx.abfe import boresch
from lazydock.gmx.abfe import fep_analysis, gather_results
from lazydock.gmx.abfe import tools

logger = logging.getLogger(__name__)

__all__ = ["run_abfe", "prepare_out_dir_and_config"]


###########################################################################
# Small helpers
###########################################################################

def _mol_basename(ligand_definition: dict) -> str:
    return Path(ligand_definition["conf"]).name


def _state_sort_key(path) -> int:
    # ``.../simulation/{lambda_type}.{state}/prod/prod.xvg`` -> state int
    return int(Path(path).parts[-3].split(".")[-1])


def _get_mdp_temperature(mdp_file) -> float:
    params = mdp.MDP().from_file(mdp_file).parameters
    if "ref-t" in params:
        return float(params["ref-t"].split()[0])
    elif "ref_t" in params:
        return float(params["ref_t"].split()[0])
    raise ValueError(f"Could not find temperature (ref-t/ref_t) in {mdp_file}")


def _first_mdp_path(sim_dir):
    """Return the path of the first available prod.mdp inside a fep
    simulation directory (for temperature extraction)."""
    prod_mdps = list(Path(sim_dir).glob("*/prod/prod.mdp"))
    if not prod_mdps:
        raise FileNotFoundError(f"No prod.mdp found under {sim_dir}")
    return prod_mdps[0]


def _link_or_copy(src, dst):
    """Copy a file (or symlink when possible) into a topology directory."""
    dst = Path(dst)
    if dst.exists():
        dst.unlink()
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy(src, dst)


###########################################################################
# Config preparation (equivalent to flow_builder.approach_flow)
###########################################################################

def _update_nwindows_config(config: dict) -> dict:
    nwindows_default = {
        "ligand": {"vdw": 11, "coul": 11},
        "complex": {"vdw": 21, "coul": 11, "bonded": 11},
    }
    nwindows = config.get("nwindows") or {}
    for key in ["ligand", "complex"]:
        if key in nwindows:
            nwindows_default[key].update(nwindows[key])
    config["nwindows"] = nwindows_default
    return config


def prepare_out_dir_and_config(global_config: dict) -> dict:
    """
    Replicate the directory creation + config assembly done in
    ``flow_builder.approach_flow`` (minus the Snakefile/scheduler parts).

    It mutates and returns ``global_config`` with the keys: ``lambdas``,
    ``complex_type``, ``mdp`` (user overrides), ``job_prefix`` ...

    The input `input/{complex,ligand}` directory produced by prepare-abfe is
    linked into ``out_approach_path/{ligand_name}/input`` (BindFlow layout).
    """
    import numpy as np

    global_config = copy.deepcopy(global_config)
    out_path = Path(global_config["out_approach_path"])
    out_path.mkdir(exist_ok=True, parents=True)

    if global_config["calculation_type"] == "fep":
        _update_nwindows_config(global_config)
        global_config["lambdas"] = {
            "ligand": {
                "vdw": list(np.round(np.linspace(0, 1, global_config["nwindows"]["ligand"]["vdw"]), 2)),
                "coul": list(np.round(np.linspace(0, 1, global_config["nwindows"]["ligand"]["coul"]), 2)),
            },
            "complex": {
                "vdw": list(np.round(np.linspace(0, 1, global_config["nwindows"]["complex"]["vdw"]), 2)),
                "coul": list(np.round(np.linspace(0, 1, global_config["nwindows"]["complex"]["coul"]), 2)),
                "bonded": list(np.round(np.linspace(0, 1, global_config["nwindows"]["complex"]["bonded"]), 2)),
            },
        }

    # Specify the complex type
    if global_config["inputs"].get("membrane"):
        global_config["complex_type"] = "membrane"
    else:
        global_config["complex_type"] = "soluble"

    # Extra mdp options if provided
    if "mdp" in global_config:
        global_config["mdp"] = global_config["mdp"]
    else:
        global_config["mdp"] = None

    # Directory structure
    input_src = Path(global_config["input_dir"])
    for i, ligand_definition in enumerate(global_config["inputs"]["ligands"]):
        input_ligand_path = Path(ligand_definition["conf"])
        # use the canonical ligand name (ligand_names), NOT the conf file stem:
        # the conf may be an intermediate file (e.g. input/ligand/ligand.gro)
        # whose stem does not match the ligand name.
        ligand_name = global_config["ligand_names"][i]
        out_ligand_path = out_path / ligand_name
        out_ligand_input_path = out_ligand_path / "input"
        (out_ligand_input_path / "complex").mkdir(exist_ok=True, parents=True)
        (out_ligand_input_path / "ligand").mkdir(exist_ok=True, parents=True)

        # Link the prepare-abfe input files (complex + ligand folders)
        for sub in ("complex", "ligand"):
            src_sub = input_src / sub
            dst_sub = out_ligand_input_path / sub
            if not dst_sub.exists() or not any(dst_sub.iterdir()):
                if not src_sub.exists():
                    raise FileNotFoundError(f"prepare-abfe input {src_sub} not found; "
                                            f"run prepare-abfe first (produces input/complex, input/ligand)")
                for f in src_sub.iterdir():
                    if f.is_file():
                        _link_or_copy(f, dst_sub / f.name)

        # Archive original files (informative only)
        with tarfile.open(out_ligand_input_path / "orig_in.tar.gz", "w:gz") as tar:
            if input_ligand_path.exists():
                tar.add(input_ligand_path, arcname=input_ligand_path.name)
            protein_conf = global_config["inputs"]["protein"]["conf"]
            tar.add(protein_conf, arcname=Path(protein_conf).name)
            if global_config["inputs"].get("cofactor"):
                cof_conf = global_config["inputs"]["cofactor"]["conf"]
                tar.add(cof_conf, arcname=Path(cof_conf).name)
            if global_config["inputs"].get("membrane"):
                mem_conf = global_config["inputs"]["membrane"]["conf"]
                tar.add(mem_conf, arcname=Path(mem_conf).name)

        # Build the replicas
        for num_replica in range(1, global_config["replicas"] + 1):
            out_replica_path = out_ligand_path / str(num_replica)
            out_replica_path.mkdir(exist_ok=True, parents=True)

    return global_config


###########################################################################
# MDP generation
###########################################################################

def _make_equi_mdps(global_config: dict, sys_type: str) -> None:
    """Equivalent of equil_setup_ligand / equil_setup_complex rules."""
    out_root = Path(global_config["out_approach_path"])
    if sys_type == "ligand":
        template_dir = templates.TemplatePath.ligand.equi
    else:
        if global_config["complex_type"] == "membrane":
            template_dir = templates.TemplatePath.complex.membrane.equi
        else:
            template_dir = templates.TemplatePath.complex.soluble.equi
    steps = sorted(step.stem for step in tools.list_if_file(template_dir, ext=".mdp"))

    for ligand_name in global_config["ligand_names"]:
        for replica in map(str, range(1, global_config["replicas"] + 1)):
            sim_dir = out_root / ligand_name / replica / sys_type / "equil-mdsim"
            mdp_template = mdp.StepMDP(step_path=template_dir)
            for step in steps:
                (sim_dir / step).mkdir(exist_ok=True, parents=True)
                output_mdp = sim_dir / step / f"{step}.mdp"
                mdp_template.set_new_step(step)
                # 平衡阶段积分器已直接在模板中设为 md+v-rescale (2026-09, 见模板注释)。
                # --equi-integrator sd 提供回退: 用户想用模板旧行为时可显式改回 sd。
                # membrane 模板本来就是 md+v-rescale, 传 sd 无意义, 忽略。
                if (step != '00_min' and global_config.get("equi_integrator", "md") == "sd"
                        and mdp_template.parameters.get("integrator", "").strip() == "md"
                        and "System" in mdp_template.parameters.get("tc-grps", "")):
                    mdp_template.set_parameters(integrator="sd")
                    # sd 需要 Langevin 恒温, 去掉 v-rescale 的 tcoupl 不合法? 
                    # v-rescale 与 sd 的 tcoupl 键均可存在(grompp 检查一致); 直接设 sd 即可,
                    # tcoupl 保留 v-rescale 会导致 sd+外部热浴 双恒温, 应删除 tcoupl.
                    if "tcoupl" in mdp_template.parameters:
                        mdp_template.parameters.pop("tcoupl")
                if "dt" in mdp_template.parameters:
                    if float(mdp_template.parameters["dt"].split(";")[0]) > global_config["dt_max"]:
                        mdp_template.set_parameters(dt=global_config["dt_max"])
                user_mdp = (global_config.get("mdp") or {}).get(sys_type, {}).get("equi", {}).get(step)
                if user_mdp:
                    mdp_template.set_parameters(**user_mdp)
                mdp_template.write(output_mdp)


def _make_fep_mdps(global_config: dict, sys_type: str) -> None:
    """Equivalent of fep_setup_ligand / fep_setup_complex rules.

    For complex, also copies the (boresch-augmented) topology into
    ``{rep}/complex/fep/topology/`` and the index file.
    """
    out_root = Path(global_config["out_approach_path"])
    if sys_type == "ligand":
        template_root = templates.TemplatePath.ligand.fep
        lambda_types = ["vdw", "coul"]
    else:
        if global_config["complex_type"] == "membrane":
            template_root = templates.TemplatePath.complex.membrane.fep
        else:
            template_root = templates.TemplatePath.complex.soluble.fep
        lambda_types = ["vdw", "coul", "bonded"]

    mdp_extra_kwargs = (global_config.get("mdp") or {}).get(sys_type, {}).get("fep", {})
    # Flexible ligands (peptides) keep their intramolecular non-bonded
    # interactions fully on during decoupling (couple-intramol=no) to avoid
    # chain unfolding in the fully-decoupled windows; cancels in the cycle.
    # NOTE: this does NOT apply to the "bonded" windows: there couple-lambda0/1
    # are both "vdw-q" (ligand non-bonded interactions are off for the whole
    # window, no chain can unfold), and couple-intramol=no there would make
    # GROMACS check every intramolecular excluded pair against rlist — 2.8 nm for
    # a 16-residue cyclic peptide — which fails grompp ("largest distance between
    # non-perturbed excluded atoms ... larger than the cut-off").  BindFlow runs
    # bonded windows with couple-intramol=yes.
    peptide_couple_no = ('no' if global_config.get("peptide_ligand") else 'yes')
    fep_rlist = global_config.get("fep_rlist")
    fep_cutoff = global_config.get("fep_cutoff")
    fep_dt_max = global_config.get("fep_dt_max", global_config["dt_max"])

    for ligand_name in global_config["ligand_names"]:
        for replica in map(str, range(1, global_config["replicas"] + 1)):
            sim_dir = out_root / ligand_name / replica / sys_type / "fep"
            for lambda_type in lambda_types:
                lambda_values = global_config["lambdas"][sys_type][lambda_type]
                couple_intramol = ('yes' if lambda_type == 'bonded' else peptide_couple_no)
                mdp.make_fep_dir_structure(
                    sim_dir=sim_dir,
                    template_dir=template_root,
                    lambda_values=lambda_values,
                    lambda_type=lambda_type,
                    sys_type=sys_type,
                    dt_max=fep_dt_max,
                    mdp_extra_kwargs=mdp_extra_kwargs,
                    couple_intramol=couple_intramol,
                    fep_rlist=fep_rlist,
                    fep_cutoff=fep_cutoff,
                )

            if sys_type == "complex":
                # ---- Topology directory (complex only) ----
                out_top_dir = sim_dir / "topology"
                out_top_dir.mkdir(exist_ok=True, parents=True)
                out_ligand_input = out_root / ligand_name / "input" / "complex"
                for itp_file in tools.list_if_file(out_ligand_input, ext=".itp"):
                    _link_or_copy(itp_file, out_top_dir / itp_file.name)
                for ndx_file in tools.list_if_file(out_ligand_input, ext=".ndx"):
                    _link_or_copy(ndx_file, out_top_dir / ndx_file.name)
                # Modify the main topology incorporating the boresch restraints
                complex_top = out_ligand_input / "complex.top"
                boresch_top = (out_root / ligand_name / replica / "complex"
                               / "equil-mdsim" / "boreschcalc" / "BoreschRestraint.top")
                fep_top = out_top_dir / "complex_boresch.top"
                with open(complex_top, "r") as original_top, \
                        open(boresch_top, "r") as boresch_top_f:
                    with open(fep_top, "w") as final_top:
                        final_top.write(original_top.read() + boresch_top_f.read())


###########################################################################
# Simulation steps
###########################################################################

def _run_step_chain(global_config, sys_type, ligand_name, replica, chain):
    """Run a chain of [step][gro/cpt] entries.

    chain: list of dicts with keys: step, run_dir, mdp, top, gro (input),
           cpt (input, optional), out_gro, out_cpt, output_finished,
           minimize (bool, for 00_min steps)
    """
    gmx = global_config.get("gmx")
    mdrun_extra = global_config["extra_directives"]["mdrun"][sys_type]
    nthreads = global_config["threads"]
    ntmpi = global_config.get("ntmpi")
    gpu_id = global_config.get("gpu_id")
    maxwarn = global_config.get("maxwarn", 2)
    retries = global_config.get("retries", 3)

    for entry in chain:
        run_dir = Path(entry["run_dir"])
        mdp_file = entry["mdp"]
        minimize = entry.get("minimize")
        # Skip steps already finished (crashed runs resume without wasting
        # completed equil/prod steps: .finished marker + output files exist).
        # Minimization steps have no .finished marker; use the output gro.
        out_gro_entry = entry.get("out_gro")
        if (entry.get("out_finished") and (run_dir / f"{entry['step']}.finished").exists()
                and Path(entry["out_gro"]).exists() and Path(entry["out_cpt"]).exists()) \
                or (minimize and out_gro_entry and Path(out_gro_entry).exists()):
            logger.info(f"step {entry['step']} already finished - skipping")
            continue
        # 00_min: minimization has no update/bonded/pme switches
        mdrun_extra_use = mdrun_extra
        if minimize:
            mdrun_extra_use = mdrun_extra.copy()
            for invalid_flag in ["update", "bonded", "pme"]:
                mdrun_extra_use.pop(invalid_flag, None)

        # Retry the (grompp+mdrun) pair like the original ``retries`` rule
        # directive did. ``mdrun -cpi`` continues from any partial checkpoint.
        last_exc = None
        for attempt in range(retries + 1):
            try:
                tools.gmx_runner(
                    gmx=gmx,
                    mdp=mdp_file,
                    topology=entry["top"],
                    structure=entry["gro"],
                    checkpoint=entry.get("cpt"),
                    index=entry.get("index"),
                    nthreads=nthreads,
                    run_dir=run_dir,
                    maxwarn=maxwarn,
                    minimize=minimize,
                    gpu_id=gpu_id,
                    ntmpi=ntmpi,
                    **mdrun_extra_use,
                )
                last_exc = None
                break
            except RuntimeError as exc:
                last_exc = exc
                if attempt < retries:
                    logger.warning(f"gmx step {entry['step']} failed (attempt "
                                   f"{attempt + 1}/{retries + 1}): {exc} -- retrying")
        if last_exc is not None:
            raise RuntimeError(f"gmx step {entry['step']} failed after "
                               f"{retries + 1} attempts: {last_exc}")

        if entry.get("out_finished"):
            tools.paths_exist(
                paths=[entry["out_gro"], entry["out_cpt"]],
                raise_error=True,
                out=entry["out_finished"],
            )


def _equi_chain_steps(global_config, sys_type) -> list:
    """Return the ordered list of equilibration step names for a system."""
    template_dir = None
    if sys_type == "ligand":
        template_dir = templates.TemplatePath.ligand.equi
    else:
        if global_config["complex_type"] == "membrane":
            template_dir = templates.TemplatePath.complex.membrane.equi
        else:
            template_dir = templates.TemplatePath.complex.soluble.equi
    return sorted(step.stem for step in tools.list_if_file(template_dir, ext=".mdp"))


def _run_equilibration(global_config, sys_type, ligand_name, replica):
    """Equivalent of equil_{ligand|complex}_simulation.smk"""
    out_root = Path(global_config["out_approach_path"])
    input_dir = out_root / ligand_name / "input" / sys_type
    steps = _equi_chain_steps(global_config, sys_type)
    top = input_dir / "ligand.top" if sys_type == "ligand" else input_dir / "complex.top"
    structure = input_dir / "ligand.gro" if sys_type == "ligand" else input_dir / "complex.gro"
    chain = []

    prev_gro = structure
    prev_cpt = None
    for idx, step in enumerate(steps):
        run_dir = out_root / ligand_name / replica / sys_type / "equil-mdsim" / step
        mdp_file = run_dir / f"{step}.mdp"
        is_min = (step == steps[0])
        out_gro = run_dir / f"{step}.gro"
        entry = {
            "step": step,
            "run_dir": run_dir,
            "mdp": mdp_file,
            "top": top,
            # chain: 00_min starts from the input structure, the rest from the
            # output gro/cpt of the previous step
            "gro": prev_gro,
            "cpt": prev_cpt,
            "out_gro": out_gro,
        }
        if not is_min:
            out_cpt = run_dir / f"{step}.cpt"
            entry["out_cpt"] = out_cpt
            entry["out_finished"] = run_dir / f"{step}.finished"
            prev_cpt = out_cpt
        else:
            entry["minimize"] = True
        # Membrane complex requires the index file for every step
        if sys_type == "complex" and global_config["complex_type"] == "membrane":
            entry["index"] = input_dir / "index.ndx"
        chain.append(entry)
        prev_gro = out_gro

    _run_step_chain(global_config, sys_type, ligand_name, replica, chain)


def _boresch_finished(run_dir: Path) -> bool:
    """Whether a previous Boresch restraint generation completed.

    ``BoreschRestraint.top`` + ``dG_off.dat`` are the two artifacts consumed
    downstream (topology for FEP windows, analytical dG_off for analysis).
    Requiring both avoids resuming on a half-written directory.
    """
    return (run_dir / "BoreschRestraint.top").exists() and (run_dir / "dG_off.dat").exists()


def _run_boresch(global_config, ligand_name, replica):
    """Equivalent of equi/boresch.smk (FEP only)."""
    out_root = Path(global_config["out_approach_path"])
    complex_equil_dir = (out_root / ligand_name / replica / "complex" / "equil-mdsim")
    prod_dir = complex_equil_dir / "prod"
    run_dir = complex_equil_dir / "boreschcalc"
    run_dir.mkdir(exist_ok=True, parents=True)

    # Boresch restraint generation is deterministic given the prod trajectory
    # (center -> FindBoreschRestraint picks the lowest-variance combination).
    # Once BoreschRestraint.top + dG_off.dat exist, skip regeneration: FEP
    # topologies were already built from this exact restraint, and re-running
    # could select a different anchor set, silently breaking consistency
    # between existing FEP windows and the analysis correction.
    if _boresch_finished(run_dir):
        logger.info("Boresch restraints already generated - skipping")
        return

    # Fix trajectory.
    tools.center_xtc(
        tpr=prod_dir / "prod.tpr",
        xtc=prod_dir / "prod.xtc",
        run_dir=run_dir,
        host_name=global_config.get("host_name", "Protein"),
        gmx=global_config.get("gmx"),
    )

    # Getting Boresch restraints (this is the ClosestRestraintFrame producer)
    temperature = _get_mdp_temperature(prod_dir / "prod.mdp")
    boresch.gen_restraint(
        topology=prod_dir / "prod.tpr",
        trajectory=run_dir / "center.xtc",
        outpath=run_dir,
        temperature=temperature,
        host_selection=global_config.get("host_selection", "protein and name CA and not moltype LIG"),
    )
    # Clean
    (run_dir / "center.xtc").unlink()


def _collect_fep_chains(global_config, sys_type, ligand_name, replica):
    """Return a list of (chain, window_label) for all FEP windows of one leg.

    Each chain is the full step chain of one lambda window
    (00_min -> 01_nvt -> 02_npt -> 03_npt_norest -> prod), i.e. the natural
    atomic unit for task-level parallelism (windows are independent).
    """
    out_root = Path(global_config["out_approach_path"])
    if sys_type == "ligand":
        top = out_root / ligand_name / "input" / "ligand" / "ligand.top"
        start_gro = (out_root / ligand_name / replica / "ligand" / "equil-mdsim"
                     / "prod" / "prod.gro")
        lambda_types = ["vdw", "coul"]
    else:
        top = (out_root / ligand_name / replica / "complex" / "fep"
               / "topology" / "complex_boresch.top")
        start_gro = (out_root / ligand_name / replica / "complex" / "equil-mdsim"
                     / "boreschcalc" / "ClosestRestraintFrame.gro")
        lambda_types = ["vdw", "coul", "bonded"]

    chains = []
    for lambda_type in lambda_types:
        lambda_values = global_config["lambdas"][sys_type][lambda_type]
        for state in range(len(lambda_values)):
            sim_root = (out_root / ligand_name / replica / sys_type / "fep"
                        / "simulation" / f"{lambda_type}.{state}")
            step_dirs = sorted(p.name for p in tools.list_if_dir(sim_root))
            chain = []
            prev_gro = start_gro
            prev_cpt = None
            for idx, step in enumerate(step_dirs):
                run_dir = sim_root / step
                mdp_file = run_dir / f"{step}.mdp"
                is_min = (idx == 0)
                out_gro = run_dir / f"{step}.gro"
                entry = {
                    "step": step,
                    "run_dir": run_dir,
                    "mdp": mdp_file,
                    "top": top,
                    "gro": prev_gro,
                    "cpt": prev_cpt,  # None for 00_min and 01_nvt (gen-vel)
                    "out_gro": out_gro,
                }
                if not is_min:
                    out_cpt = run_dir / f"{step}.cpt"
                    entry["out_cpt"] = out_cpt
                    entry["out_finished"] = run_dir / f"{step}.finished"
                    prev_cpt = out_cpt
                else:
                    entry["minimize"] = True
                if sys_type == "complex" and global_config["complex_type"] == "membrane":
                    # membrane complex uses the index file copied to fep/topology
                    entry["index"] = (out_root / ligand_name / replica / "complex"
                                      / "fep" / "topology" / "index.ndx")
                chain.append(entry)
                prev_gro = out_gro
            chains.append((chain, f"{lambda_type}.{state}"))
    return chains


def _run_fep_simulation(global_config, sys_type, ligand_name, replica):
    """Equivalent of fep_{ligand|complex}_simulation.smk (sequential)."""
    for chain, label in _collect_fep_chains(global_config, sys_type, ligand_name, replica):
        _run_step_chain(global_config, sys_type, ligand_name, replica, chain)


def _run_fep_simulation_parallel(global_config, ligand_name, replica):
    """Run all FEP windows (ligand + complex legs) concurrently.

    Two-level scheduler: a pool of ``n_parallel`` window workers; each worker
    runs one full window chain (00_min -> ... -> prod) pinned to a GPU from
    ``gpu_ids`` (round-robin).  Windows that already finished are skipped
    inside ``_run_step_chain`` via the ``.finished`` marker, so re-runs only
    execute the remaining windows.  Equilibration/Boresch are still sequential
    (they must complete before any FEP window).
    """
    n_parallel = int(global_config.get("n_parallel", 1))
    gpu_ids = list(global_config.get("gpu_ids") or [global_config.get("gpu_id", 0)])
    n_gpu = len(gpu_ids) if gpu_ids else 1
    # threads per window: keep a sane fraction of the requested total threads
    nthreads = int(global_config.get("threads", 12))

    jobs = []
    for sys_type in ("ligand", "complex"):
        for chain, label in _collect_fep_chains(global_config, sys_type, ligand_name, replica):
            jobs.append((sys_type, chain, label))
    logger.info(f"FEP parallel: {len(jobs)} windows, {n_parallel} concurrent, "
                f"GPUs={gpu_ids}, threads/window={nthreads}")

    def _worker(task):
        sys_type, chain, label = task
        # pin this worker to one GPU (round-robin over gpu_ids).  The shared
        # gmx object is safe: _run_step_chain -> tools._gmx_for_run_dir creates
        # a per-run-dir Gromacs that inherits gpu_ids, and mdrun -gpu_id is set
        # from this worker's gpu_id override.
        wid = jobs.index(task)
        gpu = gpu_ids[wid % n_gpu]
        wcfg = dict(global_config)   # shallow: only gpu_id/threads overridden
        wcfg["gpu_id"] = gpu
        # Divide the requested total threads across concurrent windows so we do
        # not oversubscribe the CPU (GROMACS OpenMP scales poorly beyond a few
        # threads per tiny window anyway).  Users can still force a value with
        # --mdrun-args "nt=N" (overrides mdrun_extra).
        wcfg["threads"] = max(1, nthreads // max(1, n_parallel))
        try:
            _run_step_chain(wcfg, sys_type, ligand_name, replica, chain)
            return label, None
        except Exception as exc:  # noqa: BLE001  surface worker errors to caller
            return label, exc

    results = {}
    with ThreadPoolExecutor(max_workers=n_parallel) as pool:
        for label, exc in pool.map(_worker, jobs):
            results[label] = exc
            if exc is not None:
                logger.error(f"FEP window {label} failed: {exc}")
    failed = {k: v for k, v in results.items() if v is not None}
    if failed:
        raise RuntimeError(f"{len(failed)} FEP window(s) failed: "
                           + "; ".join(f"{k}: {v}" for k, v in failed.items()))



def _run_fep_analysis(global_config, ligand_name, replica):
    """Equivalent of fep ana.smk + calculate_result.smk"""
    out_root = Path(global_config["out_approach_path"])
    ligand_fep_dir = out_root / ligand_name / replica / "ligand" / "fep"
    complex_fep_dir = out_root / ligand_name / replica / "complex" / "fep"

    # --- ligand contributions ---
    lig_ana = ligand_fep_dir / "ana"
    lig_ana.mkdir(exist_ok=True, parents=True)
    xvg_vdw = sorted(
        (ligand_fep_dir / "simulation").glob("vdw.*/prod/prod.xvg"),
        key=_state_sort_key,
    )
    xvg_coul = sorted(
        (ligand_fep_dir / "simulation").glob("coul.*/prod/prod.xvg"),
        key=_state_sort_key,
    )
    temperature = _get_mdp_temperature(_first_mdp_path(ligand_fep_dir / "simulation"))
    fep_analysis.get_dG_contributions(
        boresch_data=None,
        out_json_path=lig_ana / "dg_ligand_contributions.json",
        lower=None,
        upper=None,
        min_samples=500,
        temperature=temperature,
        convergency_plots_prefix=None,
        vdw=xvg_vdw,
        coul=xvg_coul,
    )

    # --- complex contributions ---
    com_ana = complex_fep_dir / "ana"
    com_ana.mkdir(exist_ok=True, parents=True)
    xvg_vdw_c = sorted(
        (complex_fep_dir / "simulation").glob("vdw.*/prod/prod.xvg"),
        key=_state_sort_key,
    )
    xvg_coul_c = sorted(
        (complex_fep_dir / "simulation").glob("coul.*/prod/prod.xvg"),
        key=_state_sort_key,
    )
    xvg_bonded = sorted(
        (complex_fep_dir / "simulation").glob("bonded.*/prod/prod.xvg"),
        key=_state_sort_key,
    )
    # Boresch dG_off.dat file
    dg_off = (out_root / ligand_name / replica / "complex" / "equil-mdsim"
              / "boreschcalc" / "dG_off.dat")
    temperature_c = _get_mdp_temperature(_first_mdp_path(complex_fep_dir / "simulation"))

    fep_analysis.get_dG_contributions(
        boresch_data=dg_off,
        out_json_path=com_ana / "dg_complex_contributions.json",
        lower=None,
        upper=None,
        min_samples=500,
        temperature=temperature_c,
        convergency_plots_prefix=None,
        vdw=xvg_vdw_c,
        coul=xvg_coul_c,
        bonded=xvg_bonded,
    )

    # --- cycle ---
    fep_analysis.get_dg_cycle(
        ligand_contributions=lig_ana / "dg_ligand_contributions.json",
        complex_contributions=com_ana / "dg_complex_contributions.json",
        out_csv=out_root / ligand_name / replica / "dG_results.csv",
    )


###########################################################################
# Full pipeline
###########################################################################

def run_abfe(
    global_config: dict,
    only_build: bool = False,
):
    """Run the complete ABFE workflow sequentially (no Snakemake).

    Parameters
    ----------
    global_config : dict
        Config containing at least: calculation_type, out_approach_path,
        input_dir, inputs, water_model, host_name, host_selection,
        ligand_names, replicas, hmr_factor, custom_ff_path, threads,
        extra_directives, retries, dt_max, lambdas, complex_type.
    only_build : bool, optional
        If True, only generate the mdp structure and the directories,
        do not run any simulation (useful for debugging). By default False.

    Returns
    -------
    None
    """
    cfg = prepare_out_dir_and_config(global_config)

    # --- MDP generation (upfront, like the setup rules) ---
    _make_equi_mdps(cfg, "ligand")
    _make_equi_mdps(cfg, "complex")

    if only_build:
        logger.info("only_build=True: MDP structure generated; "
                    "no simulation launched")
        return

    # --- Equilibrations (Ligand first, then complex) ---
    for ligand_name in cfg["ligand_names"]:
        for replica in map(str, range(1, cfg["replicas"] + 1)):
            _run_equilibration(cfg, "ligand", ligand_name, replica)
            _run_equilibration(cfg, "complex", ligand_name, replica)
            _run_boresch(cfg, ligand_name, replica)

    # --- FEP setup (needs Boresch restraints generated above) ---
    _make_fep_mdps(cfg, "ligand")
    _make_fep_mdps(cfg, "complex")

    # --- FEP simulations ---
    n_parallel = int(cfg.get("n_parallel", 1))
    for ligand_name in cfg["ligand_names"]:
        for replica in map(str, range(1, cfg["replicas"] + 1)):
            if n_parallel > 1:
                # two-level scheduler: window workers run concurrently,
                # each pinned to a GPU from gpu_ids
                _run_fep_simulation_parallel(cfg, ligand_name, replica)
            else:
                _run_fep_simulation(cfg, "ligand", ligand_name, replica)
                _run_fep_simulation(cfg, "complex", ligand_name, replica)
            _run_fep_analysis(cfg, ligand_name, replica)

    # --- Gather ---
    gather_results.get_all_fep_dgs(
        root_folder_path=cfg["out_approach_path"],
        out_csv=Path(cfg["out_approach_path"]) / "fep_results.csv",
    )
    gather_results.get_raw_fep_data(
        root_folder_path=cfg["out_approach_path"],
        out_csv=Path(cfg["out_approach_path"]) / "fep_results_raw.csv",
    )