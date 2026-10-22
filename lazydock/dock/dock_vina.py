import os
import tempfile
from typing import List, Optional, Union

from mbapy_lite.file import opts_file

from lazydock.dock.adfr import (_prepare_ligand_receptor,
                                prepare_ligand_receptor_grid)
from lazydock.dock.config import (MACROCYCLE_MAX_RING_SIZE,
                                  PREPARE_BACKEND_DEFAULT)
from lazydock.dock.flex_utils import FlexResidue
from lazydock.dock.meeko_utils import prepare_flex_docking_setup


def calc_vina_score(score_name: str, grid_buffer: float,
                    complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str], cpu: int = 1,
                    spacing: float = 0.375,
                    macrocycle: bool = False, max_ring_size: int = MACROCYCLE_MAX_RING_SIZE,
                    prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """"total", "lig_inter", "flex_inter", "other_inter", "flex_intra", "lig_intra", "torsions", "-lig_intra"
    calculate vina score of a complex pdb file.
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box padding in ANGSTROM, added on every side of the ligand
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param cpu: cpu number
    :param spacing: grid spacing in ANGSTROM, must match the value used for docking
    :param macrocycle: prepare the ligand with meeko so that a large ring is broken and
        stays flexible, matching how it was prepared for docking.  Without it the
        ligand is re-prepared by ADT, whose rigid ring yields a different rotor count
        and therefore a `total` that can not be compared with the dock score.
    :param max_ring_size: largest ring meeko may break
    :param prepare_backend: 'meeko' or 'adfr', the same value that prepared the docked
        pose MUST be passed here, otherwise the receptor/ligand pdbqt differ from the
        docked ones and the rescored `total` is not comparable with the dock score
    
    :return: score list, if None, means failed
    """
    from vina import Vina
    v = Vina(sf_name=score_name, cpu=cpu, verbosity=False)    
    with tempfile.TemporaryDirectory() as tmpdir:
        # keyword args beyond the grid: this helper's parameter order is not the
        # caller's, and a positional slip here silently swaps macrocycle for
        # max_ring_size (measured: it made every plain dock take the macrocycle path)
        v = prepare_ligand_receptor_grid(tmpdir, v, complex_pdbstr, ligand_chain, receptor_chain,
                                          grid_buffer, spacing=spacing,
                                          macrocycle=macrocycle, max_ring_size=max_ring_size,
                                          prepare_backend=prepare_backend)
        if not v:
            return None
        # compute and retrun score
        return v[0].score().tolist()


def perform_vina_dock(score_name: str, grid_buffer: float,
                      complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                      exhaustiveness: int, n_poses: int, min_rmsd: float = 1, max_evals: int = 0,
                      seed: int = None, energy_range: float = None, cpu: Optional[int] = None,
                      spacing: float = 0.375,
                      macrocycle: bool = False, max_ring_size: int = MACROCYCLE_MAX_RING_SIZE,
                      prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """
    perform docking with vina
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box padding in ANGSTROM, added on every side of the ligand
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param exhaustiveness: exhaustiveness, equal to cpu number
    :param n_poses: n_poses
    :param min_rmsd: min_rmsd
    :param max_evals: max_evals
    :param seed: random seed for vina, 0 means vina picks one at random
        (non-reproducible).  Added for the redock pipeline, which promises
        reproducibility (its pipeline.json stores a seed); without it the same
        input gave -8.62 and -7.36 kcal/mol on two consecutive runs.  Kept
        optional so that every existing caller keeps its old behaviour.
        NOTE: None is accepted as a synonym of 0 and is converted below, because
        vina's SWIG binding rejects None outright (`argument 3 of type 'int'`)
        and `redock` calls this function without a seed at all.
    :param energy_range: maximum energy difference (kcal/mol) a pose may have
        from the best one to still be written.  None means vina's own default
        (3.0), which SILENTLY DROPS poses: requesting n_poses=20 returned only
        14 in a measured run, quietly breaking the pipeline's 160-pose-per-round
        contract.  Pass a large value (e.g. 100) to actually get `n_poses`.
    :param cpu: threads for vina, None keeps the historical behaviour of tying
        the thread count to `exhaustiveness`.  Those are different things and
        tying them bites hard: a caller asking for a cheap search
        (exhaustiveness=2) silently got a 2-thread dock, measured at ~10 min per
        task where 8 threads needed ~2 min.  The redock pipeline passes its own
        `vina_cpu` here (see `redock/loop._vina_dock`); every existing caller is
        unaffected because the default reproduces the old value.
    :param spacing: grid spacing in ANGSTROM, default 0.375
    :param macrocycle: prepare the ligand with meeko, breaking a large ring so that
        its conformation can adapt during the search (requires vina >= 1.2).  Off by
        default so that existing results stay reproducible.
    :param max_ring_size: largest ring meeko may break, see MACROCYCLE_MAX_RING_SIZE
    :param prepare_backend: 'meeko' or 'adfr', which toolkit prepares the receptor and
        ligand pdbqt.  The two type a few nitrogens/polar hydrogens differently, so
        scores are NOT comparable across backends (see `PREPARE_BACKENDS`)
    
    :return: (result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr,
        ligand_pdbqtstr), None if failed
    """
    from vina import Vina

    # vina's SWIG binding is typed, so a None here is a hard TypeError rather than
    # a silent default; `redock` relies on this default and passes no seed at all
    if seed is None:
        seed = 0
    # historically the thread count was tied to `exhaustiveness`; keep that as
    # the default so existing callers behave exactly as before, but let the
    # caller override (the redock pipeline sizes its own thread budget)
    if cpu is None:
        cpu = exhaustiveness
    v = Vina(sf_name=score_name, cpu=cpu, seed=seed, verbosity=False)    
    with tempfile.TemporaryDirectory() as tmpdir:
        # keyword args beyond the grid: see the note in `calc_vina_score`
        v = prepare_ligand_receptor_grid(tmpdir, v, complex_pdbstr, ligand_chain, receptor_chain,
                                          grid_buffer, spacing=spacing,
                                          macrocycle=macrocycle, max_ring_size=max_ring_size,
                                          prepare_backend=prepare_backend)
        if not v:
            return None
        v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = v
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=min_rmsd, max_evals=max_evals)
        output_path = os.path.join(tmpdir, 'poses.pdbqt')
        write_kwargs = {} if energy_range is None else {'energy_range': float(energy_range)}
        v.write_poses(output_path, n_poses, overwrite=True, **write_kwargs)
        result_str = opts_file(output_path, 'r', way='str')
        os.remove(output_path)
        return result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def perform_vina_flex_dock(score_name: str, grid_buffer: float,
                           complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                           flex_residues: Union[str, List[str], List[FlexResidue]] = None,
                           flex_distance: float = None, flex_interaction: bool = False,
                           flex_interaction_mode: Union[str, List[str]] = 'all', flex_interaction_cutoff: float = 4.,
                           exhaustiveness: int = 32, n_poses: int = 20,
                           min_rmsd: float = 1, max_evals: int = 0,
                           seed: int = 0, cpu: int = 0, energy_range: float = None,
                           spacing: float = 0.375,
                           macrocycle: bool = False, max_ring_size: int = MACROCYCLE_MAX_RING_SIZE,
                           prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """
    perform docking with vina, receptor flexible residues are supported.
    
    flexible residues could be provided in three ways, and they could be combined:
        - flex_residues: explicit residue list, such as 'A:LYS68,A:LYS50' or [(chain, resname, resnum), ...]
        - flex_distance: select receptor residues within this distance of any ligand heavy atom
        - flex_interaction: select receptor residues which have interaction with ligand by fplip
    
    The receptor is split into rigid / flexible parts by meeko, then handed to
    `Vina.set_receptor(rigid, flex)`.  Unlike the AutoDock-GPU backend no affinity
    maps files are involved: vina computes its own maps internally from the rigid part.
    
    NOTICE: vina performs the search on CPU and flexible sidechains add many degrees
    of freedom, so this is much slower than `perform_autodock_gpu_flex_dock`.
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box padding in ANGSTROM, added on every side of the box
        that covers the ligand and every flexible residue
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param flex_residues: explicit flexible residues
    :param flex_distance: distance cutoff, select residues within this distance of ligand heavy atoms
    :param flex_interaction: FLAG, whether to select flexible residues by fplip interaction
    :param flex_interaction_mode: interaction types for fplip, 'all' to include all modes
    :param flex_interaction_cutoff: interaction distance cutoff for fplip
    :param exhaustiveness: exhaustiveness, default is 32
    :param n_poses: n_poses
    :param min_rmsd: min_rmsd
    :param max_evals: max_evals, 0 means use vina heuristics
    :param seed: random seed
    :param cpu: cpu number for vina, 0 means use all cpus
    :param energy_range: max energy difference (kcal/mol) from the best pose for
        a pose to be written; None keeps vina's default (3.0), which silently
        drops poses (see `perform_vina_dock`)
    :param spacing: grid spacing in ANGSTROM, default 0.375
    :param macrocycle: prepare the ligand with meeko, breaking a large ring so that
        its conformation can adapt during the search (requires vina >= 1.2).  Off by
        default so that existing results stay reproducible.
    :param max_ring_size: largest ring meeko may break, see MACROCYCLE_MAX_RING_SIZE
    :param prepare_backend: 'meeko' or 'adfr', which toolkit prepares the LIGAND pdbqt.
        The receptor is always split by meeko here, because ADT's `prepare_receptor`
        can not produce a rigid part plus flexible sidechains at all, so this switch
        only ever changes the ligand.  Use 'adfr' only when the ligand must match an
        ADT-prepared run: it makes ADFRsuite's `prepare_ligand` a hard requirement, and
        the two backends are not score-comparable (see `PREPARE_BACKENDS`)
    
    :return: (result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr),
        None if failed. `result_str` is a vina-format pdbqt string (NOT a dlg),
        in which every MODEL contains the ligand AND the flexible sidechains.
        A pose string can be turned into a complex by `build_flex_pose_complex_pdbstr`,
        or parsed by `lazydock.pml.autodock_utils.DlgFile` (it understands the vina
        `REMARK VINA RESULT` score line).
    """
    from vina import Vina
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor pdb string, and ligand pdbqt
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain,
                                               macrocycle, max_ring_size, prepare_backend)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, _, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        ligand_path = os.path.join(tmpdir, 'ligand.pdbqt')
        opts_file(ligand_path, 'w', way='str', data=ligand_pdbqtstr)
        # flexible residues, ADFRsuite bin dir, grid box (covering the ligand and
        # every flexible residue) and the rigid/flex receptor split; shared with
        # the AutoDock-GPU flex path (see `prepare_flex_docking_setup`)
        setup = prepare_flex_docking_setup(
            tmpdir, complex_pdbstr, ligand_chain, receptor_chain,
            receptor_pdbstr, ligand_pdbstr, grid_buffer,
            flex_residues, flex_distance, flex_interaction,
            flex_interaction_mode, flex_interaction_cutoff)
        if not setup:
            return None
        rigid_pdbqt_path, flex_pdbqt_path = setup.rigid_pdbqt_path, setup.flex_pdbqt_path
        # vina scores the flexible sidechains itself, no gpf/autogrid needed.
        # grid_center/grid_size are in ANGSTROM, which is what compute_vina_maps
        # wants; `spacing` is the separate knob that defines the voxel size
        v = Vina(sf_name=score_name, cpu=cpu, seed=seed, verbosity=False)
        v.set_receptor(rigid_pdbqt_path, flex_pdbqt_path)
        v.set_ligand_from_file(ligand_path)
        v.compute_vina_maps(center=setup.grid_center, box_size=setup.grid_size, spacing=spacing)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses,
               min_rmsd=min_rmsd, max_evals=max_evals)
        output_path = os.path.join(tmpdir, 'poses.pdbqt')
        # energy_range default (3.0 kcal/mol) silently drops poses, see
        # `perform_vina_dock`; the flex path must honour n_poses as well
        v.write_poses(output_path, n_poses, overwrite=True,
                      **({} if energy_range is None else {'energy_range': float(energy_range)}))
        result_str = opts_file(output_path, 'r', way='str')
        rigid_pdbqtstr = opts_file(rigid_pdbqt_path, 'r', way='str')
        flex_pdbqtstr = opts_file(flex_pdbqt_path, 'r', way='str')
        return result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr