'''
Date: 2026-09-01
Description: run ABFE (absolute binding free energy) simulations + analysis.
             Consumes the input/ directory produced by prepare-abfe.
             Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
import argparse
from pathlib import Path
from typing import List

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension

from lazydock.gmx.run import Gromacs
from lazydock.scripts._script_utils_ import Command, make_args_and_excute, process_batch_dir_lst


def _parse_mdp_extra(s: str) -> dict:
    """Parse a user mdp override string like:
       'ligand.equi.00_min.nsteps=1000;complex.fep.prod.nsteps=100'
    into the nested dict used by the engine config['mdp'].
    """
    result = {}
    if not s:
        return result
    for token in s.split(';'):
        token = token.strip()
        if not token:
            continue
        if '=' not in token:
            raise ValueError(f"Invalid mdp override token: {token!r} (expected a=b=c=value)")
        key_path, value = token.rsplit('=', 1)
        parts = [p for p in key_path.split('.') if p]
        if len(parts) < 3:
            raise ValueError(f"Invalid mdp override key: {key_path!r} "
                             "(expected <system>.<stage>.<step>.<param>, e.g. ligand.equi.00_min.nsteps)")
        d = result
        for p in parts[:-1]:
            d = d.setdefault(p, {})
        d[parts[-1]] = value
    return result


class _RunBase(Command):
    """Common base of the equi/fep/run subcommands: the shared argument
    table, the batch discovery (scan the prepare-abfe convention file
    complex.pdb, -n to override) and the engine config builder.
    Not registered as a subcommand itself.
    """
    HELP = """(base class, not a subcommand)"""

    def __init__(self, args, printf=print):
        super().__init__(args, printf)

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, nargs='+', default=['.'],
                          help='batch directory containing complex.pdb subdirs '
                               '(or a single case directory with input/). Default %(default)s.')
        args.add_argument('-n', '--name', type=str, default='complex.pdb',
                          help='complex structure file name in each case dir; '
                               'must match prepare-abfe -n/--name. Default %(default)s.')
        args.add_argument('--replicas', type=int, default=1, help='number of replicas. Default %(default)s.')
        args.add_argument('--threads', type=int, default=12, help='mdrun -nt. Default %(default)s.')
        args.add_argument('--ntmpi', type=int, default=None,
                          help='mdrun -ntmpi (MPI ranks). For peptide FEP windows with large '
                               'rlist (e.g. 3.1 nm) the box is too small for multi-rank domain '
                               'decomposition: set --ntmpi 1 (auto DD may fatal with '
                               '"no domain decomposition for N ranks"). Default %(default)s.')
        args.add_argument('--gpus', type=int, nargs='+', default=None,
                          help='gpu ids for mdrun (Gromacs gpu_ids). Default %(default)s.')
        args.add_argument('--mdrun-args', type=str, default='',
                          help='extra mdrun kwargs as "key=value key2=value2" '
                               '(e.g. "update=cpu pme=gpu"). Flags: "flag=1" -> bare flag. '
                               'Passed to every mdrun. Default empty.')
        args.add_argument('--maxwarn', type=int, default=2, help='grompp maxwarn. Default %(default)s.')
        args.add_argument('--dt-max', type=float, default=0.008,
                          help='max dt (ns) for mdp templates. Default %(default)s.')
        args.add_argument('--equi-integrator', type=str, default='md', choices=['md', 'sd'],
                          help='integrator for equilibration MD steps (01_nvt/02_nvt/03_npt/'
                               '04_npt/prod). Default md (v-rescale thermostat; faster, better '
                               'GPU offload). sd (Langevin, template default) is slower with '
                               'poor GPU offload; keep md for production. FEP windows always '
                               'use sd regardless (required for rigorous canonical sampling).')
        args.add_argument('--dt-max-fep', type=float, default=0.004,
                          help='max dt (ns) for FEP step mdp templates. BindFlow default 0.004 '
                               '(4 fs, relies on HMR+all-bonds); use 0.002 for safer 2 fs FEP '
                               'windows. Default %(default)s.')
        args.add_argument('--fep-rlist', type=float, default=None,
                          help='rlist (nm) for FEP windows in peptide mode (couple-intramol=no). '
                               'Default 3.1 (covers cyclic-peptide diameter ~2.6 nm + margin; '
                               'must be < half shortest box vector).')
        args.add_argument('--fep-rlist-ligand', type=float, default=None,
                          help='rlist (nm) for the LIGAND-leg FEP windows only (overrides '
                               '--fep-rlist for the ligand leg). Small-molecule ligands '
                               '(couple-intramol=yes) need only the template 1.2; large '
                               'non-standard amino-acid peptides processed through the '
                               'small-molecule path may need more. Peptide ligands: do NOT '
                               'lower below 3.1 (peptide can stretch in solvent).')
        args.add_argument('--fep-cutoff', type=float, default=None,
                          help='rcoulomb/rvdw (nm) for FEP windows in peptide mode '
                               '(couple-intramol=no). grompp checks non-perturbed excluded '
                               'pairs against the max(rlist, rvdw, rcoulomb); flexible/cyclic '
                               'peptide 1-4 pairs can stretch beyond the standard 1.0 nm, so '
                               'these are raised together. Default 1.4 when --peptide, else '
                               'template value (1.0). Pass an explicit value to override.')
        args.add_argument('--nwindows', type=str, default=None,
                          help='window counts override as "ligand.vdw=N;ligand.coul=N;'
                               'complex.vdw=N;complex.coul=N;complex.bonded=N" (default BindFlow: '
                               '11/11/21/11/11).')
        args.add_argument('--mdp-extra', type=str, default=None,
                          help='user mdp overrides (semicolon list of dotted keys): '
                               '"ligand.equi.00_min.nsteps=100;complex.fep.prod.nsteps=100".')
        args.add_argument('--fep-prod-ns', type=float, default=None,
                          help='set the FEP production (prod) step length in ns for ALL FEP '
                               'windows (ligand vdw/coul + complex vdw/coul/bonded). '
                               'Overrides the template nsteps; e.g. --fep-prod-ns 1 = 1 ns '
                               '(250000 steps @ 4 fs) per window. Use a short value (0.5-1) '
                               'to validate the pipeline quickly, or omit for the full 10 ns '
                               'template length. Takes precedence over the template but is '
                               'overridden by an explicit --mdp-extra nsteps for the same step.')
        args.add_argument('--equi-prod-ns', type=float, default=None,
                          help='set the equilibration prod step length in ns (ligand and '
                               'complex equil-mdsim). Overrides the template nsteps; e.g. '
                               '--equi-prod-ns 1 = 1 ns. Default keeps the template length.')
        args.add_argument('--n-parallel', type=int, default=1,
                          help='number of concurrent tasks/windows. For FEP: FEP windows to '
                               'run concurrently (two-level scheduler: a pool of window workers, '
                               'each running its full step chain 00_min->...->prod on one GPU). '
                               'For equil: cross-run equilibration tasks (case x leg x replica) '
                               'to run concurrently. Windows/tasks are distributed across --gpus '
                               'round-robin. Default 1 (sequential, BindFlow behavior). '
                               'Use e.g. --n-parallel 8 with 4 GPUs (2 windows per GPU, '
                               'nt per window = threads/n-parallel*ngpus ... see docs).')
        args.add_argument('--retries', type=int, default=3, help='gmx retries per step. Default %(default)s.')
        args.add_argument('--only-build', action='store_true', default=False,
                          help='only create directories + mdp structure, no simulation.')
        args.add_argument('--host-name', type=str, default='Protein',
                          help='host group name for trjconv centering. Default %(default)s.')
        args.add_argument('--host-selection', type=str, default='protein and name CA and not moltype LIG',
                          help='MDAnalysis selection for Boresch host (default excludes peptide ligand moltype). Default %(default)s.')
        args.add_argument('--peptide', action='store_true', default=False,
                          help='peptide ligand (flexible chain): use couple-intramol=no in FEP mdps '
                               'to prevent unfolding in decoupled windows. '
                               'Must match prepare-abfe --peptide. Default %(default)s.')
        return args

    def process_args(self):
        self.args.dir = process_batch_dir_lst(self.args.dir)

    def _build_global_config(self, wdir: Path) -> dict:
        """Assemble the engine config from CLI args (BindFlow global_config schema)."""
        input_dir = wdir / 'input'
        if not input_dir.exists():
            put_err(f'input/ not found in {wdir}; run prepare-abfe first.', _exit=True)
        complex_dir, ligand_dir = input_dir / 'complex', input_dir / 'ligand'
        if not (complex_dir / 'complex.gro').exists():
            put_err(f'{complex_dir}/complex.gro not found.', _exit=True)
        if not (ligand_dir / 'ligand.top').exists() and not (ligand_dir / 'ligand.gro').exists():
            put_err(f'ligand files not found in {ligand_dir}.', _exit=True)

        # out_root is the working directory itself (LazyDock convention
        # 2026-09-14: no separate --out, no {ligand_name}/ layer)
        out_root = wdir
        ligand_pdb = ligand_dir / 'ligand.pdb'
        if not ligand_pdb.exists():
            ligand_pdb = ligand_dir / 'ligand.gro'

        # nwindows override
        nwindows = None
        if self.args.nwindows:
            nwindows = {}
            for token in self.args.nwindows.split(';'):
                token = token.strip()
                if not token:
                    continue
                key, value = token.split('=')
                k1, k2 = key.strip().split('.')
                nwindows.setdefault(k1, {})[k2] = int(value.strip())

        # mdrun extra args ("key=value flag=1")
        mdrun_extra = {}
        if self.args.mdrun_args.strip():
            for token in self.args.mdrun_args.split():
                if '=' not in token:
                    raise ValueError(f'Invalid --mdrun-args token: {token!r}')
                key, value = token.split('=', 1)
                if value in ('1', 'true', 'True'):
                    mdrun_extra[key] = True
                elif value in ('0', 'false', 'False'):
                    mdrun_extra[key] = False
                else:
                    mdrun_extra[key] = value
        # BindFlow always runs mdrun with -cpi (checkpoint continuation); this
        # is harmless when no checkpoint exists yet and enables resuming after
        # a crash.  Users can disable it with --mdrun-args "cpi=0".
        mdrun_extra.setdefault('cpi', True)

        # Build mdp overrides: user --mdp-extra takes precedence over the
        # convenient --fep-prod-ns / --equi-prod-ns shortcuts (shortcuts only
        # fill in when the user did not already set nsteps for that step).
        mdp_overrides = _parse_mdp_extra(self.args.mdp_extra)
        dt_fep = self.args.dt_max_fep if self.args.dt_max_fep else self.args.dt_max
        if self.args.fep_prod_ns is not None:
            nsteps = int(round(self.args.fep_prod_ns * 1000 / dt_fep))  # ns -> steps @ dt(ns)
            # ligand: vdw/coul; complex: vdw/coul/bonded
            for sys_type, lam_types in (('ligand', ['vdw', 'coul']),
                                        ('complex', ['vdw', 'coul', 'bonded'])):
                for lam in lam_types:
                    prod = mdp_overrides.setdefault(sys_type, {}).setdefault('fep', {}) \
                        .setdefault(lam, {}).setdefault('prod', {})
                    prod.setdefault('nsteps', str(nsteps))
        if self.args.equi_prod_ns is not None:
            # equil 各步模板 dt 均 <= 0.004 (配 HMR), 且 _make_equi_mdps 只把 dt
            # 钳制到 <= dt_max, 即实际生效 dt = min(dt_max, 0.004). 换算基准必须用
            # 该实际生效 dt, 否则 --dt-max 默认 0.008 时 10ns 会被算成 5ns.
            eff_dt = min(self.args.dt_max, 0.004)
            nsteps = int(round(self.args.equi_prod_ns * 1000 / eff_dt))
            for sys_type in ('ligand', 'complex'):
                prod = mdp_overrides.setdefault(sys_type, {}).setdefault('equi', {}) \
                    .setdefault('prod', {})
                prod.setdefault('nsteps', str(nsteps))

        gpus = self.args.gpus or [0]
        gmx = Gromacs(working_dir=str(out_root), gpu_ids=gpus)

        config = {
            'calculation_type': 'fep',
            'out_approach_path': str(out_root),
            'input_dir': str(input_dir),
            'inputs': {
                'protein': {'conf': str(wdir / 'receptor.pdb')} if (wdir / 'receptor.pdb').exists()
                else {'conf': str(complex_dir / 'complex.gro')},
                'ligands': [{'conf': str(ligand_pdb)}],
            },
            'replicas': self.args.replicas,
            'water_model': None,  # not used in run stage
            'host_name': self.args.host_name,
            'host_selection': self.args.host_selection,
            'hmr_factor': None,
            'custom_ff_path': None,
            'fix_protein': False,
            'solv_d': 1.5, 'solv_bt': 'dodecahedron', 'solv_rmin': 1, 'solv_ion_conc': 150e-3,
            'threads': self.args.threads,
            'ntmpi': self.args.ntmpi,
            'dt_max': self.args.dt_max,
            'equi_integrator': self.args.equi_integrator,
            'fep_dt_max': self.args.dt_max_fep,
            'fep_rlist': self.args.fep_rlist,
            'fep_rlist_ligand': self.args.fep_rlist_ligand,
            'fep_cutoff': self.args.fep_cutoff,
            'retries': self.args.retries,
            'maxwarn': self.args.maxwarn,
            # Serial stages (equil / boresch / sequential FEP) use the first
            # available GPU explicitly, instead of letting mdrun pick the
            # default device. Parallel FEP workers override gpu_id per window.
            'gpu_id': gpus[0],
            'gpu_ids': gpus,
            'n_parallel': self.args.n_parallel,
            'nwindows': nwindows,
            'mdp': mdp_overrides,
            'extra_directives': {
                'mdrun': {'ligand': dict(mdrun_extra), 'complex': dict(mdrun_extra)},
            },
            'gmx': gmx,
            'peptide_ligand': self.args.peptide,
        }
        return config

    def _batch_configs(self) -> 'list[dict]':
        """Scan the batch dir(s) for the ABFE cases: any directory containing
        the same convention file prepare-abfe consumed (complex.pdb by default,
        see prepare-abfe -n/--name) whose input/ exists.  Also accepts a
        directory that directly contains input/ (a single case).  Returns one
        engine config per case.
        """
        cases: list[Path] = []
        for d in self.args.dir:
            d = Path(d)
            if (d / 'input').exists():
                cases.append(d)
                continue
            if not d.is_dir():
                put_err(f'dir argument should be a directory: {d}', _exit=True)
            found = get_paths_with_extension(d, [], name_substr=self.args.name, sort='natsort')
            if not found:
                put_log(f'no {self.args.name} found under {d} - nothing to do', head='ABFE')
                continue
            seen = set()
            for p in found:
                case_dir = Path(p).parent
                if case_dir in seen:
                    continue
                seen.add(case_dir)
                if (case_dir / 'input').exists():
                    cases.append(case_dir)
                else:
                    put_log(f'{case_dir}: no input/ (prepare-abfe not run?) - skipped',
                            head='ABFE')
        if not cases:
            put_err(f'no ABFE case ({self.args.name} + input/) found under {self.args.dir}',
                    _exit=True)
        return [self._build_global_config(c) for c in cases]

    @staticmethod
    def _get_engine():
        from lazydock.gmx.abfe import engine
        return engine


class Equil(_RunBase):
    """Run only the equilibration stage: cross-run task pool + Boresch (serial).

    Reuses the shared argument table, batch scan and config builder of _RunBase.
    """
    HELP = """
    run only the ABFE equilibration stage: cross-run task pool (task =
    (case, leg, replica)) -> Boresch restraints (serial).

    INPUT:
        same batch directory convention as `run` (cases = subdirs with
        complex.pdb and input/).
    """
    def main_process(self):
        configs = self._batch_configs()
        engine_ = self._get_engine()
        engine_.run_abfe_equil(configs, only_build=self.args.only_build)
        for cfg in configs:
            put_log(f'ABFE equilibration completed. Output in: {cfg["out_approach_path"]}')


class Fep(_RunBase):
    """Run only the FEP stage: windows -> analysis -> gather.

    Reuses the shared argument table, batch scan and config builder of _RunBase.
    """
    HELP = """
    run only the ABFE FEP stage: check Boresch restraints -> window
    simulation (parallel or serial) -> analysis -> gather.

    INPUT:
        same batch directory convention as `run` (cases = subdirs with
        complex.pdb and input/); requires `run-abfe equil` (or `run`) first.
    """
    def main_process(self):
        configs = self._batch_configs()
        engine_ = self._get_engine()
        engine_.run_abfe_fep(configs, only_build=self.args.only_build)
        for cfg in configs:
            put_log(f'ABFE FEP completed. Output in: {cfg["out_approach_path"]}')


class Status(Command):
    HELP = """
    show the progress of an ABFE run: per-leg (ligand/complex x vdw/coul/bonded)
    finished windows / total, running count and elapsed wall time.
    Read-only: only scans .finished markers and mdp/log metadata.

    INPUT:
        the working directory of an abfe run ({wdir}/replica_{N}/...),
        or a replica directory itself.
    """
    def __init__(self, args, printf=print):
        # NOTE: must NOT use iter_run_arg=['dir']: the generic Command.excute()
        # replaces args.dir with a bare string (one value per iteration), which
        # breaks `for d in dir` in main_process (iterates characters instead).
        super().__init__(args, printf, iter_run_arg=[])

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, nargs='+', default=['.'],
                          help='ABFE run output directory (default %(default)s).')

    def main_process(self):
        from lazydock.gmx.abfe.status import render_status, scan_abfe_progress
        for d in self.args.dir:
            rows = scan_abfe_progress(d)
            self.printf(render_status(rows))


class Run(_RunBase):
    """Full ABFE workflow subcommand (placed last: it simply composes the
    equil -> fep stages described above).
    """
    HELP = """
    run ABFE: equilibration -> Boresch restraints -> FEP simulations -> analysis -> gather.

    INPUT:
        a batch directory containing the prepare-abfe output layout:
            <case>/{input/complex, input/ligand}   (one case per subdir with complex.pdb)
        The engine builds the LazyDock layout in the SAME case directory (no
        copied input/, no {ligand_name}/ layer):
            <case>/replica_{N}/{ligand,complex}/(equil-mdsim|fep)
    """
    def main_process(self):
        configs = self._batch_configs()
        engine_ = self._get_engine()
        for cfg in configs:
            put_log(f'processing ABFE run in: {cfg["out_approach_path"]}', head='ABFE')
        if len(configs) == 1:
            engine_.run_abfe(configs[0], only_build=self.args.only_build)
        else:
            engine_.run_abfe(configs, only_build=self.args.only_build)
        for cfg in configs:
            put_log(f'ABFE run completed. Output in: {cfg["out_approach_path"]}')


_str2func = {
    'run': Run,
    'equil': Equil,
    'fep': Fep,
    'status': Status,
}


def main(sys_args: List[str] = None):
    make_args_and_excute('tools for ABFE (absolute binding free energy) run.', _str2func, sys_args)


if __name__ == '__main__':
    main()