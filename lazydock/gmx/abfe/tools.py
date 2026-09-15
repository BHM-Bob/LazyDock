'''
Date: 2026-09-01
Description: utilities for ABFE run engine, migrated from BindFlow bindflow.utils.tools.
             The gmx command layer is replaced with LazyDock's Gromacs class
             (expect-based, no subprocess piping).
             Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
import os
import re
from math import sqrt
from pathlib import Path
from typing import Iterable, List, Optional, Union

PathLike = Union[os.PathLike, str, bytes]


def paths_exist(paths: List, raise_error: bool = False, out: Union[str, None] = None) -> None:
    """Check that the paths exist.

    If ``out`` is provided and all paths exist, an empty marker file is created
    (used to mark a step as finished).
    """
    check = True
    for path in paths:
        if not Path(path).exists():
            check = False
            msg = f"Missing path/file: {path}"
            if raise_error:
                raise RuntimeError(msg)
            else:
                print(msg)
    if out and check:
        open(out, "w").close()


def list_if_dir(path: PathLike = '.') -> List[Path]:
    return [p for p in Path(path).iterdir() if p.is_dir()]


def list_if_file(path: PathLike = '.', ext: str = None) -> List[Path]:
    files = [p for p in Path(path).iterdir() if p.is_file()]
    if ext:
        files = [file for file in files if file.suffix == ext]
    return files


def sum_uncertainty_propagation(errors: Iterable[float], coefficients: Optional[Iterable[float]] = None) -> float:
    """sqrt( sum (c_i * sigma_i)^2 ) - standard uncertainty propagation."""
    if coefficients is None:
        coefficients = [1.0] * len(errors)
    else:
        coefficients = list(coefficients)
    if len(coefficients) != len(errors):
        raise ValueError("`coefficients` must have the same length as `errors`.")
    return sqrt(sum((c * e) ** 2 for c, e in zip(coefficients, errors)))


def _gmx_for_run_dir(gmx, run_dir: Path) -> 'Gromacs':
    """Return a NEW Gromacs object whose working_dir == run_dir.

    The shared ``gmx`` (created once at the output root, carrying the GPU
    selection) is only used as the ``gpu_ids`` source here; it is NEVER
    returned/reused directly: ``run_gmx_with_expect`` prepends
    ``cd <working_dir>``, so mdrun's relative outputs (-deffnm) would land in
    the output root instead of the simulation step directory.  Always creating
    a fresh per-run-dir instance also keeps the parallel task pool free of any
    shared-object state (e.g. caching) that a future Gromacs upgrade might
    introduce.  BindFlow's ``gmx_runner`` does ``os.chdir(run_dir)`` before
    launching gmx; we emulate that by building a per-run-dir Gromacs that
    inherits the gpu_ids.
    """
    from lazydock.gmx.run import Gromacs
    run_dir = Path(run_dir)
    gpu_ids = list(gmx.gpu_ids) if (gmx is not None and getattr(gmx, 'gpu_ids', None)) else [0]
    return Gromacs(working_dir=str(run_dir), gpu_ids=gpu_ids)


def gmx_runner(gmx, mdp: PathLike, topology: PathLike, structure: PathLike, checkpoint: PathLike = None,
               index: PathLike = None, nthreads: int = 12, run_dir: PathLike = '.',
               maxwarn: int = 2, minimize: bool = False, gpu_id: Optional[int] = None,
               ntmpi: int = None, **mdrun_extra):
    """Create the tpr file from mdp/structure/topology and run mdrun (LazyDock Gromacs class).

    Default commands (same as BindFlow ``tools.gmx_runner``):
        gmx grompp -f {mdp} -c {structure} -r {structure} -p {topology} -o {name}.tpr
                   (-t checkpoint) (-n index) -maxwarn {maxwarn}
        gmx mdrun -nt {nthreads} -deffnm {name} (-cpi) ...mdrun_extra

    ``ntmpi``: number of MPI ranks (mdrun -ntmpi).  Default None -> GROMACS
    auto domain decomposition (can fail with large rlist on small boxes:
    "no domain decomposition for N ranks ... minimum cell size").  Peptide
    FEP windows with rlist ~3 nm need ntmpi=1.

    All gmx invocations run inside ``run_dir`` (a per-step Gromacs object is
    created if the shared one points elsewhere).  ``gpu_id`` selects the GPU
    slot of the Gromacs object (only used when the object was created with
    several gpu_ids).
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(exist_ok=True, parents=True)
    _gmx = _gmx_for_run_dir(gmx, run_dir)
    name = Path(mdp).stem

    grompp_kwargs = dict(f=str(mdp), c=str(structure), r=str(structure),
                         p=str(topology), o=f"{name}.tpr", maxwarn=maxwarn)
    if checkpoint:
        grompp_kwargs['t'] = str(checkpoint)
    if index:
        grompp_kwargs['n'] = str(index)
    ret = _gmx.run_gmx_with_expect('grompp', **grompp_kwargs)
    if ret != 0:
        raise RuntimeError(f"gmx grompp failed (exit {ret}) for {mdp}")

    mdrun_kwargs = dict(nt=nthreads, deffnm=name)
    if ntmpi is not None:
        mdrun_kwargs['ntmpi'] = ntmpi
    if mdrun_extra:
        mdrun_kwargs.update(mdrun_extra)
    if gpu_id is not None:
        mdrun_kwargs['gpu_id'] = gpu_id
    ret = _gmx.run_gmx_with_expect('mdrun', **mdrun_kwargs)
    if ret != 0:
        raise RuntimeError(f"gmx mdrun failed (exit {ret}) for {name}")


def center_xtc(tpr: PathLike, xtc: PathLike, run_dir: PathLike, host_name: str = 'Protein',
               gmx=None) -> PathLike:
    """Center an xtc file: whole -> nojump -> center (pbc mol, ur compact).

    The third trjconv call asks for two groups: the group to center (host,
    e.g. 'Protein') and the group to output ('System'); the expect script
    answers with the two names.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(exist_ok=True, parents=True)
    _gmx = _gmx_for_run_dir(gmx, run_dir)

    # whole
    _gmx.run_gmx_with_expect('trjconv', s=str(tpr), f=str(xtc),
                             o=str(run_dir / "whole.xtc"), pbc="whole",
                             expect_actions=[{'Select group': 'System\r'}, {'>': 'q\r'}])
    # nojump
    _gmx.run_gmx_with_expect('trjconv', s=str(tpr), f=str(run_dir / "whole.xtc"),
                             o=str(run_dir / "nojump.xtc"), pbc="nojump",
                             expect_actions=[{'Select group': 'System\r'}, {'>': 'q\r'}])
    # center (ask for group to center = host_name, group to output = System)
    _gmx.run_gmx_with_expect('trjconv', s=str(tpr), f=str(run_dir / "nojump.xtc"),
                             o=str(run_dir / "center.xtc"), pbc="mol", center=True, ur="compact",
                             expect_actions=[{'Select group': f'{host_name}\r'},
                                             {'Select group': 'System\r'}])
    # Clean
    (run_dir / "whole.xtc").unlink()
    (run_dir / "nojump.xtc").unlink()

    return f"{run_dir}/center.xtc"


if __name__ == '__main__':
    pass