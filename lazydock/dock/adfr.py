import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import opts_file
from pymol import cmd

from lazydock.config import GlobalConfig
from lazydock.dock.config import (MACROCYCLE_MAX_RING_SIZE,
                                  PREPARE_BACKEND_DEFAULT, PREPARE_BACKENDS)
from lazydock.dock.meeko_utils import (prepare_ligand_with_meeko,
                                       prepare_receptor_with_meeko)
from lazydock.pml.utils import normalize_receptor_chain


_ADFR_BIN_MISS_LOGGED = False


def prepend_adfr_bin_to_path(gpf_caller: Optional[str] = None) -> Optional[str]:
    """
    resolve prepare_gpf.py and prepend its dir (ADFRsuite bin) to $PATH.

    ADFRsuite's `prepare_ligand`, `prepare_receptor`, `pythonsh` and `autogrid4`
    are all called by bare name, so their dir has to be on $PATH.  This function
    is idempotent: a dir already present is left alone, and one that is not is
    prepended once.

    When nothing can be resolved the function stays silent after the FIRST miss
    (see `_ADFR_BIN_MISS_LOGGED`), because every caller may reach this on a normal
    meeko-only run.  It does not raise: the caller is expected to fall back on the
    user's own $PATH, and the failure (if any) surfaces as the shell's own
    "command not found" when the preparation is actually attempted.

    :param gpf_caller: prepare_gpf.py path, or the ADFRsuite bin dir itself; None
        means read from GlobalConfig.named_paths["prepare_gpf"] or $PREPARE_GPF
    :return: the resolved absolute prepare_gpf.py path, None if none was resolved
    """
    global _ADFR_BIN_MISS_LOGGED
    if not gpf_caller:
        gpf_caller = GlobalConfig.named_paths["prepare_gpf"] or os.environ.get("PREPARE_GPF")
    if not gpf_caller:
        if not _ADFR_BIN_MISS_LOGGED:
            _ADFR_BIN_MISS_LOGGED = True
            put_log('prepare_gpf is not configured and can not be found with '
                    '~/.lazydock/lazydock_config.json or $PREPARE_GPF; assuming the '
                    'ADFRsuite bin dir is already on $PATH (this is logged once)')
        return None
    gpf_caller = os.path.expanduser(gpf_caller)
    bin_dir = os.path.dirname(os.path.abspath(gpf_caller))
    if bin_dir in os.environ.get('PATH', '').split(os.pathsep):
        return gpf_caller
    if not Path(gpf_caller).exists():
        return put_err(f'given gpf_caller does not exist: {gpf_caller}')
    if not Path(bin_dir).is_dir():
        return put_err(f'given bin_dir {bin_dir} is not a directory')
    os.environ['PATH'] = f'{bin_dir}{os.pathsep}{os.environ.get("PATH", "")}'
    return gpf_caller


def get_pdbqtstr_from_pdbstr(pdbstr: str, caller: str = 'prepare_ligand', entry_param: str = '-l',
                             options: List[Tuple['str', str]] = None):
    """
    call ADT to prepare a ligand or receptor pdbqt string.

    This is the ADT (AutoDockTools) route: `prepare_ligand` / `prepare_receptor` are
    external executables from ADFRsuite, so they are called through a shell by bare
    name and `$PATH` must contain the ADFRsuite bin dir (see
    `prepend_adfr_bin_to_path`).  meeko, the default backend, does the same job
    in-process and never needs this function.

    The ADT- and meeko-prepared pdbqt are NOT score-equivalent: the two type a few
    nitrogen acceptors and polar hydrogens differently (see `PREPARE_BACKENDS`).

    :param pdbstr: pdb string of the molecule to prepare
    :param caller: executable name, 'prepare_ligand' or 'prepare_receptor', or an
        absolute path to it
    :param entry_param: flag that tells `caller` which kind of input it gets,
        '-l' for ligand and '-r' for receptor
    :param options: extra options for `caller`, as (flag, value) pairs.  None means
        [('-A', 'hydrogens')], i.e. add hydrogens (a bare '-A' flag with no value)
    :return: (FLAG, pdbqt string), FLAG is True when the executable succeeded AND
        produced an output file.  On failure FLAG is False and the second item is the
        combined stdout+stderr, not a pdbqt
    """
    # ADFRsuite's prepare_ligand / prepare_receptor are called by bare name, so the
    # bin dir must be on $PATH.  A failure is logged once and otherwise ignored: the
    # user may have exported it himself, and if not, the shell reports the missing
    # command in `output_str` below.
    prepend_adfr_bin_to_path()
    
    # prepare command line options
    options = options or [('-A', 'hydrogens')]
    cmd_options = []
    [cmd_options.extend(list(i)) for i in options]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        opts_file(os.path.join(tmpdir, 'mol.pdb'), 'w', way='str', data=pdbstr)
        cmd = f"cd {tmpdir} && {caller} {entry_param} mol.pdb -o mol.pdbqt {' '.join(cmd_options)}"
        proc = subprocess.Popen(cmd, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = proc.communicate()
        output_str = stdout + stderr
        if proc.returncode != 0 or not Path(f'{tmpdir}/mol.pdbqt').exists():
            return False, output_str
        return True, opts_file(f'{tmpdir}/mol.pdbqt', way='str')


def prepare_receptor_pdbqtstr(receptor_pdbstr: str, prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """
    prepare a rigid receptor pdbqt string with the requested backend.

    Both backends add the polar hydrogens themselves, so a hydrogen-free pdb string is
    accepted.  The two are NOT interchangeable score-wise (see `PREPARE_BACKENDS` in
    `dock.config`).

    :param receptor_pdbstr: receptor pdb string
    :param prepare_backend: 'meeko' (in-process, the default) or 'adfr' (external
        `prepare_receptor`, requires ADFRsuite on $PATH -- see
        `prepend_adfr_bin_to_path`)
    :return: pdbqt string, None if failed
    """
    if prepare_backend == 'meeko':
        prepared = prepare_receptor_with_meeko(receptor_pdbstr)
        if not prepared:
            return None
        # rigid-only preparation never yields a flexible part
        return prepared[0]
    if prepare_backend == 'adfr':
        success, receptor_pdbqtstr = get_pdbqtstr_from_pdbstr(receptor_pdbstr, 'prepare_receptor', '-r')
        if not success:
            return put_err(f"Failed to prepare receptor_pdbqtstr: {receptor_pdbqtstr}")
        return receptor_pdbqtstr
    return put_err(f'unknown prepare_backend "{prepare_backend}", '
                   f'expected one of {PREPARE_BACKENDS}')


def prepare_ligand_pdbqtstr(ligand_pdbstr: str, prepare_backend: str = PREPARE_BACKEND_DEFAULT,
                            macrocycle: bool = False,
                            max_ring_size: int = MACROCYCLE_MAX_RING_SIZE):
    """
    prepare a ligand pdbqt string with the requested backend.

    the twin of `prepare_receptor_pdbqtstr`: same contract (a pdbqt string back, None
    on failure after reporting the error), but the ligand has one extra degree of
    freedom -- whether a large ring is broken open.

    Both backends add the hydrogens themselves.  The two are NOT interchangeable
    score-wise (see `PREPARE_BACKENDS` in `dock.config`).

    :param ligand_pdbstr: ligand pdb string
    :param prepare_backend: 'meeko' (in-process, the default) or 'adfr' (external
        `prepare_ligand`, requires ADFRsuite on $PATH).  IGNORED when `macrocycle` is
        True, because meeko is the only backend that can break a ring (see `macrocycle`)
    :param macrocycle: break one bond of a large ring and attach glue pseudo-atoms,
        which vina >= 1.2 pulls back closed during the search, so the ring becomes a
        normal torsion tree and its conformation can adapt to the pocket
    :param max_ring_size: largest ring meeko may break, ignored when `macrocycle` is
        False
    :return: pdbqt string, None if failed
    """
    if macrocycle:
        # ADT's prepare_ligand has NO ring/macrocycle option at all (its option list
        # holds nothing matching ring|cycl|glue|macro), so a cyclic ligand prepared
        # by it has every ring bond frozen -- the ring conformation is then literally
        # unable to move during docking.  meeko breaks one ring bond and adds glue
        # pseudo-atoms instead, which vina >= 1.2 pulls back closed during the search.
        prepared = prepare_ligand_with_meeko(ligand_pdbstr, max_ring_size)
        if not prepared:
            return None
        ligand_pdbqtstr, macro_info = prepared
        if not macro_info['ring_flexible']:
            # meeko silently ignores rings larger than `max_ring_size`, so this is
            # the only signal that the requested flexibility was not delivered
            put_err(f'macrocycle was requested but no ring was opened '
                    f'(max_ring_size={max_ring_size}, the ring may be larger); '
                    f'docking would silently run with a rigid ring.')
        else:
            put_log(f'macrocycle prepared by meeko: {macro_info["n_glue_pairs"]} ring bond(s) '
                    f'opened, {macro_info["n_branch"]} BRANCH records, '
                    f'TORSDOF={macro_info["torsdof"]}')
        return ligand_pdbqtstr
    if prepare_backend == 'meeko':
        # `break_rings=False` is mandatory here: meeko's ring breaking is on by
        # default over the whole 7..33 atom window, so a plain medium ring would be
        # opened and the resulting CG0/G0 atoms have no grid in an AutoDock map
        prepared = prepare_ligand_with_meeko(ligand_pdbstr, break_rings=False)
        if not prepared:
            return None
        return prepared[0]
    elif prepare_backend == 'adfr':
        success, ligand_pdbqtstr = get_pdbqtstr_from_pdbstr(ligand_pdbstr, 'prepare_ligand', '-l')
        if not success:
            return put_err(f"Failed to prepare ligand_pdbqtstr: {ligand_pdbqtstr}")
        return ligand_pdbqtstr
    else:
        put_err(f'unknown prepare_backend "{prepare_backend}", '
                f'expected one of {PREPARE_BACKENDS}')


def _prepare_ligand_receptor(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                             macrocycle: bool = False,
                             max_ring_size: int = MACROCYCLE_MAX_RING_SIZE,
                             prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    # the complex must be in the pymol session before any selection runs.
    # `normalize_receptor_chain` loads it only when it has to derive the receptor
    # chain itself, so an explicit `-rc` would leave the session empty and every
    # `cmd.select('complex and ...')` below would silently match nothing.
    cmd.reinitialize()
    cmd.read_pdbstr(complex_pdbstr, 'complex')
    # receptor
    receptor_chain = normalize_receptor_chain(complex_pdbstr, ligand_chain, receptor_chain,
                                             reinitialize=False)
    if not receptor_chain:
        put_err(f'there is no other chain exclude ligand chain {ligand_chain}.')
        return None
    # set receptor
    if cmd.select('receptor', f'complex and (' + ' or '.join([f'chain {c}' for c in receptor_chain]) + ')') == 0:
        put_err(f'can not find receptor with chain {receptor_chain}')
        return None
    receptor_pdbstr = cmd.get_pdbstr('receptor')
    receptor_pdbqtstr = prepare_receptor_pdbqtstr(receptor_pdbstr, prepare_backend)
    if not receptor_pdbqtstr:
        return None
    # ligand
    if cmd.select('ligand', f'complex and chain {ligand_chain}') == 0:
        put_err(f'can not find ligand with chain {ligand_chain}')
        return None
    ligand_pdbstr = cmd.get_pdbstr(f'chain {ligand_chain}')
    ligand_pdbqtstr = prepare_ligand_pdbqtstr(ligand_pdbstr, prepare_backend,
                                              macrocycle, max_ring_size)
    if not ligand_pdbqtstr:
        return None
    return receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def get_grid(lig_pdbstr: str, grid_buffer: float):
    """compute the grid box around the ligand alone, in angstrom.

    This is the rigid-receptor box: it wraps the ligand and nothing else.  The
    flexible-receptor paths need a box that also covers every flexible sidechain and
    use `flex_utils.get_grid_with_flexres` instead.  The two must not be mixed up: a
    ligand-only box would leave the flexible sidechains outside the grid.

    :param lig_pdbstr: ligand pdb string
    :param grid_buffer: box padding in ANGSTROM, added on every side
    :return: (grid_center, grid_size), both in ANGSTROM.  Consumers that want a
        grid point count instead (prepare_gpf.py's `npts`) divide by the spacing
        themselves -- see `prepare_gpf_and_map`, which is the only such consumer.
    """
    cmd.reinitialize()
    cmd.read_pdbstr(lig_pdbstr, 'ligand')
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent('ligand')
    grid_center = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    grid_size = [maxX - minX + grid_buffer*2, maxY - minY + grid_buffer*2, maxZ - minZ + grid_buffer*2]
    return grid_center, grid_size


def prepare_ligand_receptor_grid(tmpdir: str, v, complex_pdbstr: str,
                                  ligand_chain: str, receptor_chain: List[str],
                                  grid_buffer: float, spacing: float = 0.375,
                                  macrocycle: bool = False,
                                  max_ring_size: int = MACROCYCLE_MAX_RING_SIZE,
                                  prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """
    prepare ligand and receptor pdbqt file and grid
    
    both `grid_center` and `grid_size` handed to vina are in ANGSTROM, which is
    what `compute_vina_maps` documents (`spacing` is a separate knob there).
    
    :param grid_buffer: box padding in ANGSTROM
    :param spacing: grid spacing in ANGSTROM
    return: v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr
    """
    pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain,
                                           macrocycle, max_ring_size, prepare_backend)
    if not pdbstr_pack:
        return None
    receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
    # set receptor
    rec_path = os.path.join(tmpdir, 'receptor.pdbqt')
    opts_file(rec_path, 'w', way='str', data=receptor_pdbqtstr)
    v.set_receptor(rec_path)
    # set ligand
    v.set_ligand_from_string(ligand_pdbqtstr)
    # set grid; not bounding box, just box in axis
    grid_center, grid_size = get_grid(ligand_pdbstr, grid_buffer)
    v.compute_vina_maps(center=grid_center, box_size=grid_size, spacing=spacing)
    return v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def run_autogrid4(gpf_path: str, workdir: str, gpf_caller: str = None, glg_name: Optional[str] = None):
    """run autogrid4 on an existing gpf file.

    This is the single place that executes autogrid4, so both paths go through it:

        - `prepare_gpf_and_map`, whose gpf is freshly written by prepare_gpf.py
        - the meeko-based flexible-receptor path, where meeko wrote the gpf itself

    It resolves the ADFRsuite bin dir itself (see `prepend_adfr_bin_to_path`), so
    the caller does not have to arrange $PATH first.

    :param gpf_path: gpf file to run, absolute or relative to `workdir`
    :param workdir: dir to run in, where the maps are written
    :param gpf_caller: prepare_gpf.py path, None means read from config or $PREPARE_GPF
    :param glg_name: log file name, None derives it from `gpf_path`
    """
    if glg_name is None:
        glg_name = Path(os.path.basename(gpf_path)).with_suffix('.glg').name
    # prefer the autogrid4 sitting next to prepare_gpf.py, fall back on $PATH
    autogrid_caller = 'autogrid4'
    gpf_caller = prepend_adfr_bin_to_path(gpf_caller)
    if gpf_caller:
        candidate = os.path.join(os.path.dirname(gpf_caller), 'autogrid4')
        if Path(candidate).exists():
            autogrid_caller = candidate
    os.system(f'cd {workdir} && {autogrid_caller} -p {os.path.basename(gpf_path)} -l {glg_name}')


def prepare_gpf_and_map(docking_dir: str, gpf_caller: str, grid_center: list[float],
                        grid_size: list[float], spacing: float = 0.375):
    """write the gpf for the rigid receptor and run autogrid4.

    `grid_size` is a box size in ANGSTROM (what `get_grid` returns), while
    prepare_gpf.py's `npts` is a count of grid points, so the conversion happens
    here: npts = floor(size / spacing).  This function is the ONLY consumer of the
    grid box that wants point counts; vina and meeko both take angstrom.

    `spacing` defaults to 0.375, which is also prepare_gpf.py's own default, so the
    npts computed here and the spacing the gpf records always agree.  Leave it at the
    default for the AutoDock-GPU paths: this helper does not forward the value to the
    gpf, so a different spacing would size the box by one value and describe it with
    another.  The callers reject a non-default spacing up front
    (`scripts.dock.check_gpu_spacing`).

    NOTE the rounding differs from meeko's gpf writer, which forces an EVEN point count
    on every axis (`meeko/gridbox.py`: `2 * int(size / (2 * spacing))`) while this one
    floors.  A box size that is not an exact multiple of the spacing therefore loses up
    to one point per axis here -- a property of the two writers, not a bug.

    :param docking_dir: dir holding receptor.pdbqt, where the gpf and maps are written
    :param gpf_caller: prepare_gpf.py path; it also locates autogrid4 and the bin dir
        that goes on $PATH (see `prepend_adfr_bin_to_path`)
    :param grid_center: box center in angstrom
    :param grid_size: box size in angstrom
    :param spacing: grid spacing in angstrom, must be 0.375 (prepare_gpf.py's default)
    """
    grid_center_str = ','.join(map(lambda x: f'{x:.3f}', grid_center))
    grid_size_str = ','.join(map(lambda x: f'{int(x / spacing)}', grid_size))
    # pythonsh is called by bare name and lives next to prepare_gpf.py
    gpf_caller = prepend_adfr_bin_to_path(gpf_caller)
    os.system(f'cd {docking_dir} && pythonsh {gpf_caller} -r receptor.pdbqt -o receptor.gpf -p gridcenter={grid_center_str} -p npts={grid_size_str} -p ligand_types="A,C,OA,N,NA,SA,HD"')
    run_autogrid4('receptor.gpf', docking_dir, gpf_caller, 'receptor.glg')

