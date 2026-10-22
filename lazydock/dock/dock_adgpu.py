import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file

from lazydock.config import GlobalConfig
from lazydock.dock.adfr import (_prepare_ligand_receptor, get_grid,
                                prepare_gpf_and_map, prepend_adfr_bin_to_path,
                                run_autogrid4)
from lazydock.dock.config import PREPARE_BACKEND_DEFAULT
from lazydock.dock.flex_utils import FlexResidue
from lazydock.dock.meeko_utils import prepare_flex_docking_setup


def perform_autodock_gpu_dock_core(docking_dir: str, caller: str, ffile: str, lfile: str,
                                   flexfile: str = None,
                                   seed: int = 0, nrun: int = 5, nev: int = 500000,
                                   stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    """
    perform docking with AutoDock-GPU, optionally with flexible receptor sidechains

    The rigid and the flexible-receptor runs differ by exactly one command-line flag
    (`--flexres`), so they share this one implementation: pass `flexfile` to enable
    flexible docking, leave it None for a rigid receptor.

    :param docking_dir: docking dir, must be unique for parallel docking
    :param caller: caller, such as autodock_gpu_64wi_R128
    :param ffile: ffile, such as receptor.maps.fld or receptor_rigid.maps.fld
    :param lfile: lfile, such as ligand.pdbqt
    :param flexfile: flexible residue pdbqt file, such as receptor_flex.pdbqt.
        None (the default) runs a rigid-receptor dock
    :param seed, nrun, nev, stopstd, heurmax, gpu_id: AutoDock-GPU options

    :return: dlg file path, None if failed
    """
    cmd_list = [
        'cd', docking_dir, '&&',
        caller,
        '--ffile', ffile,
        '--lfile', lfile,
    ]
    if flexfile:
        cmd_list += ['--flexres', flexfile]
    cmd_list += [
        '--seed', str(seed),
        '--resnam', 'dock',
        '--nrun', str(nrun),
        '--nev', str(nev),
        '--stopstd', str(stopstd),
        '--heurmax', str(heurmax),
        '--devnum', str(gpu_id + 1),
    ]
    proc = subprocess.Popen(' '.join(cmd_list), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = proc.communicate()
    output_str = stdout + stderr
    dlg_path = os.path.join(docking_dir, 'dock.dlg')
    if proc.returncode != 0 or not Path(dlg_path).exists():
        put_err(f'docking failed, output_str: {output_str}')
        return None
    return dlg_path


def perform_autodock_gpu_dock(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                              grid_buffer: float,
                              autodock_caller: Optional[str] = None, gpf_caller: Optional[str] = None,
                              seed: int = 0, nrun: int = 5, nev: int = 500000,
                              stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0,
                              prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """
    perform docking with AutoDock-GPU

    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param grid_buffer: grid padding in ANGSTROM, added on every side of the box.
        Spacing is not a parameter here: this path's maps come from `prepare_gpf.py`
        + `autogrid4`, and the spacing that AutoDock-GPU then reads back from the
        fld is prepare_gpf.py's own 0.375 default (`meeko/gridbox.py` fixes the same
        value on the flex path), so there is no way to honour another value.  The
        CLI rejects a non-default `--spacing` up front
        (`scripts.dock.check_gpu_spacing`).
    :param autodock_caller: autodock_gpu executable path, if None, will be read from config or $AUTODOCK_GPU
    :param gpf_caller: prepare_gpf.py path, if None, will be read from config or $PREPARE_GPF
    :param seed, nrun, nev, stopstd, heurmax, gpu_id: AutoDock-GPU options
    :param prepare_backend: 'meeko' or 'adfr', which toolkit prepares the receptor and
        ligand pdbqt.  Meeko was verified to work end-to-end on this path (its pdbqt
        feeds prepare_gpf.py and autogrid4 directly, 9 maps written), but the two
        backends type some nitrogens/polar hydrogens differently, so their scores are
        NOT comparable (see `PREPARE_BACKENDS`).
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain,
                                               prepare_backend=prepare_backend)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        opts_file(os.path.join(tmpdir, 'receptor.pdbqt'), 'w', way='str', data=receptor_pdbqtstr)
        opts_file(os.path.join(tmpdir, 'ligand.pdbqt'), 'w', way='str', data=ligand_pdbqtstr)        
        # get grid center and size, both in ANGSTROM
        grid_center, grid_size = get_grid(ligand_pdbstr, grid_buffer)
        # setup docking files
        gpf_caller = prepend_adfr_bin_to_path(gpf_caller)
        if not gpf_caller:
            return 
        # this consumer wants grid point counts, so it converts with the spacing;
        # the default is used on purpose, see the docstring
        prepare_gpf_and_map(tmpdir, gpf_caller, grid_center, grid_size)
        # perform docking
        if not autodock_caller:
            autodock_caller = GlobalConfig.named_paths["autodock_gpu"] or os.environ.get("AUTODOCK_GPU")
        if not autodock_caller:
            return put_err('autodock_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $AUTODOCK_GPU')
        # `flexfile=None` is passed by keyword on purpose: the core also takes
        # `flexfile` as its 5th parameter, so leaving it out shifts every option
        # below it by one -- measured, the old positional call put the seed into
        # `flexfile` and then died joining the argv (TypeError: sequence item 9:
        # expected str instance, int found)
        docking_path = perform_autodock_gpu_dock_core(
            tmpdir, autodock_caller, 'receptor.maps.fld', 'ligand.pdbqt',
            flexfile=None, seed=seed, nrun=nrun, nev=nev,
            stopstd=stopstd, heurmax=heurmax, gpu_id=gpu_id)
        if not docking_path:
            return None
        # return docking path
        result_str = opts_file(docking_path, 'r', way='str')
        return result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def perform_autodock_gpu_flex_dock(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str], grid_buffer: float,
                                   flex_residues: Union[str, List[str], List[FlexResidue]] = None,
                                   flex_distance: float = None, flex_interaction: bool = False,
                                   flex_interaction_mode: Union[str, List[str]] = 'all', flex_interaction_cutoff: float = 4.,
                                   autodock_caller: Optional[str] = None, gpf_caller: Optional[str] = None,
                                   seed: int = 0, nrun: int = 5, nev: int = 500000,
                                   stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0,
                                   prepare_backend: str = PREPARE_BACKEND_DEFAULT):
    """
    perform docking with AutoDock-GPU with flexible receptor sidechains.
    
    flexible residues could be provided in three ways, and they could be combined:
        - flex_residues: explicit residue list, such as 'A:LYS68,A:LYS50' or [(chain, resname, resnum), ...]
        - flex_distance: select receptor residues within this distance of any ligand heavy atom
        - flex_interaction: select receptor residues which have interaction with ligand by fplip
    
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param grid_buffer: grid padding in ANGSTROM, added on every side of the box
        that covers the ligand and every flexible residue.  Spacing is not exposed on
        this path: meeko writes the rigid gpf with its own fixed 0.375 A spacing, and
        AutoDock-GPU reads that spacing from the gpf, so the two always agree.
    :param flex_residues: explicit flexible residues
    :param flex_distance: distance cutoff, select residues within this distance of ligand heavy atoms
    :param flex_interaction: FLAG, whether to select flexible residues by fplip interaction
    :param flex_interaction_mode: interaction types for fplip, 'all' to include all modes
    :param flex_interaction_cutoff: interaction distance cutoff for fplip
    :param autodock_caller: autodock_gpu executable path, if None, will be read from config or $AUTODOCK_GPU
    :param gpf_caller: prepare_gpf.py path, if None, will be read from config or
        $PREPARE_GPF.  It locates autogrid4, which is run on the gpf meeko wrote
    :param prepare_backend: 'meeko' or 'adfr', which toolkit prepares the LIGAND pdbqt.
        The receptor is always split by meeko here (ADT can not produce a rigid part
        plus flexible sidechains at all), so this switch only ever changes the ligand.
        Use 'adfr' only when the ligand must match an ADT-prepared run: it makes
        ADFRsuite's `prepare_ligand` a hard requirement, and the two backends are not
        score-comparable (see `PREPARE_BACKENDS`)
    
    :return: (result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr),
        None if failed. The 6th item is the ligand chain pdb string.
        Use `build_flex_complex_pdbstr` to build the complex pdb from the first 3 items
        and one pose pdb string from the dlg.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor pdb string, and ligand pdbqt
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain,
                                               prepare_backend=prepare_backend)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, _, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        opts_file(os.path.join(tmpdir, 'ligand.pdbqt'), 'w', way='str', data=ligand_pdbqtstr)
        # flexible residues, ADFRsuite bin dir, grid box and the rigid/flex receptor
        # split; shared with the vina flex path (see `prepare_flex_docking_setup`)
        setup = prepare_flex_docking_setup(
            tmpdir, complex_pdbstr, ligand_chain, receptor_chain,
            receptor_pdbstr, ligand_pdbstr, grid_buffer,
            flex_residues, flex_distance, flex_interaction,
            flex_interaction_mode, flex_interaction_cutoff)
        if not setup:
            return None
        rigid_pdbqt_path, flex_pdbqt_path = setup.rigid_pdbqt_path, setup.flex_pdbqt_path
        # maps are computed only from the rigid part.  meeko already wrote the gpf
        # for that rigid part (see `prepare_receptor_with_flexres`), so this path
        # must NOT call prepare_gpf_and_map: that would rebuild the gpf from
        # prepare_gpf.py and silently drop the flexres box.
        run_autogrid4(setup.rigid_gpf_path, tmpdir, gpf_caller, 'receptor_rigid.glg')
        # perform docking with flexible residues
        if not autodock_caller:
            autodock_caller = GlobalConfig.named_paths["autodock_gpu"] or os.environ.get("AUTODOCK_GPU")
        if not autodock_caller:
            return put_err('autodock_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $AUTODOCK_GPU')
        # options by keyword for the same reason as in `perform_autodock_gpu_dock`
        docking_path = perform_autodock_gpu_dock_core(
            tmpdir, autodock_caller, 'receptor_rigid.maps.fld', 'ligand.pdbqt',
            flexfile=os.path.basename(flex_pdbqt_path), seed=seed, nrun=nrun, nev=nev,
            stopstd=stopstd, heurmax=heurmax, gpu_id=gpu_id)
        if not docking_path:
            return None
        result_str = opts_file(docking_path, 'r', way='str')
        rigid_pdbqtstr = opts_file(rigid_pdbqt_path, 'r', way='str')
        flex_pdbqtstr = opts_file(flex_pdbqt_path, 'r', way='str')
        return result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr