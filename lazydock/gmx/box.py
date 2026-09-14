'''
Date: 2026-09-08
Description: box / solvate / ion 相关自由函数。
             Extracted from lazydock.scripts.run_gmx.py (simple_protein.make_box),
             parameterized (no self.args dependency) for reuse by ABFE prepare.
             run_gmx.py is intentionally NOT refactored to call these yet (backward compat).
'''
import os
import re
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import opts_file
from pymol import cmd

from lazydock.gmx.run import Gromacs
from lazydock.utils import uuid4


def get_box(mol_path: Path, padding: Union[float, List[float]], sele: str = 'not resn SOL') -> Tuple[List[float], List[float]]:
    """Return (box_center, box_size) of a molecule using PyMOL get_extent.

    box_center in nm; box_size in nm (extent + padding on both sides).
    """
    cmd.reinitialize()
    name = uuid4()
    cmd.load(str(mol_path), mol_path.name)
    cmd.select(name, sele)
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(name)
    box_center = list(map(lambda x: x / 20, [maxX + minX, maxY + minY, maxZ + minZ]))
    if isinstance(padding, float):
        padding = [padding] * 6
    padding1, padding2 = padding[:3], padding[3:]
    box_size = list(map(lambda x: x[0] / 10 + x[1] + x[2],
                        zip([maxX - minX, maxY - minY, maxZ - minZ], padding1, padding2)))
    cmd.reinitialize()
    return box_center, box_size


def apply_box_center_shift(box_center: List[float], shift: List[float]) -> List[float]:
    return [x + s for x, s in zip(box_center, shift)]


def measure_span(mol_path: Path, sele: str = None) -> float:
    """Return the maximum x/y/z extent (nm) of a molecule via PyMOL get_extent.

    Used by ABFE to derive the box padding that guarantees
    half-shortest-box-vector >= rlist (see phase3_box_auto.md).
    """
    cmd.reinitialize()
    name = uuid4()
    cmd.load(str(mol_path), mol_path.name)
    if sele:
        cmd.select(name, sele)
        ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(name)
    else:
        ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent(mol_path.name)
    cmd.reinitialize()
    return max((maxX - minX), (maxY - minY), (maxZ - minZ)) / 10  # Angstrom -> nm


def make_box(protein_path: Path, main_name: str, gmx: Gromacs,
             editconf_args: str = "-c -d 1.2 -bt dodecahedron",
             solvate_args: str = "-cs spc216.gro",
             auto_box: bool = False,
             auto_box_padding: Union[float, List[float]] = 1.2,
             auto_box_shift: List[float] = None,
             mono_lock=None) -> None:
    """Build the simulation box (editconf) and solvate.

    Mirrors run_gmx.simple_protein.make_box; parameterized for reuse.

    Parameters
    ----------
    protein_path : Path
        path to the solute gro (must have a topol.top in the same dir).
    main_name : str
        basename of the output files (e.g. 'protein' -> protein_newbox.gro ...).
    gmx : Gromacs
    editconf_args : str
        args for editconf; when auto_box, the '-d ... -bt dodecahedron' part is
        replaced by an explicit -box generated from PyMOL get_extent.
    solvate_args : str
        args for solvate.
    auto_box : bool
        use PyMOL contour box (anisotropic) instead of -d.
    auto_box_padding : float | [6]
        padding for the contour box (nm).
    auto_box_shift : [x,y,z]
        extra shift applied to the box center.
    mono_lock : threading.Lock or None
        serialize PyMOL access (pymol not thread-safe); pass a lock when used
        from multiple threads.
    """
    auto_box_shift = auto_box_shift or [0, 0, 0]
    if auto_box:
        # get shift from first editconf
        if mono_lock is not None:
            with mono_lock:
                _, box_size = get_box(protein_path, auto_box_padding)
        else:
            _, box_size = get_box(protein_path, auto_box_padding)
        manual_box_cmd = f'-box {" ".join(map(lambda x: f"{x:.2f}", box_size))}'
        editconf_args = editconf_args.replace('-d 1.2 -bt dodecahedron', ' ') + manual_box_cmd
        _, log_path = gmx.run_gmx_with_expect(f'editconf {editconf_args}', f=f'{main_name}.gro',
                                              o=f'{main_name}_newbox_tmp.gro', enable_log=True)
        shift_line = list(filter(lambda x: 'new center' in x.strip(), opts_file(log_path, way='lines')))[0]
        shift = list(map(float, re.findall(r'[\d\-\.]+', shift_line)))
        # get solvated box from first solvate
        shutil.copy(protein_path.parent / 'topol.top', protein_path.parent / 'topol_tmp.top')
        gmx.run_gmx_with_expect(f'solvate {solvate_args}', cp=f'{main_name}_newbox_tmp.gro',
                                o=f'{main_name}_solv_tmp.gro', p='topol_tmp.top')
        if mono_lock is not None:
            with mono_lock:
                solv_center, solv_size = get_box(protein_path.parent / f'{main_name}_solv_tmp.gro',
                                                 auto_box_padding, 'resn SOL')
                prot_center, _ = get_box(protein_path.parent / f'{main_name}_newbox_tmp.gro', auto_box_padding)
        else:
            solv_center, solv_size = get_box(protein_path.parent / f'{main_name}_solv_tmp.gro',
                                             auto_box_padding, 'resn SOL')
            prot_center, _ = get_box(protein_path.parent / f'{main_name}_newbox_tmp.gro', auto_box_padding)
        put_log(f'protein box size: {box_size}, tmp solvated box size: {solv_size}, '
                f'protein center: {prot_center}, tmp solvated center: {solv_center}, shift: {shift}')
        # calculate new box center
        box_center = [s + (x1 - x2) + s2 for s, x1, x2, s2 in zip(shift, solv_center, prot_center, auto_box_shift)]
        # run editconf with new box center and size
        editconf_args += f' -center {" ".join(map(lambda x: f"{x:.2f}", box_center))}'
        gmx.run_gmx_with_expect(f'editconf {editconf_args}', f=f'{main_name}.gro', o=f'{main_name}_newbox.gro')
    else:
        gmx.run_gmx_with_expect(f'editconf {editconf_args}', f=f'{main_name}.gro', o=f'{main_name}_newbox.gro')
    # solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top
    gmx.run_gmx_with_expect(f'solvate {solvate_args}', cp=f'{main_name}_newbox.gro',
                            o=f'{main_name}_solv.gro', p='topol.top')


def add_ions(protein_path: Path, main_name: str, gmx: Gromacs,
             ion_mdp: str, genion_args: str = "-pname NA -nname CL -neutral",
             genion_groups: str = "SOL", maxwarn: int = 0) -> None:
    """Grompp an ions.tpr and run genion (replace water by ions).

    Mirrors the 'grompp ions + genion' part of run_gmx.simple_protein.make_box.
    `ion_mdp` is an mdp file path used to build ions.tpr.
    """
    # grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr
    gmx.run_gmx_with_expect('grompp', f=ion_mdp, c=f'{main_name}_solv.gro', p='topol.top',
                            o='ions.tpr', maxwarn=maxwarn)
    g_groups = genion_groups
    if genion_groups == 'SOL':
        groups = gmx.get_groups('ions.tpr')
        if 'SOL' not in groups:
            put_err(f'can not find SOL group in ions.tpr, skip.')
        else:
            g_groups = groups['SOL']
    gmx.run_gmx_with_expect(f'genion {genion_args}', s='ions.tpr', o=f'{main_name}_solv_ions.gro',
                            p='topol.top', expect_actions=[{'Select a group:': f'{g_groups}\r',
                                                            'No ions to add': '', '\\timeout': ''}],
                            expect_settings={'start_timeout': 600})


#############################################################################
# 底层单步函数 (Phase 3B "选项 B" 统一: Solvate / run_gmx 共享调用)
#############################################################################

def run_editconf(gmx: Gromacs, f: str, o: str = None, bt: str = 'dodecahedron',
                 d: float = None, box: list = None, angles: list = None,
                 c: bool = False) -> None:
    """Thin wrapper over gmx editconf covering all box options."""
    kwargs: dict = dict(f=f, bt=bt)
    if o:
        kwargs['o'] = o
    if d is not None:
        kwargs['d'] = d
    if box:
        kwargs['box'] = ' '.join(map(str, box))
    if angles:
        kwargs['angles'] = ' '.join(map(str, angles))
    if c:
        kwargs['c'] = True
    gmx.run_gmx_with_expect('editconf', **kwargs)  # type: ignore[arg-type]


def run_solvate(gmx: Gromacs, cp: str, p: str, cs: str, o: str = None) -> None:
    """Thin wrapper over gmx solvate (custom water model gro supported)."""
    kwargs = dict(cp=cp, p=p, cs=cs)
    if o:
        kwargs['o'] = o
    gmx.run_gmx_with_expect('solvate', **kwargs)


def run_grompp_ions(gmx: Gromacs, f: str, c: str, p: str, o: str = 'ions.tpr',
                    maxwarn: int = 0) -> None:
    """grompp the ions.tpr used before genion."""
    gmx.run_gmx_with_expect('grompp', f=f, c=c, p=p, o=o, maxwarn=maxwarn)


def run_genion(gmx: Gromacs, s: str, p: str, o: str, pname: str = 'NA',
               nname: str = 'CL', conc: float = 150e-3, neutral: bool = True,
               rmin: float = 1.0, groups: str = 'SOL', start_timeout: int = 600) -> None:
    """genion with SOL-group selection via expect; supports conc/neutral/rmin."""
    g_groups = groups
    if groups == 'SOL':
        grps = gmx.get_groups(s)
        if 'SOL' not in grps:
            put_err(f'can not find SOL group in {s}, skip.')
        else:
            g_groups = grps['SOL']
    gmx.run_gmx_with_expect('genion', s=s, p=p, o=o, pname=pname, nname=nname,
                            conc=conc, neutral=neutral, rmin=rmin,
                            expect_actions=[{'Select a group:': f'{g_groups}\r',
                                             'No ions to add': '', '\\timeout': ''}],
                            expect_settings={'start_timeout': start_timeout})


if __name__ == '__main__':
    pass