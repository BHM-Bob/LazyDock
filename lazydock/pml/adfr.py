import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

from mbapy_lite.file import opts_file


def get_pdbqtstr_from_pdbstr(pdbstr: str, caller: str = 'prepare_ligand', entry_param: str = '-l',
                             options: List[Tuple['str', str]] = None):
    """
    call ADFR programe to prepare ligand pdbqt file or receptor pdbqt file
    
    :param pdbstr: pdb string
    :param caller: ADFR programe name, could be 'prepare_ligand' or 'prepare_receptor', or absolute path
    :param options: options for ADFR programe, if None, will be [('-A', 'hydrogens')]
    :return: FLAG, pdbqt string, FLAG is True if success
    """
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


def prepare_gpf(caller: str, grid_center: list[float], grid_size: list[int]):
    grid_center = ','.join(map(lambda x: f'{x:.3f}', grid_center))
    grid_size = ','.join(map(lambda x: f'{x:.0f}', grid_size))
    os.system(f'pythonsh {caller} -r receptor.pdbqt -o receptor.gpf -p gridcenter={grid_center} -p npts={grid_size} -p ligand_types="A,C,OA,N,NA,SA,HD"')


if __name__ == '__main__':
    pdbqtstr = get_pdbqtstr_from_pdbstr("")
