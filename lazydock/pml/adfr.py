import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Tuple

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file
from pymol import cmd


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


def _prepare_ligand_receptor_grid(tmpdir: str, v, complex_pdbstr: str,
                                  ligand_chain: str, receptor_chain: List[str], grid_buffer: float):
    """
    prepare ligand and receptor pdbqt file and grid
    """
    cmd.reinitialize()
    cmd.read_pdbstr(complex_pdbstr, 'complex')
    if not receptor_chain:
        receptor_chain = list(set(cmd.get_chains()) - set([ligand_chain]))
    if not receptor_chain:
        put_err(f'there is no other chain exclude ligand chain {ligand_chain}.')
        return None
    # set receptor
    if cmd.select('receptor', f'complex and (' + ' or '.join([f'chain {c}' for c in receptor_chain]) + ')') == 0:
        put_err(f'can not find receptor with chain {receptor_chain}')
        return None
    rec_path = os.path.join(tmpdir, 'receptor.pdbqt')
    success, receptor_pdbqtstr = get_pdbqtstr_from_pdbstr(cmd.get_pdbstr('receptor'), 'prepare_receptor', '-r')
    if not success:
        put_err(f"Failed to prepare receptor_pdbqtstr: {receptor_pdbqtstr}")
        return None
    opts_file(rec_path, 'w', way='str', data=receptor_pdbqtstr)
    v.set_receptor(rec_path)
    # set ligand
    if cmd.select('ligand', f'complex and chain {ligand_chain}') == 0:
        put_err(f'can not find ligand with chain {ligand_chain}')
        return None
    success, ligand_pdbqtstr = get_pdbqtstr_from_pdbstr(cmd.get_pdbstr(f'chain {ligand_chain}'), 'prepare_ligand', '-l')
    if not success:
        put_err(f"Failed to prepare ligand_pdbqtstr: {ligand_pdbqtstr}")
        return None
    v.set_ligand_from_string(ligand_pdbqtstr)
    # set grid; not bounding box, just box in axis
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent('ligand')
    grid_center = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    grid_size = [maxX - minX + grid_buffer*2, maxY - minY + grid_buffer*2, maxZ - minZ + grid_buffer*2]
    v.compute_vina_maps(center=grid_center, box_size=grid_size, spacing=0.375)
    return v


def calc_vina_score(score_name: str, grid_buffer: float,
                    complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str], cpu: int = 1):
    """"total", "lig_inter", "flex_inter", "other_inter", "flex_intra", "lig_intra", "torsions", "-lig_intra"
    calculate vina score of a complex pdb file.
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box extend, box will calculate with pymol.cmd.extend without align to axis.
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param cpu: cpu number
    
    :return: score list, if None, means failed
    """
    from vina import Vina
    v = Vina(sf_name=score_name, cpu=cpu, verbosity=False)    
    with tempfile.TemporaryDirectory() as tmpdir:
        v = _prepare_ligand_receptor_grid(tmpdir, v, complex_pdbstr, ligand_chain, receptor_chain, grid_buffer)
        if not v:
            return None
        # compute and retrun score
        return v.score().tolist()


def perform_dock(score_name: str, grid_buffer: float,
                complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                output_path: str,
                exhaustiveness: int, n_poses: int, min_rmsd: float = 1, max_evals: int = 0):
    """
    perform docking with vina
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box extend, box will calculate with pymol.cmd.extend without align to axis.
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param output_path: output path
    :param exhaustiveness: exhaustiveness, equal to cpu number
    :param n_poses: n_poses
    :param min_rmsd: min_rmsd
    :param max_evals: max_evals
    
    :return: True if success, False otherwise
    """
    from vina import Vina
    v = Vina(sf_name=score_name, cpu=exhaustiveness, verbosity=False)    
    with tempfile.TemporaryDirectory() as tmpdir:
        v = _prepare_ligand_receptor_grid(tmpdir, v, complex_pdbstr, ligand_chain, receptor_chain, grid_buffer)
        if not v:
            return False
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=min_rmsd, max_evals=max_evals)
        v.write_poses(output_path, n_poses, overwrite=True)
        return True

if __name__ == '__main__':
    pdbqtstr = get_pdbqtstr_from_pdbstr("")
