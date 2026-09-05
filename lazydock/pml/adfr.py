import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file
from pymol import cmd

from lazydock.config import GlobalConfig
from lazydock.pml.utils import new_pml_context


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


def _prepare_ligand_receptor(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str]):
    cmd.reinitialize()
    cmd.read_pdbstr(complex_pdbstr, 'complex')
    # receptor
    if not receptor_chain:
        receptor_chain = list(set(cmd.get_chains()) - set([ligand_chain]))
    if not receptor_chain:
        put_err(f'there is no other chain exclude ligand chain {ligand_chain}.')
        return None
    # set receptor
    if cmd.select('receptor', f'complex and (' + ' or '.join([f'chain {c}' for c in receptor_chain]) + ')') == 0:
        put_err(f'can not find receptor with chain {receptor_chain}')
        return None
    receptor_pdbstr = cmd.get_pdbstr('receptor')
    success, receptor_pdbqtstr = get_pdbqtstr_from_pdbstr(receptor_pdbstr, 'prepare_receptor', '-r')
    if not success:
        put_err(f"Failed to prepare receptor_pdbqtstr: {receptor_pdbqtstr}")
        return None
    # ligand
    if cmd.select('ligand', f'complex and chain {ligand_chain}') == 0:
        put_err(f'can not find ligand with chain {ligand_chain}')
        return None
    ligand_pdbstr = cmd.get_pdbstr(f'chain {ligand_chain}')
    success, ligand_pdbqtstr = get_pdbqtstr_from_pdbstr(ligand_pdbstr, 'prepare_ligand', '-l')
    if not success:
        put_err(f"Failed to prepare ligand_pdbqtstr: {ligand_pdbqtstr}")
        return None
    return receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def _get_grid(lig_pdbstr: str, grid_buffer: float):
    cmd.reinitialize()
    cmd.read_pdbstr(lig_pdbstr, 'ligand')
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent('ligand')
    grid_center = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    grid_size = [maxX - minX + grid_buffer*2, maxY - minY + grid_buffer*2, maxZ - minZ + grid_buffer*2]
    return grid_center, grid_size

def _prepare_ligand_receptor_grid(tmpdir: str, v, complex_pdbstr: str,
                                  ligand_chain: str, receptor_chain: List[str], grid_buffer: float):
    """
    prepare ligand and receptor pdbqt file and grid
    
    return: v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr
    """
    pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain)
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
    grid_center, grid_size = _get_grid(ligand_pdbstr, grid_buffer)
    v.compute_vina_maps(center=grid_center, box_size=grid_size, spacing=0.375)
    return v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


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
        return v[0].score().tolist()


def perform_vina_dock(score_name: str, grid_buffer: float,
                      complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                      exhaustiveness: int, n_poses: int, min_rmsd: float = 1, max_evals: int = 0):
    """
    perform docking with vina
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box extend, box will calculate with pymol.cmd.extend without align to axis.
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
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
        v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = v
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=min_rmsd, max_evals=max_evals)
        output_path = os.path.join(tmpdir, 'poses.pdbqt')
        v.write_poses(output_path, n_poses, overwrite=True)
        result_str = opts_file(output_path, 'r', way='str')
        os.remove(output_path)
        return result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def prepare_gpf_and_map(docking_dir: str, caller: str, grid_center: list[float], grid_size: list[int]):
    grid_center = ','.join(map(lambda x: f'{x:.3f}', grid_center))
    grid_size = ','.join(map(lambda x: f'{x:.0f}', grid_size))
    os.system(f'cd {docking_dir} && pythonsh {caller} -r receptor.pdbqt -o receptor.gpf -p gridcenter={grid_center} -p npts={grid_size} -p ligand_types="A,C,OA,N,NA,SA,HD"')
    os.system(f'cd {docking_dir} && autogrid4 -p receptor.gpf -l receptor.glg')


def perform_autodock_gpu_dock_core(docking_dir: str, caller: str, ffile: str, lfile: str,
                                   seed: int = 0, nrun: int = 5, nev: int = 500000,
                                   stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    """
    peform docking with AutoDock-GPU
    :param docking_dir: docking dir, must be unique for parallel docking
    :param caller: caller, such as autodock_gpu_64wi_R128
    :param ffile: ffile, such as receptor.maps.fld
    :param lfile: lfile, such as ligand.pdbqt
    :param seed: seed
    :param nrun: nrun
    :param nev: nev
    :param stopstd: stopstd
    :param heurmax: heurmax
    :param gpu_id: gpu id
    
    :return: True if success, False otherwise
    """
    cmd = [
        'cd', docking_dir, '&&',
        caller,
        '--ffile', ffile,
        '--lfile', lfile,
        '--seed', str(seed),
        '--resnam', f'dock',
        '--nrun', str(nrun),
        '--nev', str(nev),
        '--stopstd', str(stopstd),
        '--heurmax', str(heurmax),
        '--devnum', str(gpu_id+1),
    ]
    proc = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = proc.communicate()
    output_str = stdout + stderr
    # return path
    if proc.returncode != 0 or not (docking_dir / Path(f'dock.dlg')).exists():
        print(f'Error: docking failed, output_str: {output_str}')
    return os.path.join(docking_dir, 'dock.dlg')


def perform_autodock_gpu_dock(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],grid_buffer: float,
                              autodock_caller: Optional[str] = None, gpf_caller: Optional[str] = None,
                              seed: int = 0, nrun: int = 5, nev: int = 500000,
                              stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        opts_file(os.path.join(tmpdir, 'receptor.pdbqt'), 'w', way='str', data=receptor_pdbqtstr)
        opts_file(os.path.join(tmpdir, 'ligand.pdbqt'), 'w', way='str', data=ligand_pdbqtstr)        
        # get grid center and size
        grid_center, grid_size = _get_grid(ligand_pdbstr, grid_buffer)
        # setup docking files
        if not gpf_caller:
            gpf_caller = GlobalConfig.named_paths["prepare_gpf"] or os.environ.get("PREPARE_GPF")
        if not gpf_caller:
            return put_err('gpf_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $PREPARE_GPF')
        prepare_gpf_and_map(tmpdir, gpf_caller, grid_center, grid_size)
        # perform docking
        if not autodock_caller:
            autodock_caller = GlobalConfig.named_paths["autodock_gpu"] or os.environ.get("AUTODOCK_GPU")
        if not autodock_caller:
            return put_err('autodock_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $AUTODOCK_GPU')
        docking_path = perform_autodock_gpu_dock_core(tmpdir, autodock_caller, 'receptor.maps.fld', 'ligand.pdbqt',
                                                      seed, nrun, nev, stopstd, heurmax, gpu_id)
        if not docking_path:
            return None
        # return docking path
        result_str = opts_file(docking_path, 'r', way='str')
        return result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


if __name__ == '__main__':
    pdbqtstr = get_pdbqtstr_from_pdbstr("")
    complex_pdbstr = opts_file('data_tmp/pdb/complex_cycpep.pdb')
    result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr \
        = perform_autodock_gpu_dock(complex_pdbstr, 'P', ['A'], 4)
