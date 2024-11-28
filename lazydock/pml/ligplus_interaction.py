import os
import platform
import re
import tempfile
import time
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from mbapy_lite.file import opts_file
from mbapy_lite.web import TaskPool
from pymol import cmd
from tqdm import tqdm

if __name__ == '__main__':
    from lazydock import config
    from lazydock.pml.interaction_utils import sort_func
    from lazydock.utils import uuid4
else:
    from .. import config
    from ..utils import uuid4
    from .interaction_utils import sort_func
    
    
def parse_ligplus_output(result_path: str) -> Dict[str, List[Tuple[Tuple[int, str, str], Tuple[int, str, str], float]]]:
    """
    parse single LigPlus output files (ligplus.hhb or ligplus.nnb) to get atom-level interactions between receptor and ligand.
    
    file format:
    RRR.C.NNNNI.AAAA.....RRR.C.NNNNI.AAAA....DDDD
    
    Parameters:
        - result_path (str): LigPlus output file path, either ligplus.hhb or ligplus.nnb
    Returns:
        interactions (list): ((receptor_resi, receptor_resn, receptor_chain), (ligand_resi, ligand_resn, ligand_chain), distance)
    """
    lines = re.findall(r'([A-Z]+) ([A-Z]) + ([0-9]+) +(\w+) +([A-Z]+) ([A-Z]) + ([0-9]+) +(\w+) +([0-9\.]+)\n', opts_file(result_path))
    return [((int(line[2]), line[0], line[1]), (int(line[6]), line[4], line[5]), float(line[-1])) for line in lines]

    
def _run_ligplus_for_complex(ligplus_dir: str, complex_pdbstr: str,
                             receptor_chain: str, ligand_chain: str, mode: List[str], cutoff: float,
                             force_cwd: bool = False, **kwargs) -> Dict[str, List[Tuple[Tuple[int, str, str], Tuple[int, str, str], float]]]:
    """
    run LigPlus for a complex, assume comonents.cif is in the LigPlus params directory
    command lines refer to https://github.com/eachanjohnson/dimpyplot/blob/master/dimpyplot.py
    """
    # set ligplus lib path
    pwd_cmd = 'pwd'
    if platform.platform().startswith('Windows'):
        ligplus_lib = os.path.join(ligplus_dir, 'lib', 'exe_win')
        pwd_cmd = 'cd'
    elif platform.platform().startswith('Linux'):
        ligplus_lib = os.path.join(ligplus_dir, 'lib', 'exe_linux')
    else:
        raise ValueError(f'Unsupported platform: {platform.platform()}')
    ligplus_params = os.path.join(ligplus_dir, 'lib', 'params')
    # run LigPlus in temporary directory
    w_dir_handle = tempfile.TemporaryDirectory()
    w_dir = w_dir_handle.name + '/' # because hbplus just append file name to w_dir
    if force_cwd:
        now_wd = os.getcwd()
        os.chdir(w_dir)
        print(f'force change working directory to {w_dir}')
    # save complex_pdbstr into w_dir
    complex_path = os.path.join(w_dir, 'complex.pdb')
    opts_file(complex_path, 'w', data=complex_pdbstr)
    # Run HBadd
    os.system(' '.join([f'cd "{w_dir}" && {pwd_cmd} &&', os.path.join(ligplus_lib, "hbadd"), complex_path, os.path.join(ligplus_params, 'components.cif'), '-wkdir', w_dir]))
    os.system(' '.join([f'cd "{w_dir}" && {pwd_cmd} &&', os.path.join(ligplus_lib, "hbplus"), '-L', '-h', '2.90', '-d', '3.90', '-N', complex_path, '-wkdir', w_dir]))
    os.system(' '.join([f'cd "{w_dir}" && {pwd_cmd} &&', os.path.join(ligplus_lib, "hbplus"), '-L', '-h', '2.70', '-d', '3.35', complex_path, '-wkdir', w_dir]))
    os.system(' '.join([f'cd "{w_dir}" && {pwd_cmd} &&', os.path.join(ligplus_lib, "dimer"), complex_path, receptor_chain, ligand_chain]))
    os.system(' '.join([f'cd "{w_dir}" && {pwd_cmd} &&', os.path.join(ligplus_lib, "dimhtml"), 'none', '-dimp', '-dir', w_dir, '-flip', '-ctype', '1']))
    os.system(' '.join([f'cd "{w_dir}" && {pwd_cmd} &&',
        os.path.join(ligplus_lib, "ligplot"), os.path.join(w_dir, 'dimplot.pdb'), '-wkdir', w_dir,
        '-prm', os.path.join(ligplus_params, 'dimplot.prm'), '-ctype', '1'
        ]))
    # parse output
    try:
        hhb_lines = parse_ligplus_output(os.path.join(w_dir, 'ligplot.hhb'))
        nnb_lines = parse_ligplus_output(os.path.join(w_dir, 'ligplot.nnb'))
    except FileNotFoundError:
        if kwargs.get('is_retry', False):
            print(f'LigPlus output file not found in {w_dir}, retry once')
            return _run_ligplus_for_complex(ligplus_dir, complex_pdbstr, receptor_chain, ligand_chain, mode, cutoff, is_retry=True)
        else:
            print(f'LigPlus output file not found in {w_dir}, set result to empty')
        hhb_lines, nnb_lines = [], []
    # try to close temporary directory
    try:
        w_dir_handle.cleanup()
    except Exception as e:
        print(f'Error when cleaning up temporary directory: {e}, skip')
    if force_cwd:
        os.chdir(now_wd)
        print(f'change working directory back to {now_wd}')
    hhb_lines = [(line[1], line[0], line[2]) if line[0][2] == ligand_chain else line for line in hhb_lines if line[2] <= cutoff]
    nnb_lines = [(line[1], line[0], line[2]) if line[0][2] == ligand_chain else line for line in nnb_lines if line[2] <= cutoff]
    interactions = {'Hydrogen Bonds': hhb_lines, 'Non-bonded Interactions': nnb_lines}
    return {k:v for k,v in interactions.items() if k in mode}

    
def run_ligplus(ligplus_dir: str, receptor: str = None, ligand: str = None,
                complex: str = None, receptor_chain: str = None, ligand_chain: str = None,
                mode: List[str] = None, cutoff: float = 4., taskpool: TaskPool = None, force_cwd: bool = False) -> Dict[str, List[Tuple[Tuple[int, str, str], Tuple[int, str, str], float]]]:
    """
    run LigPlus to get atom-level interactions between receptor and ligand.
    Parameters:
        - ligplus_dir (str): LigPlus directory, which contains LigPlus.jar
        - receptor (str): receptor name in pymol, or pdb file path, default None
        - ligand (str): ligand name in pymol, or pdb file path, default None
        - complex (str): complex name in pymol, or pdb file path, default None
        - receptor_chain (str): receptor chain id, default None, only used when complex is provided
        - ligand_chain (str): ligand chain id, default None, only used when complex is provided
        - taskpool (TaskPool): TaskPool object for parallelization, default None
        - cutoff (float): cutoff distance, default 4.
    Returns:
        interactions (dict):
            - key: interaction type, includes: Hydrogen Bonds, Non-bonded Interactions
            - value: list of interaction info: ((receptor_resi, receptor_resn, receptor_chain), (ligand_resi, ligand_resn, ligand_chain), distance)
    """
    def load_mol(mol):
        if os.path.exists(mol) and os.path.isfile(mol):
            tmp_name = f'LAZYDOCK_TMP_OBJ_{uuid4()}'
            cmd.load(mol, tmp_name)
            chain = cmd.get_chains(tmp_name)[0]
            return tmp_name, chain
        elif mol in cmd.get_names():
            return mol, cmd.get_chains(mol)[0]
        else:
            raise ValueError(f'{mol} not found in pymol or file path')
    # load receptor and ligand, make complex
    if complex is not None:
        if receptor_chain is None or ligand_chain is None:
            raise ValueError('receptor_chain and ligand_chain must be provided when complex is provided')
        if (os.path.exists(complex) and os.path.isfile(complex)) or (complex in cmd.get_names()):
            complex_pdbstr = opts_file(complex_name)
        else:
            raise ValueError(f'{complex} not found in pymol or file path')
    elif receptor is not None and ligand is not None:
        receptor_name, receptor_chain = load_mol(receptor)
        ligand_name, ligand_chain = load_mol(ligand)
        complex_name = f'LAZYDOCK_TMP_OBJ_{uuid4()}'
        cmd.select(complex_name, f'{receptor_name} or {ligand_name}')
        complex_pdbstr = cmd.get_pdbstr(complex_name)
        cmd.delete(complex_name) # this do not delete receptor and ligand
    else:
        raise ValueError('either complex or receptor and ligand must be provided')
    # delete temporary objects if loaded from file
    for mol in [receptor_name, ligand_name]:
        if mol.startswith('LAZYDOCK_TMP_OBJ_'):
            cmd.delete(mol)
    # run LigPlus
    if taskpool is None:
        return _run_ligplus_for_complex(ligplus_dir, complex_pdbstr, receptor_chain, ligand_chain, mode, cutoff, force_cwd)
    else:
        return taskpool.add_task(None, _run_ligplus_for_complex, ligplus_dir, complex_pdbstr, receptor_chain, ligand_chain, mode, cutoff, force_cwd)


def merge_interaction_df(interaction: Dict[str, List[Tuple[Tuple[int, str, str], Tuple[int, str, str], float]]],
                         interaction_df: pd.DataFrame, distance_cutoff: float):
    """merge the interactions returned by calcu_atom_level_interactions to interaction_df."""
    # index format: CHAIN_ID:RESI:RESN
    for interaction_type, values in interaction.items():
        for single_inter in values:
            # single_inter: ((217, 'VAL', 'A'), (10, 'PHE', 'Z'), 3.71)
            receptor_res = f'{single_inter[0][2]}:{single_inter[0][0]}:{single_inter[0][1]}'
            ligand_res = f'{single_inter[1][2]}:{single_inter[1][0]}:{single_inter[1][1]}'
            points = distance_cutoff - single_inter[2]
            if ligand_res not in interaction_df.index or receptor_res not in interaction_df.columns:
                interaction_df.loc[ligand_res, receptor_res] = points
            elif np.isnan(interaction_df.loc[ligand_res, receptor_res]):
                interaction_df.loc[ligand_res, receptor_res] = points
            else:
                interaction_df.loc[ligand_res, receptor_res] += points
    return interaction_df


SUPPORTED_MODE = ['Hydrogen Bonds', 'Non-bonded Interactions']


def calcu_receptor_poses_interaction(receptor: str, poses: List[str], ligplus_dir: Optional[str] = None,
                                     mode: Union[str, List[str]] = 'Hydrogen Bonds', cutoff: float = 4.,
                                     taskpool: TaskPool = None, verbose: bool = False, force_cwd: bool = False, **kwargs):
    """
    calcu interactions between one receptor and one ligand with many poses using LigPlus.
    Parameters:
        - receptor (str): receptor pymol name
        - poses (List[str]):  ligand pymol names
        - ligplus_dir (str): LigPlus directory, which contains LigPlus.jar
        - mode (str or List[str]): supported modes: Hydrogen Bonds, Non-bonded Interactions, default Hydrogen Bonds, can be set to 'all' to include all modes.
        - taskpool (TaskPool): TaskPool object for parallelization, default None
        - cutoff (float): cutoff distance, default 4
        
    Returns:
        interactions (dict):
            - key: ligand name
            - value: interactions dict between receptor and ligand
                - key: interaction type, includes: Hydrophobic Interactions, Hydrogen Bonds, Water Bridges, Salt Bridges, pi-Stacking, pi-Cation Interactions, Halogen Bonds, Metal Complexes
                - value: list of interaction info: ((receptor_resi, receptor_resn, receptor_chain), (ligand_resi, ligand_resn, ligand_chain), distance)
        
        interaction_df (pd.DataFrame): , interactions between receptor and ligand, in the format of ligand-residue-residue matrix, with the value of each cell is the interaction score between two atoms.
            interaction_df.loc[ligand_res, receptor_res] = score
    """
    # set mode
    if mode == 'all':
        mode = SUPPORTED_MODE
    # set ligplus_dir
    if ligplus_dir is None:
        ligplus_dir = config.configs.named_paths['ligplus_dir']
        if ligplus_dir is None:
            raise ValueError('LigPlus directory not provided and not found in config.json')
    # prepare interactions
    all_interactions, interaction_df = {}, pd.DataFrame()
    # calcu for each ligand
    for ligand in tqdm(poses, desc=f'LigPlus', disable=not verbose):
        # calcu interaction
        all_interactions[ligand] = run_ligplus(ligplus_dir, receptor, ligand, mode=mode, cutoff=cutoff, taskpool=taskpool, force_cwd=force_cwd)
        # wait for taskpool
        while taskpool is not None and taskpool.count_waiting_tasks() > 0:
            time.sleep(0.1)
    # merge interactions by res
    for ligand in all_interactions:
        if taskpool is not None:
            all_interactions[ligand] = taskpool.query_task(all_interactions[ligand], True, 10**8)
        # merge interactions by res
        merge_interaction_df(all_interactions[ligand], interaction_df, cutoff)
    if not interaction_df.empty:
        # sort res
        interaction_df.sort_index(axis=0, inplace=True, key=sort_func)
        interaction_df.sort_index(axis=1, inplace=True, key=sort_func)
        interaction_df.fillna(0, inplace=True)
    else:
        return None, None
    return all_interactions, interaction_df


if __name__ == '__main__':
    # dev code
    from lazydock.pml.autodock_utils import DlgFile
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'RECEPTOR')
    dlg = DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=True, parse2std=True)
    dlg.sort_pose()
    pose_lst = []
    for i, pose in enumerate(dlg.pose_lst[:10]):
        pose_lst.append(f'ligand_{i}')
        cmd.read_pdbstr(pose.as_pdb_string(), pose_lst[-1])
        cmd.alter(f'ligand_{i}', 'type="HETATM"')
    result = run_ligplus(config.configs.named_paths['ligplus_dir'],
                         receptor='RECEPTOR', ligand=pose_lst[0], mode='Hydrogen Bonds')
    taskpool = TaskPool('threads', n_worker=4).start()
    interactions, interaction_df = calcu_receptor_poses_interaction('RECEPTOR', pose_lst, cutoff=4, taskpool=None)
    print(interactions)