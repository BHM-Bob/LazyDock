'''
Date: 2024-08-18 12:56:06
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-25 21:29:04
Description: 
'''
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from mbapy.base import put_err
from pymol import CmdException, cmd

if __name__ == '__main__':
    from lazydock.utils import uuid4
else:
    from ..utils import uuid4

RES_SEPRATOR = '::'
NULL_CHAIN = "''"


def get_h_bond_cutoff():
    """ Get the distance cutoff (center and edge) for hydrogen bonds."""
     # ideal geometry and minimally acceptable geometry
    return cmd.get('h_bond_cutoff_center'), cmd.get('h_bond_cutoff_edge')


def set_h_bond_cutoff(center: float, edge: float):
    """ Set the distance cutoff (center and edge) for hydrogen bonds."""
    cmd.set('h_bond_cutoff_center', center)
    cmd.set('h_bond_cutoff_edge', edge)
    

def get_distance_info(dist_name: str, state = 1, selection = 'all', xyz2atom = None):
    '''
    Params:
        - dist_name : str, name of distance object
        - state : int, state of distance object
        - selection : str, selection of atoms to compute distances for
        - xyz2atom : dict, mapping of xyz coordinates to atom info: {(x,y,z): (model, chain, resn, resi, elem, index)}
        
    Returns:
        - list of tuples, each tuple contains two atom info and distance value, in the format of ((model, chain, resn, resi, elem, index), (model, chain, resn, resi, elem, index), distance)
    
    Notes:
        - the order of two sele in the return is not guaranteed to be the same as the order in the input selection
        
    modified from http://pymolwiki.org/index.php/get_raw_distances
    '''
    from chempy import cpv

    state = int(state)
    if state < 1:
        state = cmd.get_state()

    # Possible return values are

    # "object:molecule"
    # "object:map"
    # "object:mesh"
    # "object:slice"
    # "object:surface"
    # "object:measurement"
    # "object:cgo"
    # "object:group"
    # "object:volume"
    # "selection"
    if dist_name not in cmd.get_names_of_type('object:measurement'):
        return put_err('no such distance object: ' + dist_name)

    raw_objects = cmd.get_session(dist_name, 1, 1, 0, 0)['names']

    if xyz2atom is None:
        xyz2atom = {}
        cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model,chain,resn,resi,elem,index)',
                        space=locals())

    r = []
    for obj in raw_objects:
        try:
            points = obj[5][2][state - 1][1]
            if points is None:
                raise ValueError
        except (KeyError, ValueError):
            continue
        for i in range(0, len(points), 6):
            xyz1 = tuple(points[i:i + 3])
            xyz2 = tuple(points[i + 3:i + 6])
            try:
                r.append([xyz2atom[xyz1], xyz2atom[xyz2], cpv.distance(xyz1, xyz2)])
            except Exception as e:
                put_err(f'Error processing distance object {obj[0]}: {e}')
    return r


def calcu_atom_level_interactions(sele1: str, sele2: str, mode = 0, cutoff=3.6, xyz2atom: Dict = None):
    """
    calcu distances between atoms set sele1 and sele2 using cmd.distance, and get the distances and atom info.
    
    Parameters:
        - sele1: str, selection 1 of atoms
        - sele2: str, selection 2 of atoms
        - mode: int, mode of cmd.distance
        - cutoff: float, cutoff of cmd.distance
        - xyz2atom: dict, mapping of xyz coordinates to atom info: {(x,y,z): (model, chain, resn, resi, elem, index)}
        
    Returns:
        - atoms: dict, atom info for each selection, in the format of 
            {'s1_a': [(model, chain, resn, resi, elem, index),...],
             's1_d': [(model, chain, resn, resi, elem, index),...],
             's2_a': [(model, chain, resn, resi, elem, index),...],
             's2_d': [(model, chain, resn, resi, elem, index),...]}
        - xyz2atom: dict, mapping of xyz coordinates to atom info: {(x,y,z): (model, chain, resn, resi, elem, index)}
        - residues: dict, residue info for each selection, in the format of {'s1_a': {model: {chain: {resn-resi: elem}}...},'s1_d': {model: {chain: {resn-resi: elem}}...},'s2_a': {model: {chain: {resn-resi: elem}}...},'s2_d': {model: {chain: {resn-resi: elem}}...}}
        - interactions: dict, interactions between atoms in each selection, in the format of {'aa': [(('receptor', 'A', 'LYS', '108', 'O', 459), ('ligand', '', 'TYR', '1', 'N', 48), 3.4605595828662383),...], 'ad': [(('receptor', 'A', 'LYS', '108', 'O', 459), ('ligand', '', 'TYR', '1', 'N', 48), 3.4605595828662383),...], 'da': [(('receptor', 'A', 'LYS', '108', 'O', 459), ('ligand', '', 'TYR', '1', 'N',
    """
    interactions = {}
    rand_id = uuid4()
    atom_selections = {
        's1_a': f'{sele1}_a_{rand_id}',
        's1_d': f'{sele1}_d_{rand_id}',
        's2_a': f'{sele2}_a_{rand_id}',
        's2_d': f'{sele2}_d_{rand_id}',
    }
    dist_names = {
        'aa': f'aa_{rand_id}',
        'ad': f'ad_{rand_id}',
        'da': f'da_{rand_id}',
        'dd': f'dd_{rand_id}',
    }
    atoms = {K:[] for K in atom_selections.keys()}
    residues = {K:{} for K in atom_selections.keys()} # model:chain:resn-resi
    if xyz2atom is None:
        xyz2atom = {}
        cmd.iterate_state(1, 'all', 'xyz2atom[x,y,z] = (model,chain,resn,resi,elem,index)', space = locals())
        
    # make selections for each atom type in sele1 and sele2
    for atom_sele in atom_selections.values():
        sele = sele1 if atom_sele.startswith(sele1) else sele2
        atom_type1 = 'donor' if atom_sele == f'{sele}_d_{rand_id}' else 'acceptor'
        atom_type2 = 'acceptor' if atom_sele == f'{sele}_d_{rand_id}' else 'donor'
        cmd.select(atom_sele, f'{sele} and ({atom_type1} and !{atom_type2})')
    # get atom info in each selection
    for atom_type in atom_selections:
        cmd.iterate(atom_selections[atom_type], lambda atom: atoms[atom_type].append(atom))
    # ganther atoms to residues
    for atom_type in atom_selections:
        for atom in atoms[atom_type]:
            # if chain is null string, set it to ''
            residues[atom_type].setdefault(RES_SEPRATOR.join([atom.model, atom.chain or NULL_CHAIN, atom.resn, atom.resi]), atom.elem)
    
    # compute contacts between each pair of atoms in each selection
    for s1_atoms, s2_atoms, dist_name in [(atom_selections['s1_a'], atom_selections['s2_d'], dist_names['ad']),
                            (atom_selections['s1_d'], atom_selections['s2_a'], dist_names['da']),
                            (atom_selections['s1_a'], atom_selections['s2_a'], dist_names['aa']),
                            (atom_selections['s1_d'], atom_selections['s2_d'], dist_names['dd'])]:
        cmd.distance(dist_name, s1_atoms, s2_atoms, cutoff=cutoff, mode=mode)
        # get distances and build interactions
        interactions[dist_name[:2]] = get_distance_info(dist_name, xyz2atom=xyz2atom)
        
    # delete atom selctions
    for atom_sele in atom_selections.values():
        cmd.delete(atom_sele)
    # delete distance objects
    for dist_name in dist_names.values():
        cmd.delete(dist_name)
    
    return atoms, xyz2atom, residues, interactions


def sort_atom_level_interactions(interactions: Dict[str, List[Tuple[Tuple[str, str, str, str, str, float],
                                                         Tuple[str, str, str, str, str, float], float]]],
                                model1: str, model2: str):
    """sort the interactions returned by calcu_atom_level_interactions, sort the order in tuple by model1 and model2."""
    for ty in list(interactions.keys()):
        values = interactions[ty]
        for i in range(len(values)):
            if values[i][0][0] == model2:
                interactions[ty][i][0], interactions[ty][i][1] = interactions[ty][i][1], interactions[ty][i][0]
            elif values[i][0][0] != model1:
                return put_err(f'{values[i][0][0]} is not in {model1} and {model2}, abort sort', interactions)
    return interactions



def merge_interaction_df(interaction: Dict[str, List[Tuple[Tuple[str, str, str, str, str, float],
                                                            Tuple[str, str, str, str, str, float], float]]],
                         interaction_df: pd.DataFrame,
                         distance_cutoff: float, nagetive_factor: float):
    """merge the interactions returned by calcu_atom_level_interactions to interaction_df."""
    # index format: CHAIN_ID:RESI:RESN
    def set_points(ty: str, points: float, nagetive_factor: float):
        points = distance_cutoff - points
        if ty in {'ad', 'da'}:
            return points
        return nagetive_factor * points
    for interaction_type, values in interaction.items():
        for single_inter in values:
            # single_inter: (('receptor', 'A', 'LYS', '108', 'O', 459), ('ligand', '', 'TYR', '1', 'N', 48), 3.4605595828662383)
            receptor_res = f'{single_inter[0][1]}:{single_inter[0][3]}:{single_inter[0][2]}'
            ligand_res = f'{single_inter[1][1]}:{single_inter[1][3]}:{single_inter[1][2]}'
            points = set_points(interaction_type, single_inter[2], nagetive_factor)
            if ligand_res not in interaction_df.index or receptor_res not in interaction_df.columns:
                interaction_df.loc[ligand_res, receptor_res] = points
            elif np.isnan(interaction_df.loc[ligand_res, receptor_res]):
                interaction_df.loc[ligand_res, receptor_res] = points
            else:
                interaction_df.loc[ligand_res, receptor_res] += points
    return interaction_df


def calcu_receptor_poses_interaction(receptor: str, poses: List[str], mode: int = 0,
                                     cutoff: float = 4., nagetive_factor: float = -1.):
    """
    calcu interactions between one receptor and one ligand with many poses.
    
    Parameters:
        - receptor: str, receptor pymol name
        - poses: list of str, ligand pymol names
        - mode: int, mode of cmd.distance
        - cutoff: float, cutoff of cmd.distance
        - nagetive_factor: float, factor to multiply the distance value for interations between acceptor and acceptor, and donor and donor.
        
    Returns:
        - interactions: dict, interactions between receptor and ligand, in the format of {'ligand': [xyz2atom, residues, interactions]}, where xyz2atom is a dict, residues is a dict, interactions is a dict.
        - interaction_df: pd.DataFrame, interactions between receptor and ligand, in the format of ligand-residue-residue matrix, with the value of each cell is the distance between two atoms.
    """
    def sort_func(index: pd.Index):
        index = index.str.split(':')
        return list(map(lambda x: (x[0], int(x[1]), x[2]), index))
    # prepare interactions
    all_interactions, interaction_df = {}, pd.DataFrame()
    # select receptor
    sele_receptor = uuid4()
    cmd.select(sele_receptor, receptor)
    # calcu for each ligand
    for ligand in poses:
        # calcu interaction
        sele_ligand = uuid4()
        cmd.select(sele_ligand, ligand)
        _, xyz2atom, residues, interactions = calcu_atom_level_interactions(sele_receptor, sele_ligand,
                                                                            mode, cutoff)
        interactions = sort_atom_level_interactions(interactions, receptor, ligand)
        # NOTE: do not save atoms, because pymol.editing._AtomProxy class can't be pickled
        all_interactions[ligand] = [xyz2atom, residues, interactions]
        cmd.delete(sele_ligand)
        # merge interactions by res
        merge_interaction_df(interactions, interaction_df, cutoff, nagetive_factor)
    cmd.delete(sele_receptor)
    if not interaction_df.empty:
        # sort res
        interaction_df.sort_index(axis=0, inplace=True, key=sort_func)
        interaction_df.sort_index(axis=1, inplace=True, key=sort_func)
        interaction_df.fillna(0, inplace=True)
    else:
        return None, None
    return interactions, interaction_df


def filter_interaction_df(interaction_df: pd.DataFrame, colum_axis_min: float = None,
                          row_axis_min: float = None, inplace: bool = False):
    """
    filter the interaction_df by the minimum value of each axis.
    
    Parameters:
        - interaction_df : pd.DataFrame, interactions between receptor and ligand, in the format of ligand-residue-residue matrix, with the value of each cell is the distance between two atoms.
        - colum_axis_min : float, minimum value of each column axis, None means no filter.
        - row_axis_min : float, minimum value of each row axis, None means no filter.
        - inplace : bool, whether to filter the original interaction_df or return a deep copy.
        
    Returns:
        - interaction_df : pd.DataFrame, filtered interactions between receptor and ligand,
                    in the format of ligand-residue-residue matrix,
                    with the value of each cell is the distance between two atoms.
    """
    if inplace:
        tmp_interaction_df = interaction_df
    else:
        tmp_interaction_df = interaction_df.copy(deep = True)
    # filter ligand
    if row_axis_min is not None:
        for ligand_res in tmp_interaction_df.index:
            if tmp_interaction_df.loc[ligand_res].abs().max() < row_axis_min:
                tmp_interaction_df.drop(ligand_res, inplace=True)
    # filter receptor
    if colum_axis_min is not None:
        for receptor_res in tmp_interaction_df.columns:
            if tmp_interaction_df[receptor_res].abs().max() < colum_axis_min:
                tmp_interaction_df.drop(receptor_res, axis=1, inplace=True)
    return tmp_interaction_df



if __name__ == '__main__':
    # dev code
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'RECEPTOR')
    cmd.load('data_tmp/pdb/LIGAND.pdb', 'LIGAND')
    cmd.select('sele1', 'RECEPTOR')
    cmd.select('sele2', 'LIGAND')
    atoms, xyz2atom, residues, interactions = calcu_atom_level_interactions('sele1', 'sele2')
    calcu_receptor_poses_interaction('RECEPTOR', ['LIGAND'])