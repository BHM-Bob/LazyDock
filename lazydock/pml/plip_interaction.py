'''
Date: 2024-10-11 10:33:10
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-05 20:08:43
Description: 
'''
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex
from pymol import cmd

if __name__ == '__main__':
    from lazydock.pml.interaction_utils import sort_func
    from lazydock.utils import uuid4
else:
    from ..utils import uuid4
    from .interaction_utils import sort_func
    
    
def get_atom_level_interactions(mol, receptor_chain: str, ligand_chain: str, mode: List[str], cutoff: float = 4.):
    """
    """
    interactions = {}
    for ligand, interaction in mol.interaction_sets.items():
        info = BindingSiteReport(interaction)
        for name, feat, values in [['Hydrophobic Interactions', info.hydrophobic_features, info.hydrophobic_info],
                                  ['Hydrogen Bonds', info.hbond_features, info.hbond_info],
                                  ['Water Bridges', info.waterbridge_features, info.waterbridge_info],
                                  ['Salt Bridges', info.saltbridge_features, info.saltbridge_info],
                                  ['pi-Stacking', info.pistacking_features, info.pistacking_info],
                                  ['pi-Cation Interactions', info.pication_features, info.pication_info],
                                  ['Halogen Bonds', info.halogen_features, info.halogen_info],
                                  ['Metal Complexes', info.metal_features, info.metal_info]]:
            interactions.setdefault(name, [])
            for value in values:
                find_idx_fn = lambda x: feat.index(list(filter(lambda y: x == y, feat))[0])
                rec_idx, lig_idx = find_idx_fn('RESNR'), find_idx_fn('RESNR_LIG')
                dist_term = 'DIST'
                if name == 'Hydrogen Bonds':
                    dist_term = 'DIST_H-A'
                elif name == 'Water Bridges':
                    dist_term = 'DIST_A-W'
                elif name == 'pi-Stacking':
                    dist_term = 'CENTDIST'
                dist_idx = find_idx_fn(dist_term)
                rec_res, lig_res, dist = value[rec_idx:rec_idx+3], value[lig_idx:lig_idx+3], float(value[dist_idx])
                if dist <= cutoff and rec_res[-1] == receptor_chain and lig_res[-1] == ligand_chain:
                    interactions[name].append((rec_res, lig_res, dist))
    return interactions


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


SUPPORTED_MODE = ['Hydrophobic Interactions', 'Hydrogen Bonds', 'Water Bridges', 'Salt Bridges', 'pi-Stacking', 'pi-Cation Interactions', 'Halogen Bonds', 'Metal Complexes']


def calcu_receptor_poses_interaction(receptor: str, poses: List[str], mode: Union[str, List[str]] = 'all', cutoff: float = 4., **kwargs):
    """
    calcu interactions between one receptor and one ligand with many poses using PLIP-python.
    Parameters:
        - receptor (str): receptor pymol name
        - poses (List[str]):  ligand pymol names
        - mode (Union[str, List[str]]): interaction types to be calculated, default 'all' to include all modes., 
            supported: Hydrophobic Interactions, Hydrogen Bonds, Water Bridges, Salt Bridges, pi-Stacking, pi-Cation Interactions, Halogen Bonds, Metal Complexes
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
    elif isinstance(mode, str) and mode in SUPPORTED_MODE:
        mode = [mode]
    elif any(m not in SUPPORTED_MODE for m in mode):
        raise ValueError(f'Unsupported mode: {mode}, supported: {SUPPORTED_MODE}')
    # else: mode = mode
    # prepare interactions
    receptor_chain = cmd.get_chains(receptor)[0]
    all_interactions, interaction_df = {}, pd.DataFrame()
    # calcu for each ligand
    for ligand in poses:
        ligand_chain = cmd.get_chains(ligand)[0]
        # calcu interaction
        sele_complex = uuid4()
        cmd.select(sele_complex, f'{ligand} or {receptor}')
        mol = PDBComplex()
        mol.load_pdb(cmd.get_pdbstr(sele_complex), as_string=True)
        try:
            mol.analyze()
            all_interactions[ligand] = get_atom_level_interactions(mol, receptor_chain, ligand_chain, mode, cutoff)
            # merge interactions by res
            merge_interaction_df(all_interactions[ligand], interaction_df, mode, cutoff)
        except Exception as e:
            print(f'Error in {ligand}: {e}, skip this pose.')
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
    calcu_receptor_poses_interaction('RECEPTOR', pose_lst)