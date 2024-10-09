'''
Date: 2024-09-30 19:28:57
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-10-09 09:29:34
Description: RRCS calculation in PyMOL, RRCS is from article "Common activation mechanism of class A GPCRs": https://github.com/elifesciences-publications/RRCS/blob/master/RRCS.py
'''
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import scipy
from pymol import cmd

def calcu_RRCS(model: str):
    """
    calcu RRCS score for each residue in the model.
    
    Parameters:
        - model: str, pymol name of the model to calculate RRCS score, must be loaded in PyMOL.
        
    Returns:
    contact_df: pd.DataFrame, contact_df[ires][jres] = contact_score
    """
    dict_coord = {} # dict to store coordinates. dict_coord[res][atom] = (x, y, z, occupancy, name)
    cmd.iterate_state(1, model, 'dict_coord.setdefault(f"{chain}:{resi}:{resn}", {}).setdefault(index, (x, y, z, q, name, resi))', space={'dict_coord': dict_coord})
    arr_coord = np.array([dict_coord[ires][iatom][:3] for ires in dict_coord for iatom in dict_coord[ires]])
    arr_n2info = np.array([[ires, iatom[-3], iatom[-2], iatom[-1]] for ires in dict_coord for iatom in dict_coord[ires].values()])
    q2_mat = arr_n2info[:, 1].astype(float).reshape(-1, 1) * arr_n2info[:, 1].astype(float).reshape(1, -1)
    # calcu d1 and d2 for each atom pair
    d1_mat = scipy.spatial.distance.cdist(arr_coord, arr_coord, metric = 'cityblock')
    d2_mat = scipy.spatial.distance.cdist(arr_coord, arr_coord, metric = 'euclidean') ** 2
    # make condition mat
    ## find close atom, if atom is close to other atom with d1 < 4.63
    is_close_atom_mat = (d1_mat < 4.63).astype(int)
    ## find residue pair, if atom is in same residue with other atom
    is_res_pair_mat = (arr_n2info[:, 0].reshape(-1, 1) != arr_n2info[:, 0].reshape(1, -1)).astype(int)
    ## find close res, if is, use backbone atom check: abs(ires_num - jres_num) < 5
    is_close_res = (np.abs(arr_n2info[:, 3].astype(int) - arr_n2info[:, 3].astype(int).reshape(1, -1)) < 5).astype(int)
    is_close_res_mat = (is_close_res.reshape(-1, 1) & is_close_res.reshape(1, -1)).astype(int)
    is_far_res_mat = 1 - is_close_res_mat
    ## find side atom, for backbone atom check, N, CA, C, O are backbone atoms, others are side atoms
    is_bb_atom = (arr_n2info[:, 2] == 'N') | (arr_n2info[:, 2] == 'CA') | (arr_n2info[:, 2] == 'C') | (arr_n2info[:, 2] == 'O')
    is_side_atom_mat = 1 - (is_bb_atom.reshape(-1, 1) & is_bb_atom.reshape(1, -1)).astype(int)
    ## calcu distance mat
    is_middle_d2 = (d2_mat >= 10.4329) & (d2_mat <= 21.4369).astype(int)
    is_small_d2 = (d2_mat < 10.4329).astype(int)
    # calcu score mat
    condition_mat = is_close_atom_mat * is_res_pair_mat * (is_close_res_mat*is_side_atom_mat + is_far_res_mat)
    score_mat = condition_mat * (is_middle_d2*1.0 + is_small_d2*(1-(d2_mat**0.5 - 3.23)/1.4)) * q2_mat
    # score_mat = np.triu(score_mat, k=0) # keep half of the score mat
    # convert dict to DataFrame
    contact_df = pd.DataFrame(data=score_mat,
                              index=list(arr_n2info[:, 0]), columns=list(arr_n2info[:, 0]))
    return contact_df

if __name__ == '__main__':
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    from mbapy.base import TimeCosts
    @TimeCosts(3)
    def test_calcu(idx):
        calcu_RRCS('receptor')
    test_calcu()