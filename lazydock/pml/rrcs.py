'''
Date: 2024-09-30 19:28:57
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-14 17:22:05
Description: RRCS calculation in PyMOL, RRCS is from article "Common activation mechanism of class A GPCRs": https://github.com/elifesciences-publications/RRCS/blob/master/RRCS.py
'''
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from pymol import cmd


def _test_close(dict_coord: Dict[str, Dict[int, Tuple[float, float, float, float]]],
                ires: str, jres: str):
    for iatom in dict_coord[ires]:
        (ix, iy, iz, iocc) = dict_coord[ires][iatom]
        for jatom in dict_coord[jres]:                  
            (jx, jy, jz, jocc) = dict_coord[jres][jatom]
            dx = abs(ix-jx)
            dy = abs(iy-jy)
            dz = abs(iz-jz)
            if dx < 4.63 and dy < 4.63 and dz < 4.63:
                return True
    return False

def _calcu_score(dict_coord: Dict[str, Dict[int, Tuple[float, float, float, float]]],
                 atomnum2name: Dict[int, str], ires: str, jres: str,
                 check_hetatm: bool, score_count: int):
    total_score = 0
    for iatom in dict_coord[ires]:
        if check_hetatm and atomnum2name[iatom] in ['N', 'CA', 'C', 'O']:
            continue
        (ix, iy, iz, iocc) = dict_coord[ires][iatom]
        for jatom in dict_coord[jres]:
            if check_hetatm and atomnum2name[jatom] in ['N', 'CA', 'C', 'O']:
                continue
            (jx, jy, jz, jocc) = dict_coord[jres][jatom]
            d2 = (ix-jx)**2 + (iy-jy)**2 + (iz-jz)**2
            if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                score = 0
            elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                score = 1.0*iocc*jocc
            else:
                score = (1-(d2**0.5 - 3.23)/1.4)*iocc*jocc
            total_score += score
            score_count[0] += 1
    return total_score


def calcu_RRCS_from_dict(dict_coord: Dict[str, Dict[int, Tuple[float, float, float, float]]],
                        atomnum2name: Dict[int, str]):
    """
    Parameters:
        - dict_coord: dict of coordinates, dict_coord[res][atom] = (x, y, z, occupancy)
        - atomnum2name: map atom number to atom name, in order to find N, CA, C, O
    Returns:
        contact_df: DataFrame of contact scores, index and columns are residue names.
    """
    contact_score = {} # dict to store final results. contact_score[ires][jres] = contact_score.
    score_count = [0] # 66320
    # calcu RRCS score for each residue pair
    for ires in dict_coord:
        ires_num = int(ires[ires.find(':')+1:ires.rfind(':')])
        contact_score[ires] = {}
        # find jres if it has any atom close to ires
        for jres in dict_coord:
            jres_num = int(jres[jres.find(':')+1:jres.rfind(':')])
            contact_score[ires][jres] = 0
            # skip because alreadly calculated
            if jres_num <= ires_num:
                continue
            # skip if jres has no atom close to ires
            if not _test_close(dict_coord, ires, jres):
                continue
            # calculate RRCS score for ires and jres
            contact_score[ires][jres] = _calcu_score(dict_coord, atomnum2name, ires, jres,
                                                     check_hetatm=abs(ires_num - jres_num) < 5, score_count=score_count)
    # convert dict to DataFrame
    contact_df = pd.DataFrame(data={k:list(v.values()) for k,v in contact_score.items()},
                              index=contact_score.keys(), columns=contact_score.keys())
    return contact_df


def calcu_RRCS_from_array(names: np.ndarray, resis: np.ndarray, positions: np.ndarray, occupancies: np.ndarray):
    """
    Parameters:
        - names: [N, ], np.ndarray of atom names
        - resns: [N, ], np.ndarray of residue names
        - resis: [N, ], np.ndarray of residue numbers
        - positions: [N, 3], np.ndarray of atom positions, shape (n_atoms, 3)
        - occupancies: [N, ], np.ndarray of atom occupancies
        
    Returns:
        contact_df: DataFrame of contact scores, index and columns are residue names.
    """
    # calcu mask
    check_atm_mask = np.abs(resis[:, None] - resis[None, :]) < 5
    atm_mask = np.in1d(names, np.array(['N', 'CA', 'C', 'O']))
    atm_mask = atm_mask[:, None] | atm_mask[None, :]
    ij_less_mask = (resis[:, None] > resis[None, :]).astype(np.int32)
    # calcu matrix
    d1_mat = np.abs(positions[:, None, :] - positions[None, :, :])# [N, N, 3]
    d2_mat = ((positions[:, None, :] - positions[None, :, :])**2).sum(axis=-1) # [N, N]
    close_mask = (d1_mat < 4.63).all(axis=-1, keepdims=False).astype(np.int32) # [N, N]
    ooc_mat = occupancies[:, None] * occupancies[None, :]
    # calculate RRCS score for each atom pair
    score_mat = np.where(d2_mat < 10.4329, ooc_mat, 0.0)
    score_mat = np.where((d2_mat >= 10.4329) & (d2_mat < 21.4369), (1-(d2_mat**0.5 - 3.23)/1.4) * ooc_mat, score_mat)
    score_mat = np.where(check_atm_mask & atm_mask, 0, score_mat * close_mask * ij_less_mask)
    # 获取分组的唯一值及其索引, row_labels is unique index
    unique_values, row_labels = np.unique(resis, return_inverse=True)
    K = len(unique_values)
    idx_mat = (row_labels[:, None] * K + row_labels[None, :]).reshape(-1)
    matrix_flat = score_mat.ravel()
    result_flat = np.bincount(idx_mat, weights=matrix_flat, minlength=K**2)
    # 重新排布为 K x K 矩阵
    result = result_flat.reshape(K, K)
    contact_df = pd.DataFrame(result, index=unique_values, columns=unique_values)
    # contact_df.to_excel('RRCS.xlsx')
    return contact_df


def calcu_RRCS_from_tensor(names: np.ndarray, resis: np.ndarray, positions: np.ndarray, occupancies: np.ndarray,
                           device: str = 'cpu'):
    """
    Parameters:
        - names: [N, ], np.ndarray of atom names
        - resns: [N, ], torch.tensor of residue names
        - resis: [N, ], torch.tensor of residue numbers
        - positions: [N, 3], torch.tensor of atom positions, shape (n_atoms, 3)
        - occupancies: [N, ], torch.tensor of atom occupancies
        - device: 'cpu' or 'cuda'
        
    Returns:
        contact_df: DataFrame of contact scores, index and columns are residue names.
    """
    import torch
    # calcu mask
    check_atm_mask = torch.abs(resis[:, None] - resis[None, :]) < 5
    atm_mask = torch.tensor(np.in1d(names, np.array(['N', 'CA', 'C', 'O'])), device=device)
    atm_mask = atm_mask[:, None] | atm_mask[None, :]
    ij_less_mask = (resis[:, None] > resis[None, :]).to(dtype=torch.int8)
    # calcu matrix
    d1_mat = torch.abs(positions[:, None, :] - positions[None, :, :])# [N, N, 3]
    d2_mat = ((positions[:, None, :] - positions[None, :, :])**2).sum(dim=-1) # [N, N]
    close_mask = (d1_mat < 4.63).all(dim=-1, keepdims=False).to(dtype=torch.int8) # [N, N]
    ooc_mat = occupancies[:, None] * occupancies[None, :]
    # calculate RRCS score for each atom pair
    score_mat = torch.where(d2_mat < 10.4329, ooc_mat, 0.0)
    score_mat = torch.where((d2_mat >= 10.4329) & (d2_mat < 21.4369), (1-(d2_mat**0.5 - 3.23)/1.4) * ooc_mat, score_mat)
    score_mat = torch.where(check_atm_mask & atm_mask, 0, score_mat * close_mask * ij_less_mask)
    # 获取分组的唯一值及其索引, row_labels is unique index
    unique_values, row_labels = torch.unique(resis, return_inverse=True)
    K = len(unique_values)
    idx_mat = (row_labels[:, None] * K + row_labels[None, :]).reshape(-1)
    matrix_flat = score_mat.ravel()
    result_flat = torch.bincount(idx_mat, weights=matrix_flat, minlength=K**2)
    # 重新排布为 K x K 矩阵
    result = result_flat.reshape(K, K)
    contact_df = pd.DataFrame(result.cpu().numpy(), index=unique_values.cpu().numpy(), columns=unique_values.cpu().numpy())
    # contact_df.to_excel('RRCS_tensor.xlsx')
    return contact_df


def calcu_RRCS(model: str, backend: str = 'numpy', _cmd = None, device: str = 'cpu'):
    """
    Parameters:
        - model: molecular name loaded in pymol
        - _cmd: pymol command object, default is cmd.

    Returns:
        contact_df: DataFrame of contact scores, index and columns are residue names.
    """
    _cmd = _cmd or cmd
    if backend in {'numpy', 'torch'}:
        # Prepare PyMOL data
        resis, names, positions, occupancies = [], [], [], []
        _cmd.iterate_state(1, model, 'resis.append(resi); names.append(name); '
                                'positions.append((x, y, z)); occupancies.append(q)',
                        space={'resis': resis, 'names': names, 'positions': positions, 
                                'occupancies': occupancies})
        resis = np.array(resis, dtype=int)
        if backend == 'torch':
            import torch
            resis = torch.tensor(resis, dtype=torch.int16, device=device)
            names = np.array(names)
            positions = torch.tensor(positions, dtype=torch.float32, device=device)
            occupancies = torch.tensor(occupancies, dtype=torch.float32, device=device)
            return calcu_RRCS_from_tensor(names, resis, positions, occupancies, device=device)
        else:
            names = np.array(names)
            positions = np.array(positions)
            occupancies = np.array(occupancies)
            return calcu_RRCS_from_array(names, resis, positions, occupancies)
    elif backend == 'dict':
        dict_coord = {} # dict to store coordinates. dict_coord[res][atom] = (x, y, z, occupancy)
        _cmd.iterate_state(1, model, 'dict_coord.setdefault(f"{chain}:{resi}:{resn}", {}).setdefault(index, (x, y, z, q))', space={'dict_coord': dict_coord})
        atomnum2name = {} # map atom number to atom name, in order to find N, CA, C, O
        _cmd.iterate(model, 'atomnum2name[index] = name', space={'atomnum2name': atomnum2name})
        return calcu_RRCS_from_dict(dict_coord, atomnum2name)
    else:
        raise ValueError('backend must be numpy or dict')


if __name__ == '__main__':
    import time
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    # calcu_RRCS('receptor', backend='torch')
    from mbapy_lite.base import TimeCosts
    @TimeCosts(6)
    def test_calcu(idx):
        calcu_RRCS('receptor', backend='numpy', device='cuda')
    test_calcu()