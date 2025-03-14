'''
Date: 2025-03-14 17:15:36
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-14 18:03:26
Description: 
'''
from typing import Union

import numpy as np
from MDAnalysis.analysis.data.filenames import Janin_ref, Rama_ref


def count_rama_favor(angles: np.ndarray, ref = 'rama'):
    """
    count thr outer, marginally allowed, and favored Ramachandran angles
    
    Parameters:
        angles: numpy.ndarray, shape (n_frames, n_res, 2), the phi and psi angles of each frame
        ref: str or numpy.ndarray, the reference data for counting, can be 'rama' or 'janin' or a numpy.ndarray
    
    Returns:
        counts: numpy.ndarray, shape (n_frames, 3), the counts of outer, marginally allowed, and favored Ramachandran angles for each frame
    """
    # 将角度转换为参考数据的索引
    phi = angles[..., 0]
    psi = angles[..., 1]
    # 计算对应参考网格的索引
    if isinstance(ref, str):
        if ref == 'rama':
            phi_idx = ((phi + 180) / 4).astype(int)
            psi_idx = ((psi + 180) / 4).astype(int)
            ref = np.load(Rama_ref)
            level = [1, 17, 15000] # from MDAnalysis/analysis/dihedrals.py
        elif ref == 'janin':
            phi_idx = (phi / 6).astype(int)
            psi_idx = (psi / 6).astype(int)
            ref = np.load(Janin_ref)
            level = [1, 6, 600] # from MDAnalysis/analysis/dihedrals.py
        else:
            raise ValueError(f'ref must be "rama" or "janin", but got {ref}')
    
    # 确保索引在合法范围内 [0, 89]
    phi_idx = np.clip(phi_idx, 0, 89)
    psi_idx = np.clip(psi_idx, 0, 89)
    
    # 注意psi在前，phi在后
    density = ref[psi_idx, phi_idx]
    
    # 创建分类掩码
    outer = (density < level[0])
    marginally_allowed = (density >= level[0]) & (density <= level[1])
    favored = (density > level[1])
    
    # 统计每个frame的分类数量
    counts = np.stack([
        outer.sum(axis=1),
        marginally_allowed.sum(axis=1),
        favored.sum(axis=1)
    ], axis=1)
    
    return counts