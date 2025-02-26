'''
Date: 2025-02-20 10:49:33
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-27 16:07:13
Description: 
'''
import numpy as np
from tqdm import tqdm


def inner_product(coords1, coords2, weights=None, backend = None):
    """Calculate the weighted inner product matrix and the E0 value.计算加权内积矩阵和E0值

    Parameters:
        coords1 (np.ndarray): The first set of coordinates, shape (n_atoms, 3).
        coords2 (np.ndarray): The second set of coordinates, shape (n_atoms, 3).
        weights (np.ndarray, optional): Weights for each atom, shape (n_atoms,).
        backend (module, optional): The backend to use for calculations, default is numpy.

    Returns:
        tuple: A flattened inner product matrix A and the E0 value.
    """
    backend = backend or np
    if weights is not None:
        w = weights[:, np.newaxis]
        A = (coords1 * w).T @ coords2
        G1 = backend.sum((coords1 * w) * coords1)
        G2 = backend.sum(weights * np.sum(coords2**2, axis=1))
    else:
        A = coords1.T @ coords2
        G1 = backend.sum(coords1**2)
        G2 = backend.sum(coords2**2)
    return A.ravel(), 0.5 * (G1 + G2)

def fast_calc_rmsd_rotation(rot, A_flat, E0, N):
    """Calculate RMSD and rotation matrix based on the inner product matrix.
    基于内积矩阵快速计算RMSD和旋转矩阵

    Parameters:
        rot (np.ndarray or None): Output array for the flattened rotation matrix, shape (9,). If None, only return RMSD.
        A_flat (np.ndarray): Flattened inner product matrix, shape (9,).
        E0 (float): Precomputed value E0.
        N (int): Number of atoms.

    Returns:
        float or tuple: If rot is None, return only the RMSD value. Otherwise, return a tuple of (RMSD, flattened rotation matrix).
    """
    # 构造4x4关键矩阵
    Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz = A_flat
    K = np.array([
        [Sxx + Syy + Szz, Syz - Szy,      Szx - Sxz,      Sxy - Syx],
        [Syz - Szy,      Sxx - Syy - Szz, Sxy + Syx,      Szx + Sxz],
        [Szx - Sxz,      Sxy + Syx,      -Sxx + Syy - Szz, Syz + Szy],
        [Sxy - Syx,      Szx + Sxz,       Syz + Szy,      -Sxx - Syy + Szz]
    ])
    
    # 计算最大特征值和对应特征向量
    eigenvalues, eigenvectors = np.linalg.eigh(K)
    max_idx = np.argmax(eigenvalues)
    max_eigen = eigenvalues[max_idx]
    quat = eigenvectors[:, max_idx]
    
    # 计算RMSD
    rmsd = np.sqrt(max(0.0, 2.0 * (E0 - max_eigen) / N))
    
    if rot is None:
        return rmsd
    
    # 四元数转旋转矩阵
    q1, q2, q3, q4 = quat / np.linalg.norm(quat)
    rot_matrix = np.array([
        [q1**2 + q2**2 - q3**2 - q4**2, 2*(q2*q3 - q1*q4),     2*(q2*q4 + q1*q3)],
        [2*(q2*q3 + q1*q4),     q1**2 - q2**2 + q3**2 - q4**2, 2*(q3*q4 - q1*q2)],
        [2*(q2*q4 - q1*q3),     2*(q3*q4 + q1*q2),     q1**2 - q2**2 - q3**2 + q4**2]
    ])
    
    # 展平旋转矩阵到输出数组
    rot[:] = rot_matrix.ravel()
    return rmsd, rot

def calc_rms_rotational_matrix(ref, conf, rot=None, weights=None):
    """Calculate RMSD and rotation matrix between two coordinate sets.
    
    Parameters:
        ref (np.ndarray): Reference coordinates, shape (n_atoms, 3)
        conf (np.ndarray): Target coordinates, shape (n_atoms, 3)
        rot (np.ndarray, optional): Output array for rotation matrix, shape (9,)
        weights (np.ndarray, optional): Weights for each atom, shape (n_atoms,)
        
    Returns:
        tuple or float: If rot is None, returns RMSD only. Otherwise returns (RMSD, rotation matrix)
    """
    A, E0 = inner_product(ref, conf, weights)
    return fast_calc_rmsd_rotation(rot, A, E0, ref.shape[0])


def batch_inner_product(batch_coords1, batch_coords2, weights=None):
    """Batch calculation of inner product matrix and E0 value for multiple frames.
    
    Args:
        batch_coords1 (np.ndarray): Reference coordinates [n_frames, n_atoms, 3]
        batch_coords2 (np.ndarray): Target coordinates [n_frames, n_atoms, 3]
        weights (np.ndarray, optional): Weights [n_atoms,] or [n_frames, n_atoms]
        
    Returns:
        tuple: (A_flat, E0) where:
            - A_flat: Flattened inner product matrix [n_frames, 9]
            - E0: Precomputed value [n_frames,]
    """
    if weights is not None:
        weights = np.asarray(weights)
        if weights.ndim == 1:
            weights = weights[np.newaxis, :, np.newaxis]  # 广播到所有帧
        else:
            weights = weights[:, :, np.newaxis]
        
        # 加权内积计算
        weighted_coords1 = batch_coords1 * weights
        G1 = np.sum(weighted_coords1 * batch_coords1, axis=(1,2))
        G2 = np.sum(weights * np.sum(batch_coords2**2, axis=2), axis=1)
        A = np.einsum('fai,faj->fij', weighted_coords1, batch_coords2)
    else:
        G1 = np.sum(batch_coords1**2, axis=(1,2))
        G2 = np.sum(batch_coords2**2, axis=(1,2))
        A = np.einsum('fai,faj->fij', batch_coords1, batch_coords2)
    
    A_flat = A.reshape(A.shape[0], 9)
    E0 = 0.5 * (G1 + G2)
    return A_flat, E0

def batch_fast_calc_rmsd(batch_rot, A_flat, E0, n_atoms):
    """
    批量快速计算RMSD和旋转矩阵的核心算法
    :param batch_rot: 输出旋转矩阵 [n_frames, 9]
    :param A_flat: 内积矩阵 [n_frames, 9]
    :param E0: 预计算值 [n_frames,]
    :param n_atoms: 原子数
    :return: RMSD数组 [n_frames,]
    """
    n_frames = A_flat.shape[0]
    S = A_flat.reshape(n_frames, 3, 3)
    
    # 构造4x4关键矩阵K [n_frames, 4, 4]
    K = np.zeros((n_frames, 4, 4))
    K[:, 0, 0] = S[:, 0, 0] + S[:, 1, 1] + S[:, 2, 2]
    K[:, 0, 1] = S[:, 1, 2] - S[:, 2, 1]
    K[:, 0, 2] = S[:, 2, 0] - S[:, 0, 2]
    K[:, 0, 3] = S[:, 0, 1] - S[:, 1, 0]
    K[:, 1, 0] = K[:, 0, 1]
    K[:, 1, 1] = S[:, 0, 0] - S[:, 1, 1] - S[:, 2, 2]
    K[:, 1, 2] = S[:, 0, 1] + S[:, 1, 0]
    K[:, 1, 3] = S[:, 2, 0] + S[:, 0, 2]
    K[:, 2, 0] = K[:, 0, 2]
    K[:, 2, 1] = K[:, 1, 2]
    K[:, 2, 2] = -S[:, 0, 0] + S[:, 1, 1] - S[:, 2, 2]
    K[:, 2, 3] = S[:, 1, 2] + S[:, 2, 1]
    K[:, 3, 0] = K[:, 0, 3]
    K[:, 3, 1] = K[:, 1, 3]
    K[:, 3, 2] = K[:, 2, 3]
    K[:, 3, 3] = -S[:, 0, 0] - S[:, 1, 1] + S[:, 2, 2]
    
    # 批量特征值分解
    eigenvalues, eigenvectors = np.linalg.eigh(K)
    max_eigenvalues = eigenvalues[:, -1]  # 取最大特征值
    quaternions = eigenvectors[:, :, -1]  # 对应特征向量
    
    # 计算RMSD
    rmsd = np.sqrt(np.maximum(0.0, 2.0 * (E0 - max_eigenvalues) / n_atoms))
    
    if batch_rot is None:
        return rmsd
    
    # 批量四元数转旋转矩阵
    q = quaternions / np.linalg.norm(quaternions, axis=1, keepdims=True)
    q0, q1, q2, q3 = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    
    batch_rot[:, 0] = q0**2 + q1**2 - q2**2 - q3**2
    batch_rot[:, 1] = 2*(q1*q2 - q0*q3)
    batch_rot[:, 2] = 2*(q1*q3 + q0*q2)
    batch_rot[:, 3] = 2*(q1*q2 + q0*q3)
    batch_rot[:, 4] = q0**2 - q1**2 + q2**2 - q3**2
    batch_rot[:, 5] = 2*(q2*q3 - q0*q1)
    batch_rot[:, 6] = 2*(q1*q3 - q0*q2)
    batch_rot[:, 7] = 2*(q2*q3 + q0*q1)
    batch_rot[:, 8] = q0**2 - q1**2 - q2**2 + q3**2
    
    return rmsd

def batch_calc_rmsd(batch_ref, batch_conf, batch_rot=None, weights=None):
    """
    批量计算RMSD和旋转矩阵
    :param batch_ref: 参考坐标 [n_frames, n_atoms, 3]
    :param batch_conf: 目标坐标 [n_frames, n_atoms, 3]
    :param batch_rot: 输出旋转矩阵 [n_frames, 9]
    :param weights: 权重 [n_atoms,] 或 [n_frames, n_atoms]
    :return: RMSD数组 [n_frames,]
    """
    A_flat, E0 = batch_inner_product(batch_ref, batch_conf, weights)
    return batch_fast_calc_rmsd(batch_rot, A_flat, E0, batch_ref.shape[1])


def fit_to(mobile_coordinates, ref_coordinates, mobile_com, ref_com, weights=None):
    r"""Perform an rmsd-fitting to determine rotation matrix and align atoms

    Parameters
    ----------
    mobile_coordinates : ndarray
        Coordinates of atoms to be aligned
    ref_coordinates : ndarray
        Coordinates of atoms to be fit against
    mobile_com: ndarray
        array of xyz coordinate of mobile center of mass
    ref_com : ndarray
        array of xyz coordinate of reference center of mass
    weights : array_like (optional)
       choose weights. With ``None`` weigh each atom equally. If a float array
       of the same length as `mobile_coordinates` is provided, use each element
       of the `array_like` as a weight for the corresponding atom in
       `mobile_coordinates`.

    Returns
    -------
    mobile_atoms : AtomGroup
        AtomGroup of translated and rotated atoms
    min_rmsd : float
        Minimum rmsd of coordinates
    """
    rot = np.zeros(9, dtype=np.float64)
    min_rmsd, R = calc_rms_rotational_matrix(ref_coordinates, mobile_coordinates,
                                             rot,  weights=weights)
    mobile_coordinates = mobile_coordinates.copy() - mobile_com
    mobile_coordinates = np.matmul(mobile_coordinates, R.reshape(3, 3))
    mobile_coordinates += ref_com
    return mobile_coordinates, min_rmsd


def rmsd(a, b, weights=None, center=False, superposition=False, backend: str = 'numpy'):
    r"""Returns RMSD between two coordinate sets `a` and `b`.

    `a` and `b` are arrays of the coordinates of N atoms of shape
    :math:`N times 3` as generated by, e.g.,
    :meth:`MDAnalysis.core.groups.AtomGroup.positions`.

    Parameters
    ----------
    a : array_like
        coordinates to align to `b`
    b : array_like
        coordinates to align to (same shape as `a`)
    weights : array_like (optional)
        1D array with weights, use to compute weighted average
    center : bool (optional)
        subtract center of geometry before calculation. With weights given
        compute weighted average as center.
    superposition : bool (optional)
        perform a rotational and translational superposition with the fast QCP
        algorithm [Theobald2005]_ before calculating the RMSD; implies
        ``center=True``.

    Returns
    -------
    rmsd : float
        RMSD between `a` and `b`
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    N = b.shape[0]
    if a.shape != b.shape:
        raise ValueError('a and b must have same shape')

    # superposition only works if structures are centered
    if center or superposition:
        # make copies (do not change the user data!)
        # weights=None is equivalent to all weights 1
        a = a - np.average(a, axis=0, weights=weights)
        b = b - np.average(b, axis=0, weights=weights)

    if weights is not None:
        if len(weights) != len(a):
            raise ValueError('weights must have same length as a and b')
        # weights are constructed as relative to the mean
        weights = np.asarray(weights, dtype=np.float64) / np.mean(weights)

    if superposition:
        rot = np.zeros(9, dtype=np.float64)
        return calc_rms_rotational_matrix(a, b, rot, weights)[0] #calc_rmsd_and_rotation(a, b, weights)
    else:
        if weights is not None:
            return np.sqrt(np.sum(weights[:, np.newaxis]
                                  * ((a - b) ** 2)) / N)
        else:
            return np.sqrt(np.sum((a - b) ** 2) / N)


def pairwise_rmsd(traj: np.ndarray, traj2: np.ndarray = None, block_size: int = 100,
                  backend: str = 'numpy', verbose: bool = False):
    """Calculate pairwise RMSD between all frames in a trajectory.
    
    Args:
        traj (np.ndarray): Trajectory data with shape (n_frames, n_atoms, 3)
        traj2 (np.ndarray, optional): Second trajectory for cross-comparison, if None, use traj itself
        block_size (int): Number of frames to process in each block
        backend (str): Computation backend ('numpy', 'torch', or 'cuda')
        verbose (bool): Whether to show progress bar
        
    Returns:
        np.ndarray: Symmetric RMSD matrix with shape (n_frames, n_frames)
    """
    traj2 = traj2 or traj
    n_frames, n_atoms, _ = traj.shape
    K = n_frames // block_size
    if backend == 'numpy':
        rmsd_matrix = np.zeros((n_frames, n_frames), dtype=np.float32)
    elif backend in {'torch', 'cuda'}:
        import torch
        traj = torch.from_numpy(traj)
        traj2 = torch.tensor(traj2)
        rmsd_matrix = torch.zeros((n_frames, n_frames), dtype=torch.float32)
    # calcu rmsd for each block
    for i in tqdm(range(K), total=K, desc='Calculating RMSD matrix', leave=False, disable=not verbose):
        # prepare block i
        start_i = i * block_size
        end_i = (i + 1) * block_size if i < K - 1 else n_frames
        block_i = traj[start_i:end_i]
        if backend == 'cuda':
            block_i = block_i.cuda()
        # calcu rmsd for block-i series
        for j in range(K):
            # prepare block j
            start_j = j * block_size
            end_j = (j + 1) * block_size if j < K - 1 else n_frames
            block_j = traj2[start_j:end_j]
            # calculate RMSD
            if backend == 'numpy':
                diff = block_i[:, np.newaxis] - block_j[np.newaxis]
                rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=-1), axis=-1))
            elif backend == 'torch':
                diff = block_i[:, None] - block_j[None]
                rmsd = torch.sqrt(torch.mean(torch.sum(diff ** 2, dim=-1), dim=-1))
            elif backend == 'cuda':
                diff = block_i[:, None] - block_j[None].cuda()
                rmsd = torch.sqrt(torch.mean(torch.sum(diff ** 2, dim=-1), dim=-1)).cpu()
            # fill rmsd matrix
            rmsd_matrix[start_i:end_i, start_j:end_j] = rmsd

    if backend in {'torch', 'cuda'}:
        rmsd_matrix = rmsd_matrix.numpy()
    return rmsd_matrix