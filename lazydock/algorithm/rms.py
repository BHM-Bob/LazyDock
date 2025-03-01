'''
Date: 2025-02-20 10:49:33
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-01 11:10:14
Description: 
'''
import numpy as np
from tqdm import tqdm


def inner_product(coords1, coords2, weights=None, backend: str = 'numpy'):
    """Calculate the weighted inner product matrix and the E0 value.计算加权内积矩阵和E0值

    Parameters:
        coords1 (np.ndarray): The first set of coordinates, shape (n_atoms, 3).
        coords2 (np.ndarray): The second set of coordinates, shape (n_atoms, 3).
        weights (np.ndarray, optional): Weights for each atom, shape (n_atoms,).
        backend (module, optional): The backend to use for calculations, default is numpy, support torch and cuda.

    Returns:
        tuple: A flattened inner product matrix A and the E0 value.
    """
    # dermine the backend
    if backend in {'torch', 'cuda'}:
        import torch
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # calculate the inner product matrix and E0 value
    if weights is not None:
        w = weights[:, np.newaxis]
        A = (coords1 * w).T @ coords2
        G1 = _backend.sum((coords1 * w) * coords1)
        G2 = _backend.sum(weights * np.sum(coords2**2, axis=1))
    else:
        A = coords1.T @ coords2
        G1 = _backend.sum(coords1**2)
        G2 = _backend.sum(coords2**2)
    return A.ravel(), 0.5 * (G1 + G2)

def fast_calc_rmsd_rotation(rot, A_flat, E0, N, backend: str = 'numpy'):
    """Calculate RMSD and rotation matrix based on the inner product matrix.
    基于内积矩阵快速计算RMSD和旋转矩阵

    Parameters:
        rot (np.ndarray or None): Output array for the flattened rotation matrix, shape (9,). If None, only return RMSD.
        A_flat (np.ndarray): Flattened inner product matrix, shape (9,).
        E0 (float): Precomputed value E0.
        N (int): Number of atoms.
        backend (module, optional): The backend to use for calculations, default is numpy, support torch and cuda.

    Returns:
        float or tuple: If rot is None, return only the RMSD value. Otherwise, return a tuple of (RMSD, flattened rotation matrix).
    """
    # dermine the backend
    if backend in {'torch', 'cuda'}:
        import torch
    arr_fn = torch.tensor if backend  in {'torch', 'cuda'} else np.array
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # 构造4x4关键矩阵
    Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz = A_flat
    K = arr_fn([
        [Sxx + Syy + Szz, Syz - Szy,      Szx - Sxz,      Sxy - Syx],
        [Syz - Szy,      Sxx - Syy - Szz, Sxy + Syx,      Szx + Sxz],
        [Szx - Sxz,      Sxy + Syx,      -Sxx + Syy - Szz, Syz + Szy],
        [Sxy - Syx,      Szx + Sxz,       Syz + Szy,      -Sxx - Syy + Szz]
    ])
    
    # 计算最大特征值和对应特征向量
    eigenvalues, eigenvectors = _backend.linalg.eigh(K)
    max_idx = _backend.argmax(eigenvalues)
    max_eigen = eigenvalues[max_idx]
    quat = eigenvectors[:, max_idx]
    
    # 计算RMSD
    rmsd = _backend.sqrt(max(0.0, 2.0 * (E0 - max_eigen) / N))
    
    if rot is None:
        return rmsd
    
    # 四元数转旋转矩阵
    q1, q2, q3, q4 = quat / _backend.linalg.norm(quat)
    rot_matrix = arr_fn([
        [q1**2 + q2**2 - q3**2 - q4**2, 2*(q2*q3 - q1*q4),     2*(q2*q4 + q1*q3)],
        [2*(q2*q3 + q1*q4),     q1**2 - q2**2 + q3**2 - q4**2, 2*(q3*q4 - q1*q2)],
        [2*(q2*q4 - q1*q3),     2*(q3*q4 + q1*q2),     q1**2 - q2**2 - q3**2 + q4**2]
    ])
    
    # 展平旋转矩阵到输出数组
    rot[:] = rot_matrix.ravel()
    return rmsd, rot

def calc_rms_rotational_matrix(ref, conf, rot=None, weights=None, backend: str = 'numpy'):
    """Calculate RMSD and rotation matrix between two coordinate sets.
    
    Parameters:
        ref (np.ndarray): Reference coordinates, shape (n_atoms, 3)
        conf (np.ndarray): Target coordinates, shape (n_atoms, 3)
        rot (np.ndarray, optional): Output array for rotation matrix, shape (9,)
        weights (np.ndarray, optional): Weights for each atom, shape (n_atoms,)
        backend (str, optional): Backend to use for computation, support 'numpy', 'torch', 'cuda'. Default is 'numpy'.
        
    Returns:
        tuple or float: If rot is None, returns RMSD only. Otherwise returns (RMSD, rotation matrix)
    """
    A, E0 = inner_product(ref, conf, weights, backend)
    return fast_calc_rmsd_rotation(rot, A, E0, ref.shape[0], backend)


def batch_inner_product(batch_coords1, batch_coords2, weights=None, backend: str = 'numpy'):
    """Batch calculation of inner product matrix and E0 value for multiple frames.
    
    Args:
        batch_coords1 (np.ndarray): Reference coordinates [n_frames, n_atoms, 3]
        batch_coords2 (np.ndarray): Target coordinates [n_frames, n_atoms, 3]
        weights (np.ndarray, optional): Weights [n_atoms,] or [n_frames, n_atoms]
        backend (str, optional): Backend to use for computation, support 'numpy', 'torch', 'cuda'. Default is 'numpy'.
        
    Returns:
        tuple: (A_flat, E0) where:
            - A_flat: Flattened inner product matrix [n_frames, 9]
            - E0: Precomputed value [n_frames,]
    """
    # dermine the backend
    if backend in {'torch', 'cuda'}:
        import torch
    arr_fn = torch.tensor if backend  in {'torch', 'cuda'} else np.array
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # calculate the inner product matrix and E0 value
    if weights is not None:
        if weights.ndim == 1:
            weights = weights[None, :, None]  # 广播到所有帧
        else:
            weights = weights[:, :, None]
        
        # 加权内积计算
        weighted_coords1 = batch_coords1 * weights
        G1 = _backend.sum(weighted_coords1 * batch_coords1, (1,2))
        G2 = _backend.sum(weights * _backend.sum(batch_coords2**2, 2), 1)
        A = _backend.einsum('fai,faj->fij', weighted_coords1, batch_coords2)
    else:
        G1 = _backend.sum(batch_coords1**2, (1,2))
        G2 = _backend.sum(batch_coords2**2, (1,2))
        A = _backend.einsum('fai,faj->fij', batch_coords1, batch_coords2)
    
    A_flat = A.reshape(A.shape[0], 9)
    E0 = 0.5 * (G1 + G2)
    return A_flat, E0

def batch_fast_calc_rmsd(batch_rot, A_flat, E0, n_atoms, backend: str = 'numpy'):
    """Batch calculation of RMSD and rotation matrix based on the inner product matrix.
    
    Parameters:
        batch_rot (np.ndarray or None): Output array for the flattened rotation matrix, shape (n_frames, 9). If None, only return RMSD.
        A_flat (np.ndarray): Flattened inner product matrix, shape (n_frames, 9).
        E0 (np.ndarray): Precomputed value E0, shape (n_frames,).
        n_atoms (int): Number of atoms.
        backend (str): Backend to use for computation, support 'numpy', 'torch', 'cuda'. Default is 'numpy'.
        
    Returns:
        np.ndarray or tuple: If batch_rot is None, returns RMSD only. Otherwise returns (RMSD, rotation matrix)
    """
    # dermine the backend
    if backend in {'torch', 'cuda'}:
        import torch
    arr_fn = torch.tensor if backend  in {'torch', 'cuda'} else np.array
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # 展平旋转矩阵到输出数组
    n_frames = A_flat.shape[0]
    S = A_flat.reshape(n_frames, 3, 3)
    # 构造4x4关键矩阵K [n_frames, 4, 4]
    if backend in {'torch', 'cuda'}:
        K = torch.zeros((n_frames, 4, 4), device=A_flat.device, dtype=A_flat.dtype)
    else:
        K = np.zeros((n_frames, 4, 4), dtype=A_flat.dtype)
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
    eigenvalues, eigenvectors = _backend.linalg.eigh(K)
    max_eigenvalues = eigenvalues[:, -1]  # 取最大特征值
    quaternions = eigenvectors[:, :, -1]  # 对应特征向量
    # 计算RMSD
    zero = _backend.zeros(1, device=A_flat.device, dtype=A_flat.dtype) if backend in {'torch', 'cuda'} else np.zeros(1, dtype=A_flat.dtype)
    rmsd = _backend.sqrt(_backend.clip(2.0 * (E0 - max_eigenvalues) / n_atoms, zero, None))
    if batch_rot is None:
        return rmsd
    # 批量四元数转旋转矩阵
    q = quaternions / _backend.linalg.norm(quaternions, axis=1, keepdims=True)
    q0, q1, q2, q3 = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    # 展平旋转矩阵到输出数组
    batch_rot[:, 0] = q0**2 + q1**2 - q2**2 - q3**2
    batch_rot[:, 1] = 2*(q1*q2 - q0*q3)
    batch_rot[:, 2] = 2*(q1*q3 + q0*q2)
    batch_rot[:, 3] = 2*(q1*q2 + q0*q3)
    batch_rot[:, 4] = q0**2 - q1**2 + q2**2 - q3**2
    batch_rot[:, 5] = 2*(q2*q3 - q0*q1)
    batch_rot[:, 6] = 2*(q1*q3 - q0*q2)
    batch_rot[:, 7] = 2*(q2*q3 + q0*q1)
    batch_rot[:, 8] = q0**2 - q1**2 - q2**2 + q3**2
    return rmsd, batch_rot

def batch_calc_rmsd_rotational_matrix(batch_ref, batch_conf, batch_rot=None, weights=None, backend: str = 'numpy'):
    """perform RMSD calculation between two sets of coordinates.
    
    Parameters:
        batch_ref (np.ndarray): Reference coordinates, shape (n_frames, n_atoms, 3)
        batch_conf (np.ndarray): Target coordinates, shape (n_frames, n_atoms, 3)
        batch_rot (np.ndarray, optional): Output array for rotation matrix, shape (n_frames, 9). Defaults to None.
        weights (np.ndarray, optional): Weights for each atom, shape (n_atoms,) or (n_frames, n_atoms). Defaults to None.
        backend (str, optional): The backend to use for calculations, default is 'numpy', supports 'numpy' and 'cuda'. Defaults to 'numpy'.
        
    Returns:
        np.ndarray or tuple: If batch_rot is None, returns RMSD only. Otherwise returns (RMSD, rotation matrix)
    """
    A_flat, E0 = batch_inner_product(batch_ref, batch_conf, weights, backend)
    return batch_fast_calc_rmsd(batch_rot, A_flat, E0, batch_ref.shape[1], backend)


def fit_to(mobile_coordinates, ref_coordinates, mobile_com, ref_com, weights=None, backend='numpy'):
    """Perform an rmsd-fitting to determine rotation matrix and align atoms

    Parameters
        mobile_coordinates, ref_coordinates : ndarray: [n_atoms, 3]
            Coordinates of atoms to be aligned
        mobile_com, ref_com: ndarray: [3,]
            array of xyz coordinate of center of mass
        weights : array_like (optional)
            choose weights. With ``None`` weigh each atom equally. If a float array
            of the same length as `mobile_coordinates` is provided, use each element
            of the `array_like` as a weight for the corresponding atom in
            `mobile_coordinates`.
        backend : str (optional)
            The backend to use for calculations, default is 'numpy', supports 'numpy' and 'cuda'. Defaults to 'numpy'.

    Returns
        mobile_coords : ndarray: [n_atoms, 3]
            AtomGroup of translated and rotated atoms
        min_rmsd : float
            Minimum rmsd of coordinates
    """
    # determine backend
    if backend in {'torch', 'cuda'}:
        import torch
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # calculate rmsd and rotation matrix
    if backend == 'cuda':
        rot = torch.zeros(9, dtype=torch.float64, device='cuda')
    else:
        rot = _backend.zeros(9, dtype=_backend.float64)
    min_rmsd, R = calc_rms_rotational_matrix(ref_coordinates, mobile_coordinates,
                                             rot,  weights=weights, backend=backend)
    # apply rotation matrix to mobile_coordinates
    mobile_coordinates = mobile_coordinates.copy() - mobile_com
    mobile_coordinates = _backend.matmul(mobile_coordinates, R.reshape(3, 3))
    mobile_coordinates += ref_com
    return mobile_coordinates, min_rmsd


def rmsd(a, b, weights=None, center=True, superposition=True, backend: str = 'numpy'):
    """Returns RMSD between two coordinate sets `a` and `b`.

    Parameters
        a, b : array_like: [n_atoms, 3]
            coordinates to align, a is the reference, b is the mobile.
        weights : array_like (optional)
            1D array with weights, use to compute weighted average
        center : bool (optional)
            subtract center of geometry before calculation. With weights given
            compute weighted average as center.
        superposition : bool (optional)
            perform a rotational and translational superposition with the fast QCP
            algorithm [Theobald2005]_ before calculating the RMSD; implies
            ``center=True``.
        backend : str (optional)
            backend to use, default is 'numpy', support 'torch' and 'cuda'

    Returns
        rmsd : float
            RMSD between `a` and `b`
    """
    # determine backend
    if backend in {'torch', 'cuda'}:
        import torch
    arr_fn = torch.tensor if backend  in {'torch', 'cuda'} else np.array
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # check input
    N = b.shape[0]
    if a.shape != b.shape:
        raise ValueError('a and b must have same shape')
    # superposition only works if structures are centered
    if center or superposition:
        if backend == 'numpy':
            a = a - np.average(a, 0, weights=weights)
            b = b - np.average(b, 0, weights=weights)
        else:
            if weights is None:
                a = a - torch.mean(a, 0)
                b = b - torch.mean(b, 0)
            else:
                a = a*weights[:, None] - torch.mean(a, 0, keepdim=True)
                b = b*weights[:, None] - torch.mean(b, 0, keepdim=True)
    # check weights
    if weights is not None:
        if len(weights) != len(a):
            raise ValueError('weights must have same length as a and b')
        # weights are constructed as relative to the mean
        weights = arr_fn(weights, dtype=np.float64) / np.mean(weights)
    # calculate RMSD
    if superposition:
        if backend == 'cuda':
            rot = torch.zeros(9, dtype=torch.float64, device='cuda')
        else:
            rot = _backend.zeros(9, dtype=_backend.float64)
        return calc_rms_rotational_matrix(a, b, rot, weights, backend=backend)[0] #calc_rmsd_and_rotation(a, b, weights)
    else:
        if weights is not None:
            return _backend.sqrt(_backend.sum(weights[:, None]
                                  * ((a - b) ** 2)) / N)
        else:
            return _backend.sqrt(_backend.sum((a - b) ** 2) / N)
        
        
def batch_rmsd(a: np.ndarray, b: np.ndarray, backend: str = 'numpy'):
    """
    batch calculate RMSD between two coordinate sets `a` and `b`, with center and superposition.
    
    Parameters
    ----------
    a, b : array_like: [n_frames, n_atoms, 3]
        coordinates to align, a is the reference, b is the mobile.
    backend : str, optional
        backend to use, default is 'numpy', support 'torch' and 'cuda'
    
    Returns
    -------
    rmsd : array_like: [n_frames]
        RMSD between `a` and `b`
    """
    # determine backend
    if backend in {'torch', 'cuda'}:
        import torch
    _backend = torch if backend  in {'torch', 'cuda'} else np
    # check input
    if a.shape != b.shape:
        raise ValueError('a and b must have same shape')
    # center
    keepdim_kwg = {'keepdims' if backend == 'numpy' else 'keepdim': True}
    a = a - _backend.mean(a, 1, **keepdim_kwg)
    b = b - _backend.mean(b, 1, **keepdim_kwg)
    # perform superposition and calc RMSD
    if backend in {'torch', 'cuda'}:
        rot = torch.zeros(a.shape[0], 9, dtype=torch.float64, device=a.device)
    else:
        rot = np.zeros((a.shape[0], 9), dtype=np.float64)
    return batch_calc_rmsd_rotational_matrix(a, b, rot, weights=None, backend=backend)[0]


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