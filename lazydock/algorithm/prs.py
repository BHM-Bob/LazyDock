
# Perform PRS calculations given and MD trajectory and a final state
# co-ordinate file

# Author: David Penkler
# Date: 17-11-2016

# speed optimzed by LazyDock: BHM-Bob G
# Date: 02-03-2025

from math import ceil, floor, log10

import numpy as np
from MDAnalysis import Universe
from tqdm import tqdm

if __name__ == "__main__":
    from lazydock.algorithm.rms import batch_fit_to, fit_to, rmsd
else:
    from .rms import batch_fit_to, fit_to, rmsd


def mean(a: np.ndarray, axis: int, keepdim: bool, backend: str = 'numpy'):
    if backend == 'numpy':
        return np.mean(a, axis=axis, keepdims=keepdim)
    else:
        return a.mean(dim=axis, keepdim=keepdim)


def repeat(a: np.ndarray, repeats: int, axis: int, backend: str = 'numpy'):
    if backend == 'numpy':
        return np.repeat(a, repeats, axis=axis)
    else:
        n_repeats = [1]*len(a.shape)
        n_repeats[axis] = repeats
        return a.repeat(*n_repeats)


def run_prs(coords: np.ndarray, initial_pose: np.ndarray, final_pose: np.ndarray,
            perturbations: int = 250, backend: str = 'numpy'):
    sum_frames, n_residues, _ = coords.shape
    if backend!= 'numpy':
        import torch as _backend
    else:
        _backend = np
    # center coords
    initial_pose -= mean(initial_pose, 0, True, backend)
    final_pose -= mean(final_pose, 0, True, backend)
    coords -= mean(coords, 1, True, backend)
    batch_first_frame = repeat(coords[0:1, :, :], len(coords), 0, backend)
    # align frames to first frame and calcu average pose
    coords, _ = batch_fit_to(coords, batch_first_frame, backend=backend)
    average_pose = coords.mean(0)
    # align frames to average pose
    stop_align = False
    for _ in range(10):
        batch_average_pose = repeat(average_pose[None, :, :], len(coords), 0, backend)
        coords, _ = batch_fit_to(coords, batch_average_pose, backend=backend)
        if stop_align:
            break
        new_average_pose = coords.mean(0)
        now_rmsd = rmsd(average_pose, new_average_pose, None, True, True, backend)
        average_pose = new_average_pose
        if now_rmsd <= 0.000001:
            stop_align = True
    # calculate corr_mat: [N*3, N*3]
    R_mat = coords.reshape(sum_frames, n_residues*3)
    corr_mat = (R_mat.T @ R_mat) / (sum_frames-1)
    # Implementing perturbations sequentially: diffP: [N*3, N*3, perturbations]
    diffP = _backend.zeros((n_residues, n_residues*3, perturbations), dtype=coords.dtype)
    if backend != 'numpy':
        diffP = diffP.to(diffE.device)
    batch_initial_pose = repeat(initial_pose.reshape(1, -1), n_residues, 0, backend)
    for s in tqdm(range(perturbations), total=perturbations, desc='perform perturbations', leave=False):
        # make random noise: rand_mat: [N, N*3]
        rand_mat = []
        for _ in range(3):
            r_i = _backend.zeros((n_residues, n_residues))
            if backend != 'numpy':
                r_i = r_i.to(coords.device)
            np.fill_diagonal(r_i, 2*np.random.rand(n_residues)-1)
            rand_mat.append(r_i)
        rand_mat = np.stack(rand_mat, axis=0).transpose(1, 2, 0).reshape(n_residues, 3*n_residues)
        # apply noise to initial_trans: perturbed_pose = [N, N*3] @ [N*3, N*3] + [N*3, 1] => [N, N*3]
        perturbed_pose = rand_mat @ corr_mat + batch_initial_pose
        # fit perturbed_pose to initial_pose
        diffP[:, :, s] = batch_fit_to(perturbed_pose.reshape(n_residues, n_residues, 3),
                                      batch_initial_pose.reshape(n_residues, n_residues, 3),
                                      backend=backend)[0].reshape(n_residues, -1)
    # calculate experimental difference： diffE： [N, 3]
    final_pose, _ = fit_to(final_pose, initial_pose, final_pose.mean(0), initial_pose.mean(0), backend=backend)
    diffE = (final_pose - initial_pose).reshape(n_residues, 3)
    # calculate pearson's coefficient
    ## 计算DTarget的向量化版本
    DTarget = np.linalg.norm(diffE, axis=1)
    ## 计算DIFF的向量化版本
    diffP_reshaped = diffP.reshape(n_residues, n_residues, 3, perturbations)
    DIFF = np.linalg.norm(diffP_reshaped, axis=2).transpose(1, 0, 2)
    # 计算RHO的向量化版本
    ## 重组DIFF为二维矩阵便于批量计算
    reshaped_diff = DIFF.transpose(1, 2, 0).reshape(-1, n_residues)
    dt_centered = DTarget - DTarget.mean()
    ## 批量计算协方差和标准差
    diff_centered = reshaped_diff - reshaped_diff.mean(axis=1, keepdims=True)
    covariances = (diff_centered @ dt_centered) / (n_residues - 1)
    std_devs = diff_centered.std(axis=1, ddof=1) * dt_centered.std(ddof=1)
    ## 避免除以零（假设数据无零标准差情况）
    max_RHO: np.ndarray = (covariances / std_devs).reshape(n_residues, perturbations).max(axis=-1)
    return max_RHO



if __name__ == "__main__":
    # extract coords
    from pathlib import Path
    from mbapy_lite.file import opts_file
    
    paths = opts_file(Path(__file__).parent.parent.parent / 'data_tmp/test_config.json', way='json')['test_rms']
    u = Universe(paths['top_path'], paths['traj_path'])
    start, stop, step = 0, 10000, 1
    ag = u.select_atoms('protein and name CA')
    stop = stop or len(u.trajectory)
    sum_frames = ceil((stop - start) / step)
    coords = np.zeros((sum_frames, len(ag), 3), dtype=np.float64)
    for current, _ in enumerate(tqdm(u.trajectory[start:stop:step],
                                     desc='Gathering coordinates', total=sum_frames, leave=False)):
        coords[current] = ag.positions.astype(np.float64)
    run_prs(coords, coords[0].copy(), coords[-1].copy(), 250)
