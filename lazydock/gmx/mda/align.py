'''
Date: 2025-03-05 21:56:36
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-04-02 19:59:50
Description: 
'''
import numpy as np
from MDAnalysis import AtomGroup, Universe
from tqdm import tqdm

if __name__ == '__main__':
    from lazydock.algorithm.rms import batch_fit_to
else:
    from ...algorithm.rms import batch_fit_to


def get_aligned_coords(u: Universe, ag: AtomGroup, start: int, step: int, stop: int,
                       ref_coords: np.ndarray = None, backend: str = "numpy",
                       verbose: bool = True, return_rmsd: bool = False) -> np.ndarray:
    """
    Align the coordinates of an AtomGroup to a reference coordinates array.
    
    Parameters
    ----------
    u : Universe
        The Universe object.
    ag : AtomGroup
        The AtomGroup to align.
    start, step, stop : int
        The start, step, and stop indices for the trajectory.
    ref_coords : np.ndarray: [N, 3] or [1, N, 3], optional, default=None
        The reference coordinates array, if None, the first frame of the trajectory will be used.
    backend : str, optional, default="numpy"
        The backend to use for the alignment. Options are "numpy", "torch", "cuda"
    verbose : bool, optional, default=True
        Whether to print progress information.
        
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]: ([n_frames, n_atoms, 3], [n_frames, n_atoms, 3])
        The original and aligned coordinates array.
    OR 
    Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray]]: ([n_frames, n_atoms, 3], ([n_frames, n_atoms, 3], [n_frames,]))
    """
    sum_frames = (len(u.trajectory) if stop is None else stop) - start
    coords = []
    for _ in tqdm(u.trajectory[start:stop:step], total=sum_frames//step,
                  desc='Gathering coordinates', leave=False, disable=not verbose):
        coords.append(ag.positions.copy().astype(np.float64))
    coords = np.array(coords, dtype=np.float64)
    if backend != 'numpy':
        import torch
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        coords = torch.tensor(coords, device=device, dtype=torch.float64)
    # do center as MDAnalysis.align.AlignTraj._single_frame
    coords: np.ndarray = coords - coords.mean(axis=1, keepdims=True)
    # copy original coords
    ori_coords = coords.copy() if backend == 'numpy' else coords.clone()
    # do center as MDAnalysis.align.AlignTraj._prepare, but coords[0] already centered
    if ref_coords is None:
        ref_coords = coords[0][None, :, :].copy() if backend == 'numpy' else coords[0][None, :, :].clone()
    # check ref_coords's shape
    if len(ref_coords.shape) == 2:
        if ref_coords.shape[0] == coords.shape[1] and ref_coords.shape[1] == 3:
            ref_coords = ref_coords[None, ]
        else:
            raise ValueError(f'Invalid ref_coords shape: {ref_coords.shape}')
    if backend == 'numpy':
        refs = np.repeat(ref_coords, coords.shape[0], axis=0)
    else:
        refs = ref_coords.repeat(coords.shape[0], 1, 1)
    # do fit as MDAnalysis.align.AlignTraj._single_frame, it use _fit_to
    if return_rmsd:
        return ori_coords, batch_fit_to(coords, refs, backend=backend)
    return ori_coords, batch_fit_to(coords, refs, backend=backend)[0]