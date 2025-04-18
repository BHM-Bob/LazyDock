'''
Date: 2025-02-20 22:02:45
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-17 18:48:16
Description: 
'''
import numpy as np
from mbapy_lite.base import put_err
from MDAnalysis.core.groups import AtomGroup


def svd_hermitian(matrix):
    """calculate svd of a hermitian matrix with torch"""
    import torch

    # 特征分解
    eigenvalues, Q = torch.linalg.eigh(matrix)
    # 按绝对值降序排列
    sorted_abs, indices = torch.sort(eigenvalues.abs(), descending=True)
    indices = indices.to(eigenvalues.device)
    eigenvalues_sorted = eigenvalues[indices]
    Q_sorted = Q[:, indices]
    # 构造左奇异矩阵
    signs = torch.sign(eigenvalues_sorted)
    U = Q_sorted * signs.reshape(1, -1)
    # 奇异值
    S = eigenvalues_sorted.abs()
    # 右奇异矩阵的共轭转置
    Vh = Q_sorted.conj().T
    
    return U, S, Vh


def generate_ordered_pairs(positions: np.ndarray, cutoff: float, backend: str = 'numpy'):
    """
    Generate all ordered pairs of atoms within a cutoff distance using NumPy operations.

    Parameters
    ----------
    positions : ndarray
        Atom coordinates as an array of shape (n_atoms, 3)
    cutoff : float
        Distance threshold

    Returns
    -------
    list of tuples
        Pairs of atom indices (i, j) where distance is less than cutoff
    """
    # determine backend
    if backend == 'numpy':
        _backend = np
    else:
        import torch as _backend
    # Compute pairwise squared distances using broadcasting
    diff = positions[:, None, :] - positions[None, :, :]
    distance_sq = (diff ** 2).sum(-1)
    # Create a mask for distances below cutoff squared
    mask = distance_sq < cutoff ** 2
    # Extract the indices where the mask is True
    i, j = _backend.where(mask)
    return  _backend.concatenate([i[None, :], j[None, :]], 0).T


def generate_valid_paris(positions: np.ndarray, cutoff: float, backend: str = 'numpy'):
    cutoff_sq = cutoff ** 2
    # Generate all pairs from neighbour_generator
    pairs = generate_ordered_pairs(positions, cutoff, backend)
    if pairs.shape[0] == 0:
        return None
    i = pairs[:, 0]
    j = pairs[:, 1]
    # Filter pairs where i < j to avoid duplicates and excluding self-pairs
    mask = j > i
    i_filtered = i[mask]
    j_filtered = j[mask]
    # Calculate squared distances using NumPy's vectorized operations
    a = positions[i_filtered]
    b = positions[j_filtered]
    distance_squared = ((a - b) ** 2).sum(1)
    # Apply cutoff and get valid indices
    valid = distance_squared < cutoff_sq
    return i_filtered[valid], j_filtered[valid]
    

def generate_matrix(positions: np.ndarray, cutoff: float, backend: str = 'numpy'):
    # determine backend
    if backend == 'numpy':
        _backend = np
    else:
        import torch as _backend
    # compute pair index
    natoms = positions.shape[0]
    valid_pair = generate_valid_paris(positions, cutoff, backend)
    if backend == 'numpy':
        matrix = np.zeros((natoms, natoms), dtype=np.float64)
    else:
        matrix = _backend.zeros((natoms, natoms), dtype=_backend.float64, device=positions.device)
    if valid_pair is None:
        return matrix
    i_filtered, j_filtered = valid_pair
    # Create matrix and set symmetric entries
    matrix[i_filtered, j_filtered] = -1.0
    matrix[j_filtered, i_filtered] = -1.0
    # Calculate diagonal entries as the count of neighbors
    row_counts = (matrix < 0).sum(1)
    if backend == 'numpy':
        np.fill_diagonal(matrix, row_counts)
    else:
        matrix = _backend.diagonal_scatter(matrix, row_counts.to(dtype=_backend.float64))
    return matrix

        
def calcu_GNMAnalysis(positions: np.ndarray, cutoff: float = 7,
                      gen_matrix_fn = None, backend: str = 'numpy', **kwargs):
    """Generate the Kirchhoff matrix of contacts.

    This generates the neighbour matrix by generating a grid of
    near-neighbours and then calculating which are are within
    the cutoff.

    Returns
    -------
        eigenvectors
        eigenvalues
    """
    # determine backend
    if backend == 'numpy':
        _backend = np
    else:
        import torch as _backend
    # calculate matrix
    gen_matrix_fn = gen_matrix_fn or generate_matrix
    matrix = gen_matrix_fn(positions, cutoff, backend=backend, **kwargs)
    # calculate eigenvectors and eigenvalues
    try:
        if backend == 'numpy':
            _, w, v = np.linalg.svd(matrix, hermitian=True)
        else:
            _, w, v = svd_hermitian(matrix)
    except Exception as e:
        return put_err(f"SVD with cutoff {cutoff} failed to converge: {e}, return None")
    list_map = _backend.argsort(w)
    if backend == 'numpy':
        w = w[list_map[1]]
        v = v[list_map[1]]
    else:
        w = w[list_map[1]].cpu().numpy()
        v = v[list_map[1]].cpu().numpy()
    return w, v


def generate_close_matrix(positions: np.ndarray, cutoff,
                          atom2residue: np.ndarray, residue_size: np.ndarray,
                          n_residue: int, weights="size", backend: str = 'numpy'):
    """Generate the Kirchhoff matrix of closeContactGNMAnalysis contacts.

    This generates the neighbour matrix by generating a grid of
    near-neighbours and then calculating which are are within
    the cutoff.

    Returns
    -------
    array
            the resulting Kirchhoff matrix
    """
    # determine backend
    if backend == 'numpy':
        _backend = np
    else:
        import torch as _backend
    # Compute residue sizes
    if weights == 'size':
        inv_sqrt_res_sizes = 1.0 / _backend.sqrt(residue_size)
    else:
        if backend == 'numpy':
            inv_sqrt_res_sizes = np.ones(n_residue, dtype=np.float64)
        else:
            inv_sqrt_res_sizes = _backend.ones(n_residue, dtype=_backend.float64, device=positions.device)
    if backend == 'numpy':
        matrix = np.zeros((n_residue, n_residue), dtype=np.float64)
    else:
        matrix = _backend.zeros((n_residue, n_residue), dtype=_backend.float64, device=positions.device)
    # Generate all atom pairs within cutoff
    # Note: Using previous generate_ordered_pairs function (adjusted for pairs)
    valid_pair = generate_valid_paris(positions, cutoff, backend)
    if valid_pair is None:
        return matrix
    i_filtered, j_filtered = valid_pair
    # Get valid residue indices
    iresidues = atom2residue[i_filtered]
    jresidues = atom2residue[j_filtered]
    # Compute contact values
    contact = inv_sqrt_res_sizes[iresidues] * inv_sqrt_res_sizes[jresidues]
    # calculate the Kirkhoff matrix
    # because current index is res_idx from atom_idx, so it can be repeat for i or j
    # thus need to use index add to accumulate the number of contacts for each residue pair
    if backend == 'numpy':
        np.add.at(matrix, (iresidues, jresidues), -contact)
        np.add.at(matrix, (jresidues, iresidues), -contact)
        np.add.at(matrix, (iresidues, iresidues), contact)
        np.add.at(matrix, (jresidues, jresidues), contact)
    else:
        matrix.index_put_((iresidues, jresidues), -contact, accumulate=True)
        matrix.index_put_((jresidues, iresidues), -contact, accumulate=True)
        matrix.index_put_((iresidues, iresidues), contact, accumulate=True)
        matrix.index_put_((jresidues, jresidues), contact, accumulate=True)
    return matrix


def genarate_atom2residue(atoms: AtomGroup):
    """
    return
        - a 1d array where each element is the residue index of the atom
        - a 1d array where each element is the number of atoms in the residue
    """
    return atoms.resindices.copy(), np.array([r.atoms.n_atoms for r in atoms.residues])


def calcu_closeContactGNMAnalysis(positions: np.ndarray, cutoff: float, atom2residue: np.ndarray,
                                  residue_size: np.ndarray, n_residue: int, weights="size",
                                  backend: str = 'numpy'):
    """Generate the Kirchhoff matrix of contacts.

    This generates the neighbour matrix by generating a grid of
    near-neighbours and then calculating which are are within
    the cutoff.

    Returns
    -------
        eigenvectors
        eigenvalues
    """
    return calcu_GNMAnalysis(positions, cutoff, gen_matrix_fn=generate_close_matrix, backend=backend,
                             atom2residue=atom2residue, residue_size=residue_size,
                             n_residue=n_residue, weights=weights)


if __name__ == '__main__':
    # dev code
    import torch
    
    coords = torch.randn(1000, 3).cuda()
    # coords = np.random.randn(1000, 3)
    for _ in range(5):
        calcu_GNMAnalysis(coords, backend='cuda')