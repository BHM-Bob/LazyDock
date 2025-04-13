'''
Date: 2025-02-20 22:02:45
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-16 20:41:43
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


def generate_valid_paris(positions, cutoff):
    positions = np.asarray(positions)
    cutoff_sq = cutoff ** 2

    # Generate all pairs from neighbour_generator
    all_pairs = generate_ordered_pairs(positions, cutoff)
    
    if not all_pairs:
        return None
    
    pairs = np.array(all_pairs)
    i = pairs[:, 0]
    j = pairs[:, 1]

    # Filter pairs where i < j to avoid duplicates
    mask = j > i
    i_filtered = i[mask]
    j_filtered = j[mask]

    # Calculate squared distances using NumPy's vectorized operations
    a = positions[i_filtered]
    b = positions[j_filtered]
    distance_squared = np.sum((a - b) ** 2, axis=1)

    # Apply cutoff and get valid indices
    valid = distance_squared < cutoff_sq
    return i_filtered[valid], j_filtered[valid]
    

def generate_matrix(positions, cutoff):
    natoms = positions.shape[0]
    valid_pair = generate_valid_paris(positions, cutoff)
    if valid_pair is None:
        return np.zeros((natoms, natoms), dtype=np.float64)
    i_filtered, j_filtered = valid_pair
    # Create matrix and set symmetric entries
    matrix = np.zeros((natoms, natoms), dtype=np.float64)
    matrix[i_filtered, j_filtered] = -1.0
    matrix[j_filtered, i_filtered] = -1.0

    # Calculate diagonal entries as the count of neighbors
    row_counts = np.sum(matrix < 0, axis=1)
    np.fill_diagonal(matrix, row_counts)

    return matrix

        
def calcu_GNMAnalysis(positions: np.ndarray, cutoff: float = 7,
                      gen_matrix_fn = None, **kwargs):
    """Generate the Kirchhoff matrix of contacts.

    This generates the neighbour matrix by generating a grid of
    near-neighbours and then calculating which are are within
    the cutoff.

    Returns
    -------
        eigenvectors
        eigenvalues
    """
    gen_matrix_fn = gen_matrix_fn or generate_matrix
    matrix = gen_matrix_fn(positions, cutoff, **kwargs)
    try:
        _, w, v = np.linalg.svd(matrix)
    except np.linalg.LinAlgError:
        return put_err(f"SVD with cutoff {cutoff} failed to converge, return None")
    list_map = np.argsort(w)
    return w[list_map[1]], v[list_map[1]]


def generate_close_matrix(positions: np.ndarray, cutoff,
                          atom2residue: np.ndarray, residue_size: np.ndarray,
                          n_residue: int, weights="size"):
    """Generate the Kirchhoff matrix of closeContactGNMAnalysis contacts.

    This generates the neighbour matrix by generating a grid of
    near-neighbours and then calculating which are are within
    the cutoff.

    Returns
    -------
    array
            the resulting Kirchhoff matrix
    """
    # Compute residue sizes
    if weights == 'size':
        inv_sqrt_res_sizes = 1.0 / np.sqrt(residue_size)
    else:
        inv_sqrt_res_sizes = np.ones(n_residue, dtype=np.float64)

    # Generate all atom pairs within cutoff
    # Note: Using previous generate_ordered_pairs function (adjusted for pairs)
    valid_pair = generate_valid_paris(positions, cutoff)
    if valid_pair is None:
        return np.zeros((n_residue, n_residue), dtype=np.float64)
    i_filtered, j_filtered = valid_pair

    # Get valid residue indices
    iresidues = atom2residue[i_filtered]
    jresidues = atom2residue[j_filtered]

    # Compute contact values
    contact = inv_sqrt_res_sizes[iresidues] * inv_sqrt_res_sizes[jresidues]

    # Initialize Kirkhoff matrix
    matrix = np.zeros((n_residue, n_residue), dtype=np.float64)

    # Update symmetric pairs
    matrix[iresidues, jresidues] -= contact
    matrix[jresidues, iresidues] -= contact
    matrix[iresidues, iresidues] += contact
    matrix[jresidues, jresidues] += contact

    # # Update diagonal elements
    # for res in range(n_residue):
    #     diagonal_contacts = contact[(iresidues == res) | (jresidues == res)]
    #     matrix[res, res] += np.sum(diagonal_contacts)
    
    # # Update diagonal elements using bincount
    # # Combine iresidues and jresidues and concatenate the contact twice
    # ire_jre_concat = np.concatenate((iresidues, jresidues))
    # contact_concat = np.concatenate((contact, contact))

    # # Compute the bincount for the combined residues
    # bincounts = np.bincount(ire_jre_concat, weights=contact_concat, minlength=n_residue)

    # # Add the bincounts to the diagonal of the matrix
    # np.fill_diagonal(matrix, bincounts)

    return matrix


def genarate_atom2residue(atoms: AtomGroup):
    """
    return
        - a 1d array where each element is the residue index of the atom
        - a 1d array where each element is the number of atoms in the residue
    """
    return atoms.resindices.copy(), np.array([r.atoms.n_atoms for r in atoms.residues])


def calcu_closeContactGNMAnalysis(positions: np.ndarray, cutoff: float, atom2residue: np.ndarray,
                                  residue_size: np.ndarray, n_residue: int, weights="size"):
    """Generate the Kirchhoff matrix of contacts.

    This generates the neighbour matrix by generating a grid of
    near-neighbours and then calculating which are are within
    the cutoff.

    Returns
    -------
        eigenvectors
        eigenvalues
    """
    return calcu_GNMAnalysis(positions, cutoff, gen_matrix_fn=generate_close_matrix,
                             atom2residue=atom2residue, residue_size=residue_size,
                             n_residue=n_residue, weights=weights)
