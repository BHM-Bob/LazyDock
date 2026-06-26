'''
Date: 2026-06-19
LastEditors: BHM-Bob
Description: GaMD reweighting algorithms for energetic and kinetic analysis.

Implements:
- parse_gamd_log: Parse Amber GaMD log file
- calc_anharmonicity: Assess anharmonicity of boost potential distribution
- reweight_1d: 1D PMF via cumulant expansion to 2nd order (custom)
- reweight_2d: 2D PMF via cumulant expansion to 2nd order (custom)
- reweight_1d_pyreweighting: 1D PMF via PyReweighting (standard tool, subprocess)
- reweight_1d_ce: 1D PMF via cumulant expansion to 3rd order (PyReweighting algorithm, embedded)
- reweight_2d_pyreweighting: 2D PMF via PyReweighting (standard tool, subprocess)
- reweight_vdr: 1D/2D PMF via gamdvdr VDR method (standard tool, subprocess)
- calc_kramers_rate: Kinetic reweighting via Kramers rate theory
- calc_diffusion_coefficient: Estimate apparent diffusion coefficient from CV MSD
- find_barriers: Identify barriers and wells in 1D PMF

Standard tool backends:
- PyReweighting: MiaoLab official reweighting tool (cumulant/Maclaurin/exponential)
  https://github.com/MiaoLab20/pyreweighting
- gamdvdr: Variable Density Reweighting (improved PyReweighting)
  https://github.com/SC-Turner/GaMD_Variable_Density_Reweighting

References:
- Miao et al., J. Chem. Theory Comput. 2014, 10, 10, 4667-4677
- Miao et al., J. Chem. Phys. 2015, 143, 124110
- Turner et al., gamdvdr, https://pypi.org/project/gamdvdr/
'''
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

# Physical constants
KB = 0.001987  # Boltzmann constant in kcal/(mol*K)


def parse_gamd_log(gamd_log_path: Path) -> Optional[pd.DataFrame]:
    """Parse Amber GaMD log file (gamd.log) into a DataFrame.

    gamd.log format (8 columns per line):
    ntwx  total_nstep  Unboosted-Potential-Energy  Unboosted-Dihedral-Energy
    Total-Force-Weight  Dihedral-Force-Weight  Boost-Energy-Potential  Boost-Energy-Dihedral

    Parameters
    ----------
    gamd_log_path : Path
        Path to the gamd.log file.

    Returns
    -------
    pd.DataFrame or None
        DataFrame with 8 columns, or None if file not found or empty.
    """
    if not gamd_log_path.exists():
        return None

    try:
        df = pd.read_csv(
            gamd_log_path, sep=r'\s+', comment='#',
            names=['ntwx', 'total_nstep', 'unboosted_potential',
                   'unboosted_dihedral', 'total_force_weight',
                   'dihedral_force_weight', 'boost_energy_potential',
                   'boost_energy_dihedral'])
        df = df.dropna()
        if len(df) == 0:
            return None
        return df
    except Exception:
        return None


def calc_anharmonicity(boost_potential: np.ndarray,
                       boost_dihedral: np.ndarray,
                       temp: float = 300.0) -> Dict:
    """Assess anharmonicity of the boost potential distribution.

    Computes the anharmonicity index gamma (entropy-based, per Miao et al.,
    JCTC 2015), as well as skewness and excess kurtosis, for the total boost
    potential delta_V = delta_V_P + delta_V_D and each component separately.

    The anharmonicity index gamma is defined as the entropy difference between
    a Gaussian with the same variance and the actual distribution:

        gamma = S_max - S_dV
              = 0.5 * ln(2*pi*e*sigma_dV^2) + integral(p(dV)*ln(p(dV)) ddV)

    where dV is dimensionless (divided by kBT), sigma_dV is the std of the
    dimensionless dV, and S_max is the maximum (Gaussian) entropy.

    When gamma = 0, dV follows exact Gaussian distribution and cumulant
    expansion to 2nd order recovers the original free energy exactly.
    As gamma increases, the dV distribution becomes less harmonic and
    reweighting accuracy degrades.

    Quality criteria (based on typical values in GaMD literature):
    - gamma < 0.01: Excellent (e.g., alanine dipeptide ~0.002)
    - 0.01 <= gamma < 0.05: Acceptable (e.g., M3 receptor ~0.007)
    - gamma >= 0.05: Poor (reweighting may be unreliable)

    Note: The previous implementation used skewness as gamma, which is NOT
    the same as the entropy-based anharmonicity defined in the original paper.
    Skewness is still reported as a supplementary statistic.

    Parameters
    ----------
    boost_potential : np.ndarray
        Boost energy for total potential (delta_V_P), shape (N,).
    boost_dihedral : np.ndarray
        Boost energy for dihedral potential (delta_V_D), shape (N,).
    temp : float
        Temperature in Kelvin (used to make dV dimensionless via kBT).

    Returns
    -------
    dict
        Keys: 'total', 'potential', 'dihedral', each containing
        {'mean', 'std', 'anharmonicity', 'skewness', 'excess_kurtosis',
         'quality'}.
    """
    delta_V = boost_potential + boost_dihedral
    kBT = KB * temp  # kcal/mol

    def _stats(data):
        mu = np.mean(data)
        sigma = np.std(data)

        # Skewness and excess kurtosis (supplementary statistics)
        if sigma > 0:
            skewness = float(np.mean(((data - mu) / sigma) ** 3))
            excess_kurtosis = float(np.mean(((data - mu) / sigma) ** 4) - 3)
        else:
            skewness = 0.0
            excess_kurtosis = 0.0

        # Anharmonicity index gamma = S_max - S_dV (entropy-based)
        # Following Miao et al., JCTC 2015, Eq. in Appendix B
        if sigma > 0 and len(data) > 1:
            # Make dV dimensionless by dividing by kBT
            dV_dimless = data / kBT
            sigma_dimless = np.std(dV_dimless)

            # S_max = 0.5 * ln(2*pi*e*sigma^2)  (Gaussian entropy)
            S_max = 0.5 * np.log(2 * np.pi * np.e * sigma_dimless ** 2)

            # S_dV = -integral(p(x)*ln(p(x)) dx)  (actual distribution entropy)
            # Estimate via histogram binning
            n_bins = min(200, max(50, len(data) // 100))
            hist, bin_edges = np.histogram(dV_dimless, bins=n_bins, density=True)
            bin_width = bin_edges[1] - bin_edges[0]
            # Only consider bins with non-zero probability
            mask = hist > 0
            S_dV = -np.sum(hist[mask] * np.log(hist[mask]) * bin_width)

            gamma = S_max - S_dV
            # Numerical safety: gamma should be >= 0 for any distribution
            gamma = max(0.0, float(gamma))
        else:
            gamma = 0.0

        # Quality assessment based on entropy-based gamma
        if gamma < 0.01:
            quality = "Excellent"
        elif gamma < 0.05:
            quality = "Acceptable"
        else:
            quality = "Poor"

        return {'mean': mu, 'std': sigma, 'anharmonicity': gamma,
                'skewness': skewness, 'excess_kurtosis': excess_kurtosis,
                'quality': quality}

    return {
        'total': _stats(delta_V),
        'potential': _stats(boost_potential),
        'dihedral': _stats(boost_dihedral),
    }


def reweight_1d(cv: np.ndarray, delta_V: np.ndarray,
                temp: float = 300.0, cutoff: int = 10,
                disc: Optional[float] = None) -> Optional[Tuple]:
    """1D PMF via cumulant expansion to 2nd order.

    Computes the reweighted PMF along a single collective variable using
    the cumulant expansion method (Miao et al., JCTC 2014).

    The reweighted free energy for each bin is:
        F(CV) = -kBT * [ln(N(CV)) - beta*C1(CV) - (beta^2/2)*C2(CV)]

    where:
        C1(CV) = <delta_V>_CV  (1st cumulant = mean of delta_V in that bin)
        C2(CV) = <delta_V^2>_CV - <delta_V>^2_CV  (2nd cumulant = variance)
        beta = 1/(kBT)

    Parameters
    ----------
    cv : np.ndarray
        Collective variable values, shape (N,).
    delta_V : np.ndarray
        Total boost potential delta_V = delta_V_P + delta_V_D, shape (N,).
    temp : float
        Temperature in Kelvin.
    cutoff : int
        Minimum number of frames per bin for reweighting.
    disc : float or None
        Bin size for CV. Auto-determined if None.

    Returns
    -------
    tuple or None
        (bin_centers, pmf_values, counts, c1_values, c2_values) or None on error.
    """
    kBT = KB * temp
    beta = 1.0 / kBT

    # Determine binning
    cv_range = np.max(cv) - np.min(cv)
    if disc is None:
        disc = cv_range / (int(np.sqrt(len(cv))) + 1)
    if disc <= 0 or cv_range <= 0: # type: ignore
        return None

    n_bins = max(10, int(cv_range / disc))
    bin_edges = np.linspace(np.min(cv), np.max(cv), n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_indices = np.digitize(cv, bin_edges) - 1

    # Clip to valid range
    valid_mask = (bin_indices >= 0) & (bin_indices < n_bins)
    bin_indices_clipped = np.clip(bin_indices, 0, n_bins - 1)

    # Compute cumulants per bin
    counts = np.zeros(n_bins, dtype=int)
    c1_values = np.full(n_bins, np.nan)
    c2_values = np.full(n_bins, np.nan)

    for i in range(n_bins):
        mask = valid_mask & (bin_indices_clipped == i)
        n_in_bin = np.sum(mask)
        counts[i] = n_in_bin
        if n_in_bin >= cutoff:
            dV_in_bin = delta_V[mask]
            c1_values[i] = np.mean(dV_in_bin)
            c2_values[i] = np.var(dV_in_bin)

    # Compute PMF via cumulant expansion
    valid_counts = counts >= cutoff
    pmf_values = np.full(n_bins, np.nan)

    if np.sum(valid_counts) == 0:
        return None

    ln_counts = np.log(counts[valid_counts].astype(float))
    c1_valid = c1_values[valid_counts]
    c2_valid = c2_values[valid_counts]

    raw_pmf = -kBT * (ln_counts - beta * c1_valid - 0.5 * beta**2 * c2_valid)

    # Shift so minimum is 0
    pmf_values[valid_counts] = raw_pmf - np.nanmin(raw_pmf)

    return bin_centers, pmf_values, counts, c1_values, c2_values


def reweight_2d(cv1: np.ndarray, cv2: np.ndarray, delta_V: np.ndarray,
                temp: float = 300.0, cutoff: int = 10,
                disc: Optional[float] = None) -> Optional[Tuple]:
    """2D PMF via cumulant expansion to 2nd order.

    Parameters
    ----------
    cv1, cv2 : np.ndarray
        Two collective variables, shape (N,).
    delta_V : np.ndarray
        Total boost potential, shape (N,).
    temp : float
        Temperature in Kelvin.
    cutoff : int
        Minimum number of frames per bin.
    disc : float or None
        Bin size for both CVs. Auto-determined if None.

    Returns
    -------
    tuple or None
        (x_centers, y_centers, pmf_grid, counts_grid) or None on error.
    """
    kBT = KB * temp
    beta = 1.0 / kBT

    # Determine binning
    if disc is None:
        disc1 = (np.max(cv1) - np.min(cv1)) / (int(np.sqrt(len(cv1))) + 1)
        disc2 = (np.max(cv2) - np.min(cv2)) / (int(np.sqrt(len(cv2))) + 1)
    else:
        disc1 = disc2 = disc

    n_bins_x = max(10, int((np.max(cv1) - np.min(cv1)) / disc1))
    n_bins_y = max(10, int((np.max(cv2) - np.min(cv2)) / disc2))

    x_edges = np.linspace(np.min(cv1), np.max(cv1), n_bins_x + 1)
    y_edges = np.linspace(np.min(cv2), np.max(cv2), n_bins_y + 1)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    # Digitize
    ix = np.clip(np.digitize(cv1, x_edges) - 1, 0, n_bins_x - 1)
    iy = np.clip(np.digitize(cv2, y_edges) - 1, 0, n_bins_y - 1)

    # Compute cumulants per 2D bin
    counts_grid = np.zeros((n_bins_x, n_bins_y), dtype=int)
    c1_grid = np.full((n_bins_x, n_bins_y), np.nan)
    c2_grid = np.full((n_bins_x, n_bins_y), np.nan)

    for i in range(n_bins_x):
        for j in range(n_bins_y):
            mask = (ix == i) & (iy == j)
            n_in_bin = np.sum(mask)
            counts_grid[i, j] = n_in_bin
            if n_in_bin >= cutoff:
                dV_in_bin = delta_V[mask]
                c1_grid[i, j] = np.mean(dV_in_bin)
                c2_grid[i, j] = np.var(dV_in_bin)

    # Compute PMF
    valid = counts_grid >= cutoff
    pmf_grid = np.full((n_bins_x, n_bins_y), np.nan)

    if np.sum(valid) == 0:
        return None

    ln_counts = np.log(counts_grid[valid].astype(float))
    raw_pmf = -kBT * (ln_counts - beta * c1_grid[valid] - 0.5 * beta**2 * c2_grid[valid])
    pmf_grid[valid] = raw_pmf - np.nanmin(raw_pmf)

    return x_centers, y_centers, pmf_grid, counts_grid


def calc_diffusion_coefficient(cv: np.ndarray, max_lag: Optional[int] = None) -> float:
    """Estimate apparent diffusion coefficient from CV mean-squared displacement.

    Uses the initial slope of the MSD vs lag plot:
        D = slope / 2  (for 1D diffusion)

    Parameters
    ----------
    cv : np.ndarray
        Collective variable time series, shape (N,).
    max_lag : int or None
        Maximum lag for MSD calculation. Default: min(1000, N//10).

    Returns
    -------
    float
        Apparent diffusion coefficient D. Returns 0.0 if insufficient data.
    """
    if len(cv) <= 100:
        return 0.0

    if max_lag is None:
        max_lag = min(1000, len(cv) // 10)

    msd_values = np.array([np.mean((cv[lag:] - cv[:-lag])**2) for lag in range(1, max_lag)])
    lags = np.arange(1, max_lag)

    # Linear fit to first 10% of lags
    n_fit = max(3, max_lag // 10)
    coeffs = np.polyfit(lags[:n_fit], msd_values[:n_fit], 1)
    return coeffs[0] / 2.0


def find_barriers(bin_centers: np.ndarray, pmf_values: np.ndarray,
                  smooth_size: Optional[int] = None) -> Optional[Dict]:
    """Identify barriers and wells in a 1D PMF profile.

    Parameters
    ----------
    bin_centers : np.ndarray
        CV bin center values.
    pmf_values : np.ndarray
        PMF values (may contain NaN for under-sampled bins).
    smooth_size : int or None
        Window size for smoothing PMF before barrier detection.
        Default: min(7, len(pmf_valid)//3).

    Returns
    -------
    dict or None
        Dictionary with keys:
        - 'well1_idx', 'well1_cv', 'well1_pmf': deepest well
        - 'well2_idx', 'well2_cv', 'well2_pmf': second deepest well
        - 'barrier_idx', 'barrier_cv', 'barrier_pmf': barrier between wells
        - 'delta_G_well1': barrier height from well1 (kcal/mol)
        - 'delta_G_well2': barrier height from well2 (kcal/mol)
        - 'curvature_well': d2F/dx2 at well1
        - 'curvature_barrier': d2F/dx2 at barrier
        Returns None if no barrier found.
    """
    from scipy.signal import argrelextrema

    valid = ~np.isnan(pmf_values)
    if np.sum(valid) < 5:
        return None

    bc_valid = bin_centers[valid]
    pmf_valid = pmf_values[valid]

    # Smooth PMF for barrier detection
    if smooth_size is None:
        smooth_size = min(7, len(pmf_valid) // 3) if len(pmf_valid) > 7 else 1

    if smooth_size > 1 and len(pmf_valid) > smooth_size:
        from scipy.ndimage import uniform_filter1d
        pmf_smooth = uniform_filter1d(pmf_valid, size=smooth_size)
    else:
        pmf_smooth = pmf_valid

    order = max(2, len(pmf_smooth) // 20)
    min_indices = argrelextrema(pmf_smooth, np.less, order=order)[0]
    max_indices = argrelextrema(pmf_smooth, np.greater, order=order)[0]

    # If no maxima found but we have minima, find barrier between minima directly
    if len(min_indices) >= 2 and len(max_indices) == 0:
        sorted_min_idx = min_indices[np.argsort(pmf_smooth[min_indices])]
        well1_idx = sorted_min_idx[0]
        well2_idx = sorted_min_idx[1]
        lo, hi = min(well1_idx, well2_idx), max(well1_idx, well2_idx)
        barrier_idx = lo + np.argmax(pmf_smooth[lo:hi+1])
        max_indices = np.array([barrier_idx])

    if len(min_indices) < 2 or len(max_indices) < 1:
        return None

    # Find the two deepest minima
    sorted_min_idx = min_indices[np.argsort(pmf_smooth[min_indices])]
    well1_idx = sorted_min_idx[0]
    well2_idx = sorted_min_idx[1] if len(sorted_min_idx) > 1 else sorted_min_idx[0]

    # Find barrier between well1 and well2
    lo, hi = min(well1_idx, well2_idx), max(well1_idx, well2_idx)
    barrier_candidates = max_indices[(max_indices > lo) & (max_indices < hi)]

    if len(barrier_candidates) == 0:
        return None

    barrier_idx = barrier_candidates[np.argmax(pmf_smooth[barrier_candidates])]

    # Compute curvatures (second derivative)
    dx = bc_valid[1] - bc_valid[0] if len(bc_valid) > 1 else 1.0

    def _curvature(idx, arr):
        if 0 < idx < len(arr) - 1:
            return (arr[idx+1] - 2*arr[idx] + arr[idx-1]) / (dx**2)
        return 0.0

    curv_well = _curvature(well1_idx, pmf_valid)
    curv_barrier = _curvature(barrier_idx, pmf_valid)

    delta_G = pmf_valid[barrier_idx] - pmf_valid[well1_idx]
    delta_G2 = pmf_valid[barrier_idx] - pmf_valid[well2_idx]

    return {
        'well1_idx': well1_idx, 'well1_cv': float(bc_valid[well1_idx]),
        'well1_pmf': float(pmf_valid[well1_idx]),
        'well2_idx': well2_idx, 'well2_cv': float(bc_valid[well2_idx]),
        'well2_pmf': float(pmf_valid[well2_idx]),
        'barrier_idx': barrier_idx, 'barrier_cv': float(bc_valid[barrier_idx]),
        'barrier_pmf': float(pmf_valid[barrier_idx]),
        'delta_G_well1': float(delta_G),
        'delta_G_well2': float(delta_G2),
        'curvature_well': float(curv_well),
        'curvature_barrier': float(curv_barrier),
    }


def calc_kramers_rate(barrier_info: Dict, D_apparent: float,
                      temp: float = 300.0) -> Dict:
    """Estimate Kramers rate from barrier info and diffusion coefficient.

    k = (omega_min * omega_barrier * D) / (2*pi*kBT) * exp(-delta_G/kBT)

    Parameters
    ----------
    barrier_info : dict
        Output from find_barriers(), must contain 'curvature_well',
        'curvature_barrier', 'delta_G_well1'.
    D_apparent : float
        Apparent diffusion coefficient from calc_diffusion_coefficient().
    temp : float
        Temperature in Kelvin.

    Returns
    -------
    dict
        Keys: 'omega_min', 'omega_barrier', 'k_gamd', 'acceleration_factor',
        'k_cMD_estimate'. Values are NaN if estimation not possible.
    """
    kBT = KB * temp
    beta = 1.0 / kBT # type: ignore

    curv_well = barrier_info['curvature_well']
    curv_barrier = barrier_info['curvature_barrier']
    delta_G = barrier_info['delta_G_well1'] # type: ignore

    result = {
        'omega_min': np.nan,
        'omega_barrier': np.nan,
        'k_gamd': np.nan,
        'acceleration_factor': np.nan,
        'k_cMD_estimate': np.nan,
    }

    # Well curvature should be positive (concave up), barrier negative (concave down)
    if curv_well <= 0 or curv_barrier >= 0 or D_apparent <= 0:
        return result

    omega_min = np.sqrt(abs(curv_well))
    omega_barrier = np.sqrt(abs(curv_barrier))

    # Kramers rate
    k_gamd = (omega_min * omega_barrier * D_apparent) / (2 * np.pi * kBT)

    result['omega_min'] = float(omega_min)
    result['omega_barrier'] = float(omega_barrier)
    result['k_gamd'] = float(k_gamd)

    return result


def calc_acceleration_factor(delta_V: np.ndarray, cv: np.ndarray,
                             barrier_cv: float, well_cv: float,
                             temp: float = 300.0) -> float:
    """Estimate GaMD acceleration factor from boost potential at barrier vs well.

    The acceleration factor is approximately:
        accel = exp(beta * (delta_V_at_barrier - delta_V_at_well))

    Parameters
    ----------
    delta_V : np.ndarray
        Total boost potential, shape (N,).
    cv : np.ndarray
        CV values, shape (N,).
    barrier_cv : float
        CV value at the barrier top.
    well_cv : float
        CV value at the well minimum.
    temp : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Acceleration factor.
    """
    kBT = KB * temp
    beta = 1.0 / kBT

    dV_at_barrier = delta_V[np.argmin(np.abs(cv - barrier_cv))]
    dV_at_well = delta_V[np.argmin(np.abs(cv - well_cv))]
    delta_V_effect = dV_at_barrier - dV_at_well

    if abs(delta_V_effect) > 0:
        return float(np.exp(beta * delta_V_effect))
    return 1.0


# ============================================================================
# Standard Tool Interfaces: PyReweighting and gamdvdr
# ============================================================================

def _check_command_available(cmd: str) -> bool:
    """Check if a command-line tool is available in PATH."""
    return shutil.which(cmd) is not None


def _prepare_pyreweighting_weights(delta_V: np.ndarray, temp: float = 300.0,
                                   workdir: Optional[str] = None) -> str:
    """Prepare weights.dat file for PyReweighting.

    Format: column 1 = dV/kBT, column 2 = timestep, column 3 = dV in kcal/mol
    Reference: https://github.com/MiaoLab20/pyreweighting

    Parameters
    ----------
    delta_V : np.ndarray
        Total boost potential in kcal/mol, shape (N,).
    temp : float
        Temperature in Kelvin.
    workdir : str or None
        Directory to write the file. Uses tempdir if None.

    Returns
    -------
    str
        Path to the weights.dat file.
    """
    kBT = KB * temp
    n_frames = len(delta_V)
    dV_kBT = delta_V / kBT
    timesteps = np.arange(1, n_frames + 1)
    data = np.column_stack([dV_kBT, timesteps, delta_V])
    filepath = os.path.join(workdir or tempfile.mkdtemp(), 'weights.dat')
    np.savetxt(filepath, data, fmt='%.6f %d %.6f')
    return filepath


def _prepare_pyreweighting_cv(cv: np.ndarray, workdir: Optional[str] = None,
                              filename: str = 'cv.dat') -> str:
    """Prepare CV data file for PyReweighting.

    Format: one column of CV values per line.

    Parameters
    ----------
    cv : np.ndarray
        Collective variable values, shape (N,).
    workdir : str or None
        Directory to write the file.
    filename : str
        Output filename.

    Returns
    -------
    str
        Path to the CV data file.
    """
    filepath = os.path.join(workdir or tempfile.mkdtemp(), filename)
    np.savetxt(filepath, cv, fmt='%.6f')
    return filepath


def _parse_xvg_pmf(xvg_path: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """Parse PyReweighting output .xvg PMF file.

    Format: two columns (RC, PMF in kcal/mol), with comment lines starting
    with '#' or '@'.

    Parameters
    ----------
    xvg_path : str
        Path to the .xvg file.

    Returns
    -------
    tuple or None
        (bin_centers, pmf_values) or None if file not found.
    """
    if not os.path.exists(xvg_path):
        return None
    data = []
    with open(xvg_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    data.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    continue
    if not data:
        return None
    arr = np.array(data)
    return arr[:, 0], arr[:, 1]


def reweight_1d_pyreweighting(cv: np.ndarray, delta_V: np.ndarray,
                              temp: float = 300.0, cutoff: int = 10,
                              disc: float = 6.0, emax: float = 20.0,
                              job: str = 'amdweight_CE',
                              pyreweighting_script: Optional[str] = None,
                              workdir: Optional[str] = None,
                              cleanup: bool = True) -> Optional[Dict]:
    """1D PMF reweighting via PyReweighting command-line tool.

    Calls PyReweighting-1D.py as a subprocess. This is the standard reference
    implementation from MiaoLab for GaMD reweighting.

    Install: clone from https://github.com/MiaoLab20/pyreweighting

    Parameters
    ----------
    cv : np.ndarray
        Collective variable values, shape (N,).
    delta_V : np.ndarray
        Total boost potential in kcal/mol, shape (N,).
    temp : float
        Temperature in Kelvin.
    cutoff : int
        Minimum number of frames per bin.
    disc : float
        Bin size (discretization) for the CV.
    emax : float
        Maximum free energy value for unsampled bins (kcal/mol).
    job : str
        Reweighting method: 'amdweight_CE' (cumulant expansion),
        'amdweight_MC' (Maclaurin series), 'amdweight' (exponential average).
    pyreweighting_script : str or None
        Path to PyReweighting-1D.py. Auto-detected if None.
    workdir : str or None
        Working directory for input/output files. Uses tempdir if None.
    cleanup : bool
        Whether to remove temporary files after completion.

    Returns
    -------
    dict or None
        Dictionary with keys depending on job type:
        - 'pmf_c2': (bin_centers, pmf_values) for CE2 (recommended)
        - 'pmf_c1': (bin_centers, pmf_values) for CE1
        - 'pmf_c3': (bin_centers, pmf_values) for CE3
        - 'pmf': (bin_centers, pmf_values) for MC/exponential
        Returns None if PyReweighting is not available or fails.
    """
    # Find PyReweighting-1D.py
    script = pyreweighting_script
    if script is None:
        script = shutil.which('PyReweighting-1D.py')
    if script is None:
        # Try common locations
        for candidate in ['PyReweighting-1D.py']:
            if os.path.exists(candidate):
                script = candidate
                break
    if script is None:
        return None

    # Create workdir
    own_workdir = workdir is None
    if own_workdir:
        workdir = tempfile.mkdtemp(prefix='pyreweighting_')
    try:
        # Prepare input files
        weights_file = _prepare_pyreweighting_weights(delta_V, temp, workdir)
        cv_file = _prepare_pyreweighting_cv(cv, workdir, 'cv.dat')

        # Build command
        cmd = [
            'python', script,
            '-input', cv_file,
            '-T', str(temp),
            '-cutoff', str(cutoff),
            '-disc', str(disc),
            '-Emax', str(emax),
            '-job', job,
            '-weight', weights_file,
        ]

        # Run PyReweighting
        result = subprocess.run(
            cmd, capture_output=True, text=True, cwd=workdir, timeout=300)

        if result.returncode != 0:
            return None

        # Parse output files
        output = {}
        if job == 'amdweight_CE':
            for order, suffix in [(1, 'c1'), (2, 'c2'), (3, 'c3')]: # type: ignore
                xvg_path = os.path.join(workdir, f'pmf-{suffix}-cv.dat.xvg')
                parsed = _parse_xvg_pmf(xvg_path)
                if parsed is not None:
                    output[f'pmf_{suffix}'] = parsed
        elif job in ('amdweight_MC', 'amdweight'):
            xvg_path = os.path.join(workdir, 'pmf-cv.dat.xvg')
            parsed = _parse_xvg_pmf(xvg_path)
            if parsed is not None:
                output['pmf'] = parsed

        return output if output else None

    finally:
        if cleanup and own_workdir:
            shutil.rmtree(workdir, ignore_errors=True)


def reweight_1d_ce(cv: np.ndarray, delta_V: np.ndarray,
                   temp: float = 300.0, cutoff: int = 10,
                   disc: Optional[float] = None,
                   emax: float = 8.0) -> Optional[Dict]:
    """1D PMF via cumulant expansion to 3rd order (PyReweighting algorithm).

    This is an embedded implementation of the PyReweighting reweight_CE
    algorithm (Miao et al., JCTC 2014), supporting cumulant expansion up to
    3rd order. It produces the same results as running
    ``PyReweighting-1D.py -job amdweight_CE`` without requiring the external
    script.

    The reweighted PMF at each order:
        PMF_C1 = -kBT * [ln(N) - beta*C1]
        PMF_C2 = -kBT * [ln(N) - beta*C1 - (beta^2/2)*C2]  (recommended)
        PMF_C3 = -kBT * [ln(N) - beta*C1 - (beta^2/2)*C2 - (beta^3/6)*C3]

    where:
        C1 = <dV>_bin  (1st cumulant)
        C2 = <dV^2>_bin - <dV>^2_bin  (2nd cumulant = variance)
        C3 = <dV^3>_bin - 3*<dV^2>*<dV> + 2*<dV>^3  (3rd cumulant)

    Parameters
    ----------
    cv : np.ndarray
        Collective variable values, shape (N,).
    delta_V : np.ndarray
        Total boost potential in kcal/mol, shape (N,).
    temp : float
        Temperature in Kelvin.
    cutoff : int
        Minimum number of frames per bin.
    disc : float or None
        Bin size for CV. Auto-determined if None (same as PyReweighting default).
    emax : float
        Maximum free energy value for unsampled bins (kcal/mol).

    Returns
    -------
    dict or None
        Dictionary with keys:
        - 'bin_centers': np.ndarray of CV bin center values
        - 'counts': np.ndarray of frame counts per bin
        - 'c1': np.ndarray of 1st cumulant per bin (in kBT units)
        - 'c2': np.ndarray of 2nd cumulant per bin (in kBT units)
        - 'c3': np.ndarray of 3rd cumulant per bin (in kBT units)
        - 'pmf_c1': PMF with 1st order correction (kcal/mol)
        - 'pmf_c2': PMF with 2nd order correction (kcal/mol, recommended)
        - 'pmf_c3': PMF with 3rd order correction (kcal/mol)
        Returns None on error.
    """
    kBT = KB * temp
    beta = 1.0 / kBT

    # Determine binning (same logic as PyReweighting)
    if disc is None:
        disc = 6.0  # PyReweighting default
    cv_min, cv_max = np.min(cv), np.max(cv)
    bins = np.arange(
        disc * (int(cv_min / disc) - 1),
        disc * (int(cv_max / disc) + 1) + disc,
        disc
    )
    n_bins = len(bins) - 1
    if n_bins < 1:
        return None
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    # Assign frames to bins
    counts = np.zeros(n_bins, dtype=int)
    dV_sums = np.zeros(n_bins)
    dV2_sums = np.zeros(n_bins)
    dV3_sums = np.zeros(n_bins)

    for i in range(len(cv)):
        j = int((cv[i] - bins[0]) / disc)
        if j >= n_bins:
            j = n_bins - 1
        if j < 0:
            j = 0
        counts[j] += 1
        dV_sums[j] += delta_V[i]
        dV2_sums[j] += delta_V[i] ** 2
        dV3_sums[j] += delta_V[i] ** 3

    # Compute cumulants per bin
    c1 = np.zeros(n_bins)
    c2 = np.zeros(n_bins)
    c3 = np.zeros(n_bins)

    for j in range(n_bins):
        if counts[j] >= cutoff:
            n = counts[j]
            dV_avg = dV_sums[j] / n
            dV2_avg = dV2_sums[j] / n
            dV3_avg = dV3_sums[j] / n
            dV_std = np.sqrt(dV2_avg - dV_avg ** 2)
            c1[j] = beta * dV_avg
            c2[j] = 0.5 * beta ** 2 * dV_std ** 2
            c3[j] = (1.0 / 6.0) * beta ** 3 * (
                dV3_avg - 3.0 * dV2_avg * dV_avg + 2.0 * dV_avg ** 3
            )

    # Compute PMF at each cumulant order
    # Base: -kBT * ln(counts)
    hist = counts.astype(float).copy()
    hist[hist == 0] = 1e-15  # avoid log(0)
    base_pmf = kBT * np.log(hist)
    base_pmf = np.max(base_pmf) - base_pmf  # normalize so min = 0

    # PMF with cumulant corrections
    c1_kcal = -c1 / beta  # convert back from kBT to kcal/mol
    c2_kcal = -c2 / beta
    c3_kcal = -c3 / beta

    pmf_c1 = base_pmf + c1_kcal
    pmf_c2 = pmf_c1 + c2_kcal
    pmf_c3 = pmf_c2 + c3_kcal

    # Normalize each PMF so minimum is 0, cap at emax
    def _normalize_pmf(pmf, emax_val):
        valid = counts >= cutoff
        if not np.any(valid):
            return pmf
        pmf_min = np.min(pmf[valid])
        pmf = pmf - pmf_min
        pmf[~valid] = emax_val
        pmf[pmf > emax_val] = emax_val
        return pmf

    pmf_c1 = _normalize_pmf(pmf_c1, emax)
    pmf_c2 = _normalize_pmf(pmf_c2, emax)
    pmf_c3 = _normalize_pmf(pmf_c3, emax)

    return {
        'bin_centers': bin_centers,
        'counts': counts,
        'c1': c1,
        'c2': c2,
        'c3': c3,
        'pmf_c1': pmf_c1,
        'pmf_c2': pmf_c2,
        'pmf_c3': pmf_c3,
    }


def reweight_2d_pyreweighting(cv1: np.ndarray, cv2: np.ndarray,
                              delta_V: np.ndarray,
                              temp: float = 300.0, cutoff: int = 10,
                              disc_x: float = 6.0, disc_y: float = 6.0,
                              emax: float = 20.0,
                              job: str = 'amdweight_CE',
                              pyreweighting_script: Optional[str] = None,
                              workdir: Optional[str] = None,
                              cleanup: bool = True) -> Optional[Dict]:
    """2D PMF reweighting via PyReweighting command-line tool.

    Calls PyReweighting-2D.py as a subprocess.

    Install: clone from https://github.com/MiaoLab20/pyreweighting

    Parameters
    ----------
    cv1, cv2 : np.ndarray
        Two collective variables, shape (N,).
    delta_V : np.ndarray
        Total boost potential in kcal/mol, shape (N,).
    temp : float
        Temperature in Kelvin.
    cutoff : int
        Minimum number of frames per bin.
    disc_x, disc_y : float
        Bin size for each CV dimension.
    emax : float
        Maximum free energy value for unsampled bins.
    job : str
        Reweighting method: 'amdweight_CE', 'amdweight_MC', or 'amdweight'.
    pyreweighting_script : str or None
        Path to PyReweighting-2D.py. Auto-detected if None.
    workdir : str or None
        Working directory for input/output files.
    cleanup : bool
        Whether to remove temporary files after completion.

    Returns
    -------
    dict or None
        Dictionary with keys depending on job type (similar to reweight_1d_pyreweighting).
        Returns None if PyReweighting is not available or fails.
    """
    script = pyreweighting_script
    if script is None:
        script = shutil.which('PyReweighting-2D.py')
    if script is None:
        return None

    own_workdir = workdir is None
    if own_workdir:
        workdir = tempfile.mkdtemp(prefix='pyreweighting_2d_')
    try:
        # Prepare input files
        weights_file = _prepare_pyreweighting_weights(delta_V, temp, workdir)
        cv_data = np.column_stack([cv1, cv2])
        cv_file = os.path.join(workdir, 'cv2d.dat')
        np.savetxt(cv_file, cv_data, fmt='%.6f %.6f')

        cmd = [
            'python', script,
            '-input', cv_file,
            '-T', str(temp),
            '-cutoff', str(cutoff),
            '-discX', str(disc_x),
            '-discY', str(disc_y),
            '-Emax', str(emax),
            '-job', job,
            '-weight', weights_file,
        ]

        result = subprocess.run(
            cmd, capture_output=True, text=True, cwd=workdir, timeout=600)

        if result.returncode != 0:
            return None

        output = {}
        if job == 'amdweight_CE':
            for order, suffix in [(1, 'c1'), (2, 'c2'), (3, 'c3')]: # type: ignore
                xvg_path = os.path.join(workdir, f'pmf-{suffix}-cv2d.dat.xvg')
                parsed = _parse_xvg_pmf(xvg_path)
                if parsed is not None:
                    output[f'pmf_{suffix}'] = parsed
        elif job in ('amdweight_MC', 'amdweight'):
            xvg_path = os.path.join(workdir, 'pmf-cv2d.dat.xvg')
            parsed = _parse_xvg_pmf(xvg_path)
            if parsed is not None:
                output['pmf'] = parsed

        return output if output else None

    finally:
        if cleanup and own_workdir:
            shutil.rmtree(workdir, ignore_errors=True)


def reweight_vdr(gamd_log_path: Union[str, Path],
                 cv_data_path: Union[str, Path],
                 conv_points: int = 9500,
                 mode: str = 'single',
                 emax: float = 8.0,
                 pbc: bool = False,
                 cores: int = 1,
                 output_dir: Optional[str] = None,
                 cleanup: bool = True) -> Optional[Dict]:
    """1D/2D PMF reweighting via gamdvdr VDR method (subprocess).

    gamdvdr (Variable Density Reweighting) is an improved reweighting method
    over PyReweighting, optimizing for large systems with variable density
    sampling.

    Install: ``pip install gamdvdr``
    Reference: https://github.com/SC-Turner/GaMD_Variable_Density_Reweighting

    Parameters
    ----------
    gamd_log_path : str or Path
        Path to the gamd.log file from GaMD simulation.
    cv_data_path : str or Path
        Path to CV data file. Format: whitespace-delimited with columns
        CV1 [CV2] frame_number (1D: 2 cols, 2D: 3 cols).
    conv_points : int
        VDR cut-off value for segmentation. Controls the density threshold.
    mode : str
        'single' for single PMF, 'convergence' for convergence analysis.
    emax : float
        Maximum free energy for unsampled regions (kcal/mol).
    pbc : bool
        Whether to handle periodic boundary conditions (e.g., dihedral angles).
    cores : int
        Number of CPU cores for parallel VDR computation.
    output_dir : str or None
        Output directory. Uses tempdir if None.
    cleanup : bool
        Whether to remove temporary files after completion.

    Returns
    -------
    dict or None
        Dictionary with keys:
        - 'pmf_path': path to the output PMF file
        - 'output_dir': path to the output directory (if not cleaned up)
        Returns None if gamdvdr is not available or fails.
    """
    if not _check_command_available('VDR'):
        return None

    gamd_log_path = Path(gamd_log_path)
    cv_data_path = Path(cv_data_path)

    if not gamd_log_path.exists() or not cv_data_path.exists():
        return None

    own_output = output_dir is None
    if own_output:
        output_dir = tempfile.mkdtemp(prefix='vdr_output_')
    try:
        cmd = [
            'VDR',
            '--gamd', str(gamd_log_path),
            '--data', str(cv_data_path),
            '--mode', mode,
            '--conv_points', str(conv_points),
            '--emax', str(emax),
            '--pbc', str(pbc),
            '--cores', str(cores),
            '--output', output_dir,
        ]

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            return None

        # Find output PMF files
        output = {'output_dir': output_dir}
        for f in os.listdir(output_dir):
            if f.endswith('.png') and 'PMF' in f:
                output['pmf_path'] = os.path.join(output_dir, f)
                break

        return output if output else None

    finally:
        if cleanup and own_output:
            shutil.rmtree(output_dir, ignore_errors=True)


def vdr_param(frames: int, anharm_tol: float = 0.01,
              std_error: float = 0.02, clusters: int = 100) -> Optional[float]:
    """Calculate optimal GaMD standard deviation limits via VDR_param.

    This helps determine the best sigma0P/sigma0D parameters for GaMD
    simulations before running them.

    Install: ``pip install gamdvdr``

    Parameters
    ----------
    frames : int
        Number of frames to be saved in the GaMD trajectory.
    anharm_tol : float
        Anharmonicity tolerance (default 0.01).
    std_error : float
        Standard error tolerance in kcal/mol (default 0.02).
    clusters : int
        Number of local clusters to generate (default 100).

    Returns
    -------
    float or None
        Maximum recommended sum of sigma0P + sigma0D, or None if VDR_param
        is not available.
    """
    if not _check_command_available('VDR_param'):
        return None

    cmd = [
        'VDR_param',
        '--frames', str(frames),
    ]
    # VDR_param doesn't support all optional params via CLI in all versions
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)

    if result.returncode != 0:
        return None

    # Parse output for the sigma limit value
    for line in result.stdout.strip().split('\n'):
        line = line.strip()
        try:
            return float(line)
        except ValueError:
            continue
    return None
