#!/usr/bin/env python3
"""
Amber Trajectory Analysis Script
Similar to ana_gmx.py but for Amber MD simulations

Supports:
- trajconv: Trajectory centering and PBC handling
- simple: Basic analysis (RMSD, RMSF, Rg, hydrogen bonds, SASA, PCA)
"""

import argparse
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web import TaskPool
from MDAnalysis import transformations as trans
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.rdf import InterRDF
from tqdm import tqdm

from lazydock.algorithm.gamd_reweight import (
    calc_acceleration_factor, calc_anharmonicity, calc_diffusion_coefficient,
    calc_kramers_rate, find_barriers, parse_gamd_log, reweight_1d, reweight_1d_ce,
    reweight_1d_pyreweighting, reweight_2d, reweight_2d_pyreweighting, reweight_vdr,
    vdr_param,
)
from lazydock.scripts._script_utils_ import (Command, check_file_num_paried,
                                             process_batch_dir_lst)


class trajconv(Command):
    """
    Trajectory centering and PBC handling for Amber simulations.
    
    This command performs the standard PBC processing workflow:
    1. unwrap - Make molecules whole across periodic boundaries
    2. center_in_box - Center the protein in the simulation box
    3. wrap - Wrap solvent/ions back into the box
    """
    HELP = """Trajectory centering and PBC handling for Amber simulations"""
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="Directory containing simulation sub-folders. Each sub-folder should contain prmtop and trajectory files. Default is %(default)s.")
        args.add_argument('-p', '--prmtop-name', type=str, default='*.prmtop',
                          help='Topology file name pattern. Default is %(default)s.')
        args.add_argument('-t', '--traj-name', type=str, default='step06_md.nc',
                          help='Trajectory file name. Default is %(default)s.')
        args.add_argument('-c', '--center-group', type=str, default='protein',
                          help='Atom selection for centering. Default is %(default)s.')
        args.add_argument('-o', '--output-suffix', type=str, default='_center',
                          help='Suffix for output trajectory. Default is %(default)s.')
        args.add_argument('--unwrap', default=False, action='store_true',
                          help='Enable unwrap for systems with broken molecules across PBC. '
                               'Warning: This is very slow (~4s per frame). '
                               'Default is fast center+wrap only (~0.02s per frame).')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='Number of parallel workers. Default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='Force re-processing even if output exists.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='Delete existing output files before processing.')
        
    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        
    def check_prmtop_traj(self, bdir) -> Tuple[List[str], List[str]]:
        """Check prmtop and trajectory file pairing for a single batch directory"""
        prmtop_paths = []
        traj_paths = []
        
        # Find prmtop files - extract extension from pattern
        # Handle patterns like "*.prmtop" -> ["prmtop"]
        prmtop_ext = self.args.prmtop_name.replace('*.', '').split('.')[-1]
        p_paths = get_paths_with_extension(bdir, [prmtop_ext])
        
        for p_path in p_paths:
            p_path = Path(p_path)
            working_dir = p_path.parent
            
            # Find corresponding trajectory
            t_path = working_dir / self.args.traj_name
            if not t_path.exists():
                # Try glob pattern
                t_paths = list(working_dir.glob(self.args.traj_name))
                if t_paths:
                    t_path = t_paths[0]
                else:
                    put_err(f"No trajectory found for {p_path}, skipping.")
                    continue
            
            prmtop_paths.append(str(p_path))
            traj_paths.append(str(t_path))
        
        # Check pairing
        invalid_roots = check_file_num_paried(prmtop_paths, traj_paths)
        if invalid_roots:
            put_err(f"Unpaired prmtop/traj files found:\n{invalid_roots}", _exit=True)
        
        return prmtop_paths, traj_paths
    
    def find_tasks(self, bdir) -> List[Tuple[str, str, Path]]:
        """Find all tasks to process in a single batch directory"""
        prmtop_paths, traj_paths = self.check_prmtop_traj(bdir)
        
        tasks = []
        for p_path, t_path in zip(prmtop_paths, traj_paths):
            working_dir = Path(p_path).parent
            tasks.append((p_path, t_path, working_dir))
        
        return tasks
    
    @staticmethod
    def process_single_trajectory(prmtop: str, traj: str, output: str,
                                   center_group: str = 'protein',
                                   unwrap: bool = False,
                                   force: bool = False, delete: bool = False) -> bool:
        """
        Process a single trajectory with PBC handling using in-memory processing.
        
        Two modes available:
        1. FAST mode (default): center + wrap only (~0.02s per frame)
           Suitable for: single proteins, stable protein-ligand complexes
        
        2. SLOW mode (--unwrap): unwrap + center + wrap (~4s per frame)
           Suitable for: multi-domain proteins, flexible ligands, broken molecules
        
        Note: Uses in_memory=True for much faster processing.
        """
        output_path = Path(output)
        
        # Check existing output
        if output_path.exists():
            if delete:
                output_path.unlink()
                put_log(f"Deleted existing: {output}")
            elif not force:
                put_log(f"Output exists, skipping: {output}")
                return True
        
        try:
            put_log(f"Processing: {Path(traj).name}")
            put_log(f"  Mode: {'SLOW (unwrap+center+wrap)' if unwrap else 'FAST (center+wrap only)'}")
            
            # Load universe with in-memory processing for speed
            u = mda.Universe(prmtop, traj, in_memory=True)
            
            # Select atoms
            center_atoms = u.select_atoms(center_group)
            other_atoms = u.select_atoms(f'not ({center_group})')
            
            put_log(f"  Total atoms: {len(u.atoms)}, Center group: {len(center_atoms)}, Other: {len(other_atoms)}")
            put_log(f"  Total frames: {len(u.trajectory)}")
            put_log(f"  Applying PBC transformations...")
            
            # Apply transformations frame by frame
            for ts in tqdm(u.trajectory, desc="Transforming frames", leave=False):
                if unwrap:
                    # SLOW mode: unwrap molecules across PBC
                    trans.unwrap(u.atoms)(ts)
                
                # Center the selected group to box center
                trans.center_in_box(center_atoms, center='mass')(ts)
                
                # Wrap other atoms back into box
                # Note: compound='residues' is removed for better performance
                trans.wrap(other_atoms)(ts)
            
            # Write output trajectory
            put_log(f"  Writing output: {output}")
            with mda.Writer(output, n_atoms=len(u.atoms)) as W:
                for ts in tqdm(u.trajectory, desc="Writing frames", leave=False):
                    W.write(u.atoms)
            
            put_log(f"  Done: {output}")
            return True
            
        except Exception as e:
            put_err(f"Error processing {traj}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    def main_process(self):
        """Main processing loop - processes a single batch_dir"""
        # Get current batch_dir (set by excute() when using iter_run_arg)
        bdir = self.args.batch_dir
        
        # Find tasks
        tasks = self.find_tasks(bdir)
        put_log(f"Found {len(tasks)} task(s) to process in {bdir}")
        
        if not tasks:
            put_log("No tasks found. Exiting.")
            return
        
        # Process tasks
        if self.args.n_workers > 1:
            # Parallel processing
            pool = TaskPool('process', self.args.n_workers).start()
            for prmtop, traj, working_dir in tasks:
                output = working_dir / f"{Path(traj).stem}{self.args.output_suffix}.nc"
                pool.add_task(str(working_dir), self.process_single_trajectory,
                             prmtop, traj, str(output), self.args.center_group,
                             self.args.unwrap, self.args.force, self.args.delete)
            pool.close()
        else:
            # Sequential processing
            for prmtop, traj, working_dir in tqdm(tasks, desc="Processing trajectories"):
                output = working_dir / f"{Path(traj).stem}{self.args.output_suffix}.nc"
                self.process_single_trajectory(prmtop, traj, str(output),
                                               self.args.center_group,
                                               self.args.unwrap,
                                               self.args.force, self.args.delete)


class simple(Command):
    """
    Simple analysis for Amber MD simulations.
    
    Supported analyses:
    - rmsd: Root Mean Square Deviation
    - rmsf: Root Mean Square Fluctuation
    - gyrate: Radius of gyration
    - hbond: Hydrogen bond analysis
    - sasa: Solvent Accessible Surface Area
    - rdf: Radial Distribution Function
    - covar: Covariance matrix and PCA
    """
    HELP = """
    Simple analysis for Amber MD simulations
    
    1. RMSD calculation
    2. RMSF calculation  
    3. Radius of gyration
    4. Hydrogen bond analysis
    5. SASA calculation
    6. RDF calculation
    7. PCA (covariance analysis)
    """
    
    SUPPORT_METHODS = ['rmsd', 'rmsf', 'gyrate', 'hbond', 'sasa', 'rdf', 'covar']
    RECOMMEND_METHODS = list(set(SUPPORT_METHODS) - {'rdf', 'covar'})
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="Directory containing simulation sub-folders. Default is %(default)s.")
        args.add_argument('-p', '--prmtop-name', type=str, default='.prmtop',
                          help='Topology file name pattern. Default is %(default)s.')
        args.add_argument('-t', '--traj-name', type=str, default='step06_md_center.nc',
                          help='Trajectory file name. Default is %(default)s.')
        args.add_argument('--methods', type=str, nargs='+', default=simple.RECOMMEND_METHODS,
                          choices=simple.SUPPORT_METHODS,
                          help='Analysis methods to run. Default is all.')
        args.add_argument('-rg', '--rms-group', type=str, default='protein',
                          help='Atom selection for RMSD/RMSF/Rg. Default is %(default)s.')
        args.add_argument('-sg', '--sasa-group', type=str, default='protein',
                          help='Atom selection for SASA. Default is %(default)s.')
        args.add_argument('-hg', '--hbond-groups', type=str, nargs=2, default=['protein', 'protein'],
                          help='Two atom selections for H-bond analysis (donor, acceptor). Default is %(default)s.')
        args.add_argument('-rdfg', '--rdf-groups', type=str, nargs=2, default=['protein', 'resname WAT'],
                          help='Two atom selections for RDF. Default is %(default)s.')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='Number of parallel workers. Default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='Force re-analysis even if output exists.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='Delete existing analysis results.')
        args.add_argument('--time-offset', type=float, default=None,
                          help='Time offset in ps to subtract from trajectory times. '
                               'Auto-detect (subtract first frame time) if not set. '
                               'Set to 0 to keep original times.')
        
    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        
    def check_prmtop_traj(self, batch_dirs) -> Tuple[List[str], List[str]]:
        """Check prmtop and trajectory file pairing across all batch directories"""
        prmtop_paths = []
        traj_paths = []
        
        for bdir in batch_dirs:
            # Find prmtop files - extract extension from pattern
            prmtop_ext = self.args.prmtop_name.split('.')[-1]
            p_paths = get_paths_with_extension(bdir, [prmtop_ext])
            
            for p_path in p_paths:
                p_path = Path(p_path)
                working_dir = p_path.parent
                
                # Determine trajectory name
                traj_name = self.args.traj_name
                
                # Find corresponding trajectory
                t_path = working_dir / traj_name
                if not t_path.exists():
                    # Try glob pattern
                    t_paths = list(working_dir.glob(traj_name))
                    if t_paths:
                        t_path = t_paths[0]
                    else:
                        put_err(f"No trajectory found for {p_path}: {t_path}, skipping.")
                        continue
                
                prmtop_paths.append(str(p_path))
                traj_paths.append(str(t_path))
        
        # Check pairing
        invalid_roots = check_file_num_paried(prmtop_paths, traj_paths)
        if invalid_roots:
            put_err(f"Unpaired prmtop/traj files found:\n{invalid_roots}", _exit=True)
        
        return prmtop_paths, traj_paths
    
    def find_tasks(self, batch_dirs) -> List[Tuple[str, str, Path]]:
        """Find all tasks to analyze across all batch directories"""
        prmtop_paths, traj_paths = self.check_prmtop_traj(batch_dirs)
        
        tasks = []
        for p_path, t_path in zip(prmtop_paths, traj_paths):
            working_dir = Path(p_path).parent
            tasks.append((p_path, t_path, working_dir))
        
        return tasks
    
    @staticmethod
    def rmsd(prmtop: str, traj: str, working_dir: Path, 
             selection: str = 'protein', force: bool = False, delete: bool = False,
             u: Optional[mda.Universe] = None, time_offset: float = 0.0) -> bool:
        """Calculate RMSD"""
        output_csv = working_dir / f'rmsd.csv'
        output_png = working_dir / f'rmsd.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"RMSD output exists, skipping: {output_csv}")
                return True
        
        try:
            if u is None:
                u = mda.Universe(prmtop, traj)
            ref = mda.Universe(prmtop, traj)
            
            # Select atoms
            mobile = u.select_atoms(selection)
            reference = ref.select_atoms(selection)
            
            # Calculate RMSD
            rmsd_analysis = rms.RMSD(mobile, reference, select=selection)
            rmsd_analysis.run()
            
            # Get results
            times = rmsd_analysis.results.rmsd[:, 1] - time_offset  # time in ps, offset corrected
            rmsd_values = rmsd_analysis.results.rmsd[:, 2]  # RMSD in Å
            
            # Save to CSV
            df = pd.DataFrame({'Time (ps)': times, 'RMSD (Å)': rmsd_values})
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(10, 6))
            plt.plot(times / 1000, rmsd_values)  # Convert to ns
            plt.xlabel('Time (ns)')
            plt.ylabel('RMSD (Å)')
            plt.title(f'RMSD - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"RMSD saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"RMSD error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    @staticmethod
    def rmsf(prmtop: str, traj: str, working_dir: Path,
             selection: str = 'protein', force: bool = False, delete: bool = False,
             u: Optional[mda.Universe] = None) -> bool:
        """Calculate RMSF (per residue)"""
        output_csv = working_dir / f'rmsf.csv'
        output_png = working_dir / f'rmsf.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"RMSF output exists, skipping: {output_csv}")
                return True
        
        try:
            if u is None:
                u = mda.Universe(prmtop, traj)
            
            # Select atoms
            atoms = u.select_atoms(selection)
            
            # Calculate RMSF
            rmsf_analysis = rms.RMSF(atoms)
            rmsf_analysis.run()
            
            # Get residue IDs and RMSF values
            residues = atoms.residues
            residue_ids = residues.resids
            
            # Average RMSF per residue
            # Build mapping from global atom index to local index in selection
            index_map = {idx: i for i, idx in enumerate(atoms.indices)}
            rmsf_by_residue = []
            for residue in residues:
                res_atoms = residue.atoms
                local_indices = [index_map[idx] for idx in res_atoms.indices if idx in index_map]
                if local_indices:
                    rmsf_by_residue.append(np.mean(rmsf_analysis.results.rmsf[local_indices]))
                else:
                    rmsf_by_residue.append(0.0)
            
            # Save to CSV
            df = pd.DataFrame({'Residue ID': residue_ids, 'RMSF (Å)': rmsf_by_residue})
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(12, 6))
            plt.plot(residue_ids, rmsf_by_residue)
            plt.xlabel('Residue ID')
            plt.ylabel('RMSF (Å)')
            plt.title(f'RMSF - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"RMSF saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"RMSF error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    @staticmethod
    def gyrate(prmtop: str, traj: str, working_dir: Path,
               selection: str = 'protein', force: bool = False, delete: bool = False,
               u: Optional[mda.Universe] = None, time_offset: float = 0.0) -> bool:
        """Calculate radius of gyration"""
        output_csv = working_dir / f'gyrate.csv'
        output_png = working_dir / f'gyrate.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"Rg output exists, skipping: {output_csv}")
                return True
        
        try:
            if u is None:
                u = mda.Universe(prmtop, traj)
            
            # Select atoms
            atoms = u.select_atoms(selection)
            
            # Calculate Rg for each frame
            times = []
            rg_values = []
            
            for ts in u.trajectory:
                times.append(ts.time - time_offset)
                rg = atoms.radius_of_gyration()
                rg_values.append(rg)
            
            # Save to CSV
            df = pd.DataFrame({'Time (ps)': times, 'Rg (Å)': rg_values})
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(10, 6))
            plt.plot(np.array(times) / 1000, rg_values)  # Convert to ns
            plt.xlabel('Time (ns)')
            plt.ylabel('Radius of Gyration (Å)')
            plt.title(f'Rg - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"Rg saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"Rg error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    @staticmethod
    def hbond(prmtop: str, traj: str, working_dir: Path,
              groups: List[str] = ['protein', 'protein'], 
              force: bool = False, delete: bool = False,
              u: Optional[mda.Universe] = None, time_offset: float = 0.0) -> bool:
        """Calculate hydrogen bonds using MDAnalysis HydrogenBondAnalysis"""
        output_csv = working_dir / f'hbond.csv'
        output_png = working_dir / f'hbond.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"H-bond output exists, skipping: {output_csv}")
                return True
        
        try:
            # Load universe if not provided
            if u is None:
                u = mda.Universe(prmtop, traj)
            
            # H-bond analysis
            hbonds = HydrogenBondAnalysis(
                universe=u,
                donors_sel=groups[0],
                hydrogens_sel='name H*',
                acceptors_sel=groups[1],
                d_a_cutoff=3.5,
                d_h_a_angle_cutoff=150
            )
            hbonds.run()
            
            # Get times from trajectory
            times = [ts.time - time_offset for ts in u.trajectory]
            n_frames = len(times)
            
            # Count H-bonds per frame from results.hbonds array
            # results.hbonds shape: (N, 6) with [frame, donor, hydrogen, acceptor, distance, angle]
            counts_per_frame = [0] * n_frames
            if hbonds.results.hbonds is not None and len(hbonds.results.hbonds) > 0:
                for hbond in hbonds.results.hbonds:
                    frame_idx = int(hbond[0])
                    if 0 <= frame_idx < n_frames:
                        counts_per_frame[frame_idx] += 1
            
            # Save to CSV
            df = pd.DataFrame({'Time (ps)': times, 'H-bond Count': counts_per_frame})
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(10, 6))
            plt.plot(np.array(times) / 1000, counts_per_frame)
            plt.xlabel('Time (ns)')
            plt.ylabel('H-bond Count')
            plt.title(f'H-bonds - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"H-bonds saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"H-bond error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    @staticmethod
    def sasa(prmtop: str, traj: str, working_dir: Path,
             selection: str = 'protein', force: bool = False, delete: bool = False,
             u: Optional[mda.Universe] = None, time_offset: float = 0.0) -> bool:
        """Calculate SASA using MDTraj ShrakeRupley algorithm"""
        output_csv = working_dir / f'sasa.csv'
        output_png = working_dir / f'sasa.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"SASA output exists, skipping: {output_csv}")
                return True
        
        try:
            # Use MDTraj for SASA calculation (more reliable than MDAnalysis)
            # MDTraj expects coordinates in nm, probe_radius in nm (0.14 nm = 1.4 Å)
            put_log("Loading trajectory with MDTraj for SASA calculation...")
            t = md.load(traj, top=prmtop)
            
            # Select atoms using MDAnalysis for selection string support
            if u is None:
                u = mda.Universe(prmtop, traj)
            atoms = u.select_atoms(selection)
            atom_indices = atoms.indices  # 0-based indices
            
            # Get times from MDTraj trajectory
            times = t.time - time_offset  # in ps, offset corrected
            
            # Calculate SASA for each frame
            # mode='atom' returns shape (n_frames, n_selected_atoms)
            # We sum over atoms to get total SASA
            put_log("Calculating SASA (this may take a while)...")
            
            # Handle unknown atom elements (e.g., virtual sites 'VS', extra points 'EP')
            # MDTraj's shrake_rupley supports change_radii parameter
            change_radii = {
                'VS': 0.0,  # Virtual site, zero radius (doesn't contribute)
                'EP': 0.0,  # Extra point, zero radius
            }
            
            sasa_per_atom = md.shrake_rupley(t, probe_radius=0.14, n_sphere_points=960, 
                                             mode='atom', change_radii=change_radii,
                                             atom_indices=atom_indices)
            
            # Sum over atoms to get total SASA per frame (convert from nm² to Å²)
            sasa_total = np.sum(sasa_per_atom, axis=1) * 100  # nm² to Å²
            
            # Save to CSV
            df = pd.DataFrame({'Time (ps)': times, 'SASA (Å²)': sasa_total})
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(10, 6))
            plt.plot(times / 1000, sasa_total)
            plt.xlabel('Time (ns)')
            plt.ylabel('SASA (Å²)')
            plt.title(f'SASA - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"SASA saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"SASA error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    @staticmethod
    def rdf(prmtop: str, traj: str, working_dir: Path,
            groups: List[str] = ['protein', 'resname WAT'],
            force: bool = False, delete: bool = False,
            u: Optional[mda.Universe] = None) -> bool:
        """Calculate RDF between two groups"""
        output_csv = working_dir / f'rdf.csv'
        output_png = working_dir / f'rdf.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"RDF output exists, skipping: {output_csv}")
                return True
        
        try:
            if u is None:
                u = mda.Universe(prmtop, traj)
            
            # Select groups
            g1 = u.select_atoms(groups[0])
            g2 = u.select_atoms(groups[1])
            
            # Calculate RDF
            rdf = InterRDF(g1, g2, nbins=75, range=(0.0, 15.0))
            rdf.run()
            
            # Save to CSV
            df = pd.DataFrame({'r (Å)': rdf.results.bins, 'g(r)': rdf.results.rdf})
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(10, 6))
            plt.plot(rdf.results.bins, rdf.results.rdf)
            plt.xlabel('r (Å)')
            plt.ylabel('g(r)')
            plt.title(f'RDF - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"RDF saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"RDF error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    @staticmethod
    def covar(prmtop: str, traj: str, working_dir: Path,
              selection: str = 'protein', force: bool = False, delete: bool = False,
              u: Optional[mda.Universe] = None) -> bool:
        """Calculate covariance matrix and PCA"""
        output_csv = working_dir / f'pca_eigenval.csv'
        output_png = working_dir / f'pca_eigenval.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"PCA output exists, skipping: {output_csv}")
                return True
        
        try:
            from MDAnalysis.analysis.pca import PCA
            
            if u is None:
                u = mda.Universe(prmtop, traj)
            
            # Select atoms and align trajectory
            u.select_atoms(selection)
            
            # Align trajectory
            align.AlignTraj(u, u, select=selection, in_memory=True).run()
            
            # PCA
            pca = PCA(u, select=selection)
            pca.run()
            
            # Get eigenvalues
            eigenvalues = pca.results.eigenvalues
            n_eigen = min(10, len(eigenvalues))
            
            # Save to CSV
            df = pd.DataFrame({
                'Component': range(1, n_eigen + 1),
                'Eigenvalue': eigenvalues[:n_eigen],
                'Variance (%)': pca.results.variance[:n_eigen] * 100
            })
            df.to_csv(output_csv, index=False)
            
            # Plot
            plt.figure(figsize=(10, 6))
            plt.bar(range(1, n_eigen + 1), eigenvalues[:n_eigen])
            plt.xlabel('Principal Component')
            plt.ylabel('Eigenvalue')
            plt.title(f'PCA Eigenvalues - {working_dir.name}')
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"PCA saved: {output_csv}, {output_png}")
            return True
            
        except Exception as e:
            put_err(f"PCA error for {working_dir}: {e}")
            import traceback
            put_err(traceback.format_exc())
            return False
    
    def main_process(self):
        """Main analysis loop"""
        # Find tasks across all batch directories
        tasks = self.find_tasks(self.args.batch_dir)
        put_log(f"Found {len(tasks)} task(s) to analyze")
        
        if not tasks:
            put_log("No tasks found. Exiting.")
            return
        
        # Run analyses for each task
        method_map = {
            'rmsd': self.rmsd,
            'rmsf': self.rmsf,
            'gyrate': self.gyrate,
            'hbond': self.hbond,
            'sasa': self.sasa,
            'rdf': self.rdf,
            'covar': self.covar,
        }
        
        for prmtop, traj, working_dir in tqdm(tasks, desc="Analyzing trajectories"):
            put_log(f"\nAnalyzing: {working_dir.name}")
            
            # Pre-load universe once for methods that use it (except SASA which uses MDTraj)
            u = None
            if 'sasa' not in self.args.methods:
                try:
                    put_log("  Loading universe (shared across analyses)...")
                    u = mda.Universe(prmtop, traj)
                except Exception as e:
                    put_err(f"  Failed to load universe: {e}")
                    continue
            
            # Determine time offset
            time_offset = 0.0
            if self.args.time_offset is None and u is not None:
                # Auto-detect: subtract first frame time so times start from 0
                first_time = u.trajectory[0].time
                if first_time > 0:
                    time_offset = first_time
                    put_log(f"  Auto-detected time offset: {time_offset:.1f} ps "
                            f"(first frame time corrected to 0)")
            elif self.args.time_offset is not None:
                time_offset = self.args.time_offset
                if time_offset > 0:
                    put_log(f"  Using specified time offset: {time_offset:.1f} ps")
            
            for method in self.args.methods:
                if method not in method_map:
                    continue
                
                put_log(f"  Running: {method}")
                
                # Prepare arguments
                kwargs = {
                    'prmtop': prmtop,
                    'traj': traj,
                    'working_dir': working_dir,
                    'force': self.args.force,
                    'delete': self.args.delete,
                }
                
                # Add universe for methods that need it (except SASA uses MDTraj)
                if method != 'sasa':
                    kwargs['u'] = u
                
                # Add method-specific arguments
                if method in ['rmsd', 'rmsf', 'gyrate']:
                    kwargs['selection'] = self.args.rms_group
                    kwargs['time_offset'] = time_offset
                elif method == 'hbond':
                    kwargs['groups'] = self.args.hbond_groups
                    kwargs['time_offset'] = time_offset
                elif method == 'sasa':
                    kwargs['selection'] = self.args.sasa_group
                    kwargs['time_offset'] = time_offset
                elif method == 'rdf':
                    kwargs['groups'] = self.args.rdf_groups
                elif method == 'covar':
                    kwargs['selection'] = self.args.rms_group
                
                # Run analysis
                try:
                    method_map[method](**kwargs)
                except Exception as e:
                    put_err(f"    Error in {method}: {e}")


class gamd(Command):
    """
    GaMD-specific analysis for Amber simulations.
    
    Supported analyses:
    - anharmonicity: Check anharmonicity of boost potential ΔV distribution
    - reweight: Energetic reweighting (1D/2D PMF via cumulant expansion)
    - kinetic: Kinetic reweighting via Kramers rate theory
    """
    HELP = """
    GaMD-specific analysis for Amber simulations
    
    1. anharmonicity: Check anharmonicity of boost potential ΔV distribution
       - Computes skewness (γ) of ΔV distribution
       - Fits Gaussian and compares with actual distribution
       - γ < 0.1: Excellent; 0.1-0.3: Acceptable; > 0.3: Poor
    
    2. reweight: Energetic reweighting via cumulant expansion
       - 1D PMF along a single CV
       - 2D PMF along two CVs
       - Requires CV data file (e.g., from simple analysis)
    
    3. kinetic: Kinetic reweighting via Kramers rate theory
       - Estimates barrier crossing rates from reweighted PMF
       - Computes acceleration factor
    """
    
    SUPPORT_METHODS = ['anharmonicity', 'reweight', 'kinetic']
    RECOMMEND_METHODS = ['anharmonicity', 'reweight']
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="Directory containing simulation sub-folders. Default is %(default)s.")
        args.add_argument('--gamd-log', type=str, default='gamd.log',
                          help='GaMD log file name in each working directory. Default is %(default)s.')
        args.add_argument('--methods', type=str, nargs='+', default=gamd.RECOMMEND_METHODS,
                          choices=gamd.SUPPORT_METHODS,
                          help='GaMD analysis methods to run. Default is %(default)s.')
        args.add_argument('--cv-file', type=str, default=None,
                          help='CV data file for reweighting (CSV with Time + CV columns). '
                               'If not set, will look for rmsd.csv in working directory.')
        args.add_argument('--cv-cols', type=str, nargs='+', default=None,
                          help='Column names in CV file to use as reaction coordinates. '
                               '1 column for 1D PMF, 2 columns for 2D PMF. '
                               'If not set, uses the first non-Time column(s).')
        args.add_argument('--temp', type=float, default=300.0,
                          help='Simulation temperature in Kelvin. Default is %(default)s.')
        args.add_argument('--Emax', type=float, default=8.0,
                          help='Max free energy for PMF display (kcal/mol). Default is %(default)s.')
        args.add_argument('--cutoff', type=int, default=10,
                          help='Minimum number of frames per bin for reweighting. Default is %(default)s.')
        args.add_argument('--disc', type=float, default=None,
                          help='Bin size for CV. Auto-determined if not set.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='Force re-analysis even if output exists.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='Delete existing analysis results.')
    
    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
    
    def find_tasks(self, batch_dirs) -> List[Tuple[Path]]:
        """Find all working directories containing gamd.log"""
        tasks = []
        for bdir in batch_dirs:
            bdir = Path(bdir)
            for gamd_log in bdir.rglob(self.args.gamd_log):
                working_dir = gamd_log.parent
                tasks.append((working_dir,))
        return tasks
    
    @staticmethod
    def anharmonicity(working_dir: Path, gamd_log_name: str = 'gamd.log',
                      force: bool = False, delete: bool = False,
                      temp: float = 300.0) -> bool:
        """Check anharmonicity of boost potential ΔV distribution"""
        output_csv = working_dir / 'gamd_anharmonicity.csv'
        output_png = working_dir / 'gamd_anharmonicity.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"Anharmonicity output exists, skipping: {output_csv}")
                return True
        
        gamd_log_path = working_dir / gamd_log_name
        df = parse_gamd_log(gamd_log_path)
        if df is None or len(df) == 0:
            put_err(f"No gamd.log data found in {working_dir}")
            return False
        
        # Call algorithm
        result = calc_anharmonicity(
            df['boost_energy_potential'].values,
            df['boost_energy_dihedral'].values,
            temp=temp
        )
        
        quality = result['total']['quality']
        gamma = result['total']['anharmonicity']
        skewness = result['total']['skewness']
        kappa = result['total']['excess_kurtosis']
        mu = result['total']['mean']
        sigma = result['total']['std']
        put_log(f"  Anharmonicity: γ={gamma:.6f}, skewness={skewness:.4f}, κ={kappa:.4f}, "
                f"μ={mu:.2f}, σ={sigma:.2f} → {quality}")
        
        # Save CSV
        result_df = pd.DataFrame({
            'Component': ['Total (ΔV_P + ΔV_D)', 'Potential (ΔV_P)', 'Dihedral (ΔV_D)'],
            'Mean': [result['total']['mean'], result['potential']['mean'], result['dihedral']['mean']],
            'Std': [result['total']['std'], result['potential']['std'], result['dihedral']['std']],
            'Anharmonicity (γ)': [result['total']['anharmonicity'], result['potential']['anharmonicity'], result['dihedral']['anharmonicity']],
            'Skewness': [result['total']['skewness'], result['potential']['skewness'], result['dihedral']['skewness']],
            'Excess_Kurtosis (κ)': [result['total']['excess_kurtosis'], result['potential']['excess_kurtosis'], result['dihedral']['excess_kurtosis']],
            'Quality': [result['total']['quality'], result['potential']['quality'], result['dihedral']['quality']],
        })
        result_df.to_csv(output_csv, index=False)
        
        # Plot: ΔV distribution vs Gaussian fit
        from scipy.stats import norm
        
        delta_V = df['boost_energy_potential'].values + df['boost_energy_dihedral'].values
        dV_P = df['boost_energy_potential'].values
        dV_D = df['boost_energy_dihedral'].values
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        for ax, data, label, color, comp_key in [
            (axes[0], delta_V, 'Total ΔV (ΔV_P + ΔV_D)', 'steelblue', 'total'),
            (axes[1], dV_P, 'Potential ΔV_P', 'coral', 'potential'),
            (axes[2], dV_D, 'Dihedral ΔV_D', 'seagreen', 'dihedral'),
        ]:
            data_mu = np.mean(data)
            data_sigma = np.std(data)
            comp_gamma = result[comp_key]['anharmonicity']
            comp_skew = result[comp_key]['skewness']
            
            n_bins = min(100, max(30, len(data) // 50))
            ax.hist(data, bins=n_bins, density=True, alpha=0.6, color=color, label='Actual')
            
            x = np.linspace(data_mu - 4 * data_sigma, data_mu + 4 * data_sigma, 300)
            ax.plot(x, norm.pdf(x, data_mu, data_sigma), 'k--', linewidth=2, label='Gaussian fit')
            
            ax.set_xlabel(f'{label} (kcal/mol)')
            ax.set_ylabel('Probability Density')
            ax.set_title(f'{label}\nγ={comp_gamma:.6f}, skew={comp_skew:.4f}')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        fig.suptitle(f'GaMD Anharmonicity Assessment - {working_dir.name}\n'
                     f'Overall Quality: {quality}', fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_png, dpi=150)
        plt.close()
        
        put_log(f"Anharmonicity saved: {output_csv}, {output_png}")
        return True
    
    @staticmethod
    def reweight(working_dir: Path, gamd_log_name: str = 'gamd.log',
                 cv_file: Optional[str] = None, cv_cols: Optional[List[str]] = None,
                 temp: float = 300.0, Emax: float = 8.0,
                 cutoff: int = 10, disc: Optional[float] = None,
                 force: bool = False, delete: bool = False) -> bool:
        """Energetic reweighting via cumulant expansion to 2nd order"""
        is_2d = cv_cols is not None and len(cv_cols) >= 2
        suffix = '2d' if is_2d else '1d'
        output_csv = working_dir / f'gamd_pmf_{suffix}.csv'
        output_png = working_dir / f'gamd_pmf_{suffix}.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"Reweighting output exists, skipping: {output_csv}")
                return True
        
        # Parse gamd.log
        gamd_log_path = working_dir / gamd_log_name
        gamd_df = parse_gamd_log(gamd_log_path)
        if gamd_df is None or len(gamd_df) == 0:
            put_err(f"No gamd.log data found in {working_dir}")
            return False
        
        delta_V = gamd_df['boost_energy_potential'].values + gamd_df['boost_energy_dihedral'].values
        n_frames_gamd = len(delta_V)
        
        # Load CV data
        cv_path = cv_file
        if cv_path is None:
            cv_path = working_dir / 'rmsd.csv'
        else:
            cv_path = working_dir / cv_path if not Path(cv_file).is_absolute() else Path(cv_file)
        
        if not cv_path.exists():
            put_err(f"CV file not found: {cv_path}. "
                    f"Run 'simple --methods rmsd' first or specify --cv-file.")
            return False
        
        cv_df = pd.read_csv(cv_path)
        n_frames_cv = len(cv_df)
        
        # Determine CV columns
        if cv_cols is None:
            time_cols = [c for c in cv_df.columns if 'time' in c.lower()]
            non_time_cols = [c for c in cv_df.columns if c not in time_cols]
            if not non_time_cols:
                put_err(f"No non-time columns found in {cv_path}")
                return False
            cv_cols = non_time_cols[:2] if is_2d else non_time_cols[:1]
        
        # Align data
        n_frames = min(n_frames_gamd, n_frames_cv)
        if n_frames_gamd != n_frames_cv:
            put_log(f"  Frame count mismatch: gamd.log={n_frames_gamd}, CV file={n_frames_cv}. "
                    f"Using first {n_frames} frames.")
        
        delta_V = delta_V[:n_frames]
        cv_data = [cv_df[col].values[:n_frames] for col in cv_cols]
        
        if len(cv_cols) == 1:
            # 1D reweighting
            pmf_result = reweight_1d(cv_data[0], delta_V, temp=temp, cutoff=cutoff, disc=disc)
            if pmf_result is None:
                put_err("1D reweighting failed")
                return False
            
            bin_centers, pmf_values, counts, c1_values, c2_values = pmf_result
            
            result_df = pd.DataFrame({
                cv_cols[0]: bin_centers,
                'PMF (kcal/mol)': pmf_values,
                'Frame_Count': counts,
                'C1 (kcal/mol)': c1_values,
                'C2 (kcal²/mol²)': c2_values,
            })
            result_df.to_csv(output_csv, index=False)
            
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[3, 1])
            valid = (~np.isnan(pmf_values)) & (pmf_values <= Emax)
            ax1.plot(bin_centers[valid], pmf_values[valid], 'b-', linewidth=1.5)
            ax1.set_xlabel(cv_cols[0])
            ax1.set_ylabel('PMF (kcal/mol)')
            ax1.set_title(f'GaMD Reweighted PMF - {working_dir.name}')
            ax1.set_ylim(bottom=np.nanmin(pmf_values[valid]) - 0.5, top=Emax + 0.5)
            ax1.grid(True, alpha=0.3)
            
            ax2.bar(bin_centers, counts, width=np.diff(np.append(bin_centers, 2*bin_centers[-1]-bin_centers[-2]))[0] * 0.8)
            ax2.set_xlabel(cv_cols[0])
            ax2.set_ylabel('Frame Count')
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"1D PMF saved: {output_csv}, {output_png}")
            return True
        else:
            # 2D reweighting
            pmf_result = reweight_2d(cv_data[0], cv_data[1], delta_V, temp=temp, cutoff=cutoff, disc=disc)
            if pmf_result is None:
                put_err("2D reweighting failed")
                return False
            
            x_centers, y_centers, pmf_grid, counts_grid = pmf_result
            
            xx, yy = np.meshgrid(x_centers, y_centers, indexing='ij')
            result_df = pd.DataFrame({
                cv_cols[0]: xx.flatten(),
                cv_cols[1]: yy.flatten(),
                'PMF (kcal/mol)': pmf_grid.flatten(),
                'Frame_Count': counts_grid.flatten(),
            })
            result_df.to_csv(output_csv, index=False)
            
            fig, ax = plt.subplots(figsize=(10, 8))
            pmf_plot = pmf_grid.copy()
            pmf_plot[pmf_plot > Emax] = np.nan
            im = ax.pcolormesh(x_centers, y_centers, pmf_plot.T, shading='auto', cmap='jet')
            plt.colorbar(im, ax=ax, label='PMF (kcal/mol)')
            ax.set_xlabel(cv_cols[0])
            ax.set_ylabel(cv_cols[1])
            ax.set_title(f'GaMD Reweighted 2D PMF - {working_dir.name}')
            plt.tight_layout()
            plt.savefig(output_png, dpi=150)
            plt.close()
            
            put_log(f"2D PMF saved: {output_csv}, {output_png}")
            return True
    
    @staticmethod
    def kinetic(working_dir: Path, gamd_log_name: str = 'gamd.log',
                cv_file: Optional[str] = None, cv_cols: Optional[List[str]] = None,
                temp: float = 300.0, cutoff: int = 10, disc: Optional[float] = None,
                force: bool = False, delete: bool = False) -> bool:
        """Kinetic reweighting via Kramers rate theory"""
        output_csv = working_dir / 'gamd_kinetic.csv'
        output_png = working_dir / 'gamd_kinetic.png'
        
        if output_csv.exists():
            if delete:
                output_csv.unlink()
                output_png.unlink(missing_ok=True)
            elif not force:
                put_log(f"Kinetic output exists, skipping: {output_csv}")
                return True
        
        # Parse gamd.log
        gamd_log_path = working_dir / gamd_log_name
        gamd_df = parse_gamd_log(gamd_log_path)
        if gamd_df is None or len(gamd_df) == 0:
            put_err(f"No gamd.log data found in {working_dir}")
            return False
        
        delta_V = gamd_df['boost_energy_potential'].values + gamd_df['boost_energy_dihedral'].values
        n_frames_gamd = len(delta_V)
        
        # Load CV data
        cv_path = cv_file
        if cv_path is None:
            cv_path = working_dir / 'rmsd.csv'
        else:
            cv_path = working_dir / cv_file if not Path(cv_file).is_absolute() else Path(cv_file)
        
        if not cv_path.exists():
            put_err(f"CV file not found: {cv_path}. "
                    f"Run 'simple --methods rmsd' first or specify --cv-file.")
            return False
        
        cv_df = pd.read_csv(cv_path)
        n_frames_cv = len(cv_df)
        
        # Determine CV column
        if cv_cols is None:
            time_cols = [c for c in cv_df.columns if 'time' in c.lower()]
            non_time_cols = [c for c in cv_df.columns if c not in time_cols]
            if not non_time_cols:
                put_err(f"No non-time columns found in {cv_path}")
                return False
            cv_col = non_time_cols[0]
        else:
            cv_col = cv_cols[0]
        
        # Align data
        n_frames = min(n_frames_gamd, n_frames_cv)
        if n_frames_gamd != n_frames_cv:
            put_log(f"  Frame count mismatch: gamd.log={n_frames_gamd}, CV file={n_frames_cv}. "
                    f"Using first {n_frames} frames.")
        
        delta_V = delta_V[:n_frames]
        cv = cv_df[cv_col].values[:n_frames]
        
        # Step 1: Compute reweighted 1D PMF
        pmf_result = reweight_1d(cv, delta_V, temp=temp, cutoff=cutoff, disc=disc)
        if pmf_result is None:
            put_err("Reweighting failed for kinetic analysis")
            return False
        
        bin_centers, pmf_values, counts, _, _ = pmf_result
        
        # Step 2: Find barriers
        barrier_info = find_barriers(bin_centers, pmf_values)
        
        if barrier_info is None:
            put_log("  No clear barrier found in PMF. Single-well potential detected.")
            valid = ~np.isnan(pmf_values)
            result_df = pd.DataFrame({
                'Property': ['Barrier detected', 'Min PMF', 'Max PMF'],
                'Value': ['No', float(np.nanmin(pmf_values[valid])), float(np.nanmax(pmf_values[valid]))],
            })
            result_df.to_csv(output_csv, index=False)
            return True
        
        # Step 3: Compute diffusion coefficient
        D_apparent = calc_diffusion_coefficient(cv)
        
        # Step 4: Compute Kramers rate
        kramers_result = calc_kramers_rate(barrier_info, D_apparent, temp=temp)
        
        # Step 5: Compute acceleration factor
        accel_factor = calc_acceleration_factor(
            delta_V, cv,
            barrier_cv=barrier_info['barrier_cv'],
            well_cv=barrier_info['well1_cv'],
            temp=temp
        )
        kramers_result['acceleration_factor'] = accel_factor
        if accel_factor > 0 and not np.isnan(kramers_result['k_gamd']):
            kramers_result['k_cMD_estimate'] = kramers_result['k_gamd'] / accel_factor
        
        # Build results dict
        results = {
            'cv_col': cv_col,
            'well1_cv': barrier_info['well1_cv'],
            'well1_pmf': barrier_info['well1_pmf'],
            'well2_cv': barrier_info['well2_cv'],
            'well2_pmf': barrier_info['well2_pmf'],
            'barrier_cv': barrier_info['barrier_cv'],
            'barrier_pmf': barrier_info['barrier_pmf'],
            'delta_G_well1_kcal': barrier_info['delta_G_well1'],
            'delta_G_well2_kcal': barrier_info['delta_G_well2'],
            'curvature_well': barrier_info['curvature_well'],
            'curvature_barrier': barrier_info['curvature_barrier'],
            'D_apparent': D_apparent,
            'omega_min': kramers_result['omega_min'],
            'omega_barrier': kramers_result['omega_barrier'],
            'k_gamd_arbitrary': kramers_result['k_gamd'],
            'acceleration_factor': kramers_result['acceleration_factor'],
            'k_cMD_estimate': kramers_result.get('k_cMD_estimate', np.nan),
        }
        
        delta_G = barrier_info['delta_G_well1']
        if not np.isnan(kramers_result['k_gamd']):
            put_log(f"  Barrier: ΔG‡={delta_G:.2f} kcal/mol at CV={barrier_info['barrier_cv']:.2f}")
            put_log(f"  Acceleration factor: {accel_factor:.2f}")
        else:
            put_log(f"  Barrier: ΔG‡={delta_G:.2f} kcal/mol, but curvature/D not suitable for Kramers estimate")
        
        # Save CSV
        result_df = pd.DataFrame({
            'Property': list(results.keys()),
            'Value': list(results.values()),
        })
        result_df.to_csv(output_csv, index=False)
        
        # Plot PMF with annotated barriers
        valid = ~np.isnan(pmf_values)
        bc_valid = bin_centers[valid]
        pmf_valid = pmf_values[valid]
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), height_ratios=[3, 1])
        
        ax1.plot(bc_valid, pmf_valid, 'b-', linewidth=1.5)
        ax1.plot(barrier_info['well1_cv'], barrier_info['well1_pmf'], 'go', markersize=10,
                 label=f'Well 1 ({barrier_info["well1_cv"]:.2f})')
        ax1.plot(barrier_info['well2_cv'], barrier_info['well2_pmf'], 'bs', markersize=10,
                 label=f'Well 2 ({barrier_info["well2_cv"]:.2f})')
        ax1.plot(barrier_info['barrier_cv'], barrier_info['barrier_pmf'], 'r^', markersize=12,
                 label=f'Barrier ({delta_G:.2f} kcal/mol)')
        ax1.vlines(barrier_info['barrier_cv'], barrier_info['well1_pmf'], barrier_info['barrier_pmf'],
                    colors='red', linestyles='dashed', alpha=0.7)
        ax1.set_xlabel(cv_col)
        ax1.set_ylabel('PMF (kcal/mol)')
        ax1.set_title(f'GaMD Kinetic Analysis - {working_dir.name}')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        n_show = min(10000, n_frames)
        time_ps = cv_df[cv_df.columns[0]].values[:n_show] if len(cv_df.columns) > 0 else np.arange(n_show)
        ax2.plot(time_ps, cv[:n_show], 'gray', linewidth=0.3, alpha=0.7)
        ax2.set_xlabel('Time (ps)')
        ax2.set_ylabel(cv_col)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_png, dpi=150)
        plt.close()
        
        put_log(f"Kinetic analysis saved: {output_csv}, {output_png}")
        return True
    
    def main_process(self):
        """Main analysis loop"""
        tasks = self.find_tasks(self.args.batch_dir)
        put_log(f"Found {len(tasks)} task(s) to analyze")
        
        if not tasks:
            put_log("No tasks found. Exiting.")
            return
        
        method_map = {
            'anharmonicity': self.anharmonicity,
            'reweight': self.reweight,
            'kinetic': self.kinetic,
        }
        
        for (working_dir,) in tqdm(tasks, desc="GaMD Analysis"):
            put_log(f"\nAnalyzing: {working_dir.name}")
            
            for method in self.args.methods:
                if method not in method_map:
                    continue
                
                put_log(f"  Running: {method}")
                
                kwargs = {
                    'working_dir': working_dir,
                    'gamd_log_name': self.args.gamd_log,
                    'force': self.args.force,
                    'delete': self.args.delete,
                    'temp': self.args.temp,
                }
                
                if method in ['reweight', 'kinetic']:
                    kwargs['cv_file'] = self.args.cv_file
                    kwargs['cv_cols'] = self.args.cv_cols
                    kwargs['cutoff'] = self.args.cutoff
                    kwargs['disc'] = self.args.disc
                
                if method == 'reweight':
                    kwargs['Emax'] = self.args.Emax
                
                try:
                    method_map[method](**kwargs)
                except Exception as e:
                    put_err(f"    Error in {method}: {e}")
                    import traceback
                    put_err(traceback.format_exc())


def main(sys_args: List[str] = None):
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='Amber Trajectory Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # trajconv command
    trajconv_parser = subparsers.add_parser('trajconv', help='Trajectory centering and PBC handling')
    trajconv.make_args(trajconv_parser)
    
    # simple command
    simple_parser = subparsers.add_parser('simple', help='Simple analysis (RMSD, RMSF, Rg, etc.)')
    simple.make_args(simple_parser)
    
    # gamd command
    gamd_parser = subparsers.add_parser('gamd', help='GaMD analysis (anharmonicity, reweighting, kinetics)')
    gamd.make_args(gamd_parser)
    
    args = parser.parse_args(sys_args)
    
    if args.command == 'trajconv':
        cmd = trajconv(args)
        cmd.excute()
    elif args.command == 'simple':
        cmd = simple(args)
        cmd.excute()
    elif args.command == 'gamd':
        cmd = gamd(args)
        cmd.excute()
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
