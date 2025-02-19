'''
Date: 2025-01-16 10:08:37
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-19 22:26:33
Description: 
'''
import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.plot import save_show
from MDAnalysis.analysis import dihedrals, gnm, helix_analysis, diffusionmap, align, rms
from tqdm import tqdm

from lazydock.scripts._script_utils_ import excute_command
from lazydock.scripts.ana_gmx import mmpbsa


def smv(arr: np.ndarray, w: int = 50):
    return np.convolve(arr, np.ones(w), "valid") / w


class simple(mmpbsa):
    HELP = """
    simple analysis collections for MDAnalysis.
    input is centered trajectory file, such as md_center.xtc.
    
    1. elastic network， using a Gaussian network model with only close contacts, analyzed the thermal fluctuation behavior of the system
    2. average twist of the helix
    3. Radial distribution function of specific residue(s)
    4. Ramachandran plot
    5. Janin plot 
    6. Inter-frame pairwise RMSD
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('-top', '--top-name', type = str, required=True,
                          help='topology file name in each sub-directory, such as md.tpr.')
        args.add_argument('-traj', '--traj-name', type = str, required=True,
                          help='trajectory file name in each sub-directory, such as md_center.xtc.')
        args.add_argument('--elastic-select', type = str, default='name CA',
                          help='elastic select group, default is %(default)s.')
        args.add_argument('--twist-select', type = str, nargs='+', default = [],
                          help='twist select group, can be multiple MDAnalysis selection string, default is %(default)s.')
        args.add_argument('--radial-select', type = str, nargs='+', default = [],
                          help='radial select group, can be multiple MDAnalysis selection string, default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='force to re-run the analysis, default is %(default)s.')
        return args
        
    @staticmethod
    def elastic(u: mda.Universe, w_dir: Path, elastic_select: str, force: bool = False):
        if os.path.exists(w_dir / 'elastic.xlsx') and not force:
            return put_log('Elastic analysis already calculated, use -F to re-run.')
        put_log('Calculating elastic network')
        nma1 = gnm.closeContactGNMAnalysis(u, select=elastic_select, cutoff=7.0, weights='size')
        nma1.run()
        df = pd.DataFrame(nma1.results)
        df.to_excel(w_dir / 'elastic.xlsx', index=False)
        fig = plt.figure(figsize=(10, 8))
        sns.lineplot(data=df, x='times', y='eigenvalues', alpha=0.6, ax=fig.gca())
        sns.lineplot(x=df['times'][50:], y=smv(df['eigenvalues']), label='SMV', ax=fig.gca())
        save_show(w_dir / 'elastic.png', 600, show=False)
        plt.close(fig=fig)
        del nma1
        
    @staticmethod
    def twist(u: mda.Universe, w_dir: Path, twist_select: List[str]):
        dfs = {}
        for i, twist_sel in enumerate(twist_select):
            h = helix_analysis.HELANAL(u, select=twist_sel,
                                    ref_axis=[0, 0, 1]).run()
            dfs[twist_sel] = h.results.local_twists
            twi_i = dfs[twist_sel].mean(axis=1)
            fig = plt.figure(figsize=(10, 8))
            sns.lineplot(x=np.arange(len(twi_i)), y=twi_i, alpha=0.6, ax=fig.gca())
            sns.lineplot(x=np.arange(50, len(twi_i-50)), y=smv(twi_i), label='SMV', ax=fig.gca())
            save_show(w_dir / f'twist_{i}.png', 600, show=False)
            plt.close(fig=fig)
        opts_file(w_dir / 'twist.pkl', 'wb', way='pkl', data=dfs)
            
    @staticmethod
    def radial_dist(u: mda.Universe, w_dir: Path, radial_select: List[str]):
        dfs = {}
        for i, radial_sel in enumerate(radial_select):
            r = helix_analysis.RadialDistribution(u, select=radial_sel, bins=100).run()
            dfs[radial_sel] = r.results.radial_distribution
            fig = plt.figure(figsize=(10, 8))
            sns.lineplot(x=r.results.bins, y=r.results.radial_distribution, alpha=0.6, ax=fig.gca())
            save_show(w_dir / f'radial_{i}.png', 600, show=False)
            plt.close(fig=fig)
        opts_file(w_dir / 'radial.pkl', 'wb', way='pkl', data=dfs)
        
    @staticmethod
    def ramachandran(u: mda.Universe, w_dir: Path, force: bool = False):
        if os.path.exists(w_dir / 'ramachandran.npz') and not force:
            return put_log('Ramachandran plot already calculated, use -F to re-run.')
        protein = gnm.GNMAnalysis(u, select='protein', cutoff=7.0)
        rama = dihedrals.Ramachandran(protein).run()
        np.savez_compressed(w_dir / 'ramachandran.npz', angles = rama.results.angles)
        fig, ax = plt.subplots(figsize=(10, 8))
        rama.plot(color='black', marker='.', ref=True, ax=ax)
        save_show(w_dir / 'ramachandran.png', 600, show=False)
        plt.close(fig=fig)
    
    @staticmethod
    def janin(u: mda.Universe, w_dir: Path, force: bool = False):
        if os.path.exists(w_dir / 'janin.npz') and not force:
            return put_log('Janin plot already calculated, use -F to re-run.')
        protein = gnm.GNMAnalysis(u, select='protein', cutoff=7.0)
        janin = dihedrals.Janin(protein).run()
        np.savez_compressed(w_dir / 'janin.npz', angles = janin.results.angles)
        fig, ax = plt.subplots(figsize=(10, 8))
        janin.plot(color='black', marker='.', ref=True, ax=ax)
        save_show(w_dir / 'janin.png', 600, show=False)
        plt.close(fig=fig)
    
    @staticmethod
    def inter_frame_rmsd(u: mda.Universe, w_dir: Path, force: bool = False):
        if os.path.exists(w_dir / 'inter_frame_rmsd.pkl') and not force:
            return put_log('Inter-frame RMSD already calculated, use -F to re-run.')
        put_log('Aligning inter-frame RMSD')
        aligner = align.AlignTraj(u, u, select='name CA', in_memory=True).run()
        put_log('Calculating inter-frame RMSD')
        matrix = diffusionmap.DistanceMatrix(u, select='name CA').run()
        opts_file(w_dir / 'inter_frame_rmsd.pkl', 'wb', way='pkl', data=matrix.results.dist_matrix)
        plt.imshow(matrix.results.dist_matrix, cmap='viridis')
        plt.xlabel('Frame')
        plt.ylabel('Frame')
        plt.colorbar(label=r'RMSD ($\AA$)')
        save_show(w_dir / 'inter_frame_rmsd.png', 600, show=False)
        del aligner, matrix
        
    def main_process(self):
        print(f'find {len(self.tasks)} tasks.')
        # process each task
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            complex_path = Path(complex_path).resolve()
            u = mda.Universe(top_path, traj_path)
            # elastic network
            self.elastic(u, wdir, self.args.elastic_select, self.args.force)
            # twist
            self.twist(u, wdir, self.args.twist_select)
            # radial distribution
            self.radial_dist(u, wdir, self.args.radial_select)
            # ramachandran
            self.ramachandran(u, wdir, self.args.force)
            # janin
            self.janin(u, wdir, self.args.force)
            # inter-frame RMSD
            self.inter_frame_rmsd(u, wdir, self.args.force)
            bar.update(1)
        bar.close()


_str2func = {
    'simple': simple,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for MDAnalysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    main()