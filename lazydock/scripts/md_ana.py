'''
Date: 2025-01-16 10:08:37
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-21 20:37:52
Description: 
'''
import argparse
import os
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Union

from Bio import BiopythonDeprecationWarning

warnings.simplefilter('ignore', BiopythonDeprecationWarning)
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
from lazydock.gmx.mda.gnm import (calcu_closeContactGNMAnalysis,
                                  calcu_GNMAnalysis, genarate_atom2residue)
from lazydock.gmx.mda.rms import pairwise_rmsd
from lazydock.scripts._script_utils_ import (clean_path, excute_command,
                                             process_batch_dir_lst)
from lazydock.scripts.ana_gmx import mmpbsa
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.plot import save_show
from mbapy_lite.web_utils.task import TaskPool
from MDAnalysis.analysis import (align, diffusionmap, dihedrals, gnm,
                                 helix_analysis, rms)
from tqdm import tqdm


def smv(arr: np.ndarray, w: int = 50):
    return np.convolve(arr, np.ones(w), "valid") / w


def make_args(args: argparse.ArgumentParser):
    args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                        help="dir which contains many sub-folders, each sub-folder contains docking result files. Default is %(default)s.")
    args.add_argument('-top', '--top-name', type = str, default='md.tpr',
                        help='topology file name in each sub-directory, such as md.tpr. Default is %(default)s.')
    args.add_argument('-traj', '--traj-name', type = str, default='md_center.xtc',
                        help='trajectory file name in each sub-directory, such as md_center.xtc. Default is %(default)s.')
    args.add_argument('-b', '--begin-frame', type=int, default=0,
                        help='First frame to start the analysis. Default is %(default)s.')
    args.add_argument('-e', '--end-frame', type=int, default=None,
                        help='First frame to start the analysis. Default is %(default)s.')
    args.add_argument('-step', '--traj-step', type=int, default=1,
                        help='Step while reading trajectory. Default is %(default)s.')
    args.add_argument('-F', '--force', default=False, action='store_true',
                        help='force to re-run the analysis, default is %(default)s.')
    return args


class elastic(mmpbsa):
    HELP = """
    simple analysis collections for MDAnalysis.
    input is centered trajectory file, such as md_center.xtc.
    
    elastic network, using a Gaussian network model with only close contacts,
    analyzed the thermal fluctuation behavior of the system
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('-sele', '--select', type = str, default='protein and name CA',
                            help='selection for analysis.')
        args.add_argument('-close', '--close', action='store_true', default=False,
                            help='Using a Gaussian network model with only close contacts.')
        args.add_argument('-fast', '--fast', action='store_true', default=False,
                            help='Using fast computation implement.')
        args.add_argument('--cutoff', type = float, default=7,
                            help='elastic network neighber cutoff, default is %(default)s.')
        args.add_argument('--backend', type = str, default='numpy', choices=['numpy', 'torch', 'cuda'],
                            help='calculation backend, default is %(default)s.')
        args.add_argument('--block-size', type = int, default=100,
                            help='calculation block size, default is %(default)s.')
        args.add_argument('-np', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        return args
    
    def fast_calcu(self, u: mda.Universe, args: argparse.ArgumentParser):
        start, step, stop = args.begin_frame, args.traj_step, args.end_frame
        pool = TaskPool('process', args.n_workers).start()
        put_log(f'Calculating elastic network using lazydock.gmx.mda.gnm.calcu_{"closeContact" if args.close else ""}GNMAnalysis')
        ag, v, t = u.select_atoms(args.select), [], []
        sum_frames = (len(u.trajectory) if stop is None else stop) - start
        for frame in tqdm(u.trajectory[start:stop:step], total=sum_frames//step, desc='Calculating frames', leave=False):
            t.append(frame.time)
            if args.close:
                atom2res, res_size = genarate_atom2residue(ag)
                pool.add_task(frame.time, calcu_closeContactGNMAnalysis,
                              ag.positions.copy(), args.cutoff,
                              atom2res, res_size, ag.n_residues, 'size')
            else:
                pool.add_task(frame.time, calcu_GNMAnalysis, ag.positions.copy(), args.cutoff)
            pool.wait_till(lambda : pool.count_waiting_tasks() == 0, wait_each_loop=0.001, update_result_queue=False)
        # gether results
        for t_i in t:
            v.append(pool.query_task(t_i, True, 999)[0])
        pool.close()
        return np.array(t), np.array(v)
    
    def calcu(self, u: mda.Universe, args: argparse.ArgumentParser):
        put_log(f'Calculating elastic network using MDAnalysis.analysis.gnm.{"closeContact" if args.close else ""}GNMAnalysis')
        if args.close:
            nma = gnm.closeContactGNMAnalysis(u, select=args.select, cutoff=args.cutoff, weights='size')
        else:
            nma = gnm.GNMAnalysis(u, select=args.select, cutoff=args.cutoff)
        nma.run(start=args.begin_frame, step=args.traj_step, stop=args.end_frame, verbose=True)
        t, v = np.array(nma.results['times']).copy(), np.array(nma.results['eigenvalues']).copy()
        del nma
        return t, v

    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        file_name = 'elastic' + ('_close' if args.close else '')
        if os.path.exists(w_dir / f'{file_name}.xlsx') and not args.force:
            return put_log('Elastic analysis already calculated, use -F to re-run.')
        if args.fast:
            t, v = self.fast_calcu(u, args)
        else:
            t, v = self.calcu(u, args)
        df = pd.DataFrame({'times': t, 'eigenvalues': v})
        df.to_excel(w_dir / f'{file_name}.xlsx', index=False)
        fig = plt.figure(figsize=(10, 8))
        sns.lineplot(data=df, x='times', y='eigenvalues', alpha=0.6, ax=fig.gca())
        sns.lineplot(x=df['times'][49:], y=smv(df['eigenvalues']), label='SMV', ax=fig.gca())
        save_show(w_dir / f'{file_name}.png', 600, show=False)
        plt.close(fig=fig)
        
    def main_process(self):
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # process each task
        bar = tqdm(total=len(self.tasks), desc='Calculating')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            u = mda.Universe(top_path, traj_path)
            wdir = Path(wdir).resolve()
            self.analysis(u, wdir, self.args)
            bar.update(1)
        bar.close()


class rmsd(elastic):
    HELP = """Pairwise RMSD of a trajectory to itself"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('-sele', '--select', type = str, default='protein and name CA',
                            help='selection for analysis.')
        args.add_argument('--backend', type = str, default='numpy', choices=['numpy', 'torch', 'cuda'],
                            help='calculation backend, default is %(default)s.')
        args.add_argument('--block-size', type = int, default=100,
                            help='calculation block size, default is %(default)s.')
        return args
        
    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.traj_step, args.end_frame
        if os.path.exists(w_dir / 'inter_frame_rmsd.pkl') and not force:
            return put_log('Inter-frame RMSD already calculated, use -F to re-run.')
        put_log('Aligning inter-frame RMSD')
        aligner = align.AlignTraj(u, u, select=args.select, in_memory=True).run(verbose=True, start=start, step=step, stop=stop)
        # calcu interaction for each frame
        ag, coords = u.select_atoms(args.select), []
        sum_frames = (len(u.trajectory) if stop is None else stop) - start
        for _ in tqdm(u.trajectory[start:stop:step], total=sum_frames//step, desc='Gathering coordinates', leave=False):
            coords.append(ag.positions.copy())
        coords = np.array(coords)
        matrix: np.ndarray = pairwise_rmsd(coords, block_size=args.block_size, backend=args.backend, verbose=True)
        np.savez_compressed(w_dir / 'inter_frame_rmsd.npz', matrix=matrix.astype(np.float16))
        plt.imshow(matrix, cmap='viridis')
        plt.xlabel('Frame')
        plt.ylabel('Frame')
        plt.colorbar(label=r'RMSD ($\AA$)')
        save_show(w_dir / 'inter_frame_rmsd.png', 600, show=False)
        plt.close()
        del aligner, matrix
        

class rama(elastic):
    HELP = """plot Ramachandran and Janin"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('-rstep', '--rama-step', type=int, default=100,
                          help='Step while reading trajectory for plotting amachandran plot and Janin plot. Default is %(default)s.')
        args.add_argument('-alpha', '--alpha', type=float, default=0.2,
                          help='Scatter alpha for plotting amachandran plot and Janin plot. Default is %(default)s.')
        args.add_argument('-size', '--size', type=float, default=80,
                          help='Scatter size for plotting amachandran plot and Janin plot. Default is %(default)s.')

    def ramachandran(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.rama_step, args.end_frame
        if os.path.exists(w_dir / 'ramachandran.npz') and not force:
            return put_log('Ramachandran plot already calculated, use -F to re-run.')
        protein = u.select_atoms('protein')
        rama = dihedrals.Ramachandran(protein).run(start=start, step=step, stop=stop)
        np.savez_compressed(w_dir / 'ramachandran.npz', angles = rama.results.angles)
        fig, ax = plt.subplots(figsize=(10, 8))
        rama.plot(color='black', marker='.', ref=True, ax=ax, alpha=args.alpha, s=args.size)
        save_show(w_dir / 'ramachandran.png', 600, show=False)
        plt.close(fig=fig)

    def janin(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.rama_step, args.end_frame
        if os.path.exists(w_dir / 'janin.npz') and not force:
            return put_log('Janin plot already calculated, use -F to re-run.')
        protein = u.select_atoms('protein')
        janin = dihedrals.Janin(protein).run(start=start, step=step, stop=stop)
        np.savez_compressed(w_dir / 'janin.npz', angles = janin.results.angles)
        fig, ax = plt.subplots(figsize=(10, 8))
        janin.plot(color='black', marker='.', ref=True, ax=ax, alpha=args.alpha, s=args.size)
        save_show(w_dir / 'janin.png', 600, show=False)
        plt.close(fig=fig)

    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        self.ramachandran(u, w_dir, args)
        self.janin(u, w_dir, args)


class sele_ana(elastic):
    HELP = """
    sele analysis collections for MDAnalysis.
    input is centered trajectory file, such as md_center.xtc.
    
    1. average twist of the helix
    2. Radial distribution function of specific residue(s)
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('--twist-select', type = str, nargs='+', default = [],
                          help='twist select group, can be multiple MDAnalysis selection string, default is %(default)s.')
        args.add_argument('--radial-select', type = str, nargs='+', default = [],
                          help='radial select group, can be multiple MDAnalysis selection string, default is %(default)s.')
        args.add_argument('-np', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')

    def twist(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        twist_select, nw, start, step, stop = args.twist_select, args.n_workers, args.begin_frame, args.traj_step, args.end_frame
        dfs = {}
        for i, twist_sel in enumerate(twist_select):
            h = helix_analysis.HELANAL(u, select=twist_sel,
                                    ref_axis=[0, 0, 1]).run(n_workers=nw, backend='multiprocessing', start=start, step=step, stop=stop)
            dfs[twist_sel] = h.results.local_twists
            twi_i = dfs[twist_sel].mean(axis=1)
            fig = plt.figure(figsize=(10, 8))
            sns.lineplot(x=np.arange(len(twi_i)), y=twi_i, alpha=0.6, ax=fig.gca())
            sns.lineplot(x=np.arange(50, len(twi_i-50)), y=smv(twi_i), label='SMV', ax=fig.gca())
            save_show(w_dir / f'twist_{i}.png', 600, show=False)
            plt.close(fig=fig)
        opts_file(w_dir / 'twist.pkl', 'wb', way='pkl', data=dfs)

    def radial_dist(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        radial_select, nw, start, step, stop = args.radial_select, args.n_workers, args.begin_frame, args.traj_step, args.end_frame
        dfs = {}
        for i, radial_sel in enumerate(radial_select):
            r = helix_analysis.RadialDistribution(u, select=radial_sel, bins=100).run(n_workers=nw, backend='multiprocessing', start=start, step=step, stop=stop)
            dfs[radial_sel] = r.results.radial_distribution
            fig = plt.figure(figsize=(10, 8))
            sns.lineplot(x=r.results.bins, y=r.results.radial_distribution, alpha=0.6, ax=fig.gca())
            save_show(w_dir / f'radial_{i}.png', 600, show=False)
            plt.close(fig=fig)
        opts_file(w_dir / 'radial.pkl', 'wb', way='pkl', data=dfs)

    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        self.twist(u, w_dir, args)
        self.radial_dist(u, w_dir, args)
        
        
class pair_rmsd(elastic):
    HELP = """
    sele analysis collections for MDAnalysis.
    input is centered trajectory file, such as md_center.xtc.
    
    Pairwise RMSD between two trajectories
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-rd', '--ref-dir', type = str, required=True,
                          help="dir which contains ref input files, default is %(default)s.")
        make_args(args)
        args.add_argument('--ref-chain-name', type = str, required=True,
                          help='receptor chain name, such as "A".')
        args.add_argument('--chain-name', type = str, required=True,
                          help='receptor chain name, such as "A".')

    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        self.args.ref_dir = clean_path(self.args.ref_dir)

    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.traj_step, args.end_frame
        if os.path.exists(w_dir / 'inter_frame_rmsd.pkl') and not force:
            return put_log('Inter-frame RMSD already calculated, use -F to re-run.')
        put_log('Aligning inter-frame RMSD')
        ag, coords = u.select_atoms('name CA'), []
        aligner = align.AlignTraj(ag, self.ref_ag, select='name CA', in_memory=True).run(verbose=True, start=start, step=step, stop=stop)
        # calcu interaction for each frame
        sum_frames = (len(u.trajectory) if stop is None else stop) - start
        for _ in tqdm(u.trajectory[start:stop:step], total=sum_frames//step, desc='Gathering coordinates', leave=False):
            coords.append(ag.positions.copy())
        coords = np.array(coords)
        matrix = pairwise_rmsd(coords, block_size=50, backend='cuda', verbose=True)
        np.savez_compressed(w_dir / 'inter_frame_rmsd.npz', matrix=matrix)
        plt.imshow(matrix, cmap='viridis')
        plt.xlabel('Frame')
        plt.ylabel('Frame')
        plt.colorbar(label=r'RMSD ($\AA$)')
        save_show(w_dir / 'inter_frame_rmsd.png', 600, show=False)
        plt.close()
        del aligner, matrix
        
    def main_process(self):
        self.ref_top_path, self.ref_traj_paths = self.check_top_traj(bdir = self.args.ref_dir)
        self.ref_u = mda.Universe(self.ref_top_path, self.ref_traj_path, in_memory=True)
        self.ref_ag = self.ref_u.atoms[self.ref_u.atoms.chainIDs == self.args.ref_chain_name]
        super().main_process()


_str2func = {
    'elastic': elastic,
    'rmsd': rmsd,
    'rama': rama,
    'sele': sele_ana,
    'pari-rmsd': pair_rmsd,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for MDAnalysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    main()