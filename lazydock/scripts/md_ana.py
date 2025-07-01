'''
Date: 2025-01-16 10:08:37
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-06-02 12:21:20
Description: 
'''
import argparse
import os
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Union

from Bio import BiopythonDeprecationWarning
from matplotlib.ticker import FuncFormatter

warnings.simplefilter('ignore', BiopythonDeprecationWarning)
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
import seaborn as sns
from lazydock.algorithm.rms import batch_fit_to, batch_rmsd, pairwise_rmsd
from lazydock.gmx.mda.align import get_aligned_coords
from lazydock.gmx.mda.gnm import (calcu_closeContactGNMAnalysis,
                                  calcu_GNMAnalysis, genarate_atom2residue)
from lazydock.gmx.mda.utils import filter_atoms_by_chains
from lazydock.scripts._script_utils_ import (clean_path, excute_command,
                                             process_batch_dir_lst)
from lazydock.scripts.ana_gmx import mmpbsa
from matplotlib.cm import ScalarMappable
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import opts_file
from mbapy_lite.plot import save_show
from mbapy_lite.web_utils.task import TaskPool
from MDAnalysis.analysis import align, dihedrals, gnm, helix_analysis
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
        self.pool: TaskPool = None
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('-sele', '--select', type = str, default='protein and name CA',
                            help='selection for analysis.')
        args.add_argument('-c', '--chains', type = str, nargs='+', default=None,
                          help='chain of molecular to be included into calculation. Default is %(default)s.')
        args.add_argument('-close', '--close', action='store_true', default=False,
                            help='Using a Gaussian network model with only close contacts.')
        args.add_argument('-fast', '--fast', action='store_true', default=False,
                            help='Using fast computation implement.')
        args.add_argument('--cutoff', type = float, default=7,
                            help='elastic network neighber cutoff, default is %(default)s.')
        args.add_argument('--backend', type = str, default='numpy', choices=['numpy', 'torch', 'cuda'],
                            help='calculation backend, default is %(default)s.')
        args.add_argument('-nw', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        return args
    
    def fast_calcu(self, u: mda.Universe, args: argparse.ArgumentParser):
        start, step, stop = args.begin_frame, args.traj_step, args.end_frame
        ag, v, t = u.select_atoms(args.select), [], []
        if args.chains is not None:
            ag = filter_atoms_by_chains(ag, args.chains)
        sum_frames = (len(u.trajectory) if stop is None else stop) - start
        for frame in tqdm(u.trajectory[start:stop:step], total=sum_frames//step, desc='Calculating frames', leave=False):
            t.append(frame.time)
            if args.backend in {'torch', 'cuda'}:
                import torch
                device = 'cuda' if args.backend == 'cuda' else 'cpu'
                positions = torch.tensor(ag.positions, dtype=torch.float64, device=device)
            else:
                positions = ag.positions.copy()
            if args.close:
                atom2res, res_size = genarate_atom2residue(ag)
                if args.backend in {'torch', 'cuda'}:
                    atom2res = torch.tensor(atom2res, dtype=torch.int64, device=device)
                    res_size = torch.tensor(res_size, dtype=torch.float64, device=device)
                self.pool.add_task(frame.time, calcu_closeContactGNMAnalysis,
                                    positions, args.cutoff,
                                    atom2res, res_size, ag.n_residues, 'size', backend=args.backend)
            else:
                self.pool.add_task(frame.time, calcu_GNMAnalysis, positions, args.cutoff, backend=args.backend)
            self.pool.wait_till(lambda : self.pool.count_waiting_tasks() == 0, wait_each_loop=0.001, update_result_queue=False)
        # gether results
        for t_i in t:
            v.append(self.pool.query_task(t_i, True, 999)[0])
        return np.array(t), np.array(v)
    
    def calcu(self, u: mda.Universe, args: argparse.ArgumentParser):
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
        if self.__class__.__name__ == 'elastic':
            self.pool = TaskPool('process', self.args.n_workers).start()
        # process each task
        bar = tqdm(total=len(self.tasks), desc='Calculating')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir.replace(str(self.args.batch_dir), '')}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            u = mda.Universe(top_path, traj_path)
            wdir = Path(wdir).resolve()
            self.analysis(u, wdir, self.args)
            bar.update(1)
        if self.pool:
            self.pool.close(1)
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
        args.add_argument('-f', '-fit', '--fit', action='store_true', default=False,
                            help='Using pair-wise centering and comformation fitting before the RMSD, else use batch_fit_to to align mobile to the ref\'s first frame.')
        args.add_argument('--backend', type = str, default='numpy', choices=['numpy', 'torch', 'cuda'],
                            help='calculation backend, default is %(default)s.')
        args.add_argument('--block-size', type = int, default=100,
                            help='calculation block size, default is %(default)s.')
        return args
        
    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.traj_step, args.end_frame
        if os.path.exists(w_dir / 'inter_frame_rmsd.pkl') and not force:
            return put_log('Inter-frame RMSD already calculated, use -F to re-run.')
        # calcu interaction for each frame
        ag, coords = u.select_atoms(args.select), []
        sum_frames = (len(u.trajectory) if stop is None else stop) - start
        for _ in tqdm(u.trajectory[start:stop:step], total=sum_frames//step, desc='Gathering coordinates', leave=False):
            coords.append(ag.positions.copy().astype(np.float64))
        coords = np.array(coords)
        if args.backend != 'numpy':
            import torch
            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
            coords = torch.tensor(coords, device=device, dtype=torch.float64)
        # do center as MDAnalysis.align.AlignTraj._single_frame
        coords: np.ndarray = coords - coords.mean(axis=1, keepdims=True)
        # do fit and calcu RMSD depend on backend
        if args.fit:
            # NOTE: if the frames already aligned to the first frame,
            #       the RMSD between the i-j will also like pair-wise-align RMSD.
            #       Because BOTH are aligned to the same first frame, like it aligned to each other!
            #       So the result wit and without fit will be the same!
            matrix = np.zeros((coords.shape[0], coords.shape[0])) if args.backend == 'numpy' else torch.zeros((coords.shape[0], coords.shape[0]), dtype=torch.float64, device=device)
            # calcu rmsd for each block
            K = int(np.ceil(coords.shape[0] / args.block_size))
            for i in tqdm(range(K), total=K, desc='Calculating RMSD matrix', leave=False):
                # prepare block i
                start_i = i * args.block_size
                end_i = (i + 1) * args.block_size if i < K - 1 else coords.shape[0]
                # do center as MDAnalysis.align.AlignTraj._prepare, but coords[0] already centered
                if args.backend == 'numpy':
                    refs = np.repeat(coords[start_i:end_i, :, :], coords.shape[0], axis=0)
                    mobile = np.repeat(coords[:, :, :], args.block_size, axis=0)
                else:
                    refs = coords[start_i:end_i, :, :].repeat(coords.shape[0], 1, 1)
                    mobile = coords[:, :, :].repeat(args.block_size, 1, 1)
                # do fit as MDAnalysis.align.AlignTraj._single_frame, it use _fit_to
                matrix[start_i:end_i, :] = batch_fit_to(mobile, refs, backend=args.backend)[1].reshape(args.block_size, coords.shape[0])
            if args.backend != 'numpy':
                matrix = matrix.cpu().numpy()
        else:
            # do center as MDAnalysis.align.AlignTraj._prepare, but coords[0] already centered
            if args.backend == 'numpy':
                refs = np.repeat(coords[0][None, :, :], coords.shape[0], axis=0)
            else:
                refs = coords[0][None, :, :].repeat(coords.shape[0], 1, 1)
            # do fit as MDAnalysis.align.AlignTraj._single_frame, it use _fit_to
            coords = batch_fit_to(coords, refs, backend=args.backend)[0]
            matrix: np.ndarray = pairwise_rmsd(coords, block_size=args.block_size, backend=args.backend, verbose=True)
        np.savez_compressed(w_dir / 'inter_frame_rmsd.npz', matrix=matrix.astype(np.float16))
        plt.imshow(matrix, cmap='viridis')
        plt.xlabel('Frame', fontsize=16, weight='bold')
        plt.ylabel('Frame', fontsize=16, weight='bold')
        plt.gca().tick_params(labelsize=14)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('RMSD ($\AA$)', fontsize=16, weight='bold')
        cbar.ax.tick_params(labelsize=14)
        save_show(w_dir / 'inter_frame_rmsd.png', 600, show=False)
        plt.close()
        

class rama(elastic):
    HELP = """plot Ramachandran and Janin"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('-rstep', '--rama-step', type=int, default=100,
                          help='Step while reading trajectory for plotting amachandran plot and Janin plot. Default is %(default)s.')
        args.add_argument('-np', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        args.add_argument('-alpha', '--alpha', type=float, default=0.2,
                          help='Scatter alpha for plotting amachandran plot and Janin plot. Default is %(default)s.')
        args.add_argument('-size', '--size', type=float, default=80,
                          help='Scatter size for plotting amachandran plot and Janin plot. Default is %(default)s.')
        args.add_argument('-c', '--chains', type = str, nargs='+', default=None,
                          help='chain of molecular to be included into calculation. Default is %(default)s.')
        
    def plot_matrix(self, result, file_name: str, w_dir: Path, args: argparse.ArgumentParser):
        fig, ax = plt.subplots(figsize=(10, 8))
        result.results.angles = result.results.angles[::args.rama_step]
        result.plot(color='black', marker='.', ref=True, ax=ax, alpha=args.alpha, s=args.size)
        save_show(w_dir / file_name, 600, show=False)
        plt.close(fig=fig)

    def ramachandran(self, ag: mda.AtomGroup, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.traj_step, args.end_frame
        if os.path.exists(w_dir / 'ramachandran.npz') and not force:
            return put_log('Ramachandran plot already calculated, use -F to re-run.')
        rama = dihedrals.Ramachandran(ag).run(start=start, step=step, stop=stop,
                                              n_workers=args.n_workers, backend='multiprocessing')
        np.savez_compressed(w_dir / 'ramachandran.npz', angles = rama.results.angles)
        self.plot_matrix(rama, 'ramachandran.png', w_dir, args)

    def janin(self, ag: mda.AtomGroup, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.traj_step, args.end_frame
        if os.path.exists(w_dir / 'janin.npz') and not force:
            return put_log('Janin plot already calculated, use -F to re-run.')
        janin = dihedrals.Janin(ag).run(start=start, step=step, stop=stop,
                                        n_workers=args.n_workers, backend='multiprocessing')
        np.savez_compressed(w_dir / 'janin.npz', angles = janin.results.angles)
        self.plot_matrix(janin, 'janin.png', w_dir, args)

    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        ag = u.select_atoms('protein')
        if self.args.chains is not None:
            ag = filter_atoms_by_chains(ag, self.args.chains)
        put_log(f'get {len(ag)} atoms for analysis.')
        self.ramachandran(ag, w_dir, args)
        self.janin(ag, w_dir, args)


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
        
        
class pca(elastic):
    HELP = """
    Principal component analysis of the trajectory
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        make_args(args)
        args.add_argument('-c', '--chains', type = str, nargs='+', default=None,
                          help='chain of molecular to be included into calculation. Default is %(default)s.')
        args.add_argument('-sele', '--select', type = str, default='backbone',
                            help='selection for analysis.')
        args.add_argument('--backend', type = str, default='numpy', choices=['numpy', 'torch', 'cuda'],
                            help='calculation backend, default is %(default)s.')
        args.add_argument('--n-components', type = int, default=None,
                            help='number of components to calculate, default is %(default)s.')
        return args
    
    def plot_PCA(self, df: pd.DataFrame, w_dir: Path, args: argparse.ArgumentParser):
        # 配置参数
        pc_columns = ['PC1', 'PC2', 'PC3']
        n = len(pc_columns)
        time = df['Time (ps)'].values
        # 创建3x3子图
        fig, axes = plt.subplots(n, n, figsize=(12, 12), 
                                gridspec_kw={'wspace':0.05, 'hspace':0.05})
        # 颜色条配置
        cmap = plt.get_cmap('viridis')
        sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(time.min(), time.max()))
        # 遍历所有子图组合
        for i in range(n):
            for j in range(n):
                ax = axes[i, j]
                if i <=1:
                    ax.set_xticks([])
                if j > 0:
                    ax.set_yticks([])
                # 绘制散点（核心优化参数）
                sc = ax.scatter(df[pc_columns[j]], df[pc_columns[i]],
                                c=time, cmap=cmap, marker='.', alpha=0.8, linewidths=0,
                                rasterized=True) # 启用栅格化加速
                # 添加坐标标签
                if i == n-1:
                    ax.set_xlabel(pc_columns[j], labelpad=5, fontsize=16, weight='bold')
                if j == 0:
                    ax.set_ylabel(pc_columns[i], labelpad=5, fontsize=16, weight='bold')
                ax.tick_params(labelsize=14)
        # 添加颜色条
        cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(sm, cax=cax, label='Time (ps)')
        cax.tick_params(labelsize=14)
        cax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x/1000:.0f}'))
        cax.set_ylabel('Time (ns)', fontsize=16)
        save_show(w_dir / 'PCA.png', 600, show=False)
        plt.close()
        
    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        force, start, step, stop = args.force, args.begin_frame, args.traj_step, args.end_frame
        if os.path.exists(w_dir / 'PCA.csv') and not force:
            return put_log('PCA already calculated, use -F to re-run.')
        # determine backend
        if args.backend != 'numpy':
            import torch as _backend
        else:
            _backend = np
        # get aligned coords: [n_frames, n_atoms, 3]
        ag = u.select_atoms(args.select)
        if args.chains is not None:
            ag = filter_atoms_by_chains(ag, args.chains)
        ori_coords, coords = get_aligned_coords(u, ag, start, step, stop, backend=args.backend, verbose=True)
        # calculate correlation matrix
        ## Center the data
        if args.backend != 'numpy':
            coords -= _backend.mean(ori_coords, dim=0, keepdim=True)
        else:
            coords -= _backend.mean(ori_coords, axis=0, keepdims=True)
        # eigenvalue decomposition
        ## Reshape to [n_frames, 3N] and calculate covariance
        coords = coords.reshape(coords.shape[0], -1)
        if args.backend != 'numpy':
            self.cov = _backend.cov(coords.T)
            e_vals, e_vects = _backend.linalg.eigh(self.cov)
            sort_idx = _backend.argsort(e_vals, descending=True)
        else:
            self.cov = _backend.cov(coords, rowvar=False)
            e_vals, e_vects = _backend.linalg.eig(self.cov)
            sort_idx = _backend.argsort(e_vals)[::-1]
        self.variance = e_vals[sort_idx]
        self._p_components = e_vects[:, sort_idx]
        # set n components
        if args.n_components is None:
            args.n_components = len(self.variance)
        self.variance = self.variance[:args.n_components]
        self.cumulated_variance = (_backend.cumsum(self.variance, 0) / _backend.sum(self.variance))[:args.n_components]
        self.p_components = self._p_components[:, :args.n_components]
        # do the transform
        dot = coords @ self.p_components[:, :3]
        # save result as dataframe
        if args.backend != 'numpy':
            dot = dot.cpu().numpy()
        df = pd.DataFrame(dot, columns=['PC{}'.format(i+1) for i in range(3)])
        df['Time (ps)'] = df.index * u.trajectory.dt
        df.to_csv(w_dir / 'PCA.csv', index=False)
        # plot result
        self.plot_PCA(df, w_dir, args)
        
        
class show_chain(elastic):
    HELP = """
    PPrint chain info of a tpr file
    """
    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        print() # new line after tqdm
        for i, chain_i in enumerate(np.unique(u.atoms.chainIDs)):
            idx = u.atoms.chainIDs == chain_i
            print(f'chain {i}: {chain_i}: {len(u.atoms[idx])} atoms')    


_str2func = {
    'elastic': elastic,
    'rmsd': rmsd,
    'rama': rama,
    'sele': sele_ana,
    'pari-rmsd': pair_rmsd,
    'pca': pca,
    'show-chain': show_chain,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for MDAnalysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    main()