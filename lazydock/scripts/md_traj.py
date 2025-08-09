'''
Date: 2025-02-27 22:08:05
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-07 19:02:11
Description: 
'''

import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Union

import matplotlib as mpl
import MDAnalysis as mda
from matplotlib.ticker import FuncFormatter
import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
from lazydock.scripts._script_utils_ import excute_command
from lazydock.scripts.ana_gmx import RRCS
from matplotlib import pyplot as plt
from mbapy_lite.stats import pca
from mbapy_lite.web_utils.task import TaskPool
from mbapy_lite.base import put_err, put_log
from mbapy_lite.plot import save_show
from tqdm import tqdm


def run_md_traj_simple_analysis(top_path: str, traj_path: str, args: argparse.Namespace):
    root = Path(top_path).resolve().parent
    main_name = Path(traj_path).stem
    if (root / f'{main_name}_nematic_order.csv').exists() and not args.force:
        return put_log(f'{root} already calculated, skip.')
    def save_df(data: np.ndarray, name: str, cols: str = None):
        df = pd.DataFrame(data, columns=cols or [name])
        df.to_csv(root / f'{main_name}_{name}.csv')
    # load trajectory
    u = mda.Universe(top_path, traj_path)
    # select atoms
    if args.chains is not None and len(args.chains):
        idx = u.atoms.chainIDs == args.chains[0]
        for chain_i in args.chains[1:]:
            idx = idx | (u.atoms.chainIDs == chain_i)
        ag = u.atoms[idx]
    else:
        ag = u.atoms
    # load trajectory
    t = md.load(traj_path, top=str(root / args.gro_name))[args.begin_frame:args.end_frame:args.traj_step]
    mol = t.atom_slice(ag.ids)
    # directors
    save_df(md.compute_directors(mol)[:, 0, :], 'directors', ['x', 'y', 'z'])
    # order
    save_df(md.compute_nematic_order(mol), 'nematic_order')
    # density
    save_df(md.density(mol), 'density')
    # isothermal_compressability_kappa_T
    compressability = []
    for i in range(0, len(u.trajectory), 100):
        compressability.append(md.isothermal_compressability_kappa_T(t[i*100: (i+1)*100], 300))
    save_df(compressability, 'isothermal_compressability_kappa_T')
    # interia
    inters = md.compute_inertia_tensor(mol)
    save_df(np.stack([inters[:, 0, 0], inters[:, 1, 1], inters[:, 2, 2]], axis=1), 'inertia_tensor', ['xx', 'yy', 'zz'])
    # kabsch_sander
    csr_matrixs = md.kabsch_sander(mol)
    np.savez_compressed(root / f'{main_name}_kabsch_sander.npz', csr_matrixs=csr_matrixs)
    mean_mat = np.zeros_like(csr_matrixs[0].toarray())
    frame_avg = np.zeros(len(mol.time))
    for i, mat_i in enumerate(csr_matrixs):
        mat_i = mat_i.toarray()
        mean_mat += mat_i
        frame_avg[i] = mat_i.mean().mean()
    mean_mat /= len(csr_matrixs)
    np.savez_compressed(root / f'{main_name}_kabsch_sander_mean_mat.npz', csr_matrixs=mean_mat)
    save_df(frame_avg, 'kabsch_sander_frame_avg')
    # shape
    save_df(md.principal_moments(mol), 'principal_moments', ['x', 'y', 'z'])
    save_df(md.asphericity(mol), 'asphericity')
    save_df(md.acylindricity(mol), 'acylindricity')
    save_df(md.relative_shape_antisotropy(mol), 'relative_shape_antisotropy')
    save_df(md.density(mol), 'density')
    # bond angle and dihedrals
    phi = md.compute_phi(mol, False)
    psi = md.compute_psi(mol, False)
    phi_, psi_ = pca(phi[1], 10), pca(psi[1], 10)
    scatter = plt.scatter(phi_[:, 0], psi_[:, 0], c=t.time, alpha=0.1, cmap='viridis')
    sm = mpl.cm.ScalarMappable(norm=scatter.norm, cmap=scatter.cmap)
    sm.set_array(t.time)
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x/1000:.0f}'))
    cbar.ax.tick_params(labelsize=14)
    cbar.ax.set_ylabel('Time (ns)', fontsize=16, weight='bold')
    plt.gca().tick_params(labelsize=14)
    plt.xlabel('PCA-0 of $\phi$', fontsize=16, weight='bold')
    plt.ylabel('PCA-0 of $\psi$', fontsize=16, weight='bold')
    save_show(str(root / f'{main_name}_phi_psi.png'), 600)
    plt.close()


class simple_analysis(RRCS):
    HELP = """
    network analysis for MD-TASK
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    def main_process(self):
        # load origin dfs from data file
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # process each complex
        pool, tasks = TaskPool('process', self.args.n_workers).start(), []
        for top_path, traj_path in tqdm(self.tasks, desc='process tasks'):
            # perform analysis
            tasks.append(pool.add_task(None, run_md_traj_simple_analysis, top_path, traj_path, self.args))
            pool.wait_till(lambda: pool.count_waiting_tasks() == 0, 0.01, update_result_queue=False)
        pool.wait_till_tasks_done(tasks)
        pool.close(1)
        

_str2func = {
    'simple-analysis': simple_analysis,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS analysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    # dev code
    
    main()