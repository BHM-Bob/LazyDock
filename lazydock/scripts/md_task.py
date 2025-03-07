'''
Date: 2025-02-01 11:07:08
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-07 16:20:42
Description: 
'''
import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

import MDAnalysis as mda
import networkx as nx
import numpy as np
import pandas as pd
from lazydock.scripts._script_utils_ import excute_command
from lazydock.scripts.ana_gmx import mmpbsa
from lazydock_md_task.scripts.calc_correlation import plot_map
from lazydock_md_task.scripts.contact_map_v2 import (calculate_contacts,
                                                     load_and_preprocess_traj,
                                                     plot_network,
                                                     save_network_data)
from lazydock_md_task.scripts.prs_v2 import main as prs_main
from matplotlib import pyplot as plt
from mbapy.web_utils.task import TaskPool
from mbapy_lite.base import put_err, put_log
from mbapy_lite.plot import save_show
from tqdm import tqdm


def construct_graph(frame, atoms_inices: np.ndarray, threshold=6.7):
    # calculate distance matrix and build edges
    nodes = range(0, len(atoms_inices))
    dist = np.sqrt(np.sum((frame.positions[atoms_inices][:, None] - frame.positions[atoms_inices][None, :])**2, axis=-1))
    indices = np.where(dist < threshold)
    edges = list(filter(lambda x: x[0] < x[1], zip(indices[0], indices[1])))
    # build networkx graph
    protein_graph = nx.Graph()
    protein_graph.add_nodes_from(nodes)
    protein_graph.add_edges_from(edges)
    return protein_graph

def calcu_nextwork_from_frame(g):
    try:
        bc = nx.betweenness_centrality(g, normalized=False)
        bc = np.asarray(list(bc.values())).reshape(-1)
        path_dict = dict(nx.all_pairs_shortest_path_length(g))
        path_df = pd.DataFrame(path_dict)
        path_df.sort_index(axis=0, inplace=True)
        path_df.sort_index(axis=1, inplace=True)
        return bc, path_df.values
    except Exception as ex:
        return put_err(f'Error calculating BC for frame: {ex}')


class network(mmpbsa):
    HELP = """
    network analysis for MD-TASK
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.result_suffix = 'network.npz'
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('-top', '--top-name', type = str, default='md.tpr',
                            help='topology file name in each sub-directory, such as md.tpr. Default is %(default)s.')
        args.add_argument('-traj', '--traj-name', type = str, default='md_center.xtc',
                            help='trajectory file name in each sub-directory, such as md_center.xtc. Default is %(default)s.')
        args.add_argument('-c', '--chain', type = str, default=None,
                          help='chain of molecular to be included into calculation. Default is %(default)s.')
        args.add_argument("--threshold", type=float, default=6.7,
                          help="Maximum distance threshold in Angstroms when constructing graph (default: %(default)s)")
        args.add_argument('-b', '--begin-frame', type=int, default=0,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-e', '--end-frame', type=int, default=None,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-step', '--traj-step', type=int, default=1,
                          help='Step while reading trajectory. Default is %(default)s.')
        args.add_argument('-np', "--n-workers", type=int, default=4,
                          help="Number of workers to parallelize the calculation, default: %(default)s")
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='force to re-run the analysis, default is %(default)s.')
        return args
    
    def calcu_network(self, topol_path: Path, traj_path: Path):
        # prepare trajectory and topology
        u = mda.Universe(str(topol_path), str(traj_path))
        atoms = u.select_atoms("(name CB and protein) or (name CA and resname GLY and protein)")
        if self.args.chain is not None:
            atoms = atoms[atoms.chainIDs == self.args.chain]
        # prepare and run parallel calculation
        sum_frames = (len(u.trajectory) if self.args.end_frame is None else self.args.end_frame) - self.args.begin_frame
        total_bc, total_dj_path = [None] * sum_frames, [None] * sum_frames
        for current, frame in enumerate(tqdm(u.trajectory[self.args.begin_frame:self.args.end_frame:self.args.traj_step],
                                             desc='Calculating network', total=sum_frames, leave=False)):
            pg = construct_graph(frame, atoms.indices, self.args.threshold)
            self.pool.add_task(current, calcu_nextwork_from_frame, pg)
            self.pool.wait_till(lambda: self.pool.count_waiting_tasks() == 0, 0.001, update_result_queue=False)
        # gather results
        self.pool.wait_till(lambda: self.pool.count_done_tasks() == sum_frames, 0)
        for i in tqdm(list(self.pool.tasks.keys()), total=sum_frames, desc='Gathering results', leave=False):
            total_bc[i], total_dj_path[i] = self.pool.query_task(i, True, 999)
        self.pool.clear()
        return np.asarray(total_bc), np.asarray(total_dj_path)
        
    def save_results(self, top_path: Path, total_bc: np.ndarray, total_dj_path: np.ndarray):
        total_bc, total_dj_path = total_bc.mean(axis=0), total_dj_path.mean(axis=0)
        # plot Average Shortest Path figure
        plt.plot(total_bc)
        plt.xlabel('Residue (aa)', fontsize=16, weight='bold')
        plt.ylabel('Betweenness Centrality', fontsize=16, weight='bold')
        save_show(top_path.parent / f'{top_path.stem}_Betweenness Centrality.png', 600, show=False)
        plt.close()
        # plot Average Shortest Path figure
        plt.imshow(total_dj_path, cmap='viridis')
        plt.xlabel('Residue (aa)', fontsize=16, weight='bold')
        plt.ylabel('Residue (aa)', fontsize=16, weight='bold')
        plt.colorbar(label=r'Average Shortest Path')
        save_show(top_path.parent / f'{top_path.stem}_Average Shortest Path.png', 600, show=False)
        plt.close()
        # save outputs, float32 and int16 for saving space
        np.savez(top_path.parent / f'{top_path.stem}_network.npz', total_bc=total_bc.astype(np.float32),
                 total_dj_path=total_dj_path.astype(np.int16))
    
    def main_process(self):
        self.pool = TaskPool('process', self.args.n_workers).start()
        # load origin dfs from data file
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # process each task
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            top_path, traj_path = Path(top_path), Path(traj_path)
            if os.path.exists(os.path.join(wdir, f'{top_path.stem}_{self.result_suffix}')) and not self.args.force:
                put_log(f'{top_path.stem}_{self.result_suffix} already exists, skip.')
                continue
            results = self.calcu_network(Path(top_path), Path(traj_path))
            self.save_results(top_path, *results)
            bar.update(1)
        self.pool.close()
        
        
class correlation(network):
    HELP = """
    correlation analysis for MD-TASK
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.result_suffix = 'corr_matrix.npz'
        
    def correlate(self, coords):
        # residues shape: (n_traj_frame, n_res, 3)
        n_traj, n_res, _ = coords.shape
        # centerlize coords for each residue
        mean = np.mean(coords, axis=0, keepdims=True)
        delta = coords - mean  # (n_traj, n_res, 3)
        # calculate magnitude for each residue
        dot_products = np.sum(delta ** 2, axis=2)  # (n_traj, n_res, 3) => (n_traj, n_res)
        mean_dots = np.mean(dot_products, axis=0)  # (n_traj, n_res) => (n_res,)
        magnitudes = np.sqrt(mean_dots)  # (n_res,)
        
        # 计算所有残基对的协方差矩阵
        # 使用einsum计算每个残基对在所有时间点的点积之和，再除以时间数得到平均值
        cov_matrix = np.einsum('tix,tjx->ij', delta, delta) / n_traj
        # calculate correlation matrix
        corr_matrix = cov_matrix / (magnitudes[:, None] * magnitudes[None, :])
        
        return corr_matrix
    
    def save_results(self, top_path: Path, corr_matrix: np.ndarray):
        np.savez(top_path.parent / f'{top_path.stem}_corr_matrix.npz', corr_matrix=corr_matrix.astype(np.float32))
        plot_map(corr_matrix, top_path.stem, top_path.parent / f'{top_path.stem}_corr_matrix')
    
    def calcu_network(self, topol_path: Path, traj_path: Path):
        # prepare trajectory and topology
        u = mda.Universe(str(topol_path), str(traj_path))
        atoms = u.select_atoms("(name CA and protein)")
        if self.args.chain is not None:
            atoms = atoms[atoms.chainIDs == self.args.chain]
        # extract coords
        sum_frames = (len(u.trajectory) if self.args.end_frame is None else self.args.end_frame) - self.args.begin_frame
        coords = np.zeros((len(u.trajectory), len(atoms), 3), dtype=np.float64)
        for current, _ in enumerate(tqdm(u.trajectory[self.args.begin_frame:self.args.end_frame:self.args.traj_step],
                                         desc='Gathering coordinates', total=sum_frames, leave=False)):
            coords[current] = atoms.positions
        # calculate correlation matrix and save, show
        sorted_residx = np.argsort(atoms.resids)
        corr_matrix = self.correlate(coords[:, sorted_residx, :])
        return [corr_matrix]
     
        
class prs(network):
    HELP = """
    prs analysis for MD-TASK
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.result_suffix = 'PRS.png'
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        network.make_args(args)
        args.add_argument("--perturbations", type=int, default=250,
                          help="number of perturbations, default: %(default)s.")

    def calcu_network(self, topol_path: Path, traj_path: Path):
        return [prs_main(top_path=str(topol_path), traj_path=str(traj_path), chains=[self.args.chain],
                         start=self.args.begin_frame, stop=self.args.end_frame, step=self.args.traj_step,
                         perturbations=self.args.perturbations, n_worker=self.args.n_workers)]
    
    def save_results(self, top_path: Path, max_RHO: np.ndarray):
        # save to csv
        pd.DataFrame({'Residue': np.arange(len(max_RHO)), 'max_RHO': max_RHO}).to_csv(top_path.parent / f'{top_path.stem}_PRS.csv', index=False)
        # plot figure
        plt.plot(max_RHO)
        plt.xlabel('Residue (aa)', fontsize=16, weight='bold')
        plt.ylabel('Correlation coefficient', fontsize=16, weight='bold')
        plt.gca().tick_params(labelsize=14, axis='both')
        save_show(top_path.parent / f'{top_path.stem}_PRS.png', 600, show=False)
        plt.close()



class contact_map(network):
    HELP = """
    contact map analysis for MD-TASK
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.result_suffix = 'contact_map.png'
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        network.make_args(args)
        args.add_argument("--residue", type=str, required=True, help="traget residue name, such as LYS111")
        args.add_argument("--nodesize", type=int, default=2900, help="node size in drawing, default: %(default)s.")
        args.add_argument("--nodefontsize", type=float, default=9.5, help="node font size in drawing, default: %(default)s.")
        args.add_argument("--edgewidthfactor", type=float, default=10.0, help="edge width factor in drawing, default: %(default)s.")
        args.add_argument("--edgelabelfontsize", type=float, default=8.0, help="edge label font size in drawing, default: %(default)s.")
        return args
    
    def process_args(self):
        super().process_args()
        self.args.residue = self.args.residue.upper()
        self.args.chain = self.args.chain.upper()
        self.prefix = self.args.residue.split(".")[0] if "." in self.args.residue else self.args.residue
        
    def calcu_network(self, traj_path: Path, topol_path: Path):
        # 1. load trajectory and topology
        traj = load_and_preprocess_traj(str(traj_path), str(topol_path), self.args.step)
        # 2. calculate contacts
        contacts, n_frames = calculate_contacts(
            traj, self.args.residue, self.args.chain, self.args.threshold/10
        )
        # 3. generate edges list
        center_node = f"{self.args.residue}.{self.args.chain}"
        edges_list = [[center_node, edge[1], count/n_frames] for edge, count in contacts.items()]
        # 4. create graph object
        contact_graph = nx.Graph()
        contact_graph.add_weighted_edges_from(edges_list)
        # 5. save output results
        output_csv = f"{self.prefix}_chain{self.args.chain}_network.csv"
        save_network_data(edges_list, output_csv)
        # 6. generate visualization graph
        output_png = f"{self.prefix}_chain{self.args.chain}_contact_map.png"
        plot_network(contact_graph, edges_list, output_png,
                     node_size=self.args.nodesize, node_fontsize=self.args.nodefontsize,
                     edgewidth_factor=self.args.edgewidthfactor, edgelabel_fontsize=self.args.edgelabelfontsize)


_str2func = {
    'network': network,
    'correlation': correlation,
    'prs': prs,
    'contact-map': contact_map,
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