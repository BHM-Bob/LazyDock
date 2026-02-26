'''
Date: 2026-02-19
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2026-02-26 21:59:15
Description: steps most from http://www.mdtutorials.com/gmx/umbrella
'''
import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from MDAnalysis import Universe
from pymol import cmd
from tqdm import tqdm

from lazydock.gmx.prepare_ff import insert_content
from lazydock.gmx.run import Gromacs
from lazydock.scripts._script_utils_ import make_args_and_excute
from lazydock.scripts.run_gmx import simple_complex as _simple_complex
from lazydock.scripts.run_gmx import simple_protein as _simple_protein
from lazydock.utils import uuid4


class pull(_simple_protein):
    HELP = """
    run pull process for umbrella sampling, steps most from http://www.mdtutorials.com/gmx/umbrella
    Note: 
        1. pull-chain will be named as PULL_Chain in pull.ndx, and restrain-chain will be named as RESTRAIN_Chain.
        2. if use posres, you should add define=-DXXX in each mdp file(pull, nvt, npt, pull_md).
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.indexs = {}
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--protein-name', type = str, required=True,
                          help='protein name in each sub-directory, such as protein.gro.')
        args.add_argument('-pc', '--pull-chain', type=str, required=True,
                          help='chain name for pull, such as A.')
        args.add_argument('-rc', '--restrain-chain', type=str, required=True,
                          help='chain name for position restrain, such as B.')
        args.add_argument('-cca', '--calcu-center-atom', type=int, default=None,
                          help='add pull-group1-pbcatom = #number_of_your_central_atom in mdp file, will replace group idx with restrain group cca value, default is %(default)s.')
        args.add_argument('-posres', '--position-restrain', type=str, nargs='+', default=None,
                          help='add position restrain define in topol-itp file, format is CHAIN_ID RESTRAIN_DEF, example is A POSRES_RESTRAIN_CHAIN. Default is %(default)s.')
        args.add_argument('-po', '--pull-only', action='store_true', default=False,
                          help='only run pull process, default is %(default)s.')
        args.add_argument('-sd', '--start-distance', type=float, default=None,
                          help='start distance (nm) from pull trajectory for sampling MD run, default is %(default)s.')
        args.add_argument('-ed', '--end-distance', type=float, default=None,
                          help='end distance (nm) from pull trajectory for sampling MD run, default is %(default)s.')
        args.add_argument('-di', '--distance-interval', type=float, default=0.2,
                          help='distance (nm) interval from pull trajectory for sampling MD run, default is %(default)s.')
        args.add_argument('--pull-mdp', type = str, required=True,
                          help='pull mdp file, a file-path, will copy to working directory.')
        args.add_argument('--nvt-mdp', type = str, required=True,
                          help='nvt mdp file, if is a file-path, will copy to working directory.')
        args.add_argument('--npt-mdp', type = str, required=True,
                          help='npt mdp file, if is a file-path, will copy to working directory.')
        args.add_argument('--md-mdp', type = str, required=True,
                          help='production md mdp file, if is a file-path, will copy to working directory.')
        args.add_argument('--maxwarn', type=int, default=0,
                          help='maxwarn for em,nvt,npt,md gmx grompp command, default is %(default)s.')
        
    def process_args(self):
        super().process_args()
        if self.args.position_restrain is not None:
            if not len(self.args.position_restrain) % 2 == 0:
                put_err("position_restrain must be even number of arguments, input as CHAIN_ID_1 RESTRAIN_DEF_1 CHAIN_ID_2 RESTRAIN_DEF_2. exit")
                exit(1)
        
    def get_complex_atoms_index(self, u: Universe):
        pull_idx = u.atoms.chainIDs == self.args.pull_chain
        res_idx = u.atoms.chainIDs == self.args.restrain_chain
        put_log(f"pull atoms: {pull_idx.sum()}, restrain atoms: {res_idx.sum()}.")
        return pull_idx, res_idx
    
    def get_index_range(self, idx: np.ndarray):
        return idx.argmax(), idx.shape[0] - idx[::-1].argmax()
        
    def make_restrain_idx(self, protein_path: Path, gmx: Gromacs):
        # get resi range for restrain chain and pull chain
        u = Universe(str(protein_path.parent / 'npt.tpr'), str(protein_path.parent / 'npt.xtc'))
        pull_idx, res_idx = self.get_complex_atoms_index(u)
        pull_min, pull_max = self.get_index_range(pull_idx)
        res_min, res_max = self.get_index_range(res_idx)
        pull_range_str, res_range_str = f"{pull_min+1}-{pull_max}", f"{res_min+1}-{res_max}"
        # make restrain index file for pull and restrain chain
        groups = gmx.get_groups('npt.tpr')
        gmx.run_gmx_with_expect('make_ndx', f=str(protein_path.parent / 'npt.tpr'), o='pull.ndx',
                                expect_actions=[{'>': f'a {pull_range_str}\r'}, {'>': f'name {len(groups)} PULL_Chain\r'},
                                                {'>': f'a {res_range_str}\r'}, {'>': f'name {len(groups)+1} RESTRAIN_Chain\r'},
                                                {'>': 'q\r'}])
        
    def apply_restrain_def(self, gmx: Gromacs):
        """add position restrain define in topol_Protein_chain_A.itp file"""
        to_insert_content = '\n\n#ifdef POSERES_DEF\n#include "posre_Protein_chain_CHAIN_ID.itp"\n#endif\n\n'
        for chain_id, restrain_def in zip(self.args.position_restrain[::2], self.args.position_restrain[1::2]):
            itp_path = str(gmx.wdir / f'topol_Protein_chain_{chain_id}.itp')
            if os.path.exists(itp_path):
                insert_content(itp_path, '#endif', to_insert_content.replace('CHAIN_ID', chain_id).replace('POSERES_DEF', restrain_def))
        
    def calcu_com_com_distance(self, u: Universe):
        res_idx, pull_idx = self.get_complex_atoms_index(u)
        res_ag, pull_ag = u.atoms[res_idx], u.atoms[pull_idx]
        com_dis_df = pd.DataFrame(columns=['frame_idx', 'com_com_dist'])
        for frame_idx, _ in tqdm(enumerate(u.trajectory), total=u.trajectory.n_frames,
                                 desc=f'calculate com-com distance for chain {self.args.restrain_chain} and {self.args.pull_chain}', leave=False):
            com_com_dist = np.linalg.norm(res_ag.center_of_mass() - pull_ag.center_of_mass(), axis=0) / 10 # convert to nm
            com_dis_df.loc[frame_idx] = [frame_idx, com_com_dist]
        return com_dis_df
    
    def apply_center_atom(self, gmx: Gromacs, u: Universe, mdps: Dict[str, str]):
        _, res_idx = self.get_complex_atoms_index(u)
        res_ag = u.atoms[res_idx]
        center_pos = res_ag.center_of_mass()
        center_com_dist = np.linalg.norm(center_pos - res_ag.positions, axis=0)
        res_center_atom_idx = np.argmin(center_com_dist).astype(int)
        for mdp_name in mdps:
            mdp_config = opts_file(os.path.join(gmx.working_dir, mdps[mdp_name]))
            if 'pull-group1-pbcatom' in mdp_config:
                put_err(f"mdp file {mdp_name} already has pull-group1-pbcatom, skip.")
                continue
            new_mdp_config = mdp_config + f'\n\npull-group1-pbcatom = {res_center_atom_idx}'
            opts_file(os.path.join(gmx.working_dir, mdp_name), 'w', data=new_mdp_config)
    
    def run_sample(self, sample_gro: str, sample_main_name: str, gmx: Gromacs, mdps: Dict[str, str]):
        # STEP 6.1: run nvt equlibration
        gmx.run_gmx_with_expect('grompp', f=mdps['nvt'], c=sample_gro,
                                p='topol.top', r=sample_gro, n='pull.ndx',
                                o=f'{sample_main_name}_nvt.tpr', maxwarn=self.args.maxwarn)
        gmx.run_gmx_with_expect('mdrun -v', deffnm=f'{sample_main_name}_nvt')
        # STEP 6.2: run npt equlibration
        gmx.run_gmx_with_expect('grompp', f=mdps['npt'], c=f'{sample_main_name}_nvt.gro',
                                p='topol.top', r=f'{sample_main_name}_nvt.gro', t=f'{sample_main_name}_nvt.cpt', n='pull.ndx',
                                o=f'{sample_main_name}_npt.tpr', maxwarn=self.args.maxwarn)
        gmx.run_gmx_with_expect('mdrun -v', deffnm=f'{sample_main_name}_npt')
        # STEP 6.3: run md simulation
        gmx.run_gmx_with_expect('grompp', f=mdps['md'], c=f'{sample_main_name}_npt.gro',
                                t=f'{sample_main_name}_npt.cpt', p='topol.top', r=f'{sample_main_name}_npt.gro', n='pull.ndx',
                                o=f'{sample_main_name}_md.tpr', imd=f'{sample_main_name}_md.gro', maxwarn=self.args.maxwarn)
        gmx.run_gmx_with_expect('mdrun -v', deffnm=f'{sample_main_name}_md')

    @staticmethod
    def plot_hist(xvg_path: str, png_path: str, title: str = None):
        # 读取histograms.xvg（假设有多列，第一列是位置，后面每列是一个窗口的直方图）
        data = np.loadtxt(xvg_path, comments=['#', '@'])
        x = data[:, 0]
        histograms = data[:, 1:]
        np.save(png_path.replace('.png', '.npy'), histograms)
        plt.figure(figsize=(10, 6))
        for i in range(histograms.shape[1]):
            plt.plot(x, histograms[:, i], alpha=0.5, label=f'Window {i}')
        plt.grid(True)
        plt.xlabel('COM-COM Distance (nm)')
        plt.ylabel('Probability Density')
        plt.title('Umbrella Sampling Histograms')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(png_path, dpi=600)

    def main_process(self):
        # get protein paths
        if os.path.isdir(self.args.batch_dir):
            proteins_path = get_paths_with_extension(self.args.batch_dir, [], name_substr=self.args.protein_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.batch_dir}, exit.', _exit=True)
        put_log(f'get {len(proteins_path)} protein(s)')
        # check mdp files
        mdp_names = ['pull', 'nvt', 'npt', 'md']
        mdp_exist = list(map(lambda x: os.path.isfile(getattr(self.args, f'{x}_mdp')), mdp_names))
        if not all(mdp_exist):
            missing_names = [n for n, e in zip(mdp_names, mdp_exist) if not e]
            put_log(f'Warning: can not find mdp files in abspath: {", ".join(missing_names)}, skip.')
        # process each complex
        for protein_path in tqdm(proteins_path, total=len(proteins_path)):
            protein_path = Path(protein_path).resolve()
            main_name = protein_path.stem
            gmx = Gromacs(working_dir=str(protein_path.parent))
            # prepare gmx env and mdp files
            mdps = self.get_mdp(protein_path.parent, mdp_names)
            # check if md.tpr exists, if yes, skip STEP 1 ~ 3
            if os.path.exists(protein_path.parent / 'pull.tpr'):
                put_log(f'{protein_path} already done with pull.tpr, skip pulling run.')
            else:
                # STEP 1: make restraints index file
                self.make_restrain_idx(protein_path, gmx)
                # STEP 2: gmx grompp -f md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
                gmx.run_gmx_with_expect('grompp', f=mdps['pull'], c='npt.gro', p='topol.top', r='npt.gro', n='pull.ndx',
                                            t='npt.cpt', o='pull.tpr', maxwarn=self.args.maxwarn)
                # STEP 3: gmx mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg
                ## check apply pos restrain def
                if self.args.position_restrain:
                    self.apply_restrain_def(gmx)
                gmx.run_gmx_with_expect('mdrun -v', deffnm='pull', pf='pullf.xvg', px='pullx.xvg')
                gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f pullx.xvg -o pullx.png -t "Pull X of {main_name}" -csv pullx.csv -ns')
                gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f pullf.xvg -o pullf.png -t "Pull F of {main_name}" -csv pullf.csv -ns')
            # STEP 4: calculate com-com distance v.s. time curve
            u = Universe(str(protein_path.parent / 'pull.tpr'), str(protein_path.parent / 'pull.xtc'))
            com_dis_df = self.calcu_com_com_distance(u)
            com_dis_df.to_csv(protein_path.parent / 'pull_com_com_dist.csv', index=False)
            fig, ax = plt.subplots(figsize=(9, 6))
            ax.plot(com_dis_df['frame_idx'], com_dis_df['com_com_dist'])
            ax.set_xlabel('Frame Index', fontsize=14)
            ax.set_ylabel('COM-COM Distance (nm)', fontsize=14)
            ax.grid(linestyle='--')
            fig.savefig(protein_path.parent / 'pull_com_com_dist.png', dpi=600)
            ## calcu center atoms for each chain and add to mdp file
            if self.args.calcu_center_atom:
                self.apply_center_atom(gmx, u, mdps)
            ## check if pull-only
            if self.args.pull_only:
                put_log(f'pull only option enabled, skip sampling run.')
                continue
            # STEP 5: sample MD trajectory frame from pull trajectory
            start_dist = self.args.start_distance or com_dis_df['com_com_dist'].values[0]
            end_dist = self.args.end_distance or com_dis_df['com_com_dist'].max()
            sample_df = pd.DataFrame(columns=['dist', 'real_dist', 'frame_idx', 'sample_name', 'main_name'])
            for dist_i in tqdm(np.arange(start_dist, end_dist, self.args.distance_interval),
                                  desc=f'sample MD trajectory from pull trajectory {main_name}', leave=False):
                frame_idx = np.argmin(np.abs(com_dis_df['com_com_dist'] - dist_i))
                u.trajectory[frame_idx]
                real_dist = com_dis_df['com_com_dist'].values[frame_idx]
                sample_df.loc[len(sample_df)] = [dist_i, real_dist, frame_idx, f'pull_sample_{frame_idx}.gro', f'pull_sample_{frame_idx}']
                print(f'\n\nsample {dist_i:.2f} nm (real {real_dist:.2f} nm) at frame {frame_idx}')
                u.atoms.write(protein_path.parent / sample_df.loc[len(sample_df)-1, 'sample_name'], reindex=False)
            sample_df.to_csv(protein_path.parent / 'pull_sample.csv', index=False)
            # STEP 6: run MD simulation for each sample
            for sample_gro in tqdm(sample_df['sample_name'], desc=f'run MD simulation for each sample', leave=False):
                sample_main_name = sample_gro.split(".")[0]
                self.run_sample(sample_gro, sample_main_name, gmx, mdps)
            # STEP 7: perform WHAM analysis
            opts_file(str(protein_path.parent / 'pull_wham_tpr_lst.dat'), 'w',
                      data='\n'.join([f'{sample_main_name}_md.tpr' for sample_main_name in sample_df['main_name']]))
            opts_file(str(protein_path.parent / 'pull_wham_pullf.dat'), 'w',
                      data='\n'.join([f'{sample_main_name}_md_pullf.xvg' for sample_main_name in sample_df['main_name']]))
            gmx.run_gmx_with_expect('wham',  it='pull_wham_tpr_lst.dat', if_='pull_wham_pullf.dat',
                                    o='pull_wham_pme.xvg', hist='pull_wham_hist.xvg', unit='kCal')
            # STEP 8: plot WHAM histogram and curve
            gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f pull_wham_pme.xvg -o pull_wham_pme.png -t "WHAM PME of {main_name}" -csv pull_wham_pme.csv')
            self.plot_hist('pull_wham_hist.xvg', 'pull_wham_pme.png')


class any_sample(pull):
    HELP = """
    run any sample GROMACS simulation.
    This command is designed for running supplement pull sample while the first pull-sample is not enough.
    This command can be used for running any com-com distance sample.
    This command is a one-command-one-task command, which means it is not for batch processing.
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.iter_run_arg = [] # no batch processing, so just skip
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-n', '--protein-name', type = str, required=True,
                          help='protein name in each sub-directory, such as protein.gro.')
        args.add_argument('-s', '--suffix', type = str, required=True,
                          help='suffix for result name (exclude sample GROMACS name).')
        args.add_argument('-dis', '--distance', type=float, nargs='+', required=True,
                          help='start distance (nm) from pull trajectory for sampling MD run.')
        args.add_argument('--nvt-mdp', type = str, required=True,
                          help='nvt mdp file, if is a file-path, will copy to working directory.')
        args.add_argument('--npt-mdp', type = str, required=True,
                          help='npt mdp file, if is a file-path, will copy to working directory.')
        args.add_argument('--md-mdp', type = str, required=True,
                          help='production md mdp file, if is a file-path, will copy to working directory.')
        args.add_argument('--maxwarn', type=int, default=0,
                          help='maxwarn for em,nvt,npt,md gmx grompp command, default is %(default)s.')

    def process_args(self):
        # no batch processing, so just skip
        pass
    
    def main_process(self):
        # check mdp files
        mdp_names = ['nvt', 'npt', 'md']
        mdp_exist = list(map(lambda x: os.path.isfile(getattr(self.args, f'{x}_mdp')), mdp_names))
        if not all(mdp_exist):
            missing_names = [n for n, e in zip(mdp_names, mdp_exist) if not e]
            put_log(f'Warning: can not find mdp files in abspath: {", ".join(missing_names)}, skip.')
        # run sample
        suffix = self.args.suffix
        protein_path = Path(self.args.protein_name).resolve()
        main_name = protein_path.stem
        gmx = Gromacs(working_dir=str(protein_path.parent))
        # prepare gmx env and mdp files
        mdps = self.get_mdp(protein_path.parent, mdp_names)
        # check if pull.tpr exists, if not, exit
        if not os.path.exists(protein_path.parent / 'pull.tpr'):
            return put_err(f'can not find pull.tpr, exit.')
        # STEP 1: load com-com distance, 
        u = Universe(str(protein_path.parent / 'pull.tpr'), str(protein_path.parent / 'pull.xtc'))
        com_dis_df = pd.read_csv(protein_path.parent / 'pull_com_com_dist.csv')
        sample_df = pd.read_csv(protein_path.parent / 'pull_sample.csv')
        start_run_idx = len(sample_df)
        # STEP 2: extract sample gro files from pull trajectory
        for dist_i in tqdm(self.args.distance, desc=f'sample MD trajectory from pull trajectory {main_name}', leave=False):
            nearst_idx = np.argmin(np.abs(sample_df['dist'] - dist_i))
            frame_idx = np.argmin(np.abs(com_dis_df['com_com_dist'] - dist_i))
            if frame_idx in sample_df['frame_idx']:
                print(f'frame {frame_idx} has been sampled, skip.')
                continue
            u.trajectory[frame_idx]
            real_dist = com_dis_df['com_com_dist'].values[frame_idx]
            sample_df.loc[len(sample_df)] = [dist_i, real_dist, frame_idx, f'pull_sample_{frame_idx}.gro', f'pull_sample_{frame_idx}']
            print(f'\n\nsample {dist_i:.2f} nm (real {real_dist:.2f} nm) at frame {frame_idx}')
            print(f'nearst sample in previous run is {sample_df.loc[nearst_idx, "dist"]:.2f} nm')
            u.atoms.write(protein_path.parent / sample_df.loc[len(sample_df)-1, 'sample_name'], reindex=False)
        sample_df.to_csv(protein_path.parent / f'pull_sample_{suffix}.csv', index=False)
        # STEP 6: run MD simulation for each sample
        for sample_gro in tqdm(sample_df['sample_name'][start_run_idx:], desc=f'run MD simulation for each sample', leave=False):
            sample_main_name = sample_gro.split(".")[0]
            self.run_sample(sample_gro, sample_main_name, gmx, mdps)
        # STEP 7: perform WHAM analysis
        sample_df = sample_df.sort_values(by='dist')
        opts_file(str(protein_path.parent / f'pull_wham_tpr_lst_{suffix}.dat'), 'w',
                    data='\n'.join([f'{sample_main_name}_md.tpr' for sample_main_name in sample_df['main_name']]))
        opts_file(str(protein_path.parent / f'pull_wham_pullf_{suffix}.dat'), 'w',
                    data='\n'.join([f'{sample_main_name}_md_pullf.xvg' for sample_main_name in sample_df['main_name']]))
        gmx.run_gmx_with_expect('wham',  it=f'pull_wham_tpr_lst_{suffix}.dat', if_=f'pull_wham_pullf_{suffix}.dat',
                                o=f'pull_wham_pme_{suffix}.xvg', hist=f'pull_wham_hist_{suffix}.xvg', unit='kCal')
        # STEP 8: plot WHAM histogram and curve
        gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f pull_wham_pme_{suffix}.xvg -o pull_wham_pme_{suffix}.png -t "WHAM PME of {main_name}" -csv pull_wham_pme_{suffix}.csv')
        self.plot_hist(f'pull_wham_hist_{suffix}.xvg', f'pull_wham_pme_{suffix}.png')


_str2func = {
    'pull': pull,
    'any-sample': any_sample,
}


def main(sys_args: List[str] = None):
    make_args_and_excute('tools for GROMACS.', _str2func, sys_args)


if __name__ == "__main__":
    main()