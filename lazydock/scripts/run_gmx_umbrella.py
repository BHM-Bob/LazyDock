'''
Date: 2026-02-19
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2026-02-20 00:19:50
Description: steps most from http://www.mdtutorials.com/gmx/umbrella
'''
import argparse
import os
import re
import shutil
import time
from datetime import datetime
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

from lazydock.gmx.run import Gromacs
from lazydock.scripts._script_utils_ import Command, process_batch_dir_lst
from lazydock.scripts.run_gmx import simple_complex as _simple_complex
from lazydock.scripts.run_gmx import simple_protein as _simple_protein
from lazydock.utils import uuid4


class pull(_simple_protein):
    HELP = """
    run pull process for umbrella sampling, steps most from http://www.mdtutorials.com/gmx/umbrella
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
        args.add_argument('-posres', '--position-restrain', action='store_true', default=False,
                          help='add position restrain define in topol-itp file for restrain chain, default is %(default)s.')
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
        
    def calcu_com_com_distance(self, u: Universe):
        res_idx, pull_idx = self.get_complex_atoms_index(u)
        res_ag, pull_ag = u.atoms[res_idx], u.atoms[pull_idx]
        com_dis_df = pd.DataFrame(columns=['frame_idx', 'com_com_dist'])
        for frame_idx, _ in tqdm(enumerate(u.trajectory), total=u.trajectory.n_frames,
                                 desc=f'calculate com-com distance for chain {self.args.restrain_chain} and {self.args.pull_chain}', leave=False):
            com_com_dist = np.linalg.norm(res_ag.center_of_mass() - pull_ag.center_of_mass(), axis=0) / 10 # convert to nm
            com_dis_df.loc[frame_idx] = [frame_idx, com_com_dist]
        return com_dis_df

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
            # STEP 5: sample MD trajectory frame from pull trajectory
            start_dist = self.args.start_distance or com_dis_df['com_com_dist'].values[0]
            end_dist = self.args.end_distance or com_dis_df['com_com_dist'].max()
            sample_df = pd.DataFrame(columns=['dist_idx', 'real_dist', 'frame_idx', 'sample_name'])
            for dist_i in tqdm(np.arange(start_dist, end_dist, self.args.distance_interval),
                                  desc=f'sample MD trajectory from pull trajectory {main_name}', leave=False):
                frame_idx = np.argmin(np.abs(com_dis_df['com_com_dist'] - dist_i))
                u.trajectory[frame_idx]
                real_dist = com_dis_df['com_com_dist'].values[frame_idx]
                sample_df.loc[len(sample_df)] = [dist_i, real_dist, frame_idx, f'pull_sample_{frame_idx}.gro']
                print(f'\n\nsample {dist_i:.2f} nm (real {real_dist:.2f} nm) at frame {frame_idx}')
                u.atoms.write(protein_path.parent / sample_df.loc[len(sample_df)-1, 'sample_name'], reindex=False)
            sample_df.to_csv(protein_path.parent / 'pull_sample.csv', index=False)
            # STEP 6: run MD simulation for each sample
            for sample_gro in tqdm(sample_df['sample_name'], desc=f'run MD simulation for each sample', leave=False):
                sample_main_name = sample_gro.split(".")[0]
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
                                        t=f'{sample_main_name}_npt.cpt', p='topol.top', n='pull.ndx',
                                        o=f'{sample_main_name}_md.tpr', imd=f'{sample_main_name}_md.gro', maxwarn=self.args.maxwarn)
                gmx.run_gmx_with_expect('mdrun -v', deffnm=sample_main_name)
            

class simple_protein(_simple_protein):
    HELP = """
    run simple protein GROMACS simulation
    1. gmx editconf -f protein.gro -o protein_newbox.gro -c -d 1.0 -bt cubic
    2. gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top
    3. gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr
    4. gmx genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
    5. gmx grompp -f minim.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr
    6. gmx mdrun -v -deffnm em
    7. gmx energy -f em.edr -o potential.xvg # At the prompt, type "10 0" to select Potential (10); zero (0) terminates input.
    8. gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    9. gmx mdrun -deffnm nvt
    10. gmx energy -f nvt.edr -o temperature.xvg # Type "16 0" at the prompt to select the temperature of the system and exit.
    11. gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    12. gmx mdrun -deffnm npt
    13. gmx energy -f npt.edr -o pressure.xvg # Type "18 0" at the prompt to select the pressure of the system and exit.
    14. gmx energy -f npt.edr -o density.xvg # using energy and entering "24 0" at the prompt.
    15. gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
    16. gmx mdrun -v -ntomp 4 -deffnm md -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu
    
    if step 16 terminated, you can use gmx mdrun -s md.tpr -cpi md.cpt -v -ntomp 4 -deffnm md -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.indexs = {}
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        _simple_protein.make_args(args)
        return args
        
    def energy_minimization(self, protein_path: Path, main_name: str, gmx: Gromacs, mdps: Dict[str, str]):
        # STEP 5: grompp -f minim.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr
        gmx.run_gmx_with_expect('grompp', f=mdps['em'], c=f'{main_name}_solv_ions.gro',
                                    p='topol.top', o='em.tpr', maxwarn=self.args.maxwarn)
        # STEP 6: mdrun -v -deffnm em
        gmx.run_gmx_with_expect(f'mdrun {self.args.em_args}', deffnm='em')
        # STEP 7: energy -f em.edr -o potential.xvg
        gmx.run_gmx_with_expect('energy', f='em.edr', o='potential.xvg',
                                    expect_actions=[{'line or a zero.': f'{self.args.potential_groups}\r'}])
        os.system(f'cd "{protein_path.parent}" && dit xvg_compare -c 1 -f potential.xvg -o potential.png -t "EM Potential of {main_name}" -csv {main_name}_potential.csv -ns')
        
    def equilibration(self, protein_path: Path, main_name: str, gmx: Gromacs, mdps: Dict[str, str]):
        # STEP 8: grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
        gmx.run_gmx_with_expect('grompp', f=mdps['nvt'], c='em.gro', r='em.gro', p='topol.top',
                                    o='nvt.tpr', n=self.indexs.get('nvt', None), maxwarn=self.args.maxwarn)
        # STEP 9: mdrun -deffnm nvt
        gmx.run_gmx_with_expect('mdrun', deffnm='nvt')
        # STEP 10: energy -f nvt.edr -o temperature.xvg
        gmx.run_gmx_with_expect('energy', f='nvt.edr', o='temperature.xvg',
                                    expect_actions=[{'line or a zero.': f'{self.args.temperature_groups}\r'}])
        os.system(f'cd "{protein_path.parent}" && dit xvg_compare -c 1 -f temperature.xvg -o temperature.png -smv -ws 10 -t "NVT Temperature of {main_name}" -csv {main_name}_temperature.csv -ns')
        # STEP 11: grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
        gmx.run_gmx_with_expect('grompp', f=mdps['npt'], c='nvt.gro', r='nvt.gro', t='nvt.cpt',
                                    p='topol.top', o='npt.tpr', n=self.indexs.get('npt', None), maxwarn=self.args.maxwarn)
        # STEP 12: mdrun -deffnm npt
        gmx.run_gmx_with_expect('mdrun', deffnm='npt')
        # STEP 13: energy -f npt.edr -o pressure.xvg
        gmx.run_gmx_with_expect('energy', f='npt.edr', o='pressure.xvg',
                                    expect_actions=[{'line or a zero.': f'{self.args.pressure_groups}\r'}])
        os.system(f'cd "{protein_path.parent}" && dit xvg_compare -c 1 -f pressure.xvg -o pressure.png -smv -ws 10 -t "NPT Pressure of {main_name}" -csv {main_name}_pressure.csv -ns')
        # STEP 14: energy -f npt.edr -o density.xvg
        gmx.run_gmx_with_expect('energy', f='npt.edr', o='density.xvg',
                                    expect_actions=[{'line or a zero.': f'{self.args.density_groups}\r'}])
        os.system(f'cd "{protein_path.parent}" && dit xvg_compare -c 1 -f density.xvg -o density.png -smv -ws 10 -t "NPT Density of {main_name}" -csv {main_name}_density.csv -ns')
        
    def production_md(self, protein_path: Path, main_name: str, gmx: Gromacs, mdps: Dict[str, str]):
        # STEP 15: grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
        gmx.run_gmx_with_expect('grompp', f=mdps['md'], c='npt.gro', t='npt.cpt', p='topol.top',
                                    o='md.tpr', imd='md.gro', n=self.indexs.get('md', None), maxwarn=self.args.maxwarn)
        # STEP 16: mdrun -v -ntomp 4 -deffnm md -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu
        gmx.run_gmx_with_expect(f'mdrun {self.args.mdrun_args}', deffnm='md')

    def main_process(self):
        # get protein paths
        if os.path.isdir(self.args.batch_dir):
            proteins_path = get_paths_with_extension(self.args.batch_dir, [], name_substr=self.args.protein_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.batch_dir}, exit.', _exit=True)
        put_log(f'get {len(proteins_path)} protein(s)')
        # check mdp files
        mdp_names = ['ion', 'em', 'nvt', 'npt','md']
        mdp_exist = list(map(lambda x: os.path.isfile(getattr(self.args, f'{x}_mdp')), mdp_names))
        if not all(mdp_exist):
            missing_names = [n for n, e in zip(mdp_names, mdp_exist) if not e]
            put_log(f'Warning: can not find mdp files in abspath: {", ".join(missing_names)}, skip.')
        # sleep until start time
        self.sleep_until_start_time()
        # process each complex
        for protein_path in tqdm(proteins_path, total=len(proteins_path)):
            protein_path = Path(protein_path).resolve()
            main_name = protein_path.stem
            # check if md.tpr exists, if yes, skip
            if os.path.exists(protein_path.parent / 'md.tpr'):
                put_log(f'{protein_path} already done with md.tpr, skip.')
                continue
            # prepare gmx env and mdp files
            gmx = Gromacs(working_dir=str(protein_path.parent))
            mdps = self.get_mdp(protein_path.parent)
            # STEP 1 ~ 4: make box, solvate, ions
            self.make_box(protein_path, main_name, gmx, mdps)
            # STEP 5 ~ 7: energy minimization
            self.energy_minimization(protein_path, main_name, gmx, mdps)
            # STEP 8 ~ 14: equilibration
            self.equilibration(protein_path, main_name, gmx, mdps)
            exit()
            # STEP 15 ~ 16: production md
            self.production_md(protein_path, main_name, gmx, mdps)


class simple_complex(_simple_complex):
    HELP = _simple_complex.HELP.replace('protein', 'complex')
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = _simple_complex.make_args(args)
        return args
        
    def equilibration(self, protein_path: Path, main_name: str, gmx: Gromacs, mdps: Dict[str, str]):
        from lazydock.scripts.prepare_gmx import ligand
        lig_name = Path(self.args.ligand_name).stem
        # STEP 1: make restraints file for ligand
        # make_ndx -f jz4.gro -o index_jz4.ndx
        gmx.run_gmx_with_expect('make_ndx', f=f'{lig_name}.gro', o=f'index_{lig_name}.ndx',
                                    expect_actions=[{'>': '0 & ! a H*\r'}, {'>': 'q\r'}])
        # gmx genrestr -f jz4.gro -n index_jz4.ndx -o posre_jz4.itp -fc 1000 1000 1000
        gmx.run_gmx_with_expect('genrestr', f=f'{lig_name}.gro', n=f'index_{lig_name}.ndx', o=f'posre_{lig_name}.itp', fc='1000 1000 1000',
                                     expect_actions=[{'Select a group:': '3\r'}])
        # STEP 2: add restraints info into topol.top
        res_info = f'\n; Ligand position restraints\n#ifdef {self.args.lig_posres}\n#include "posre_{lig_name}.itp"\n#endif\n\n'
        ligand.insert_content(protein_path.parent / 'topol.top', f'#include "{lig_name}.itp"\n', res_info)
        # STEP 3: make tc-grps index file
        # gmx make_ndx -f em.gro -o index.ndx
        tc_groups = self.args.tc_groups
        if self.args.tc_groups == 'auto':
            groups = gmx.get_groups('em.tpr')
            if 'Protein' not in groups and 'LIG' not in groups:
                put_err(f'can not find Protein or LIG group in em.tpr, skip.')
            else:
                tc_groups = f'{groups["Protein"]} | {groups["LIG"]}'
        gmx.run_gmx_with_expect('make_ndx', f='em.gro', o='tc_index.ndx',
                                    expect_actions=[{'>': f'{tc_groups}\r'}, {'>': 'q\r'}])
        for k in ['nvt', 'npt','md']:
            self.indexs[k] = 'tc_index.ndx'
        super().equilibration(protein_path, main_name, gmx, mdps)


_str2func = {
    'pull': pull,
    'simple-protein': simple_protein,
    'simple-complex': simple_complex,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for n, func in _str2func.items():
        func.make_args(subparsers.add_parser(n, description=func.HELP))

    args = args_paser.parse_args(sys_args)
    if args.sub_command in _str2func:
        _str2func[args.sub_command](args).excute()


if __name__ == "__main__":
    # pass
    # main(r'complex -d data_tmp/gmx/complex -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir data_tmp/gmx/charmm36-jul2022.ff'.split())
    
    main()