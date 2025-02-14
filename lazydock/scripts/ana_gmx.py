import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

if 'MBAPY_PLT_AGG' in os.environ:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
from tqdm import tqdm
from MDAnalysis import Universe
from mbapy_lite.base import put_err, put_log
from mbapy.plot import save_show
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.web import TaskPool
from lazydock.gmx.run import Gromacs
from lazydock.gmx.mda.convert import PDBConverter
from lazydock.pml.interaction_utils import calcu_pdbstr_interaction
from lazydock.pml.plip_interaction import run_plip_analysis
from lazydock.pml.rrcs import calcu_RRCS_from_tensor, calcu_RRCS_from_array
from lazydock.scripts.ana_interaction import simple_analysis, pml_mode, plip_mode
from lazydock.scripts._script_utils_ import Command, clean_path, excute_command


class simple(Command):
    HELP = """
    simple analysis for GROMACS simulation
    0. gmx_mpi trjconv -s md.tpr -f md.xtc -o md_center.xtc -pbc mol -center
    
    1. gmx_mpi rms -s md.tpr -f md_center.xtc -o rmsd.xvg -tu ns 
    2. gmx_mpi rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg
    3. gmx_mpi gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg
    4. gmx_mpi hbond -s md.tpr -f md_center.xtc -num -dt 10
    
    5. gmx_mpi sasa -s md.tpr -f md_center.xtc -o sasa_total.xvg -or sasa_res.xvg -tu ns 
    6. gmx_mpi covar -s md.tpr -f md_center.xtc -o eigenval.xvg -tu ns 
    
    7. free energy landscape from rmsd and gyrate by MD-DaVis
    8. Probability Density Function from rmsd and gyrate
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str,
                          help='directory to store the prepared files')
        args.add_argument('-n', '--main-name', type = str,
                          help='main name in each sub-directory, such as md.tpr.')
        args.add_argument('-cg', '--center-group', type = str, default='1',
                          help='group to center the trajectory, default is %(default)s.')
        args.add_argument('-rg', '--rms-group', type = str, default='4',
                          help='group to calculate rmsd, rmsf, and gyrate, default is %(default)s.')
        args.add_argument('-hg', '--hbond-group', type = str, default='1',
                          help='group to calculate hbond, default is %(default)s.')
        args.add_argument('-sg', '--sasa-group', type = str, default='4',
                          help='group to calculate sasa, default is %(default)s.')
        args.add_argument('-eg', '--eigenval-group', type = str, default='4',
                          help='group to calculate eigenval, default is %(default)s.')
        args.add_argument('-xmax', '--eigenval-xmax', type = int, default=15,
                          help='max value of eigenval, default is %(default)s.')
        return args

    @staticmethod
    def trjconv(gmx: Gromacs, main_name: str, center_group: str = '1', **kwargs):
        gmx.run_gmx_with_expect('trjconv', s=f'{main_name}.tpr', f=f'{main_name}.xtc', o=f'{main_name}_center.xtc', pbc='mol', center=True,
                                    expect_actions=[{'Select a group:': f'{center_group}\r', '\\timeout': f'{center_group}\r'},
                                                    {'Select a group:': '0\r', '\\timeout': '0\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        
    @staticmethod
    def rms(gmx: Gromacs, main_name: str, group: str = '4', **kwargs):
        gmx.run_gmx_with_expect('rms', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'rmsd.xvg', tu='ns',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f rmsd.xvg -o rmsd.png -smv -ws 10 -t "RMSD of {main_name}" -csv {main_name}_rmsd.csv -ns')
        
    @staticmethod
    def rmsf(gmx: Gromacs, main_name: str, group: str = '4', res: bool = True, **kwargs):
        gmx.run_gmx_with_expect('rmsf', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'rmsf.xvg', res=res,
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f rmsf.xvg -o rmsf.png -t "RMSF of {main_name}" -csv {main_name}_rmsf.csv -ns')
        
    @staticmethod
    def gyrate(gmx: Gromacs, main_name: str, group: str = '4', **kwargs):
        gmx.run_gmx_with_expect('gyrate', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'gyrate.xvg',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f gyrate.xvg -o gyrate.png -smv -ws 10 -t "Gyrate of {main_name}" -csv {main_name}_gyrate.csv -ns')
        
    @staticmethod
    def hbond(gmx: Gromacs, main_name: str, group: str = '1', dt=10, **kwargs):
        gmx.run_gmx_with_expect('hbond', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc',
                                    num=f'{main_name}_hbond_num.xvg', dist=f'{main_name}_hbond_dist.xvg',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f {main_name}_hbond_num.xvg -o hbond_num.png -smv -ws 10 -t "H-bond num of {main_name}" -csv {main_name}_hbond_num.csv -ns')
        gmx.run_command_with_expect(f'dit xvg_show -f {main_name}_hbond_dist.xvg -o hbond_dist.png -ns')

    @staticmethod
    def sasa(gmx: Gromacs, main_name: str, group: str = '4', **kwargs):
        gmx.run_gmx_with_expect('sasa -or sasa_res.xvg', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc',
                                    o=f'sasa_total.xvg', odg=f'sasa_dg.xvg', tv='sasa_tv.xvg', tu='ns',
                                    expect_actions=[{'>': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        for ty in ['total', 'res', 'dg', 'tv']:
            gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f sasa_{ty}.xvg -o sasa_{ty}.png -smv -ws 10 -t "SASA {ty} of {main_name}" -csv {main_name}_sasa_{ty}.csv -ns')

    @staticmethod
    def covar(gmx: Gromacs, main_name: str, group: str = '4', xmax: int = 15, **kwargs):
        gmx.run_gmx_with_expect('covar', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'eigenval.xvg', tu='ns',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_command_with_expect(f'dit xvg_compare -c 1 -f eigenval.xvg -o eigenval.png -xmin 0 -xmax {xmax} -smv -ws 10 -t "Eigenval of {main_name}" -csv {main_name}_eigenval.csv -ns')
    
    @staticmethod
    def free_energy_landscape(gmx: Gromacs, main_name: str, **kwargs):
        # MD-DaVis
        gmx.run_command_with_expect(f'md-davis landscape_xvg -c -T 300 -x rmsd.xvg -y gyrate.xvg -o FEL.html -n FEL -l "RMSD-Rg" --axis_labels "dict(x=\'RMSD (in nm)\', y=\'Rg (in nm)\', z=\'Free Energy (kJ mol<sup>-1</sup>)<br>\')"')
        # gmx and dit
        gmx.run_command_with_expect(f'dit xvg_combine -f rmsd.xvg gyrate.xvg -c 0,1 1 -l RMSD Gyrate -o rmsd_gyrate.xvg -x "Time (ps)"')
        gmx.run_command_with_expect(f'sham -f rmsd_gyrate.xvg -ls sham.xpm')
        gmx.run_command_with_expect(f'it xpm_show -f sham.xpm -m 3d --x_precision 1 --y_precision 1 --z_precision 1 -cmap jet --colorbar_location right -o rmsd_gyrate.png -ns')
        
    
    @staticmethod
    def plot_PDF(gmx: Gromacs, main_name: str, **kwargs):
        """idea from https://pymolwiki.org/index.php/Geo_Measures_Plugin"""
        # read data and calculate density
        x = pd.read_csv(f'{gmx.working_dir}/{main_name}_rmsd.csv').values[:, -1]
        y = pd.read_csv(f'{gmx.working_dir}/{main_name}_gyrate.csv').values[:, -1]
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
        # plot scatter
        fig, ax = plt.subplots()
        pdf = ax.scatter(x, y, c=z, s=50, edgecolor="none", cmap=plt.cm.jet)
        # Hide right and top spines
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")
        # Set x and y limits
        plt.xlim(x.min() - 1, x.max() + 1)
        plt.ylim(y.min() - 1, y.max() + 1)
        # Set x and y labels
        plt.xlabel('RMSD (in nm)')
        plt.ylabel('Rg (in nm)')
        # Adding the color bar
        colbar = plt.colorbar(pdf)
        colbar.set_label("Probability Density Function")
        save_show(os.path.join(gmx.working_dir, f'{main_name}_PDF.png'), 600, show=False)
        plt.close(fig)
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        if not os.path.isdir(self.args.dir):
            put_err(f'dir argument should be a directory: {self.args.dir}, exit.', _exit=True)
        
    def main_process(self):
        # get complex paths
        complexs_path = get_paths_with_extension(self.args.dir, [], name_substr=self.args.main_name)
        put_log(f'get {len(complexs_path)} task(s)')
        # process each complex
        for complex_path in tqdm(complexs_path, total=len(complexs_path)):
            complex_path = Path(complex_path).resolve()
            gmx = Gromacs(working_dir=str(complex_path.parent))
            # perform trjconv
            self.trjconv(gmx, main_name=complex_path.stem, center_group=self.args.center_group)
            # perform analysis
            self.rms(gmx, main_name=complex_path.stem, group=self.args.rms_group)
            self.rmsf(gmx, main_name=complex_path.stem, group=self.args.rms_group)
            self.gyrate(gmx, main_name=complex_path.stem, group=self.args.rms_group)
            self.hbond(gmx, main_name=complex_path.stem, group=self.args.hbond_group)
            self.sasa(gmx, main_name=complex_path.stem, group=self.args.sasa_group)
            self.covar(gmx, main_name=complex_path.stem, group=self.args.eigenval_group, xmax=self.args.eigenval_xmax)
            # perform free energy landscape by MD-DaVis
            self.free_energy_landscape(gmx, main_name=complex_path.stem)
            # plot PDF
            self.plot_PDF(gmx, main_name=complex_path.stem)
            
            
class mmpbsa(simple_analysis):
    HELP = """
    mmpbsa analysis for GROMACS simulation
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser, mmpbsa_args: bool = True):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, required=True,
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        if mmpbsa_args:
            args.add_argument('-i', '--input', type = str, required=True,
                              help=f"gmx_MMPBSA input file name in each sub-folder, such as mmpbsa.in")
            args.add_argument('-np', '--np', type = int, required=True,
                              help=f"npi np argument for gmx_MMPBSA")
        args.add_argument('-top', '--top-name', type = str, required=True,
                          help=f"topology file name in each sub-folder.")
        args.add_argument('-traj', '--traj-name', type = str, required=True,
                          help=f"trajectory file name in each sub-folder.")
        args.add_argument('--receptor-chain-name', type = str, required=True,
                          help='receptor chain name, such as "A".')
        args.add_argument('--ligand-chain-name', type = str, required=True,
                          help='ligand chain name, such as "LIG".')
        return args
        
    def get_complex_atoms_index(self, u: Universe):
        rec_idx = u.atoms.chainIDs == self.args.receptor_chain_name
        lig_idx = u.atoms.chainIDs == self.args.ligand_chain_name
        put_log(f"receptor atoms: {rec_idx.sum()}, ligand atoms: {lig_idx.sum()}.")
        return rec_idx, lig_idx
    
    def check_top_traj(self):
        top_paths = get_paths_with_extension(self.args.batch_dir, [os.path.split(self.args.top_name)[-1]], name_substr=self.args.top_name)
        traj_paths = get_paths_with_extension(self.args.batch_dir, [os.path.split(self.args.traj_name)[-1]], name_substr=self.args.traj_name)
        if not self.check_file_num_paried(top_paths, traj_paths):
            put_err(f"The number of top and traj files is not equal, please check the input files.\ninvalid roots:{self.invalid_roots}", _exit=True)
        return top_paths, traj_paths
    
    def process_args(self):
        self.args.batch_dir = clean_path(self.args.batch_dir)
        if not os.path.isdir(self.args.batch_dir):
            put_err(f'batch_dir argument should be a directory: {self.args.batch_dir}, exit.', _exit=True)
        # load origin dfs from data file
        self.top_paths, self.traj_paths = self.check_top_traj()
        for r_path, l_path in zip(self.top_paths, self.traj_paths):
            self.tasks.append((r_path, l_path))
    
    def main_process(self):
        print(f'find {len(self.tasks)} tasks.')
        # run tasks
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            # get receptor and ligand atoms index range
            u = Universe(top_path, traj_path)
            rec_idx, lig_idx = self.get_complex_atoms_index(u)
            rec_range_str, lig_range_str = f"{rec_idx.min()}-{rec_idx.max()}", f"{lig_idx.min()}-{lig_idx.max()}"
            # make index file for receptor and ligand
            gmx = Gromacs(working_dir=wdir)
            gmx.run_gmx_with_expect('make_ndx', f=os.path.basename(top_path), o='mmpbsa_tmp.ndx',
                                        expect_actions=[{'>': 'q\r'}])
            sum_groups = opts_file(os.path.join(gmx.working_dir, 'mmpbsa_tmp.ndx')).count(']')
            gmx.run_gmx_with_expect('make_ndx', f=os.path.basename(top_path), o='mmpbsa.ndx',
                                        expect_actions=[{'>': f'a {rec_range_str}\r'}, {'>': f'name {sum_groups+1} MMPBSA_Receptor\r'},
                                                        {'>': f'a {lig_range_str}\r'}, {'>': f'name {sum_groups+2} MMPBSA_Ligand\r'},
                                                        {'>': 'q\r'}])
            # call gmx_MMPBSA
            cmd_str = f'gmx_MMPBSA -O -i {self.args.input} -cs {self.args.top_name} -ct {self.args.traj_name} -ci mmpbsa.in -cg {sum_groups+1} {sum_groups+2} -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv'
            os.system(f'cd "{gmx.working_dir}" && mpirun -np {self.args.np} {cmd_str}')
            bar.update(1)
    
    
def run_pdbstr_interaction_analysis(pdbstr: str, receptor_chain: str, ligand_chain: str,
                                    method: str, mode: str, cutoff: float, hydrogen_atom_only: bool, **kwargs):
    if method == 'pymol':
        inter = calcu_pdbstr_interaction(f'chain {receptor_chain}', f'chain {ligand_chain}', pdbstr, mode, cutoff, hydrogen_atom_only)
    elif method == 'plip':
       inter = run_plip_analysis(pdbstr, receptor_chain, ligand_chain, mode, cutoff)
    else:
        return put_err(f"method {method} not supported, return None.")
    return inter


class interaction(mmpbsa):
    HELP = """
    interaction analysis for GROMACS simulation
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.alter_chain = None
        self.alter_res = None
        self.alter_atm = None

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        mmpbsa.make_args(args, mmpbsa_args=False)
        args.add_argument('--alter-ligand-chain', type = str, default=None,
                          help='alter ligand chain name from topology to user-define, such as "Z".')
        args.add_argument('--alter-ligand-res', type = str, default=None,
                          help='alter ligand res name from topology to user-define, such as "UNK".')
        args.add_argument('--alter-ligand-atm', type = str, default=None,
                          help='alter ligand atom type from topology to user-define, such as "HETATM".')
        args.add_argument('--method', type = str, default='pymol', choices=['pymol', 'plip'],
                          help='interaction method, default is %(default)s.')
        args.add_argument('--mode', type = str, default='all',
                          help=f'interaction mode, multple modes can be separated by comma, all method support `\'all\'` model.\npymol: {",".join(pml_mode)}\nplip: {",".join(plip_mode)}')
        args.add_argument('--cutoff', type = float, default=4,
                          help='distance cutoff for interaction calculation, default is %(default)s.')
        args.add_argument('--hydrogen-atom-only', default=False, action='store_true',
                          help='only consider hydrogen bond acceptor and donor atoms, this only works when method is pymol, default is %(default)s.')
        args.add_argument('--output-style', type = str, default='receptor', choices=['receptor'],
                          help='output style\n receptor: resn resi distance')
        args.add_argument('--ref-res', type = str, default='',
                          help='reference residue name, input string shuld be like GLY300,ASP330, also support a text file contains this format string as a line.')
        args.add_argument('-np', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        args.add_argument('-b', '--begin-frame', type=int, default=0,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-e', '--end-frame', type=int, default=None,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-step', '--traj-step', type=int, default=1,
                          help='Step while reading trajectory. Default is %(default)s.')
    
    def process_args(self):
        # self.args.alter_ligand_chain will passed to final interaction calcu function
        if self.args.alter_ligand_chain is not None:
            self.alter_chain = {self.args.ligand_chain_name: self.args.alter_ligand_chain}
        else:
            self.args.alter_ligand_chain = self.args.ligand_chain_name
        # in lazydock, default ligand res name is the chain name too, so alter chain name to alter-ligand-res
        if self.args.alter_ligand_res is not None:
            self.alter_res = {self.args.ligand_chain_name: self.args.alter_ligand_res}
        # set self.alter_atm
        if self.args.alter_ligand_atm is not None:
            self.alter_atm = {self.args.ligand_chain_name: self.args.alter_ligand_atm}
        # output formater
        self.output_formater = getattr(self, f'output_fromater_{self.args.output_style}')
        # check batch dir AND top and traj
        mmpbsa.process_args(self)
        # check method and mode AND load ref_res
        self.args.batch_dir = str(self.args.batch_dir)
        simple_analysis.process_args(self)
        
    def main_process(self):
        print(f'find {len(self.tasks)} tasks.')
        # run tasks
        pool = TaskPool('process', self.args.n_workers).start()
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            # load pdbstr from traj
            u = Universe(top_path, traj_path)
            rec_idx, lig_idx = self.get_complex_atoms_index(u)
            complex_ag = u.atoms[rec_idx | lig_idx]
            # calcu interaction for each frame
            sum_frames = (len(u.trajectory) if self.args.end_frame is None else self.args.end_frame) - self.args.begin_frame
            for frame in tqdm(u.trajectory[self.args.begin_frame:self.args.end_frame:self.args.traj_step],
                              total=sum_frames//self.args.traj_step, desc='Calculation frames', leave=False):
                pdbstr = PDBConverter(complex_ag).fast_convert(alter_chain=self.alter_chain, alter_res=self.alter_res, alter_atm=self.alter_atm)
                pool.add_task(frame.time, run_pdbstr_interaction_analysis, pdbstr,
                              self.args.receptor_chain_name, self.args.alter_ligand_chain,
                              self.args.method, self.args.mode, self.args.cutoff, self.args.hydrogen_atom_only)
                pool.wait_till(lambda: pool.count_waiting_tasks() == 0, 0.001, update_result_queue=False)
            # merge interactions
            interactions, df = {}, pd.DataFrame()
            for k in list(pool.tasks.keys()):
                i = len(df)
                interactions[k] = pool.query_task(k, True, 10)
                df.loc[i, 'time'] = k
                df.loc[i, 'ref_res'] = ''
                for inter_mode, inter_value in interactions[k].items():
                    fmt_string = self.output_formater(inter_value, self.args.method)
                    df.loc[i, inter_mode] = fmt_string
                    for r in self.args.ref_res:
                        if r in fmt_string and not r in df.loc[i,'ref_res']:
                            df.loc[i,'ref_res'] += f'{r},'
            # save result
            top_path = Path(top_path).resolve()
            df.to_csv(str(top_path.parent / f'{top_path.stem}_{self.args.method}_interactions.csv'), index=False)
            opts_file(str(top_path.parent / f'{top_path.stem}_{self.args.method}_interactions.pkl'), 'wb', way='pkl', data=interactions)
            # other things
            pool.task = {}
            bar.update(1)
        pool.close()


def _calcu_RRCS(resis: np.ndarray, names: np.ndarray, positions: np.ndarray,
                occupancies: np.ndarray, backend: str):
    if backend == 'numpy':
        return calcu_RRCS_from_array(names, resis, positions, occupancies)
    elif backend in {'torch', 'cuda'}:
        device = 'cuda' if backend == 'cuda' else 'cpu'
        import torch
        resis = torch.tensor(resis, dtype=torch.int32, device=device)
        sort_idx = torch.argsort(resis)
        resis = resis[sort_idx]
        names = np.array(names)[sort_idx.cpu().numpy()]
        positions = torch.tensor(positions, dtype=torch.float32, device=device)[sort_idx]
        occupancies = torch.tensor(occupancies, dtype=torch.float32, device=device)[sort_idx]
        return calcu_RRCS_from_tensor(names, resis, positions, occupancies, device=device)


class RRCS(mmpbsa):
    HELP = """
    RRCS analysis for GROMACS simulation
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, required=True,
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('-top', '--top-name', type = str, required=True,
                          help=f"topology file name in each sub-folder.")
        args.add_argument('-traj', '--traj-name', type = str, required=True,
                          help=f"trajectory file name in each sub-folder.")
        args.add_argument('-c', '--chains', type = str, default=None, nargs='+',
                          help='chain of molecular to be included into calculation. Default is %(default)s.')
        args.add_argument('-np', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        args.add_argument('-b', '--begin-frame', type=int, default=0,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-e', '--end-frame', type=int, default=None,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-step', '--traj-step', type=int, default=1,
                          help='Step while reading trajectory. Default is %(default)s.')
        args.add_argument('--backend', type=str, default='numpy', choices=['numpy', 'torch', 'cuda'],
                          help='backend for RRCS calculation. Default is %(default)s.')
    
    def main_process(self):
        print(f'find {len(self.tasks)} tasks.')
        # run tasks
        pool = TaskPool('process', self.args.n_workers).start()
        bar = tqdm(total=len(self.tasks), desc='Calculating RRCS')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            bar.set_description(f"{wdir}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            # load pdbstr from traj
            u = Universe(top_path, traj_path)
            if self.args.chains is not None and len(self.args.chains):
                idx = u.atoms.chainIDs == self.args.chains[0]
                for chain_i in self.args.chains[1:]:
                    idx = idx | (u.atoms.chainIDs == chain_i)
                ag = u.atoms[idx]
            else:
                ag = u.atoms
            print(f'find {np.unique(u.atoms.chainIDs)} chains, {len(ag)} atoms in {self.args.chains}')
            # calcu interaction for each frame
            sum_frames = (len(u.trajectory) if self.args.end_frame is None else self.args.end_frame) - self.args.begin_frame
            for frame in tqdm(u.trajectory[self.args.begin_frame:self.args.end_frame:self.args.traj_step],
                              total=sum_frames//self.args.traj_step, desc='Calculation frames', leave=False):
                if not hasattr(ag, 'occupancies'):
                    occupancies = np.ones(len(ag))
                else:
                    occupancies = ag.occupancies
                pool.add_task(frame.time, _calcu_RRCS, ag.resids.copy(), ag.names.copy(),
                              ag.positions.copy(), occupancies.copy(), self.args.backend)
                pool.wait_till(lambda: pool.count_waiting_tasks() == 0, 0.01, update_result_queue=False)
            # merge result
            scores, frames = {}, list(pool.tasks.keys())
            for k in frames:
                df = pool.query_task(k, True, 10)
                scores[k] = df.values
            # save result
            top_path = Path(top_path).resolve()
            np.savez(str(top_path.parent / f'{top_path.stem}_RRCS.npz'),
                     scores=scores, frames=np.array(frames), resis=ag.resids, resns=ag.resnames, chains=ag.chainIDs)
            # other things
            pool.task = {}
            bar.update(1)
        pool.close()
        

_str2func = {
    'simple': simple,
    'mmpbsa': mmpbsa,
    'interaction': interaction,
    'rrcs': RRCS,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS analysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    # dev code
    # main('interaction -bd data_tmp/gmx/run1 -top md.tpr -traj md_center.xtc -np 8 --receptor-chain-name A --ligand-chain-name LIG --alter-ligand-chain Z --alter-ligand-res UNK --method plip --mode all'.split(' '))
    # main('rrcs -d data_tmp/gmx/run1 -top md.tpr -traj md_center.xtc -c A Z -np 2 --backend cuda'.split(' '))
    
    main()