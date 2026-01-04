import argparse
import os
import shutil
from pathlib import Path
import time
from typing import Dict, List, Tuple, Union

if 'MBAPY_PLT_AGG' in os.environ:
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lazydock.algorithm.utils import vectorized_sliding_average
from lazydock.gmx.mda.convert import PDBConverter, FakeAtomGroup
from lazydock.gmx.run import Gromacs
from lazydock.pml.interaction_utils import calcu_pdbstr_interaction
from lazydock.pml.plip_interaction import check_support_mode, run_plip_analysis
from lazydock.pml.rrcs import calcu_RRCS_from_array, calcu_RRCS_from_tensor
from lazydock.pml.thirdparty.modevectors import modevectors
from lazydock.scripts._script_utils_ import (Command, check_file_num_paried,
                                             excute_command,
                                             process_batch_dir_lst)
from lazydock.scripts.ana_interaction import (plip_mode, pml_mode,
                                              simple_analysis)
from matplotlib.ticker import FuncFormatter
from mbapy_lite.plot import save_show
from mbapy_lite.base import put_err, put_log, split_list
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.web import TaskPool
from MDAnalysis import AtomGroup, Universe
from scipy.stats import gaussian_kde
from tqdm import tqdm


class trjconv(Command):
    HELP = """"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--main-name', type=str, default='md.tpr',
                          help='main name in each sub-directory, such as md.tpr, default is %(default)s.')
        args.add_argument('-g', '--groups', type=str, nargs='+', default=['1', '0'],
                          help='groups for gmx trjconv, default is %(default)s.')
        args.add_argument('-ndx', '--index', type=str, default=None,
                          help='index file name in each sub-directory, such as tc_index.ndx, default is %(default)s.')
        args.add_argument('-pbc', type=str, default='mol', choices=['mol', 'atom', 'res', 'whole', 'cluster', 'nojump'],
                          help='pbc option for gmx trjconv, default is %(default)s.')
        args.add_argument('-ur', type=str, default='compact', choices=['rect', 'tric', 'compact'],
                          help='ur option for gmx trjconv, default is %(default)s.')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='number of workers to parallel. Default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='force to re-run the analysis, default is %(default)s.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='delete the exist analysis result, default is %(default)s.')

    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        
    def run_gmx_cmd(self, working_dir, *args, **kwargs):
        gmx = Gromacs(working_dir=working_dir)
        gmx.run_gmx_with_expect(*args, **kwargs)
        
    def main_process(self):
        # get complex paths
        complexs_path = get_paths_with_extension(self.args.batch_dir, [], name_substr=self.args.main_name)
        put_log(f'get {len(complexs_path)} task(s)')
        pool = TaskPool('threads', self.args.n_workers).start()
        exp_acts = []
        for g in self.args.groups:
            exp_acts.append({'Select a group:': f'{g}\r', '\\timeout': f'{g}\r'})
        # process each complex
        for complex_path in tqdm(complexs_path, total=len(complexs_path)):
            complex_path = Path(complex_path).resolve()
            main_name = complex_path.stem
            # check trjconv result exists
            if (complex_path.parent / f'{main_name}_center.xtc').exists():
                if self.args.delete:
                    os.remove(str(complex_path.parent / f'{main_name}_center.xtc'))
                    put_log(f'{main_name}_center.xtc deleted.')
                elif not self.args.force:
                    put_log(f'{main_name}_center.xtc already exists, skip trjconv.')
                    continue
            # perform trjconv, to avoid threads data conflict, just pass a nameless instance to pool
            pool.add_task(str(complex_path), self.run_gmx_cmd, complex_path.parent,
                          'trjconv', s=f'{main_name}.tpr', f=f'{main_name}.xtc', o=f'{main_name}_center.xtc',
                          n=self.args.index, pbc=self.args.pbc, ur=self.args.ur, center=True,
                          expect_actions=exp_acts, expect_settings={'timeout': 10})
            time.sleep(5) # delay for 5 seconds, avoid `expect` conflicts
            pool.wait_till(lambda: pool.count_waiting_tasks() == 0, 1)
        pool.close(1)


class make_ndx(trjconv):
    HELP = """"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-f', '--main-name', type=str, default='md.tpr',
                          help='main name in each sub-directory, such as md.tpr, default is %(default)s.')
        args.add_argument('-g', '--groups', type=str, nargs='+', default=['1', '0'],
                          help='groups option for make_ndx, will send each string to gmx make_ndx and automatically add a "q" in the last, default is %(default)s.')
        args.add_argument('-o', '--output', type=str, default='ana_index.ndx',
                          help='output index file name in each sub-directory, such as ana_index.ndx, default is %(default)s.')
        args.add_argument('-n', '--index', type=str, default=None,
                          help='index file name in each sub-directory, such as tc_index.ndx, default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='force to re-run the analysis, default is %(default)s.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='delete the exist analysis result, default is %(default)s.')
        
    def main_process(self):
        # get complex paths
        complexs_path = get_paths_with_extension(self.args.batch_dir, [], name_substr=self.args.main_name)
        put_log(f'get {len(complexs_path)} task(s)')
        # process each complex
        for complex_path in tqdm(complexs_path, total=len(complexs_path)):
            complex_path = Path(complex_path).resolve()
            gmx = Gromacs(working_dir=str(complex_path.parent))
            main_name = complex_path.stem
            # check result exists
            if os.path.exists(os.path.join(gmx.working_dir, self.args.output)):
                if self.args.delete:
                    os.remove(os.path.join(gmx.working_dir, self.args.output))
                    put_log(f'{self.args.output} deleted.')
                elif not self.args.force:
                    put_log(f'{self.args.output} already exists, skip.')
                    continue
            # perform make_ndx
            exp_acts = []
            for g in self.args.groups:
                exp_acts.append({'>': f'{g}\r', '\\timeout': f'{g}\r'})
            exp_acts.append({'>': 'q\r'})
            gmx.run_gmx_with_expect('make_ndx', f=f'{main_name}.tpr', n=self.args.index, o=self.args.output,
                                    expect_actions=exp_acts, expect_settings={'timeout': 5})


class simple(trjconv):
    HELP = """
    simple analysis for GROMACS simulation
    
    1. gmx_mpi rms -s md.tpr -f md_center.xtc -o rmsd.xvg -tu ns 
    2. gmx_mpi rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg
    3. gmx_mpi gyrate -s md.tpr -f md_center.xtc -o gyrate.xvg
    4. gmx_mpi hbond -s md.tpr -f md_center.xtc -num -dt 10
    
    5. gmx_mpi sasa -s md.tpr -f md_center.xtc -o sasa_total.xvg -or sasa_res.xvg -tu ns 
    6. gmx_mpi covar -s md.tpr -f md_center.xtc -o eigenval.xvg -tu ns 
    
    7. free energy landscape from rmsd and gyrate by MD-DaVis
    8. Probability Density Function from rmsd and gyrate
    """
    SUPPORT_METHODS = ['rms', 'rmsf', 'gyrate', 'hbond', 'sasa', 'covar', 'dssp', 'FEL', 'PDF']
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--main-name', type = str, default='md.tpr',
                          help='main name in each sub-directory, default is %(default)s.')
        args.add_argument('-t', '--main-type', type = str, nargs='+', default=[],
                          help='main name in each sub-directory, default is %(default)s.')
        args.add_argument('-ndx', '--index', type=str, default=None,
                          help='index file name in each sub-directory, such as ana_index.ndx, default is %(default)s.')
        args.add_argument('--methods', type = str, nargs='+', default=simple.SUPPORT_METHODS, choices=simple.SUPPORT_METHODS,
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('--dit-style', type = str, default=None,
                          help='DIT.mplstyle style file path, default is %(default)s.')
        args.add_argument('-rg', '--rms-group', type = str, default='4',
                          help='group to calculate rmsd, rmsf, and gyrate, default is %(default)s.')
        args.add_argument('-hg', '--hbond-group', type = int, nargs='+', default=[1, 1],
                          help='group to calculate hbond, default is %(default)s.')
        args.add_argument('-sg', '--sasa-group', type = str, default='4',
                          help='group to calculate sasa, default is %(default)s.')
        args.add_argument('-eg', '--eigenval-group', type = str, default='4',
                          help='group to calculate eigenval, default is %(default)s.')
        args.add_argument('-dg', '--dssp-group', type = str, default='1',
                          help='group to calculate DSSP, default is %(default)s.')
        args.add_argument('--dssp-num', action='store_true', default=False,
                          help='wheter to calculate DSSP number, default is %(default)s.')
        args.add_argument('--dssp-clear', action='store_true', default=False,
                          help='wheter to send --clear arg to gmx dssp, default is %(default)s.')
        args.add_argument('-xmax', '--eigenval-xmax', type = int, default=15,
                          help='max value of eigenval, default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='force to re-run the analysis, default is %(default)s.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='delete the exist analysis result, default is %(default)s.')
        args.add_argument('--task-suffix', type = str, default='',
                          help='suffix of task, default is %(default)s.')
        return args
        
    @staticmethod
    def rms(gmx: Gromacs, main_name: str, index: str = None, group: str = '4', force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_rmsd{gmx.task_uid}.csv')):
            if delete:
                (gmx.wdir / f'{main_name}_rmsd.csv').unlink(missing_ok=True)
                (gmx.wdir / f'rmsd.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'rmsd.png').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_rmsd.csv already exists, skip.')
        gmx.run_gmx_with_expect('rms', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'rmsd{gmx.task_uid}.xvg', tu='ns', n=index,
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_cmd_with_expect(f'dit xvg_compare -c 1 -f rmsd{gmx.task_uid}.xvg -o rmsd{gmx.task_uid}.png -smv -ws 10 -t "RMSD of {main_name}" -csv {main_name}_rmsd{gmx.task_uid}.csv -ns')
        
    @staticmethod
    def rmsf(gmx: Gromacs, main_name: str, index: str = None, group: str = '4', res: bool = True, force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_rmsf{gmx.task_uid}.csv')):
            if delete:
                (gmx.wdir / f'{main_name}_rmsf.csv').unlink(missing_ok=True)
                (gmx.wdir / f'rmsf.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'rmsf.png').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_rmsf.csv already exists, skip.')
        gmx.run_gmx_with_expect('rmsf', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'rmsf{gmx.task_uid}.xvg', res=res, n=index,
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_cmd_with_expect(f'dit xvg_compare -c 1 -f rmsf{gmx.task_uid}.xvg -o rmsf{gmx.task_uid}.png -t "RMSF of {main_name}" -csv {main_name}_rmsf{gmx.task_uid}.csv -ns')
        
    @staticmethod
    def gyrate(gmx: Gromacs, main_name: str, index: str = None, group: str = '4', force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_gyrate{gmx.task_uid}.csv')):
            if delete:
                (gmx.wdir / f'{main_name}_gyrate.csv').unlink(missing_ok=True)
                (gmx.wdir / f'gyrate.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'gyrate.png').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_gyrate.csv already exists, skip.')
        gmx.run_gmx_with_expect('gyrate', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'gyrate{gmx.task_uid}.xvg', n=index,
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_cmd_with_expect(f'dit xvg_compare -c 1 -f gyrate{gmx.task_uid}.xvg -o gyrate{gmx.task_uid}.png -smv -ws 10 -t "Gyrate of {main_name}" -csv {main_name}_gyrate{gmx.task_uid}.csv -ns')
        
    @staticmethod
    def hbond(gmx: Gromacs, main_name: str, index: str = None, group: Tuple[int, int] = (1, 1), dt=10, force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_hbond_num{gmx.task_uid}.csv')):
            if delete:
                (gmx.wdir / f'{main_name}_hbond_dist.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_hbond_num.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_hbond_num.csv').unlink(missing_ok=True)
                (gmx.wdir / f'hbond_dist.png').unlink(missing_ok=True)
                (gmx.wdir / f'hbond_num.png').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_hbond_num.csv already exists, skip.')
        gmx.run_gmx_with_expect('hbond', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc',
                                    num=f'{main_name}_hbond_num{gmx.task_uid}.xvg', dist=f'{main_name}_hbond_dist{gmx.task_uid}.xvg', n=index,
                                    expect_actions=[{'Select a group:': f'{group[0]}\r', '\\timeout': f'{group[0]}\r'},
                                                    {'Select a group:': f'{group[1]}\r', '\\timeout': f'{group[1]}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_cmd_with_expect(f'dit xvg_compare -c 1 -f {main_name}_hbond_num{gmx.task_uid}.xvg -o hbond_num{gmx.task_uid}.png -smv -ws 10 -t "H-bond num of {main_name}" -csv {main_name}_hbond_num{gmx.task_uid}.csv -ns')
        gmx.run_cmd_with_expect(f'dit xvg_show -f {main_name}_hbond_dist{gmx.task_uid}.xvg -o hbond_dist{gmx.task_uid}.png -ns')

    @staticmethod
    def sasa(gmx: Gromacs, main_name: str, index: str = None, group: str = '4', force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_sasa_tv{gmx.task_uid}.csv')):
            if delete:
                for ty in ['total', 'res', 'dg', 'tv']:
                    (gmx.wdir / f'{main_name}_sasa_{ty}{gmx.task_uid}.xvg').unlink(missing_ok=True)
                    (gmx.wdir / f'{main_name}_sasa_{ty}{gmx.task_uid}.png').unlink(missing_ok=True)
                    (gmx.wdir / f'{main_name}_sasa_{ty}{gmx.task_uid}.csv').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_sasa_tv.csv already exists, skip.')
        gmx.run_gmx_with_expect(f'sasa -or {main_name}_sasa_res{gmx.task_uid}.xvg', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc',
                                    o=f'{main_name}_sasa_total{gmx.task_uid}.xvg', odg=f'{main_name}_sasa_dg{gmx.task_uid}.xvg', tv=f'{main_name}_sasa_tv{gmx.task_uid}.xvg', tu='ns', n=index,
                                    expect_actions=[{'>': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        for ty in ['total', 'res', 'dg', 'tv']:
            gmx.run_cmd_with_expect(f'dit xvg_compare -c 1 -f {main_name}_sasa_{ty}{gmx.task_uid}.xvg -o {main_name}_sasa_{ty}{gmx.task_uid}.png -smv -ws 10 -t "SASA {ty} of {main_name}" -csv {main_name}_sasa_{ty}{gmx.task_uid}.csv -ns')

    @staticmethod
    def covar(gmx: Gromacs, main_name: str, index: str = None, group: str = '4', xmax: int = 15, force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_eigenval{gmx.task_uid}.csv')):
            if delete:
                (gmx.wdir / f'{main_name}_eigenval{gmx.task_uid}.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_eigenval{gmx.task_uid}.png').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_eigenval{gmx.task_uid}.csv').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_eigenval.csv already exists, skip.')
        gmx.run_gmx_with_expect('covar', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'{main_name}_eigenval{gmx.task_uid}.xvg', tu='ns', n=index,
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        gmx.run_cmd_with_expect(f'dit xvg_compare -c 1 -f {main_name}_eigenval{gmx.task_uid}.xvg -o {main_name}_eigenval{gmx.task_uid}.png -xmin 0 -xmax {xmax} -t "Eigenval of {main_name}" -csv {main_name}_eigenval{gmx.task_uid}.csv -ns')
    
    @staticmethod
    def dssp(gmx: Gromacs, main_name: str, index: str = None, group: str = None, num: bool = False, clear: bool = False,
             force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_dssp_mat{gmx.task_uid}.dat')):
            if delete:
                (gmx.wdir / f'{main_name}_dssp_num{gmx.task_uid}.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_dssp_mat{gmx.task_uid}.dat').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_dssp_mat{gmx.task_uid}.xpm').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_dssp_mat{gmx.task_uid}.png').unlink(missing_ok=True)
                (gmx.wdir / f'{main_name}_dssp_num{gmx.task_uid}.png').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_dssp_mat{gmx.task_uid}.dat already exists, skip.')
        kwgs = {}
        if num:
            kwgs = {'num': f'{main_name}_dssp_num{gmx.task_uid}.xvg'}
        if clear:
            kwgs.update({'_clear': True})
        gmx.run_gmx_with_expect('dssp', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'{main_name}_dssp_mat{gmx.task_uid}.dat',
                                sel=group, n=index, _hmode='dssp', tu='ns', **kwgs)
        gmx.run_cmd_with_expect(f'dit dssp -f {main_name}_dssp_mat{gmx.task_uid}.dat -o {main_name}_dssp_mat{gmx.task_uid}.xpm')
        gmx.run_cmd_with_expect(f'dit xpm_show -f {main_name}_dssp_mat{gmx.task_uid}.xpm -o {main_name}_dssp_mat{gmx.task_uid}.png -xs 0.01 --x_precision 0 -x "Time (ns)" -y "Residues (aa)" -ns')
        if num:
            gmx.run_cmd_with_expect(f'dit xvg_compare -c 1-10 -f {main_name}_dssp_num{gmx.task_uid}.xvg -o {main_name}_dssp_num{gmx.task_uid}.png -t "DSSP number of {main_name}" -csv {main_name}_dssp_num{gmx.task_uid}.csv -ns')
    
    @staticmethod
    def free_energy_landscape(gmx: Gromacs, main_name: str, force: bool = False, delete: bool = False, **kwargs):
        # MD-DaVis
        if os.path.exists(os.path.join(gmx.working_dir, f'FEL{gmx.task_uid}.html')):
            if delete:
                (gmx.wdir / f'FEL{gmx.task_uid}.html').unlink(missing_ok=True)
            if not force:
                return put_log(f'FEL{gmx.task_uid}.html already exists, skip.')
        gmx.run_cmd_with_expect(f'md-davis landscape_xvg -c -T 300 -x rmsd.xvg -y gyrate.xvg -o FEL{gmx.task_uid}.html -n FEL -l "RMSD-Rg" --axis_labels "dict(x=\'RMSD (in nm)\', y=\'Rg (in nm)\', z=\'Free Energy (kJ mol<sup>-1</sup>)<br>\')"')
        # gmx and dit
        if os.path.exists(os.path.join(gmx.working_dir, f'rmsd_gyrate{gmx.task_uid}.png')):
            if delete:
                (gmx.wdir / f'rmsd_gyrate{gmx.task_uid}.png').unlink(missing_ok=True)
                (gmx.wdir / f'rmsd_gyrate{gmx.task_uid}.xvg').unlink(missing_ok=True)
                (gmx.wdir / f'sham{gmx.task_uid}.xpm').unlink(missing_ok=True)
            if not force:
                return put_log(f'rmsd_gyrate{gmx.task_uid}.png already exists, skip.')
        gmx.run_cmd_with_expect(f'dit xvg_combine -f rmsd{gmx.task_uid}.xvg gyrate{gmx.task_uid}.xvg -c 0,1 1 -l RMSD Gyrate -o rmsd_gyrate{gmx.task_uid}.xvg -x "Time (ps)"')
        gmx.run_gmx_with_expect(f'sham -f rmsd_gyrate{gmx.task_uid}.xvg -ls sham{gmx.task_uid}.xpm --ngrid 100')
        gmx.run_cmd_with_expect(f'dit xpm_show -f sham{gmx.task_uid}.xpm --x_precision 2 --y_precision 2 -x "RMSD (nm)" -y "Rg (nm)" -cmap jet --colorbar_location right -ip gaussian -o rmsd_gyrate{gmx.task_uid}.png -ns')
        
    @staticmethod
    def plot_PDF(gmx: Gromacs, main_name: str, force: bool = False, delete: bool = False, **kwargs):
        if os.path.exists(os.path.join(gmx.working_dir, f'{main_name}_PDF.png')):
            if delete:
                (gmx.wdir / f'{main_name}_PDF.png').unlink(missing_ok=True)
            if not force:
                return put_log(f'{main_name}_PDF.png already exists, skip rms.')
        """idea from https://pymolwiki.org/index.php/Geo_Measures_Plugin"""
        # read data and calculate density
        if not os.path.exists(f'{gmx.working_dir}/{main_name}_rmsd.csv') or not os.path.exists(f'{gmx.working_dir}/{main_name}_gyrate.csv'):
            return put_log(f'{main_name}_rmsd.csv or {main_name}_gyrate.csv not exists, skip rms.')
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
        plt.xlabel('RMSD (in nm)', fontsize=16, weight='bold')
        plt.ylabel('Rg (in nm)', fontsize=16, weight='bold')
        ax.tick_params(labelsize=14)
        # Adding the color bar
        cbar = plt.colorbar(pdf)
        cbar.set_label("Probability Density Function", fontsize=16, weight='bold')
        cbar.ax.tick_params(labelsize=14)
        save_show(os.path.join(gmx.working_dir, f'{main_name}_PDF.png'), 600, show=False)
        plt.close(fig)
        
    def main_process(self):
        # get complex paths
        complexs_path = get_paths_with_extension(self.args.batch_dir, self.args.main_type, name_substr=self.args.main_name)
        put_log(f'get {len(complexs_path)} task(s)')
        # process each complex
        for complex_path in tqdm(complexs_path, total=len(complexs_path)):
            complex_path = Path(complex_path).resolve()
            gmx = Gromacs(working_dir=str(complex_path.parent))
            gmx.task_uid = self.args.task_suffix
            main_name = complex_path.stem
            if (complex_path.parent / f'{main_name}.tpr').exists() and (complex_path.parent / f'{main_name}_center.xtc').exists():
                put_log(f'Perform analysis for {main_name}.tpr and {main_name}_center.xtc.')
            else:
                put_err(f'{main_name}.tpr or {main_name}_center.xtc not exists in {complex_path.parent}, skip.')
                continue
            # copy DIT.mplstyle file to working directory
            if self.args.dit_style and os.path.exists(self.args.dit_style):
                shutil.copy(self.args.dit_style, str(complex_path.parent))
            # perform analysis
            if 'rms' in self.args.methods:
                self.rms(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.rms_group,
                        force=self.args.force, delete=self.args.delete)
            if 'rmsf' in self.args.methods:
                self.rmsf(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.rms_group,
                        force=self.args.force, delete=self.args.delete)
            if 'gyrate' in self.args.methods:
                self.gyrate(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.rms_group,
                            force=self.args.force, delete=self.args.delete)
            if 'hbond' in self.args.methods:
                self.hbond(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.hbond_group,
                        force=self.args.force, delete=self.args.delete)
            if 'sasa' in self.args.methods:
                self.sasa(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.sasa_group,
                        force=self.args.force, delete=self.args.delete)
            if 'covar' in self.args.methods:
                self.covar(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.eigenval_group,
                        xmax=self.args.eigenval_xmax, force=self.args.force, delete=self.args.delete)
            if 'dssp' in self.args.methods:
                self.dssp(gmx, main_name=complex_path.stem, index=self.args.index, group=self.args.dssp_group,
                        num=self.args.dssp_num, clear=self.args.dssp_clear, force=self.args.force, delete=self.args.delete)
            # perform free energy landscape by MD-DaVis
            if 'FEL' in self.args.methods:
                self.free_energy_landscape(gmx, main_name=complex_path.stem, force=self.args.force, delete=self.args.delete)
            # plot PDF
            if 'PDF' in self.args.methods:
                self.plot_PDF(gmx, main_name=complex_path.stem, force=self.args.force, delete=self.args.delete)
            
            
class mmpbsa(simple):
    HELP = """
    mmpbsa analysis for GROMACS simulation
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser, mmpbsa_args: bool = True):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        if mmpbsa_args:
            args.add_argument('-i', '--input', type = str, required=True,
                              help=f"gmx_MMPBSA input file name in each sub-folder, such as mmpbsa.in")
            args.add_argument('-o', '--output', type = str, default='MMPBSA_FINAL_RESULTS',
                              help=f"gmx_MMPBSA output file name, such as MMPBSA_FINAL_RESULTS")
            args.add_argument('-np', '--np', type = int, required=True,
                              help=f"npi np argument for gmx_MMPBSA")
        args.add_argument('-top', '--top-name', type = str, default='md.tpr',
                          help="topology file name in each sub-folder, default is %(default)s.")
        args.add_argument('-traj', '--traj-name', type = str, default='md_center.xtc',
                          help="trajectory file name in each sub-folder, default is %(default)s.")
        args.add_argument('--receptor-chain-name', type = str, required=True,
                          help='receptor chain name, such as "A".')
        args.add_argument('--ligand-chain-name', type = str, required=True,
                          help='ligand chain name, such as "LIG".')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='force to re-run the analysis, default is %(default)s.')
        return args
        
    def get_complex_atoms_index(self, u: Universe):
        rec_idx = u.atoms.chainIDs == self.args.receptor_chain_name
        lig_idx = u.atoms.chainIDs == self.args.ligand_chain_name
        put_log(f"receptor atoms: {rec_idx.sum()}, ligand atoms: {lig_idx.sum()}.")
        return rec_idx, lig_idx
    
    def get_index_range(self, idx: np.ndarray):
        return idx.argmax(), idx.shape[0] - idx[::-1].argmax()
    
    def check_top_traj(self, bdir = None):
        bdir = bdir or self.args.batch_dir
        top_paths = get_paths_with_extension(bdir, [os.path.split(self.args.top_name)[-1]], name_substr=self.args.top_name)
        traj_paths = get_paths_with_extension(bdir, [os.path.split(self.args.traj_name)[-1]], name_substr=self.args.traj_name)
        invalid_roots = check_file_num_paried(top_paths, traj_paths)
        if invalid_roots:
            put_err(f"The number of top and traj files is not equal, please check the input files.\ninvalid roots:{invalid_roots}", _exit=True)
        return top_paths, traj_paths
        
    def find_tasks(self):
        tasks = []
        for r_path, l_path in zip(self.top_paths, self.traj_paths):
            tasks.append((r_path, l_path))
        return tasks
    
    def main_process(self):
        # load origin dfs from data file
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # run tasks
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            wdir_repr = os.path.relpath(wdir, self.args.batch_dir) # relative path to batch_dir, shorter
            bar.set_description(f"{wdir_repr}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            # check results
            if os.path.exists(os.path.join(wdir, self.args.output+'.csv')) and not self.args.force:
                put_log(f"{self.args.output}.csv already exists, skip.")
                bar.update(1)
                continue
            # get receptor and ligand atoms index range
            u = Universe(top_path, traj_path)
            rec_idx, lig_idx = self.get_complex_atoms_index(u)
            rec_min, rec_max = self.get_index_range(rec_idx)
            lig_min, lig_max = self.get_index_range(lig_idx)
            rec_range_str, lig_range_str = f"{rec_min+1}-{rec_max}", f"{lig_min+1}-{lig_max}"
            # make index file for receptor and ligand
            gmx = Gromacs(working_dir=wdir)
            gmx.run_gmx_with_expect('make_ndx', f=os.path.basename(top_path), o='mmpbsa_tmp.ndx',
                                    expect_actions=[{'>': 'q\r'}])
            sum_groups = opts_file(os.path.join(gmx.working_dir, 'mmpbsa_tmp.ndx')).count(']')
            gmx.run_gmx_with_expect('make_ndx', f=os.path.basename(top_path), o='mmpbsa.ndx',
                                    expect_actions=[{'>': f'a {rec_range_str}\r'}, {'>': f'name {sum_groups} MMPBSA_Receptor\r'},
                                                    {'>': f'a {lig_range_str}\r'}, {'>': f'name {sum_groups+1} MMPBSA_Ligand\r'},
                                                    {'>': 'q\r'}])
            # check MMPBSA parameters input file
            if not os.path.exists(os.path.join(wdir, self.args.input)):
                if os.path.exists(self.args.input):
                    input_name = os.path.basename(self.args.input)
                    shutil.copy(self.args.input, os.path.join(wdir, input_name))
                else:
                    put_err(f"input file {self.args.input} not exists, skip.")
                    continue
            else:
                input_name = self.args.input
            # call gmx_MMPBSA
            cmd_str = f'gmx_MMPBSA -O -i {input_name} -cs {self.args.top_name} -ct {self.args.traj_name} -ci mmpbsa.ndx -cg {sum_groups} {sum_groups+1} -cp topol.top -o {self.args.output}.dat -eo {self.args.output}.csv -nogui'
            os.system(f'cd "{gmx.working_dir}" && mpirun -np {self.args.np} {cmd_str}')
            bar.update(1)
    
    
def run_pdbstr_interaction_analysis(fake_ag: FakeAtomGroup, receptor_chain: str, ligand_chain: str,
                                    method: str, mode: str, cutoff: float, hydrogen_atom_only: bool,
                                    alter_chain: Dict[str,str] = None, alter_res: Dict[str,str] = None, alter_atm: Dict[str,str] = None):
    pdbstr = PDBConverter(fake_ag).fast_convert(alter_chain=alter_chain, alter_res=alter_res, alter_atm=alter_atm)
    if method == 'pymol':
        inter = calcu_pdbstr_interaction(f'chain {receptor_chain}', f'chain {ligand_chain}', pdbstr, mode, cutoff, hydrogen_atom_only)
    elif method == 'plip':
        mode = check_support_mode(mode)
        inter = run_plip_analysis(pdbstr, receptor_chain, ligand_chain, mode, cutoff)
    else:
        return put_err(f"method {method} not supported, return None.")
    return inter


class interaction(simple_analysis, mmpbsa):
    HELP = """
    interaction analysis for GROMACS simulation
    """
    def __init__(self, args, printf=print):
        Command.__init__(self, args, printf, ['batch_dir'])
        self.alter_chain = {}
        self.alter_res = None
        self.alter_atm = None

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        mmpbsa.make_args(args, mmpbsa_args=False)
        args.add_argument('-gro', '--gro-name', type = str, default='md.gro',
                          help=f"gro file name in each sub-folder.")
        args.add_argument('--suffix', type = str, default='',
                          help='suffix for output file name, default is %(default)s.')
        args.add_argument('--alter-receptor-chain', type = str, default=None,
                          help='alter receptor chain name from topology to user-define, such as "A".')
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
        args.add_argument('--max-plot', type=int, default=None,
                          help='max res to plot, filter by (max, mean) interaction freq. Default is %(default)s.')
        args.add_argument('--skip-plot', action='store_true', default=False,
                          help='skip plot. Default is %(default)s.')
        args.add_argument('--ref-res', type = str, default='',
                          help='reference residue name, input string shuld be like GLY300,ASP330, also support a text file contains this format string as a line.')
        args.add_argument('-nw', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        args.add_argument('-b', '--begin-frame', type=int, default=1,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-e', '--end-frame', type=int, default=None,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-step', '--traj-step', type=int, default=1,
                          help='Step while reading trajectory. Default is %(default)s.')
        args.add_argument('--plot-time-unit', type=int, default=100,
                          help='time unit for plot in X-Axis, default is %(default)s.')
        args.add_argument('--yticks-interval', type=int, default=10,
                          help='interval for y axis ticks, default is %(default)s.')
        args.add_argument('--fig-size', type=int, nargs=2, default=[9, 6],
                          help='figure size, default is %(default)s.')
    
    def process_args(self):
        # self.args.alter_ligand_chain will passed to final interaction calcu function
        if self.args.alter_ligand_chain is not None:
            self.alter_chain[self.args.ligand_chain_name] = self.args.alter_ligand_chain
        else:
            self.args.alter_ligand_chain = self.args.ligand_chain_name
        # in lazydock, default ligand res name is the chain name too, so alter chain name to alter-ligand-res
        if self.args.alter_ligand_res is not None:
            self.alter_res = {self.args.ligand_chain_name: self.args.alter_ligand_res}
        # set self.alter_atm
        if self.args.alter_ligand_atm is not None:
            self.alter_atm = {self.args.alter_ligand_chain: self.args.alter_ligand_atm}
        # set alter receptor chain
        if self.args.alter_receptor_chain is not None:
            self.alter_chain[self.args.receptor_chain_name] = self.args.alter_receptor_chain
        else:
            self.args.alter_receptor_chain = self.args.receptor_chain_name
        # output formater
        self.output_formater = getattr(self, f'output_fromater_{self.args.output_style}')
        # check batch dir， check method and mode AND load ref_res
        simple_analysis.process_args(self)
        
    def calcu_interaction(self, top_path: str, gro_path: str, traj_path: str, pool: TaskPool):
        # load pdbstr from traj
        u, u2 = Universe(top_path, traj_path), Universe(gro_path)
        u.atoms.residues.resids = u2.atoms.residues.resids
        rec_idx, lig_idx = self.get_complex_atoms_index(u)
        if rec_idx.sum() == 0 or lig_idx.sum() == 0:
            return put_err(f"no atoms found in receptor or ligand, skip.", (None, None))
        complex_ag = u.atoms[rec_idx | lig_idx]
        # calcu interaction for each frame
        sum_frames = (len(u.trajectory) if self.args.end_frame is None else self.args.end_frame) - self.args.begin_frame
        for frame in tqdm(u.trajectory[self.args.begin_frame:self.args.end_frame:self.args.traj_step],
                        total=sum_frames//self.args.traj_step, desc='Calculating frames', leave=False):
            fake_ag = FakeAtomGroup(complex_ag)
            pool.add_task(frame.time, run_pdbstr_interaction_analysis, fake_ag,
                        self.args.alter_receptor_chain, self.args.alter_ligand_chain,
                        self.args.method, self.args.mode, self.args.cutoff, self.args.hydrogen_atom_only,
                        alter_chain=self.alter_chain, alter_res=self.alter_res, alter_atm=self.alter_atm)
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
        return interactions, df
        
    def main_process(self):
        # load origin dfs from data file
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        suffix = self.args.suffix
        print(f'find {len(self.tasks)} tasks.')
        # run tasks
        pool = TaskPool('process', self.args.n_workers).start()
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            wdir_repr = os.path.relpath(wdir, self.args.batch_dir) # relative path to batch_dir, shorter
            bar.set_description(f"{wdir_repr}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            # calcu interaction and save to file OR load results if have been calculated before and not force recalculate
            top_path = Path(top_path).resolve()
            gro_path = str(top_path.parent / self.args.gro_name)
            csv_path = str(top_path.parent / f'{top_path.stem}_{self.args.method}_interactions{suffix}.csv')
            pkl_path = str(top_path.parent / f'{top_path.stem}_{self.args.method}_interactions{suffix}.pkl')
            if (not os.path.exists(csv_path) or not os.path.exists(pkl_path)) or self.args.force:
                interactions, df = self.calcu_interaction(str(top_path), gro_path, traj_path, pool)
                if interactions is None:
                    put_log(f'no interaction found, skip.')
                    bar.update(1)
                    continue
                df.to_csv(csv_path, index=False)
                opts_file(pkl_path, 'wb', way='pkl', data=interactions)
            else:
                interactions = opts_file(pkl_path, 'rb', way='pkl')
                # ckeck whether to plot
                if self.args.skip_plot:
                    pool.task = {}
                    bar.update(1)
                    continue
            # transform interaction df to plot matrix df
            times = split_list(list(interactions.keys()), self.args.plot_time_unit)
            plot_df = pd.DataFrame()
            for i, time_u in tqdm(enumerate(times), desc='Gathering interactions', total=len(times)):
                for time_i in time_u:
                    lst = [single_v for v in interactions[time_i].values() if v for single_v in v if single_v]
                    for single_inter in lst:
                        receptor_res = f'{single_inter[0][1]}{single_inter[0][0]}'
                        if i not in plot_df.index or receptor_res not in plot_df.columns:
                            plot_df.loc[i, receptor_res] = 1
                        elif np.isnan(plot_df.loc[i, receptor_res]):
                            plot_df.loc[i, receptor_res] = 1
                        else:
                            plot_df.loc[i, receptor_res] += 1
                if i not in plot_df.index:
                    plot_df.loc[i, :] = np.nan
                plot_df.loc[i, :] /= len(time_u)
            # filter plot df by (max, mean) inter value
            if self.args.max_plot is not None and len(plot_df.columns) > self.args.max_plot:
                sort_val = {k: (v.max(), v.mean()) for k, v in plot_df.items()}
                sorted_val = sorted(list(sort_val.keys()), key=lambda k: sort_val[k])
                to_del = sorted_val[:len(sort_val)-self.args.max_plot]
                put_log(f'delete {to_del} from plot_df')
                plot_df.drop(to_del, axis=1, inplace=True)
            # sort residue index
            plot_df = plot_df[sorted(list(plot_df.columns), key=lambda x: int(x[3:]))]
            # save to csv and plot
            plot_df.to_csv(str(top_path.parent / f'{top_path.stem}_{self.args.method}_plot_df.csv'), index=False)
            if not plot_df.empty:
                fig, ax = plt.subplots(figsize=self.args.fig_size)
                sns.heatmap(plot_df, xticklabels=list(plot_df.columns),
                            cmap='viridis', cbar_kws={'label': 'Interaction frequency'}, ax=ax)
                y_ticks = list(range(0, len(plot_df)+1, self.args.yticks_interval))
                ax.set_yticks(y_ticks, list(map(str, y_ticks)))
                ax.tick_params(labelsize=14, axis='both')
                plt.xlabel('Residues (aa)', fontsize=16, weight='bold')
                plt.ylabel('Time (ns)', fontsize=16, weight='bold')
                cbar = ax.collections[0].colorbar
                cbar.ax.tick_params(labelsize=14)
                cbar.ax.set_ylabel('Interaction frequency', fontsize=16)
                save_show(str(top_path.parent / f'{top_path.stem}_{self.args.method}_interactions{suffix}.png'), 600, show=False)
                plt.close(fig)
            # other things
            pool.clear()
            bar.update(1)
        pool.close(timeout=1)


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
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('-top', '--top-name', type = str, default='md.tpr',
                          help="topology file name in each sub-folder, default is %(default)s.")
        args.add_argument('-gro', '--gro-name', type = str, default='md.gro',
                          help="gro file name in each sub-folder, default is %(default)s.")
        args.add_argument('-traj', '--traj-name', type = str, default='md_center.xtc',
                          help="trajectory file name in each sub-folder, default is %(default)s.")
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
        
    @staticmethod
    def plot_average_heatmap(scores: np.ndarray, top_path: Path):
        plt.imshow(scores.mean(axis=0), cmap='viridis')
        plt.xlabel('Residue Index', fontsize=16, weight='bold')
        plt.ylabel('Residue Index', fontsize=16, weight='bold')
        plt.gca().tick_params(labelsize=14)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.set_ylabel('Average RRCS', fontsize=16, weight='bold')
        save_show(str(top_path.parent / f'{top_path.stem}_RRCS.png'), 600, show=False)
        plt.close()
        
        
    @staticmethod
    def plot_diag_vs_frame(scores: np.ndarray, top_path: Path):
        # determine the n
        diags = np.stack([np.diag(score, k=-1) for score in scores], axis=0).T
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.heatmap(diags, cmap='viridis', cbar_kws={'label': r'RRCS'}, xticklabels=1000, ax=ax)
        ax.set_aspect('auto')
        plt.xlabel('Frames', fontsize=16, weight='bold')
        plt.ylabel('Residue Index', fontsize=16, weight='bold')
        cbar = ax.collections[0].colorbar
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.set_ylabel('RRCS', fontsize=16, weight='bold')
        save_show(str(top_path.parent / f'{top_path.stem}_RRCS_vs_Frame.png'), 600, show=False)
        plt.close()
        
    @staticmethod
    def plot_lines(scores: np.ndarray, top_path: Path, ag_gro: AtomGroup):
        def _set_vertical_title(ax, text, pos):
            """设置垂直Y轴标题的辅助函数"""
            t = ax.yaxis.set_label_coords(*pos)
            ax.set_ylabel(text, 
                        rotation=90,
                        labelpad=25,
                        fontsize=16,
                        verticalalignment='center',
                        fontweight='bold')
        n_frames, N, _ = scores.shape
        
        # 创建3行1列的竖向布局，设置高度比例和画布尺寸
        fig = plt.figure(figsize=(9, 8))
        gs = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1, 1, 1], hspace=0)
        
        # 预定义颜色循环
        colors = plt.cm.tab20(np.linspace(0, 1, 10))
        
        # 公共参数设置
        plot_params = {
            'linewidth': 2.5,
            'alpha': 0.7,
        }
        style_params = {
            'xlabel': 'Time (ns)',
            'ylabel_pos': (-0.12, 0.5),  # Y轴标题位置参数
            'smooth_window': 100
        }

        # 第一个子图：全时段均值最大的top10
        ax1 = plt.subplot(gs[0])
        flat_full_mean = np.mean(scores, axis=0).ravel()
        top_indices = np.argpartition(flat_full_mean, -10)[-10:]
        i1, j1 = np.unravel_index(top_indices, (N, N))
        labels1 = [f'{ag_gro.residues.resnames[i]}{i}-{ag_gro.residues.resnames[j]}{j}' 
                for i, j in zip(i1, j1)]
        data1 = vectorized_sliding_average(scores[:, i1, j1], style_params['smooth_window'])
        for idx in range(10):
            ax1.plot(data1[:, idx], color=colors[idx], **plot_params)
        _set_vertical_title(ax1, "Top Average Contacts", style_params['ylabel_pos'])
        
        # 第二个子图：差异最大的top10（后增前减）
        ax2 = plt.subplot(gs[1], sharex=ax1)
        front_mean = np.mean(scores[:max(1, int(n_frames*0.25))], axis=0)
        back_mean = np.mean(scores[-max(1, int(n_frames*0.25)):], axis=0)
        diff = back_mean - front_mean
        top_diff_indices = np.argpartition(diff.ravel(), -10)[-10:]
        i2, j2 = np.unravel_index(top_diff_indices, (N, N))
        labels2 = [f'{ag_gro.residues.resnames[i]}{i}-{ag_gro.residues.resnames[j]}{j}' 
                for i, j in zip(i2, j2)]
        data2 = vectorized_sliding_average(scores[:, i2, j2], style_params['smooth_window'])
        for idx in range(10):
            ax2.plot(data2[:, idx], color=colors[idx], **plot_params)
        _set_vertical_title(ax2, "Increasing Contacts", style_params['ylabel_pos'])
        
        # 第三个子图：差异最大的top10（前增后减）
        ax3 = plt.subplot(gs[2], sharex=ax1)
        diff = front_mean - back_mean
        top_dec_indices = np.argpartition(diff.ravel(), -10)[-10:]
        i3, j3 = np.unravel_index(top_dec_indices, (N, N))
        labels3 = [f'{ag_gro.residues.resnames[i]}{i}-{ag_gro.residues.resnames[j]}{j}' 
                for i, j in zip(i3, j3)]
        data3 = vectorized_sliding_average(scores[:, i3, j3], style_params['smooth_window'])
        for idx in range(10):
            ax3.plot(data3[:, idx], color=colors[idx], **plot_params)
        _set_vertical_title(ax3, "Decreasing Contacts", style_params['ylabel_pos'])

        # 设置图例和XY轴范围
        max_val = max(np.max(data1), np.max(data2), np.max(data3))
        for ax, labels, top_info in zip([ax1, ax2, ax3], [labels1, labels2, labels3],
                                        ['Top 10 Average', 'Top 10 Increasing', 'Top 10 Decreasing']):
            ax.set_xlim((0, n_frames))
            ax.set_ylim((-5, max_val+5))
            leg = ax.legend(labels, title=f'{top_info} Residue Pairs', title_fontsize=16,
                            loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=14,
                            frameon=False, labelspacing=0.1, ncol=1)
            leg._legend_box.align = "left"
            for text, color in zip(leg.get_texts(), colors):
                text.set_color(color)
        
        # 调整坐标轴显示
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax3.set_xlabel(style_params['xlabel'], fontsize=16, weight='bold')
        ax3.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x/100:.0f}'))
        
        # 统一设置坐标轴样式
        for ax in [ax1, ax2, ax3]:
            ax.tick_params(axis='both', which='major', labelsize=14)
            ax.grid(True, linestyle='--', alpha=0.6)
        
        # 优化布局并保存
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # 为图例留出右侧空间
        save_show(str(top_path.parent / f'{top_path.stem}_RRCS_vertical.png'), 600, show=False)
        plt.close()

    def main_process(self):
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # run tasks
        pool = TaskPool('process', self.args.n_workers).start()
        bar = tqdm(total=len(self.tasks), desc='Calculating RRCS')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            top_path = Path(top_path).resolve()
            wdir_repr = os.path.relpath(wdir, self.args.batch_dir) # relative path to batch_dir, shorter
            bar.set_description(f"{wdir_repr}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            if os.path.exists(os.path.join(wdir, f'{top_path.stem}_RRCS.npz')) and not self.args.force:
                put_log(f'{top_path.stem}_RRCS.npz already exists, skip.')
                continue
            # load pdbstr from traj
            u = Universe(str(top_path), traj_path)
            u_gro = Universe(str(top_path.parent / self.args.gro_name))
            if self.args.chains is not None and len(self.args.chains):
                idx = u.atoms.chainIDs == self.args.chains[0]
                for chain_i in self.args.chains[1:]:
                    idx = idx | (u.atoms.chainIDs == chain_i)
                ag = u.atoms[idx]
                ag_gro = u_gro.atoms[idx]
            else:
                ag = u.atoms
                ag_gro = u_gro.atoms
            print(f'find {np.unique(u.atoms.chainIDs)} chains, {len(ag)} atoms in {self.args.chains}')
            # calcu interaction for each frame
            sum_frames = (len(u.trajectory) if self.args.end_frame is None else self.args.end_frame) - self.args.begin_frame
            for frame in tqdm(u.trajectory[self.args.begin_frame:self.args.end_frame:self.args.traj_step],
                              total=sum_frames//self.args.traj_step, desc='Calculating frames', leave=False):
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
            scores = np.stack(list(scores.values()), axis=0)
            # save result
            np.savez_compressed(str(top_path.parent / f'{top_path.stem}_RRCS.npz'),
                                scores=scores, frames=np.array(frames), resis=ag.resids,
                                resns=ag.resnames, chains=ag.chainIDs)
            # plot
            self.plot_average_heatmap(scores, top_path)
            self.plot_diag_vs_frame(scores, top_path)
            self.plot_lines(scores, top_path, ag_gro)
            # other things
            pool.task = {}
            bar.update(1)
        pool.close()


def run_porcupine(top_path: str, traj_path: str, start: int = 0, stop: int = None, step: int = 1):
    from pymol import cmd
    cmd.reinitialize()
    u = Universe(top_path, traj_path)
    sum_frames = len(u.trajectory)
    put_log(f'{traj_path} has {len(u.atoms)} atoms with {sum_frames} frames')
    del u
    # load top and traj
    stop = stop or sum_frames
    sum_frames = (stop - start)//step
    cmd.load(top_path, 'mol')
    cmd.load_traj(traj_path, 'mol', interval=step, start=start, stop=stop)
    cmd.split_states('mol')
    # run modevectors
    modevectors(f'mol_{1:0>4d}', f'mol_{sum_frames:0>4d}',
                cutoff=.0, head_length=0.3, head=0.2, headrgb="1.0,0.2,0.1", tailrgb="1.0,0.2,0.1", notail=0)
    cmd.set('cartoon_trace', 0)
    cmd.set('cartoon_tube_radius', 0.3)
    cmd.disable('all')
    cmd.enable(f'mol_{1:0>4d}', 1)
    cmd.enable('modevectors', 1)
    cmd.set('ray_shadow', 0)
    top_path: Path = Path(top_path)
    cmd.save(str(top_path.parent / f'{top_path.stem}_porcupine.pse'), f'mol_{1:0>4d} or modevectors', format='pse')


class porcupine(mmpbsa):
    HELP = """
    porcupine plot analysis for GROMACS simulation
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-top', '--top-name', type = str, default='md.tpr',
                          help="topology file name in each sub-folder, default is %(default)s.")
        args.add_argument('-traj', '--traj-name', type = str, default='md_center.xtc',
                          help="trajectory file name in each sub-folder, default is %(default)s.")
        args.add_argument('-b', '--begin-frame', type=int, default=0,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-e', '--end-frame', type=int, default=None,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-step', '--traj-step', type=int, default=1,
                          help='Step while reading trajectory. Default is %(default)s.')
        args.add_argument('-nw', '--n-workers', type=int, default=4,
                          help='number of workers to parallel. Default is %(default)s.')
        args.add_argument('-D', '--delete', default=False, action='store_true',
                          help='delete the exist analysis result, default is %(default)s.')
        
    def main_process(self):
        # load origin dfs from data file
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # process each complex
        pool, tasks = TaskPool('process', self.args.n_workers).start(), []
        for top_path, traj_path in self.tasks:
            top_path = Path(top_path).resolve()
            # check result exists
            result_path = top_path.parent / f'{top_path.stem}_porcupine.pse'
            if result_path.exists():
                if self.args.delete:
                    result_path.unlink()
                    put_log(f'{top_path.stem}_porcupine.pse deleted.')
                elif not self.args.force:
                    put_log(f'{top_path.stem}_porcupine.pse already exists, skip.')
                    continue
            # perform porcupine
            tasks.append(pool.add_task(None, run_porcupine, str(top_path), traj_path,
                                       self.args.begin_frame, self.args.end_frame, self.args.traj_step))
            pool.wait_till(lambda: pool.count_waiting_tasks() == 0, 0.01, update_result_queue=False)
        pool.wait_till_tasks_done(tasks)
        pool.close()


_str2func = {
    'trjconv': trjconv,
    'make_ndx': make_ndx,
    'simple': simple,
    'mmpbsa': mmpbsa,
    'interaction': interaction,
    'rrcs': RRCS,
    'porcupine': porcupine,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS analysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    # dev code
    # main('rrcs -d data_tmp/gmx/run1 -top md.tpr -traj md_center.xtc -c A Z -np 2 --backend cuda'.split(' '))
    
    main()