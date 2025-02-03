import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from pymol import cmd
from tqdm import tqdm

from lazydock.gmx.run import Gromacs
from lazydock.scripts._script_utils_ import Command, clean_path
from lazydock.utils import uuid4


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
        gmx.run_command_with_expect('trjconv', s=f'{main_name}.tpr', f=f'{main_name}.xtc', o=f'{main_name}_center.xtc', pbc='mol', center=True,
                                    expect_actions=[{'Select a group:': f'{center_group}\r', '\\timeout': f'{center_group}\r'},
                                                    {'Select a group:': '0\r', '\\timeout': '0\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        
    @staticmethod
    def rms(gmx: Gromacs, main_name: str, group: str = '4', **kwargs):
        gmx.run_command_with_expect('rms', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'rmsd.xvg', tu='ns',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        os.system(f'cd "{gmx.working_dir}" && dit xvg_compare -c 1 -f rmsd.xvg -o rmsd.png -smv -t "RMSD of {main_name}" -csv {main_name}_rmsd.csv -ns')
        
    @staticmethod
    def rmsf(gmx: Gromacs, main_name: str, group: str = '4', res: bool = True, **kwargs):
        gmx.run_command_with_expect('rmsf', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'rmsf.xvg', res=res,
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        os.system(f'cd "{gmx.working_dir}" && dit xvg_compare -c 1 -f rmsf.xvg -o rmsf.png -t "RMSF of {main_name}" -csv {main_name}_rmsf.csv -ns')
        
    @staticmethod
    def gyrate(gmx: Gromacs, main_name: str, group: str = '4', **kwargs):
        gmx.run_command_with_expect('gyrate', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'gyrate.xvg',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        os.system(f'cd "{gmx.working_dir}" && dit xvg_compare -c 1 -f gyrate.xvg -o gyrate.png -smv -t "Gyrate of {main_name}" -csv {main_name}_gyrate.csv -ns')
        
    @staticmethod
    def hbond(gmx: Gromacs, main_name: str, group: str = '1', dt=10, **kwargs):
        gmx.run_command_with_expect('hbond', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc',
                                    num=f'{main_name}_hbond_num.xvg', dist=f'{main_name}_hbond_dist.xvg',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        os.system(f'cd "{gmx.working_dir}" && dit xvg_compare -c 1 -f {main_name}_hbond_num.xvg -o hbond_num.png -smv -t "H-bond num of {main_name}" -csv {main_name}_hbond_num.csv -ns')
        os.system(f'cd "{gmx.working_dir}" && dit xvg_show -f {main_name}_hbond_dist.xvg -o hbond_dist.png -ns')

    @staticmethod
    def sasa(gmx: Gromacs, main_name: str, group: str = '4', **kwargs):
        gmx.run_command_with_expect('sasa -or sasa_res.xvg', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc',
                                    o=f'sasa_total.xvg', odg=f'sasa_dg.xvg', tv='sasa_tv.xvg', tu='ns',
                                    expect_actions=[{'>': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        for ty in ['total', 'res', 'dg', 'tv']:
            os.system(f'cd "{gmx.working_dir}" && dit xvg_compare -c 1 -f sasa_{ty}.xvg -o sasa_{ty}.png -smv -t "SASA {ty} of {main_name}" -csv {main_name}_sasa_{ty}.csv -ns')

    @staticmethod
    def covar(gmx: Gromacs, main_name: str, group: str = '4', xmax: int = 15, **kwargs):
        gmx.run_command_with_expect('covar', s=f'{main_name}.tpr', f=f'{main_name}_center.xtc', o=f'eigenval.xvg', tu='ns',
                                    expect_actions=[{'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'},
                                                    {'Select a group:': f'{group}\r', '\\timeout': f'{group}\r'}],
                                    expect_settings={'timeout': 10}, **kwargs)
        os.system(f'cd "{gmx.working_dir}" && dit xvg_compare -c 1 -f eigenval.xvg -o eigenval.png -xmin 0 -xmax {xmax} -smv -t "Eigenval of {main_name}" -csv {main_name}_eigenval.csv -ns')
    
    @staticmethod
    def free_energy_landscape(gmx: Gromacs, main_name: str, **kwargs):
        os.system(f'cd "{gmx.working_dir}" && md-davis landscape_xvg -c -T 300 -x rmsd.xvg -y gyrate.xvg -o FEL.html -n FEL -l "RMSD-Rg" --axis_labels "dict(x=\'RMSD (in nm)\', y=\'Rg (in nm)\', z=\'Free Energy (kJ mol<sup>-1</sup>)<br>\')"')
    
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


_str2func = {
    'simple': simple,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS analysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    simple_args = simple.make_args(subparsers.add_parser('simple', description=simple.HELP))

    args = args_paser.parse_args(sys_args)
    if args.sub_command in _str2func:
        _str2func[args.sub_command](args).excute()


if __name__ == '__main__':
    main()