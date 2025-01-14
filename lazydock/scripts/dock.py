'''
Date: 2024-12-04 20:58:39
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-01-14 19:30:23
Description: 
'''

import argparse
import os
import time
from functools import wraps
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from mbapy_lite.base import Configs, put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file, write_sheets
from mbapy_lite.plot import save_show
from mbapy_lite.web import TaskPool, random_sleep
from tqdm import tqdm

from lazydock.pml.autodock_utils import DlgFile
from lazydock.scripts._script_utils_ import Command, clean_path, excute_command
from lazydock.web.dinc import run_dock_on_DINC_ensemble
from lazydock.web.hdock import run_dock_on_HDOCK, run_dock_on_HPEPDOCK


class vina(Command):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.taskpool = None
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                                help='config file directory contains config files, each sub-directory is a task. Default is %(default)s.')
        args.add_argument('-c', '--config-name', type = str, default='config.txt',
                                help='config file name, each config is a task. Default is %(default)s.')
        args.add_argument('-v', '--vina-name', type = str, default='vina',
                                help='vina executable name to call. Default is %(default)s.')
        args.add_argument('--vina-args', type = str, default='--log ./log.txt',
                                help='args for vina executable. Default is %(default)s.')
        args.add_argument('-n', '--n-workers', type=int, default=1,
                                help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.taskpool = TaskPool('threads', self.args.n_workers).start()
        self.args.vina_args = '' if self.args.vina_args in {"''", '""', ''} else self.args.vina_args
        
    @staticmethod
    def run_vina(config_path: Path, vina_name: str, vina_args: str):
        print(f'current: {config_path}')
        config = opts_file(config_path, way='lines')
        out_name = list(filter(lambda x: x.startswith('out'), config))[0].split('=')[1].strip()
        if (config_path.parent / out_name).exists():
            print(f'{config_path.parent} has done, skip')
            return 
        cmd_string = f'cd "{config_path.parent}" && {vina_name} --config ./{config_path.name} {vina_args}'
        print(f'running: {cmd_string}')
        os.system(cmd_string)
        
    def main_process(self):
        if os.path.isdir(self.args.dir):
            configs_path = get_paths_with_extension(self.args.dir, ['.txt'], name_substr=self.args.config_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.dir}, exit.', _exit=True)
        print(f'get {len(configs_path)} config(s) for docking')
        tasks = []
        for config_path in tqdm(configs_path, total=len(configs_path)):
            tasks.append(self.taskpool.add_task(None, self.run_vina, Path(config_path), self.args.vina_name, self.args.vina_args))
            while self.taskpool.count_waiting_tasks() > 1:
                time.sleep(1)
        self.taskpool.wait_till_tasks_done(tasks)
        self.taskpool.close()
        
        
def hdock_run_fn_warpper(result_prefix: str = 'HDOCK', result_name: str = 'HDOCK_all_results.tar.gz'):
    def ret_warpper(func):
        def core_wrapper(*args, **kwargs):
            config_path = args[0] if len(args) > 0 else kwargs.get('config_path', None)
            if config_path is None:
                put_err('config_path is required, exit.', _exit=True)
            # get parameters from config file
            if isinstance(config_path, Path):
                root = config_path.parent
                parameters = hdock.get_paramthers_from_config(config_path)
                parameters['receptor_path'] = config_path.parent / parameters['receptor']
                parameters['ligand_path'] = config_path.parent / parameters['ligand']
            # OR get parameters from receptor and ligand path in tuple
            elif isinstance(config_path, tuple):
                root = Path(config_path[0]).parent
                parameters = {'receptor_path': config_path[0], 'ligand_path': config_path[1]}
            # un-expected type
            else:
                put_err(f'config_path type not support: {type(config_path)}, exit.', _exit=True)
            # check if done
            if (root / f'{result_prefix}_all_results.tar.gz').exists() or (root / result_name).exists():
                return print(f'{root} has done, skip')
            # perform docking
            print(f'current: {root}')
            ret =  func(*args, parameters=parameters, **kwargs)
            return ret
        return core_wrapper
    return ret_warpper


class hdock(Command):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                          help='vina config file directory contains config files (named "config.txt"), each sub-directory is a task. Default is %(default)s.')
        args.add_argument('-r', '--receptor', type = str, default=None,
                          help="receptor pdb file name, optional. If provided, will ignore config.txt.")
        args.add_argument('-l', '--ligand', type = str, default=None,
                          help="ligand pdb file name, optional. If provided, will ignore config.txt.")
        args.add_argument('--config-name', type=str, default=None,
                          help='Vina config for specifing parameters. Default is %(default)s.')
        args.add_argument('-m', '--method', type = str, default='web', choices=['web'],
                          help='docking method. Currently support "web". Default is %(default)s.')
        args.add_argument('--email', type = str, default=None,
                          help='email address for HDOCK web server. Default is %(default)s.')
        args.add_argument('-gui', '--gui', action='store_true', default=False,
                          help='show browser GUI. Default is %(default)s.')
        return args
    
    @staticmethod
    def get_paramthers_from_config(config_path: Path) -> Dict:
        config_lines = config_path.read_text().split('\n')
        paramthers = {line.split('=')[0].strip():line.split('=')[1].strip() for line in config_lines if not line.startswith('#')}
        return {k:v[:v.find('#')] for k,v in paramthers.items()}

    @staticmethod
    @hdock_run_fn_warpper('HDOCK')
    def run_hdock_web(config_path: Union[Path, Tuple[str, str]], parameters: Dict[str, str] = None, email=None, browser=None,
                      args: argparse.ArgumentParser = None):
        w_dir = config_path.parent if isinstance(config_path, Path) else Path(config_path[0]).parent
        run_dock_on_HDOCK(receptor_path=parameters['receptor_path'], ligand_path=parameters['ligand_path'],
                          w_dir=w_dir, email=email, browser=browser)

    @staticmethod
    @hdock_run_fn_warpper('HDOCK')
    def run_hdock_local(config_path: Union[Path, Tuple[str, str]], parameters: Dict[str, str] = None, **kwargs):
        raise NotImplementedError('local docking not implemented yet.')

    def main_process(self):
        if not os.path.isdir(self.args.dir):
            put_err(f'dir argument should be a directory: {self.args.config}, exit.', _exit=True)
        if self.args.receptor is not None and self.args.ligand is not None:
            r_paths = get_paths_with_extension(self.args.dir, [], name_substr=self.args.receptor)
            l_paths = get_paths_with_extension(self.args.dir, [], name_substr=self.args.ligand)
            if len(r_paths) != len(l_paths):
                r_roots = [os.path.dirname(p) for p in r_paths]
                l_roots = [os.path.dirname(p) for p in l_paths]
                roots_count = {root: r_roots.count(root)+l_roots.count(root) for root in (set(r_roots) | set(l_roots))}
                invalid_roots = '\n'.join([root for root, count in roots_count.items() if count != 2])
                return put_err(f"The number of receptor and ligand files is not equal, please check the input files.\ninvalid roots:\n{invalid_roots}")
            configs_path = [(r, l) for r, l in zip(r_paths, l_paths)]
        elif self.args.config_name is not None:
            configs_path = get_paths_with_extension(self.args.dir, ['.txt'], name_substr=self.args.config_name)
        else:
            return put_err('config_name or receptor and ligand should be provided, exit.', _exit=True)
        print(f'get {len(configs_path)} config(s) for docking')
        # allow browser gui
        if self.args.gui:
            from mbapy_lite.web import Browser
            browser = Browser(options=[f"--user-agent={Configs.web.chrome_driver_path}"])
        else:
            browser = None
        # perform docking
        dock_fn = getattr(self, f'run_hdock_{self.args.method}')
        for config_path in tqdm(configs_path, total=len(configs_path)):
            dock_fn(Path(config_path).resolve() if isinstance(config_path, str) else config_path,
                    email=self.args.email, browser=browser, args=self.args)
            random_sleep(300, 180) # sleep 3~5 minutes to avoid overloading the server


class hpepdock(hdock):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = hdock.make_args(args)
        # make sure hpepdock only support web method
        method_arg_idx = args._actions.index(list(filter(lambda x: x.dest =='method', args._actions))[0])
        args._actions[method_arg_idx].choices = ['web']
        args._actions[method_arg_idx].default = 'web'
        args._actions[method_arg_idx].help = 'docking method. Currently support "web". Default is %(default)s.'
        return args
    
    @staticmethod
    @hdock_run_fn_warpper('HPEPDOCK')
    def run_hdock_web(config_path: Path, parameters: Dict[str, str] = None, email=None, browser=None,
                      args: argparse.ArgumentParser = None):
        w_dir = config_path.parent if isinstance(config_path, Path) else Path(config_path[0]).resolve().parent
        run_dock_on_HPEPDOCK(receptor_path=parameters['receptor_path'], ligand_path=parameters['ligand_path'],
                             w_dir=w_dir, email=email, browser=browser)


class dinc_ensemble(hdock):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = hdock.make_args(args)
        args.add_argument('--use-config-box', action='store_true', default=False,
                          help='FLAG, whether to use grid box from config file. Default is %(default)s.')
        # make sure only support web method
        method_arg_idx = args._actions.index(list(filter(lambda x: x.dest =='method', args._actions))[0])
        args._actions[method_arg_idx].choices = ['web']
        args._actions[method_arg_idx].default = 'web'
        args._actions[method_arg_idx].help = 'docking method. Currently support "web". Default is %(default)s.'
        # make sure eamil is required
        email_arg_idx = args._actions.index(list(filter(lambda x: x.dest == 'email', args._actions))[0])
        args._actions[email_arg_idx].required = True
        args._actions[email_arg_idx].help = 'email address for HDOCK web server. Required.'
        return args
    
    def process_args(self):
        if self.args.receptor is not None and self.args.ligand is not None and self.args.config_name is not None:
            return put_err('receptor, ligand and config_name cannot be provided at the same time, exit.', _exit=True)
        return super().process_args()
    
    @staticmethod
    @hdock_run_fn_warpper('', 'DINC-Ensemble_result.zip')
    def run_hdock_web(config_path: Path, parameters: Dict[str, str] = None, email=None, browser=None,
                      args: argparse.ArgumentParser = None):
        w_dir = config_path.parent if isinstance(config_path, Path) else Path(config_path[0]).resolve().parent
        put_log(f'working in {w_dir}')
        if args.use_config_box:
            box_center = {k:parameters[f'center_{k}'] for k in ['x', 'y', 'z']}
            box_size = {k:parameters[f'size_{k}'] for k in ['x', 'y', 'z']}
        else:
            box_center, box_size = 'receptor', 'ligand'
        put_log(f'parameters: {parameters}, box_center: {box_center}, box_size: {box_size}')
        run_dock_on_DINC_ensemble(receptor_path=parameters['receptor_path'], ligand_path=parameters['ligand_path'],
                                  email=email, box_center=box_center, box_size=box_size, w_dir=str(w_dir), browser=browser)


def convert_result_run_convert(input_path: Path, output_path: Path, method: str):
    if method in {'lazydock', 'obabel'}:
        getattr(convert_result, f'run_convert_{method}')(input_path, output_path)
    else:
        put_err(f'unsupported convert method: {method}, exit.', _exit=True)


class convert_result(Command):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.taskpool = None

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                                help='input directory. Default is %(default)s.')
        args.add_argument('-n', '--name', type = str, default='',
                          help='input file name. Default is %(default)s.')
        args.add_argument('-i', '--input-type', type = str, default='pdbqt,dlg',
                          help='input file type. Default is %(default)s.')
        args.add_argument('-o', '--output-type', type = str, default='pdb', choices=['pdb'],
                          help='input file type. Default is %(default)s.')
        args.add_argument('-s', '--suffix', type = str, default='',
                          help='output file suffix. Default is %(default)s.')
        args.add_argument('-m', '--method', type = str, default='lazydock', choices=['lazydock', 'obabel'],
                          help='convert tools to use. Currently support "lazydock, obabel". Default is %(default)s.')
        args.add_argument('--n-workers', type=int, default=4,
                          help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        self.args.input_type = self.args.input_type.split(',')
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.taskpool = TaskPool('process', self.args.n_workers).start()
        
    @staticmethod
    def run_convert_lazydock(input_path: Path, output_path: Path):
        dlg = DlgFile(input_path)
        pdb_string = '\n'.join([f'MODEL {i+1}\n'+pose.as_pdb_string().replace('\r\n', '\n')+'ENDMDL' for i, pose in enumerate(dlg.pose_lst)])
        opts_file(output_path, 'w', data=pdb_string)

    @staticmethod
    def run_convert_obabel(input_path: Path, output_path: Path):
        ty1, ty2 = input_path.suffix.lower()[1:], output_path.suffix.lower()[1:]
        os.system(f'obabel -i{ty1} "{str(input_path)}" -o{ty2} -O "{str(output_path)}"')

    def main_process(self):
        input_paths = get_paths_with_extension(self.args.dir, self.args.input_type, name_substr=self.args.name)
        print(f'get {len(input_paths)} input(s) for convert:\n', '\n'.join([f'{i+1}. {x}' for i, x in enumerate(input_paths)]))
        if input('start convert? (y/n) ').lower() != 'y':
            return 
        tasks = []
        for input_path in tqdm(input_paths, total=len(input_paths)):
            input_path = Path(input_path)
            output_path = input_path.parent / f'{input_path.stem}{self.args.suffix}.{self.args.output_type}'
            tasks.append(self.taskpool.add_task(None, convert_result_run_convert, input_path, output_path, self.args.method))
            while self.taskpool.count_waiting_tasks() > 0:
                time.sleep(1)
        self.taskpool.wait_till_tasks_done(tasks)
        self.taskpool.close()


class cluster_result(Command):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.taskpool = None

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                                help='input directory. Default is %(default)s.')
        args.add_argument('-n', '--name', type = str, default='',
                          help='input file name, such as "dock.pdbqt". Default is %(default)s.')
        args.add_argument('-m', '--method', type = str, default='pose', choices=['pose', 'interaction'],
                          help='cluster method to use. Currently support "pose, interaction". Default is %(default)s.')
        args.add_argument('--range', type = str, default='2,7',
                          help='range of cluster groups, input as "min,max". Default is %(default)s.')
        args.add_argument('--repeat', type=int, default=3,
                          help='number of repeats for clustering. Default is %(default)s.')
        args.add_argument('--n-workers', type=int, default=4,
                          help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        self.args.range = [int(x) for x in self.args.range.split(',')]
        if len(self.args.range) != 2 or self.args.range[0] >= self.args.range[1]:
            put_err(f'range should be "min,max", and min < max, got {self.args.range}, exit.', _exit=True)
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.taskpool = TaskPool('process', self.args.n_workers).start()
        
    @staticmethod
    def cluster_on_pose(input_path: Path, _range: Tuple[int, int], repeat: int):
        wdir = input_path.parent / 'cluster_pose'
        os.makedirs(str(wdir), exist_ok=True)
        dlg = DlgFile(input_path)
        dlg.sort_pose(inplace=True)
        rmsd = dlg.rmsd('numpy')
        sns.clustermap(rmsd)
        save_show(str(wdir / 'RMSD_cluster.png'), 600, show=False)
        rs_df = pd.DataFrame(columns=['repeat', 'k', 'sse','ssr', 'r'])
        dfs = {'rs_df': rs_df}
        for i in range(repeat):
            for k in range(_range[0], _range[1]+1):
                groups_idx = dlg.rmsd_cluster(k) # [N, ]
                groups_idx_lst = list(groups_idx)
                group_sizes = {idx: groups_idx_lst.count(idx) for idx in range(k)}
                sse, ssr, r = dlg.calcu_SSE_SSR(rmsd, groups_idx)
                rs_df.loc[len(rs_df)] = [i, k, sse, ssr, r]
                dfs[f'repeat_{i}_k_{k}'] = pd.DataFrame(data={'pose_idx': range(len(dlg.pose_lst)),
                                                              'groups_idx': groups_idx,
                                                              'group_size': [group_sizes[groups_idx[idx]] for idx in range(len(dlg))],
                                                              'energy': list(map(lambda x: x.energy, dlg.pose_lst))})
        fig = plt.figure()
        sns.lineplot(data=rs_df, x='k', y='r', hue='repeat', ax=fig.gca())
        save_show(str(wdir / 'SSR-SSE_lines.png'), 600, show=False)
        plt.close(fig)
        write_sheets(str(wdir / 'cluster.xlsx'), dfs)
        put_log(f'cluster result saved to {wdir}.')

    @staticmethod
    def cluster_on_interaction(input_path: Path):
        raise NotImplementedError('interaction clustering not implemented yet.')

    def main_process(self):
        input_paths = get_paths_with_extension(self.args.dir, [], name_substr=self.args.name)
        tasks = []
        for input_path in tqdm(input_paths, total=len(input_paths)):
            input_path = Path(input_path)
            tasks.append(self.taskpool.add_task(None, getattr(self, f'cluster_on_{self.args.method}'),
                                                input_path, self.args.range, self.args.repeat))
            while self.taskpool.count_waiting_tasks() > 0:
                time.sleep(1)
        self.taskpool.wait_till_tasks_done(tasks)
        self.taskpool.close()


_str2func = {
    'vina': vina,
    'hdock': hdock,
    'hpepdock': hpepdock,
    'dinc-ensemble': dinc_ensemble,
    'convert-result': convert_result,
    'cluster-result': cluster_result,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'perform docking with Vina or other docking software.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    vina_args = vina.make_args(subparsers.add_parser('vina', description='perform vina molecular docking.'))
    hdock_args = hdock.make_args(subparsers.add_parser('hdock', description='perform HDOCK molecular docking.'))
    hpepdock_args = hpepdock.make_args(subparsers.add_parser('hpepdock', description='perform HPEPDOCK molecular docking.'))
    dinc_ensemble_args = dinc_ensemble.make_args(subparsers.add_parser('dinc-ensemble', description='perform DINC-Ensemble molecular docking.'))
    convert_result_args = convert_result.make_args(subparsers.add_parser('convert-result', description='convert docking result file to pdb format.'))
    cluster_result_args = cluster_result.make_args(subparsers.add_parser('cluster-result', description='cluster docking result.'))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    # main('convert-result -d data_tmp/docking/ligand1 -m lazydock --n-workers 1'.split())
    # main('dinc-ensemble -d data_tmp/docking/ligand1 -m web --email 2262029386@qq.com --config-name config.txt --use-config-box'.split())
    # main('cluster-result -d data_tmp/docking/ligand1 -n dock.pdbqt'.split())
    
    main()