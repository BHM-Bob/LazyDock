'''
Date: 2024-12-04 20:58:39
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-04 21:44:00
Description: 
'''

import argparse
import os
import time
from pathlib import Path
from typing import Dict, List

from mbapy_lite.base import put_err
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web import TaskPool
from tqdm import tqdm

if __name__ == '__main__':
    from lazydock.scripts._script_utils_ import (Command, clean_path,
                                                 excute_command)
else:
    from ._script_utils_ import Command, clean_path, excute_command
    

class vina(Command):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.taskpool = None
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-c', '--config', type = str, default='.',
                                help='config file path or directory contains config files (named "config.txt"), each sub-directory is a task. Default is %(default)s.')
        args.add_argument('-v', '--vina-name', type = str, default='vina',
                                help='vina executable name to call. Default is %(default)s.')
        args.add_argument('-n', '--n-workers', type=int, default=1,
                                help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.config = clean_path(self.args.config)
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.taskpool = TaskPool('threads', self.args.n_workers).start()
        
    @staticmethod
    def run_vina(config_path: Path, vina_name: str):
        print(f'current: {config_path}')
        if (config_path.parent / 'dock.pdbqt').exists():
            print(f'{config_path.parent} has done, skip')
            return 
        os.system(f'cd "{config_path.absolute().parent}" && {vina_name} --config ./config.txt --log ./log.txt')
        
    def main_process(self):
        if os.path.isfile(self.args.config):
            configs_path = [self.args.config]
        elif os.path.isdir(self.args.config):
            configs_path = get_paths_with_extension(self.args.config, ['.txt'], name_substr='config')
        else:
            put_err(f'config file or directory not found: {self.args.config}, exit.', _exit=True)
        print(f'get {len(configs_path)} config(s) for docking')
        tasks = []
        for config_path in tqdm(configs_path, total=len(configs_path)):
            tasks.append(self.taskpool.add_task(None, self.run_vina, Path(config_path)))
            while self.taskpool.count_waiting_tasks() > 1:
                time.sleep(1)
        self.taskpool.wait_till_tasks_done(tasks)
        self.taskpool.close()


_str2func = {
    'vina': vina,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'perform docking with Vina or other docking software.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    vina_args = vina.make_args(subparsers.add_parser('vina', description='perform vina molecular docking.'))
    
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main(['vina'])