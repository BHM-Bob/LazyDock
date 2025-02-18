'''
Date: 2024-11-23 19:53:42
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-18 15:39:50
Description:
'''
import argparse
import os
import traceback
from pathlib import Path
from typing import Any, Dict, List, Union

from mbapy_lite.base import (check_parameters_path, parameter_checker, put_err,
                             put_log)
from mbapy_lite.file import opts_file


def clean_path(path: str):
    return Path(path.replace('"', '').replace("'", '')).resolve()

def _print(content: str, f, verbose = True):
    if f is not None:
        f.write(content+'\n')
    if verbose:
        print(content)

def show_args(args, args_name: List[str], printf = print):
    printf('')
    for arg_name in args_name:
        printf(f'get arg: {arg_name}: {getattr(args, arg_name)}')
    printf('')

class Command:
    def __init__(self, args: argparse.Namespace, printf = print,
                 iter_run_arg: List[str] = None) -> None:
        self.args = args
        self.printf = printf
        self._pickle_except_list = []
        self.iter_run_arg = iter_run_arg or []
        
    def process_args(self):
        pass
    
    def main_process(self):
        pass
        
    def excute(self):
        self.process_args()
        show_args(self.args, list(self.args.__dict__.keys()), self.printf)
        if self.iter_run_arg:
            # [arg1_values: List, arg2_values: List, ...]
            iter_args = [getattr(self.args, n).copy() for n in self.iter_run_arg]
            sum_run = min(len(i) for i in iter_args)
            for i, args in enumerate(zip(*iter_args)):
                cur_arg_info = []
                # args: [arg1_value1, arg2_value2, ...]
                for n, v in zip(self.iter_run_arg, args):
                    setattr(self.args, n, v)
                    cur_arg_info.append(f'{n}={v}')
                put_log(f'running iter[{i+1}/{sum_run}] for args: {", ".join(cur_arg_info)}')
                self.main_process()
        else:
            return self.main_process()
    
    def save_session(self, module_name: str, module_path: str = 'mbapy.scripts', path: str = os.curdir):
        if not Path(path).parent.exists():
            os.makedirs(Path(path).parent, exist_ok=True)
        session = {'__module_path__': module_path, '__module_name__': module_name, '__cmd_name__': self.__class__.__name__}
        for k,v in self.args.__dict__.items():
            if k not in self._pickle_except_list:
                session[k] = v
        opts_file(path, 'wb', way = 'pkl', data = session)
        
    @parameter_checker(path = check_parameters_path)
    def load_session(self, path: str):
        return opts_file(path, 'rb', way = 'pkl')
        
    def exec_from_session(self, session: Union[str, Dict[str, Any]]):
        if isinstance(session, str) and check_parameters_path(session):
            session = self.load_session(session)
        
    
def excute_command(args_paser: argparse.ArgumentParser, sys_args: List[str],
                   _str2func: Dict[str, callable]):
    args = args_paser.parse_args(sys_args)
    
    if args.sub_command in _str2func:
        cmd = _str2func[args.sub_command]
        try:
            # if is inherite from Command, excute it
            if isinstance(cmd, type) and issubclass(cmd, Command):
                cmd(args).excute()
            # if is a function, excute it
            elif callable(cmd):
                cmd(args)
            # else, just try
            else:
                cmd(args)
        except Exception as e:
            traceback.print_exception(type(e), e, e.__traceback__)
    else:
        put_err(f'no such sub commmand: {args.sub_command}')
    
        
        
