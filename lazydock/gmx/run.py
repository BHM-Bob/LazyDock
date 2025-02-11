'''
Date: 2024-12-18 10:48:32
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-11 15:07:12
Description:
'''
import os
from typing import Any, Dict, List

from mbapy_lite.base import get_fmt_time, put_log
from mbapy_lite.file import opts_file
from mbapy_lite.game import BaseInfo

if __name__ == '__main__':
    from lazydock.utils import get_storage_path
else:
    from ..utils import get_storage_path


class Gromacs(BaseInfo):
    def __init__(self, call_name: str = 'gmx', working_dir: str = '.'):
        super().__init__()
        self.call_name = call_name
        self.working_dir = os.path.abspath(working_dir)
        
    def kwargs2cmd(self, kwargs: Dict[str, str]):
        return ' '.join([f'{"--" if k.startswith("_") else "-"}{k} {v}' for k, v in kwargs.items()])
    
    def gen_command(self, sub_commmand: str, **kwargs):
        """
        Run gromacs command.
        
        Parameters:
            - sub_commmand: str, sub command of gromacs, such as grompp, mdrun, etc.
            - **kwargs: dict, keyword arguments for gromacs command.
                - if the key starts with '_', it will be added with '--' before the key, else with '-' before the key.
                - if the value is bool, it will be added with no value.
                - if the value is list or tuple, it will be joined with space.
                - if the value is not str, it will be converted to str directly via str().
        """
        for k in list(kwargs.keys()):
            if kwargs[k] is None:
                del kwargs[k]
            elif isinstance(kwargs[k], bool):
                if kwargs[k]:
                    kwargs[k] = ''
                else:
                    kwargs[f'no{k}'] = ''
            elif isinstance(kwargs[k], (list, tuple)):
                kwargs[k] =' '.join(map(str, kwargs[k]))
            elif not isinstance(kwargs[k], str):
                kwargs[k] = str(kwargs[k])
        return f'cd "{self.working_dir}" && {self.call_name} {sub_commmand} {self.kwargs2cmd(kwargs)}'
    
    def run_command_with_expect(self, cmd: str, expect_actions: List[Dict[str, str]] = None,
                                expect_settings: Dict[str, Any] = None, enable_log: bool = False):
        """
        Run gromacs command with expect script.
        
        Parameters: 
            - sub_commmand: str, sub command of gromacs, such as grompp, mdrun, etc.
            - expect_actions: List[Dict[str, str]], the expect script actions.
                - key: the string to match the output of the command. if key is '\\timeout', it will be treated as a timeout.
                - value: the value to send
            - expect_settings: dict, the expect script settings.
                - timeout: int, default is -1, the timeout for expect script to start.
            - enable_log: bool, default is False, whether to enable log file.
            - **kwargs: dict, keyword arguments for gromacs command, generate by gen_command() method.
        """
        # set name and paths
        scripts_dir = os.path.join(self.working_dir, 'LazyDock_gmx_scripts')
        os.makedirs(scripts_dir, exist_ok=True)
        scripts_name = f'{get_fmt_time("%Y-%m-%d-%H-%M-%S.%f")}'
        # generate command
        if enable_log:
            log_path = os.path.join(scripts_dir, f'{scripts_name}.log')
            cmd = f'{cmd} > ./LazyDock_gmx_scripts/{scripts_name}.log'.replace('&&', '\n')
        put_log(f'Get command:\n{cmd}', head='LazyDock')
        # just run the command if no expect actions
        if expect_actions is None or not expect_actions:
            ret_val = os.system(cmd)
            if enable_log:
                ret_val = (ret_val, log_path)
            return ret_val
        # save cmd to bash file
        bash_path = os.path.join(scripts_dir, f'{scripts_name}.sh')
        opts_file(bash_path, 'w', data=cmd)
        # create expect script
        expect_settings = expect_settings or {}
        expect_lines = []
        expect_lines.append(f'set timeout {expect_settings.get("timeout", -1)}')
        expect_lines.append(f'spawn bash "{bash_path}"')
        for action in expect_actions:
            expect_lines.append('expect {')
            for key, value in action.items():
                if key == '\\timeout':
                    expect_lines.append(f'    timeout {{\n        puts "===TIMEOUT==="\n        send "{value}"}}')
                else:
                    expect_lines.append(f'    "{key}" {{\n        send "{value}"}}')
            expect_lines.append('}\n')
        expect_lines.append('interact')
        expect_script = '\n'.join(expect_lines)
        # save expect script to file and run it
        script_path = os.path.join(scripts_dir, f'{scripts_name}.exp')
        opts_file(script_path, 'w', data=expect_script)
        put_log(f'Running expect script: {script_path}', head='LazyDock')
        ret_val = os.system(f'cd "{self.working_dir}" && expect "{script_path}"')
        if enable_log:
            ret_val = (ret_val, log_path)
        return ret_val
    
    def run_gmx_with_expect(self, sub_commmand: str, expect_actions: List[Dict[str, str]] = None,
                            expect_settings: Dict[str, Any] = None, enable_log: bool = False, **kwargs):
        """
        Run gromacs command with expect script.
        """
        cmd = self.gen_command(sub_commmand, **kwargs)
        self.run_command_with_expect(cmd, expect_actions, expect_settings, enable_log)
        
    def run_normal_cmd_with_expect(self, cmd: str, expect_actions: List[Dict[str, str]] = None,
                                  expect_settings: Dict[str, Any] = None, enable_log: bool = False):
        """
        Run normal command with expect script.
        """
        self.run_command_with_expect(cmd, expect_actions, expect_settings, enable_log)
        
        
        
if __name__ == '__main__':
    gmx = Gromacs()
    gmx.run_gmx_with_expect('grompp', f='topol.top', c='conf.gro', p='topol.top', o='tpr', maxwarn=1)
    gmx.run_gmx_with_expect('pdb2gmx -f receptor.pdb -o processed.gro -ter -ignh', [{')': '1\r'}, {'None': '1\r'}, {'None': '1\r'}, {'None': '1\r'}])