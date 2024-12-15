'''
Date: 2024-12-15 19:25:42
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-15 19:35:39
Description: 
'''

from typing import List, Union

from pymol import cmd


def start_lazydock_server(host: str = 'localhost', port: int = 8085, quiet: int = 1):
    from lazydock.pml.server import VServer
    print(f'Starting LazyDock server on {host}:{port}, quiet={quiet}')
    VServer(host, port, not bool(quiet))

cmd.extend('start_lazydock_server', start_lazydock_server)


def align_pose_to_axis_warp(pml_name: str, fixed: Union[List[float], str] = 'center', state: int = 0):
    from lazydock.pml.align_to_axis import align_pose_to_axis
    align_pose_to_axis(pml_name, fixed, state, warp=True)

cmd.extend('align_pose_to_axis', align_pose_to_axis_warp)