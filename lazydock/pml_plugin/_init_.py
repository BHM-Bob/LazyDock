'''
Date: 2024-12-15 19:25:42
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-22 20:46:09
Description: 
'''
import os
from typing import List, Union

from pymol import cmd


def start_lazydock_server(host: str = 'localhost', port: int = 8085, quiet: int = 1):
    from lazydock.pml.server import VServer
    print(f'Starting LazyDock server on {host}:{port}, quiet={quiet}')
    VServer(host, port, not bool(quiet))

cmd.extend('start_lazydock_server', start_lazydock_server)


def align_pose_to_axis_warp(pml_name: str, move_name: str = None, fixed: Union[List[float], str] = 'center', state: int = 0, move_method: str = 'rotate', dss: int = 1, quiet: int = 0):
    from lazydock.pml.align_to_axis import align_pose_to_axis
    print('try drag coords or matrix manually before running this command will get secondary structure remained.')
    align_pose_to_axis(pml_name, move_name, fixed, state, move_method, dss, quiet)

cmd.extend('align_pose_to_axis', align_pose_to_axis_warp)


def open_vina_config_as_box(config_path: str, spacing: float = 1.0):
    from mbapy_lite.file import opts_file

    from lazydock.pml.thirdparty.draw_bounding_box import draw_box
    if not os.path.exists(config_path):
        return print(f'Config file {config_path} not found, skip.')
    cfg = opts_file(config_path, way='lines')
    get_line = lambda n: [line for line in cfg if line.startswith(n)][0]
    center = {line.split('=')[0].strip(): float(line.split('=')[1].strip()) for line in map(get_line, ['center_x', 'center_y', 'center_z'])}
    size = {line.split('=')[0].strip(): float(line.split('=')[1].strip()) for line in map(get_line, ['size_x', 'size_y', 'size_z'])}
    print(f'center: {center}, size: {size}')
    minz, miny, minz = [(center[f'center_{k}'] - size[f'size_{k}'] / 2) * spacing for k in ['x', 'y', 'z']]
    maxz, maxy, maxz = [(center[f'center_{k}'] + size[f'size_{k}'] / 2) * spacing for k in ['x', 'y', 'z']]
    draw_box(minz, miny, minz, maxz, maxy, maxz)
    
cmd.extend('open_vina_config_as_box', open_vina_config_as_box)


print('LazyDock plugin loaded.')
print('''
Commands (python API):
    start_lazydock_server(host='localhost', port=8085, quiet=1)
    align_pose_to_axis(pml_name, move_name='', fixed='center', state=0, move_method='rotate', dss=1, quite=0)
    open_vina_config_as_box(config_path, spacing=1.0)
''')