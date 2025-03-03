'''
Date: 2024-12-15 19:25:42
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-01-07 12:14:35
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


def align_pose_to_axis_warp(pml_name: str, move_name: str = None, fixed: Union[List[float], str] = 'center', state: int = 0, move_method: str = 'rotate', dss: int = 0, quiet: int = 0):
    from lazydock.pml.align_to_axis import align_pose_to_axis
    align_pose_to_axis(pml_name, move_name, fixed, state, move_method, dss, quiet)

cmd.extend('align_pose_to_axis', align_pose_to_axis_warp)


def open_vina_config_as_box(config_path: str, spacing: float = 1.0, linewidth: float = 2.0, r: float = 1.0, g: float = 1.0, b: float = 1.0):
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
    draw_box(minz, miny, minz, maxz, maxy, maxz, linewidth=linewidth, r=r, g=g, b=b)
    
cmd.extend('open_vina_config_as_box', open_vina_config_as_box)



def calcu_RRCS(model: str):
    from lazydock.pml.rrcs import calcu_RRCS
    df = calcu_RRCS(model)
    path = os.path.abspath(f'{model}_RRCS.xlsx')
    df.to_excel(path)
    print(f'RRCS saved to {path}')

cmd.extend('calcu_RRCS', calcu_RRCS)


def apply_shader_from_interaction_df(df_path: str, obj_name: str, cmap: str = 'coolwarm', alpha_mode: str = None,
                                     show_cbar: bool = False):
    from lazydock.pml.shader import Shader, ShaderValues
    values = ShaderValues().from_interaction_df(df_path, obj_name)
    shader = Shader(cmap)
    shader.create_colors_in_pml(values)
    shader.apply_shader_values(values, alpha_mode=alpha_mode)
    if show_cbar:
        shader.show_cbar(show=True)
    
cmd.extend('apply_shader_from_interaction_df', apply_shader_from_interaction_df)


def apply_shader_from_df(df_path: str, chain_col: str, resi_col: str, c_value_col: str,
                         obj_name: str, cmap: str = 'coolwarm', alpha_mode: str = None,
                         save_cbar: bool = False):
    from lazydock.pml.shader import Shader, ShaderValues
    values = ShaderValues().from_cols_df(df_path, obj_name, chain_col, resi_col, c_value_col)
    shader = Shader(cmap)
    shader.create_colors_in_pml(values)
    shader.apply_shader_values(values, alpha_mode=alpha_mode)
    if save_cbar:
        shader.plor_cbar(save=True)
    
cmd.extend('apply_shader_from_df', apply_shader_from_df)


print('LazyDock plugin loaded.')
print('''
Commands (python API):
    start_lazydock_server(host='localhost', port=8085, quiet=1)
    align_pose_to_axis(pml_name, move_name='', fixed='center', state=0, move_method='rotate', dss=1, quite=0)
    open_vina_config_as_box(config_path, spacing=1.0)
    calcu_RRCS(model: str)
    apply_shader_from_interaction_df(df_path: str, obj_name: str, cmap: str = 'coolwarm', alpha_mode: str ='cartoon_transparency')
    apply_shader_from_df(df_path: str, chain_col: str, resi_col: str, c_value_col: str, obj_name: str, cmap: str = 'coolwarm', alpha_mode: str ='cartoon_transparency', save_cbar: bool = False)
''')