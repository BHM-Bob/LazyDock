'''
Date: 2024-12-15 19:25:42
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-07-05 19:06:27
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
    from lazydock.pml.thirdparty.draw_bounding_box import draw_box
    from mbapy_lite.file import opts_file
    if not os.path.exists(config_path):
        return print(f'Config file {config_path} not found, skip.')
    spacing = float(spacing)
    cfg = opts_file(config_path, way='lines')
    get_line = lambda n: [line for line in cfg if line.startswith(n)][0]
    center = {line.split('=')[0].strip(): float(line.split('=')[1].strip()) for line in map(get_line, ['center_x', 'center_y', 'center_z'])}
    size = {line.split('=')[0].strip(): float(line.split('=')[1].strip()) for line in map(get_line, ['size_x', 'size_y', 'size_z'])}
    print(f'center: {center}, size: {size}')
    minx, miny, minz = [(center[f'center_{k}'] - size[f'size_{k}'] / 2) * spacing for k in ['x', 'y', 'z']]
    maxx, maxy, maxz = [(center[f'center_{k}'] + size[f'size_{k}'] / 2) * spacing for k in ['x', 'y', 'z']]
    draw_box(minx, miny, minz, maxx, maxy, maxz, linewidth=linewidth, r=r, g=g, b=b)
    
cmd.extend('open_vina_config_as_box', open_vina_config_as_box)


def align_by_sele(ref_sele: str, move_sele: str, move_all: str = None, n_cycle: int = 1):
    """Align move_sele to ref_sele, move_all is None, then move move_sele only, else move move_all by align matrix.
    1. extract coords from ref_sele and move_sele.
    2. calculate align matrix.
    3. move move_sele by align matrix.
    4. if move_all is not None, move move_all by align matrix.
    
    Parameters:
        ref_sele: Reference selection
        move_sele: Selection to be aligned
        move_all: Optional selection to apply the same transformation (defaults to move_sele)
        n_cycle: Number of alignment iterations to perform (default: 1)
    """
    move_all = move_all or move_sele
    # Import required APIs inside the function as requested
    from lazydock.algorithm.rms import calc_rms_rotational_matrix
    from lazydock.pml.utils import rotation_matrix_to_pymol_matrix
    import numpy as np
    
    # Extract coordinates from selections using PyMOL API
    ref_coords = np.array(cmd.get_coords(ref_sele))
    move_coords = np.array(cmd.get_coords(move_sele))
    
    if ref_coords.shape[0] != move_coords.shape[0]:
        print(f"Error: Number of atoms don't match - ref_sele has {ref_coords.shape[0]} atoms, move_sele has {move_coords.shape[0]} atoms")
        return
    
    # Store the first RMSD value
    first_rmsd = None
    
    # Perform alignment iterations
    for cycle in range(n_cycle):
        # Calculate rotation matrix using RMS API
        rot_matrix_flat = np.zeros(9)
        rmsd, rot_matrix = calc_rms_rotational_matrix(ref_coords, move_coords, rot_matrix_flat)
        print(f'Cycle {cycle}: RMSD = {rmsd:.4f}, Rotation Matrix:\n{rot_matrix}')
        # Store the first RMSD value
        if first_rmsd is None:
            first_rmsd = rmsd
        
        # Reshape to 3x3 matrix
        rot_matrix = rot_matrix.reshape(3, 3)
        
        # Calculate centers for translation
        ref_center = np.mean(ref_coords, axis=0)
        move_center = np.mean(move_coords, axis=0)
        
        # Convert to PyMOL matrix format using the utility function
        pymol_matrix = rotation_matrix_to_pymol_matrix(rot_matrix, ref_center, move_center)
        
        # Apply transformation to move_sele
        cmd.transform_selection(move_all, pymol_matrix, homogenous=0)
        
        # Update move_coords for next iteration
        if cycle < n_cycle - 1:  # Don't need to update after the last iteration
            move_coords = np.array(cmd.get_coords(move_sele))
    
    print(f"Aligned {move_sele} to {ref_sele} with initial RMSD: {first_rmsd:.4f} after {n_cycle} cycle(s)")

cmd.extend('align_by_sele', align_by_sele)


def calcu_RRCS(model: str):
    from lazydock.pml.rrcs import calcu_RRCS
    df = calcu_RRCS(model)
    path = os.path.abspath(f'{model}_RRCS.xlsx')
    df.to_excel(path)
    print(f'RRCS saved to {path}')

cmd.extend('calcu_RRCS', calcu_RRCS)


def apply_shader_from_interaction_df(df_path: str, obj_name: str, sum_axis: str = 0, cmap: str = 'coolwarm', alpha_mode: str = None,
                                     show_cbar: bool = False):
    from lazydock.pml.shader import Shader, ShaderValues
    values = ShaderValues().from_interaction_df(df_path, obj_name, int(sum_axis))
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
    align_by_sele(ref_sele, move_sele, move_all=None, n_cycle=1)
    calcu_RRCS(model: str)
    apply_shader_from_interaction_df(df_path: str, obj_name: str, cmap: str = 'coolwarm', alpha_mode: str ='cartoon_transparency')
    apply_shader_from_df(df_path: str, chain_col: str, resi_col: str, c_value_col: str, obj_name: str, cmap: str = 'coolwarm', alpha_mode: str ='cartoon_transparency', save_cbar: bool = False)
''')