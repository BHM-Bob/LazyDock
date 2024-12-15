from typing import List, Union

import numpy as np
from compas.geometry import oriented_bounding_box_numpy
from mbapy.base import put_err, split_list
from pymol import cmd, editing


def sort_vertices(vertices):
    """
    对顶点进行排序，以匹配plot_bounding_box中的顺序。
    
    参数:
    vertices (numpy.ndarray): 8个顶点的坐标，形状为 (8, 3)。
    
    返回:
    list: 排序后的顶点列表。
    """
    for i, j in zip(range(3), [8, 4, 2]):
        tmp = []
        for sub in split_list(vertices, j):
            tmp.extend(sorted(sub, key=lambda v: v[i]))
        vertices = tmp
    return vertices


def calcu_bounding_box(pml_name: str = None, coords: np.ndarray = None, state: int = 0):
    index2coords = {}
    if coords is None:
        cmd.iterate_state(state, pml_name, 'index2coords[index] = [x, y, z]', space=locals())
        coords = np.array(list(index2coords.values()))
        index2coords = {index: i for i, index in enumerate(list(index2coords.keys()))}
    if not coords.any():
        return put_err('No coordinates found in selection, return None.', None)
    bounding_box_vertices = oriented_bounding_box_numpy(coords)
    sorted_vertices = sort_vertices(bounding_box_vertices)
    return coords, index2coords, sorted_vertices


def align_bounding_box_to_axis(coords: np.ndarray, bounding_box_vertices: np.ndarray,
                               fixed_coords: Union[List[float], str] = 'center'):
    """
    以指定的点为不动点对点云进行旋转，使其包围盒的各边与坐标轴对齐。
    
    参数:
    coords (numpy.ndarray): 点云数据，形状为 (N, 3)。
    fixed_coords (numpy.ndarray): 不动点坐标，形状为 (3,)。
    
    返回:
    numpy.ndarray: 旋转后的点云数据。
    """
    if isinstance(fixed_coords, str):
        if fixed_coords == 'center':
            fixed_coords = coords.mean(axis=0)
        else:
            return put_err(f'Unsupported fixed_point value: {fixed_coords}, return original vertices.', bounding_box_vertices)
    bounding_box_vertices = np.array(bounding_box_vertices)
    # Select the first vertex and its three adjacent vertices
    v0 = bounding_box_vertices[0]
    v1 = bounding_box_vertices[1]
    v2 = bounding_box_vertices[2]
    v4 = bounding_box_vertices[4]
    # Compute the three principal axes from the bounding box edges
    axis_x = v1 - v0
    axis_y = v2 - v0
    axis_z = v4 - v0
    # Normalize the axes to get unit vectors
    axis_x = axis_x / np.linalg.norm(axis_x)
    axis_y = axis_y / np.linalg.norm(axis_y)
    axis_z = axis_z / np.linalg.norm(axis_z)
    # Construct the rotation matrix from the principal axes
    rotation_matrix = np.column_stack((axis_x, axis_y, axis_z))
    # Rotate the point cloud around the fixed point
    rotated_points = np.dot(coords - fixed_coords, rotation_matrix) + fixed_coords
    rotated_box = np.dot(bounding_box_vertices - fixed_coords, rotation_matrix) + fixed_coords
    return rotated_points, rotated_box, rotation_matrix, fixed_coords


def align_pose_to_axis(pml_name: str, fixed: Union[List[float], str] = 'center', state: int = 0, move_method: str = 'transform', dss: bool = True):
    """
    TODO: RMS is rigth(=0), but second structure all trun to raandom coil.
    """
    # get coords
    coords, index2coords, sorted_vertices = calcu_bounding_box(pml_name, state=state)
    if isinstance(fixed, str) and fixed != 'center':
        fixed_coords = []
        cmd.iterate_state(state, fixed, 'fixed_coords.append([x, y, z])', space=locals())
        fixed_coords = np.array(fixed_coords).mean(axis=0)
    else:
        fixed_coords = fixed
    # align bounding box
    aligned_coords, aligned_box, rotation_matrix, fixed_coords = align_bounding_box_to_axis(coords, sorted_vertices, fixed_coords=fixed_coords)
    # create pymol rotation matrix
    pml_mat = np.zeros((4, 4))
    pml_mat[:3, :3] = rotation_matrix.T # np is matmul, but pymol is dot product
    pml_mat[-1, :] = list(-fixed_coords) + [1]
    pml_mat[:, -1] = list(fixed_coords) + [1]
    # move to aligned position
    if move_method == 'transform':
        cmd.transform_selection(pml_name, pml_mat.flatten().tolist(), homogenous=1)
    elif move_method == 'alter':
        cmd.alter_state(state, pml_name, 'x = aligned_coords[index2coords[index], 0]', space=locals())
        cmd.alter_state(state, pml_name, 'y = aligned_coords[index2coords[index], 1]', space=locals())
        cmd.alter_state(state, pml_name, 'z = aligned_coords[index2coords[index], 2]', space=locals())
        cmd.sort(pml_name)
    else:
        put_err(f'Unsupported move_method: {move_method}, only support transform and alter, skip transform.')
    if dss:
        editing.dss(pml_name)
    cmd.rebuild(pml_name)
    print(f'pymol transform matrix: {pml_mat.flatten().tolist()}')
    return aligned_coords, aligned_box, rotation_matrix, fixed_coords

def set_axes_equal(ax):
    """
    设置坐标轴的比例相等。
    
    参数:
    ax (Axes3D): 3D坐标轴。
    """
    xlim = ax.get_xlim3d()
    ylim = ax.get_ylim3d()
    zlim = ax.get_zlim3d()

    x_range = abs(xlim[1] - xlim[0])
    x_middle = np.mean(xlim)
    y_range = abs(ylim[1] - ylim[0])
    y_middle = np.mean(ylim)
    z_range = abs(zlim[1] - zlim[0])
    z_middle = np.mean(zlim)

    # 找到最大的范围
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def plot_bounding_box(ax, vertices, color='r'):
    """
    绘制长方体的边。
    
    参数:
    ax (Axes3D): 3D坐标轴。
    vertices (list): 长方体的8个顶点坐标。
    """
    # 定义长方体的12条边
    edges = [
        (0, 1), (1, 3), (3, 2), (2, 0),  # 底部面
        (4, 5), (5, 7), (7, 6), (6, 4),  # 顶部面
        (0, 4), (1, 5), (2, 6), (3, 7)   # 连接上下两个面
    ]
    
    # 绘制每条边
    for edge in edges:
        ax.plot(*zip(vertices[edge[0]], vertices[edge[1]]), linestyle='-', color=color)


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    from lazydock.pml.thirdparty.draw_bounding_box import draw_bounding_box
    
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    cmd.select('receptor_CA', 'name CA and receptor')
    coords, _, vertics = calcu_bounding_box('receptor_CA')
    aligned_coords, aligned_box, rotation_matrix, fixed_coords = align_bounding_box_to_axis(coords, vertics)
    _, _, aligned_vertics = calcu_bounding_box(coords=aligned_coords)
    
    aligned_coords, aligned_box, rotation_matrix, fixed_coords = align_pose_to_axis('receptor', move_method='alter')
    draw_bounding_box('receptor')
    cmd.save('data_tmp/pdb/RECEPTOR_bounding_box.pse')
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(aligned_coords[:, 0], aligned_coords[:, 1], aligned_coords[:, 2], c='r', alpha=0.5, marker='o', s=10, label='Original Points')
    plot_bounding_box(ax, aligned_box, color='r')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    set_axes_equal(ax)
    ax.legend()
    plt.show()
    