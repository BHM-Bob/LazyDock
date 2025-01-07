import itertools
from typing import List, Union

import numpy as np
from compas.geometry import oriented_bounding_box_numpy
from mbapy.base import put_err, split_list
from pymol import cmd, editing


def sort_vertices(vertices, sort_orders=[(0, 1, 2), (1, 2, 0), (2, 0, 1)], tolerance=0.1):
    """
    对顶点进行排序，以匹配plot_bounding_box中的顺序。
    
    参数:
    vertices (numpy.ndarray): 8个顶点的坐标，形状为 (8, 3)。
    sort_orders (list of tuples): 排序顺序的列表，每个元组包含三个索引，表示XYZ的排序顺序。
    
    返回:
    list: [P0, P0-x, ...] 或 None（如果所有排序顺序都失败）。
    """
    neighbors = {
        0: [1, 2, 4], 1: [0, 3, 5], 2: [0, 3, 6], 3: [1, 2, 7], 4: [0, 5, 6], 5: [1, 4, 7], 6: [2, 4, 7], 7: [3, 5, 6]
    }
    def is_orthogonal(v1, v2, tolerance=1e-6):
        """检查两个向量是否正交，包含指定的误差容许"""
        _v1, _v2 = v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)
        return abs(np.dot(_v1, _v2)) < tolerance

    def check_orthogonality(points, tolerance=1e-6):
        """检查所有顶点的三条棱是否正交"""
        for i in range(8):
            vectors = []
            for j in neighbors[i]:
                vectors.append(np.array(points[j]) - np.array(points[i]))
            for v1, v2 in itertools.combinations(vectors, 2):
                if not is_orthogonal(v1, v2, tolerance):
                    return False
        return True

    for sort_order in sort_orders:
        vertices_copy = vertices.copy()
        for i, j in zip(sort_order, [4, 2, 1]):
            tmp = []
            for sub in split_list(vertices_copy, j):
                tmp.extend(sorted(sub, key=lambda v: v[i]))
            vertices_copy = tmp
        
        if check_orthogonality(vertices_copy, tolerance):
            return vertices_copy
    
    raise ValueError("No valid sorting order found for the given vertices.")


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


def normalize(v):
    """归一化向量"""
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def angle_between(v1, v2):
    """计算两个向量之间的夹角（弧度）"""
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.arccos(np.clip(cos_theta, -1.0, 1.0))

def rotation_matrix_from_vectors(v1, v2):
    """计算从向量v1旋转到向量v2的旋转矩阵"""
    v1 = normalize(v1)
    v2 = normalize(v2)
    v = np.cross(v1, v2)
    s = np.linalg.norm(v)
    c = np.dot(v1, v2)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    if s == 0:
        # 如果v1和v2同向或反向，不需要旋转
        return np.eye(3)
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
    return R

def rotation_matrix_x(angle):
    """绕X轴旋转的旋转矩阵"""
    c = np.cos(np.radians(angle))
    s = np.sin(np.radians(angle))
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])

def rotation_matrix_y(angle):
    """绕Y轴旋转的旋转矩阵"""
    c = np.cos(np.radians(angle))
    s = np.sin(np.radians(angle))
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])

def rotation_matrix_z(angle):
    """绕Z轴旋转的旋转矩阵"""
    c = np.cos(np.radians(angle))
    s = np.sin(np.radians(angle))
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])

def apply_rotation(coords, angles):
    """应用旋转到顶点"""
    # 构建旋转矩阵
    Rx = rotation_matrix_x(angles[0])
    Ry = rotation_matrix_y(angles[1])
    Rz = rotation_matrix_z(angles[2])
    
    # 组合旋转矩阵
    R = np.dot(Rz, np.dot(Ry, Rx))
    
    # 应用旋转矩阵到每个顶点
    rotated_vertices = np.dot(coords, R.T)
    return rotated_vertices

def align_edge_with_x_axis(v0, v1):
    """计算将边v0到v1与X轴对齐所需的旋转角度"""
    edge1 = np.array(v1) - np.array(v0)
    x_axis = np.array([1, 0, 0])
    R = rotation_matrix_from_vectors(edge1, x_axis)
    # 计算欧拉角
    sy = np.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0])
    singular = sy < 1e-6
    if not singular:
        x = np.arctan2(R[2, 1], R[2, 2])
        y = np.arctan2(-R[2, 0], sy)
        z = np.arctan2(R[1, 0], R[0, 0])
    else:
        x = np.arctan2(-R[1, 2], R[1, 1])
        y = np.arctan2(-R[2, 0], sy)
        z = 0
    # 将弧度转换为度
    euler_angles = np.degrees([x, y, z])
    return euler_angles

def align_edge_with_y_axis_after_x_alignment(v0, v2):
    """计算在v0-v1与X轴对齐后，将v0-v2与Y轴对齐所需的绕X轴旋转角度"""
    edge2 = np.array(v2) - np.array(v0)
    # 只考虑Y和Z分量，因为X分量在与X轴对齐后为0
    yz_projection = edge2[1:]  # 取Y和Z分量
    y_axis = np.array([1, 0])  # Y轴在YZ平面上的投影
    angle = np.arctan2(yz_projection[1], yz_projection[0])  # 计算绕X轴旋转角度
    return np.degrees(angle)

def rotate_bounding_box_to_axis(bounding_box_vertices: np.ndarray):
    """
    计算包围盒绕轴旋转的角度，使包围盒的各边与坐标轴对齐。
    
    参数:
    bounding_box_vertices： 包围盒的8个顶点坐标，形状为 (8, 3)。
    
    返回:
    numpy.ndarray: 旋转角度。
    """
    bounding_box_vertices = np.array(bounding_box_vertices)
    angles1 = align_edge_with_x_axis(bounding_box_vertices[0], bounding_box_vertices[1])
    vertices_rotated = apply_rotation(bounding_box_vertices, angles1)
    angle2 = align_edge_with_y_axis_after_x_alignment(vertices_rotated[0], vertices_rotated[2])
    return angles1.tolist(), [angle2, 0, 0]


def align_pose_to_axis(pml_name: str, move_name: str = None, fixed: Union[List[float], str] = 'center', state: int = 0,
                       move_method: str = 'rotate', dss: bool = False, quite: int = 1):
    """
    Parameters:
        - pml_name (str): pymol object name to calculate minimum bounding box and tranformation matrix.
        - move_name (str): pymol object name to move to aligned position.
        - fixed (Union[List[float], str]): fixed point to align with, not works when move_method is 'rotate'.
        - state (int): pml state to calculate minimum bounding box.
        - move_method (str): method to move to aligned position, support 'transform', 'alter' and 'rotate'.
        - dss (bool): whether to apply dss to pml object.
        - quite (int): whether to print log information.
    
    Returns:
        - aligned_coords (numpy.ndarray): aligned coordinates.
        - aligned_box (numpy.ndarray): aligned bounding box.
        - rotation_matrix (numpy.ndarray): rotation matrix to align with.
        - fixed_coords (numpy.ndarray): fixed point to align with.

    Notes:
        - TODO: RMS is rigth(=0), but second structure all trun to raandom coil, fixed by rotate method.
        - TODO: in rotate method, apply angles1 gets right result, but apply angles2 gets wrong result, pymol is fine.
    """
    move_name = move_name or pml_name
    # get coords
    coords, index2coords, sorted_vertices = calcu_bounding_box(pml_name, state=state)
    if isinstance(fixed, str) and fixed != 'center':
        fixed_coords = []
        cmd.iterate_state(state, fixed, 'fixed_coords.append([x, y, z])', space=locals())
        fixed_coords = np.array(fixed_coords).mean(axis=0)
    else:
        fixed_coords = fixed
    if move_method in {'transform', 'alter'}:
        # align bounding box
        aligned_coords, aligned_box, rotation_matrix, fixed_coords = align_bounding_box_to_axis(coords, sorted_vertices, fixed_coords=fixed_coords)
        # create pymol rotation matrix
        pml_mat = np.zeros((4, 4))
        pml_mat[:3, :3] = rotation_matrix.T # np is matmul, but pymol is dot product
        pml_mat[-1, :] = list(-fixed_coords) + [1]
        pml_mat[:, -1] = list(fixed_coords) + [1]
        if not quite:
            print(f'pymol transform matrix: {pml_mat.flatten().tolist()}')
    elif move_method == 'rotate':
        # calculate rotation degree
        angles1, angles2 = rotate_bounding_box_to_axis(sorted_vertices)
        if not quite:
            print(f'Rotate angles1: {angles1}, angles2: {angles2}')
    else:
        return put_err(f'Unsupported move_method: {move_method}, only support transform and alter, skip transform.')
    # move to aligned position
    if move_method == 'transform':
        cmd.transform_selection(move_name, pml_mat.flatten().tolist(), homogenous=0)
    elif move_method == 'alter':
        cmd.alter_state(state, move_name, 'x = aligned_coords[index2coords[index], 0]', space=locals())
        cmd.alter_state(state, move_name, 'y = aligned_coords[index2coords[index], 1]', space=locals())
        cmd.alter_state(state, move_name, 'z = aligned_coords[index2coords[index], 2]', space=locals())
        cmd.sort(move_name)
    elif move_method == 'rotate':
        for angles in [angles1, angles2]:
            for axis, angle in zip(['x', 'y', 'z'], angles):
                if angle != 0:
                    cmd.rotate(axis, angle, move_name)
            aligned_coords = apply_rotation(coords, angles)
            aligned_box = apply_rotation(sorted_vertices, angles)
        rotation_matrix = None
    if (isinstance(dss, str) and int(dss)) or (isinstance(dss, int) and dss):
        print('Applying DSS to aligned object, it may has bad effect on some complex structure.')
        editing.dss(move_name)
    cmd.rebuild(move_name)
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
    # 标记每个点
    for i in range(8):
        ax.text(vertices[i][0], vertices[i][1], vertices[i][2], str(i))


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    from lazydock.pml.thirdparty.draw_bounding_box import draw_bounding_box
    
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    cmd.select('receptor_CA', 'name CA and receptor')
    coords, _, vertics = calcu_bounding_box('receptor_CA')
    aligned_coords, aligned_box, rotation_matrix, fixed_coords = align_bounding_box_to_axis(coords, vertics)
    _, _, aligned_vertics = calcu_bounding_box(coords=aligned_coords)
    
    aligned_coords, aligned_box, rotation_matrix, fixed_coords = align_pose_to_axis('receptor', move_method='rotate')
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
    