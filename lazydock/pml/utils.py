'''
Date: 2024-12-18 18:35:52
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-18 18:43:31
Description: 
'''
import os

import numpy as np
from pymol import cmd

if __name__ == '__main__':
    from lazydock.utils import uuid4
else:
    from ..utils import uuid4

def get_seq(pose: str, fasta: bool = False):
    if os.path.isfile(pose):
        if fasta:
            pose_name = os.path.basename(pose).split('.')[0]
        else:
            pose_name = uuid4()
        cmd.load(pose, pose_name)
        pose = pose_name
    seq = cmd.get_fastastr(pose)
    if fasta:
        return seq
    return ''.join(seq.split('\n')[1:])


def rotation_matrix_to_pymol_matrix(rotation_matrix, ref_center, move_center):
    """
    Convert a 3x3 rotation matrix to PyMOL's 4x4 transformation matrix format.
    
    Parameters:
        rotation_matrix (np.ndarray): 3x3 rotation matrix
        ref_center (np.ndarray): Center of reference coordinates [x, y, z]
        move_center (np.ndarray): Center of moving coordinates [x, y, z]
        
    Returns:
        list: PyMOL-specific 4x4 matrix in the format expected by cmd.transform_selection
    """
    # Calculate translation needed after rotation
    # We want: ref_center = rotation_matrix * (move_center + pre_translation) + post_translation
    # For simplicity, we'll use post_translation only
    post_translation = ref_center - np.dot(rotation_matrix, move_center)
    
    # Create PyMOL-specific 4x4 matrix (homogenous=0 format)
    # [m0 m1 m2 m3]
    # [m4 m5 m6 m7]
    # [m8 m9 m10 m11]
    # [m12 m13 m14 m15]
    pymol_matrix = [
        rotation_matrix[0,0], rotation_matrix[0,1], rotation_matrix[0,2], post_translation[0],
        rotation_matrix[1,0], rotation_matrix[1,1], rotation_matrix[1,2], post_translation[1],
        rotation_matrix[2,0], rotation_matrix[2,1], rotation_matrix[2,2], post_translation[2],
        0, 0, 0, 1
    ]
    
    return pymol_matrix


if __name__ == '__main__':
    print(get_seq('data_tmp/pdb/RECEPTOR.pdb', False))