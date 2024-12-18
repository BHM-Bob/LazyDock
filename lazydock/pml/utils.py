'''
Date: 2024-12-18 18:35:52
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-18 18:43:31
Description: 
'''
import os

from pymol import cmd

if __name__ == '__main__':
    from lazydock.utils import uuid4
else:
    from .utils import uuid4

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


if __name__ == '__main__':
    print(get_seq('data_tmp/pdb/RECEPTOR.pdb', False))