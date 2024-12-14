'''
Date: 2024-12-13 23:07:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-14 22:08:11
Description: 
'''
import argparse
from pathlib import Path
from typing import List

from mbapy_lite.file import get_paths_with_extension
from pymol import cmd
from tqdm import tqdm

if __name__ == '__main__':
    from lazydock.pml.align_to_axis import align_pose_to_axis
    from lazydock.scripts._script_utils_ import clean_path, show_args
else:
    from ..pml.align_to_axis import align_pose_to_axis
    from ._script_utils_ import clean_path, show_args
    

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'Get pocket residues from ProteinPlus web server.')
    args_paser.add_argument('-d', '--dir', type = str, default='.',
                            help='input dir. Default is %(default)s.')
    args_paser.add_argument('-n', '--name', type = str, default='',
                            help='sub string in file name. Default is %(default)s.')
    args_paser.add_argument('-t', '--type', type = str, default='pdb',
                            help='file type. Default is %(default)s.')
    args_paser.add_argument('-s', '--suffix', type=str, default='aligned_to_axis',
                            help='output suffix. Default is %(default)s.')
    args_paser.add_argument('-r', '--recursive', action='store_true',
                            help='search subdirectories recursively. Default is %(default)s.')
    args = args_paser.parse_args(sys_args)
    # process IO path
    args.dir = clean_path(args.dir)
    args.type = args.type.split('.') if ',' in args.type else [args.type]
    paths = get_paths_with_extension(args.dir, args.type, recursive=args.recursive, name_substr=args.name)
    if not paths:
        raise ValueError(f'No file found in {args.receptor}.')
    # show args
    show_args(args, ['dir', 'type', 'suffix', 'recursive'])
    # sort bonds
    for path in tqdm(paths, total=len(paths)):
        path = Path(path)
        cmd.reinitialize()
        cmd.load(str(path), 'mol')
        align_pose_to_axis('mol')
        cmd.save(f'{path.parent}/{path.stem}_{args.suffix}.{path.suffix}', 'mol')


if __name__ == "__main__":
    main()