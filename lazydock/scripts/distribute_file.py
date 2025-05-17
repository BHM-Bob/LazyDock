'''
Date: 2024-12-13 20:18:59
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-05-17 12:23:43
Description: steps most from http://www.mdtutorials.com/gmx
'''

import argparse
import os
import shutil
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import get_dir, get_paths_with_extension
from tqdm import tqdm

from lazydock.scripts._script_utils_ import clean_path, show_args


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'Get pocket residues from ProteinPlus web server.')
    args_paser.add_argument('-sr', '--source-root', type = str,
                            help='root to the files or directory to be distributed.')
    args_paser.add_argument('-sn', '--source-name', type = str,
                            help='sub string in source files or directory name.')
    args_paser.add_argument('-st', '--source-type', type = str, nargs='+', default=[],
                            help='file type in source directory.')
    args_paser.add_argument('-m', '--mode', type = str, choices=['file', 'dir'], default='file',
                            help='source search mode, file or directory. Default is %(default)s.')
    args_paser.add_argument('-dr', '--dist-root', type = str,
                            help='path to the dist root directory.')
    args_paser.add_argument('-dn', '--dist-name', type = str,
                            help='sub string in dist directory name.')
    args_paser.add_argument('-dt', '--dist-type', type = str, nargs='+', default=[],
                            help='file type in dist directory. Default is %(default)s.')
    args_paser.add_argument('-in', '--item-name', type = str, default='',
                            help='sub string in dist directory item files name. Default is %(default)s.')
    args_paser.add_argument('-it', '--item-type', type = str, default='',
                            help='file type in dist directory. Default is %(default)s.')
    args_paser.add_argument('-n', '--new-name', type = str, default=None,
                            help='sub string in file name. Default is %(default)s.')
    args = args_paser.parse_args(sys_args)
    # process IO path
    args.source_root = clean_path(args.source_root)
    args.dist_root = clean_path(args.dist_root)
    # show args
    show_args(args, ['source_root', 'source_name', 'source_type', 'mode', 'dist_root',
                     'dist_name', 'dist_type', 'item_name', 'new_name'])
    # search source files or directories
    if args.mode == 'file':
        sources = get_paths_with_extension(args.source_root, args.source_type, name_substr=args.source_name)
    else:
        sources = get_dir(args.source_root, 0, None, args.source_type, True, args.source_name)
    if not sources:
        return put_err(f'No dist directory found in {args.source_root}, exit.')
    # search dist directories
    dist_dirs = get_dir(args.dist_root, file_extensions=args.dist_type, recursive=True, dir_name_substr=args.dist_name)
    if not dist_dirs:
        return put_err(f'No dist directory found in {args.dist_root}, exit.')
    # copy file or directory to dist directory
    for dist_dir in tqdm(dist_dirs, total=len(dist_dirs)):
        for source in sources:
            dist_name = os.path.basename(source) if args.new_name is None else args.new_name
            dist_path = os.path.join(dist_dir, dist_name)
            if args.mode == 'file':
                shutil.copy(source, dist_path)
            else:
                shutil.copytree(source, dist_path)


if __name__ == "__main__":
    main()