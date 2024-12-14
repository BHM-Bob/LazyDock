'''
Date: 2024-11-23 20:03:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-23 20:28:09
Description: 
'''
import argparse
import os
from typing import Dict, List

from mbapy_lite.file import get_paths_with_extension, opts_file

if __name__ == '__main__':
    from lazydock.gmx.thirdparty.sort_mol2_bonds import sort_bonds
    from lazydock.scripts._script_utils_ import clean_path, show_args
else:
    from ..gmx.thirdparty.sort_mol2_bonds import sort_bonds
    from ._script_utils_ import clean_path, show_args
    

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'Sort mol2 bonds by atom index.')
    args_paser.add_argument('-i', '--input', type = str, default='.',
                            help='input mol2 file or dir path. Default is %(default)s.')
    args_paser.add_argument('-s', '--suffix', type = str, default='_sorted',
                            help='suffix of sorted mol2 file name. Default is %(default)s.')
    args_paser.add_argument('-r', '--recursive', action='store_true', default=False,
                            help='FLAG, recursive search. Default is %(default)s.')
    args = args_paser.parse_args(sys_args)
    # process IO path
    args.input = clean_path(args.input)
    if os.path.isfile(args.input):
        paths = [args.input]
    elif os.path.isdir(args.input):
        paths = get_paths_with_extension(args.input, ['.mol2'], recursive=args.recursive)
    else:
        assert False, f'Input path {args.input} is not a file or dir.'
    # show args
    show_args(args, ['input', 'suffix', 'recursive'])
    # sort bonds
    for path in paths:
        print(f'Sorting bonds of {path} ...')
        sorted_mol2 = sort_bonds(opts_file(path))
        sorted_path = os.path.join(os.path.dirname(path), os.path.basename(path).split('.')[0] + args.suffix + '.mol2')
        opts_file(sorted_path, 'w', data=sorted_mol2)
        print(f'Sorted bonds saved to {sorted_path}.')


if __name__ == "__main__":
    main()