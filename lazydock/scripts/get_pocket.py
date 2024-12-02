'''
Date: 2024-11-27 17:01:54
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-02 17:00:54
Description: 
'''
import argparse
import os
import shutil
from pathlib import Path
from typing import Dict, List

from mbapy_lite.file import get_paths_with_extension, opts_file
from pymol import cmd

if __name__ == '__main__':
    from lazydock.scripts._script_utils_ import clean_path, show_args
    from lazydock.web.gen_pocket import (get_pocket_box_from_ProteinPlus,
                                         parse_pocket_box_from_ProteinPlus)
else:
    from ..web.gen_pocket import (get_pocket_box_from_ProteinPlus,
                                  parse_pocket_box_from_ProteinPlus)
    from ._script_utils_ import clean_path, show_args
    

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'Get pocket residues from ProteinPlus web server.')
    args_paser.add_argument('-r', '--receptor', type = str, default='.',
                            help='input mol2 file or dir path. Default is %(default)s.')
    args_paser.add_argument('-l', '--ligand', type = str, default=None,
                            help='ligand sdf file path. Default is %(default)s.')
    args_paser.add_argument('-o', '--output', type=str, default=None,
                            help='output dir. Default is %(default)s.')
    args_paser.add_argument('-m', '--method', type=str, default='extend', choices=['extend', 'mean'],
                            help='pocket calculation method. Default is %(default)s.')
    args = args_paser.parse_args(sys_args)
    # process IO path
    args.receptor = clean_path(args.receptor)
    if os.path.isfile(args.receptor):
        receptor_paths = [args.receptor]
    elif os.path.isdir(args.receptor):
        receptor_paths = get_paths_with_extension(args.receptor, ['.pdb'], recursive=args.recursive)
    else:
        raise ValueError(f'Input path {args.receptor} is not a pdb file or directory.')
    if not receptor_paths:
        raise ValueError(f'No pdb file found in {args.receptor}.')
    # show args
    show_args(args, ['receptor', 'ligand', 'output'])
    # sort bonds
    for receptor_path in receptor_paths:
        print(f'Getting pocket box from {receptor_path} ...')
        output_dir = Path(receptor_path).parent / f'{Path(receptor_path).stem}_pocket_box'
        if output_dir.exists():
            print(f'Output dir {output_dir} already exists, skip.')
            continue
        output_dir.mkdir(exist_ok=True, parents=True)
        shutil.copy(receptor_path, output_dir) # so that the zip file is in the same dir with pdb file
        # get pocket box
        get_pocket_box_from_ProteinPlus(os.path.join(output_dir, Path(receptor_path).name), ligand_path=args.ligand)
        # parse pocket box
        zip_path = get_paths_with_extension(output_dir, ['.zip'], recursive=False)[0]
        for index in [[0], [1], [2], [0, 1], [0, 1, 2]]:
            idx_str = ','.join(map(str, index))
            if args.method == 'extend':
                pocket = parse_pocket_box_from_ProteinPlus(zip_path, index, True, method=args.method)
            elif args.method == 'mean':
                df, pocket = parse_pocket_box_from_ProteinPlus(zip_path, index, True, method=args.method)
                df.to_excel(output_dir / f"pocket_{idx_str}_box.xlsx", index=False)
            # save pocket residues
            opts_file(output_dir / f"pocket_{idx_str}_residues.json", 'w', way='json', data=pocket)
            cmd.save(output_dir / f"pocket_{idx_str}_box.pse")
            print(f'Pocket residues saved to {output_dir}.')


if __name__ == "__main__":
    main()