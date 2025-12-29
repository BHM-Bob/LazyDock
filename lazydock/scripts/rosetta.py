'''
Date: 2025-02-20 10:00:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-20 10:00:00
Description: Prepare ligand files for docking
'''
import argparse
import os
from pathlib import Path
from typing import List, Optional

from tqdm import tqdm
import pandas as pd
import pyrosetta
from mbapy_lite.base import put_err
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web import TaskPool

from lazydock.pyrt.energy_utils import calcu_interface_energy
from lazydock.scripts._script_utils_ import Command, excute_command


class cacl_energy(Command):
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
        # init pyrosetta
        pyrosetta.init(silent=False)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--batch-dir', nargs='+', type=str, default=['.'],
                          help='Batch directory path. default: current directory.')
        args.add_argument('-n', '--name', type=str, default='',
                          help='pdb file name substring. default: empty.')
        args.add_argument('-r', '--receptor-chain', type=str, required=True,
                          help='Receptor chain ID. required.')
        args.add_argument('-l', '--ligand-chain', type=str, required=True,
                          help='Ligand chain ID. required.')
        args.add_argument('-sfxn', '--scorefxn-name', type=str, default='ref2015',
                          help='Score function name. default: ref2015.')
        args.add_argument('-o', '--output', type=str, default='pyrosetta_energy.csv',
                          help='Output CSV file path. default: pyrosetta_energy.csv')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='Number of workers. default: 1.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
    
    def main_process(self):
        pdb_paths = get_paths_with_extension(self.args.batch_dir, ['.pdb'], name_substr=self.args.name)
        df = pd.DataFrame(columns=['pdb_path', 'energy'])
        # parallel
        if self.args.n_workers > 1:
            pool = TaskPool('process', self.args.n_workers, report_error=True).start()
        for pdb_path in tqdm(pdb_paths, desc='Calculating interface energy'):
            if self.args.n_workers > 1:
                pool.add_task(pdb_path, calcu_interface_energy, pdb_path,
                              self.args.receptor_chain, self.args.ligand_chain, self.args.scorefxn_name)
                pool.wait_till_free()
            else:
                energy = calcu_interface_energy(pdb_path, self.args.receptor_chain,
                                                self.args.ligand_chain, self.args.scorefxn_name)
                df.loc[len(df)] = [pdb_path, energy]
        if self.args.n_workers > 1:
            for pdb_path in tqdm(list(pool.tasks.keys()), desc='Querying results from TaskPool'):
                df.loc[len(df)] = [pdb_path, pool.query_task(pdb_path, block=True, timeout=30)]
            pool.close(1)
        df.to_csv(self.args.output, index=False)



_str2func = {
    'calc-energy': cacl_energy,
    # 'relax': relax,
}


def main(sys_args: Optional[List[str]] = None):
    args_paser = argparse.ArgumentParser(description='Prepare ligand files for docking')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))
    
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main()