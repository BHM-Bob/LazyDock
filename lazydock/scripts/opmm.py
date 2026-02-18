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

import pandas as pd
from mbapy_lite.base import put_err
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web import TaskPool
from tqdm import tqdm

from lazydock.opmm.relax import ForceFieldMinimizer
from lazydock.scripts._script_utils_ import Command, excute_command
from openmm import app as openmm_app


def _relax_worker(pdb_path: str, output_path: str, chain: str, stiffness: float, 
                 max_iter: int, tolerance: int, platform: str, constraints: str, 
                 restrain_backbone: bool):
    # 映射constraints字符串到OpenMM对象
    if constraints == 'hbond':
        constraints_obj = openmm_app.HBonds
    elif constraints == 'none' or constraints is None:
        constraints_obj = None
    else:
        constraints_obj = openmm_app.HBonds  # 默认值
    
    relaxer = ForceFieldMinimizer(
        stiffness=stiffness,
        max_iterations=max_iter,
        tolerance=tolerance,
        platform=platform,
        constraints=constraints_obj  # pyright: ignore[reportArgumentType]
    )
    
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    result_pdb, ret = relaxer(
        pdb_str, 
        output_path, 
        restrain_chain=chain,
        restrain_backbone=restrain_backbone,
        return_info=True
    )
    
    with open(output_path, 'w') as f:
        f.write(result_pdb)
    
    return pdb_path, output_path


class relax(Command):
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--batch-dir', nargs='+', type=str, default=['.'],
                          help='Batch directory path. default: current directory.')
        args.add_argument('-n', '--name', type=str, default='',
                          help='pdb file name substring. default: empty.')
        args.add_argument('-rc', '--restrain-chain', type=str, nargs='+', default=[],
                          help='Chain ID to restrain. default: empty.')
        args.add_argument('-it', '--max-iter', type=int, default=1000,
                          help='Maximum iterations for relaxation. default: 1000.')
        args.add_argument('-s', '--stiffness', type=float, default=10**10,
                          help='Stiffness for restraints. default: 10^10.')
        args.add_argument('-t', '--tolerance', type=int, default=10,
                          help='Tolerance for minimization. default: 10.')
        args.add_argument('-p', '--platform', type=str, default='CUDA',
                          choices=['CUDA', 'CPU'],
                          help='Platform for computation. default: CUDA.')
        args.add_argument('-c', '--constraints', type=str, default='hbond',
                          choices=['hbond', 'none'],
                          help='Constraints type: hbond (HBonds) or none (None). default: hbond.')
        args.add_argument('-rb', '--restrain-backbone', action='store_true', default=False,
                          help='Restrain backbone atoms. default: False.')
        args.add_argument('-o', '--output-suffix', type=str, default='_relaxed',
                          help='Output PDB file suffix. default: _relaxed.')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='Number of workers. default: 1.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
    
    def main_process(self):
        # 处理多个目录
        pdb_paths = []
        for batch_dir in self.args.batch_dir:
            dir_paths = get_paths_with_extension(batch_dir, ['.pdb'], name_substr=self.args.name)
            pdb_paths.extend(dir_paths)
        
        if not pdb_paths:
            self.printf("Warning: No PDB files found in specified directories.")
            return
        
        # 检查是否有约束链
        if not self.args.restrain_chain:
            self.printf("Warning: No restrain chain specified. Using empty list.")
        
        # parallel
        if self.args.n_workers > 1:
            pool = TaskPool('process', self.args.n_workers, report_error=True).start()
        
        # Process each PDB file
        for pdb_path in tqdm(pdb_paths, desc='Relaxing structures'):
            output_path = pdb_path.replace('.pdb', f'{self.args.output_suffix}.pdb')
            
            if self.args.n_workers > 1:
                pool.add_task(pdb_path, _relax_worker, pdb_path, output_path,
                                                     self.args.restrain_chain, self.args.stiffness,
                                                     self.args.max_iter, self.args.tolerance,
                                                     self.args.platform, self.args.constraints,
                                                     self.args.restrain_backbone)
                pool.wait_till_free()
            else:
                # 单线程处理
                try:
                    _relax_worker(pdb_path, output_path, self.args.restrain_chain, 
                                 self.args.stiffness, self.args.max_iter, self.args.tolerance,
                                 self.args.platform, self.args.constraints, self.args.restrain_backbone)
                    self.printf(f"Successfully relaxed: {pdb_path} -> {output_path}")
                except Exception as e:
                    self.printf(f"Error processing {pdb_path}: {str(e)}")
        
        if self.args.n_workers > 1:
            pool.close(1)


_str2func = {
    'relax': relax,
}


def main(sys_args: Optional[List[str]] = None):
    args_paser = argparse.ArgumentParser(description='OpenMM-based structure relaxation')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))
    
    excute_command(args_paser, sys_args, _str2func)  # pyright: ignore[reportArgumentType]


if __name__ == "__main__":
    main()