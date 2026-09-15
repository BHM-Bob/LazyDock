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
from pymol import cmd
from tqdm import tqdm

from lazydock.opmm.relax import ForceFieldMinimizer
from lazydock.scripts._script_utils_ import Command, excute_command
from openmm import app as openmm_app


def _relax_worker(pdb_path: str, output_path: str, chain: str, stiffness: float, 
                 max_iter: int, tolerance: int, platform: str, constraints: str, 
                 restrain_backbone: bool, cyclic_chains: list = None, device_index=None,
                 disulfide_chains: list = None):
    # 映射constraints字符串到OpenMM对象
    if constraints == 'hbond':
        constraints_obj = openmm_app.HBonds
    elif constraints == 'all':
        constraints_obj = openmm_app.AllBonds
    elif constraints == 'none' or constraints is None:
        constraints_obj = None
    else:
        constraints_obj = openmm_app.HBonds  # 默认值
    
    relaxer = ForceFieldMinimizer(
        stiffness=stiffness,
        max_iterations=max_iter,
        tolerance=tolerance,
        platform=platform,
        constraints=constraints_obj,  # pyright: ignore[reportArgumentType]
        cyclic_chains=cyclic_chains,
        disulfide_chains=disulfide_chains,
        device_index=device_index,
    )
    
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    result_pdb, _ = relaxer(
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
        args.add_argument('-relax', '--relax-chain', type=str, nargs='+', default=None,
                          help='Relaxed chain. If given, will set restrain-chain to others.')
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
                          choices=['hbond', 'all', 'none'],
                          help='Constraints type: hbond (HBonds), all (All Bonds) or none (None). default: hbond.')
        args.add_argument('-rb', '--restrain-backbone', action='store_true', default=False,
                          help='Restrain backbone atoms. default: False.')
        args.add_argument('-o', '--output-suffix', type=str, default='_relaxed',
                          help='Output PDB file suffix. default: _relaxed.')
        args.add_argument('--cyclic-chains', type=str, nargs='+', default=None,
                          help='Cyclic peptide chains to keep the head-tail peptide bond during relaxation, '
                               'e.g. --cyclic-chains P, default is %(default)s.')
        args.add_argument('--disulfide-chains', type=str, nargs='+', default=None,
                          help='Chains containing disulfide bonds to protect during relaxation '
                               '(SG-SG bond via CHARMM36 DISU patch or CustomBondForce), '
                               'e.g. --disulfide-chains P, default is %(default)s.')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='Number of workers. default: 1.')
        args.add_argument('--gpus', type=int, nargs='+', default=None,
                          help='GPU device ids to use for parallel relaxation, e.g. --gpus 0 1, '
                               'default is %(default)s (use default device).')
        args.add_argument('--n-task-per-gpu', type=int, default=1,
                          help='max number of concurrent tasks on one GPU, default is %(default)s.')
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
        if not self.args.restrain_chain and not self.args.relax_chain:
            self.printf("Warning: No restrain chain specified. Using empty list.")
        
        # parallel
        pool = TaskPool('process', self.args.n_workers, report_error=True).start()
        
        # GPU 槽位: 把 --gpus 展开成 [g0 x n_task_per_gpu, g1 x n_task_per_gpu, ...],
        # 任务按 index 轮询分配, 保证单卡可承载多个并发任务
        gpu_slots = []
        if self.args.gpus:
            for g in self.args.gpus:
                gpu_slots.extend([g] * max(1, self.args.n_task_per_gpu))
        
        # Process each PDB file
        for task_i, pdb_path in enumerate(tqdm(pdb_paths, desc='Relaxing structures')):
            output_path = pdb_path.replace('.pdb', f'{self.args.output_suffix}.pdb')
            if self.args.relax_chain:
                cmd.reinitialize()
                cmd.load(pdb_path)
                all_chains = cmd.get_chains('all')
                self.args.restrain_chain = list(set(all_chains) - set(self.args.relax_chain))
            device_index = gpu_slots[task_i % len(gpu_slots)] if gpu_slots else None
            pool.add_task(pdb_path, _relax_worker, pdb_path, output_path,
                                                    self.args.restrain_chain, self.args.stiffness,
                                                    self.args.max_iter, self.args.tolerance,
                                                    self.args.platform, self.args.constraints,
                                                    self.args.restrain_backbone,
                                                    self.args.cyclic_chains,
                                                    device_index,
                                                    self.args.disulfide_chains)
            pool.wait_till_free()
        pool.wait_till_all_done()
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