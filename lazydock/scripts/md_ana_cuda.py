import argparse
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import MDAnalysis as mda
import numpy as np
import torch
from mbapy_lite.base import put_err, put_log
from tqdm import tqdm

from lazydock.gmx.mda.align import get_aligned_coords
from lazydock.gmx.mda.utils import filter_atoms_by_chains
from lazydock.scripts._script_utils_ import (clean_path, excute_command,
                                             process_batch_dir_lst)
from lazydock.scripts.ana_gmx import mmpbsa


class trjalign(mmpbsa):
    """
    对MD轨迹进行对齐并输出到新文件。
    接受align链和apply链参数，基于align链进行对齐，将变换应用到apply链，输出对齐后的轨迹。
    使用CUDA加速计算旋转矩阵。
    """
    HELP = """
    Align MD trajectory based on specified chains and output to a new file.

    This tool performs trajectory alignment using specified align chains as reference,
    applies the transformation to apply chains, and writes the aligned trajectory
    to a new file. Uses CUDA for accelerated computation of rotation matrix.

    Example:
        python md_ana_cuda.py trjalign -d ./simulations -align A B -apply C D -o aligned.xtc
    """

    def __init__(self, args, printf=print):
        super().__init__(args, printf)

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files. Default is %(default)s.")
        args.add_argument('-top', '--top-name', type=str, default='md.tpr',
                          help='topology file name in each sub-directory, such as md.tpr. Default is %(default)s.')
        args.add_argument('-traj', '--traj-name', type=str, default='md_center.xtc',
                          help='trajectory file name in each sub-directory, such as md_center.xtc. Default is %(default)s.')
        args.add_argument('-align', '--align-chains', type=str, nargs='+', required=True,
                          help='Chain(s) to use for alignment reference. Multiple chains can be specified. Required.')
        args.add_argument('-apply', '--apply-chains', type=str, nargs='+', required=True,
                          help='Chain(s) to apply the alignment transformation to. Multiple chains can be specified. Required.')
        args.add_argument('-sele', '--select', type=str, default='protein',
                          help='Selection string for alignment atoms. Default is %(default)s.')
        args.add_argument('-o', '--output', type=str, default='aligned.xtc',
                          help='Output trajectory file name. Default is %(default)s.')
        args.add_argument('-b', '--begin-frame', type=int, default=0,
                          help='First frame to start the analysis. Default is %(default)s.')
        args.add_argument('-e', '--end-frame', type=int, default=None,
                          help='Last frame to end the analysis. Default is %(default)s.')
        args.add_argument('-step', '--traj-step', type=int, default=1,
                          help='Step while reading trajectory. Default is %(default)s.')
        args.add_argument('-F', '--force', default=False, action='store_true',
                          help='Force to re-run the analysis even if output exists. Default is %(default)s.')
        return args

    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)

    def align(self, u: mda.Universe, args: argparse.ArgumentParser,
              align_select: Optional[str] = None, align_chains: Optional[List[str]] = None):
        """
        Align the trajectory based on specified chains.
        """
        # Check CUDA availability
        device_str = 'cuda' if torch.cuda.is_available() else 'cpu'
        backend_str = 'cuda' if device_str == 'cuda' else 'torch'
        if device_str == 'cpu':
            put_log('Warning: CUDA not available, falling back to CPU.')
        else:
            put_log(f'Using CUDA device: {torch.cuda.get_device_name(0)}')

        # Get align atom group (for calculating alignment transformation)
        if not align_select and not args.align_chains:
            align_ag = u.atoms
        else:
            align_ag = u.select_atoms(align_select or args.select)
        if args.align_chains is not None or align_chains is not None:
            align_ag = filter_atoms_by_chains(align_ag, align_chains or args.align_chains)
        put_log(f'Align selection: {len(align_ag)} atoms from chains {args.align_chains}')
        
        # Process each frame
        start, step, stop = args.begin_frame, args.traj_step, args.end_frame
        ori_trj, aligned_trj, rmsd, rot = get_aligned_coords(u, align_ag, start, step, stop,
                                                             backend=backend_str, return_rmsd=True, return_rot=True)
        return align_ag, ori_trj, aligned_trj, rmsd, rot
    
    def analysis(self, u: mda.Universe, w_dir: Path, args: argparse.ArgumentParser):
        output_path = w_dir / args.output

        # Check if output already exists
        if output_path.exists() and not args.force:
            return put_log(f'Aligned trajectory {output_path} already exists, use -F to re-run.')
        
        start, step, stop = args.begin_frame, args.traj_step, args.end_frame
        sum_frames = (len(u.trajectory) if stop is None else stop) - start

        # Get apply atom group (atoms to be written to output)
        if args.apply_chains is not None:
            apply_ag = filter_atoms_by_chains(u.atoms, args.apply_chains)
        else:
            apply_ag = u.atoms
        put_log(f'Apply selection: {len(apply_ag)} atoms from chains {args.apply_chains}')
        
        # Save PDB file for the apply selection (first frame)
        for ftype in ['.pdb', '.gro']:
            pdb_path = w_dir / output_path.with_suffix(ftype).name
            apply_ag.write(str(pdb_path))
            put_log(f'Structure saved to {pdb_path}')

            # Prepare output trajectory writer
            writer = mda.Writer(str(output_path), n_atoms=len(apply_ag))
            # Write first frame, since it is the reference frame
            writer.write(apply_ag)
        
        align_ag, _, _, _, rot = self.align(u, args, align_select=args.select)
        Rs = rot.cpu().numpy()
        
        # Get reference coordinates from first frame
        u.trajectory[args.begin_frame]
        ref_coords = align_ag.positions.copy().astype(np.float64)
        ref_com = ref_coords.mean(axis=0)
        
        # Process each frame
        for i, _ in tqdm(enumerate(u.trajectory[start:stop:step]), total=sum_frames//step,
                        desc='Rotating trajectory', leave=False):
            apply_coords = apply_ag.positions.copy()
            apply_ag.positions = (apply_coords - apply_coords.mean(axis=0, keepdims=True)) @ Rs[i] + ref_com
            writer.write(apply_ag)

        writer.close()
        put_log(f'Aligned trajectory saved to {output_path}')

    def main_process(self):
        self.top_paths, self.traj_paths = self.check_top_traj()
        self.tasks = self.find_tasks()
        print(f'find {len(self.tasks)} tasks.')
        # process each task
        bar = tqdm(total=len(self.tasks), desc='Aligning trajectories')
        for top_path, traj_path in self.tasks:
            wdir = os.path.dirname(top_path)
            wdir_repr = os.path.relpath(wdir, self.args.batch_dir)
            bar.set_description(f"{wdir_repr}: {os.path.basename(top_path)} and {os.path.basename(traj_path)}")
            print(f'Loading {traj_path}...')
            u = mda.Universe(top_path, traj_path)
            print(f'Loading {traj_path} done. {len(u.trajectory)} frames.')
            wdir = Path(wdir).resolve()
            self.analysis(u, wdir, self.args)
            bar.update(1)
        bar.close()


_str2func = {
    'trjalign': trjalign,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description='CUDA-accelerated tools for MDAnalysis.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == '__main__':
    # dev code

    main()
