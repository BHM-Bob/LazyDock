#!/usr/bin/env python3
"""
Date: 2025-06-13
Description: Amber (pmemd) 分子动力学模拟运行脚本
参考: amber_workflow.md

STEPS (黄金标准5步法):
1. Minimization Step 1: Restraining the solute (约束溶质)
2. Minimization Step 2: No restraints (全系统无约束)
3. Heating: 0K -> 300K with restraints (NVT加热，约束主链)
4. NPT Equilibration with weak restraints (弱约束NPT平衡)
5. NPT Equilibration without restraints (无约束NPT平衡)
6. Production MD (生产模拟)
"""
import argparse
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from tqdm import tqdm

from lazydock.scripts._script_utils_ import Command, clean_path, process_batch_dir_lst


class AmberRunner(Command):
    """Amber模拟运行基类"""
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.amber_home = os.environ.get('AMBERHOME', '')
        if not self.amber_home:
            put_err('AMBERHOME environment variable not set, exit.', _exit=True)
        
        # 检测可用的pmemd
        self.pmemd_bin = self._detect_pmemd()
        put_log(f'using pmemd binary: {self.pmemd_bin}')
    
    def _detect_pmemd(self) -> str:
        """检测可用的pmemd版本"""
        # 优先检测GPU版本
        for binary in ['pmemd.cuda_SPFP', 'pmemd.cuda_DPFP', 'pmemd.cuda', 'pmemd']:
            result = subprocess.run(['which', binary], capture_output=True, text=True)
            if result.returncode == 0:
                return binary
        return 'pmemd.cuda'  # 默认
    
    def run_pmemd(self, mdin_file: Path, prmtop: Path, inpcrd: Path,
                  output_prefix: str, ref_crd: Optional[Path] = None,
                  working_dir: Optional[Path] = None) -> bool:
        """运行pmemd模拟"""
        if working_dir is None:
            working_dir = inpcrd.parent
        
        # 构建命令
        cmd_parts = [
            self.pmemd_bin, '-O',
            '-i', str(mdin_file),
            '-o', str(working_dir / f'{output_prefix}.out'),
            '-p', str(prmtop),
            '-c', str(inpcrd),
            '-r', str(working_dir / f'{output_prefix}.rst'),
        ]
        
        # 添加参考坐标（用于约束）
        if ref_crd is not None:
            cmd_parts.extend(['-ref', str(ref_crd)])
        
        # 添加轨迹输出（非最小化）
        mdin_content = mdin_file.read_text()
        if 'imin = 0' in mdin_content or 'imin=0' in mdin_content:
            cmd_parts.extend(['-x', str(working_dir / f'{output_prefix}.nc')])
            cmd_parts.extend(['-inf', str(working_dir / f'{output_prefix}.mdinfo')])
        
        cmd = ' '.join(cmd_parts)
        put_log(f'running: {cmd}')
        
        result = subprocess.run(cmd, shell=True, capture_output=True, 
                               text=True, cwd=str(working_dir))
        
        if result.returncode != 0:
            put_err(f'pmemd failed: {result.stderr}')
            return False
        
        return True
    
    def write_mdin(self, working_dir: Path, filename: str, 
                   content: str) -> Path:
        """写入mdin文件"""
        mdin_path = working_dir / filename
        mdin_path.write_text(content)
        return mdin_path


class simple_protein(AmberRunner):
    """
    run simple protein Amber simulation (黄金标准5步法)
    
    STEPS:
    1. Minimization with restraints (约束溶质)
    2. Minimization without restraints (全系统无约束)
    3. Heating: 0K -> 300K with restraints (NVT加热)
    4. NPT equilibration with weak restraints (弱约束NPT)
    5. NPT equilibration without restraints (无约束NPT)
    6. Production MD
    """
    
    HELP = """
    run simple protein Amber simulation (黄金标准5步法)
    
    STEPS:
    1. Minimization with restraints (约束溶质)
    2. Minimization without restraints (全系统无约束)
    3. Heating: 0K -> 300K with restraints (NVT加热)
    4. NPT equilibration with weak restraints (弱约束NPT)
    5. NPT equilibration without restraints (无约束NPT)
    6. Production MD
    """
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--prmtop-name', type=str, default='*.prmtop',
                          help='prmtop file name pattern. Default is %(default)s.')
        args.add_argument('-c', '--inpcrd-name', type=str, default='*.inpcrd',
                          help='inpcrd file name pattern. Default is %(default)s.')
        
        # Minimization parameters
        args.add_argument('--min1-maxcyc', type=int, default=5000,
                          help='max cycles for minimization step 1. Default is %(default)s.')
        args.add_argument('--min1-ncyc', type=int, default=2500,
                          help='steepest descent steps for minimization step 1. Default is %(default)s.')
        args.add_argument('--min1-restraint-wt', type=float, default=10.0,
                          help='restraint weight for minimization step 1. Default is %(default)s.')
        args.add_argument('--min1-restraint-mask', type=str, default='!@H=',
                          help='restraint mask for minimization step 1. Default is %(default)s.')
        args.add_argument('--min2-maxcyc', type=int, default=5000,
                          help='max cycles for minimization step 2. Default is %(default)s.')
        args.add_argument('--min2-ncyc', type=int, default=2500,
                          help='steepest descent steps for minimization step 2. Default is %(default)s.')
        
        # Heating parameters
        args.add_argument('--heat-nstlim', type=int, default=50000,
                          help='number of steps for heating (100ps @ 2fs). Default is %(default)s.')
        args.add_argument('--heat-temp0', type=float, default=300.0,
                          help='target temperature for heating. Default is %(default)s.')
        args.add_argument('--heat-restraint-wt', type=float, default=5.0,
                          help='restraint weight for heating. Default is %(default)s.')
        args.add_argument('--heat-restraint-mask', type=str, default='@CA,C,N',
                          help='restraint mask for heating. Default is %(default)s.')
        
        # NPT parameters
        args.add_argument('--npt1-nstlim', type=int, default=50000,
                          help='number of steps for NPT step 1 (100ps @ 2fs). Default is %(default)s.')
        args.add_argument('--npt1-restraint-wt', type=float, default=1.0,
                          help='restraint weight for NPT step 1. Default is %(default)s.')
        args.add_argument('--npt1-restraint-mask', type=str, default='@CA,C,N',
                          help='restraint mask for NPT step 1. Default is %(default)s.')
        args.add_argument('--npt2-nstlim', type=int, default=100000,
                          help='number of steps for NPT step 2 (200ps @ 2fs). Default is %(default)s.')
        
        # Production MD parameters
        args.add_argument('--md-nstlim', type=int, default=50000000,
                          help='number of steps for production MD (100ns @ 2fs). Default is %(default)s.')
        args.add_argument('--md-temp0', type=float, default=300.0,
                          help='target temperature for production MD. Default is %(default)s.')
        args.add_argument('--md-ntpr', type=int, default=5000,
                          help='energy output frequency for production MD. Default is %(default)s.')
        args.add_argument('--md-ntwx', type=int, default=5000,
                          help='coordinate output frequency for production MD. Default is %(default)s.')
        args.add_argument('--md-ntwr', type=int, default=50000,
                          help='restart file output frequency for production MD. Default is %(default)s.')
        
        # Common parameters
        args.add_argument('--cut', type=float, default=10.0,
                          help='nonbonded cutoff distance (Angstrom). Default is %(default)s.')
        args.add_argument('--dt', type=float, default=0.002,
                          help='time step (ps). Default is %(default)s.')
        args.add_argument('--barostat', type=int, default=2,
                          choices=[1, 2],
                          help='barostat: 1=Berendsen, 2=Monte Carlo (GPU推荐). Default is %(default)s.')
        args.add_argument('--start-step', type=int, default=1,
                          help='start from which step (1-6). Default is %(default)s.')
        args.add_argument('--end-step', type=int, default=6,
                          help='end at which step (1-6). Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
    
    def get_input_files(self, working_dir: Path) -> Optional[Tuple[Path, Path]]:
        """获取输入文件"""
        # 查找prmtop文件
        if '*' in self.args.prmtop_name:
            prmtop_paths = list(working_dir.glob(self.args.prmtop_name))
        else:
            prmtop_paths = [working_dir / self.args.prmtop_name]
        
        if not prmtop_paths or not prmtop_paths[0].exists():
            put_err(f'no prmtop file found in {working_dir}')
            return None
        
        prmtop = prmtop_paths[0]
        
        # 查找inpcrd文件
        if '*' in self.args.inpcrd_name:
            inpcrd_paths = list(working_dir.glob(self.args.inpcrd_name))
        else:
            inpcrd_paths = [working_dir / self.args.inpcrd_name]
        
        if not inpcrd_paths or not inpcrd_paths[0].exists():
            put_err(f'no inpcrd file found in {working_dir}')
            return None
        
        inpcrd = inpcrd_paths[0]
        
        return prmtop, inpcrd
    
    def create_min1_mdin(self) -> str:
        """创建能量最小化第1步的mdin文件"""
        return f"""Minimization Step 1: Restraining the solute
 &cntrl
   imin = 1,
   maxcyc = {self.args.min1_maxcyc},
   ncyc = {self.args.min1_ncyc},
   ntb = 1,
   cut = {self.args.cut},
   ntr = 1,
   restraint_wt = {self.args.min1_restraint_wt},
   restraintmask = '{self.args.min1_restraint_mask}',
   ntpr = 100,
 /
"""
    
    def create_min2_mdin(self) -> str:
        """创建能量最小化第2步的mdin文件"""
        return f"""Minimization Step 2: No restraints
 &cntrl
   imin = 1,
   maxcyc = {self.args.min2_maxcyc},
   ncyc = {self.args.min2_ncyc},
   ntb = 1,
   cut = {self.args.cut},
   ntr = 0,
   ntpr = 100,
 /
"""
    
    def create_heat_mdin(self) -> str:
        """创建加热步骤的mdin文件"""
        return f"""Heating from 0K to {self.args.heat_temp0}K with restraints
 &cntrl
   imin = 0,
   nstlim = {self.args.heat_nstlim},
   dt = {self.args.dt},
   ntx = 1,
   irest = 0,
   ntb = 1,
   cut = {self.args.cut},
   ntc = 2,
   ntf = 2,
   tempi = 0.0,
   temp0 = {self.args.heat_temp0},
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 1,
   restraint_wt = {self.args.heat_restraint_wt},
   restraintmask = '{self.args.heat_restraint_mask}',
   ntpr = 500,
   ntwx = 500,
   ntwr = 5000,
   ig = -1,
 /
"""
    
    def create_npt1_mdin(self) -> str:
        """创建NPT第1步的mdin文件"""
        return f"""NPT Equilibration with weak restraints
 &cntrl
   imin = 0,
   nstlim = {self.args.npt1_nstlim},
   dt = {self.args.dt},
   ntx = 5,
   irest = 1,
   ntb = 2,
   ntp = 1,
   pres0 = 1.0,
   taup = 2.0,
   barostat = {self.args.barostat},
   cut = {self.args.cut},
   ntc = 2,
   ntf = 2,
   temp0 = 300.0,
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 1,
   restraint_wt = {self.args.npt1_restraint_wt},
   restraintmask = '{self.args.npt1_restraint_mask}',
   ntpr = 500,
   ntwx = 500,
   ntwr = 5000,
   ioutfm = 1,
 /
"""
    
    def create_npt2_mdin(self) -> str:
        """创建NPT第2步的mdin文件"""
        return f"""NPT Equilibration without restraints
 &cntrl
   imin = 0,
   nstlim = {self.args.npt2_nstlim},
   dt = {self.args.dt},
   ntx = 5,
   irest = 1,
   ntb = 2,
   ntp = 1,
   pres0 = 1.0,
   taup = 2.0,
   barostat = {self.args.barostat},
   cut = {self.args.cut},
   ntc = 2,
   ntf = 2,
   temp0 = 300.0,
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 0,
   ntpr = 1000,
   ntwx = 1000,
   ntwr = 10000,
   ioutfm = 1,
 /
"""
    
    def create_md_mdin(self) -> str:
        """创建生产模拟的mdin文件"""
        return f"""Production MD
 &cntrl
   imin = 0,
   nstlim = {self.args.md_nstlim},
   dt = {self.args.dt},
   ntx = 5,
   irest = 1,
   ntb = 2,
   ntp = 1,
   pres0 = 1.0,
   taup = 2.0,
   barostat = {self.args.barostat},
   cut = {self.args.cut},
   ntc = 2,
   ntf = 2,
   temp0 = {self.args.md_temp0},
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 0,
   ntpr = {self.args.md_ntpr},
   ntwx = {self.args.md_ntwx},
   ntwr = {self.args.md_ntwr},
   ioutfm = 1,
   ig = -1,
 /
"""
    
    def run_step(self, step: int, prmtop: Path, inpcrd: Path, 
                 working_dir: Path) -> bool:
        """运行指定步骤"""
        put_log(f'running step {step}...')
        
        if step == 1:
            # Minimization step 1: with restraints
            mdin = self.write_mdin(working_dir, 'min1.in', self.create_min1_mdin())
            rst_file = working_dir / 'min1.rst'
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 1.')
                return True
            return self.run_pmemd(mdin, prmtop, inpcrd, 'min1', 
                                 ref_crd=inpcrd, working_dir=working_dir)
        
        elif step == 2:
            # Minimization step 2: without restraints
            mdin = self.write_mdin(working_dir, 'min2.in', self.create_min2_mdin())
            rst_file = working_dir / 'min2.rst'
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 2.')
                return True
            inpcrd_step2 = working_dir / 'min1.rst'
            if not inpcrd_step2.exists():
                put_err(f'min1.rst not found, cannot run step 2')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step2, 'min2',
                                 working_dir=working_dir)
        
        elif step == 3:
            # Heating: 0K -> 300K with restraints
            mdin = self.write_mdin(working_dir, 'heat.in', self.create_heat_mdin())
            rst_file = working_dir / 'heat.rst'
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 3.')
                return True
            inpcrd_step3 = working_dir / 'min2.rst'
            if not inpcrd_step3.exists():
                put_err(f'min2.rst not found, cannot run step 3')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step3, 'heat',
                                 ref_crd=inpcrd_step3, working_dir=working_dir)
        
        elif step == 4:
            # NPT step 1: with weak restraints
            mdin = self.write_mdin(working_dir, 'npt1.in', self.create_npt1_mdin())
            rst_file = working_dir / 'npt1.rst'
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 4.')
                return True
            inpcrd_step4 = working_dir / 'heat.rst'
            if not inpcrd_step4.exists():
                put_err(f'heat.rst not found, cannot run step 4')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step4, 'npt1',
                                 ref_crd=inpcrd_step4, working_dir=working_dir)
        
        elif step == 5:
            # NPT step 2: without restraints
            mdin = self.write_mdin(working_dir, 'npt2.in', self.create_npt2_mdin())
            rst_file = working_dir / 'npt2.rst'
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 5.')
                return True
            inpcrd_step5 = working_dir / 'npt1.rst'
            if not inpcrd_step5.exists():
                put_err(f'npt1.rst not found, cannot run step 5')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step5, 'npt2',
                                 working_dir=working_dir)
        
        elif step == 6:
            # Production MD
            mdin = self.write_mdin(working_dir, 'md.in', self.create_md_mdin())
            rst_file = working_dir / 'md.rst'
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 6.')
                return True
            inpcrd_step6 = working_dir / 'npt2.rst'
            if not inpcrd_step6.exists():
                put_err(f'npt2.rst not found, cannot run step 6')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step6, 'md',
                                 working_dir=working_dir)
        
        return False
    
    def main_process(self):
        # 处理每个batch目录
        for batch_dir in self.args.batch_dir:
            put_log(f'processing batch directory: {batch_dir}')
            
            # 获取所有子目录
            subdirs = [d for d in batch_dir.iterdir() if d.is_dir()]
            
            if not subdirs:
                # 如果没有子目录，尝试直接处理batch_dir
                subdirs = [batch_dir]
            
            for working_dir in tqdm(subdirs, total=len(subdirs)):
                put_log(f'processing: {working_dir}')
                
                # 获取输入文件
                input_files = self.get_input_files(working_dir)
                if input_files is None:
                    put_err(f'failed to get input files in {working_dir}, skip.')
                    continue
                
                prmtop, inpcrd = input_files
                
                # 运行指定步骤
                for step in range(self.args.start_step, self.args.end_step + 1):
                    if not self.run_step(step, prmtop, inpcrd, working_dir):
                        put_err(f'step {step} failed in {working_dir}, stop.')
                        break
                else:
                    put_log(f'successfully completed all steps in {working_dir}')


class continue_md(AmberRunner):
    """
    continue a previous Amber MD simulation
    
    从已有的.rst文件继续生产模拟
    """
    
    HELP = """
    continue a previous Amber MD simulation
    
    从已有的.rst文件继续生产模拟
    """
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, default='.',
                          help='working directory. Default is %(default)s.')
        args.add_argument('-p', '--prmtop', type=str, required=True,
                          help='prmtop file path.')
        args.add_argument('-c', '--inpcrd', type=str, required=True,
                          help='input coordinate/restart file path.')
        args.add_argument('-o', '--output-prefix', type=str, default='md_continue',
                          help='output file prefix. Default is %(default)s.')
        args.add_argument('--nstlim', type=int, default=50000000,
                          help='number of steps for continuation. Default is %(default)s.')
        args.add_argument('--temp0', type=float, default=300.0,
                          help='target temperature. Default is %(default)s.')
        args.add_argument('--cut', type=float, default=10.0,
                          help='nonbonded cutoff distance. Default is %(default)s.')
        args.add_argument('--dt', type=float, default=0.002,
                          help='time step (ps). Default is %(default)s.')
        args.add_argument('--barostat', type=int, default=2,
                          choices=[1, 2],
                          help='barostat: 1=Berendsen, 2=Monte Carlo. Default is %(default)s.')
        args.add_argument('--ntpr', type=int, default=5000,
                          help='energy output frequency. Default is %(default)s.')
        args.add_argument('--ntwx', type=int, default=5000,
                          help='coordinate output frequency. Default is %(default)s.')
        args.add_argument('--ntwr', type=int, default=50000,
                          help='restart file output frequency. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        self.args.prmtop = clean_path(self.args.prmtop)
        self.args.inpcrd = clean_path(self.args.inpcrd)
    
    def create_mdin(self) -> str:
        """创建继续模拟的mdin文件"""
        return f"""Production MD (continued)
 &cntrl
   imin = 0,
   nstlim = {self.args.nstlim},
   dt = {self.args.dt},
   ntx = 5,
   irest = 1,
   ntb = 2,
   ntp = 1,
   pres0 = 1.0,
   taup = 2.0,
   barostat = {self.args.barostat},
   cut = {self.args.cut},
   ntc = 2,
   ntf = 2,
   temp0 = {self.args.temp0},
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 0,
   ntpr = {self.args.ntpr},
   ntwx = {self.args.ntwx},
   ntwr = {self.args.ntwr},
   ioutfm = 1,
   ig = -1,
 /
"""
    
    def main_process(self):
        working_dir = self.args.dir
        
        # 检查输入文件
        prmtop = Path(self.args.prmtop)
        inpcrd = Path(self.args.inpcrd)
        
        if not prmtop.exists():
            put_err(f'prmtop file not found: {prmtop}', _exit=True)
        if not inpcrd.exists():
            put_err(f'inpcrd file not found: {inpcrd}', _exit=True)
        
        # 创建mdin文件
        mdin = self.write_mdin(working_dir, f'{self.args.output_prefix}.in', 
                               self.create_mdin())
        
        # 运行模拟
        if self.run_pmemd(mdin, prmtop, inpcrd, self.args.output_prefix,
                         working_dir=working_dir):
            put_log(f'successfully continued MD: {self.args.output_prefix}')
        else:
            put_err(f'failed to continue MD: {self.args.output_prefix}')


def main():
    parser = argparse.ArgumentParser(description='Amber MD simulation runner')
    subparsers = parser.add_subparsers(dest='command', help='commands')
    
    # simple_protein command
    simple_protein_parser = subparsers.add_parser('simple_protein', 
                                                   help='run simple protein simulation')
    simple_protein.make_args(simple_protein_parser)
    
    # continue_md command
    continue_md_parser = subparsers.add_parser('continue_md',
                                               help='continue a previous MD simulation')
    continue_md.make_args(continue_md_parser)
    
    args = parser.parse_args()
    
    if args.command == 'simple_protein':
        cmd = simple_protein(args)
    elif args.command == 'continue_md':
        cmd = continue_md(args)
    else:
        parser.print_help()
        return
    
    cmd.excute()


if __name__ == '__main__':
    main()
