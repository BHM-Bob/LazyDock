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
import re
import subprocess
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from tqdm import tqdm

from lazydock.scripts._script_utils_ import (Command, check_file_num_paried,
                                             clean_path, process_batch_dir_lst)


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
    
    def parse_mdinfo(self, mdinfo_path: Path) -> Optional[Dict]:
        """解析pmemd的.mdinfo文件，获取进度信息"""
        if not mdinfo_path.exists():
            return None
        
        try:
            content = mdinfo_path.read_text()
            info = {}
            
            # 解析总步数 (NSTEP)
            nstep_match = re.search(r'NSTEP\s*=\s*(\d+)', content)
            if nstep_match:
                info['current_step'] = int(nstep_match.group(1))
            
            # 解析时间 (TIME)
            time_match = re.search(r'TIME\s*\(PS\)\s*=\s*([\d.]+)', content)
            if time_match:
                info['time_ps'] = float(time_match.group(1))
            
            # 解析温度 (TEMP)
            temp_match = re.search(r'TEMP\s*\(K\)\s*=\s*([\d.]+)', content)
            if temp_match:
                info['temperature'] = float(temp_match.group(1))
            
            # 解析压力 (PRESS)
            press_match = re.search(r'PRESS\s*=\s*([-\d.]+)', content)
            if press_match:
                info['pressure'] = float(press_match.group(1))
            
            # 解析总能量 (ETOT)
            etot_match = re.search(r'ETOT\s*=\s*([-\d.]+)', content)
            if etot_match:
                info['total_energy'] = float(etot_match.group(1))
            
            # 解析性能信息 (ns/day)
            perf_match = re.search(r'ns/day\s*=\s*([\d.]+)', content)
            if perf_match:
                info['ns_per_day'] = float(perf_match.group(1))
            
            return info
        except Exception:
            return None
    
    def monitor_progress(self, mdinfo_path: Path, total_steps: int, 
                         interval: int = 30, desc: str = "MD"):
        """监控模拟进度，使用tqdm显示进度条"""
        # 创建tqdm进度条
        pbar = tqdm(total=total_steps, desc=desc, unit="steps", 
                    bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}] {postfix}')
        
        last_step = 0
        
        while True:
            time.sleep(interval)
            
            # 检查mdinfo文件是否存在
            if not mdinfo_path.exists():
                # 检查模拟是否已经结束（通过检查rst文件）
                rst_file = mdinfo_path.with_suffix('.rst')
                if rst_file.exists():
                    # 模拟完成，更新到100%
                    pbar.n = total_steps
                    pbar.refresh()
                    pbar.close()
                    break
                continue
            
            info = self.parse_mdinfo(mdinfo_path)
            if info is None:
                continue
            
            current_step = info.get('current_step', 0)
            if current_step == 0 or current_step <= last_step:
                continue
            
            # 更新进度条
            pbar.n = current_step
            
            # 构建postfix信息
            postfix_parts = []
            if 'temperature' in info:
                postfix_parts.append(f"T={info['temperature']:.1f}K")
            # 压力为0.0时不显示（NVT步骤）
            if 'pressure' in info and abs(info['pressure']) > 0.01:
                postfix_parts.append(f"P={info['pressure']:.1f}bar")
            if 'ns_per_day' in info:
                postfix_parts.append(f"{info['ns_per_day']:.1f}ns/day")
            
            if postfix_parts:
                pbar.set_postfix_str(" | ".join(postfix_parts))
            
            pbar.refresh()
            last_step = current_step
            
            # 如果已经完成，退出监控
            if current_step >= total_steps:
                pbar.close()
                break
    
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
        is_minimization = 'imin = 1' in mdin_content or 'imin=1' in mdin_content
        
        if not is_minimization:
            cmd_parts.extend(['-x', str(working_dir / f'{output_prefix}.nc')])
            cmd_parts.extend(['-inf', str(working_dir / f'{output_prefix}.mdinfo')])
        
        cmd = ' '.join(cmd_parts)
        put_log(f'running: {cmd}')
        
        # 对于非最小化步骤，启动进度监控线程
        monitor_thread = None
        if not is_minimization:
            # 从mdin文件中提取nstlim
            nstlim_match = re.search(r'nstlim\s*=\s*(\d+)', mdin_content)
            if nstlim_match:
                total_steps = int(nstlim_match.group(1))
                mdinfo_path = working_dir / f'{output_prefix}.mdinfo'
                # 从output_prefix提取步骤描述 (e.g., "step03_heat" -> "heat")
                step_desc = output_prefix.split('_')[-1] if '_' in output_prefix else output_prefix
                monitor_thread = threading.Thread(
                    target=self.monitor_progress,
                    args=(mdinfo_path, total_steps, 30, step_desc),
                    daemon=True
                )
                monitor_thread.start()
        
        result = subprocess.run(cmd, shell=True, capture_output=True, 
                               text=True, cwd=str(working_dir))
        
        # 等待监控线程结束
        if monitor_thread is not None:
            monitor_thread.join(timeout=5)
        
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


class simple_md(AmberRunner):
    """
    run simple Amber MD simulation (黄金标准5步法)
    
    适用于蛋白质、配体或复合物系统。
    Amber使用Langevin动力学(ntt=3)，不需要像GROMACS那样定义温度耦合组。
    
    STEPS:
    1. Minimization with restraints (约束溶质)
    2. Minimization without restraints (全系统无约束)
    3. Heating: 0K -> 300K with restraints (NVT加热)
    4. NPT equilibration with weak restraints (弱约束NPT)
    5. NPT equilibration without restraints (无约束NPT)
    6. Production MD
    """
    
    HELP = """
    run simple Amber MD simulation (黄金标准5步法)
    
    适用于蛋白质、配体或复合物系统。
    Amber使用Langevin动力学(ntt=3)，不需要像GROMACS那样定义温度耦合组。
    
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
    
    def check_prmtop_inpcrd(self, batch_dir: Path) -> Tuple[List[str], List[str]]:
        """检查prmtop和inpcrd文件是否配对，类似于ana_gmx的check_top_traj"""
        prmtop_paths = get_paths_with_extension(
            batch_dir,
            [os.path.split(self.args.prmtop_name)[-1]],
            name_substr=self.args.prmtop_name
        )
        inpcrd_paths = get_paths_with_extension(
            batch_dir,
            [os.path.split(self.args.inpcrd_name)[-1]],
            name_substr=self.args.inpcrd_name
        )
        
        invalid_roots = check_file_num_paried(prmtop_paths, inpcrd_paths)
        if invalid_roots:
            put_err(f"The number of prmtop and inpcrd files is not equal, please check the input files.\n"
                   f"invalid roots: {invalid_roots}", _exit=True)
        
        return prmtop_paths, inpcrd_paths
    
    def find_tasks(self) -> List[Tuple[Path, Path, Path]]:
        """查找所有任务，返回(prmtop_path, inpcrd_path, working_dir)列表"""
        tasks = []
        for batch_dir in self.args.batch_dir:
            prmtop_paths, inpcrd_paths = self.check_prmtop_inpcrd(batch_dir)
            for prmtop_path, inpcrd_path in zip(prmtop_paths, inpcrd_paths):
                working_dir = Path(prmtop_path).parent
                tasks.append((Path(prmtop_path), Path(inpcrd_path), working_dir))
        return tasks
    
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
    
    STEP_NAMES = {
        1: 'min1',
        2: 'min2',
        3: 'heat',
        4: 'npt1',
        5: 'npt2',
        6: 'md'
    }
    
    def get_step_filename(self, step: int, suffix: str = 'in') -> str:
        """获取步骤文件名，格式为 stepXX_name.suffix"""
        name = self.STEP_NAMES.get(step, f'step{step}')
        return f'step{step:02d}_{name}.{suffix}'
    
    def get_step_rst_file(self, step: int, working_dir: Path) -> Path:
        """获取步骤的rst文件路径"""
        return working_dir / self.get_step_filename(step, 'rst')
    
    def run_step(self, step: int, prmtop: Path, inpcrd: Path, 
                 working_dir: Path) -> bool:
        """运行指定步骤"""
        put_log(f'running step {step}...')
        
        # 获取步骤文件名
        mdin_file = self.get_step_filename(step, 'in')
        output_prefix = f'step{step:02d}_{self.STEP_NAMES.get(step, f"step{step}")}'
        
        if step == 1:
            # Minimization step 1: with restraints
            mdin = self.write_mdin(working_dir, mdin_file, self.create_min1_mdin())
            rst_file = self.get_step_rst_file(step, working_dir)
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 1.')
                return True
            return self.run_pmemd(mdin, prmtop, inpcrd, output_prefix, 
                                 ref_crd=inpcrd, working_dir=working_dir)
        
        elif step == 2:
            # Minimization step 2: without restraints
            mdin = self.write_mdin(working_dir, mdin_file, self.create_min2_mdin())
            rst_file = self.get_step_rst_file(step, working_dir)
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 2.')
                return True
            inpcrd_step2 = self.get_step_rst_file(1, working_dir)
            if not inpcrd_step2.exists():
                put_err(f'{inpcrd_step2} not found, cannot run step 2')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step2, output_prefix,
                                 working_dir=working_dir)
        
        elif step == 3:
            # Heating: 0K -> 300K with restraints
            mdin = self.write_mdin(working_dir, mdin_file, self.create_heat_mdin())
            rst_file = self.get_step_rst_file(step, working_dir)
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 3.')
                return True
            inpcrd_step3 = self.get_step_rst_file(2, working_dir)
            if not inpcrd_step3.exists():
                put_err(f'{inpcrd_step3} not found, cannot run step 3')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step3, output_prefix,
                                 ref_crd=inpcrd_step3, working_dir=working_dir)
        
        elif step == 4:
            # NPT step 1: with weak restraints
            mdin = self.write_mdin(working_dir, mdin_file, self.create_npt1_mdin())
            rst_file = self.get_step_rst_file(step, working_dir)
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 4.')
                return True
            inpcrd_step4 = self.get_step_rst_file(3, working_dir)
            if not inpcrd_step4.exists():
                put_err(f'{inpcrd_step4} not found, cannot run step 4')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step4, output_prefix,
                                 ref_crd=inpcrd_step4, working_dir=working_dir)
        
        elif step == 5:
            # NPT step 2: without restraints
            mdin = self.write_mdin(working_dir, mdin_file, self.create_npt2_mdin())
            rst_file = self.get_step_rst_file(step, working_dir)
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 5.')
                return True
            inpcrd_step5 = self.get_step_rst_file(4, working_dir)
            if not inpcrd_step5.exists():
                put_err(f'{inpcrd_step5} not found, cannot run step 5')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step5, output_prefix,
                                 working_dir=working_dir)
        
        elif step == 6:
            # Production MD
            mdin = self.write_mdin(working_dir, mdin_file, self.create_md_mdin())
            rst_file = self.get_step_rst_file(step, working_dir)
            if rst_file.exists():
                put_log(f'{rst_file} already exists, skip step 6.')
                return True
            inpcrd_step6 = self.get_step_rst_file(5, working_dir)
            if not inpcrd_step6.exists():
                put_err(f'{inpcrd_step6} not found, cannot run step 6')
                return False
            return self.run_pmemd(mdin, prmtop, inpcrd_step6, output_prefix,
                                 working_dir=working_dir)
        
        return False
    
    def main_process(self):
        """批量处理所有目录中的prmtop/inpcrd文件对"""
        # 查找所有任务
        tasks = self.find_tasks()
        put_log(f'found {len(tasks)} task(s).')
        
        if not tasks:
            put_err('no tasks found, exit.')
            return
        
        # 运行所有任务
        bar = tqdm(total=len(tasks), desc='Running MD')
        for prmtop_path, inpcrd_path, working_dir in tasks:
            wdir_repr = os.path.relpath(working_dir, self.args.batch_dir[0] if self.args.batch_dir else '.')
            bar.set_description(f"{wdir_repr}: {prmtop_path.name} and {inpcrd_path.name}")
            
            # 运行指定步骤
            for step in range(self.args.start_step, self.args.end_step + 1):
                if not self.run_step(step, prmtop_path, inpcrd_path, working_dir):
                    put_err(f'step {step} failed for {prmtop_path}, stop.')
                    break
            else:
                put_log(f'successfully completed all steps for {prmtop_path}')
            
            bar.update(1)


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

    # simple_md command (renamed from simple_protein)
    simple_md_parser = subparsers.add_parser('simple_md',
                                              help='run simple MD simulation (protein/ligand/complex)')
    simple_md.make_args(simple_md_parser)

    # continue_md command
    continue_md_parser = subparsers.add_parser('continue_md',
                                               help='continue a previous MD simulation')
    continue_md.make_args(continue_md_parser)

    args = parser.parse_args()

    if args.command == 'simple_md':
        cmd = simple_md(args)
    elif args.command == 'continue_md':
        cmd = continue_md(args)
    else:
        parser.print_help()
        return

    cmd.excute()


if __name__ == '__main__':
    main()
