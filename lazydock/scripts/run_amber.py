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
    
    def parse_gamd_log(self, gamd_log_path: Path) -> Optional[Dict]:
        """解析GaMD的gamd.log文件最后一行，获取boost能量和力权重信息
        
        gamd.log格式（每行8列）：
        ntwx  total_nstep  Unboosted-Potential-Energy  Unboosted-Dihedral-Energy
        Total-Force-Weight  Dihedral-Force-Weight  Boost-Energy-Potential  Boost-Energy-Dihedral
        """
        if not gamd_log_path.exists():
            return None
        
        try:
            # 读取文件最后几行，找到最后一个数据行
            with open(gamd_log_path, 'rb') as f:
                # 从文件末尾向前读取
                f.seek(0, 2)
                file_size = f.tell()
                if file_size == 0:
                    return None
                # 读取最后2KB足够获取最后一行
                read_size = min(file_size, 2048)
                f.seek(file_size - read_size)
                tail_content = f.read().decode('utf-8', errors='ignore')
            
            lines = tail_content.strip().split('\n')
            # 从后往前找第一个非注释、非空的数据行
            for line in reversed(lines):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        return {
                            'total_nstep': int(parts[1]),
                            'unboosted_potential': float(parts[2]),
                            'unboosted_dihedral': float(parts[3]),
                            'total_force_weight': float(parts[4]),
                            'dihedral_force_weight': float(parts[5]),
                            'boost_energy_potential': float(parts[6]),
                            'boost_energy_dihedral': float(parts[7]),
                        }
                    except (ValueError, IndexError):
                        continue
            return None
        except Exception:
            return None

    def monitor_progress(self, mdinfo_path: Path, total_steps: int, 
                         interval: int = 30, desc: str = "MD",
                         gamd_log_path: Optional[Path] = None):
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
            
            # GaMD专属指标
            if gamd_log_path is not None:
                gamd_info = self.parse_gamd_log(gamd_log_path)
                if gamd_info is not None:
                    bP = gamd_info['boost_energy_potential']
                    bD = gamd_info['boost_energy_dihedral']
                    fwP = gamd_info['total_force_weight']
                    fwD = gamd_info['dihedral_force_weight']
                    gamd_parts = []
                    if bP > 0.01:
                        gamd_parts.append(f"bP={bP:.1f}")
                    if bD > 0.01:
                        gamd_parts.append(f"bD={bD:.1f}")
                    # force weight 低于0.9时显示，提示boost力不完整
                    if fwP < 0.9:
                        gamd_parts.append(f"fwP={fwP:.2f}")
                    if fwD < 0.9:
                        gamd_parts.append(f"fwD={fwD:.2f}")
                    if gamd_parts:
                        postfix_parts.append(" ".join(gamd_parts))
            
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
                # 检测是否为GaMD模拟，如果是则监控gamd.log
                is_gamd = 'igamd' in mdin_content and re.search(r'igamd\s*=\s*[1-9]', mdin_content)
                gamd_log_path = working_dir / 'gamd.log' if is_gamd else None
                monitor_thread = threading.Thread(
                    target=self.monitor_progress,
                    args=(mdinfo_path, total_steps, 30, step_desc),
                    kwargs={'gamd_log_path': gamd_log_path},
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
        mdin_path.write_text(content.rstrip() + "\n")
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


class gamd_md(AmberRunner):
    """
    run GaMD simulation using Amber pmemd

    Supports a two-stage GaMD workflow:
    1. GaMD equilibration to collect boost statistics and ramp up the boost
    2. GaMD production with fixed boost parameters
    """

    HELP = """
    run GaMD simulation using Amber pmemd

    Supports a two-stage GaMD workflow:
    1. GaMD equilibration to collect boost statistics and ramp up the boost
    2. GaMD production with fixed boost parameters

    NOTE: The input file must be a restart containing velocities obtained from
    a standard Amber equilibration workflow such as simple_md NPT equilibration.
    """

    def __init__(self, args, printf=print):
        super().__init__(args, printf)

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type=str, nargs='+', default=['.'],
                          help='dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.')
        args.add_argument('-n', '--prmtop-name', type=str, default='*.prmtop',
                          help='prmtop file name pattern. Default is %(default)s.')
        args.add_argument('-c', '--inpcrd-name', type=str, default='*.inpcrd',
                          help='inpcrd file name pattern. Default is %(default)s.')

        args.add_argument('--eq-nstlim', type=int, default=3000000,
                          help='number of steps for GaMD equilibration. Default is %(default)s.')
        args.add_argument('--prod-nstlim', type=int, default=100000000,
                          help='number of steps for GaMD production MD. Default is %(default)s.')
        args.add_argument('--ntcmd', type=int, default=1000000,
                          help='number of steps for GaMD equilibration statistics collection. Default is %(default)s.')
        args.add_argument('--ieqdtd', type=int, default=None,
                          help='DEPRECATED: use --nteb instead. If set, mapped to nteb as fallback.')
        args.add_argument('--nteb', type=int, default=1000000,
                          help='number of steps for GaMD equilibration boost application. Default is %(default)s.')
        args.add_argument('--ntave', type=int, default=None,
                          help='GaMD averaging interval, must divide ntcmd and nteb. Default is chosen automatically.')
        args.add_argument('--ntcmdprep', type=int, default=None,
                          help='GaMD conventional MD prep steps before statistics collection. Default is ntcmd/4 when not set.')
        args.add_argument('--ntebprep', type=int, default=None,
                          help='GaMD boost prep steps before biased MD. Default is nteb/4 when not set.')
        args.add_argument('--igamd', type=int, default=3,
                          choices=[1, 2, 3, 4],
                          help='GaMD boost type: 1=dihedral, 2=dual, 3=dual hybrid. Default is %(default)s.')
        args.add_argument('--sigma0P', type=float, default=3.0,
                          help='GaMD total potential sigma threshold. Default is %(default)s.')
        args.add_argument('--sigma0D', type=float, default=3.0,
                          help='GaMD dihedral potential sigma threshold. Default is %(default)s.')
        args.add_argument('--ntp', type=int, default=1,
                          choices=[0, 1, 2],
                          help='pressure coupling flag for GaMD runs. Default is %(default)s.')
        args.add_argument('--barostat', type=int, default=2,
                          choices=[1, 2],
                          help='barostat type for NPT GaMD runs. Default is %(default)s.')
        args.add_argument('--cut', type=float, default=10.0,
                          help='nonbonded cutoff distance (Angstrom). Default is %(default)s.')
        args.add_argument('--dt', type=float, default=0.002,
                          help='time step (ps). Default is %(default)s.')
        args.add_argument('--prod-temp0', type=float, default=300.0,
                          help='target temperature for GaMD production. Default is %(default)s.')
        args.add_argument('--prod-ntpr', type=int, default=50000,
                          help='energy output frequency for GaMD production MD. Default is %(default)s.')
        args.add_argument('--prod-ntwx', type=int, default=50000,
                          help='coordinate output frequency for GaMD production MD. Default is %(default)s.')
        args.add_argument('--prod-ntwr', type=int, default=5000000,
                          help='restart file output frequency for GaMD production MD. Default is %(default)s.')
        args.add_argument('--start-phase', type=str, default='equil',
                          choices=['equil', 'prod'],
                          help='start from which GaMD phase: equil or prod. Default is %(default)s.')
        args.add_argument('--end-phase', type=str, default='prod',
                          choices=['equil', 'prod'],
                          help='end at which GaMD phase: equil or prod. Default is %(default)s.')
        return args

    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)

        if self.args.nteb is None:
            if self.args.ieqdtd is not None:
                put_log(f'--ieqdtd is deprecated, use --nteb instead. Mapping ieqdtd={self.args.ieqdtd} to nteb.')
                self.args.nteb = self.args.ieqdtd
            else:
                self.args.nteb = 1000000
        if self.args.ntave is None:
            self.args.ntave = 100000 if self.args.ntcmd >= 100000 and self.args.nteb >= 100000 else 1
        if self.args.ntcmdprep is None:
            self.args.ntcmdprep = max(1, self.args.ntcmd // 4)
        if self.args.ntebprep is None:
            self.args.ntebprep = max(1, self.args.nteb // 4)

        if self.args.ntcmd <= 0:
            put_err('ntcmd must be positive.', _exit=True)
        if self.args.nteb <= 0:
            put_err('nteb must be positive.', _exit=True)
        if self.args.ntave <= 0:
            put_err('ntave must be positive.', _exit=True)
        if self.args.ntcmd % self.args.ntave != 0:
            put_err('ntcmd must be a multiple of ntave.', _exit=True)
        if self.args.nteb % self.args.ntave != 0:
            put_err('nteb must be a multiple of ntave.', _exit=True)
        if not (0 <= self.args.ntcmdprep <= self.args.ntcmd):
            put_err('ntcmdprep must be between 0 and ntcmd.', _exit=True)
        if not (0 <= self.args.ntebprep <= self.args.nteb):
            put_err('ntebprep must be between 0 and nteb.', _exit=True)

    def check_prmtop_inpcrd(self, batch_dir: Path) -> Tuple[List[str], List[str]]:
        """检查prmtop和inpcrd文件是否配对"""
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

    def get_phase_filename(self, phase: str, suffix: str = 'in') -> str:
        return f'gamd_{phase}.{suffix}'

    def get_phase_rst_file(self, phase: str, working_dir: Path) -> Path:
        return working_dir / f'gamd_{phase}.rst'

    def create_gamd_equil_mdin(self) -> str:
        """创建GaMD预平衡mdin文件"""
        return f"""GaMD Equilibration: collect stats and ramp boost
 &cntrl
     imin = 0,
     irest = 1,
     ntx = 5,
   nstlim = {self.args.eq_nstlim},
   dt = {self.args.dt},
   ntc = 2,
   ntf = 2,
   cut = {self.args.cut},
   ntb = 2,
   ntp = {self.args.ntp},
   pres0 = 1.0,
   taup = 2.0,
   temp0 = {self.args.prod_temp0},
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 0,
   iwrap = 1,
   ntpr = 5000,
   ntwx = 5000,
   ntwr = 500000,
   ioutfm = 1,
   ig = -1,
   igamd = {self.args.igamd},
   iE = 1,
   irest_gamd = 0,
   ntcmd = {self.args.ntcmd},
   ntcmdprep = {self.args.ntcmdprep},
   nteb = {self.args.nteb},
   ntebprep = {self.args.ntebprep},
   ntave = {self.args.ntave},
   sigma0P = {self.args.sigma0P},
   sigma0D = {self.args.sigma0D},
 /"""

    def create_gamd_prod_mdin(self) -> str:
        """创建GaMD生产模拟mdin文件"""
        return f"""GaMD Production Run
 &cntrl
   imin = 0,
   irest = 1,
   ntx = 5,
   nstlim = {self.args.prod_nstlim},
   dt = {self.args.dt},
   ntc = 2,
   ntf = 2,
   cut = {self.args.cut},
   ntb = 2,
   ntp = {self.args.ntp},
   pres0 = 1.0,
   taup = 2.0,
   temp0 = {self.args.prod_temp0},
   ntt = 3,
   gamma_ln = 2.0,
   ntr = 0,
   iwrap = 1,
   ntpr = {self.args.prod_ntpr},
   ntwx = {self.args.prod_ntwx},
   ntwr = {self.args.prod_ntwr},
   ioutfm = 1,
   ig = -1,
   igamd = {self.args.igamd},
   iE = 1,
   irest_gamd = 1,
   ntcmd = 0,
   nteb = 0,
   ntave = {self.args.ntave},
   sigma0P = {self.args.sigma0P},
   sigma0D = {self.args.sigma0D},
 /"""

    def _check_velocities_with_pmemd(self, prmtop: Path, inpcrd: Path, working_dir: Path) -> Tuple[bool, str]:
        """Run a lightweight pmemd check to determine if the input restart contains velocities.
        Returns (has_velocities, pmemd_stderr).
        This writes a tiny mdin and runs pmemd; if pmemd errors with "could not find enough velocities" we return False.
        """
        check_mdin = """Check velocities
 &cntrl
   imin = 0,
   irest = 1,
   ntx = 5,
   nstlim = 0,
   ntpr = 0,
   ntwx = 0,
   ntwr = 0,
 /
"""
        mdin_path = working_dir / 'gamd_check_vel.in'
        mdin_path.write_text(check_mdin.rstrip() + "\n")

        cmd_parts = [self.pmemd_bin, '-O', '-i', str(mdin_path), '-o', str(working_dir / 'gamd_check_vel.out'),
                     '-p', str(prmtop), '-c', str(inpcrd)]
        result = subprocess.run(cmd_parts, capture_output=True, text=True, cwd=str(working_dir))
        stderr = (result.stderr or '') + (result.stdout or '')
        if result.returncode == 0:
            return True, stderr
        # detect missing velocities message
        if 'could not find enough velocities' in stderr.lower() or 'i could not find enough velocities' in stderr.lower():
            return False, stderr
        # unknown failure - return False and include stderr for diagnostics
        return False, stderr

    def has_velocities(self, prmtop: Path, inpcrd: Path, working_dir: Path) -> bool:
        """Heuristic check whether the given coordinate/restart file contains velocities.
        Uses a lightweight pmemd invocation to probe for velocities. Returns True if velocities present.
        """
        try:
            ok, stderr = self._check_velocities_with_pmemd(prmtop, inpcrd, working_dir)
            if ok:
                return True
            # if pmemd indicated missing velocities, return False
            if 'could not find enough velocities' in stderr.lower() or 'i could not find enough velocities' in stderr.lower():
                return False
            # otherwise conservatively treat as False and log message
            put_log(f'pmemd velocity-check returned non-zero exit but no explicit missing-velocity message. stderr: {stderr[:200]}')
            return False
        except Exception as e:
            put_log(f'exception while checking velocities: {e}')
            return False

    def check_gamd_equil_quality(self, working_dir: Path) -> bool:
        """Check GaMD equilibration output quality by parsing the final boost parameters.
        Returns True if quality is acceptable, False if boost is too weak (warning only, not fatal).
        """
        equil_out = working_dir / 'gamd_equil.out'
        if not equil_out.exists():
            put_log(f'GaMD equilibration output not found: {equil_out}, skip quality check.')
            return True

        try:
            import re
            content = equil_out.read_text()
            # Find the last "GaMD updated parameters" lines for kP and kD
            kp_matches = re.findall(
                r'GaMD updated parameters: step,VmaxP,VminP,VavgP,sigmaVP,k0P,kP,EthreshP\s*=\s*\d+\s+[\d\-\.\s]+\s+(\d+\.\d+)\s+(\d+\.\d+)',
                content)
            kd_matches = re.findall(
                r'GaMD updated parameters: step,VmaxD,VminD,VavgD,sigmaVD,k0D,kD,EthreshD\s*=\s*\d+\s+[\d\-\.\s]+\s+(\d+\.\d+)\s+(\d+\.\d+)',
                content)

            if kp_matches:
                last_k0P, last_kP = kp_matches[-1]
                last_k0P, last_kP = float(last_k0P), float(last_kP)
                if last_kP < 1e-4:
                    put_err(f'GaMD equilibration produced very weak total-potential boost (kP={last_kP:.6f}, k0P={last_k0P:.4f}). '
                            f'This means GaMD will have almost no acceleration on the total potential energy. '
                            f'Consider: (1) using igamd=2 (dihedral-only boost), or '
                            f'(2) adjusting sigma0P, or '
                            f'(3) increasing ntcmd for better statistics.',
                            _exit=False)
                else:
                    put_log(f'GaMD equilibration boost quality: kP={last_kP:.6f}, k0P={last_k0P:.4f}')

            if kd_matches:
                last_k0D, last_kD = kd_matches[-1]
                last_k0D, last_kD = float(last_k0D), float(last_kD)
                if last_kD < 1e-4:
                    put_err(f'GaMD equilibration produced very weak dihedral boost (kD={last_kD:.6f}, k0D={last_k0D:.4f}). '
                            f'Consider adjusting sigma0D or increasing ntcmd.',
                            _exit=False)
                else:
                    put_log(f'GaMD equilibration boost quality: kD={last_kD:.6f}, k0D={last_k0D:.4f}')

            return True
        except Exception as e:
            put_log(f'exception while checking GaMD equilibration quality: {e}')
            return True


    def run_phase(self, phase: str, prmtop: Path, inpcrd: Path, working_dir: Path) -> bool:
        mdin = self.write_mdin(working_dir, self.get_phase_filename(phase),
                               self.create_gamd_equil_mdin() if phase == 'equil' else self.create_gamd_prod_mdin())
        output_prefix = f'gamd_{phase}'
        rst_file = self.get_phase_rst_file(phase, working_dir)
        if rst_file.exists():
            put_log(f'{rst_file} already exists, skip GaMD {phase} phase.')
            return True

        if phase == 'equil':
            run_inpcrd = inpcrd
            if not self.has_velocities(prmtop, run_inpcrd, working_dir):
                put_err('Input restart file lacks velocities. Provide a restart from simple_md NPT equilibration (for example step05_npt2.rst or later).', _exit=False)
                return False
        else:
            equil_rst = self.get_phase_rst_file('equil', working_dir)
            if not equil_rst.exists():
                put_err(f'GaMD equilibration restart not found: {equil_rst}. Run equil first or place the restart file in the working directory.')
                return False
            run_inpcrd = equil_rst
            if not self.has_velocities(prmtop, run_inpcrd, working_dir):
                put_err('GaMD equilibration restart lacks velocities required for production. Provide a restart with velocities generated by simple_md NPT equilibration or a previous MD run.', _exit=False)
                return False
            # Check equilibration quality before production
            self.check_gamd_equil_quality(working_dir)

        return self.run_pmemd(mdin, prmtop, run_inpcrd, output_prefix, working_dir=working_dir)

    def main_process(self):
        tasks = self.find_tasks()
        put_log(f'found {len(tasks)} GaMD task(s).')

        if not tasks:
            put_err('no tasks found, exit.')
            return

        phase_order = ['equil', 'prod']
        start_index = phase_order.index(self.args.start_phase)
        end_index = phase_order.index(self.args.end_phase)
        if start_index > end_index:
            put_err('start-phase must be earlier than or equal to end-phase.', _exit=True)
        run_phases = phase_order[start_index:end_index + 1]

        bar = tqdm(total=len(tasks), desc='Running GaMD')
        for prmtop_path, inpcrd_path, working_dir in tasks:
            wdir_repr = os.path.relpath(working_dir, self.args.batch_dir[0] if self.args.batch_dir else '.')
            bar.set_description(f"{wdir_repr}: {prmtop_path.name} and {inpcrd_path.name}")
            for phase in run_phases:
                if not self.run_phase(phase, prmtop_path, inpcrd_path, working_dir):
                    put_err(f'GaMD {phase} phase failed for {prmtop_path}, stop.')
                    break
            else:
                put_log(f'successfully completed GaMD phases {run_phases} for {prmtop_path}')
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


def main(sys_args: List[str] = None):
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

    # gamd_md command
    gamd_md_parser = subparsers.add_parser('gamd_md',
                                           help='run GaMD simulation with Amber pmemd')
    gamd_md.make_args(gamd_md_parser)

    args = parser.parse_args(sys_args)

    if args.command == 'simple_md':
        cmd = simple_md(args)
    elif args.command == 'continue_md':
        cmd = continue_md(args)
    elif args.command == 'gamd_md':
        cmd = gamd_md(args)
    else:
        parser.print_help()
        return

    cmd.excute()


if __name__ == '__main__':
    main()
