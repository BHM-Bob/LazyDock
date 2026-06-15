#!/usr/bin/env python3
"""
Date: 2025-06-13
Description: Amber (pmemd) 系统准备脚本
参考: amber_workflow.md

STEPS for protein:
1. clean PDB file by pdb4amber
2. build system using tleap (solvation, add ions)

STEPS for ligand:
1. generate ligand parameters using antechamber
2. check missing parameters using parmchk2
3. build system using tleap

STEPS for complex:
1. prepare receptor and ligand separately
2. combine and build complex system using tleap
"""
import argparse
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from pymol import cmd
from tqdm import tqdm

from lazydock.scripts._script_utils_ import Command, clean_path


class AmberCommand(Command):
    """Amber命令基类"""
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.amber_home = os.environ.get('AMBERHOME', '')
        if not self.amber_home:
            put_err('AMBERHOME environment variable not set, exit.', _exit=True)
        
        # 检测AmberTools路径（用于tleap, antechamber等工具）
        # 优先使用CONDA_PREFIX（如果安装了ambertools）
        self.amber_tools_home = os.environ.get('CONDA_PREFIX', self.amber_home)
        
        # 验证tleap是否存在
        self.tleap_bin = self._find_executable('tleap')
        self.antechamber_bin = self._find_executable('antechamber')
        self.parmchk2_bin = self._find_executable('parmchk2')
        self.pdb4amber_bin = self._find_executable('pdb4amber')
    
    def _find_executable(self, name: str) -> str:
        """查找Amber可执行文件，优先在AmberTools路径中查找"""
        # 首先在AmberTools路径中查找
        if self.amber_tools_home:
            tool_path = Path(self.amber_tools_home) / 'bin' / name
            if tool_path.exists():
                return str(tool_path)
        
        # 然后在AMBERHOME中查找
        if self.amber_home:
            amber_path = Path(self.amber_home) / 'bin' / name
            if amber_path.exists():
                return str(amber_path)
        
        # 最后使用PATH中的命令
        return name
    
    def run_command(self, cmd: str, cwd: Optional[str] = None, 
                    check: bool = True, capture_output: bool = True) -> subprocess.CompletedProcess:
        """运行shell命令"""
        put_log(f'running: {cmd}')
        if capture_output:
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                   text=True, cwd=cwd)
        else:
            # 对于长时间运行的命令（如tleap），直接输出到终端
            result = subprocess.run(cmd, shell=True, cwd=cwd)
        if check and result.returncode != 0:
            put_err(f'command failed: {cmd}\nstderr: {result.stderr if capture_output else ""}')
        return result


class protein(AmberCommand):
    """
    prepare single protein for Amber MDS.
    
    STEPS:
    1. clean PDB file by pdb4amber
    2. build system using tleap (solvation, add ions)
    """
    
    HELP = """
    prepare single protein for Amber MDS.
    
    STEPS:
    1. clean PDB file by pdb4amber
    2. build system using tleap (solvation, add ions)
    """
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, default='.',
                          help='protein directory. Default is %(default)s.')
        args.add_argument('-n', '--protein-name', type=str,
                          help='protein file name in each sub-directory.')
        args.add_argument('--protein-ff', type=str, default='ff19SB',
                          choices=['ff19SB', 'ff14SB'],
                          help='protein force field. Default is %(default)s.')
        args.add_argument('--water-model', type=str, default='opc',
                          choices=['opc', 'tip3p', 'tip4pew'],
                          help='water model. Default is %(default)s.')
        args.add_argument('--box-size', type=float, default=12.0,
                          help='solvation box buffer size (Angstrom). Default is %(default)s.')
        args.add_argument('--box-type', type=str, default='oct',
                          choices=['oct', 'box'],
                          help='box type: oct (octahedron) or box. Default is %(default)s.')
        args.add_argument('--neutralize', action='store_true', default=True,
                          help='add ions to neutralize system. Default is %(default)s.')
        args.add_argument('--ion-conc', type=float, default=0.0,
                          help='ion concentration (M). Default is %(default)s.')
        args.add_argument('--positive-ion', type=str, default='Na+',
                          help='positive ion type. Default is %(default)s.')
        args.add_argument('--negative-ion', type=str, default='Cl-',
                          help='negative ion type. Default is %(default)s.')
        args.add_argument('--reduce', action='store_true', default=True,
                          help='add hydrogens using reduce. Default is %(default)s.')
        args.add_argument('--dry', action='store_true', default=True,
                          help='remove water from input PDB. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
    
    def get_water_source(self) -> Tuple[str, str]:
        """获取水模型source命令和solvate命令"""
        water_map = {
            'opc': ('leaprc.water.opc', 'OPCBOX'),
            'tip3p': ('leaprc.water.tip3p', 'TIP3PBOX'),
            'tip4pew': ('leaprc.water.tip4pew', 'TIP4PEWBOX'),
        }
        return water_map[self.args.water_model]
    
    def clean_pdb(self, input_pdb: Path, output_pdb: Path) -> bool:
        """使用pdb4amber清理PDB文件"""
        if output_pdb.exists():
            put_log(f'{output_pdb} already exists, skip pdb4amber.')
            return True
        
        cmd_parts = [self.pdb4amber_bin, '-i', str(input_pdb), '-o', str(output_pdb)]
        if self.args.dry:
            cmd_parts.append('--dry')
        if self.args.reduce:
            cmd_parts.append('--reduce')
        
        cmd = ' '.join(cmd_parts)
        result = self.run_command(cmd, cwd=str(input_pdb.parent), check=False)
        
        if result.returncode != 0:
            put_err(f'pdb4amber failed for {input_pdb}, try without reduce')
            # 尝试不带reduce
            cmd_parts = [self.pdb4amber_bin, '-i', str(input_pdb), '-o', str(output_pdb)]
            if self.args.dry:
                cmd_parts.append('--dry')
            cmd = ' '.join(cmd_parts)
            result = self.run_command(cmd, cwd=str(input_pdb.parent))
        
        return result.returncode == 0
    
    def build_system(self, protein_path: Path, output_dir: Path) -> bool:
        """使用tleap构建系统"""
        prmtop_path = output_dir / f'{protein_path.stem}.prmtop'
        inpcrd_path = output_dir / f'{protein_path.stem}.inpcrd'
        
        if prmtop_path.exists() and inpcrd_path.exists():
            put_log(f'{prmtop_path} and {inpcrd_path} already exist, skip tleap.')
            return True
        
        water_source, water_box = self.get_water_source()
        
        # 构建tleap输入脚本
        tleap_script = f"""source leaprc.protein.{self.args.protein_ff}
source {water_source}

receptor = loadpdb {protein_path.name}

# Solvation
solvate{self.args.box_type} receptor {water_box} {self.args.box_size}

# Add ions
"""
        if self.args.neutralize:
            tleap_script += f"addions receptor {self.args.positive_ion} 0\n"
            tleap_script += f"addions receptor {self.args.negative_ion} 0\n"
        
        if self.args.ion_conc > 0:
            tleap_script += f"addions receptor {self.args.positive_ion} {self.args.ion_conc}\n"
            tleap_script += f"addions receptor {self.args.negative_ion} {self.args.ion_conc}\n"
        
        tleap_script += f"""
# Save files
savepdb receptor {protein_path.stem}_solv.pdb
saveamberparm receptor {protein_path.stem}.prmtop {protein_path.stem}.inpcrd
quit
"""
        
        # 写入tleap脚本
        tleap_file = output_dir / 'tleap.in'
        tleap_file.write_text(tleap_script)
        
        # 运行tleap
        cmd = f'{self.tleap_bin} -f {tleap_file.name}'
        result = self.run_command(cmd, cwd=str(output_dir), check=False)
        
        if result.returncode != 0:
            put_err(f'tleap failed for {protein_path}')
            return False
        
        return True
    
    def main_process(self):
        # 获取蛋白质路径
        if os.path.isdir(self.args.dir):
            proteins_path = get_paths_with_extension(self.args.dir, [], 
                                                      name_substr=self.args.protein_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.dir}, exit.', _exit=True)
        
        put_log(f'get {len(proteins_path)} protein(s)')
        
        # 处理每个蛋白质
        for protein_path in tqdm(proteins_path, total=len(proteins_path)):
            protein_path = Path(protein_path).resolve()
            working_dir = protein_path.parent
            
            # STEP 1: clean PDB
            cleaned_pdb = working_dir / f'{protein_path.stem}_clean.pdb'
            if not self.clean_pdb(protein_path, cleaned_pdb):
                put_err(f'failed to clean PDB: {protein_path}, skip.')
                continue
            
            # STEP 2: build system with tleap
            if not self.build_system(cleaned_pdb, working_dir):
                put_err(f'failed to build system: {protein_path}, skip.')
                continue
            
            put_log(f'successfully prepared: {protein_path}')


class ligand(AmberCommand):
    """
    prepare single ligand for Amber MDS.
    
    STEPS:
    1. convert ligand to mol2 format (if needed)
    2. generate ligand parameters using antechamber
    3. check missing parameters using parmchk2
    4. build system using tleap
    """
    
    HELP = """
    prepare single ligand for Amber MDS.
    
    STEPS:
    1. convert ligand to mol2 format (if needed)
    2. generate ligand parameters using antechamber
    3. check missing parameters using parmchk2
    4. build system using tleap
    """
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, default='.',
                          help='ligand directory. Default is %(default)s.')
        args.add_argument('-n', '--ligand-name', type=str,
                          help='ligand file name in each sub-directory.')
        args.add_argument('--ligand-ff', type=str, default='gaff2',
                          choices=['gaff', 'gaff2'],
                          help='ligand force field. Default is %(default)s.')
        args.add_argument('--charge-method', type=str, default='bcc',
                          choices=['bcc', 'resp', 'gas'],
                          help='charge calculation method. Default is %(default)s.')
        args.add_argument('--net-charge', type=int, default=0,
                          help='net charge of ligand. Default is %(default)s.')
        args.add_argument('--residue-name', type=str, default='LIG',
                          help='residue name for ligand. Default is %(default)s.')
        args.add_argument('--water-model', type=str, default='opc',
                          choices=['opc', 'tip3p', 'tip4pew'],
                          help='water model. Default is %(default)s.')
        args.add_argument('--box-size', type=float, default=20.0,
                          help='solvation box buffer size (Angstrom). Default is %(default)s.')
        args.add_argument('--box-type', type=str, default='oct',
                          choices=['oct', 'box'],
                          help='box type: oct (octahedron) or box. Default is %(default)s.')
        args.add_argument('--neutralize', action='store_true', default=True,
                          help='add ions to neutralize system. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
    
    def get_water_source(self) -> Tuple[str, str]:
        """获取水模型source命令和solvate命令"""
        water_map = {
            'opc': ('leaprc.water.opc', 'OPCBOX'),
            'tip3p': ('leaprc.water.tip3p', 'TIP3PBOX'),
            'tip4pew': ('leaprc.water.tip4pew', 'TIP4PEWBOX'),
        }
        return water_map[self.args.water_model]
    
    def convert_to_mol2(self, ligand_path: Path, mol2_path: Path) -> bool:
        """转换配体为mol2格式"""
        if mol2_path.exists():
            put_log(f'{mol2_path} already exists, skip conversion.')
            return True
        
        # 使用obabel转换
        cmd = f'obabel -i pdb {ligand_path} -o mol2 -O {mol2_path}'
        result = self.run_command(cmd, check=False)
        
        if result.returncode != 0 or not mol2_path.exists():
            put_err(f'failed to convert {ligand_path} to mol2')
            return False
        
        return True
    
    def run_antechamber(self, mol2_path: Path, prep_path: Path,
                        frcmod_path: Path) -> bool:
        """运行antechamber生成配体参数"""
        if prep_path.exists() and frcmod_path.exists():
            put_log(f'{prep_path} and {frcmod_path} already exist, skip antechamber.')
            return True

        working_dir = mol2_path.parent

        # 运行antechamber
        at_type = self.args.ligand_ff.replace('gaff', 'gaff')
        cmd = (f'{self.antechamber_bin} -i {mol2_path.name} -fi mol2 '
               f'-o {prep_path.name} -fo prepi '
               f'-c {self.args.charge_method} '
               f'-at {at_type} '
               f'-rn {self.args.residue_name} '
               f'-nc {self.args.net_charge}')

        result = self.run_command(cmd, cwd=str(working_dir), check=False)

        if result.returncode != 0:
            put_err(f'antechamber failed for {mol2_path}')
            return False

        # 运行parmchk2检查缺失参数
        cmd = f'{self.parmchk2_bin} -i {prep_path.name} -f prepi -o {frcmod_path.name}'
        result = self.run_command(cmd, cwd=str(working_dir), check=False)

        if result.returncode != 0:
            put_err(f'parmchk2 failed for {mol2_path}')
            return False

        return True
    
    def build_system(self, ligand_path: Path, prep_path: Path,
                     frcmod_path: Path, output_dir: Path) -> bool:
        """使用tleap构建配体系统"""
        prmtop_path = output_dir / f'{ligand_path.stem}.prmtop'
        inpcrd_path = output_dir / f'{ligand_path.stem}.inpcrd'

        if prmtop_path.exists() and inpcrd_path.exists():
            put_log(f'{prmtop_path} and {inpcrd_path} already exist, skip tleap.')
            return True

        water_source, water_box = self.get_water_source()

        # 构建tleap输入脚本
        # loadamberprep 会创建以残基名命名的变量（如 LIG）
        # 使用残基名作为变量名
        residue_name = self.args.residue_name

        tleap_script = f"""source leaprc.{self.args.ligand_ff}
source {water_source}

# 加载配体参数文件，会自动创建名为 {residue_name} 的变量
loadamberprep {prep_path.name}
loadamberparams {frcmod_path.name}

# Solvation - 使用残基名作为变量名
solvate{self.args.box_type} {residue_name} {water_box} {self.args.box_size}

# Add ions
"""
        if self.args.neutralize:
            tleap_script += f"addions {residue_name} Na+ 0\n"
            tleap_script += f"addions {residue_name} Cl- 0\n"

        tleap_script += f"""
# Save files
savepdb {residue_name} {ligand_path.stem}_solv.pdb
saveamberparm {residue_name} {ligand_path.stem}.prmtop {ligand_path.stem}.inpcrd
quit
"""
        
        # 写入tleap脚本
        tleap_file = output_dir / 'tleap_ligand.in'
        tleap_file.write_text(tleap_script)

        # 运行tleap
        cmd = f'{self.tleap_bin} -f {tleap_file.name}'
        result = self.run_command(cmd, cwd=str(output_dir), check=False)

        if result.returncode != 0:
            put_err(f'tleap failed for {ligand_path}')
            return False

        return True
    
    def main_process(self):
        # 获取配体路径
        if os.path.isdir(self.args.dir):
            ligands_path = get_paths_with_extension(self.args.dir, ['.pdb', '.mol2', '.sdf'],
                                                     name_substr=self.args.ligand_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.dir}, exit.', _exit=True)
        
        put_log(f'get {len(ligands_path)} ligand(s)')
        
        # 处理每个配体
        for ligand_path in tqdm(ligands_path, total=len(ligands_path)):
            ligand_path = Path(ligand_path).resolve()
            working_dir = ligand_path.parent
            
            # STEP 1: convert to mol2 if needed
            mol2_path = working_dir / f'{ligand_path.stem}.mol2'
            if ligand_path.suffix.lower() != '.mol2':
                if not self.convert_to_mol2(ligand_path, mol2_path):
                    put_err(f'failed to convert {ligand_path} to mol2, skip.')
                    continue
            else:
                mol2_path = ligand_path
            
            # STEP 2: run antechamber
            prep_path = working_dir / f'{ligand_path.stem}.prep'
            frcmod_path = working_dir / f'{ligand_path.stem}.frcmod'
            if not self.run_antechamber(mol2_path, prep_path, frcmod_path):
                put_err(f'failed to run antechamber for {ligand_path}, skip.')
                continue
            
            # STEP 3: build system with tleap
            if not self.build_system(ligand_path, prep_path, frcmod_path, working_dir):
                put_err(f'failed to build system for {ligand_path}, skip.')
                continue
            
            put_log(f'successfully prepared: {ligand_path}')


class complex(AmberCommand):
    """
    prepare protein-ligand complex for Amber MDS.
    
    - input complex.pdb should have two chains, one for receptor and one for ligand.
    - complex.pdb should already add hydrogens.
    
    STEPS:
    1. extract receptor and ligand from complex.pdb by chain name
    2. prepare receptor using pdb4amber
    3. prepare ligand using antechamber
    4. combine and build complex system using tleap
    """
    
    HELP = """
    prepare protein-ligand complex for Amber MDS.
    
    - input complex.pdb should have two chains, one for receptor and one for ligand.
    - complex.pdb should already add hydrogens.
    
    STEPS:
    1. extract receptor and ligand from complex.pdb by chain name
    2. prepare receptor using pdb4amber
    3. prepare ligand using antechamber
    4. combine and build complex system using tleap
    """
    
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, default='.',
                          help='complex directory. Default is %(default)s.')
        args.add_argument('-n', '--complex-name', type=str,
                          help='complex file name in each sub-directory.')
        args.add_argument('-rc', '--receptor-chain-name', type=str,
                          help='receptor chain name.')
        args.add_argument('-lc', '--ligand-chain-name', type=str,
                          help='ligand chain name.')
        args.add_argument('--protein-ff', type=str, default='ff19SB',
                          choices=['ff19SB', 'ff14SB'],
                          help='protein force field. Default is %(default)s.')
        args.add_argument('--ligand-ff', type=str, default='gaff2',
                          choices=['gaff', 'gaff2'],
                          help='ligand force field. Default is %(default)s.')
        args.add_argument('--water-model', type=str, default='opc',
                          choices=['opc', 'tip3p', 'tip4pew'],
                          help='water model. Default is %(default)s.')
        args.add_argument('--box-size', type=float, default=12.0,
                          help='solvation box buffer size (Angstrom). Default is %(default)s.')
        args.add_argument('--box-type', type=str, default='oct',
                          choices=['oct', 'box'],
                          help='box type: oct (octahedron) or box. Default is %(default)s.')
        args.add_argument('--neutralize', action='store_true', default=True,
                          help='add ions to neutralize system. Default is %(default)s.')
        args.add_argument('--ion-conc', type=float, default=0.0,
                          help='ion concentration (M). Default is %(default)s.')
        args.add_argument('--ligand-charge-method', type=str, default='bcc',
                          choices=['bcc', 'resp', 'gas'],
                          help='charge calculation method for ligand. Default is %(default)s.')
        args.add_argument('--ligand-net-charge', type=int, default=0,
                          help='net charge of ligand. Default is %(default)s.')
        args.add_argument('--ligand-residue-name', type=str, default='LIG',
                          help='residue name for ligand. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
    
    def get_water_source(self) -> Tuple[str, str]:
        """获取水模型source命令和solvate命令"""
        water_map = {
            'opc': ('leaprc.water.opc', 'OPCBOX'),
            'tip3p': ('leaprc.water.tip3p', 'TIP3PBOX'),
            'tip4pew': ('leaprc.water.tip4pew', 'TIP4PEWBOX'),
        }
        return water_map[self.args.water_model]
    
    @staticmethod
    def extract_receptor_ligand(ipath: str, receptor_chain_name: str, 
                                ligand_chain_name: str, opath_r: str, opath_l: str) -> bool:
        """从复合物中提取受体和配体（受体移除氢原子，配体保留）"""
        cmd.reinitialize()
        cmd.load(ipath, 'complex')
        
        # 提取受体（移除氢原子，让 pdb4amber --reduce 重新添加）
        if cmd.select('receptor', f'complex and chain {receptor_chain_name}') == 0:
            put_err(f'receptor chain {receptor_chain_name} has zero atom in {ipath}, skip this complex.')
            return False
        else:
            cmd.select('receptor_noH', 'receptor and not hydrogen')
            cmd.save(opath_r, 'receptor_noH')
        
        # 提取配体（保留氢原子）
        if cmd.select('ligand', f'complex and chain {ligand_chain_name}') == 0:
            put_err(f'ligand chain {ligand_chain_name} has zero atom in {ipath}, skip this complex.')
            return False
        else:
            cmd.save(opath_l, 'ligand')
        
        return True
    
    def prepare_receptor(self, receptor_path: Path, working_dir: Path) -> Optional[Path]:
        """准备受体"""
        cleaned_pdb = working_dir / f'{receptor_path.stem}_clean.pdb'
        if cleaned_pdb.exists():
            put_log(f'{cleaned_pdb} already exists, skip pdb4amber.')
            return cleaned_pdb

        # 使用 --reduce 重新生成氢原子，确保与力场兼容
        cmd_parts = [self.pdb4amber_bin, '-i', str(receptor_path), '-o', str(cleaned_pdb), '--dry', '--reduce']
        cmd_str = ' '.join(cmd_parts)
        result = self.run_command(cmd_str, cwd=str(working_dir), check=False)

        if result.returncode != 0:
            put_err(f'pdb4amber failed for {receptor_path}')
            return None

        return cleaned_pdb
    
    def prepare_ligand(self, ligand_path: Path, working_dir: Path) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
        """准备配体，返回(prep_path, frcmod_path, mol2_path)"""
        # 转换配体为mol2
        mol2_path = working_dir / f'{ligand_path.stem}.mol2'
        if ligand_path.suffix.lower() != '.mol2':
            if not mol2_path.exists():
                cmd_str = f'obabel -i pdb {ligand_path} -o mol2 -O {mol2_path}'
                result = self.run_command(cmd_str, check=False)
                if result.returncode != 0:
                    put_err(f'failed to convert {ligand_path} to mol2')
                    return None, None, None
        else:
            mol2_path = ligand_path

        # 运行antechamber
        prep_path = working_dir / f'{ligand_path.stem}.prep'
        frcmod_path = working_dir / f'{ligand_path.stem}.frcmod'

        if prep_path.exists() and frcmod_path.exists():
            put_log(f'{prep_path} and {frcmod_path} already exist, skip antechamber.')
            return prep_path, frcmod_path, mol2_path

        # 生成GAFF格式的mol2文件（用于loadmol2，包含正确的GAFF原子类型）
        gaff_mol2_path = working_dir / f'{ligand_path.stem}_gaff.mol2'
        cmd_str = (f'{self.antechamber_bin} -i {mol2_path.name} -fi mol2 '
                   f'-o {gaff_mol2_path.name} -fo mol2 '
                   f'-c {self.args.ligand_charge_method} '
                   f'-at {self.args.ligand_ff} '
                   f'-rn {self.args.ligand_residue_name} '
                   f'-nc {self.args.ligand_net_charge}')

        result = self.run_command(cmd_str, cwd=str(working_dir), check=False)
        if result.returncode != 0:
            put_err(f'antechamber failed to generate GAFF mol2 for {ligand_path}')
            return None, None, None

        # 生成prep文件（用于备份兼容）
        cmd_str = (f'{self.antechamber_bin} -i {mol2_path.name} -fi mol2 '
                   f'-o {prep_path.name} -fo prepi '
                   f'-c {self.args.ligand_charge_method} '
                   f'-at {self.args.ligand_ff} '
                   f'-rn {self.args.ligand_residue_name} '
                   f'-nc {self.args.ligand_net_charge}')

        result = self.run_command(cmd_str, cwd=str(working_dir), check=False)
        if result.returncode != 0:
            put_err(f'antechamber failed to generate prep for {ligand_path}')
            return None, None, None

        # 运行parmchk2
        cmd_str = f'{self.parmchk2_bin} -i {prep_path.name} -f prepi -o {frcmod_path.name}'
        result = self.run_command(cmd_str, cwd=str(working_dir), check=False)
        if result.returncode != 0:
            put_err(f'parmchk2 failed for {ligand_path}')
            return None, None, None

        return prep_path, frcmod_path, gaff_mol2_path
    
    def build_complex(self, receptor_path: Path,
                      prep_path: Path, frcmod_path: Path,
                      mol2_path: Path, working_dir: Path) -> bool:
        """构建复合物系统"""
        prmtop_path = working_dir / 'complex.prmtop'
        inpcrd_path = working_dir / 'complex.inpcrd'
        
        if prmtop_path.exists() and inpcrd_path.exists():
            put_log(f'{prmtop_path} and {inpcrd_path} already exist, skip tleap.')
            return True
        
        water_source, water_box = self.get_water_source()
        
        # 构建tleap输入脚本
        # 使用教程方法：先加载所有力场，再加载结构
        tleap_script = f"""# Load force fields (按标准顺序：蛋白 -> 配体 -> 水/离子)
source leaprc.protein.{self.args.protein_ff}
source leaprc.{self.args.ligand_ff}
source {water_source}

# Load receptor
receptor = loadpdb {receptor_path.name}

# Load ligand parameters (使用loadmol2保持坐标，而不是loadamberprep)
loadamberparams {frcmod_path.name}
lig = loadmol2 {mol2_path.name}

# Combine
complex = combine {{receptor lig}}

# Save gas-phase complex (for coordinate verification)
savepdb complex complex_gas.pdb

# Solvation
solvate{self.args.box_type} complex {water_box} {self.args.box_size}

# Add ions
"""
        if self.args.neutralize:
            tleap_script += "addions complex Na+ 0\n"
            tleap_script += "addions complex Cl- 0\n"

        if self.args.ion_conc > 0:
            tleap_script += f"addions complex Na+ {self.args.ion_conc}\n"
            tleap_script += f"addions complex Cl- {self.args.ion_conc}\n"

        tleap_script += f"""
# Save files
savepdb complex complex_solv.pdb
saveamberparm complex complex.prmtop complex.inpcrd
quit
"""
        
        # 写入tleap脚本
        tleap_file = working_dir / 'tleap_complex.in'
        working_dir.mkdir(parents=True, exist_ok=True)
        tleap_file.write_text(tleap_script)
        put_log(f'written tleap script to: {tleap_file}')

        # 运行tleap - 使用完整路径，不捕获输出以支持长时间运行
        cmd_str = f'{self.tleap_bin} -f {tleap_file}'
        self.run_command(cmd_str, cwd=str(working_dir), check=False, capture_output=False)

        # 检查输出文件是否存在（tleap 可能返回非零退出码但成功生成文件）
        if not prmtop_path.exists() or not inpcrd_path.exists():
            put_err(f'tleap failed for complex in {working_dir}')
            return False

        return True
    
    def main_process(self):
        # 获取复合物文件
        if os.path.isdir(self.args.dir):
            complexes_path = get_paths_with_extension(self.args.dir, ['.pdb'],
                                                       name_substr=self.args.complex_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.dir}, exit.', _exit=True)
        
        put_log(f'get {len(complexes_path)} complex(es)')
        
        # 处理每个复合物
        for complex_path in tqdm(complexes_path, total=len(complexes_path)):
            complex_path = Path(complex_path).resolve()
            working_dir = complex_path.parent
            
            # STEP 1: extract receptor and ligand
            receptor_path = working_dir / f'1_{complex_path.stem}_receptor.pdb'
            ligand_path = working_dir / f'1_{complex_path.stem}_ligand.pdb'
            
            if receptor_path.exists() and ligand_path.exists():
                put_log(f'{receptor_path} and {ligand_path} already exist, skip extraction.')
            else:
                if not self.extract_receptor_ligand(str(complex_path),
                                                    self.args.receptor_chain_name,
                                                    self.args.ligand_chain_name,
                                                    str(receptor_path),
                                                    str(ligand_path)):
                    put_err(f'failed to extract receptor/ligand from {complex_path}, skip.')
                    continue
            
            # STEP 2: prepare receptor
            receptor_prepared = self.prepare_receptor(receptor_path, working_dir)
            if receptor_prepared is None:
                put_err(f'failed to prepare receptor for {complex_path}, skip.')
                continue
            
            # STEP 3: prepare ligand
            prep_path, frcmod_path, mol2_path = self.prepare_ligand(ligand_path, working_dir)
            if prep_path is None or frcmod_path is None or mol2_path is None:
                put_err(f'failed to prepare ligand for {complex_path}, skip.')
                continue

            # STEP 4: build complex
            if not self.build_complex(receptor_prepared,
                                      prep_path, frcmod_path, mol2_path, working_dir):
                put_err(f'failed to build complex for {complex_path}, skip.')
                continue
            
            put_log(f'successfully prepared: {complex_path}')


def main():
    parser = argparse.ArgumentParser(description='Amber system preparation')
    subparsers = parser.add_subparsers(dest='command', help='commands')
    
    # protein command
    protein_parser = subparsers.add_parser('protein', help='prepare protein')
    protein.make_args(protein_parser)
    
    # ligand command
    ligand_parser = subparsers.add_parser('ligand', help='prepare ligand')
    ligand.make_args(ligand_parser)
    
    # complex command
    complex_parser = subparsers.add_parser('complex', help='prepare complex')
    complex.make_args(complex_parser)
    
    args = parser.parse_args()
    
    if args.command == 'protein':
        cmd = protein(args)
    elif args.command == 'ligand':
        cmd = ligand(args)
    elif args.command == 'complex':
        cmd = complex(args)
    else:
        parser.print_help()
        return
    
    cmd.excute()


if __name__ == '__main__':
    main()
