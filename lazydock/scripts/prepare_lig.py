'''
Date: 2025-02-20 10:00:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-20 10:00:00
Description: Prepare ligand files for docking
'''

import argparse
import os
from pathlib import Path
from typing import List

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web_utils.task import TaskPool
from pymol import cmd
from tqdm import tqdm

from lazydock.scripts._script_utils_ import Command, excute_command, process_batch_dir_lst


class smiles2pdb(Command):
    SOURCE_REPR = 'SMILES'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        # Check if RDKit is available
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            put_err('RDKit is required for SMILES to PDB conversion. Please install it with: pip install rdkit', _exit=True)
        # Store RDKit modules for later use
        self.Chem = Chem
        self.AllChem = AllChem
        self.transfer_method = self.Chem.MolFromSmiles
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-i', '--input', type=str, required=True,
                          help=f'{smiles2pdb.SOURCE_REPR} string to convert to PDB format. Required.')
        args.add_argument('-o', '--output', type=str, required=True,
                          help=f'Output PDB file path. Required.')
        return args
    
    def process_args(self):
        
        # Process output path
        self.args.output = Path(self.args.output).resolve()
        if not self.args.output.parent.exists():
            self.args.output.parent.mkdir(parents=True, exist_ok=True)
    
    def main_process(self):
        try:
            # Convert SMILES to molecule
            mol = self.transfer_method(self.args.input)
            if mol is None:
                put_err(f'Failed to create molecule from {self.SOURCE_REPR}: {self.args.input}', _exit=True)
            
            # Add hydrogens
            mol = self.Chem.AddHs(mol, addResidueInfo=True)
            
            # Generate 3D coordinates
            self.AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # Minimize the structure
            self.AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            
            # Write to PDB file
            pdb_block = self.Chem.MolToPDBBlock(mol)
            with open(self.args.output, 'w') as f:
                f.write(pdb_block)
            
            self.printf(f'Successfully converted {self.SOURCE_REPR} to PDB: {self.args.output}')
        except Exception as e:
            put_err(f'Error during {self.SOURCE_REPR} to PDB conversion: {str(e)}', _exit=True)


class seq2pdb(smiles2pdb):
    SOURCE_REPR = 'Amino acid sequence'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.transfer_method = self.Chem.MolFromFASTA
        
        
class cif2pdb(Command):
    HELP = """"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--main-name', type=str, required=True,
                          help='file in each sub-directory, such as model.cif.')
        args.add_argument('--new-name', type=str, default=None,
                          help='new name of the pdb, such as complex.pdb, default is %(default)s.')
        args.add_argument('--suffix', type=str, default=None,
                          help='suffix of the output pdb, such as _transfer, default is %(default)s.')

    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        if self.args.new_name:
            self.new_name_fn = self.new_name
        elif self.args.suffix:
            self.new_name_fn = self.add_suffix
        else:
            self.new_name_fn = self.only_pdb
        
    @staticmethod
    def only_pdb(cif_path: Path, *args, **kwargs):
        return cif_path.with_suffix('.pdb')
    
    @staticmethod
    def add_suffix(cif_path: Path, suffix: str, *args, **kwargs):
        return cif_path.with_suffix(f'{suffix}.pdb')
    
    @staticmethod
    def new_name(cif_path: Path, new_name: str, *args, **kwargs):
        return cif_path.with_name(f'{new_name}').with_suffix('.pdb')
        
    def main_process(self):
        # get complex paths
        cif_paths = get_paths_with_extension(self.args.batch_dir, [self.args.main_name], name_substr=self.args.main_name)
        put_log(f'get {len(cif_paths)} task(s)')
        # process each
        for cif_path in tqdm(cif_paths, total=len(cif_paths)):
            cif_path = Path(cif_path).resolve()
            cmd.reinitialize()
            cmd.set('connect_mode', 4)
            cmd.set('pdb_conect_all', 'on')
            cmd.load(str(cif_path))
            cmd.save(str(self.new_name_fn(cif_path, suffix=self.args.suffix, new_name=self.args.new_name)))


class fix_cp(Command):
    HELP = """Fix cyclic peptide (头尾碳骨架环肽) topology.

case=backbone: 检测并修复首尾大环肽的拓扑错误.
    AlphaFold3 预测环肽时, 会在环化位点(C端残基羰基C与N端残基氨基N之间)留下线性C端羧基氧 OXT,
    导致 OXT 与环化键上的 N 几乎重叠 (<0.5A). 此工具:
      1. 删除环化位点 C 端残基的 OXT 原子(及其 CONECT 记录)
      2. 将 C端羰基C 与 N端氨基N 之间的 CONECT 肽键连接补上
    (仅编辑拓扑连接, 不做任何坐标移动)

future cases: e.g. disulfide bond, side-chain cyclization."""
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains input pdb files or sub-folders, default is %(default)s.")
        args.add_argument('-n', '--name', type=str, default='',
                          help="input pdb file name substring, default is %(default)s.")
        args.add_argument('-c', '--case', type=str, default='backbone', choices=['backbone'],
                          help='fix case, default is %(default)s, support: backbone(首尾环肽).')
        args.add_argument('--only-chains', type=str, nargs='+', default=None,
                          help='only fix given chains (cyclic peptide chains), e.g. --only-chains P, '
                               'default is %(default)s (auto-detect).')
        args.add_argument('--min-cn', type=float, default=2.0,
                          help='max C-N distance for treating as cyclic peptide closure site (Angstrom), '
                               'default is %(default)s.')
        args.add_argument('--suffix', type=str, default='_fix',
                          help='suffix of the output pdb file name, default is %(default)s.')
        args.add_argument('-relax', '--relax', action='store_true',
                          help='after fixing topology, run OpenMM relaxation with cyclic chain protection '
                               '(sew head-tail peptide bond + CustomBondForce), default is %(default)s.')
        args.add_argument('-it', '--max-iter', type=int, default=1000,
                          help='max iterations for relaxation, default is %(default)s.')
        args.add_argument('-t', '--tolerance', type=int, default=10,
                          help='tolerance for relaxation, default is %(default)s.')
        args.add_argument('-p', '--platform', type=str, default='CUDA', choices=['CUDA', 'CPU'],
                          help='platform for relaxation, default is %(default)s.')
        args.add_argument('--constraints', type=str, default='hbond',
                          choices=['hbond', 'all', 'none'],
                          help='constraints type for relaxation, default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
    
    @staticmethod
    def _parse_pdb(pdb_path: str):
        """解析 PDB, 返回 atoms(按serial), 原始原子行, 以及 CONECT 邻接表"""
        atoms = {}  # serial -> (chain, resi, name, px, py, pz)
        atom_lines = []  # 原始行
        serial2conect = {}  # serial -> list of neighbor serials
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    serial = int(line[6:11])
                    name = line[12:16].strip()
                    chain = line[21]
                    resi = int(line[22:26])
                    px, py, pz = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    atoms[serial] = (chain, resi, name, px, py, pz)
                    atom_lines.append(line)
                elif line.startswith('CONECT'):
                    nums = [int(line[i:i+5]) for i in range(6, len(line.rstrip('\n')), 5)]
                    if nums:
                        serial2conect.setdefault(nums[0], []).extend(nums[1:])
        return atoms, atom_lines, serial2conect
    
    def _detect_cyclic_sites(self, atoms, serial2conect):
        """检测候选环化位点: 链的 C端残基 C(=O) 与 N端残基 N 距离 < min_cn 且有 OXT"""
        from collections import defaultdict
        # 按链收集残基
        chain_residues = defaultdict(dict)  # chain -> {resi: {name: serial}}
        for serial, (chain, resi, name, *_coords) in atoms.items():
            chain_residues[chain].setdefault(resi, {})[name] = serial
        
        sites = []
        for chain, residues in sorted(chain_residues.items()):
            resis = sorted(residues)
            if len(resis) < 2:
                continue
            first, last = resis[0], resis[-1]
            r_first, r_last = residues[first], residues[last]
            if not (r_first.get('N') and r_last.get('C') and r_last.get('OXT')):
                continue
            sN, sC, sOXT = r_first['N'], r_last['C'], r_last['OXT']
            xN = atoms[sN][3:6]; xC = atoms[sC][3:6]
            d = ((xN[0]-xC[0])**2 + (xN[1]-xC[1])**2 + (xN[2]-xC[2])**2) ** 0.5
            if d < self.args.min_cn:
                sites.append({
                    'chain': chain,
                    'first_resi': first, 'last_resi': last,
                    'sN': sN, 'sC': sC, 'sOXT': sOXT,
                    'd_cn': d,
                })
        return sites
    
    def _fix_backbone(self, pdb_path: Path):
        atoms, atom_lines, serial2conect = self._parse_pdb(str(pdb_path))
        sites = self._detect_cyclic_sites(atoms, serial2conect)
        
        # 链过滤
        if self.args.only_chains:
            sites = [s for s in sites if s['chain'] in self.args.only_chains]
        
        if not sites:
            put_log(f'{pdb_path.name}: no cyclic peptide closure site found, skip.')
            return None
        self._last_fixed_chains = [s['chain'] for s in sites]
        
        oxt_serials = {s['sOXT'] for s in sites}
        # 1. 构造新的 CONECT: 删除 OXT 的键, 在 C端C 与 N端N 之间建立肽键
        new_conect = {}  # serial -> set(neighbors)
        for serial, nbrs in serial2conect.items():
            if serial in oxt_serials:
                continue  # 删除 OXT 的 CONECT 行(原子行也会删)
            nbrs = [n for n in nbrs if n not in oxt_serials]  # 从邻居中去掉 OXT
            new_conect[serial] = set(nbrs)
        
        # 2. 补肽键: C端C <-> N端N
        for s in sites:
            new_conect.setdefault(s['sC'], set()).add(s['sN'])
            new_conect.setdefault(s['sN'], set()).add(s['sC'])
            put_log(f"{pdb_path.name}: fix chain {s['chain']}: add peptide bond "
                    f"C(res{s['last_resi']},serial{s['sC']}) - N(res{s['first_resi']},serial{s['sN']}), "
                    f"C-N dist={s['d_cn']:.3f}A, removed OXT(serial{s['sOXT']})")
        
        # 3. 写输出
        lines = []
        for line in atom_lines:
            serial = int(line[6:11])
            if serial not in oxt_serials:
                lines.append(line)
        # 重新生成 CONECT(按 serial 排序): PDB CONECT 每行放 [自身 serial + 最多4个邻居], 邻居多于4个拆多行
        conect_out = []
        for serial in sorted(new_conect):
            nbrs = sorted(new_conect[serial])
            if not nbrs:
                continue
            for i in range(0, len(nbrs), 4):
                row = f'CONECT{serial:5d}' + ''.join(f'{n:5d}' for n in nbrs[i:i+4]) + '\n'
                conect_out.append(row)
        lines.extend(conect_out)
        lines.append('END\n')
        return lines
    
    def _relax_one(self, fix_path: Path, cyclic_chains: List[str]):
        """对修复后的 PDB 做 OpenMM 弛豫, 保护环化键(缝环+CustomBondForce)"""
        from lazydock.opmm.relax import ForceFieldMinimizer
        from openmm import app as openmm_app
        if self.args.constraints == 'hbond':
            constraints_obj = openmm_app.HBonds
        elif self.args.constraints == 'all':
            constraints_obj = openmm_app.AllBonds
        else:
            constraints_obj = None
        relaxer = ForceFieldMinimizer(
            stiffness=10,
            max_iterations=self.args.max_iter,
            tolerance=self.args.tolerance,
            platform=self.args.platform,
            constraints=constraints_obj,
            cyclic_chains=cyclic_chains,
        )
        result_pdb, ret = relaxer(str(fix_path), None, return_info=True)
        with open(fix_path, 'w') as f:
            f.write(result_pdb)
        self.printf(f'relaxed: {fix_path.name} (efinal={ret["efinal"]:.1f} kJ/mol)')

    def main_process(self):
        pdb_paths = [Path(p).resolve() for p in get_paths_with_extension(
            self.args.batch_dir, ['.pdb'], name_substr=self.args.name)]
        put_log(f'get {len(pdb_paths)} pdb file(s) in {self.args.batch_dir}')
        for pdb_path in tqdm(pdb_paths, total=len(pdb_paths)):
            out_lines = self._fix_backbone(pdb_path)
            if out_lines is None:
                continue
            # batch 模式下每个文件输出到同目录 + suffix, 避免固定 output 互相覆盖
            out_path = pdb_path.with_name(pdb_path.stem + self.args.suffix + '.pdb')
            with open(out_path, 'w') as f:
                f.writelines(out_lines)
            self.printf(f'fixed: {pdb_path.name} -> {out_path.name}')
            if self.args.relax:
                # 获取本次修复的环化链
                cyclic_chains = self._last_fixed_chains
                if cyclic_chains:
                    self._relax_one(out_path, cyclic_chains)


_str2func = {
    'smiles2pdb': smiles2pdb,
    'seq2pdb': seq2pdb,
    'cif2pdb': cif2pdb,
    'fix-cp': fix_cp,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description='Prepare ligand files for docking')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))
    
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main()