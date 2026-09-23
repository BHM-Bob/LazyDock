'''
Date: 2025-02-20 10:00:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-20 10:00:00
Description: Prepare ligand files for docking
'''

import argparse
import os
import traceback
from collections import defaultdict
from pathlib import Path
from typing import List

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web_utils.task import TaskPool
from pymol import cmd
from tqdm import tqdm

from lazydock.scripts._script_utils_ import (Command, excute_command,
                                             process_batch_dir_lst)


class smiles2pdb(Command):
    SOURCE_REPR = 'SMILES'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        # Check if RDKit is available
        try:
            from rdkit import Chem # type: ignore
            from rdkit.Chem import AllChem # type: ignore
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
        args.add_argument('--force', action='store_true',
                          help='force overwrite the output file, default is %(default)s.')

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
        cif_paths = get_paths_with_extension(self.args.batch_dir, ['.cif'],
                                             name_substr=self.args.main_name, sort='natsort')
        put_log(f'get {len(cif_paths)} task(s)')
        # process each
        for cif_path in tqdm(cif_paths, total=len(cif_paths)):
            cif_path = Path(cif_path).resolve()
            if self.args.force and self.new_name_fn(cif_path, _args=self.args).exists():
                put_log(f'skip: {self.new_name_fn(cif_path, _args=self.args)}')
                continue
            try:
                cmd.reinitialize()
                cmd.set('connect_mode', 4)
                cmd.set('pdb_conect_all', 'on')
                cmd.load(str(cif_path))
                cmd.save(str(self.new_name_fn(cif_path, suffix=self.args.suffix, new_name=self.args.new_name)))
            except:
                put_err(f'error with {cif_path}')
                traceback.print_exc()


def _fix_cp_relax_one(fix_path: Path, cyclic_chains: List[str], disulfide_chains: List[str],
                      constraints: str, max_iter: int, tolerance: float,
                      platform: str, gpu_index: int):
    """对修复后的 PDB 做 OpenMM 弛豫, 保护环化键/二硫键"""
    from openmm import app as openmm_app

    from lazydock.opmm.relax import ForceFieldMinimizer
    if constraints == 'hbond':
        constraints_obj = openmm_app.HBonds
    elif constraints == 'all':
        constraints_obj = openmm_app.AllBonds
    else:
        constraints_obj = None
    relaxer = ForceFieldMinimizer(
        stiffness=10,
        max_iterations=max_iter,
        tolerance=tolerance,
        platform=platform,
        constraints=constraints_obj,
        cyclic_chains=cyclic_chains,
        disulfide_chains=disulfide_chains,
        device_index=gpu_index,
    )
    result_pdb, ret = relaxer(str(fix_path), None, return_info=True)
    with open(fix_path, 'w') as f:
        f.write(result_pdb)
    put_log(f'relaxed: {fix_path.name} (efinal={ret["efinal"]:.1f} kJ/mol)')


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
        args.add_argument('-c', '--case', type=str, default='backbone', choices=['backbone', 'disulfide', 'auto'],
                          help='fix case, default is %(default)s: '
                               'backbone(仅头尾碳骨架环) / disulfide(仅二硫键环) / auto(前置检测分流).')
        args.add_argument('--only-chains', type=str, nargs='+', default=None,
                          help='only fix given chains (cyclic peptide chains), e.g. --only-chains P, '
                               'default is %(default)s (auto-detect).')
        args.add_argument('--min-cn', type=float, default=2.0,
                          help='max C-N distance for treating as cyclic peptide closure site (Angstrom), '
                               'default is %(default)s.')
        args.add_argument('--ss-max-dist', type=float, default=3.5,
                          help='max SG-SG distance for treating as disulfide bond (Angstrom), '
                               'default is %(default)s.')
        args.add_argument('--remove-outlier-atoms', type=str, nargs='+', default=None,
                          help='remove atoms too far from their residue centroid (chains to process), '
                               'e.g. --remove-outlier-atoms A P, default is %(default)s (disabled).')
        args.add_argument('--outlier-dist', type=float, default=10.0,
                          help='distance threshold (Angstrom) from residue centroid for an atom to be '
                               'considered an outlier and removed, only used when --remove-outlier-atoms '
                               'is set, default is %(default)s.')
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
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='number of workers for parallel relaxing, default is %(default)s.')
        args.add_argument('--gpus', type=int, nargs='+', default=[0],
                          help='GPU indices to use for parallel relaxing, default is %(default)s.')
        args.add_argument('--n-task-per-gpu', type=int, default=1,
                          help='max number of concurrent tasks on one GPU, default is %(default)s.')
        args.add_argument('-F', '--force', action='store_true',
                          help='force overwrite output files, default is %(default)s.')
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
    
    @staticmethod
    def _has_cyclic_mark(sN, sC, sOXT, serial2conect):
        """判断链是否带环化标记: 头尾 N-C 间已有 CONECT 键, 或链尾残基有 OXT(线性C端残留)"""
        if sOXT:
            return True
        return (sC in serial2conect.get(sN, [])) or (sN in serial2conect.get(sC, []))

    def _detect_cyclic_sites(self, atoms, serial2conect, only_chains=None):
        """检测候选环化位点: 链的 C端残基 C(=O) 与 N端残基 N 距离 < min_cn.
        OXT 非必需(上游输出可能已无 OXT 且已带头尾 CONECT), 但记录是否存在供后续删除.
        only_chains: None 检测所有链; 有值时仅检测对应链, 链不存在时报告并跳过."""
        from collections import defaultdict

        # 按链收集残基
        chain_residues = defaultdict(dict)  # chain -> {resi: {name: serial}}
        for serial, (chain, resi, name, *_coords) in atoms.items():
            chain_residues[chain].setdefault(resi, {})[name] = serial
        
        if only_chains:
            for c in only_chains:
                if c not in chain_residues:
                    put_log(f'chain {c} not found in this complex, skip.')
        
        sites = []
        for chain, residues in sorted(chain_residues.items()):
            if only_chains and chain not in only_chains:
                continue
            resis = sorted(residues)
            if len(resis) < 2:
                continue
            first, last = resis[0], resis[-1]
            r_first, r_last = residues[first], residues[last]
            if not (r_first.get('N') and r_last.get('C')):
                continue
            sN, sC = r_first['N'], r_last['C']
            sOXT = r_last.get('OXT', None)
            xN = atoms[sN][3:6]; xC = atoms[sC][3:6]
            d = ((xN[0]-xC[0])**2 + (xN[1]-xC[1])**2 + (xN[2]-xC[2])**2) ** 0.5
            if d < self.args.min_cn:
                sites.append({
                    'chain': chain,
                    'first_resi': first, 'last_resi': last,
                    'sN': sN, 'sC': sC, 'sOXT': sOXT,
                    'd_cn': d,
                })
            elif self._has_cyclic_mark(sN, sC, sOXT, serial2conect):
                # 有环化标记(头尾 CONECT/OXT)但空间未闭合: 构象问题而非拓扑问题, 报告并跳过
                put_log(f'chain {chain}: cyclic mark found but C-N dist={d:.1f}A >= min_cn, '
                        f'peptide is NOT spatially closed, skip (cannot fix topology on an open conformation)')
        return sites

    def _detect_disulfide_sites(self, atoms, serial2conect, only_chains=None):
        """检测候选二硫键位点: 链内恰 2 个 CYS, 且 SG-SG 有 CONECT 或距离 < ss_max_dist.
        返回 sites 列表(dict: chain/res1/res2/sS1/sS2/d_ss).
        only_chains: None 检测所有链; 有值时仅检测对应链, 链不存在时报告并跳过."""
        chain_residues = defaultdict(dict)  # chain -> {resi: {name: serial}}
        for serial, (chain, resi, name, *_coords) in atoms.items():
            chain_residues[chain].setdefault(resi, {})[name] = serial

        if only_chains:
            for c in only_chains:
                if c not in chain_residues:
                    put_log(f'chain {c} not found in this complex, skip.')

        sites = []
        for chain, residues in sorted(chain_residues.items()):
            if only_chains and chain not in only_chains:
                continue
            # 恰 2 个 CYS(含 SG 的残基)
            cys_resis = [r for r, d in residues.items() if 'SG' in d]
            if len(cys_resis) != 2:
                continue
            r1, r2 = cys_resis
            sS1, sS2 = residues[r1].get('SG'), residues[r2].get('SG')
            if not (sS1 and sS2):
                continue
            x1 = atoms[sS1][3:6]; x2 = atoms[sS2][3:6]
            d_ss = ((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2) ** 0.5
            has_bond = (sS2 in serial2conect.get(sS1, [])) or (sS1 in serial2conect.get(sS2, []))
            if has_bond or d_ss < self.args.ss_max_dist:
                sites.append({
                    'chain': chain,
                    'res1': r1, 'res2': r2,
                    'sS1': sS1, 'sS2': sS2,
                    'd_ss': d_ss,
                    'has_bond': has_bond,
                })
        return sites
    
    def _fix_backbone(self, pdb_path: Path):
        atoms, atom_lines, serial2conect = self._parse_pdb(str(pdb_path))
        sites = self._detect_cyclic_sites(atoms, serial2conect, self.args.only_chains)

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
            if s['sOXT']:
                put_log(f"{pdb_path.name}: fix chain {s['chain']}: add peptide bond "
                        f"C(res{s['last_resi']},serial{s['sC']}) - N(res{s['first_resi']},serial{s['sN']}), "
                        f"C-N dist={s['d_cn']:.3f}A, removed OXT(serial{s['sOXT']})")
            else:
                put_log(f"{pdb_path.name}: fix chain {s['chain']}: add peptide bond "
                        f"C(res{s['last_resi']},serial{s['sC']}) - N(res{s['first_resi']},serial{s['sN']}), "
                        f"C-N dist={s['d_cn']:.3f}A, no OXT")
        
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

    def _fix_disulfide(self, pdb_path: Path):
        """修复二硫键环: 保留全部原子(不删 OXT), 仅确保 SG-SG 间存在 CONECT 键.
        返回 (lines, sites)."""
        atoms, atom_lines, serial2conect = self._parse_pdb(str(pdb_path))
        sites = self._detect_disulfide_sites(atoms, serial2conect, self.args.only_chains)

        if not sites:
            put_log(f'{pdb_path.name}: no disulfide site found, skip.')
            return None
        self._last_fixed_chains = [s['chain'] for s in sites]

        # 1. 更新 CONECT: 确保 SG-SG 键存在(不删任何原子行)
        new_conect = {serial: set(nbrs) for serial, nbrs in serial2conect.items()}
        removed_nc = []
        for s in sites:
            if not s['has_bond']:
                new_conect.setdefault(s['sS1'], set()).add(s['sS2'])
                new_conect.setdefault(s['sS2'], set()).add(s['sS1'])
                put_log(f"{pdb_path.name}: fix chain {s['chain']}: add S-S bond "
                        f"SG(res{s['res1']},serial{s['sS1']}) - SG(res{s['res2']},serial{s['sS2']}), "
                        f"S-S dist={s['d_ss']:.3f}A, OXT kept")
            else:
                put_log(f"{pdb_path.name}: chain {s['chain']}: S-S bond exists "
                        f"(res{s['res1']}-res{s['res2']}, d={s['d_ss']:.3f}A), OXT kept")
            # 防御: 二硫键环的首尾 N 端氨基与 C 端羧基之间不应有任何键连接.
            # 若 CONECT 中有错误的首尾 N-C 键(上游距离推断引入), 将其移除
            chain_res = {ser: v for ser, v in atoms.items() if v[0] == s['chain']}
            resis = sorted({v[1] for v in chain_res.values()})
            fN_ser = next((ser for ser, v in chain_res.items()
                           if v[1] == resis[0] and v[2] == 'N'), None)
            lC_ser = next((ser for ser, v in chain_res.items()
                           if v[1] == resis[-1] and v[2] == 'C'), None)
            if fN_ser and lC_ser:
                if lC_ser in new_conect.get(fN_ser, set()):
                    new_conect[fN_ser].discard(lC_ser)
                    removed_nc.append((s['chain'], fN_ser, lC_ser))
                if fN_ser in new_conect.get(lC_ser, set()):
                    new_conect[lC_ser].discard(fN_ser)
                    removed_nc.append((s['chain'], lC_ser, fN_ser))
        for chain, s1, s2 in set(removed_nc):
            put_log(f'{pdb_path.name}: chain {chain}: removed spurious head-to-tail N-C bond '
                    f'(serial {s1}-{s2}) for disulfide ring')

        # 2. 写输出: 保留全部 ATOM 行 + 更新 CONECT
        lines = list(atom_lines)  # 二硫键 case 不删原子(保留 OXT 等)
        conect_out = []
        for serial in sorted(new_conect):
            nbrs = sorted(new_conect[serial])
            if not nbrs:
                continue
            for i in range(0, len(nbrs), 4):
                conect_out.append(f'CONECT{serial:5d}' + ''.join(f'{n:5d}' for n in nbrs[i:i+4]) + '\n')
        lines.extend(conect_out)
        lines.append('END\n')
        return lines

    def _remove_outlier_atoms_from_lines(self, lines: List[str], pdb_path: Path) -> List[str]:
        """在所有 case 修复逻辑之后、relax 之前执行(可选):
        对 --remove-outlier-atoms 指定链, 删除与所属残基质心距离 > --outlier-dist 的离群原子,
        并同步清理其 CONECT 记录. 逐链汇总报告删除的原子(按原子名统计).
        后续 relax 阶段(PDBFixer)会自动补齐缺失原子, 因此这里只做记录级清理.
        """
        if not self.args.remove_outlier_atoms:
            return lines
        repair_chains = set(self.args.remove_outlier_atoms)
        atoms = {}  # serial -> (chain, resi, name, px, py, pz)
        serial2conect = defaultdict(list)
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                serial = int(line[6:11])
                name = line[12:16].strip()
                chain = line[21]
                resi = int(line[22:26])
                px, py, pz = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atoms[serial] = (chain, resi, name, px, py, pz)
            elif line.startswith('CONECT'):
                nums = [int(line[i:i+5]) for i in range(6, len(line.rstrip('\n')), 5)]
                if nums:
                    serial2conect.setdefault(nums[0], []).extend(nums[1:])
        if not atoms:
            return lines

        # 逐链按残基计算质心, 检测离群原子
        removed = defaultdict(lambda: defaultdict(int))  # chain -> atom name -> count
        outlier_serials = set()
        for chain in repair_chains:
            chain_atoms = {ser: v for ser, v in atoms.items() if v[0] == chain}
            resis = sorted({v[1] for v in chain_atoms.values()})
            for resi in resis:
                ser_list = [ser for ser, v in chain_atoms.items() if v[1] == resi]
                n_atoms = len(ser_list)
                if n_atoms == 0:
                    continue
                cx = sum(atoms[ser][3] for ser in ser_list) / n_atoms
                cy = sum(atoms[ser][4] for ser in ser_list) / n_atoms
                cz = sum(atoms[ser][5] for ser in ser_list) / n_atoms
                for ser in ser_list:
                    v = atoms[ser]
                    d = ((v[3]-cx)**2 + (v[4]-cy)**2 + (v[5]-cz)**2) ** 0.5
                    if d > self.args.outlier_dist:
                        outlier_serials.add(ser)
                        removed[chain][v[2]] += 1

        if not outlier_serials:
            return lines

        # 重建输出: 删除离群原子行与旧 CONECT 行, 丢弃离群 serial 的键, 再追加新 CONECT
        out_lines = []
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                if int(line[6:11]) not in outlier_serials:
                    out_lines.append(line)
            elif not line.startswith(('CONECT', 'END')):
                out_lines.append(line)
        new_conect = {}
        for serial, nbrs in serial2conect.items():
            if serial in outlier_serials:
                continue
            nbrs = [n for n in nbrs if n not in outlier_serials]
            if nbrs:
                new_conect[serial] = sorted(set(nbrs))
        for serial in sorted(new_conect):
            nbrs = sorted(new_conect[serial])
            for i in range(0, len(nbrs), 4):
                out_lines.append(f'CONECT{serial:5d}' + ''.join(f'{n:5d}' for n in nbrs[i:i+4]) + '\n')
        out_lines.append('END\n')

        # 报告: 单 PDB 分链汇总
        report_lines = []
        n_total = 0
        for chain in sorted(removed):
            items = ', '.join(f'{name}: {cnt}' for name, cnt in sorted(removed[chain].items()))
            n_chain = sum(removed[chain].values())
            n_total += n_chain
            report_lines.append(f'  chain {chain}: removed {n_chain} atom(s) ({items})')
        put_log(f'{pdb_path.name}: removed {n_total} outlier atom(s) farther than '
                f'{self.args.outlier_dist}A from residue centroid:\n' + '\n'.join(report_lines))
        return out_lines

    def main_process(self):
        pdb_paths = [Path(p).resolve() for p in get_paths_with_extension(
            self.args.batch_dir, ['.pdb'], name_substr=self.args.name, sort='natsort')]
        put_log(f'get {len(pdb_paths)} pdb file(s) in {self.args.batch_dir}')
        
        # GPU 槽位: 把 --gpus 展开成 [g0 x n_task_per_gpu, g1 x n_task_per_gpu, ...],
        # 任务按 index 轮询分配, 保证单卡可承载多个并发任务
        gpu_slots = []
        if self.args.gpus:
            for g in self.args.gpus:
                gpu_slots.extend([g] * max(1, self.args.n_task_per_gpu))
        pool = TaskPool('process', self.args.n_workers, report_error=True).start()
                
        for task_i, pdb_path in tqdm(enumerate(pdb_paths), total=len(pdb_paths)):
            out_path = pdb_path.with_name(pdb_path.stem + self.args.suffix + '.pdb')
            if out_path.exists() and not self.args.force:
                continue
            
            try:
                out_lines = self._fix_one(pdb_path)
            except Exception:
                traceback.print_exc()
                continue
            if out_lines is None:
                continue
            # 所有 case 修复逻辑之后、relax 之前: 删除离群原子(可选)
            out_lines = self._remove_outlier_atoms_from_lines(out_lines, pdb_path)
            with open(out_path, 'w') as f:
                f.writelines(out_lines)
            self.printf(f'fixed: {pdb_path.name} -> {out_path.name}')
            if self.args.relax and self._last_fixed_chains:
                if self._last_fixed_case == 'disulfide':
                    pool.add_task(None, _fix_cp_relax_one, out_path, [], self._last_fixed_chains,
                                  constraints=self.args.constraints, max_iter=self.args.max_iter,
                                  tolerance=self.args.tolerance, platform=self.args.platform,
                                  gpu_index=gpu_slots[task_i % len(gpu_slots)])
                else:
                    pool.add_task(None, _fix_cp_relax_one, out_path, self._last_fixed_chains, [],
                                  constraints=self.args.constraints, max_iter=self.args.max_iter,
                                  tolerance=self.args.tolerance, platform=self.args.platform,
                                  gpu_index=gpu_slots[task_i % len(gpu_slots)])
            pool.wait_till_free()
        pool.wait_till_all_done()
        pool.close(1)

    def _fix_one(self, pdb_path: Path):
        """按 case 分发单个文件修复: backbone / disulfide / auto."""
        if self.args.case == 'disulfide':
            self._last_fixed_case = 'disulfide'
            return self._fix_disulfide(pdb_path)
        if self.args.case == 'backbone':
            self._last_fixed_case = 'backbone'
            return self._fix_backbone(pdb_path)
        # auto: 前置检测分流
        atoms, _, serial2conect = self._parse_pdb(str(pdb_path))
        disu_sites = self._detect_disulfide_sites(atoms, serial2conect, self.args.only_chains)
        if disu_sites:
            put_log(f'{pdb_path.name}: auto -> disulfide')
            self._last_fixed_case = 'disulfide'
            return self._fix_disulfide(pdb_path)
        cyc_sites = self._detect_cyclic_sites(atoms, serial2conect, self.args.only_chains)
        if cyc_sites:
            put_log(f'{pdb_path.name}: auto -> backbone')
            self._last_fixed_case = 'backbone'
            return self._fix_backbone(pdb_path)
        put_log(f'{pdb_path.name}: auto: neither disulfide nor backbone closure detected, skip.')
        return None


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