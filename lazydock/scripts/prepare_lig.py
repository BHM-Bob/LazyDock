'''
Date: 2025-02-20 10:00:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-20 10:00:00
Description: Prepare ligand files for docking
'''

import argparse
import multiprocessing
import traceback
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, TypedDict

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
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import AllChem  # type: ignore
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
        args.add_argument('--keep-duplicate-conect', action='store_true',
                          help='keep duplicated CONECT neighbours from PyMOL as-is. By default '
                               'cif2pdb rewrites the CONECT block to drop repeated neighbours '
                               '(PyMOL writes a bond twice on both sides of a non-standard '
                               'residue) and to keep a single trailing END record, '
                               'default is %(default)s.')

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

    @staticmethod
    def _dedup_conect(pdb_path: Path) -> int:
        """重写 PDB 的 CONECT 块: 去掉重复邻居、合并同一原子的多行、只保留尾部一个 END.

        PyMOL 在 connect_mode=4 + pdb_conect_all=on 下, 对跨非标准残基的肽键
        会双向各写两次(如 `CONECT 1052 1051 1053 1054 1054` 与
        `CONECT 1054 1052 1055 1052`), 该键其实只有一根, 不是双键;
        另外当输入 PDB 已带 END 时会写出两个 CONECT 块(中间夹一个 END).
        这些都不影响 OpenMM(读入时按 `(a, b)` 对去重), 但会让下游按
        CONECT 计数/解析人困惑, 且严格解析器(如部分 GROMACS 工具链)可能报错.

        返回被去除的重复邻居条目数; 0 表示原文件本身干净.
        """
        lines = pdb_path.read_text().splitlines()
        conect: Dict[int, List[int]] = {}
        order: List[int] = []
        for line in lines:
            if not line.startswith('CONECT'):
                continue
            nums = [int(line[i:i + 5]) for i in range(6, len(line.rstrip()), 5)]
            if not nums:
                continue
            if nums[0] not in conect:
                conect[nums[0]] = []
                order.append(nums[0])
            conect[nums[0]].extend(nums[1:])
        if not conect:
            return 0
        removed = sum(len(v) - len(dict.fromkeys(v)) for v in conect.values())

        body = [ln for ln in lines if not ln.startswith('CONECT')]
        while body and body[-1].strip() in ('END', ''):
            body.pop()
        for serial in order:
            neighbours = list(dict.fromkeys(conect[serial]))
            for i in range(0, len(neighbours), 4):
                body.append(f'CONECT{serial:5d}'
                            + ''.join(f'{n:5d}' for n in neighbours[i:i + 4]))
        body.append('END')
        pdb_path.write_text('\n'.join(body) + '\n')
        return removed

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
            out_path = Path(self.new_name_fn(cif_path, suffix=self.args.suffix,
                                             new_name=self.args.new_name))
            try:
                cmd.reinitialize()
                cmd.set('connect_mode', 4)
                cmd.set('pdb_conect_all', 'on')
                cmd.load(str(cif_path))
                cmd.save(str(out_path))
                if not self.args.keep_duplicate_conect:
                    n_removed = self._dedup_conect(out_path)
                    if n_removed:
                        put_log(f'{out_path.name}: removed {n_removed} duplicated CONECT '
                                f'neighbour(s) written by PyMOL')
            except:
                put_err(f'error with {cif_path}')
                traceback.print_exc()


def _fix_cp_relax_one(fix_path: Path, backbone_chains: List[str], disulfide_chains: List[str],
                      restrain_chains: Optional[List[str]], restrain_backbone: bool,
                      constraints: str, max_iter: int, tolerance: float,
                      platform: str, gpu_index: int,
                      custom_residue_xml: Optional[List[str]] = None,
                      ss_max_dist: float = 0.3) -> None:
    """对修复后的 PDB 做 OpenMM 弛豫, 保护头尾肽键(backbone 闭合)/二硫键.

    ss_max_dist: SG-SG 距离阈值, 单位 nm(与 ForceFieldMinimizer 一致),
        决定 PDBFixer 之后是否由 _sew_disulfide_bonds 重补 S-S 拓扑键.
        调用方需把 fix-cp 的 --ss-max-dist(单位 A)换算后传入."""
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
        backbone_chains=backbone_chains,
        disulfide_chains=disulfide_chains,
        ss_max_dist=ss_max_dist,
        device_index=gpu_index,
        custom_residue_xml=custom_residue_xml,
    )
    result_pdb, ret = relaxer(str(fix_path), None, restrain_chains, restrain_backbone, return_info=True)
    if not result_pdb:
        raise RuntimeError(f'relaxation returned empty structure for {fix_path.name}')
    with open(fix_path, 'w') as f:
        f.write(result_pdb)
    put_log(f'relaxed: {fix_path.name} (efinal={ret["efinal"]:.1f} kJ/mol)')


class BackboneSite(TypedDict):
    """头尾闭合(backbone 环化)候选位点. sOXT 为 None 表示链尾无 OXT(上游已无线性 C 端)."""
    chain: str
    first_resi: int
    last_resi: int
    sN: int
    sC: int
    sOXT: Optional[int]
    d_cn: float


class DisulfideSite(TypedDict):
    """二硫键候选位点. has_bond 表示输入 CONECT 中两 SG 已直接相连."""
    chain: str
    res1: int
    res2: int
    sS1: int
    sS2: int
    d_ss: float
    has_bond: bool


@dataclass
class FixResult:
    """单个 PDB 的修复结果, 由 _fix_one 返回给 main_process 调度 relax.

    backbone_chains:  存在头尾肽键闭合(backbone 环化)的链, 需要 relax 用 CustomBondForce 保护首尾键.
    disulfide_chains: 存在 S-S 闭合的链, 需要 relax 走 CHARMM36 DISU patch / 手动缝接 S-S.
    一条链只会有其一(-c disulfide 模式下 backbone_chains 为空).
    """
    lines: List[str]
    backbone_chains: List[str] = field(default_factory=list)
    disulfide_chains: List[str] = field(default_factory=list)


class fix_cp(Command):
    HELP = """Fix cyclic peptide topology.

两种环化(two cyclization modes):
  - backbone: 头尾闭合, C 端羰基 C 与 N 端氨基 N 成肽键.
  - disulfide: 侧链闭合, 两个 Cys 的 SG-SG 成二硫键.

case=backbone: 检测并修复首尾大环肽的拓扑错误.
    AlphaFold3 预测环肽时, 会在环化位点(C端残基羰基C与N端残基氨基N之间)留下线性C端羧基氧 OXT,
    导致 OXT 与环化键上的 N 几乎重叠 (<0.5A). 此工具:
      1. 删除环化位点 C 端残基的 OXT 原子(及其 CONECT 记录)
      2. 将 C端羰基C 与 N端氨基N 之间的 CONECT 肽键连接补上
    (仅编辑拓扑连接, 不做任何坐标移动)

case=disulfide: 检测并补全链内 SG-SG 的 CONECT 二硫键.

case=auto: 前置检测两类位点后分流; 同一结构可同时含受体二硫键与配体头尾闭合."""
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
        args.add_argument('-rc', '--restrain-chain', type=str, nargs='+', default=None,
                          help='restrain chain in relaxation, default is %(default)s.')
        args.add_argument('-rb', '--restrain-backbone', action='store_true',
                          help='restrain backbone for given restrain chains in relaxation, default is %(default)s.')
        args.add_argument('--min-cn', type=float, default=2.0,
                          help='max C-N distance for treating as a backbone closure site (Angstrom), '
                               'default is %(default)s.')
        args.add_argument('--ss-max-dist', type=float, default=None,
                          help='max SG-SG distance for treating as a disulfide bond (Angstrom). '
                               'default is %(default)s: disabled, i.e. only SG-SG pairs already '
                               'carrying a CONECT record are treated as disulfides, and distance-based '
                               'pairing is left to OpenMM/PDBFixer (0.3 nm threshold). '
                               'Pass an explicit value (e.g. 3.5) only when the structure contains '
                               'non-standard sulfur residues (SG3, e.g. LE1) not recognized by OpenMM, '
                               'or when the S-S distance exceeds OpenMM 0.3 nm but the bond should '
                               'still be restored. When set, the same value is also forwarded to '
                               'relax (converted to nm) so that PDBFixer will not drop such S-S bonds.')
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
                               '(backbone: sew head-tail bond + CustomBondForce; '
                               'disulfide: CHARMM36 DISU patch), default is %(default)s.')
        args.add_argument('-it', '--max-iter', type=int, default=1000,
                          help='max iterations for relaxation, default is %(default)s.')
        args.add_argument('-t', '--tolerance', type=float, default=10,
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
        args.add_argument('--ncaa-xml', type=str, nargs='+', default=None,
                          help='external non-standard amino acid template XML file(s) or dir(s) '
                               'loaded at runtime. Each XML is paired with a same-named YAML '
                               '(e.g. my_nc_aa.xml + my_nc_aa.yml) that describes atom renames. '
                               'default is %(default)s (only builtin templates).')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
    
    @staticmethod
    def _parse_pdb(pdb_path: str) -> Tuple[Dict[int, Tuple[str, int, str, float, float, float]],
                                           List[str], Dict[int, List[int]]]:
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
    def _has_backbone_mark(sN: int, sC: int, sOXT: Optional[int],
                           serial2conect: Dict[int, List[int]]) -> bool:
        """判断链是否带 backbone 闭合标记: 头尾 N-C 间已有 CONECT 键, 或链尾残基有 OXT(线性C端残留)"""
        if sOXT:
            return True
        return (sC in serial2conect.get(sN, [])) or (sN in serial2conect.get(sC, []))

    def _detect_backbone_sites(self, atoms: Dict[int, Tuple[str, int, str, float, float, float]],
                               serial2conect: Dict[int, List[int]],
                               only_chains: Optional[List[str]] = None) -> List[BackboneSite]:
        """检测候选 backbone 闭合位点(仅头尾肽键, 不含二硫键):
        链的 C端残基 C(=O) 与 N端残基 N 距离 < min_cn.
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
            elif self._has_backbone_mark(sN, sC, sOXT, serial2conect):
                # 有 backbone 闭合标记(头尾 CONECT/OXT)但空间未闭合: 构象问题而非拓扑问题, 报告并跳过
                put_log(f'chain {chain}: backbone closure mark found but C-N dist={d:.1f}A >= min_cn, '
                        f'peptide is NOT spatially closed, skip (cannot fix topology on an open conformation)')
        return sites

    def _detect_disulfide_sites(self, atoms: Dict[int, Tuple[str, int, str, float, float, float]],
                                serial2conect: Dict[int, List[int]],
                                only_chains: Optional[List[str]] = None) -> List[DisulfideSite]:
        """检测候选二硫键位点: 同链内有 SG-SG 的 CONECT 记录, 或距离 < ss_max_dist.
        判据只看原子名 SG(不校验残基是否为 CYS/是否标准), 因此非标准含硫残基
        (如 LE1 青霉胺的 SG3)也能被检出. 注意: 上游 PDB 里 LE1 的硫原子名为 SG,
        在 fix-cp 之后才被 fix_ncAA 规范化为 SG3, 所以本函数看到的是 SG.
        only_chains: None 检测所有链; 有值时仅检测对应链, 链不存在时报告并跳过.
        逐对判断: 同一链内任意两个含 SG 的残基, 距离 < ss_max_dist 或已成键
        即视为一对二硫键, 每个 SG 只归入一对, 支持一条链上多对.
        ss_max_dist 为 None 时只按 CONECT 判定, 距离配对交给
        OpenMM/PDBFixer(固定 0.3 nm), 避免本函数与 OpenMM 阈值不一致.
        """
        chain_residues = defaultdict(dict)  # chain -> {resi: {name: serial}}
        for serial, (chain, resi, name, *_coords) in atoms.items():
            chain_residues[chain].setdefault(resi, {})[name] = serial

        if only_chains:
            for c in only_chains:
                if c not in chain_residues:
                    put_log(f'chain {c} not found in this complex, skip.')

        def sq_dist(s1, s2):
            x1 = atoms[s1][3:6]
            x2 = atoms[s2][3:6]
            return sum((a - b) ** 2 for a, b in zip(x1, x2)) ** 0.5

        sites = []
        for chain, residues in sorted(chain_residues.items()):
            if only_chains and chain not in only_chains:
                continue
            cys_resis = [r for r, d in residues.items() if 'SG' in d]
            if len(cys_resis) < 2:
                continue
            # 逐对配对: 只保留距离阈值内或已成键的组合, 每个 SG 只归入一对.
            # ss_max_dist 为 None 时不做距离判定, 只认 CONECT.
            ss_max_dist = self.args.ss_max_dist
            pairs = []
            for i in range(len(cys_resis)):
                for j in range(i + 1, len(cys_resis)):
                    r1, r2 = cys_resis[i], cys_resis[j]
                    sS1, sS2 = residues[r1]['SG'], residues[r2]['SG']
                    d_ss = sq_dist(sS1, sS2)
                    has_bond = ((sS2 in serial2conect.get(sS1, []))
                                or (sS1 in serial2conect.get(sS2, [])))
                    if has_bond or (ss_max_dist is not None and d_ss < ss_max_dist):
                        pairs.append((d_ss, r1, r2, sS1, sS2, has_bond))
            used = set()
            for d_ss, r1, r2, sS1, sS2, has_bond in sorted(pairs):
                if sS1 in used or sS2 in used:
                    continue
                used.update((sS1, sS2))
                sites.append({
                    'chain': chain,
                    'res1': r1, 'res2': r2,
                    'sS1': sS1, 'sS2': sS2,
                    'd_ss': d_ss,
                    'has_bond': has_bond,
                })
        return sites
    
    def _fix_backbone(self, pdb_path: Path) -> Optional[FixResult]:
        """backbone case: 删 OXT + 补头尾肽键 CONECT, 不动坐标. 无位点返回 None."""
        atoms, atom_lines, serial2conect = self._parse_pdb(str(pdb_path))
        sites = self._detect_backbone_sites(atoms, serial2conect, self.args.only_chains)

        if not sites:
            put_log(f'{pdb_path.name}: no backbone closure site found, skip.')
            return None
        backbone_chains = list(dict.fromkeys(s['chain'] for s in sites))
        
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
        return FixResult(lines=lines, backbone_chains=backbone_chains)

    def _fix_disulfide(self, pdb_path: Path) -> Optional[FixResult]:
        """修复二硫键环: 保留全部原子(不删 OXT), 仅确保 SG-SG 间存在 CONECT 键.
        无位点返回 None."""
        atoms, atom_lines, serial2conect = self._parse_pdb(str(pdb_path))
        sites = self._detect_disulfide_sites(atoms, serial2conect, self.args.only_chains)

        if not sites:
            put_log(f'{pdb_path.name}: no disulfide site found, skip.')
            return None
        disulfide_chains = list(dict.fromkeys(s['chain'] for s in sites))

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
        return FixResult(lines=lines, disulfide_chains=disulfide_chains)

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
        # 用 spawn 启动子进程: CUDA Runtime 不是 fork-safe 的, 从已加载 CUDA/pymol 的父进程
        # fork 出的子进程一旦初始化 CUDA 就会死锁, spawn 让子进程全新解释器启动以规避该问题.
        # 注: mbapy_lite 的 TaskPool 只把 mp_pool_init_kwargs 原样转发给 multiprocessing.Pool,
        # 而 multiprocessing.Pool 本身就是"绑定默认 context"的 BaseContext.Pool, 不接受 context
        # 参数(传了会 TypeError), 因此这里改为切换默认 start method.
        multiprocessing.set_start_method('spawn', force=True)
        pool = TaskPool('process', self.args.n_workers, report_error=True).start()
        relax_tasks = {}  # task_name -> out_path, 用于事后校验 relax 是否真的成功
                
        for task_i, pdb_path in tqdm(enumerate(pdb_paths), total=len(pdb_paths)):
            out_path = pdb_path.with_name(pdb_path.stem + self.args.suffix + '.pdb')
            if out_path.exists() and not self.args.force:
                continue
            
            try:
                result = self._fix_one(pdb_path)
            except Exception:
                traceback.print_exc()
                continue
            if result is None:
                continue
            # 所有 case 修复逻辑之后、relax 之前: 删除离群原子(可选)
            out_lines = self._remove_outlier_atoms_from_lines(result.lines, pdb_path)
            # 非标准氨基酸(ncAA)原子名规范化(如 LE1: SG->SG3, C8->CG1, C9->CG2),
            # 统一为 CHARMM36/GROMACS 约定, 供后续 OpenMM 模板匹配与 GROMACS pdb2gmx 识别
            from lazydock.opmm.custom_residues import CustomResidues
            from lazydock.opmm.fix_ncAA import fix_ncAA_pdb
            custom = CustomResidues(self.args.ncaa_xml)
            out_lines = fix_ncAA_pdb(out_lines, pdb_path, custom=custom)
            with open(out_path, 'w') as f:
                f.writelines(out_lines)
            self.printf(f'fixed: {pdb_path.name} -> {out_path.name}')
            if self.args.relax and (result.backbone_chains or result.disulfide_chains):
                # 两类环化可并存(不同链): disulfide 让 relax 启用 CHARMM36 DISU patch,
                # backbone 用 CustomBondForce 保护首尾肽键.
                task_name = pool.add_task(None, _fix_cp_relax_one, out_path,
                                          result.backbone_chains, result.disulfide_chains,
                                          restrain_chains=self.args.restrain_chain, restrain_backbone=self.args.restrain_backbone,
                                          constraints=self.args.constraints, max_iter=self.args.max_iter,
                                          tolerance=self.args.tolerance, platform=self.args.platform,
                                          gpu_index=gpu_slots[task_i % len(gpu_slots)],
                                          custom_residue_xml=self.args.ncaa_xml,
                                          # --ss-max-dist 单位 A, relax 侧单位 nm; 未传(None)时用 relax 默认 0.3nm
                                          ss_max_dist=(self.args.ss_max_dist / 10.0
                                                       if self.args.ss_max_dist is not None else 0.3))
                relax_tasks[task_name] = out_path
            pool.wait_till_free()
        pool.wait_till_all_done()
        # relax 失败时磁盘上留下的仍是 pre-relax(仅拓扑修复)版本, 不能以成功流程静默结束.
        # _fix_cp_relax_one 成功时返回 None, 因此任何非 None 结果(异常对象或 TaskStatus)都视为失败.
        failed = []
        for task_name, out_path in relax_tasks.items():
            result = pool.query_task(task_name)
            if result is not None:
                failed.append((out_path, result))
        pool.close(1)
        if failed:
            for out_path, err in failed:
                put_err(f'relax FAILED for {out_path.name}: {err!r} '
                        f'-- file on disk is the pre-relax (topology-fixed only) version',
                        warning_level=3)
            put_err(f'{len(failed)}/{len(relax_tasks)} relax task(s) failed, '
                    f'outputs are NOT relaxed', warning_level=3, _exit=1)

    def _fix_one(self, pdb_path: Path) -> Optional[FixResult]:
        """按 case 分发单个文件修复: backbone / disulfide / auto. 无命中返回 None."""
        if self.args.case == 'disulfide':
            return self._fix_disulfide(pdb_path)
        if self.args.case == 'backbone':
            return self._fix_backbone(pdb_path)
        # auto: 前置检测分流.
        # 二硫键本身不需要 fix-cp 写 CONECT: OpenMM/PDBFixer 按 `distance < 0.3 nm`
        # 自动成键(标准 CYS), 非标准含硫残基则由 relax 的 disulfide_chains 缝接.
        # 因此检出二硫键的真正价值是把链名报给 relax, 而非改拓扑.
        atoms, _, serial2conect = self._parse_pdb(str(pdb_path))
        disu_sites = self._detect_disulfide_sites(atoms, serial2conect, self.args.only_chains)
        bb_sites = self._detect_backbone_sites(atoms, serial2conect, self.args.only_chains)
        disu_chains = list(dict.fromkeys(s['chain'] for s in disu_sites))
        if bb_sites:
            # 头尾闭合必须在 fix-cp 阶段补 CONECT(OpenMM 不会推断), 因此优先走 backbone:
            # OXT 需删除、首尾肽键需补上, 这两件事只有 _fix_backbone 做.
            put_log(f'{pdb_path.name}: auto -> backbone'
                    + (f' + disulfide {disu_chains}' if disu_chains else ''))
            result = self._fix_backbone(pdb_path)
            if result is not None:
                result.disulfide_chains = disu_chains
            return result
        elif disu_sites:
            put_log(f'{pdb_path.name}: auto -> disulfide {disu_chains}')
            return self._fix_disulfide(pdb_path)
        else:
            put_log(f'{pdb_path.name}: auto: neither disulfide nor backbone closure detected, skip.')
            return None


_str2func = {
    'smiles2pdb': smiles2pdb,
    'seq2pdb': seq2pdb,
    'cif2pdb': cif2pdb,
    'fix-cp': fix_cp,
}


def main(sys_args: Optional[List[str]] = None) -> None:
    args_paser = argparse.ArgumentParser(description='Prepare ligand files for docking')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))
    
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main()