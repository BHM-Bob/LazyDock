'''
Date: 2026-08-28 10:00:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2026-08-28 10:30:00
Description: Check PDB structure for atom clashes and covalent bond length errors using OpenBabel
'''

import argparse
from pathlib import Path
from typing import List

import numpy as np
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.web_utils.task import TaskPool
from openbabel import openbabel as ob
from scipy.spatial import KDTree
from tqdm import tqdm

# openbabel 读取时对元素列缺失等格式问题会刷警告, 此处静默(0=quiet)
ob.obErrorLog.SetOutputLevel(0)

from lazydock.scripts._script_utils_ import Command, excute_command

# openbabel 中未知原子序数的 VdW 半径返回 0.0，此时用这个默认值兜底
DEFAULT_VDW_RADII = 1.7


def atom_repr(atom: 'ob.OBAtom') -> str:
    """生成原子的可读描述: RESNAME:ATOMNAME(chain:resnum)"""
    residue = atom.GetResidue()
    if residue is None:
        return f'ATOM{atom.GetIdx()}({atom.GetAtomicNum()})'
    return f'{residue.GetName()}:{residue.GetAtomID(atom).strip()}({residue.GetChain()}:{residue.GetNum()})'


# ============ 模块级检测函数 ============
# 设计为仅依赖 args(argparse Namespace)与显式参数, 不依赖类实例,
# 模块级函数可被 pickle, 便于后续用 multiprocessing 做并行检测.

def report_atom(atom: 'ob.OBAtom', report_atom_nums) -> bool:
    """--report-atoms 过滤: 原子是否属于指定元素"""
    if report_atom_nums is None:
        return True
    return atom.GetAtomicNum() in report_atom_nums


def keep_atom(atom: 'ob.OBAtom', skip_chains, only_chains) -> bool:
    """链过滤: 按 --skip-chains / --only-chains 决定原子是否保留"""
    r = atom.GetResidue()
    chain = r.GetChain() if r is not None else ''
    if skip_chains and chain in skip_chains:
        return False
    if only_chains and chain not in only_chains:
        return False
    return True


def check_atom_clashes(mol: 'ob.OBMol', args, report_atom_nums=None) -> List[dict]:
    """基于 openbabel GetVdwRad/GetCovalentRad 的原子碰撞检测

    - 键连(1-2)原子对: 不在此检测, 由 bond length 检测(GetEquibLength)评估, 共价键用键长评估最专业
    - 1-3(隔一个原子)原子对: 距离 < bond_overlap_scale * (cov_r_i + cov_r_j) 视为碰撞, 严格检测
    - 非键原子对: 距离 < vdw_scale * (vdw_r_i + vdw_r_j) 视为碰撞
    - 极性 H 与受体(N/O/S/F/Cl)之间距离 > hbond_exempt 视为正常氢键, 不报为 clash
    - 每条 clash 按 overlap(limit - distance) 分级: > severe_overlap 为 severe, 其余为 mild
    - overlap < min_overlap 的轻微接触不报告
    """
    # 受体元素集合
    acceptor_nums = {ob.GetAtomicNum(s) for s in ['N', 'O', 'S', 'F', 'Cl']}
    atoms = [a for a in ob.OBMolAtomIter(mol)]
    # 链过滤: 生成保留掩码
    keep = np.array([keep_atom(a, args.skip_chains, args.only_chains) for a in atoms], dtype=bool)
    atoms = [a for i, a in enumerate(atoms) if keep[i]]
    coords = np.array([[a.GetX(), a.GetY(), a.GetZ()] for a in atoms])

    # KDTree 候选对: 依据 VdW 半径和(cov 半径和更小, 能被覆盖)
    max_r = max(float(ob.GetVdwRad(a.GetAtomicNum()) or DEFAULT_VDW_RADII) for a in atoms)
    kdtree = KDTree(coords)
    candidate_pairs = kdtree.query_pairs(2 * max_r * args.vdw_scale + 1e-3)

    clashes = []
    for i, j in candidate_pairs:
        a1, a2 = atoms[i], atoms[j]
        # --report-atoms: 至少一个原子属于指定元素才报告
        if not (report_atom(a1, report_atom_nums) or report_atom(a2, report_atom_nums)):
            continue
        if args.ignore_h and (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1):
            continue
        if args.no_hh and a1.GetAtomicNum() == 1 and a2.GetAtomicNum() == 1:
            continue
        # 键连(1-2)原子对: 由 bond length 检测(GetEquibLength)评估
        if a1.IsConnected(a2):
            continue
        dist = float(np.linalg.norm(coords[i] - coords[j]))
        # 1-3 关系(隔一个原子): 非键但受拓扑约束, 用共价半径和阈值严格检测
        if any(nb.IsConnected(a2) for nb in ob.OBAtomAtomIter(a1)):
            limit = args.bond_overlap_scale * (
                (ob.GetCovalentRad(a1.GetAtomicNum()) or 0.0) +
                (ob.GetCovalentRad(a2.GetAtomicNum()) or 0.0))
            relation = '1-3'
        else:
            limit = args.vdw_scale * (
                (ob.GetVdwRad(a1.GetAtomicNum()) or DEFAULT_VDW_RADII) +
                (ob.GetVdwRad(a2.GetAtomicNum()) or DEFAULT_VDW_RADII))
            relation = 'non-bond'
        if dist < limit:
            # 氢键豁免: 极性 H 与受体(N/O/S/F/Cl)距离 > hbond_exempt 视为正常氢键
            if args.hbond_exempt > 0:
                a1_is_h = a1.GetAtomicNum() == 1
                a2_is_h = a2.GetAtomicNum() == 1
                if a1_is_h and not a2_is_h:
                    h_atom, other = a1, a2
                elif a2_is_h and not a1_is_h:
                    h_atom, other = a2, a1
                else:
                    h_atom = other = None
                if h_atom is not None and h_atom.IsPolarHydrogen() \
                        and other.GetAtomicNum() in acceptor_nums:
                    if dist > args.hbond_exempt:
                        continue
            overlap = limit - dist
            if overlap < args.min_overlap:
                continue
            clashes.append({
                'atom1': atom_repr(a1),
                'atom2': atom_repr(a2),
                'relation': relation,
                'severity': 'severe' if overlap > args.severe_overlap else 'mild',
                'distance': round(dist, 3),
                'limit': round(limit, 3),
                'overlap': round(overlap, 3),
            })
    return clashes


def check_bond_lengths(mol: 'ob.OBMol', args, report_atom_nums=None) -> List[dict]:
    """基于 openbabel GetEquibLength 的共价键键长检测"""
    errors = []
    for bond in ob.OBMolBondIter(mol):
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if not (keep_atom(a1, args.skip_chains, args.only_chains)
                and keep_atom(a2, args.skip_chains, args.only_chains)):
            continue
        # --report-atoms: 至少一个原子属于指定元素才报告
        if not (report_atom(a1, report_atom_nums) or report_atom(a2, report_atom_nums)):
            continue
        if args.ignore_h and (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1):
            continue
        length = bond.GetLength()
        eq_length = bond.GetEquibLength()
        if eq_length <= 0 or length <= 0:
            continue
        diff = abs(length - eq_length)
        tol = max(args.bond_tol * eq_length, args.bond_tol_abs)
        if diff > tol:
            errors.append({
                'atom1': atom_repr(a1),
                'atom2': atom_repr(a2),
                'bond_order': bond.GetBondOrder(),
                'distance': round(length, 3),
                'equil_length': round(eq_length, 3),
                'diff': round(diff, 3),
                'tolerance': round(tol, 3),
            })
    return errors


def check_one_file(pdb_path, args, report_atom_nums=None) -> dict:
    """对单个 PDB 文件执行完整检测, 返回结果 dict(供并行调用)"""
    conv = ob.OBConversion()
    if not conv.SetInAndOutFormats('pdb', 'pdb'):
        put_err(f'failed to set pdb format, exit.', _exit=True)
    mol = ob.OBMol()
    if not conv.ReadFile(mol, str(pdb_path)):
        put_err(f'failed to read pdb file: {pdb_path}', _exit=True)
    if args.add_h:
        mol.AddHydrogens()  # openbabel 补全缺失的 H

    n_kept = sum(1 for a in ob.OBMolAtomIter(mol)
                 if keep_atom(a, args.skip_chains, args.only_chains))
    n_bonds = sum(1 for b in ob.OBMolBondIter(mol)
                  if keep_atom(b.GetBeginAtom(), args.skip_chains, args.only_chains)
                  and keep_atom(b.GetEndAtom(), args.skip_chains, args.only_chains))
    result = {'file': str(pdb_path), 'atoms': n_kept, 'bonds': n_bonds}
    if not args.skip_clash:
        result['clashes'] = check_atom_clashes(mol, args, report_atom_nums)
        result['n_clashes'] = len(result['clashes'])
    if not args.skip_bonds:
        result['bond_errors'] = check_bond_lengths(mol, args, report_atom_nums)
        result['n_bond_errors'] = len(result['bond_errors'])
    return result


class check_structure(Command):
    HELP = """Check PDB structure for atom clashes and covalent bond length errors using OpenBabel.

Atom clash: 
  bonded(1-2) atom pair: NOT checked here, evaluated by bond length check (GetEquibLength, the most professional way).
  1-3 atom pair: distance < bond_overlap_scale * (covalent_radius(i) + covalent_radius(j)),
  reference from OpenBabel GetCovalentRad, strictly checked.
  non-bonded atom pair: distance < vdw_scale * (vdw_radius(i) + vdw_radius(j)), reference from OpenBabel GetVdwRad.
  Polar H to acceptor (N/O/S/F/Cl) pairs with distance > hbond_exempt are normal hydrogen bonds, not clashes.
Clashes are graded by overlap (limit - distance): severe if > severe_overlap, else mild;
overlap < min_overlap contacts are not reported.
Bond length: |length - equil_length| > bond_tol * equil_length (percent) or > bond_tol_abs (Angstrom),
equil_length reference from OpenBabel OBBond.GetEquibLength.

Output modes (--summary, independent of -o):
  None   (default): print all details to console only when -o is not given
  each   : for each file, print clash/bond error counts and ratios
  detail : always print all details, even when -o is given

Filtering:
  --report-atoms C N O: only report clash/bond errors involving at least one atom of the given element symbols."""
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains input pdb files or sub-folders containing input pdb files, default is %(default)s.")
        args.add_argument('-n', '--name', type=str, default='',
                          help="input pdb file name substring, default is %(default)s.")
        args.add_argument('-o', '--output', type=str, default=None,
                          help='output json file path, default is %(default)s, which means print to console.')
        args.add_argument('--summary', type=str, choices=['each', 'detail'], default=None,
                          help='console output mode: each for per-file counts with ratio, '
                               'detail for all details, default is %(default)s (same as current behavior).')
        args.add_argument('--skip-clash', action='store_true',
                          help='skip atom clash check.')
        args.add_argument('--skip-bonds', action='store_true',
                          help='skip bond length check.')
        args.add_argument('--vdw-scale', type=float, default=0.8,
                          help='scale factor of vdw radii sum for clash detection, default is %(default)s.')
        args.add_argument('--bond-overlap-scale', type=float, default=0.7,
                          help='scale factor of covalent radii sum for bonded(1-2)/1-3 atom pair clash detection, '
                               'default is %(default)s. Bonded atoms in normal bond length (~covalent radii sum) '
                               'are always shorter than vdw radii sum, so this threshold detects severe overlaps '
                               'of bonded/1-3 atoms without false positives on normal bonds.')
        args.add_argument('--no-hh', action='store_true',
                          help='ignore H-H atom pair clashes (H-H pairs of e.g. methyl group are naturally close).')
        args.add_argument('--hbond-exempt', type=float, default=1.5,
                          help='polar H (bonded to N/O/S) to acceptor (N/O/S/F/Cl) pairs with distance greater than '
                               'this value are treated as normal hydrogen bonds and NOT reported as clash. '
                               'Set to 0 to disable, default is %(default)s Angstrom.')
        args.add_argument('--min-overlap', type=float, default=0.0,
                          help='only report clashes with overlap (limit - distance) greater than this value, '
                               'default is %(default)s.')
        args.add_argument('--severe-overlap', type=float, default=0.6,
                          help='clashes with overlap greater than this value are marked as severe, '
                               'default is %(default)s Angstrom.')
        args.add_argument('--bond-tol', type=float, default=0.15,
                          help='bond length error tolerance (fraction of equil length), default is %(default)s.')
        args.add_argument('--bond-tol-abs', type=float, default=0.3,
                          help='bond length error absolute tolerance in Angstrom, default is %(default)s.')
        args.add_argument('--ignore-h', action='store_true',
                          help='ignore hydrogen atoms in both checks.')
        args.add_argument('--add-h', action='store_true',
                          help='add missing hydrogens before check (default is False, keep PDB atoms as is).')
        args.add_argument('--only-chains', type=str, nargs='+', default=None,
                          help='only check given chains, e.g. --only-chains A B, default is %(default)s.')
        args.add_argument('--skip-chains', type=str, nargs='+', default=None,
                          help='skip given chains, default is %(default)s.')
        args.add_argument('--report-atoms', type=str, nargs='+', default=None,
                          help='only report clash/bond errors involving at least one atom of these element symbols, '
                               'e.g. --report-atoms C N O, default is %(default)s (report all).')
        args.add_argument('-nw', '--n-workers', type=int, default=1,
                          help='number of workers for parallel check, default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
        if self.args.skip_clash and self.args.skip_bonds:
            put_err('skip-clash and skip-bonds are both set, nothing to check, exit.', _exit=True)
        if self.args.output:
            self.args.output = Path(self.args.output).resolve()
            if not self.args.output.parent.exists():
                self.args.output.parent.mkdir(parents=True, exist_ok=True)
        # 校验元素符号并转为小写 set, 用于 --report-atoms 过滤
        if self.args.report_atoms:
            invalid = [s for s in self.args.report_atoms if ob.GetAtomicNum(s) == 0]
            if invalid:
                put_err(f'invalid element symbol(s): {invalid}, exit.', _exit=True)
            self.report_atom_nums = set(ob.GetAtomicNum(s) for s in self.args.report_atoms)
        else:
            self.report_atom_nums = None
    
    def _print_detail(self, result: dict):
        """打印单个文件的详细结果(标题行 + 所有 clash/bond error 详情)"""
        n_clash = result.get('n_clashes', 0)
        n_bond = result.get('n_bond_errors', 0)
        self.printf(f"\n== {Path(result['file'])}: {result['atoms']} atoms, "
                    f"{n_clash} clash(es), {n_bond} bond error(s) ==")
        for clash in result.get('clashes', []):
            self.printf(f"  clash: {clash['atom1']} -- {clash['atom2']} "
                        f"[{clash['relation']}/{clash['severity']}]: "
                        f"dist={clash['distance']}, limit={clash['limit']}, overlap={clash['overlap']}")
        for err in result.get('bond_errors', []):
            self.printf(f"  bond: {err['atom1']} -- {err['atom2']} (order {err['bond_order']}): "
                        f"len={err['distance']}, equil={err['equil_length']}, tol={err['tolerance']}")
    
    def main_process(self):
        pdb_paths = [Path(p).resolve() for p in get_paths_with_extension(
            self.args.batch_dir, ['.pdb'], name_substr=self.args.name)]
        put_log(f'get {len(pdb_paths)} pdb file(s) in {self.args.batch_dir}')
        
        tasks = []
        pool = TaskPool('process', self.args.n_workers, report_error=True).start()
        for pdb_path in tqdm(pdb_paths, total=len(pdb_paths)):
            tasks.append(pool.add_task(pdb_path, check_one_file, pdb_path, self.args, self.report_atom_nums))
            pool.wait_till_free()
            
        pool.wait_till_all_done()
        results = [pool.query_task(t_i) for t_i in tasks]
        pool.close(1)
        
        if self.args.output:
            opts_file(self.args.output, 'w', way='json', data=results)
            self.printf(f'results saved to: {self.args.output}')
        
        # console output control by --summary, independent of -o
        if self.args.summary == 'each':
            # 逐文件打印 clash/bond error 数量与比例
            for result in results:
                atoms, bonds = result['atoms'], result['bonds']
                n_clash = result.get('n_clashes', 0)
                n_bond = result.get('n_bond_errors', 0)
                n_severe = sum(1 for c in result.get('clashes', []) if c.get('severity') == 'severe')
                clash_pct = n_clash / atoms * 100 if atoms else 0
                bond_pct = n_bond / bonds * 100 if bonds else 0
                self.printf(f"{Path(result['file'])}: {atoms} atoms, {bonds} bonds, "
                            f"{n_clash} clash(es) ({clash_pct:.2f}% of atoms, {n_severe} severe), "
                            f"{n_bond} bond error(s) ({bond_pct:.2f}% of bonds)")
        elif self.args.summary == 'detail':
            # 打印所有详细信息, 与 -o 并列, 即使指定了输出文件也打印
            for result in results:
                self._print_detail(result)
        elif not self.args.output:
            # None: 当前行为, 无 -o 时打印所有详情
            for result in results:
                self._print_detail(result)
        
        # 汇总统计
        if len(results) > 1:
            total_clash = sum(r.get('n_clashes', 0) for r in results)
            total_bond = sum(r.get('n_bond_errors', 0) for r in results)
            put_log(f'total: {len(results)} file(s), {total_clash} clash(es), {total_bond} bond error(s)')


_str2func = {
    'check-structure': check_structure,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description='Check PDB structure for atom clashes and covalent bond length errors')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k, description=v.HELP))
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main()