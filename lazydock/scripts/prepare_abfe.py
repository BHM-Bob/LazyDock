'''
Date: 2026-08-28
LastEditors: BHM-Bob 2262029386@qq.com
Description: prepare ABFE systems (complex + ligand topologies).
             Supports receptor multi-chain and peptide ligands (as protein).
'''
import argparse
import logging
import os
from pathlib import Path
from typing import List, Union

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension

# system_builder imported lazily in _get_builder (avoids openff probe on CLI parse)
from lazydock.gmx.run import Gromacs
from lazydock.gmx.abfe import system_builder
from lazydock.scripts._script_utils_ import Command, make_args_and_excute, process_batch_dir_lst
from tqdm import tqdm


def _fix_ligand_for_abfe(pdb_path: Path, ligand_resname: str = 'LIG') -> Path:
    """Return the pdb path with the ligand residue renamed (in-place edit of a copy)."""
    from mbapy_lite.file import opts_file
    lines = opts_file(pdb_path, way='lines')
    new_lines = []
    for line in lines: # type: ignore
        if line.startswith(('ATOM', 'HETATM')):
            line = f'{line[:17]}{ligand_resname:<3s}{line[20:]}'
        new_lines.append(line)
    out = pdb_path.with_suffix('.ligrenamed.pdb')
    opts_file(out, 'w', way='lines', data=new_lines)
    return out


def _auto_solv_d(pdb_path: Path, rlist: float, pad: float = 1.2,
                 margin: float = 0.3) -> float:
    """Derive editconf -d (nm) for one leg such that the box
    half-shortest-vector >= rlist + margin.

    For cubic / dodecahedron / octahedron boxes (editconf -bt) the
    half-shortest-vector is a/2 with a the box edge (measured, 2026-09-11),
    and editconf -d grows the box isotropically: a = span + 2*d.
    ->  a = 2 * max(rlist + margin, span/2 + pad)
        d = a/2 - span/2
    """
    from lazydock.gmx.box import measure_span
    span = measure_span(pdb_path)
    half = max(rlist + margin, span / 2 + pad)
    d = half - span / 2
    put_log(f'auto solv-d: span={span:.2f} nm, rlist={rlist:.2f} -> '
            f'half-vector={half:.2f} nm, d={d:.2f} nm', head='ABFE')
    return d


def auto_solv_d_pair(ligand_pdb: Path, complex_leg_span_src: Path, peptide: bool,
                     fep_rlist_ligand: float = None) -> tuple:
    """Return (d_ligand, d_complex) for --solv-d-auto.

    peptide: ligand-leg rlist is fixed >= 3.1 (iron rule, never lowered);
             complex-leg rlist = same as ligand (symmetric, user decision 2).
    """
    rlist_lig = fep_rlist_ligand or (3.1 if peptide else 1.2)
    # peptide iron rule: never below 3.1
    if peptide:
        rlist_lig = max(rlist_lig, 3.1)
    # complex leg rlist = ligand rlist (symmetric), padding 1.2
    d_lig = _auto_solv_d(ligand_pdb, rlist_lig)
    d_com = _auto_solv_d(complex_leg_span_src, rlist_lig)
    put_log(f'--solv-d-auto: ligand d={d_lig:.2f}, complex d={d_com:.2f} nm', head='ABFE')
    return d_lig, d_com


def _check_solv_d_conflict(args) -> None:
    """The three box-mode flags are mutually exclusive (user review 2026-09-12):
      --solv-d  : manual isotropic editconf -d (both legs, same value)
      --solv-d-auto : auto isotropic editconf -d per leg (shape via --solv-bt)
      --auto-box    : auto anisotropic contour box (triclinic, each axis >= rlist)
    They express different box-building strategies; combining them is ambiguous,
    so error out instead of silently resolving (auto-box no longer implies
    solv-d-auto: the contour box has its own rlist handling)."""
    if args.solv_d_auto and args.solv_d is not None:
        put_err('--solv-d-auto and --solv-d are mutually exclusive: both provided. '
                'Remove one (auto is preferred for correctness).', _exit=True)
    if args.auto_box and args.solv_d_auto:
        put_err('--auto-box and --solv-d-auto are mutually exclusive: both provided. '
                'They are DIFFERENT box strategies (anisotropic contour vs isotropic '
                '-d, see help). Pick one.', _exit=True)
    if args.auto_box and args.solv_d is not None:
        put_err('--auto-box and --solv-d are mutually exclusive: both provided. '
                'The contour box derives its axes from rlist; a manual --solv-d '
                'conflicts with it.', _exit=True)


class Complex(Command):
    HELP = """
    prepare ABFE system: complex + ligand.
    
    INPUT:
        complex.pdb with two (or more) chains: receptor chain(s) + ligand chain.
        - receptor can have MULTIPLE chains (--receptor-chains A B ...)
        - ligand can be a PEPTIDE (--peptide, treated as protein via pdb2gmx)
          or a small molecule (espaloma/openff/gaff parameterization)

    OUTPUT (BindFlow-compatible dir structure):
        input/complex/  complex.gro complex.top + itps + index.ndx
        input/ligand/   ligand.gro ligand.top + itps
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf, iter_run_arg=['dir'])

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, nargs='+', default=['.'],
                          help='directory containing complex.pdb. Default is %(default)s.')
        args.add_argument('-n', '--name', type=str, default='complex.pdb',
                          help='complex structure file name in each dir. Default is %(default)s.')
        args.add_argument('-rc', '--receptor-chains', type=str, nargs='+', default=None,
                          help='receptor chain name(s), e.g. A or A B C (multi-chain supported). '
                               'Optional: if omitted, receptor chains are auto-derived as '
                               'ALL chains - ligand chain.')
        args.add_argument('-lc', '--ligand-chain', type=str, required=True,
                          help='ligand chain name, e.g. Z or P (REQUIRED, no auto-detection: '
                               'avoids ambiguity; receptor chains are derived from it).')
        args.add_argument('--peptide', action='store_true',
                          help='treat ligand as a PEPTIDE (protein-protein form, pdb2gmx).')
        args.add_argument('--ligand-ff', type=str, default='espaloma', choices=['espaloma', 'openff', 'gaff'],
                          help='small molecule force field backend (ignored if --peptide). Default %(default)s.')
        args.add_argument('--protein-ff', type=str, default='amber99sb-ildn',
                          help='force field for pdb2gmx (protein/peptide), e.g. amber99sb-ildn, '
                               'charmm36-mar2019, charmm36-jul2022. '
                               'If a PATH to an existing force field DIRECTORY (e.g. '
                               '/data/ff/charmm36-jul2022.ff) is given, it is used as an '
                               'external force field (GMXLIB points to it; the code is taken '
                               'from the directory name). Otherwise it is treated as a GROMACS '
                               'built-in force field code. Default %(default)s.')
        args.add_argument('--pdb2gmx-args', type=str, default='-ter -ignh',
                          help='args passed to pdb2gmx. Default %(default)s.')
        args.add_argument('--water-model', type=str, default='amber/tip3p',
                          help='water model. Default %(default)s.')
        args.add_argument('--solv-d', type=float, default=None,
                          help='editconf -d (nm), single value applied to BOTH ligand and '
                               'complex legs. Default 1.5. Conflicts with --solv-d-auto '
                               '(error out).')
        args.add_argument('--solv-d-auto', action='store_true',
                          help='auto-derive editconf -d per leg (ligand/complex) from the '
                               'system span and the FEP rlist, guaranteeing '
                               'half-shortest-box-vector >= rlist. '
                               'ISOTROPIC: the box grows equally in all directions '
                               '(shape set by --solv-bt: cubic/dodecahedron/octahedron). '
                               'Conflicts with --solv-d and --auto-box (both error out).')
        args.add_argument('--solv-bt', type=str, default='cubic',
                          choices=['cubic', 'dodecahedron', 'octahedron'],
                          help='box type for --solv-d-auto (isotropic editconf -d). '
                               'Default %(default)s (editconf default). '
                               'All choices share half-shortest-vector = a/2, so at the '
                               'same a the water volume ranks: cubic a^3 (most), '
                               'octahedron ~0.77*a^3, dodecahedron ~0.707*a^3 (least, '
                               'measured 2026-09-11). '
                               'NOTE: --solv-d-auto CANNOT adapt to the solute contour '
                               '(isotropic); use --auto-box for anisotropic contour boxes.')
        args.add_argument('--auto-box', action='store_true',
                          help='use PyMOL contour box (ANISOTROPIC triclinic/orthogonal, '
                               'each axis = max(span_i + 2*pad, 2*(rlist+0.3)) >= rlist) '
                               'instead of the isotropic editconf -d. For elongated '
                               'receptors (GPCR/multimers) saves water (~40%% vs cubic). '
                               'This is a *different branch* from --solv-d-auto: contour '
                               'vs isotropic. Mutually exclusive with --solv-d-auto and '
                               '--solv-d (error out if combined).')
        args.add_argument('--auto-box-pad', type=float, default=1.2,
                          help='padding (nm) for --auto-box contour (per axis, both sides). '
                               'Default %(default)s.')
        args.add_argument('--ion-conc', type=float, default=150e-3, help='ion concentration (M). Default %(default)s.')
        args.add_argument('--hmr-factor', type=float, default=2.5,
                          help='HMR factor (Hydrogen Mass Repartition). Default %(default)s '
                               'matching BindFlow FEP driver (fep_full_run.py hard-codes 2.5). '
                               'Must be used with 4 fs FEP timestep (dt_max >= 0.004). '
                               'Pass 0 to disable HMR (then use dt_max <= 0.002).')
        args.add_argument('--maxwarn', type=int, default=0, help='maxwarn for grompp. Default %(default)s.')
        args.add_argument('--builder-dir', type=str, default='builder', help='builder dir name. Default %(default)s.')
        args.add_argument('--fep-rlist-ligand', type=float, default=None,
                          help='rlist (nm) required for the ligand-leg FEP windows. Used by the '
                               'automatic box derivation to guarantee '
                               'half-shortest-box-vector >= rlist. Default: peptide 3.1; '
                               'small molecule 1.2.')
        return args

    def process_args(self):
        self.args.dir = process_batch_dir_lst(self.args.dir)

    def main_process(self):
        # --protein-ff: 单参数判定 (2026-09-08, Phase 2):
        #   * 已存在的目录路径 → 外部力场: 设置 GMXLIB (MakeInputs.custom_ff_path),
        #     力场 code 取目录名 stem (如 charmm36-jul2022.ff → charmm36-jul2022)
        #   * 其他字符串 → GROMACS 内置力场 code (custom_ff_path=None)
        ff_arg = self.args.protein_ff
        if os.path.isdir(ff_arg):
            self.args.custom_ff_path = Path(ff_arg).resolve()
            self.args.protein_ff = Path(ff_arg.rstrip('/')).stem
            put_log(f'--protein-ff is an existing directory: external force field, '
                    f'GMXLIB parent handling in MakeInputs, code = {self.args.protein_ff}')
        else:
            self.args.custom_ff_path = None
        # 盒子模式互斥检查: 必须在找 complex / 判 skip 之前 (参数合法性优先于
        # input/ 存在与否, 否则已存在 input 的目录会跳过检查而静默接受冲突 flag)
        _check_solv_d_conflict(self.args)
        # find complex files (REAL batch: process EVERY match, not just the first)
        complex_files = []
        for _dir in (self.args.dir if isinstance(self.args.dir, list) else [self.args.dir]):
            if os.path.isdir(_dir):
                complex_files.extend(get_paths_with_extension(_dir, [],
                                                              name_substr=self.args.name, sort='natsort'))
            else:
                put_err(f'dir argument should be a directory: {_dir}', _exit=True)
        if not complex_files:
            put_err(f'no {self.args.name} found in {self.args.dir}', _exit=True)
        put_log(f'found {len(complex_files)} complex file(s) to process: '
                + ', '.join(str(Path(p)) for p in complex_files))

        skipped, failed = 0, 0
        for complex_path in tqdm(map(Path, complex_files), desc='ABFE prepare', total=len(complex_files)):
            complex_path = complex_path.resolve()
            wdir = complex_path.parent
            # skip if input/ already exists (idempotent batch re-run)
            if (wdir / 'input').exists():
                put_log(f'skip {complex_path}: input/ already exists.', head='ABFE')
                skipped += 1
                continue
            try:
                self._prepare_single(complex_path, wdir)
            except NotImplementedError:
                # 明确的功能缺失 (如小分子路径) 不能静默吞掉 - 抛出让用户看到
                raise
            except Exception as e:
                put_err(f'processing {complex_path} failed: {e}')
                failed += 1
        if skipped or failed:
            put_log(f'batch done: processed {len(complex_files)}, skipped {skipped} (input exists), failed {failed}.')
        else:
            put_log(f'ABFE prepare completed for {len(complex_files)} complex file(s).')

    def _prepare_single(self, complex_path: Path, wdir: Path):
        """Run the full prepare pipeline for one complex.pdb (single batch member).

        Peptide mode (new "three-stage" architecture, 2026-09-11):
          1. pdb2gmx: 整链 complex.pdb 一次 prepare-gmx protein (-merge no 保留两链
             moleculetype Protein_chain_X, 相对坐标天然保留 - 修复 Phase 2 的
             受体/配体独立居中导致的重叠 bug);
          2. ABFE 后处理: MakeInputs 读整链 gro/top, 按原子数识别配体 moleculetype
             并 rename 为 LIG (couple-moltype 用), posres 同步;
          3. 溶剂化: ABFE Solvate (自定义水模型/HMR/rlist/auto-box 全保留)。
        Small-molecule mode: 保持原 toff 路径 (espaloma/openff/gaff), 不动。
        """
        put_log(f'processing ABFE prepare for: {complex_path}', head='ABFE')

        if not self.args.peptide:
            # 小分子路径未实现 (2026-09-11, user review #4): complex 腿需要
            # prepare-gmx complex (toff/CGenFF 拼装) 支持, 目前不存在, 直接报错
            # 而非用 protein 路径误处理含小分子的 complex.pdb。
            raise NotImplementedError(
                'small-molecule ABFE prepare is NOT implemented yet: the complex '
                'leg requires prepare-gmx complex (toff/CGenFF) support. '
                'See docs/dev/abfe/small_mol_future_prepare_gmx_complex.md. '
                'Use --peptide for peptide ligands.')

        # 链推导: ligand chain (必填, 不自动探测), receptor = ALL - ligand (可选显式覆盖)
        from pymol import cmd
        import shutil
        cmd.reinitialize()
        cmd.load(str(complex_path), 'complex')
        all_chains = list(cmd.get_chains('complex'))
        lig_chain = self.args.ligand_chain
        if lig_chain not in all_chains:
            put_err(f'ligand chain {lig_chain} not found in PDB (chains: {all_chains}).', _exit=True)
        if self.args.receptor_chains:
            rec_chains = self.args.receptor_chains
        else:
            rec_chains = [c for c in all_chains if c != lig_chain]
        if not rec_chains:
            put_err(f'no receptor chains left after removing ligand chain {lig_chain}. '
                    f'Use --receptor-chains explicitly.', _exit=True)
        put_log(f'chains: all={all_chains}, receptor={rec_chains}, ligand={lig_chain}', head='ABFE')

        # ligand.pdb: 小分子路径切出 (toff 参数化用) + peptide 路径的 ligand 腿
        # (独立 protein, 无相对坐标问题)。peptide complex 腿直接用整链 complex.pdb。
        ligand_pdb = wdir / 'ligand.pdb'
        if not ligand_pdb.exists():
            if cmd.select('ligand', f'complex and chain {lig_chain}') == 0:
                put_err(f'ligand chain {lig_chain} has zero atoms.', _exit=True)
            cmd.save(str(ligand_pdb), 'ligand')
            put_log(f'ligand saved: {ligand_pdb}', head='ABFE')
        if not self.args.peptide:
            # 小分子: 残基名强制 LIG (toff 参数化按 resname 处理)
            ligand_pdb = _fix_ligand_for_abfe(ligand_pdb, 'LIG')
            put_log(f'ligand resname renamed to LIG: {ligand_pdb}', head='ABFE')

        # -- pdb2gmx via prepare-gmx (Phase 2 复用, 2026-09-11 整链重构)
        #    import prepare_gmx.main() 传 list-str 参数; 出错会抛异常 →
        #    上层 try-catch 捕获。
        #    prepare-gmx protein 产出: {stem}.gro + topol.top (工作目录内)。
        #    力场目录沿用 prepare-gmx 的 --ff-dir (copy ff 到 cwd, 过渡期策略)。
        from lazydock.scripts.prepare_gmx import main as prepare_gmx_main

        def _run_prepare_gmx_protein(wdir_p, pdb_name, merge_val='no'):
            """Run prepare-gmx protein in wdir_p for pdb_name; skip if gro/top exist.

            - peptide complex 腿: 整链 complex.pdb, -merge no (保留两链 moleculetype)
            - peptide ligand 腿 / 小分子: 单链, -merge all 无影响
            Returns (gro_abs, top_abs, gro_H_abs): gro_abs/top_abs 为 prepare-gmx
            产出的 gro/top (工作目录, 供 MakeInputs 读取); gro_H_abs 为加氢后的
            染色体注释用的 gro/pdb (含 -ignh 去掉的 H, 供 PyMOL span 测量;
            对单链它与 gro_abs 等价)。posre.itp(s) + ff link 同步到 builder/。
            """
            gro_p = wdir_p / (Path(pdb_name).stem + '.gro')
            top_p = wdir_p / 'topol.top'
            if gro_p.exists() and top_p.exists():
                put_log(f'{gro_p.name}/topol.top already exist, skip prepare-gmx.', head='ABFE')
            else:
                # pdb2gmx-args: ABFE 强制 -water none (水由后续 solvate 加入);
                # -merge 按腿: 整链(peptide complex) no, 其余 all (prepare-gmx 默认)。
                pdb2gmx_args = self.args.pdb2gmx_args
                if '-water' not in pdb2gmx_args:
                    pdb2gmx_args += ' -water none'
                if f'-merge' not in pdb2gmx_args:
                    pdb2gmx_args += f' -merge {merge_val}'
                argv = ['protein', '-d', str(wdir_p), '-n', pdb_name,
                        '--n-term', 'auto', '--c-term', 'auto',
                        '--pdb2gmx-args', pdb2gmx_args, '--chain-num', '1']
                if self.args.custom_ff_path:
                    argv += ['--ff-dir', str(self.args.custom_ff_path)]
                put_log(f'running prepare-gmx protein in {wdir_p} for {pdb_name} '
                        f'(pdb2gmx args: {pdb2gmx_args})', head='ABFE')
                prepare_gmx_main(argv)
                if not (gro_p.exists() and top_p.exists()):
                    raise RuntimeError(f'prepare-gmx protein did not produce {gro_p}/topol.top in {wdir_p}')
            # posre.itp(s) (pdb2gmx 默认固定名, 与 topol.top 同目录) → builder 目录
            # (gmx_process/MakeInputs 在 builder/wd 里 parmed 读 topol.top, 相对
            # include 只在 top 同目录解析; 整链 per-chain topol_*.itp 也同步)
            builder_dir_p = wdir / 'builder'
            builder_dir_p.mkdir(exist_ok=True, parents=True)
            posres = sorted(wdir_p.glob('posre*.itp'))
            if posres:
                for p in posres:
                    shutil.copy(p, builder_dir_p / p.name)
                put_log(f'copied {[p.name for p in posres]} -> {builder_dir_p}', head='ABFE')
            else:
                put_log(f'no posre*.itp found in {wdir_p}, skip posre copy.', head='ABFE')
            topol_files = sorted(wdir_p.glob('topol_*.itp'))  # per-chain topologies
            if topol_files:
                for p in topol_files:
                    shutil.copy(p, builder_dir_p / p.name)
                put_log(f'copied {[p.name for p in topol_files]} -> {builder_dir_p}', head='ABFE')
            # 力场目录 → builder 符号链接 (topol.top 的 #include "charmm36-jul2022.ff/..."
            # 是相对路径, parmed 预处理器只在 top 所在目录解析, 不看 GMXLIB)
            if self.args.custom_ff_path:
                ff_link = builder_dir_p / Path(self.args.custom_ff_path).name
                if not ff_link.exists():
                    os.symlink(str(self.args.custom_ff_path), str(ff_link))
                put_log(f'linked {self.args.custom_ff_path} -> {ff_link}', head='ABFE')
            return gro_p, top_p

        # complex 腿: 整链 complex.pdb 一次 prepare-gmx protein (-merge no 保留两链
        # moleculetype, 相对坐标天然保留 - 修复 Phase 2 受体/配体独立居中重叠 bug)。
        # 用 complex_path.name (-n/--name 指定的文件名), 与 batch 扫描的输入一致。
        _complex_gro, _complex_top = _run_prepare_gmx_protein(wdir, complex_path.name, merge_val='no')

        gmx = Gromacs(working_dir=str(wdir))

        # build MakeInputs
        # peptide: 整链 gro/top (MakeInputs 读整链, 拆 sys_protein/sys_ligand + rename LIG)
        #          配体腿 = 独立 peptide protein (parmed 读它的 top/gro, 与整链原子数一致
        #          时 _peptide_n_atoms 精确匹配整链中的配体 moleculetype)
        # 小分子: 受体 gro/top + 配体 toff
        protein_def = {'conf': str(_complex_gro), 'top': str(_complex_top),
                       'ff': {'code': self.args.protein_ff}}
        peptide_def = None
        ligand_def = None
        if self.args.peptide:
            # 整链结构: 单独处理, 不走 system_combiner 拼接 (相对坐标已保留)
            ligand_pep_dir = wdir / 'ligand_pep'
            ligand_pep_dir.mkdir(exist_ok=True)
            ligand_pep_def = None
            # prepare-gmx protein 只在 wdir_p 内找 ligand.pdb → 先复制过去
            shutil.copy(ligand_pdb, ligand_pep_dir / 'ligand.pdb')
            _run_prepare_gmx_protein(ligand_pep_dir, 'ligand.pdb', merge_val='all')
            # ligand 腿: 独立 peptide protein (parmed 读单链 top/gro, 原子数 238
            # = 整链中的配体 moleculetype, _peptide_n_atoms 精确匹配)
            ligand_pep_def = {'conf': str(ligand_pep_dir / 'ligand.gro'),
                              'top': str(ligand_pep_dir / 'topol.top'),
                              'ff': {'code': self.args.protein_ff}}
            ligand_def = ligand_pep_def
            peptide_def = {'conf': str(_complex_gro), 'top': str(_complex_top),
                           'ff': {'code': self.args.protein_ff}}
        else:
            # 小分子路径: 入口 (L254) 已抛 NotImplementedError, 此处不可达;
            # 保留显式 raise 作为防御 (避免蛋白整链逻辑误处理小分子 complex.pdb)
            raise NotImplementedError(
                'small-molecule ABFE prepare is NOT implemented yet (complex leg '
                'requires prepare-gmx complex/toff support). '
                'See docs/dev/abfe/small_mol_future_prepare_gmx_complex.md.')

        out_input = wdir / 'input'
        # 有效 rlist (--auto-box 也需要, 不能只依赖 --solv-d-auto 的推导):
        #   peptide 铁律: ligand-leg rlist >= 3.1, 永不降低
        rlist_eff = self.args.fep_rlist_ligand
        if rlist_eff is None:
            if self.args.peptide:
                rlist_eff = 3.1
            else:
                rlist_eff = 1.2
        elif self.args.peptide:
            rlist_eff = max(rlist_eff, 3.1)
        put_log(f'rlist adjust from {self.args.fep_rlist_ligand} to {rlist_eff}', head='ABFE')
        # --solv-d-auto: 按腿推导 d (ligand/complex 分开), 保证半盒矢 >= rlist
        # (互斥已在 main_process 校验, 此处只选分支)
        solv_d_ligand = solv_d_complex = None
        if self.args.auto_box:
            # 轮廓盒模式: Solvate 内部根据 rlist 自适应轴 (每轴 >= 2*(rlist+0.3)),
            # 不再显式传 d; 但 rlist 必须传下去, 否则退化为纯 span 盒
            solv_d_ligand = solv_d_complex = None
            put_log(f'--auto-box: contour box derived from rlist={rlist_eff} '
                    f'inside Solvate.', head='ABFE')
        elif self.args.solv_d_auto:
            # span 测量源: ligand 腿用独立配体 gro (peptide) 或 ligand pdb (小分子),
            # complex 腿用整链 gro (pdb2gmx -ignh 已重建氢, gro 为全原子, 直接测)
            if self.args.peptide:
                ligand_meas = ligand_pep_dir / 'ligand.gro'
            else:
                ligand_meas = ligand_pdb
            complex_meas = _complex_gro
            solv_d_ligand, solv_d_complex = auto_solv_d_pair(
                ligand_meas, complex_meas, bool(self.args.peptide),
                fep_rlist_ligand=rlist_eff)
        elif self.args.solv_d is not None:
            solv_d_ligand = solv_d_complex = self.args.solv_d
        else:
            solv_d_ligand = solv_d_complex = 1.5  # 默认

        builder = system_builder.MakeInputs(
            protein=protein_def,
            host_name='Protein',
            water_model=self.args.water_model,
            custom_ff_path=self.args.custom_ff_path,
            hmr_factor=self.args.hmr_factor,
            solv_d=solv_d_ligand,   # 回退值 (两腿都用)
            solv_d_ligand=solv_d_ligand,
            solv_d_complex=solv_d_complex,
            solv_bt=self.args.solv_bt,
            solv_ion_conc=self.args.ion_conc,
            builder_dir=wdir / self.args.builder_dir,
            gmx=gmx,
            peptide_definition=peptide_def,
            maxwarn=self.args.maxwarn,
            fep_rlist_ligand=rlist_eff,
            auto_box=self.args.auto_box,
            auto_box_pad=self.args.auto_box_pad,
            ligand_chain=self.args.ligand_chain,
            whole_complex=bool(self.args.peptide),  # 整链模式仅 peptide (小分子走 toff)
        )
        with builder:
            builder(ligand_definition=ligand_def, out_dir=out_input)

        put_log(f'ABFE prepare completed. Output in: {out_input}', head='ABFE')


_str2func = {
    'complex': Complex,
}


def main(sys_args: List[str] = None):
    make_args_and_excute('tools for ABFE (absolute binding free energy) preparation.', _str2func, sys_args)


if __name__ == '__main__':
    main()
