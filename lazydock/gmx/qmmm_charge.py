#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QM/MMGBSA 电荷管线: pymol within 选择 + parmed 电荷统计
========================================================
对 gmx_MMPBSA 生成的 _GMXMMPBSA_COM.pdb 做 pymol within 选择(QM区)，
再统计 QM 区净电荷，输出 gmx_MMPBSA 输入文件的
qmcharge_rec / qmcharge_lig / qmcharge_com。

=====================================================================
残基编号的坑 (为什么不用残基号直查 prmtop)
=====================================================================
_GMXMMPBSA_COM.pdb 是多链复合物, 每条链残基号各自从 1 开始(链内号);
COM.prmtop 是 Amber 全局连续编号(从 0 起)。
   例 (rbb_d7): 链A 残基 1..443, 链B 残基 1..18
               prmtop: A 残基 0..442, B 残基 443..460 (偏移 442)
偏移 = 配体前面所有链的残基数之和 - 1, 随体系变化(四链更大), 不可 hardcode。
若拿 pymol 的"链内号"直查 prmtop.residues[i], 会取到受体残基, 电荷全错。

本脚本的稳健映射(不依赖偏移、不依赖 pymol 原子 id):
  1. pymol 只取"选中残基身份" (chain, resi, resn) —— 不用 atom.id
     (pymol 的 atom.id 不是文件序, 是内部重排索引, 不能用于映射)
  2. 解析 PDB 文件, 构建 (chain, resi) -> [serial...] (serial 连续 1..N)
  3. serial -> prmtop 原子索引: prm_idx = serial - 1
  4. 对选中残基的原子, 取 prmtop 残基对象, 残基去重后电荷求和

与 gmx_MMPBSA 一致性
--------------------
make_top.py qm_residues="within X" 分支:
  - 从 complex_trajs[0] 第一帧(dump 0) 生成 _GMXMMPBSA_COM.pdb
  - 受体残基 r 与配体残基 l 任意原子对距离 ≤ X -> 双侧残基全选
  - 全原子(含H)欧氏距离, 无 PBC 镜像
pymol 等价:
  - byres(chain A within X of chain B) + byres(chain B within X of chain A)
  - within = S1 中距 S2 任意原子 ≤ X 的原子(含H, 原子中心距离)
  - byres 展开为完整残基

qc 说明: 输出的整数满足 gmx_MMPBSA main.py 硬检查
  qmcharge_lig + qmcharge_rec == qmcharge_com
qmcharge_com 先按 rec+lig 求和取整(物理上二者无重叠)。
"""
import re
from collections import defaultdict

import parmed as pmd


def parse_pdb_residue_serials(pdb: str):
    """(chain, resi) -> [PDB serial...]. 要求 serial 连续 1..N, 否则报错."""
    res_serials = defaultdict(list)
    serials = []
    n_atom = 0
    with open(pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                n_atom += 1
                ser = int(line[6:11])
                serials.append(ser)
                key = (line[21], int(line[22:26].strip()))
                res_serials[key].append(ser)
    if serials != list(range(1, n_atom + 1)):
        raise ValueError(
            f'PDB serial 不连续(1..{n_atom}), 无法做 serial->prmtop 映射。'
            '请确认 pdb 未被重新编号/筛选。')
    return res_serials


def qm_charges(prmtop: str, res_serials, rec_keys, lig_keys):
    """统计 rec/lig 残基的净电荷(残基对象去重, round 取整)."""
    prm = pmd.load_file(prmtop)
    n_prm = len(prm.atoms)

    def charge_of(keys):
        q, seen, nres = 0.0, set(), 0
        for ch, resi, resn in keys: # type: ignore
            for ser in res_serials[(ch, resi)]:
                idx = ser - 1
                if not (0 <= idx < n_prm):
                    raise ValueError(f'serial {ser} 超出 prmtop 原子范围 1..{n_prm}')
                r = prm.atoms[idx].residue
                if r not in seen:
                    seen.add(r)
                    nres += 1
                    q += sum(a.charge for a in r.atoms)
        return q, nres

    q_rec, nr = charge_of(rec_keys)
    q_lig, nl = charge_of(lig_keys)
    # 先各自 round 再求和, 保证整数精确满足 gmx_MMPBSA 硬检查:
    #   qmcharge_com == qmcharge_rec + qmcharge_lig
    # (不能 round(q_rec+q_lig): 个别情况 round(rec)+round(lig) != round(rec+lig)
    #  会让硬检查失败, 例如 rec=0.6, lig=0.6)
    r, l = round(q_rec), round(q_lig)
    return r, l, r + l, nr, nl


def update_input_file(in_path: str, q_rec: int, q_lig: int, q_com: int):
    """把 qmcharge_* 写回 gmx_MMPBSA 输入文件的 &gb 段(就地更新, 保留注释)."""
    with open(in_path) as f:
        lines = f.readlines()
    values = {'qmcharge_com': q_com, 'qmcharge_rec': q_rec, 'qmcharge_lig': q_lig}
    hit = set()
    out = []
    for line in lines:
        stripped = line.split('#')[0].strip()          # 去注释
        matched = False
        for key in values:
            m = re.match(r'^\s*' + key + r'\s*=', stripped)
            if m:
                indent = line[:len(line) - len(line.lstrip())]
                out.append(f'{indent}{key} = {values[key]}      # [auto] pymol+parmed 计算\n')
                hit.add(key)
                matched = True
                break
        if not matched:
            out.append(line)
    miss = set(values) - hit
    if miss:
        raise SystemExit(f'输入文件 {in_path} 缺少键 {sorted(miss)}，请在 &gb 段手工添加。')
    with open(in_path, 'w') as f:
        f.writelines(out)
    return sorted(hit)