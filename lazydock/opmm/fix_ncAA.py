# -*- coding: utf-8 -*-
"""
非标准氨基酸 (ncAA) 支持: 原子名规范化 + PDBFixer 模板注册.

背景: 输入 PDB 中修饰氨基酸(如 LE1 青霉胺)的侧链原子名可能与
CHARMM36/GROMACS 官方模板不一致(如 PDB: SG/C8/C9 vs CHARMM36: SG3/CG1/CG2),
导致 OpenMM/PDBFixer/GROMACS 无法按模板匹配原子. 本模块提供:

  1. fix_ncAA_pdb: 把非标准残基原子名规范化为 CHARMM36/GROMACS 约定(fix-cp 写文件前调用).
  2. register_ncAA_templates: 把 storage/opmm 下 ncAA XML 模板注册到 PDBFixer,
     使其能为这些残基补氢(LE1 无 L 等价物可借, 否则 OpenMM createSystem
     会因拓扑缺氢/模板含氢而报错). 需在 addMissingHydrogens 之前调用.
"""

from __future__ import annotations

import os
from typing import Dict, Optional

from mbapy_lite.base import put_log
from openmm import app as openmm_app
from openmm import unit
from openmm.vec3 import Vec3

from lazydock.utils import get_storage_path

# 封端残基(ACE/NH2): 原子名与模板一致, 无规范化需求, 但需注册模板补氢.
# ACE 在 PDBFixer 内置 _standardTemplates 中(内置 ACE.pdb 无氢, 会补不上),
# 注册时需把它移出 _standardTemplates 以强制走注册模板分支.
CAPS_RESIDUES = {'ACE', 'NH2'}

# N-甲基氨基酸(SAR/MLE): PDBFixer 替换表会把 SAR->GLY、MLE->LEU(静默去甲基),
# 与 D-aa 同款问题; 但这里模板自带 CN 甲基, 需先"替换补氢再改回名"(与 D-aa 同款).
# 原子名与模板一致(CN/CA/...), 无规范化需求.
NCME_RESIDUES = {'SAR', 'MLE'}

# 主库内建模板文件(storage/opmm 相对路径). 其余 ncAA 模板由调用方通过
# custom_residue_xml 参数在运行时加载, 见 lazydock.opmm.custom_residues.
NCAA_TEMPLATE_FILES = ['opmm/charmm36_nc_aa.xml', 'opmm/charmm36_caps.xml']

# 内建残基的原子名规范化(上游 PDB 原子名 -> CHARMM36/GROMACS 原子名).
# 只有随主库分发的残基才写在这里; 外部模板的同名映射请放在配对 YAML 的
# atom_rename 字段, 由 CustomResidues 提供.
BUILTIN_ATOM_RENAME = {
    'LE1': {'SG': 'SG3', 'C8': 'CG1', 'C9': 'CG2'},  # Penicillamine 青霉胺
}


def _collect_templates(custom: Optional['CustomResidues'] = None):
    """汇总内建模板与运行时外部模板的 (文件列表, 残基名集合, 原子名映射, 截断别名).

    custom 为 None 或空时只返回内建集, 行为与改造前一致.
    """
    files = [get_storage_path(rel) for rel in NCAA_TEMPLATE_FILES]
    resnames = set(BUILTIN_ATOM_RENAME) | set(CAPS_RESIDUES) | set(NCME_RESIDUES)
    atom_rename: Dict[str, Dict[str, str]] = {
        k: dict(v) for k, v in BUILTIN_ATOM_RENAME.items()}
    aliases: Dict[str, str] = {}
    if custom:
        files = files + custom.xml_files
        resnames |= custom.residue_names
        for resname, mapping in custom.atom_rename.items():
            atom_rename.setdefault(resname, {}).update(mapping)
        aliases = custom.aliases
    return files, resnames, atom_rename, aliases


def fix_ncAA_pdb(lines, pdb_path=None, custom=None):
    """把 ATOM/HETATM 行中非标准残基的原子名规范化为 CHARMM36/GROMACS 约定.
    仅改原子名(13-16列), 坐标/元素/CONECT(按 serial)不受影响.
    custom: CustomResidues 实例, 提供外部模板的原子名映射; None 时只用内建表.
    返回新行列表; pdb_path 仅用于日志."""
    _, _, atom_rename, _ = _collect_templates(custom)
    n_renamed = 0
    out = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')) and line[17:20].strip() in atom_rename:
            ren = atom_rename[line[17:20].strip()]
            old = line[12:16].strip()
            new = ren.get(old)
            if new:
                line = line[:12] + f'{new:>4}' + line[16:]
                n_renamed += 1
        out.append(line)
    if n_renamed:
        put_log(f'fix_ncAA: renamed {n_renamed} ncAA atom name(s) in '
                f'{pdb_path.name if pdb_path else "pdb"}', bg_color='blue')
    return out


def _type_to_element(type_):
    """按 CHARMM36 原子类型首字母推断元素(类型首字母即元素, 除 Cl 外)."""
    t = (type_ or '').upper()
    if t.startswith('CL'):
        return openmm_app.element.chlorine
    return {'H': openmm_app.element.hydrogen,
            'C': openmm_app.element.carbon,
            'N': openmm_app.element.nitrogen,
            'O': openmm_app.element.oxygen,
            'S': openmm_app.element.sulfur,
            'P': openmm_app.element.phosphorus,
            'F': openmm_app.element.fluorine,
            'I': openmm_app.element.iodine,
            'BR': openmm_app.element.bromine}[t[0]]


def _parse_ncAA_xml(xml_path):
    """解析 ncAA XML 模板文件, 返回 {残基名: (topology, positions, terminal)}.
    positions 用占位坐标(全 0): PDBFixer 补氢时氢坐标由父原子位置推导并做能量优化,
    模板坐标不被使用. 模板原子(含氢)全部视为非 terminal, PDBFixer 会补全部缺失氢;
    参与二硫键时的巯基氢(HG3)删除由 relax 阶段按拓扑键判断处理."""
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_path)
    out = {}
    res_els = tree.getroot().find('Residues')
    if res_els is None:
        return out
    for res_el in res_els.findall('Residue'):
        name = res_el.get('name')
        atoms = [(a_el.get('name'), _type_to_element(a_el.get('type')))
                 for a_el in res_el.findall('Atom')]
        bonds = [(b_el.get('atomName1'), b_el.get('atomName2'))
                 for b_el in res_el.findall('Bond')]
        topology = openmm_app.Topology()
        chain = topology.addChain()
        residue = topology.addResidue(name, chain)
        atoms_by_name = {aname: topology.addAtom(aname, elem, residue)
                         for aname, elem in atoms}
        for a1, a2 in bonds:
            if a1 in atoms_by_name and a2 in atoms_by_name:
                topology.addBond(atoms_by_name[a1], atoms_by_name[a2])
        positions = unit.Quantity([Vec3(0, 0, 0)] * len(atoms), unit.nanometer)
        terminal = [False] * len(atoms)
        out[name] = (topology, positions, terminal)
    return out


def register_ncAA_templates(fixer, resnames=None, custom=None):
    """把 ncAA/封端模板(带氢)注册到 PDBFixer, 使其能为这些残基补氢.
    resnames: 只注册给定残基名(建议传拓扑中出现过的), None 则全部注册.
    custom: CustomResidues 实例, 提供运行时加载的外部模板; None 时只用内建集.
    必须在 fixer.addMissingHydrogens() 之前调用. 返回注册模板数."""
    files, builtin_resnames, _, aliases = _collect_templates(custom)
    if resnames is None:
        resnames = set(builtin_resnames)
    # 把被 PDB 4 字符限制截断的残基名一并视为待注册(如 A1A -> A1A2Y),
    # 使模板同时以截断名注册一份.
    resnames = set(resnames) | {trunc for trunc, full in aliases.items()
                                if full in resnames}
    n_registered = 0
    for path in files:
        for name, (topology, positions, terminal) in _parse_ncAA_xml(path).items():
            # 同时按原名与截断名注册(A1A2Y 与 A1A 指向同一模板)
            alias_list = [name] + [trunc for trunc, full in aliases.items() if full == name]
            for alias in alias_list:
                if alias not in resnames:
                    continue
                # ACE 在 PDBFixer 内置模板中(ACE.pdb 无氢), 用它补氢会失败;
                # 移出 _standardTemplates, 强制 _describeVariant 走注册模板分支(带氢).
                fixer._standardTemplates.discard(alias)
                fixer.registerTemplate(_rename_template(topology, alias), positions,
                                       terminal=terminal)
                n_registered += 1
    return n_registered


def _rename_template(topology, new_name):
    """复制一份模板拓扑并改残基名(用于截断名别名, 避免同名冲突)."""
    out = openmm_app.Topology()
    chain = out.addChain()
    old_res = next(iter(topology.residues()))
    residue = out.addResidue(new_name, chain, old_res.id)
    atoms = {a.name: out.addAtom(a.name, a.element, residue) for a in old_res.atoms()}
    for a1, a2 in topology.bonds():
        out.addBond(atoms[a1.name], atoms[a2.name])
    return out


def _load_template_bonds(custom=None):
    """解析 ncAA XML, 返回 {残基名: [(原子名1, 原子名2), ...]}, 含截断名别名."""
    import xml.etree.ElementTree as ET
    files, _, _, aliases = _collect_templates(custom)
    templates = {}
    for path in files:
        if not path or not os.path.exists(path):
            continue
        res_els = ET.parse(path).getroot().find('Residues')
        if res_els is None:
            continue
        for rel_el in res_els.findall('Residue'):
            templates[rel_el.get('name')] = [
                (b.get('atomName1'), b.get('atomName2')) for b in rel_el.findall('Bond')]
    # 截断名反向映射(A1A / A1A2 -> A1A2Y)
    for trunc, full in aliases.items():
        if full in templates:
            templates[trunc] = templates[full]
    return templates


def sew_ncAA_internal_bonds(topology, positions, custom=None):
    """按 ncAA 模板补全非标准残基缺失的"内部重原子键".

    背景: PDBFixer 重建拓扑时, 残基内部键只有两条来源——
      1) createStandardBonds() 仅覆盖 OpenMM Topology._standardBonds 中的标准残基;
      2) 从原拓扑继承已存在的键(且要求至少一端是非标准残基).
    因此非标准残基的内部键完全依赖上游文件里的 CONECT 记录. 但上游并不总是可靠:
    某些 CIF 直接转 PDB 时没有 CONECT; PyMOL 转换对以数字开头的残基名(如 9YT)
    也可能整条推断失败. 缺键时 OpenMM createSystem 会报
    "<|place_holder_mm_span_0442|> set of atoms matches <RES>, but the bonds are different".

    本函数在 _fix 之后、建 system 之前调用: 逐残基比对模板键, 只补"两端原子都在
    该残基内且尚未成键"的组合, 不重复添加, 也不触碰跨残基键. 对已带完整 CONECT
    的输入等价于空操作, 因此可作为标准步骤常开.
    返回 (topology, positions, 补键数)."""
    templates = _load_template_bonds(custom)
    if not templates:
        return topology, positions, 0

    modeller = openmm_app.Modeller(topology, positions)
    existing = set()
    for a1, a2 in modeller.topology.bonds():
        existing.add((a1.index, a2.index))
        existing.add((a2.index, a1.index))
    n_added = 0
    for residue in modeller.topology.residues():
        tmpl_bonds = templates.get(residue.name)
        if not tmpl_bonds:
            continue
        atoms = {a.name: a for a in residue.atoms()}
        for n1, n2 in tmpl_bonds:
            a1, a2 = atoms.get(n1), atoms.get(n2)
            if a1 is None or a2 is None:
                continue
            if (a1.index, a2.index) in existing:
                continue
            modeller.topology.addBond(a1, a2)
            existing.add((a1.index, a2.index))
            existing.add((a2.index, a1.index))
            n_added += 1
    return modeller.topology, modeller.positions, n_added