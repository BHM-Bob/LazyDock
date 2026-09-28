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

from mbapy_lite.base import put_log
from openmm import app as openmm_app
from openmm import unit
from openmm.vec3 import Vec3

from lazydock.utils import get_storage_path

# 非标准氨基酸原子名规范化: 残基名 -> {PDB原子名: CHARMM36/GROMACS 原子名}
NCAA_ATOM_RENAME = {
    'LE1': {'SG': 'SG3', 'C8': 'CG1', 'C9': 'CG2'},  # Penicillamine 青霉胺
}
NCAA_RESIDUES = set(NCAA_ATOM_RENAME)

# ncAA OpenMM XML 模板文件(storage/opmm 相对路径), 与 relax.py 加载列表保持一致
NCAA_TEMPLATE_FILES = ['opmm/charmm36_nc_aa.xml']


def fix_ncAA_pdb(lines, pdb_path=None):
    """把 ATOM/HETATM 行中非标准残基的原子名规范化为 CHARMM36/GROMACS 约定.
    仅改原子名(13-16列), 坐标/元素/CONECT(按 serial)不受影响.
    返回新行列表; pdb_path 仅用于日志."""
    n_renamed = 0
    out = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')) and line[17:20] in NCAA_ATOM_RENAME:
            ren = NCAA_ATOM_RENAME[line[17:20]]
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


def register_ncAA_templates(fixer, resnames=None):
    """把 storage/opmm 的 ncAA 模板(带氢)注册到 PDBFixer, 使其能为 ncAA 补氢.
    resnames: 只注册给定残基名(建议传拓扑中出现过的), None 则全部注册.
    必须在 fixer.addMissingHydrogens() 之前调用. 返回注册模板数."""
    if resnames is None:
        resnames = set(NCAA_RESIDUES)
    n_registered = 0
    for rel in NCAA_TEMPLATE_FILES:
        for name, (topology, positions, terminal) in _parse_ncAA_xml(get_storage_path(rel)).items():
            if name not in resnames or name in fixer._standardTemplates:
                continue
            fixer.registerTemplate(topology, positions, terminal=terminal)
            n_registered += 1
    return n_registered