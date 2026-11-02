# -*- coding: utf-8 -*-
"""外部非标准氨基酸(ncAA)模板的运行时加载.

背景
----
除主库 storage/opmm/charmm36_nc_aa.xml 中已开源的少量模板外, 项目实践中常遇到
新的非标准氨基酸. 若每加一个残基都要改代码(残基名白名单/模板文件列表/FF 加载
列表), 维护成本高.

本模块提供一套运行时加载机制: 用户把 OpenMM FF XML(残基模板 + 缺失键角参数)与
同名 YAML(原子名规范化等元数据)放在任意目录, 通过 CLI 参数或构造参数传入 XML
路径即可支持新残基, 无需修改代码.

配对约定
--------
XML 与同名 YAML 成对使用(如 my_nc_aa.xml 与 my_nc_aa.yml). 调用方只传 XML,
YAML 按同名规则自动发现; 缺失 YAML 时元数据为空, 不影响纯模板情形.

YAML 字段(均可选)
-----------------
atom_rename:    {残基名: {上游原子名: 模板原子名}}
resname_alias:  {全名: [别名, ...]}  留空则按前 3 / 前 4 字符自动派生
"""
from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Set

import yaml
from mbapy_lite.base import put_log

# 主库内建的 ncAA 模板(已开源, 随包分发). 外部模板无需包含这些残基.
BUILTIN_TEMPLATE_FILES = ['opmm/charmm36_nc_aa.xml', 'opmm/charmm36_caps.xml']


def _parse_residue_names(xml_path: str) -> List[str]:
    """从 OpenMM FF XML 中取出所有 <Residue name="..."> 的残基名."""
    import xml.etree.ElementTree as ET
    try:
        root = ET.parse(xml_path).getroot()
    except Exception as e:
        put_log(f'custom_residues: failed to parse {xml_path}: {e}', bg_color='red')
        return []
    res_els = root.find('Residues')
    if res_els is None:
        return []
    return [el.get('name') for el in res_els.findall('Residue') if el.get('name')]  # pyright: ignore[reportReturnType]


def _find_meta(xml_path: str) -> Optional[str]:
    """XML 同目录下查找同名 YAML(.yml / .yaml)."""
    stem = os.path.splitext(xml_path)[0]
    for ext in ('.yml', '.yaml'):
        candidate = stem + ext
        if os.path.exists(candidate):
            return candidate
    return None


class CustomResidues:
    """运行时自定义残基模板集.

    参数
    ----
    xml_paths:
        外部模板 XML 路径列表. 也接受目录(自动收集目录下的 *.xml)或字符串.
    meta_paths:
        显式指定 YAML 路径. 一般留空, 由 XML 同名规则自动发现.
    """

    def __init__(self, xml_paths: Optional[Iterable[str]] = None,
                 meta_paths: Optional[Iterable[str]] = None):
        self._xml_files: List[str] = []
        self._residue_names: Set[str] = set()
        self._atom_rename: Dict[str, Dict[str, str]] = {}
        self._aliases: Dict[str, str] = {}  # 截断名 -> 全名

        xml_paths = self._normalize(xml_paths)
        meta_paths = self._normalize(meta_paths)

        for path in xml_paths:
            if not os.path.exists(path):
                put_log(f'custom_residues: template file not found, skipped: {path}', bg_color='red')
                continue
            self._xml_files.append(path)
            names = _parse_residue_names(path)
            self._residue_names.update(names)
            for name in names:
                self._auto_alias(name)

        # 元数据: 显式指定的优先, 未指定则按同名自动发现
        metas = meta_paths or [m for m in (_find_meta(p) for p in self._xml_files) if m]
        for meta in metas:
            self._load_meta(meta)

    @staticmethod
    def _normalize(paths) -> List[str]:
        """把 str / 单个路径 / 目录统一成"文件路径列表"."""
        if not paths:
            return []
        if isinstance(paths, str):
            paths = [paths]
        out: List[str] = []
        for p in paths:
            p = os.path.abspath(os.path.expanduser(str(p)))
            if os.path.isdir(p):
                out.extend(sorted(os.path.join(p, f) for f in os.listdir(p)
                                  if f.lower().endswith('.xml')))
            else:
                out.append(p)
        return out

    def _auto_alias(self, name: str):
        """残基名超过 PDB 3 字符字段宽度时, 自动派生截断别名.

        两种写出端的截断长度不一致: OpenMM PDBFile.writeFile 取前 3 字符,
        PyMOL cmd.save 保留前 4 字符. 两者都要能命中模板, 否则 PDBFixer
        无法为该残基补氢.
        """
        if len(name) <= 3:
            return
        for alias in {name[:3], name[:4]}:
            if alias != name:
                self._aliases.setdefault(alias, name)

    def _load_meta(self, meta_path: str):
        try:
            # YAML 允许中文注释, 必须显式用 UTF-8 读取; 否则会退回到系统默认
            # 编码(某些环境下是 ASCII), 遇到中文注释直接解码失败.
            with open(meta_path, encoding='utf-8') as f:
                data = yaml.safe_load(f) or {}
        except Exception as e:
            put_log(f'custom_residues: failed to parse {meta_path}: {e}', bg_color='red')
            return
        if not isinstance(data, dict):
            return
        for resname, mapping in (data.get('atom_rename') or {}).items():
            if mapping:
                self._atom_rename.setdefault(resname, {}).update(
                    {str(k): str(v) for k, v in mapping.items()})
        for full, aliases in (data.get('resname_alias') or {}).items():
            if not aliases:
                continue
            if isinstance(aliases, str):
                aliases = [aliases]
            for alias in aliases:
                self._aliases[str(alias)] = str(full)

    # ------------------------------------------------------------ 只读属性

    @property
    def xml_files(self) -> List[str]:
        """供 ForceField(*files) 加载的模板文件列表."""
        return list(self._xml_files)

    @property
    def residue_names(self) -> Set[str]:
        """模板中定义的残基名(含自动派生的截断别名)."""
        return set(self._residue_names) | set(self._aliases)

    @property
    def atom_rename(self) -> Dict[str, Dict[str, str]]:
        """{残基名: {上游原子名: 模板原子名}}."""
        return {k: dict(v) for k, v in self._atom_rename.items()}

    @property
    def aliases(self) -> Dict[str, str]:
        """{截断名: 全名}."""
        return dict(self._aliases)

    def alias_of(self, resname: str) -> List[str]:
        """返回某全名对应的所有截断别名."""
        return [a for a, full in self._aliases.items() if full == resname]

    def __bool__(self) -> bool:
        return bool(self._xml_files)

    def __len__(self) -> int:
        return len(self._residue_names)

    def __repr__(self) -> str:
        return (f'CustomResidues(files={len(self._xml_files)}, '
                f'residues={len(self._residue_names)})')