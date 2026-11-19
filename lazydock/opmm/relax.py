# most of the code is from PepGLAD
import io
import os
from types import SimpleNamespace
from typing import List, Optional, Union

import openmm
import pdbfixer
from mbapy_lite.base import put_log
from openmm import app as openmm_app
from openmm import unit

from lazydock.opmm.custom_residues import CustomResidues
from lazydock.opmm.fix_ncAA import (register_ncAA_templates,
                                   sew_ncAA_internal_bonds)
from lazydock.utils import get_storage_path

# 二硫键可能出现在标准 CYS(SG) 或含硫修饰氨基酸(如 LE1 青霉胺 SG3)上, 统一识别:
#   - 缝接/保护二硫键时按原子名匹配
#   - 注意: 输入 PDB 中 LE1 的硫原子原始名为 SG, 会在 fix-cp 阶段被
#     fix_ncAA 规范化为 SG3, 故两处均需兼容
SULFUR_ATOM_NAMES = ('SG', 'SG3')

# 巯基氢原子名: 标准 CYS 为 HG, LE1(青霉胺) 为 HG3, 其余含硫修饰残基沿用 CHARMM 的 HG1.
SULFUR_HYDROGEN_NAMES = ('HG', 'HG1', 'HG3')

ENERGY = unit.kilojoules_per_mole
LENGTH = unit.nanometer  # pyright: ignore[reportAttributeAccessIssue]

# backbone 闭合(头尾肽键)保护: CustomBondForce 的平衡键长与力常数.
# 头尾肽键是标准酰胺 C-N, CHARMM36 无对应 patch 参数, 故用 CustomBondForce 兜底.
BACKBONE_BOND_LEN = 0.134  # nm, 酰胺 C-N 平衡键长 (~1.34 A)
BACKBONE_BOND_STIFFNESS = 3e5  # kJ/mol/nm^2

# 常见离子元素 -> 离子力场模板名(charmm36.xml 不含离子, 需额外加载 ions.xml)
from openmm.app.element import (calcium, chlorine, magnesium, potassium,
                                sodium, zinc)

ION_ELEMENT_TO_TEMPLATE = {
    calcium: 'CAL', magnesium: 'MG', zinc: 'ZN',
    sodium: 'SOD', potassium: 'POT', chlorine: 'CLA',
}
ION_TEMPLATE_NAMES = set(ION_ELEMENT_TO_TEMPLATE.values())

# D-氨基酸: PDB 三字码(DAL/DPN/...) -> GROMACS 约定名(DALA/DPHE/...), 需要时用于改名.
# PDBFixer 的 substitutions 表会把 D* 替换成 L* (DAL->ALA 等), 从而丢失 D 型手性;
# 我们通过"先记录位置 -> 允许补H(替换成L) -> 再把残基名改回 D"来保留 D 型.
# 映射: PDB D名 -> pdbfixer 替换目标(L名) -> GROMACS D名(实测自 charmm36-jul2022 aminoacids.rtp)
D_AMINO_ACID_PDB_TO_L = {
    # 实测自 pdbfixer.substitutions
    'DAH': 'PHE', 'DAL': 'ALA', 'DAR': 'ARG', 'DAS': 'ASP',
    'DCY': 'CYS', 'DGL': 'GLU', 'DGN': 'GLN', 'DHA': 'ALA',
    'DHI': 'HIS', 'DIL': 'ILE', 'DIV': 'VAL', 'DLE': 'LEU',
    'DLY': 'LYS', 'DNP': 'ALA', 'DPN': 'PHE', 'DPR': 'PRO',
    'DSN': 'SER', 'DSP': 'ASP', 'DTH': 'THR', 'DTR': 'TRP',
    'DTY': 'TYR', 'DVA': 'VAL',
}
# L名(按 CHARMM 组氨酸三码 HSD/HSE/HSP 处理) -> GROMACS D名(实测自 aminoacids.rtp)
D_AMINO_ACID_L_TO_GMX = {
    'ALA': 'DALA', 'ARG': 'DARG', 'ASN': 'DASN', 'ASP': 'DASP',
    'CYS': 'DCYS', 'GLN': 'DGLN', 'GLU': 'DGLU', 'GLY': 'DGLY',
    'HSD': 'DHSD', 'HSE': 'DHSE', 'HSP': 'DHSP', 'ILE': 'DILE',
    'LEU': 'DLEU', 'LYS': 'DLYS', 'MET': 'DMET', 'PHE': 'DPHE',
    'PRO': 'DPRO', 'SER': 'DSER', 'THR': 'DTHR', 'TRP': 'DTRP',
    'TYR': 'DTYR', 'VAL': 'DVAL',
}


class ForceFieldMinimizer(object):
    """PDBFixer + OpenMM 结构弛豫器(CHARMM36 力场).

    典型用法::

        minimizer = ForceFieldMinimizer(platform='CUDA',
                                        backbone_chains=['P'],
                                        disulfide_chains=['P'])
        relaxed_pdb, info = minimizer(pdb_path, out_path,
                                      restrain_chain=None, restrain_backbone=False,
                                      return_info=True)

    流程: ``_fix``(PDBFixer 补缺失原子/氢, 并保真 D-氨基酸名) -> ``_sew_caps_bonds``
    (C 端 NH2 酰胺封端) -> 按需缝接 backbone 闭合键/二硫键 -> ``createSystem``
    (恒 ``ignoreExternalBonds=False``, 以便 OpenMM 应用 DISU 等 patch) ->
    LocalEnergyMinimizer 最小化 -> 写回 PDB(补回 writeFooter 丢弃的跨残基键).

    特殊拓扑的处理职责划分:
      - D-氨基酸: PDBFixer 会将其静默替换为 L 型, 这里改回 D 名(坐标不变, 手性靠坐标保持);
      - 非标准氨基酸: 由 ``custom_residues`` 提供的模板注册给 PDBFixer 补氢, 并从模板补内部重原子键;
      - 离子: ``charmm36.xml`` 无离子模板, 探测到时联合加载 ``ions.xml``;
      - backbone 闭合(头尾肽键): 缝接首 N-尾 C 键, 并用 ``CustomBondForce`` 保护,
        因为 CHARMM36 没有对应的 patch 参数;
      - 二硫键: 分标准与非标两条路径, 但都不加额外键势(键势均来自力场 patch 参数):
          * 标准 CYS(SG-SG): 由 PDBFixer 的 ``createDisulfideBonds`` 按 **硬编码 0.3 nm**
            距离判定并补拓扑键, createSystem 时匹配 charmm36.xml 的 DISU patch
            (SG->SM, CT1->CT2), 直接生效.
            注意: PDBFixer 重建拓扑时会**丢弃**输入 CONECT 中两端均为标准残基的键,
            因此 fix-cp 即便按更大的 ``--ss-max-dist`` 把 0.32 nm 的 SG-SG 写进 CONECT,
            该键仍会在 ``_fix`` 中被丢弃(实测 0.29 nm 保留 / 0.32 nm 丢弃, 与 CONECT 有无无关).
            若需要保留这种超出 0.3 nm 的 S-S, 由本类在 PDBFixer **之后**用
            ``_sew_disulfide_bonds`` 按 ``self.ss_max_dist`` 重新补键(该参数默认 0.3 nm);
          * 非标含硫残基(如 LE1 的 SG3)与 CYS 成键: OpenMM 不认识这种组合, 由本类分三步处理
            -- ``_sew_disulfide_bonds`` 按 ``SULFUR_ATOM_NAMES`` 与 ss_max_dist 缝接 SG3-SG 拓扑键,
            ``_remove_disulfide_hydrogens`` 删掉两侧巯基氢(``SULFUR_HYDROGEN_NAMES``: HG/HG1/HG3)
            使其呈 CYX 式无氢状态, 最后靠 ``charmm36_nc_aa.xml`` 的 LE1SS patch 把 SG3/SG 都
            改为 SM 原子类型并补 CG301-SM 等键参数(SM-SM 二硫键参数沿用官方 CHARMM36).
    """

    def __init__(self, stiffness=10.0, max_iterations=0, tolerance=10, platform='CUDA', constraints=openmm_app.HBonds,
                 backbone_chains: Optional[List[str]] = None,
                 disulfide_chains: Optional[List[str]] = None,
                 ss_max_dist: float = 0.3,
                 device_index: Optional[int] = None,
                 custom_residue_xml: Optional[Union[str, List[str]]] = None):
        """初始化弛豫器.

        Args:
            stiffness (float): 位置约束(CustomExternalForce)的力常数, 单位 kJ/mol/nm^2,
                仅在使用 ``restrain_chain`` 时生效. 默认 10.
            max_iterations (int): 最小化最大迭代步数, 0 表示由 OpenMM 自行决定(收敛或达到内部上限).
                默认 0.
            tolerance (float): 最小化收敛判据, 单位 kJ/mol/nm. 默认 10.
            platform (str): OpenMM 计算平台, ``'CUDA'`` 或 ``'CPU'``. 默认 ``'CUDA'``.
            constraints (openmm.app.Constraint or None): 键长约束方式,
                ``HBonds`` 约束含氢键长, ``AllBonds`` 全约束, ``None`` 不约束. 默认 ``HBonds``.
            backbone_chains (List[str] | None): 需要保护 backbone 闭合(头尾肽键)的链 ID 列表.
                对这些链缝接首 N - 尾 C 肽键, 并加 ``CustomBondForce`` 防止最小化时断链.
                默认 None(不处理).
            disulfide_chains (List[str] | None): 含二硫键、需要保护的链 ID 列表.
                标准 CYS 的 S-S 由 OpenMM 的 DISU patch 处理; 非标含硫残基(如 LE1 的 SG3)
                由本类手工缝接二硫键并删除巯基氢. 默认 None(不处理).
            ss_max_dist (float): SG-SG 距离阈值(单位 nm), 阈值内且尚未成键的残基对判定为二硫键并补键.
                默认 0.3(即 3.0 A).
            device_index (int | None): CUDA 设备序号, 仅 ``platform='CUDA'`` 时生效.
                None 表示不设置 ``DeviceIndex``(使用 OpenMM 默认卡). 默认 None.
            custom_residue_xml (str | List[str] | None): 运行时外部非标准氨基酸模板的 XML 路径
                (任意路径), 需与同名 YAML 配对使用. 用于在不修改代码的前提下支持新的 ncAA.
                None 或空表示只用内建模板集. 默认 None.
        """
        super().__init__()
        self.stiffness = stiffness * ENERGY/(LENGTH ** 2) # coefficient for restraints force, unit: kJ*mol^-1*nm^-2  # pyright: ignore[reportOperatorIssue]
        self.max_iterations = max_iterations
        self.tolerance = tolerance * ENERGY/LENGTH  # pyright: ignore[reportOperatorIssue]
        assert platform in ('CUDA', 'CPU')
        self.platform = platform
        self.constraints = constraints
        # backbone 闭合(头尾肽键)支持: 对 backbone_chains 中的链, 在 PDBFixer 补全后缝接首尾肽键,
        # 并加 CustomBondForce 保护该键, 避免最小化时断链
        self.backbone_chains = list(backbone_chains or [])
        # 二硫键支持: 对 disulfide_chains 中的链, 补全 SG-SG 拓扑键.
        # createSystem 恒以 ignoreExternalBonds=False 建系, OpenMM 会自动应用
        # CHARMM36 DISU patch(S-S 平衡键长 0.2029nm, k=144766 kJ/mol/nm^2),
        # 无需手动加键势.
        self.disulfide_chains = list(disulfide_chains or [])
        self.ss_max_dist = ss_max_dist  # nm, SG-SG 距离阈值内判定为二硫键并补键(3.0A)
        # GPU 设备选择: 仅 CUDA 时生效. 向前兼容: None 表示用默认卡(不设置 DeviceIndex)
        self.device_index = device_index
        # 运行时外部非标准氨基酸模板(任意路径的 XML, 与同名 YAML 配对).
        # 用于在不改代码的前提下支持新的 ncAA; None/空 时行为与内建集完全一致.
        self.custom_residues = CustomResidues(custom_residue_xml)

    def _fix(self, pdb_str):
        # D-氨基酸处理: PDB 三字码(DAL/...) 会被 pdbfixer 替换成 L 型, 丢失 D 手性.
        # 策略: 记录 D 残基位置 -> 正常替换+补H(用 L 模板补 H) -> 再把残基名改回 D 名,
        # 以保留 D 型标记(坐标本就不动, 手性靠坐标保持).
        d_res_pos = {}  # (chain_id, resSeq) -> (L名, GROMACS D名)
        for ln in pdb_str.splitlines():
            if ln.startswith(('ATOM', 'HETATM')):
                rn = ln[17:20]
                if rn in D_AMINO_ACID_PDB_TO_L:
                    l_name = D_AMINO_ACID_PDB_TO_L[rn]
                    gmx_d_name = D_AMINO_ACID_L_TO_GMX.get(l_name)
                    if gmx_d_name:
                        d_res_pos[(ln[21], ln[22:26].strip())] = (l_name, gmx_d_name)

        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
        fixer.findNonstandardResidues()
        # N-甲基氨基酸(SAR/MLE)在 PDBFixer 替换表中(SAR->GLY, MLE->LEU):
        # 若走 replaceNonstandardResidues 会被替换并删除 N-CH3(CN) 原子, 丢失甲基.
        # 因此把它们从替换列表中剔除, 保留原残基(含 CN), 由注册的带氢模板补氢.
        fixer.nonstandardResidues = [
            (r, w) for r, w in fixer.nonstandardResidues if r.name not in ('SAR', 'MLE')]  # type: ignore
        fixer.replaceNonstandardResidues()

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        # 非标准氨基酸(ncAA, 如 LE1 青霉胺): 将 storage/opmm 下的 ncAA 模板
        # 注册到 PDBFixer(registerTemplate), 使 addMissingHydrogens 能为其补氢.
        # 不注册的话 LE1 补不上氢, OpenMM createSystem 会因模板含氢而缺原子报错.
        n_registered = register_ncAA_templates(
            fixer, set(res.name for res in fixer.topology.residues()),
            custom=self.custom_residues)
        if n_registered:
            put_log(f'registered {n_registered} ncAA template(s) to PDBFixer for hydrogen addition')
        fixer.addMissingAtoms(seed=0)
        fixer.addMissingHydrogens()

        # 非标准残基的内部重原子键: 上游CIF文件未必带全 CONECT(如 PyMOL 转换时
        # 以数字开头的残基名可能整条无法推断), 缺键时 createSystem 会报
        # "bonds are different". 这里按 ncAA 模板统一补全, 已成键的不重复添加.
        topology, positions, n_bonds = sew_ncAA_internal_bonds(
            fixer.topology, fixer.positions, custom=self.custom_residues)
        if n_bonds:
            put_log(f'ncAA: added {n_bonds} internal bond(s) from templates')
        fixer.topology, fixer.positions = topology, positions

        # 将 D 残基名改回 GROMACS D 名(如 DALA), 供后续 OpenMM 模板匹配和 GMX 识别.
        for res in fixer.topology.residues():
            key = (res.chain.id, str(res.id))
            if key in d_res_pos:
                l_name, gmx_d_name = d_res_pos[key]
                if res.name == l_name:
                    put_log(f'restore D-amino acid {l_name} (res {key}) -> {gmx_d_name} '
                            f'after PDBFixer replacement')
                    res.name = gmx_d_name

        # 将单原子离子残基重命名为离子模板名(CAL/MG/...), 与 ions.xml 对齐.
        #   - 读取 PDB 时优先用元素列(77-78列); 若无元素列则按原子名猜测,
        #     其中 "单原子残基且原子名/残基名以 CA 开头" 会被猜为钙(elem.calcium).
        #     (用"单原子残基"限定的原因: 蛋白的 alpha 碳原子名也叫 CA, 但其残基含多个原子,
        #      不会触发这一分支, 从而避免把 alpha 碳误判为钙离子)
        #   - 因此这里 res 中唯一的原子, 其 .element 已经是解析好的 Element 对象(如 elem.calcium).
        # 本步(lazydock 侧): 以该 Element 对象为 key 查 ION_ELEMENT_TO_TEMPLATE, 命中则把
        # 残基名/原子名改写成 charmm36 兼容的离子模板名(如 CA -> CAL), 与 ions.xml 中
        # 的 Residue name 对齐; 否则 OpenMM 在 charmm36.xml 中找不到该单原子残基模板而报错.
        for res in fixer.topology.residues():
            atoms_ = list(res.atoms())
            # 判定依据为：单原子残基 & 元素已解析(现有模板: CAL/MG/ZN/SOD/POT/CLA)
            if len(atoms_) == 1 and atoms_[0].element is not None:
                template = ION_ELEMENT_TO_TEMPLATE.get(atoms_[0].element)
                if template and res.name != template:
                    put_log(f'rename ion residue {res.name}({atoms_[0].name}) -> {template} '
                            f'for OpenMM force field')
                    res.name = template
                    atoms_[0].name = template

        # 直接返回拓扑+坐标(而非字符串), 避免 writeFile 把 4 字符残基名(DALA)截断成 3 字符
        return self._dedup_bonds(fixer.topology, fixer.positions)

    @staticmethod
    def _dedup_bonds(topology, positions):
        """删除拓扑中的重复键(同一对原子被 addBond 多次).

        根因: PDBFixer._addAtomsToTopology 在重建拓扑时先 createStandardBonds()
        按 residues.xml 建立肽键(-C to N), 随后又遍历原拓扑把所有"至少一端为
        非标准残基"的键再 addBond 一次, 该循环没有已存在检查. 因此凡是非标准
        氨基酸与其相邻残基之间的肽键(以及 ncAA 端的环化键)都会变成 2 个 bond 对象,
        PDBFile.writeFile 每遇到一个 bond 对象就写一个邻居 serial, 输出即出现
        "CONECT 1 261 261", 按 PDB 约定被读成双键(N 表现为 4 连接或双键酰胺).

        力场建系本身不受影响(createSystem 按模板匹配, 实测无重复键), 但写出的
        PDB 拓扑语义错误, 下游按 CONECT 建系的工具会读错. 此处统一去重."""
        seen = set()
        kept = []
        n_dup = 0
        for bond in topology.bonds():
            a1, a2 = bond
            key = (a1.index, a2.index) if a1.index < a2.index else (a2.index, a1.index)
            if key in seen:
                n_dup += 1
                continue
            seen.add(key)
            kept.append(bond)
        if n_dup:
            topology._bonds = kept
            put_log(f'dedup: removed {n_dup} duplicate bond(s) from topology',
                    bg_color='blue')
        return topology, positions

    @staticmethod
    def _res_id(res):
        """残基 id 转 int, 失败返回 0(避免异常)."""
        try:
            return int(res.id)
        except (TypeError, ValueError):
            return 0

    def _sew_caps_bonds(self, topology, positions):
        """C 端酰胺封端(NH2)支持, 需在 _fix 之后、缝接二硫之前无条件调用.

        背景: PDBFixer 的 _addAtomsToTopology 重建拓扑时, 只恢复"至少一端非标准"
        的跨残基键(pdbfixer.py 键恢复逻辑). 当链尾为 C 端酰胺封端(NH2)时,
        CYS 与 NH2 都在 OpenMM Topology._standardBonds(residues.xml) 中,
        CYS14:C - NH215:N 这个跨残基肽键会被丢弃; 同时链尾残基被 PDBFixer
        当作蛋白 C 端, 用 missingTerminals 补了 OXT. 更麻烦的是: NH2 因无
        外部键会被 addMissingAtoms 拆成独立链(链 id 与主链相同, 单残基),
        导致:
          - createSystem(ignoreExternalBonds=False) 时 CYS/NH2 模板匹配失败
          - NH2 残基 N 无外部键, minimize 时被 LJ 推开
        本步: 对每个 NH2 残基, 按 chain.id 内找 id 比它小的最大残基作为前驱,
        删除前驱的 OXT, 并补回 C-N 肽键(NH2 残基本身原子/氢由 caps 模板补;
        二硫缝键在之后进行, 互不影响).
        返回 (topology, positions, n_oxt_removed, n_cn_added)."""
        modeller = openmm_app.Modeller(topology, positions)
        # NH2 可能被拆成独立链(chain.id 相同), 按 chain.id 聚合全部残基再按 id 匹配
        resis_by_chain = {}
        for chain in modeller.topology.chains():
            resis_by_chain.setdefault(chain.id, []).extend(chain.residues())
        caps = []  # (prev_res, nh2_res)
        for cid, resis in resis_by_chain.items():
            nh2_reses = [r for r in resis if r.name == 'NH2']
            for nh2_res in nh2_reses:
                prevs = [r for r in resis if r is not nh2_res
                         and self._res_id(r) < self._res_id(nh2_res)]
                if prevs:
                    caps.append((max(prevs, key=self._res_id), nh2_res))
        if not caps:
            return topology, positions, 0, 0
        to_delete = [a for prev, _ in caps for a in prev.atoms() if a.name == 'OXT']
        n_oxt = len(to_delete)
        if to_delete:
            modeller.delete(to_delete)
        n_cn = 0
        # delete 后拓扑重建, 按 (chain.id, res.id) 重新定位原子再补键
        for cid, pid, nid in sorted({(p.chain.id, p.id, n.id) for p, n in caps}):
            prev = next((r for c in modeller.topology.chains()
                         for r in c.residues() if c.id == cid and r.id == pid), None)
            nh2 = next((r for c in modeller.topology.chains()
                        for r in c.residues() if c.id == cid and r.id == nid), None)
            if prev is None or nh2 is None:
                continue
            atoms_prev = {a.name: a for a in prev.atoms()}
            atoms_cap = {a.name: a for a in nh2.atoms()}
            c_prev, n_cap = atoms_prev.get('C'), atoms_cap.get('N')
            if c_prev is None or n_cap is None:
                continue
            bonded = any({c_prev, n_cap} == set(b) for b in modeller.topology.bonds())
            if not bonded:
                modeller.topology.addBond(c_prev, n_cap)
                n_cn += 1
        return modeller.topology, modeller.positions, n_oxt, n_cn

    def _sew_backbone_chains(self, topology, positions, backbone_chains: List[str]):
        """在 PDBFixer 补全后的 PDB 上缝接 backbone 闭合链的首尾肽键:
        删除每条链尾残基的 OXT 与首残基 N 端多余 H(保留1个酰胺NH), 并 addBond(首N, 尾C).
        返回 (topology, positions, bonded_chain_ids)."""
        modeller = openmm_app.Modeller(topology, positions)
        bonded_chains = []
        for chain_id in backbone_chains:
            chain = next((c for c in modeller.topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            resis = list(chain.residues())
            if len(resis) < 2:
                continue
            names_f = {a.name: a for a in resis[0].atoms()}
            names_l = {a.name: a for a in resis[-1].atoms()}
            if 'N' not in names_f or 'C' not in names_l:
                continue
            # 首删尾残基 OXT(线性C端残留) 与 首残基N端多余H.
            # 头尾成环后 N 的成键状态取决于首残基类型:
            #   PRO: N 是亚胺(N连CD+CA+前残基C), 无H -> 删全部端H
            #   其他: N 是酰胺(NH), 保留 1 个端H -> 删 2 个多余的
            to_delete = []
            if 'OXT' in names_l:
                to_delete.append(names_l['OXT'])
            if resis[0].name == 'PRO':
                h_keys = ('H', 'H2', 'H1', 'H3')
            else:
                h_keys = ('H1', 'H2', 'H3')
            for hk in h_keys:
                if hk in names_f:
                    to_delete.append(names_f[hk])
            if to_delete:
                modeller.delete(to_delete)
            bonded_chains.append(chain_id)
        # delete 后 topology 重建, 需重新定位原子再加键
        for chain_id in bonded_chains:
            chain = next((c for c in modeller.topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            resis = list(chain.residues())
            names_f = {a.name: a for a in resis[0].atoms()}
            names_l = {a.name: a for a in resis[-1].atoms()}
            if 'N' in names_f and 'C' in names_l:
                n_at, c_at = names_f['N'], names_l['C']
                # 该键可能已由 PDBFixer 的异质原子键恢复循环带入(首/尾残基为 ncAA 时),
                # 此处不再重复添加, 否则 CONECT 出现重复邻居(被读作双键).
                if any({n_at, c_at} == set(b) for b in modeller.topology.bonds()):
                    continue
                modeller.topology.addBond(n_at, c_at)
        return modeller.topology, modeller.positions, bonded_chains

    def _remove_disulfide_hydrogens(self, modeller):
        """删除参与二硫键的巯基氢: 对 disulfide_chains 链内, 如果某硫原子(SG/SG3)
        已与另一硫原子成键(二硫键), 则删除其上的巯基氢.
        巯基氢命名因残基而异: 标准 CYS 为 HG, LE1 为 HG3, 其余含硫修饰残基可能为 HG1.
        必须在 _sew_disulfide_bonds 之后调用. 返回删除数."""
        n_removed = 0
        sulfur_hydrogens = {}
        for chain_id in self.disulfide_chains:
            chain = next((c for c in modeller.topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            s_atoms = [a for r in chain.residues() for a in r.atoms() if a.name in SULFUR_ATOM_NAMES]
            for a in s_atoms:
                sulfur_hydrogens[a] = [h for h in a.residue.atoms()
                                       if h.element == openmm_app.element.hydrogen
                                       and h.name in SULFUR_HYDROGEN_NAMES]
        # 找出已互相成键的硫原子对(二硫键)
        bonded_sg = set()
        for (b1, b2) in modeller.topology.bonds():
            if b1 in sulfur_hydrogens and b2 in sulfur_hydrogens:
                bonded_sg.add(b1)
                bonded_sg.add(b2)
        to_remove = [h for a in bonded_sg for h in sulfur_hydrogens[a]]
        if to_remove:
            modeller.delete(to_remove)
            n_removed = len(to_remove)
        return n_removed

    def _sew_disulfide_bonds(self, modeller):
        """在 modeller 拓扑上补全会 S-S 键: 对 disulfide_chains 链内,
        SG/SG3(标准CYS与含硫修饰氨基酸如LE1) 尚无键且距离 < ss_max_dist 的原子对 addBond.
        返回补键数."""
        n_added = 0
        for chain_id in self.disulfide_chains:
            chain = next((c for c in modeller.topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            sulfurs = [a for r in chain.residues() for a in r.atoms() if a.name in SULFUR_ATOM_NAMES]
            for i in range(len(sulfurs)):
                for j in range(i+1, len(sulfurs)):
                    bonded_already = any((a in (sulfurs[i], sulfurs[j]) and b in (sulfurs[i], sulfurs[j]))
                                         for a, b in modeller.topology.bonds())
                    if bonded_already:
                        continue
                    x1, y1, z1 = modeller.positions[sulfurs[i].index].value_in_unit(LENGTH)
                    x2, y2, z2 = modeller.positions[sulfurs[j].index].value_in_unit(LENGTH)
                    d = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ** 0.5
                    if d < self.ss_max_dist:
                        modeller.topology.addBond(sulfurs[i], sulfurs[j])
                        n_added += 1
        return n_added

    def _ensure_backbone_conect(self, pdb_str, topology):
        """补写头尾肽键(首残基 N - 尾残基 C)的 CONECT 记录.

        PDBFile.writeFooter 只写"至少一端是非标准残基"的键以及 CYS SG-SG,
        因此当环肽首尾都是标准残基时, 该键虽然存在于
        内存拓扑且 relax 全程受 CustomBondForce 保护(键长 1.34A 正确), 但写出的
        PDB 里没有任何 CONECT 记录, 下游按 CONECT 建系会读成线性肽.
        此处按写出的 serial 手工补一条 CONECT(已有的按要求跳过, 不产生重复邻居).
        返回 (pdb_str, 补写条数)."""
        if not self.backbone_chains:
            return pdb_str, 0
        pairs = []
        for chain_id in self.backbone_chains:
            chain = next((c for c in topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            resis = list(chain.residues())
            if len(resis) < 2:
                continue
            names_f = {a.name: a for a in resis[0].atoms()}
            names_l = {a.name: a for a in resis[-1].atoms()}
            if 'N' in names_f and 'C' in names_l:
                pairs.append((names_f['N'], names_l['C']))
        if not pairs:
            return pdb_str, 0
        # 从写出的文本反查 serial(keepIds=True 时残基 id 即拓扑 id)
        serial = {}
        existing = {}
        for ln in pdb_str.splitlines():
            if ln.startswith(('ATOM', 'HETATM')):
                try:
                    serial[(ln[21], ln[22:26].strip(), ln[12:16].strip())] = int(ln[6:11])
                except ValueError:
                    continue
            elif ln.startswith('CONECT'):
                nums = [int(ln[i:i+5]) for i in range(6, len(ln)-1, 5)
                        if ln[i:i+5].strip()]
                if nums:
                    existing.setdefault(nums[0], set()).update(nums[1:])
        new_lines = []
        for n_at, c_at in pairs:
            k1 = (n_at.residue.chain.id, str(n_at.residue.id).strip(), n_at.name)
            k2 = (c_at.residue.chain.id, str(c_at.residue.id).strip(), c_at.name)
            s1, s2 = serial.get(k1), serial.get(k2)
            if s1 is None or s2 is None:
                continue
            # 已写出(任一端 CONECT 里有对方)则跳过
            if s2 in existing.get(s1, set()) or s1 in existing.get(s2, set()):
                continue
            new_lines.append(f'CONECT{s1:5d}{s2:5d}')
        if not new_lines:
            return pdb_str, 0
        # 插到 END 之前(END 可能不存在则直接追加)
        lines = pdb_str.splitlines()
        end_idx = len(lines)
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].strip() == 'END':
                end_idx = i
                break
        lines[end_idx:end_idx] = new_lines
        return '\n'.join(lines) + '\n', len(new_lines)

    def _get_pdb_string(self, topology, positions):
        with io.StringIO() as f:
            openmm_app.PDBFile.writeFile(topology, positions, f, keepIds=True)
            pdb_str = f.getvalue()
        # PDBFile.writeFile 会把 4 字符残基名(如 DALA)截断成 3 字符(DAL).
        # 这里按 (chain, resSeq) 把 D 残基名恢复为完整 4 字符, 供下游 GROMACS 识别.
        d_res = {}  # (chain_id, resSeq) -> 完整 D 残基名
        for res in topology.residues():
            if res.name in D_AMINO_ACID_L_TO_GMX.values():
                d_res[(res.chain.id, res.id)] = res.name
        if d_res:
            lines = pdb_str.splitlines()
            for i, ln in enumerate(lines):
                if ln.startswith(('ATOM', 'HETATM')):
                    key = (ln[21], ln[22:26].strip())
                    if key in d_res:
                        lines[i] = ln[:17] + d_res[key] + ln[21:]
            pdb_str = '\n'.join(lines)
        # 头尾肽键(首N-尾C)可能因两端均为标准残基而未被 writeFooter 写入, 补写.
        pdb_str, n_bb = self._ensure_backbone_conect(pdb_str, topology)
        if n_bb:
            put_log(f'backbone: wrote {n_bb} missing CONECT record(s) for '
                    f'head-to-tail peptide bond(s)', bg_color='blue')
        return pdb_str

    def _add_backbone_bond_force(self, system, topology):
        """对每条 backbone 闭合链的首N-尾C 加 CustomBondForce 键势, 保护头尾肽键.
        注意: 全局参数名用唯一前缀 k_bb/r0_bb, 避免与 restrain 的 CustomExternalForce
        中同名参数 k 的 default value 冲突(openmm 要求 system 内同名参数 default 一致)."""
        if not self.backbone_chains:
            return
        cbf = openmm.CustomBondForce('0.5*k_bb*(r-r0_bb)^2')
        cbf.addGlobalParameter('k_bb', BACKBONE_BOND_STIFFNESS)
        cbf.addGlobalParameter('r0_bb', BACKBONE_BOND_LEN)
        n_added = 0
        for chain_id in self.backbone_chains:
            chain = next((c for c in topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            resis = list(chain.residues())
            if len(resis) < 2:
                continue
            names_f = {a.name: a for a in resis[0].atoms()}
            names_l = {a.name: a for a in resis[-1].atoms()}
            if 'N' in names_f and 'C' in names_l:
                cbf.addBond(names_f['N'].index, names_l['C'].index)
                n_added += 1
        if n_added:
            system.addForce(cbf)

    def _minimize_pdb(self, pdb, restrain_chain: Optional[Union[str, List[str]]] = None,
                      restrain_backbone: bool = False):
        # charmm36 系列自定义模板列表. 无条件全量加载, 由 OpenMM 自行按残基名匹配:
        #   - ions.xml: 离子模板(CAL/MG/ZN/SOD/POT/CLA), charmm36.xml 不含离子
        #   - charmm36_d_aa.xml: D-氨基酸模板(DALA/DPHE/...)
        #   - charmm36_nc_aa.xml: 非标准氨基酸(如 LE1 青霉胺/SAR 肌氨酸/MLE N-甲基亮氨酸)
        #   - charmm36_caps.xml: 封端模板(ACE N端乙酰 / NH2 C端酰胺)
        # 无需探测拓扑中是否含对应残基, 未命中的模板不影响其他残基建系.
        EXTRA_FF_FILES = ['opmm/ions.xml', 'opmm/charmm36_d_aa.xml',
                          'opmm/charmm36_nc_aa.xml', 'opmm/charmm36_caps.xml']
        # 运行时外部 ncAA 模板(任意路径), 追加在内建文件之后. 未传时列表不变.
        ff_files = (['charmm36.xml'] + [get_storage_path(f) for f in EXTRA_FF_FILES]
                    + self.custom_residues.xml_files)
        force_field = openmm_app.ForceField(*ff_files)  # referring to http://docs.openmm.org/latest/userguide/application/02_running_sims.html
        # ignoreExternalBonds 恒为 False. 该 flag 只影响"残基模板匹配"阶段(是否要求
        # 模板声明的外部键与拓扑完全一致), 不会删除任何键 -- data.bonds 直接来自
        # topology.bonds(), 普通肽键始终参与 HarmonicBond/Angle 建系.
        # 实测(环肽 + 另一条链含 5 对二硫键, 共 251 条肽键):
        #   True : 251/251 肽键有键势, 但 CYS 被静默匹配成 CYM 模板(SG 电荷 -0.80),
        #          S-S 键在拓扑中存在却拿不到键参数(CHARMM36 无 SS-SS, 只有 SM-SM),
        #          于是既无键势也无 1-2 排斥 -> 最小化时塌缩到 0.55A.
        #   False: DISU patch 正常应用(SG -> SM, -0.08), 5/5 对 S-S 收敛到 2.03A.
        # 三种 case 的结论:
        #   - 仅二硫键: 必须 False, 否则 DISU patch 无法应用;
        #   - 仅 backbone 闭合: False 与 True 等价(实测头尾键 1.336 vs 1.340A,
        #     肽键 0/247 异常). 首/尾残基只有 1 个跨残基键, 签名与普通残基一致,
        #     无需放宽匹配; 头尾键由下方 CustomBondForce 单独保护, 不依赖该 flag;
        #   - 两者并存: 必须 False(二硫键需求压过 backbone).
        # 因此硬编码 False, 不再按 has_bb/has_disu 推导.
        system = force_field.createSystem(pdb.topology, constraints=self.constraints,
                                          ignoreExternalBonds=False)
        # 头尾肽键保护(首 N - 尾 C), 该键无标准 patch 参数, 用 CustomBondForce 兜底
        self._add_backbone_bond_force(system, pdb.topology)
        # 二硫键无需手动加键势: ignoreExternalBonds=False 时 DISU patch 会自动
        # 应用(SM-SM 有真实键/角/二面参数).

        # Add constraints to restrain_chain
        restrain_chain = restrain_chain or []
        restrain_name = ['N', 'CA', 'C', 'O'] if restrain_backbone else []
        if restrain_chain:
            force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
            force.addGlobalParameter("k", self.stiffness)
            for p in ["x0", "y0", "z0"]:
                force.addPerParticleParameter(p)

            for i, a in enumerate(pdb.topology.atoms()):
                if a.residue.chain.id in restrain_chain and (not restrain_backbone or a.name in restrain_name):
                    force.addParticle(i, pdb.positions[i])

            system.addForce(force)

        # Set up the integrator and simulation
        integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
        platform = openmm.Platform.getPlatformByName(self.platform)
        if self.platform == 'CUDA' and self.device_index is not None:
            # DeviceIndex 为进程级全局设置; 多卡并行时每个 worker 子进程只设一块卡
            platform.setPropertyDefaultValue('DeviceIndex', str(self.device_index))
        simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)

        # Perform minimization
        ret = {}
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        ret["einit"] = state.getPotentialEnergy().value_in_unit(ENERGY)
        ret["posinit"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)

        simulation.minimizeEnergy(maxIterations=self.max_iterations, tolerance=self.tolerance)

        state = simulation.context.getState(getEnergy=True, getPositions=True)
        ret["efinal"] = state.getPotentialEnergy().value_in_unit(ENERGY)
        ret["pos"] = state.getPositions(asNumpy=True).value_in_unit(LENGTH)
        ret["min_pdb"] = self._get_pdb_string(simulation.topology, state.getPositions())

        return ret['min_pdb'], ret

    def _add_energy_remarks(self, pdb_str, ret):
        pdb_lines = pdb_str.splitlines()
        pdb_lines.insert(1, "REMARK   1  FINAL ENERGY:   {:.3f} KCAL/MOL".format(ret['efinal']))
        pdb_lines.insert(1, "REMARK   1  INITIAL ENERGY: {:.3f} KCAL/MOL".format(ret['einit']))
        return "\n".join(pdb_lines)

    def __call__(self, pdb_str, out_path, restrain_chain: Optional[Union[str, List[str]]] = None, restrain_backbone: bool = False, return_info=True):
        if '\n' not in pdb_str and pdb_str.lower().endswith(".pdb"):
            with open(pdb_str) as f:
                pdb_str = f.read()

        pdb_fixed_topology, pdb_fixed_positions = self._fix(pdb_str)
        # C 端酰胺封端(NH2): PDBFixer 重建拓扑会丢 NH2.N-前残基.C 肽键并给前残基补 OXT,
        # 无条件恢复(仅对链尾确为 NH2 的链生效, 其余链无副作用).
        pdb_fixed_topology, pdb_fixed_positions, n_oxt, n_cn = self._sew_caps_bonds(
            pdb_fixed_topology, pdb_fixed_positions)
        if n_oxt or n_cn:
            put_log(f'caps: removed {n_oxt} OXT, added {n_cn} C-N bond(s) for C-terminal amide cap(NH2)')
        # 统一走内存拓扑传递(backbone 闭合/二硫键/D残基名均需完整保留, 避免字符串往返丢失)
        if self.backbone_chains or self.disulfide_chains:
            # 缝接 backbone 闭合链/二硫键
            if self.backbone_chains:
                topology, positions, _bonded = self._sew_backbone_chains(
                    pdb_fixed_topology, pdb_fixed_positions, self.backbone_chains)
                modeller = openmm_app.Modeller(topology, positions)
            else:
                modeller = openmm_app.Modeller(pdb_fixed_topology, pdb_fixed_positions)
            if self.disulfide_chains:
                self._sew_disulfide_bonds(modeller)
                # 删除参与二硫键硫上的巯基氢(如 LE1 的 HG3 / CYS 的 HG1), 拓扑上已无氢,
                # 使 OpenMM 的 LE1SS/DISU patch 能正确匹配(CYX 式无氢二硫)
                self._remove_disulfide_hydrogens(modeller)
            pdb_like = SimpleNamespace(topology=modeller.topology, positions=modeller.positions)
            pdb_min, ret = self._minimize_pdb(pdb_like, restrain_chain, restrain_backbone)
        else:
            pdb_like = SimpleNamespace(topology=pdb_fixed_topology, positions=pdb_fixed_positions)
            pdb_min, ret = self._minimize_pdb(pdb_like, restrain_chain, restrain_backbone)
        pdb_min = self._add_energy_remarks(pdb_min, ret)
        if out_path and os.path.exists(out_path):
            with open(out_path, 'w') as f:
                f.write(pdb_min)
        if return_info:
            return pdb_min, ret
        else:
            return pdb_min


if __name__ == '__main__':
    import sys
    force_field = ForceFieldMinimizer(stiffness=10**8)
    # force_field(sys.argv[1], sys.argv[2], restrain_chain=[sys.argv[3]])
    force_field('data_tmp/docking/CB1R_0.pdb', 'data_tmp/docking/CB1R_0_openmm_relax_test.pdb',
                restrain_chain=['A'], restrain_backbone=True)
