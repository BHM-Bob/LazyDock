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

from lazydock.utils import get_storage_path

ENERGY = unit.kilojoules_per_mole
LENGTH = unit.nanometer  # pyright: ignore[reportAttributeAccessIssue]

# 常见离子元素 -> 离子力场模板名(charmm36.xml 不含离子, 需额外加载 ions.xml)
from openmm.app.element import (calcium, chlorine, magnesium, potassium,
                                sodium, zinc)

ION_ELEMENT_TO_TEMPLATE = {
    calcium: 'CAL', magnesium: 'MG', zinc: 'ZN',
    sodium: 'SOD', potassium: 'POT', chlorine: 'CLA',
}
ION_TEMPLATE_NAMES = set(ION_ELEMENT_TO_TEMPLATE.values())


class ForceFieldMinimizer(object):

    def __init__(self, stiffness=10.0, max_iterations=0, tolerance=10, platform='CUDA', constraints=openmm_app.HBonds,
                 cyclic_chains: Optional[List[str]] = None,
                 cyclic_bond_len: float = 0.134, cyclic_bond_stiffness: float = 3e5,
                 disulfide_chains: Optional[List[str]] = None,
                 ss_bond_len: float = 0.2029, ss_bond_stiffness: float = 144766,
                 ss_max_dist: float = 0.3,
                 device_index: Optional[int] = None):
        super().__init__()
        self.stiffness = stiffness * ENERGY/(LENGTH ** 2) # coefficient for restraints force, unit: kJ*mol^-1*nm^-2  # pyright: ignore[reportOperatorIssue]
        self.max_iterations = max_iterations
        self.tolerance = tolerance * ENERGY/LENGTH  # pyright: ignore[reportOperatorIssue]
        assert platform in ('CUDA', 'CPU')
        self.platform = platform
        self.constraints = constraints
        # 环肽(碳骨架)支持: 对 cyclic_chains 中的链, 在 PDBFixer 补全后缝接首尾肽键,
        # 并加 CustomBondForce 保护该环化键, 避免最小化时断链
        self.cyclic_chains = list(cyclic_chains or [])
        self.cyclic_bond_len = cyclic_bond_len  # nm, 环化键平衡长(酰胺 C-N ~1.34A)
        self.cyclic_bond_stiffness = cyclic_bond_stiffness  # kJ/mol/nm^2
        # 二硫键支持: 对 disulfide_chains 中的链, 补全 SG-SG 键.
        # 当 ignoreExternalBonds=False 时 openmm 自动应用 CHARMM36 DISU patch
        # (S-S 平衡键长 0.2029nm, k=144766 kJ/mol/nm^2), 无需手动键势;
        # 当与 cyclic_chains 并存(forced ignoreExternalBonds=True)时, 用 CustomBondForce 手动保护.
        self.disulfide_chains = list(disulfide_chains or [])
        self.ss_bond_len = ss_bond_len  # nm, 二硫键 SG-SG 平衡长(CHARMM36: 2.029A)
        self.ss_bond_stiffness = ss_bond_stiffness  # kJ/mol/nm^2
        self.ss_max_dist = ss_max_dist  # nm, SG-SG 距离阈值内判定为二硫键并补键(3.0A)
        # GPU 设备选择: 仅 CUDA 时生效. 向前兼容: None 表示用默认卡(不设置 DeviceIndex)
        self.device_index = device_index

    def _fix(self, pdb_str):
        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms(seed=0)
        fixer.addMissingHydrogens()

        # 上游(OpenMM PDBFile 解析)已经把离子残基的 element 识别好了:
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

        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        return out_handle.getvalue()

    def _sew_cyclic_chains(self, pdb_str: str, cyclic_chains: List[str]):
        """在 PDBFixer 补全后的 PDB 上缝接环化链首尾肽键:
        删除每条链尾残基的 OXT 与首残基 N 端多余 H(保留1个酰胺NH), 并 addBond(首N, 尾C).
        返回 (topology, positions, bonded_chain_ids)."""
        pdb = openmm_app.PDBFile(io.StringIO(pdb_str))
        modeller = openmm_app.Modeller(pdb.topology, pdb.positions)
        bonded_chains = []
        for chain_id in cyclic_chains:
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
            # 环化后 N 的成键状态取决于首残基类型:
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
                modeller.topology.addBond(names_f['N'], names_l['C'])
        return modeller.topology, modeller.positions, bonded_chains

    def _sew_disulfide_bonds(self, modeller):
        """在 modeller 拓扑上补全会 S-S 键: 对 disulfide_chains 链内,
        SG-SG 尚无键且距离 < ss_max_dist 的原子对 addBond.
        返回补键数."""
        n_added = 0
        for chain_id in self.disulfide_chains:
            chain = next((c for c in modeller.topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            sgs = [a for r in chain.residues() for a in r.atoms() if a.name == 'SG']
            for i in range(len(sgs)):
                for j in range(i+1, len(sgs)):
                    bonded_already = any((a in (sgs[i], sgs[j]) and b in (sgs[i], sgs[j]))
                                         for a, b in modeller.topology.bonds())
                    if bonded_already:
                        continue
                    x1, y1, z1 = modeller.positions[sgs[i].index].value_in_unit(LENGTH)
                    x2, y2, z2 = modeller.positions[sgs[j].index].value_in_unit(LENGTH)
                    d = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ** 0.5
                    if d < self.ss_max_dist:
                        modeller.topology.addBond(sgs[i], sgs[j])
                        n_added += 1
        return n_added

    def _get_pdb_string(self, topology, positions):
        with io.StringIO() as f:
            openmm_app.PDBFile.writeFile(topology, positions, f, keepIds=True)
            return f.getvalue()

    def _add_cyclic_bond_force(self, system, topology):
        """对每条环化链的首N-尾C 加 CustomBondForce 键势, 保护环化键.
        注意: 全局参数名用唯一前缀 k_cyc/r0_cyc, 避免与 restrain 的 CustomExternalForce
        中同名参数 k 的 default value 冲突(openmm 要求 system 内同名参数 default 一致)."""
        if not (self.cyclic_chains and self.cyclic_bond_stiffness > 0):
            return
        cbf = openmm.CustomBondForce('0.5*k_cyc*(r-r0_cyc)^2')
        cbf.addGlobalParameter('k_cyc', self.cyclic_bond_stiffness)
        cbf.addGlobalParameter('r0_cyc', self.cyclic_bond_len)
        n_added = 0
        for chain_id in self.cyclic_chains:
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

    def _add_disulfide_bond_force(self, system, topology):
        """对 disulfide_chains 链内所有 SG-SG 加 CustomBondForce 保护.
        仅在 ignoreExternalBonds=True(与 cyclic_chains 并存, DISU patch 无法自动应用)时调用."""
        if not (self.disulfide_chains and self.ss_bond_stiffness > 0):
            return
        cbf = openmm.CustomBondForce('0.5*k_ss*(r-r0_ss)^2')
        cbf.addGlobalParameter('k_ss', self.ss_bond_stiffness)
        cbf.addGlobalParameter('r0_ss', self.ss_bond_len)
        n_added = 0
        for chain_id in self.disulfide_chains:
            chain = next((c for c in topology.chains() if c.id == chain_id), None)
            if chain is None:
                continue
            sgs = [a for r in chain.residues() for a in r.atoms() if a.name == 'SG']
            for i in range(len(sgs)):
                for j in range(i+1, len(sgs)):
                    cbf.addBond(sgs[i].index, sgs[j].index)
                    n_added += 1
        if n_added:
            system.addForce(cbf)

    def _minimize_pdb(self, pdb, restrain_chain: Optional[Union[str, List[str]]] = None,
                      restrain_backbone: bool = False):
        # charmm36.xml 不含离子模板; 若拓扑中有单原子离子残基, 追加加载 ions.xml
        ion_resids = ION_TEMPLATE_NAMES.intersection(
            res.name for res in pdb.topology.residues())
        if ion_resids:
            ions_xml = get_storage_path('opmm/ions.xml')
            force_field = openmm_app.ForceField("charmm36.xml", ions_xml)
            put_log(f'loaded ions.xml for residues: {sorted(ion_resids)}', bg_col='blue')
        else:
            force_field = openmm_app.ForceField("charmm36.xml") # referring to http://docs.openmm.org/latest/userguide/application/02_running_sims.html
        # ignoreExternalBonds 判定:
        #   - 仅二硫键: False -> openmm 自动应用 CHARMM36 DISU patch (S-S 有真实键/角/二面参数, 收敛2.03A)
        #   - 仅碳骨架环: True -> 头尾跨残基键无标准 patch, 忽略后可手动加键势保护
        #   - 两者并存: True -> DISU patch 不可用, 需手动 S-S CustomBondForce 保护
        has_cyc, has_disu = bool(self.cyclic_chains), bool(self.disulfide_chains)
        ignore_external = has_cyc or not has_disu
        system = force_field.createSystem(pdb.topology, constraints=self.constraints,
                                          ignoreExternalBonds=ignore_external)
        # 环化键保护
        self._add_cyclic_bond_force(system, pdb.topology)
        # 二硫键保护(仅 ignoreExternalBonds=True 时需手动)
        if has_disu and ignore_external:
            self._add_disulfide_bond_force(system, pdb.topology)

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

    def _minimize(self, pdb_str: str, restrain_chain: Optional[Union[str, List[str]]] = None,
                  restrain_backbone: bool = False):
        pdb = openmm_app.PDBFile(io.StringIO(pdb_str))
        return self._minimize_pdb(pdb, restrain_chain, restrain_backbone)

    def _add_energy_remarks(self, pdb_str, ret):
        pdb_lines = pdb_str.splitlines()
        pdb_lines.insert(1, "REMARK   1  FINAL ENERGY:   {:.3f} KCAL/MOL".format(ret['efinal']))
        pdb_lines.insert(1, "REMARK   1  INITIAL ENERGY: {:.3f} KCAL/MOL".format(ret['einit']))
        return "\n".join(pdb_lines)

    def __call__(self, pdb_str, out_path, restrain_chain: Optional[Union[str, List[str]]] = None, restrain_backbone: bool = False, return_info=True):
        if '\n' not in pdb_str and pdb_str.lower().endswith(".pdb"):
            with open(pdb_str) as f:
                pdb_str = f.read()

        pdb_fixed = self._fix(pdb_str)
        if self.cyclic_chains or self.disulfide_chains:
            # 缝接环化链/二硫键, 用内存拓扑传递, 避免字符串往返丢失 addBond 键
            if self.cyclic_chains:
                topology, positions, _bonded = self._sew_cyclic_chains(pdb_fixed, self.cyclic_chains)
                modeller = openmm_app.Modeller(topology, positions)
            else:
                pdb = openmm_app.PDBFile(io.StringIO(pdb_fixed))
                modeller = openmm_app.Modeller(pdb.topology, pdb.positions)
            if self.disulfide_chains:
                self._sew_disulfide_bonds(modeller)
            pdb_like = SimpleNamespace(topology=modeller.topology, positions=modeller.positions)
            pdb_min, ret = self._minimize_pdb(pdb_like, restrain_chain, restrain_backbone)
        else:
            pdb_min, ret = self._minimize(pdb_fixed, restrain_chain, restrain_backbone)
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
