#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OpenMM约束系统详解：createSystem constraints vs CustomExternalForce

这个测试对比两种不同的约束实现方式，帮助理解它们的工作原理和区别。
"""

import io
import os
import numpy as np
import pdbfixer
from openmm import app as openmm_app
from openmm import openmm, unit
import unittest


def fix_pdb(pdb_str):
    """Fix PDB structure using pdbfixer"""
    fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    
    out_handle = io.StringIO()
    openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
    return out_handle.getvalue()


class TestConstraintsVsExternalForce(unittest.TestCase):
    """对比测试：createSystem constraints vs CustomExternalForce"""
    
    @classmethod
    def setUpClass(cls):
        """加载测试数据"""
        pdb_path = 'data_tmp/pdb/pdbstr.pdb'
        with open(pdb_path, 'r') as f:
            cls.pdb_str = f.read()
        cls.fixed_pdb_str = fix_pdb(cls.pdb_str)
        cls.pdb = openmm_app.PDBFile(io.StringIO(cls.fixed_pdb_str))
    
    def test_1_constraints_none(self):
        """Test 1: createSystem constraints=None - 无约束"""
        print("\n" + "="*70)
        print("Test 1: createSystem(constraints=None) - 无任何键约束")
        print("="*70)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        # constraints=None: 不添加任何键约束
        system = force_field.createSystem(self.pdb.topology, constraints=None)
        
        # 检查系统中的力
        print("\n系统中的力:")
        for i, force in enumerate(system.getForces()):
            print(f"  {i}: {force.__class__.__name__}")
        
        # 运行最小化
        self._run_minimization(system, "无约束")
    
    def test_2_constraints_hbonds(self):
        """Test 2: createSystem constraints=HBonds - 只约束H-X键"""
        print("\n" + "="*70)
        print("Test 2: createSystem(constraints=HBonds) - 只约束涉及H的键")
        print("="*70)
        print("""
约束说明：
- HBonds: 只约束涉及氢原子的键（N-H, O-H, C-H等）
- 这在标准MD中很常见，因为H质量小，需要更小的积分步长
- 约束后，H原子与重原子之间的距离被固定
        """)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        # constraints=HBonds: 只约束涉及H的键
        system = force_field.createSystem(self.pdb.topology, constraints=openmm_app.HBonds)
        
        print("\n系统中的力:")
        for i, force in enumerate(system.getForces()):
            print(f"  {i}: {force.__class__.__name__}")
        
        result = self._run_minimization(system, "HBonds约束")
        
        # 分析涉及H的键
        h_bonds = []
        for bond in self.pdb.topology.bonds():
            if 'H' in bond[0].name or 'H' in bond[1].name:
                h_bonds.append((bond[0].name, bond[1].name))
        
        print(f"\n涉及H的键数量: {len(h_bonds)}")
        print(f"示例: {h_bonds[:5]}")
        
        return result
    
    def test_3_constraints_allbonds(self):
        """Test 3: createSystem constraints=AllBonds - 约束所有键"""
        print("\n" + "="*70)
        print("Test 3: createSystem(constraints=AllBonds) - 约束所有键")
        print("="*70)
        print("""
约束说明：
- AllBonds: 约束所有键，不仅是H参与的键
- 进一步减少系统自由度
- 用于更快的模拟，但可能影响精度
        """)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        system = force_field.createSystem(self.pdb.topology, constraints=openmm_app.AllBonds)
        
        print("\n系统中的力:")
        for i, force in enumerate(system.getForces()):
            print(f"  {i}: {force.__class__.__name__}")
        
        return self._run_minimization(system, "AllBonds约束")
    
    def test_4_constraints_hangles(self):
        """Test 4: createSystem constraints=HAngles - 约束键和H参与的键角"""
        print("\n" + "="*70)
        print("Test 4: createSystem(constraints=HAngles) - 约束键和H参与的键角")
        print("="*70)
        print("""
约束说明：
- HAngles: 约束所有键 + 涉及H的键角
- 最严格的约束选项
- 进一步减少自由度，允许更大的积分步长
        """)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        system = force_field.createSystem(self.pdb.topology, constraints=openmm_app.HAngles)
        
        print("\n系统中的力:")
        for i, force in enumerate(system.getForces()):
            print(f"  {i}: {force.__class__.__name__}")
        
        return self._run_minimization(system, "HAngles约束")
    
    def test_5_custom_external_force_high_stiffness(self):
        """Test 5: CustomExternalForce (高刚度) - 模拟刚性约束"""
        print("\n" + "="*70)
        print("Test 5: CustomExternalForce (stiffness=10^10) - 模拟刚性位置约束")
        print("="*70)
        print("""
约束说明（你的API实现）：
- CustomExternalForce: 自定义外部力
- 公式: E = 0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
- 这是一个软约束，允许原子在平衡位置附近振动
- 即使k很大，也是软约束，不是刚性约束
        """)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        system = force_field.createSystem(self.pdb.topology, constraints=None)
        
        # 添加CustomExternalForce约束（高刚度）
        force = openmm.CustomExternalForce(
            "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
        )
        force.addGlobalParameter("k", 10**10 * unit.kilojoules_per_mole / (unit.nanometer**2))
        for p in ["x0", "y0", "z0"]:
            force.addPerParticleParameter(p)
        
        # 约束配体的所有原子
        for i, atom in enumerate(self.pdb.topology.atoms()):
            if atom.residue.chain.id == 'B':  # 配体链
                force.addParticle(i, self.pdb.positions[i])
        
        system.addForce(force)
        
        print("\n系统中的力:")
        for i, force in enumerate(system.getForces()):
            print(f"  {i}: {force.__class__.__name__}")
        
        return self._run_minimization(system, "CustomExternalForce(高刚度)")
    
    def test_6_custom_external_force_low_stiffness(self):
        """Test 6: CustomExternalForce (低刚度) - 软约束"""
        print("\n" + "="*70)
        print("Test 6: CustomExternalForce (stiffness=10^4) - 软位置约束")
        print("="*70)
        print("""
约束说明：
- 较低刚度的软约束
- 允许原子有较大位移
- 常用于引导结构优化
        """)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        system = force_field.createSystem(self.pdb.topology, constraints=None)
        
        # 添加CustomExternalForce约束（低刚度）
        force = openmm.CustomExternalForce(
            "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
        )
        force.addGlobalParameter("k", 10**4 * unit.kilojoules_per_mole / (unit.nanometer**2))
        for p in ["x0", "y0", "z0"]:
            force.addPerParticleParameter(p)
        
        # 约束配体的所有原子
        for i, atom in enumerate(self.pdb.topology.atoms()):
            if atom.residue.chain.id == 'B':
                force.addParticle(i, self.pdb.positions[i])
        
        system.addForce(force)
        
        return self._run_minimization(system, "CustomExternalForce(低刚度)")
    
    def test_7_backbone_vs_sidechain(self):
        """Test 7: 你的API实现 - 约束主链原子"""
        print("\n" + "="*70)
        print("Test 7: 约束主链原子 (N, CA, C, O)")
        print("="*70)
        print("""
你的API实现：
- restrain_chain=['A'], restrain_backbone=True
- 只约束选定链的主链原子
- 使用中等刚度(10^8)
        """)
        
        force_field = openmm_app.ForceField("charmm36.xml")
        system = force_field.createSystem(self.pdb.topology, constraints=openmm_app.HBonds)
        
        # 添加主链约束
        force = openmm.CustomExternalForce(
            "0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)"
        )
        force.addGlobalParameter("k", 10**8 * unit.kilojoules_per_mole / (unit.nanometer**2))
        for p in ["x0", "y0", "z0"]:
            force.addPerParticleParameter(p)
        
        backbone_atoms = ['N', 'CA', 'C', 'O']
        constrained_count = 0
        
        for i, atom in enumerate(self.pdb.topology.atoms()):
            if atom.residue.chain.id == 'A':  # 受体链
                if atom.name in backbone_atoms and atom.element.name != 'hydrogen':
                    force.addParticle(i, self.pdb.positions[i])
                    constrained_count += 1
        
        print(f"\n约束的主链原子数: {constrained_count}")
        system.addForce(force)
        
        return self._run_minimization(system, "约束主链")
    
    def _run_minimization(self, system, name):
        """运行能量最小化并分析结果"""
        integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
        platform = openmm.Platform.getPlatformByName("CUDA")
        simulation = openmm_app.Simulation(self.pdb.topology, system, integrator, platform)
        simulation.context.setPositions(self.pdb.positions)
        
        # 获取初始状态
        state_initial = simulation.context.getState(getEnergy=True, getPositions=True)
        pos_initial = state_initial.getPositions(asNumpy=True)
        energy_initial = state_initial.getPotentialEnergy()
        
        print(f"\n初始能量: {energy_initial.value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
        
        # 执行最小化
        simulation.minimizeEnergy(maxIterations=1000, 
                                 tolerance=10 * unit.kilojoules_per_mole / unit.nanometer)
        
        # 获取最终状态
        state_final = simulation.context.getState(getEnergy=True, getPositions=True)
        pos_final = state_final.getPositions(asNumpy=True)
        energy_final = state_final.getPotentialEnergy()
        
        energy_change = (energy_final - energy_initial).value_in_unit(unit.kilocalories_per_mole)
        print(f"最终能量: {energy_final.value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
        print(f"能量变化: {energy_change:+.2f} kcal/mol")
        
        # 分析配体移动
        max_move = 0.0
        avg_move = 0.0
        moved_count = 0
        ligand_count = 0
        
        for i, atom in enumerate(self.pdb.topology.atoms()):
            if atom.residue.chain.id == 'B' and atom.element.name != 'hydrogen':
                diff = pos_final[i] - pos_initial[i]
                distance = np.linalg.norm(diff) * 10  # 转换为Å
                max_move = max(max_move, distance)
                avg_move += distance
                ligand_count += 1
                if distance > 0.01:
                    moved_count += 1
        
        avg_move /= ligand_count if ligand_count > 0 else 1
        
        print(f"\n配体B链移动分析:")
        print(f"  - 移动原子数: {moved_count}/{ligand_count}")
        print(f"  - 最大移动: {max_move:.4f} Å")
        print(f"  - 平均移动: {avg_move:.4f} Å")
        
        return {
            'name': name,
            'initial_energy': energy_initial,
            'final_energy': energy_final,
            'energy_change': energy_change,
            'max_ligand_move': max_move,
            'avg_ligand_move': avg_move,
            'pos_initial': pos_initial,
            'pos_final': pos_final
        }


class TestRootCauseAnalysis(unittest.TestCase):
    """根本原因分析测试"""
    
    def test_hbond_constraints_issue(self):
        """分析HBond约束对ML生成配体的影响"""
        print("\n" + "="*70)
        print("根本原因分析：HBond约束对ML生成配体的影响")
        print("="*70)
        print("""
问题分析：
1. 你的API使用了: createSystem(constraints=openmm_app.HBonds)
2. ML生成的配体氢键几何结构可能不符合标准力场参数
3. HBond约束会强制保持特定的键长，导致系统能量无法正常优化

对比：
- createSystem的constraints参数：约束的是化学键（N-H, O-H等）
- CustomExternalForce：约束的是原子位置（空间位置）

结论：
ML生成的配体问题不在于"位置约束"，而在于"化学键约束"！
        """)
        
        pdb_path = 'data_tmp/pdb/pdbstr.pdb'
        with open(pdb_path, 'r') as f:
            pdb_str = f.read()
        
        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms(seed=0)
        fixer.addMissingHydrogens()
        
        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        fixed_pdb = openmm_app.PDBFile(io.StringIO(out_handle.getvalue()))
        
        # 统计氢键数量
        h_bonds_in_ligand = 0
        h_bonds_in_receptor = 0
        
        for bond in fixed_pdb.topology.bonds():
            atom1, atom2 = bond
            has_h = 'H' in atom1.name or 'H' in atom2.name
            
            if atom1.residue.chain.id == 'B' or atom2.residue.chain.id == 'B':
                if has_h:
                    h_bonds_in_ligand += 1
            else:
                if has_h:
                    h_bonds_in_receptor += 1
        
        print(f"\n配体B链中涉及H的键数量: {h_bonds_in_ligand}")
        print(f"受体A链中涉及H的键数量: {h_bonds_in_receptor}")
        
        print("""
建议解决方案：
1. 使用constraints=None，让力场自动优化键长
2. 或使用更宽松的约束系统
3. 确保配体残基类型被正确定义
        """)


if __name__ == "__main__":
    print("="*70)
    print("OpenMM约束系统详解")
    print("createSystem constraints vs CustomExternalForce 对比测试")
    print("="*70)
    print("""
关键区别：
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
| 特性              | createSystem constraints | CustomExternalForce    |
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
| 约束对象          | 化学键/键角              | 原子位置               |
| 约束类型          | 刚性约束                 | 软约束                 |
| 数学实现          | 直接消除自由度           | 添加势能项             |
| 原子移动          | 完全固定                 | 允许小幅度振动         |
| 计算效率          | 更快                     | 较慢                   |
| 适用场景          | MD模拟加速               | 结构优化/引导         |
| 刚度控制          | 固定(完全约束)           | 可调(0.5*k*(x-x0)^2)   |
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

你的问题根源：
createSystem的constraints=HBonds会约束所有涉及H的化学键，
对于ML生成的配体，这可能导致不合理的键长被强制固定，
使能量最小化无法进行。
    """)
    
    unittest.main(verbosity=2)
