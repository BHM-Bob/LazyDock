#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
测试新的ForceFieldMinimizer API
验证：constraints=None + 约束受体backbone
期望结果：配体动、侧链动、主链不动
"""

import io
import os
import tempfile
import numpy as np
from openmm import app as openmm_app
from openmm import openmm, unit

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from lazydock.opmm.relax import ForceFieldMinimizer


def analyze_movements(pos_initial, pos_final, pdb, chain_id, backbone_atoms):
    """分析原子的移动情况"""
    results = {
        'backbone': {'max': 0, 'avg': 0, 'moved': 0, 'total': 0},
        'sidechain': {'max': 0, 'avg': 0, 'moved': 0, 'total': 0},
        'all': {'max': 0, 'avg': 0, 'moved': 0, 'total': 0}
    }
    
    backbone_moves = []
    sidechain_moves = []
    all_moves = []
    
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == chain_id and atom.element.name != 'hydrogen':
            diff = pos_final[i] - pos_initial[i]
            distance = np.linalg.norm(diff) * 10  # 转换为Å
            all_moves.append(distance)
            
            if atom.name in backbone_atoms:
                backbone_moves.append(distance)
            else:
                sidechain_moves.append(distance)
    
    # 汇总结果
    if backbone_moves:
        results['backbone']['total'] = len(backbone_moves)
        results['backbone']['max'] = max(backbone_moves)
        results['backbone']['avg'] = np.mean(backbone_moves)
        results['backbone']['moved'] = sum(1 for m in backbone_moves if m > 0.01)
    
    if sidechain_moves:
        results['sidechain']['total'] = len(sidechain_moves)
        results['sidechain']['max'] = max(sidechain_moves)
        results['sidechain']['avg'] = np.mean(sidechain_moves)
        results['sidechain']['moved'] = sum(1 for m in sidechain_moves if m > 0.01)
    
    if all_moves:
        results['all']['total'] = len(all_moves)
        results['all']['max'] = max(all_moves)
        results['all']['avg'] = np.mean(all_moves)
        results['all']['moved'] = sum(1 for m in all_moves if m > 0.01)
    
    return results


def test_new_api():
    """测试新的API：constraints=None + 约束受体backbone"""
    print("\n" + "="*70)
    print("测试新API: constraints=None + 约束受体backbone")
    print("="*70)
    
    # 读取PDB
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    # 创建两个minimizer进行对比
    print("\n--- 对比测试 ---")
    
    # 测试1: 旧行为 (constraints=HBonds)
    print("\n[测试1] 旧行为: constraints=HBonds (默认)")
    minimizer_old = ForceFieldMinimizer(
        stiffness=10**8, 
        max_iterations=1000, 
        tolerance=10,
        platform='CUDA',
        constraints=openmm_app.HBonds  # 默认值
    )
    
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        output_path_old = tmp.name
    
    try:
        pdb_min_old, info_old = minimizer_old(
            pdb_str, output_path_old, 
            restrain_chain=['A'], 
            restrain_backbone=True
        )
        
        print(f"  初始能量: {info_old['einit']:.2f} kcal/mol")
        print(f"  最终能量: {info_old['efinal']:.2f} kcal/mol")
        print(f"  能量变化: {info_old['efinal'] - info_old['einit']:+.2f} kcal/mol")
        
    finally:
        if os.path.exists(output_path_old):
            os.unlink(output_path_old)
    
    # 测试2: 新行为 (constraints=None)
    print("\n[测试2] 新行为: constraints=None")
    minimizer_new = ForceFieldMinimizer(
        stiffness=10**8, 
        max_iterations=1000, 
        tolerance=10,
        platform='CUDA',
        constraints=None  # 新参数
    )
    
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        output_path_new = tmp.name
    
    try:
        pdb_min_new, info_new = minimizer_new(
            pdb_str, output_path_new, 
            restrain_chain=['A'], 
            restrain_backbone=True
        )
        
        print(f"  初始能量: {info_new['einit']:.2f} kcal/mol")
        print(f"  最终能量: {info_new['efinal']:.2f} kcal/mol")
        print(f"  能量变化: {info_new['efinal'] - info_new['einit']:+.2f} kcal/mol")
        
    finally:
        if os.path.exists(output_path_new):
            os.unlink(output_path_new)
    
    # 分析移动情况
    print("\n" + "="*70)
    print("移动情况详细分析")
    print("="*70)
    
    # 使用pdbfixer处理PDB
    import pdbfixer
    fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    
    out_handle = io.StringIO()
    openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
    pdb = openmm_app.PDBFile(io.StringIO(out_handle.getvalue()))
    
    backbone_atoms = ['N', 'CA', 'C', 'O']
    
    # 分析新API的结果
    print("\n[新API: constraints=None + 约束受体backbone]")
    
    results_receptor = analyze_movements(
        info_new['posinit'], info_new['pos'], pdb, 'A', backbone_atoms
    )
    results_ligand = analyze_movements(
        info_new['posinit'], info_new['pos'], pdb, 'B', backbone_atoms
    )
    
    print(f"\n  受体A链 (约束backbone):")
    print(f"    主链:")
    print(f"      - 移动原子数: {results_receptor['backbone']['moved']}/{results_receptor['backbone']['total']}")
    print(f"      - 最大移动: {results_receptor['backbone']['max']:.6f} Å")
    print(f"      - 平均移动: {results_receptor['backbone']['avg']:.6f} Å")
    print(f"    侧链:")
    print(f"      - 移动原子数: {results_receptor['sidechain']['moved']}/{results_receptor['sidechain']['total']}")
    print(f"      - 最大移动: {results_receptor['sidechain']['max']:.6f} Å")
    print(f"      - 平均移动: {results_receptor['sidechain']['avg']:.6f} Å")
    
    print(f"\n  配体B链:")
    print(f"    - 移动原子数: {results_ligand['all']['moved']}/{results_ligand['all']['total']}")
    print(f"    - 最大移动: {results_ligand['all']['max']:.6f} Å")
    print(f"    - 平均移动: {results_ligand['all']['avg']:.6f} Å")
    
    # 验证期望的行为
    print("\n" + "="*70)
    print("期望行为验证")
    print("="*70)
    
    checks = []
    
    # 检查1: 配体应该移动
    ligand_moved = results_ligand['all']['max'] > 0.1
    checks.append(("配体移动", ligand_moved, results_ligand['all']['max']))
    
    # 检查2: 受体侧链应该移动
    sidechain_moved = results_receptor['sidechain']['max'] > 0.1
    checks.append(("受体侧链移动", sidechain_moved, results_receptor['sidechain']['max']))
    
    # 检查3: 受体主链应该保持不动（或移动很小）
    backbone_fixed = results_receptor['backbone']['max'] < 0.5
    checks.append(("受体主链不动", backbone_fixed, results_receptor['backbone']['max']))
    
    # 检查4: 能量应该降低
    energy_improved = (info_new['efinal'] - info_new['einit']) < 0
    checks.append(("能量降低", energy_improved, info_new['efinal'] - info_new['einit']))
    
    print("\n验证结果:")
    all_passed = True
    for name, passed, value in checks:
        status = "✅ 通过" if passed else "❌ 失败"
        print(f"  {name}: {status} (值: {value:.4f})")
        if not passed:
            all_passed = False
    
    print("\n" + "="*70)
    if all_passed:
        print("🎉 所有测试通过！新API工作正常。")
    else:
        print("⚠️ 部分测试失败，需要进一步调试。")
    print("="*70)
    
    return all_passed


def test_different_stiffness():
    """测试不同的刚度值"""
    print("\n" + "="*70)
    print("测试不同刚度值的影响")
    print("="*70)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    stiffness_values = [10**4, 10**6, 10**8, 10**10]
    backbone_atoms = ['N', 'CA', 'C', 'O']
    
    import pdbfixer
    fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(pdb_str))
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=0)
    fixer.addMissingHydrogens()
    
    out_handle = io.StringIO()
    openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
    pdb = openmm_app.PDBFile(io.StringIO(out_handle.getvalue()))
    
    print(f"\n{'刚度值':<12} {'能量变化(kcal/mol)':<20} {'主链最大移动(Å)':<18} {'配体最大移动(Å)':<18}")
    print("-"*70)
    
    for stiffness in stiffness_values:
        minimizer = ForceFieldMinimizer(
            stiffness=stiffness,
            max_iterations=1000,
            tolerance=10,
            platform='CUDA',
            constraints=None
        )
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            pdb_min, info = minimizer(
                pdb_str, output_path,
                restrain_chain=['A'],
                restrain_backbone=True
            )
            
            energy_change = info['efinal'] - info['einit']
            
            results_receptor = analyze_movements(
                info['posinit'], info['pos'], pdb, 'A', backbone_atoms
            )
            results_ligand = analyze_movements(
                info['posinit'], info['pos'], pdb, 'B', backbone_atoms
            )
            
            print(f"10^{int(np.log10(stiffness)):<10} {energy_change:<20.2f} {results_receptor['backbone']['max']:<18.6f} {results_ligand['all']['max']:<18.6f}")
            
        finally:
            if os.path.exists(output_path):
                os.unlink(output_path)


def test_comparison_with_leucine_constraint():
    """与test_leucine_constraint.py的测试进行对比"""
    print("\n" + "="*70)
    print("与CB1R测试对比")
    print("="*70)
    
    # 读取CB1R PDB
    cb1r_path = 'data_tmp/docking/CB1R_0.pdb'
    if not os.path.exists(cb1r_path):
        print("CB1R_0.pdb not found, skipping comparison...")
        return
    
    with open(cb1r_path, 'r') as f:
        cb1r_str = f.read()
    
    # 测试新API
    minimizer = ForceFieldMinimizer(
        stiffness=10**5,
        max_iterations=2000,
        tolerance=1,
        platform='CUDA',
        constraints=None
    )
    
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        output_path = tmp.name
    
    try:
        pdb_min, info = minimizer(
            cb1r_str, output_path,
            restrain_chain=['A'],
            restrain_backbone=True
        )
        
        print(f"\nCB1R测试结果 (constraints=None):")
        print(f"  初始能量: {info['einit']:.2f} kcal/mol")
        print(f"  最终能量: {info['efinal']:.2f} kcal/mol")
        print(f"  能量变化: {info['efinal'] - info['einit']:+.2f} kcal/mol")
        
        # 分析移动
        import pdbfixer
        fixer = pdbfixer.PDBFixer(pdbfile=io.StringIO(cb1r_str))
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms(seed=0)
        fixer.addMissingHydrogens()
        
        out_handle = io.StringIO()
        openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, out_handle, keepIds=True)
        pdb = openmm_app.PDBFile(io.StringIO(out_handle.getvalue()))
        
        backbone_atoms = ['N', 'CA', 'C', 'O']
        
        results = analyze_movements(
            info['posinit'], info['pos'], pdb, 'A', backbone_atoms
        )
        
        print(f"\n  受体A链主链移动:")
        print(f"    - 最大移动: {results['backbone']['max']:.6f} Å")
        print(f"    - 平均移动: {results['backbone']['avg']:.6f} Å")
        print(f"    - 移动原子数: {results['backbone']['moved']}/{results['backbone']['total']}")
        
        print(f"\n  受体A链侧链移动:")
        print(f"    - 最大移动: {results['sidechain']['max']:.6f} Å")
        print(f"    - 平均移动: {results['sidechain']['avg']:.6f} Å")
        
        # 验证结果
        backbone_ok = results['backbone']['max'] < 0.5
        sidechain_moved = results['sidechain']['max'] > 0.1
        
        print(f"\n  验证:")
        print(f"    - 主链固定: {'✅' if backbone_ok else '❌'}")
        print(f"    - 侧链移动: {'✅' if sidechain_moved else '❌'}")
        
    finally:
        if os.path.exists(output_path):
            os.unlink(output_path)


if __name__ == "__main__":
    print("="*70)
    print("测试新的ForceFieldMinimizer API")
    print("="*70)
    print("""
测试目标：
验证 constraints=None + 约束受体backbone 时：
1. 配体(Chain B)应该移动
2. 受体侧链应该移动  
3. 受体主链应该保持不动

这对于ML生成的配体优化非常重要！
    """)
    
    # 运行测试
    success = test_new_api()
    test_different_stiffness()
    test_comparison_with_leucine_constraint()
    
    print("\n" + "="*70)
    print("测试完成")
    print("="*70)
    
    if success:
        print("\n✅ 新API工作正常！")
        print("使用示例:")
        print('  minimizer = ForceFieldMinimizer(')
        print('      stiffness=10**8,')
        print('      max_iterations=1000,')
        print('      platform="CUDA",')
        print('      constraints=None  # 关键参数！')
        print('  )')
        print('  ')
        print('  # 约束受体主链，优化配体和受体侧链')
        print('  result, info = minimizer(pdb_str, output_path,')
        print('      restrain_chain=["A"],')
        print('      restrain_backbone=True')
        print('  )')
