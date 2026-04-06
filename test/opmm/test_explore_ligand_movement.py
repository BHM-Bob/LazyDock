#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import os
import tempfile
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


def create_system(pdb, constraints=None):
    """Create a molecular system from PDB"""
    force_field = openmm_app.ForceField("charmm36.xml")
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    return system


def add_backbone_restraints(system, pdb, chain_id, stiffness=10**8):
    """Add restraints to backbone atoms of a specific chain"""
    force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", stiffness * unit.kilojoules_per_mole / (unit.nanometer**2))
    for p in ["x0", "y0", "z0"]:
        force.addPerParticleParameter(p)
    
    backbone_atoms = ['N', 'CA', 'C', 'O']
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == chain_id:
            if atom.name in backbone_atoms and atom.element.name != 'hydrogen':
                force.addParticle(i, pdb.positions[i])
    
    system.addForce(force)
    return system


def add_ligand_restraints(system, pdb, chain_id, stiffness=10**8):
    """Add restraints to all atoms of a specific chain (ligand)"""
    force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", stiffness * unit.kilojoules_per_mole / (unit.nanometer**2))
    for p in ["x0", "y0", "z0"]:
        force.addPerParticleParameter(p)
    
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == chain_id:
            force.addParticle(i, pdb.positions[i])
    
    system.addForce(force)
    return system


def run_minimization(pdb, system, platform='CUDA', max_iterations=1000, tolerance=10):
    """Run energy minimization and return results"""
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    platform_obj = openmm.Platform.getPlatformByName(platform)
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform_obj)
    simulation.context.setPositions(pdb.positions)
    
    # Get initial state
    state_initial = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True)
    
    # Perform minimization
    simulation.minimizeEnergy(maxIterations=max_iterations, 
                             tolerance=tolerance * unit.kilojoules_per_mole / unit.nanometer)
    
    # Get final state
    state_final = simulation.context.getState(getEnergy=True, getPositions=True)
    
    return {
        'initial_energy': state_initial.getPotentialEnergy(),
        'final_energy': state_final.getPotentialEnergy(),
        'initial_positions': state_initial.getPositions(asNumpy=True),
        'final_positions': state_final.getPositions(asNumpy=True),
        'initial_forces': state_initial.getForces(asNumpy=True) if platform == 'CUDA' else None
    }


def analyze_movement(pos_initial, pos_final, pdb, chain_id):
    """Analyze atom movements for a specific chain"""
    max_move = 0.0
    avg_move = 0.0
    moved_count = 0
    atom_count = 0
    
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == chain_id and atom.element.name != 'hydrogen':
            diff = pos_final[i] - pos_initial[i]
            distance = np.linalg.norm(diff)
            max_move = max(max_move, distance)
            avg_move += distance
            atom_count += 1
            if distance > 0.01:  # 0.01 Å threshold
                moved_count += 1
    
    avg_move /= atom_count if atom_count > 0 else 1
    return max_move, avg_move, moved_count, atom_count


def test_no_constraints():
    """Test 1: Optimization without any constraints"""
    print("\n" + "="*60)
    print("Test 1: 无约束优化 (No Constraints)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    results = run_minimization(pdb, system)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze receptor movement
    max_receptor, avg_receptor, moved_receptor, total_receptor = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'A'
    )
    print(f"\n受体A链:")
    print(f"  - 移动原子数: {moved_receptor}/{total_receptor}")
    print(f"  - 最大移动: {max_receptor*10:.4f} Å")
    print(f"  - 平均移动: {avg_receptor*10:.4f} Å")
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return results


def test_receptor_backbone_fixed():
    """Test 2: Fix receptor backbone and optimize"""
    print("\n" + "="*60)
    print("Test 2: 固定受体主链 (Receptor Backbone Fixed)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    system = add_backbone_restraints(system, pdb, 'A', stiffness=10**8)
    results = run_minimization(pdb, system)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze receptor backbone
    backbone_atoms = ['N', 'CA', 'C', 'O']
    backbone_moves = []
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'A' and atom.name in backbone_atoms and atom.element.name != 'hydrogen':
            diff = results['final_positions'][i] - results['initial_positions'][i]
            distance = np.linalg.norm(diff) * 10  # Convert to Å
            backbone_moves.append(distance)
    
    if backbone_moves:
        print(f"\n受体主链原子移动:")
        print(f"  - 最大移动: {max(backbone_moves):.6f} Å")
        print(f"  - 平均移动: {np.mean(backbone_moves):.6f} Å")
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return results


def test_ligand_fixed():
    """Test 3: Fix ligand and optimize receptor"""
    print("\n" + "="*60)
    print("Test 3: 固定配体 (Ligand Fixed)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    system = add_ligand_restraints(system, pdb, 'B', stiffness=10**8)
    results = run_minimization(pdb, system)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze receptor movement
    max_receptor, avg_receptor, moved_receptor, total_receptor = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'A'
    )
    print(f"\n受体A链:")
    print(f"  - 移动原子数: {moved_receptor}/{total_receptor}")
    print(f"  - 最大移动: {max_receptor*10:.4f} Å")
    print(f"  - 平均移动: {avg_receptor*10:.4f} Å")
    
    return results


def test_no_hbond_constraints():
    """Test 4: Optimization without HBond constraints"""
    print("\n" + "="*60)
    print("Test 4: 移除氢键约束 (No HBond Constraints)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    # No constraints
    system = create_system(pdb, constraints=None)
    system = add_backbone_restraints(system, pdb, 'A', stiffness=10**8)
    results = run_minimization(pdb, system)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return results


def test_lower_stiffness():
    """Test 5: Lower constraint stiffness"""
    print("\n" + "="*60)
    print("Test 5: 降低约束刚度 (Lower Stiffness: 10^6)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    system = add_backbone_restraints(system, pdb, 'A', stiffness=10**6)  # Lower stiffness
    results = run_minimization(pdb, system)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze receptor backbone
    backbone_atoms = ['N', 'CA', 'C', 'O']
    backbone_moves = []
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'A' and atom.name in backbone_atoms and atom.element.name != 'hydrogen':
            diff = results['final_positions'][i] - results['initial_positions'][i]
            distance = np.linalg.norm(diff) * 10
            backbone_moves.append(distance)
    
    if backbone_moves:
        print(f"\n受体主链原子移动:")
        print(f"  - 最大移动: {max(backbone_moves):.6f} Å")
        print(f"  - 平均移动: {np.mean(backbone_moves):.6f} Å")
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return results


def test_steepest_descent():
    """Test 6: Use steepest descent minimizer instead of default L-BFGS"""
    print("\n" + "="*60)
    print("Test 6: 使用最速下降法 (Steepest Descent Minimizer)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    system = add_backbone_restraints(system, pdb, 'A', stiffness=10**8)
    
    # Use steepest descent
    integrator = openmm.CustomIntegrator(0.0)
    integrator.addGlobalVariable("energy", 0)
    integrator.addPerDofVariable("oldPos", 0)
    integrator.addComputeGlobal("energy", "energy")
    integrator.addComputePerDof("oldPos", "x")
    integrator.addComputePerDof("x", "x - m * f / mass")
    integrator.addComputePerDof("x", "x + 0.5 * dt * v")
    integrator.addComputePerDof("v", "v + 0.5 * dt * f / mass")
    integrator.addComputePerDof("v", "v * friction")
    integrator.addComputePerDof("x", "x + 0.5 * dt * v")
    
    platform = openmm.Platform.getPlatformByName("CUDA")
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    state_initial = simulation.context.getState(getEnergy=True, getPositions=True)
    
    # Manual steepest descent steps
    for step in range(100):
        simulation.step(1)
    
    state_final = simulation.context.getState(getEnergy=True, getPositions=True)
    
    print(f"\n初始能量: {state_initial.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {state_final.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(state_final.getPotentialEnergy() - state_initial.getPotentialEnergy()).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    pos_initial = state_initial.getPositions(asNumpy=True)
    pos_final = state_final.getPositions(asNumpy=True)
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        pos_initial, pos_final, pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return {
        'initial_energy': state_initial.getPotentialEnergy(),
        'final_energy': state_final.getPotentialEnergy(),
        'initial_positions': pos_initial,
        'final_positions': pos_final
    }


def test_receptor_sidechain_only():
    """Test 7: Fix only receptor sidechains, allow backbone movement"""
    print("\n" + "="*60)
    print("Test 7: 只固定受体侧链 (Fix Sidechains Only)")
    print("="*60)
    
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    
    # Add restraints to receptor sidechains (non-backbone atoms)
    force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", 10**8 * unit.kilojoules_per_mole / (unit.nanometer**2))
    for p in ["x0", "y0", "z0"]:
        force.addPerParticleParameter(p)
    
    backbone_atoms = ['N', 'CA', 'C', 'O']
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'A':
            if atom.name not in backbone_atoms and atom.element.name != 'hydrogen':
                force.addParticle(i, pdb.positions[i])
    
    system.addForce(force)
    results = run_minimization(pdb, system)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze receptor backbone movement
    backbone_moves = []
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'A' and atom.name in backbone_atoms and atom.element.name != 'hydrogen':
            diff = results['final_positions'][i] - results['initial_positions'][i]
            distance = np.linalg.norm(diff) * 10
            backbone_moves.append(distance)
    
    if backbone_moves:
        print(f"\n受体主链原子移动:")
        print(f"  - 最大移动: {max(backbone_moves):.6f} Å")
        print(f"  - 平均移动: {np.mean(backbone_moves):.6f} Å")
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return results


def test_compare_with_cb1r():
    """Test 8: Compare with previous CB1R results"""
    print("\n" + "="*60)
    print("Test 8: 对比CB1R测试结果 (Compare with CB1R Results)")
    print("="*60)
    
    # Re-run the same test as in test_leucine_constraint.py
    pdb_path = 'data_tmp/pdb/pdbstr.pdb'
    with open(pdb_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    system = create_system(pdb, constraints=openmm_app.HBonds)
    system = add_backbone_restraints(system, pdb, 'A', stiffness=10**5)  # Lower stiffness like CB1R test
    results = run_minimization(pdb, system, max_iterations=2000, tolerance=1)
    
    print(f"\n初始能量: {results['initial_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"最终能量: {results['final_energy'].value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    print(f"能量变化: {(results['final_energy'] - results['initial_energy']).value_in_unit(unit.kilocalories_per_mole):.2f} kcal/mol")
    
    # Analyze receptor backbone
    backbone_atoms = ['N', 'CA', 'C', 'O']
    backbone_moves = []
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'A' and atom.name in backbone_atoms and atom.element.name != 'hydrogen':
            diff = results['final_positions'][i] - results['initial_positions'][i]
            distance = np.linalg.norm(diff) * 10
            backbone_moves.append(distance)
    
    if backbone_moves:
        print(f"\n受体主链原子移动:")
        print(f"  - 最大移动: {max(backbone_moves):.6f} Å")
        print(f"  - 平均移动: {np.mean(backbone_moves):.6f} Å")
    
    # Analyze ligand movement
    max_ligand, avg_ligand, moved_ligand, total_ligand = analyze_movement(
        results['initial_positions'], results['final_positions'], pdb, 'B'
    )
    print(f"\n配体B链:")
    print(f"  - 移动原子数: {moved_ligand}/{total_ligand}")
    print(f"  - 最大移动: {max_ligand*10:.4f} Å")
    print(f"  - 平均移动: {avg_ligand*10:.4f} Å")
    
    return results


if __name__ == "__main__":
    print("="*70)
    print("探索固定受体主链后配体不移动的原因")
    print("="*70)
    
    # Run all tests
    results = []
    results.append(("无约束优化", test_no_constraints))
    results.append(("固定受体主链", test_receptor_backbone_fixed))
    results.append(("固定配体", test_ligand_fixed))
    results.append(("移除氢键约束", test_no_hbond_constraints))
    results.append(("降低约束刚度", test_lower_stiffness))
    results.append(("最速下降法", test_steepest_descent))
    results.append(("只固定侧链", test_receptor_sidechain_only))
    results.append(("对比CB1R测试", test_compare_with_cb1r))
    
    print("\n" + "="*70)
    print("测试结果总结")
    print("="*70)
    print(f"\n{'测试名称':<25} {'初始能量':<15} {'最终能量':<15} {'能量变化':<15} {'配体最大移动(Å)':<15}")
    print("-"*85)
    
    for name, test_func in results:
        try:
            result = test_func()
            initial = result['initial_energy'].value_in_unit(unit.kilocalories_per_mole)
            final = result['final_energy'].value_in_unit(unit.kilocalories_per_mole)
            change = final - initial
            
            max_ligand, _, _, _ = analyze_movement(
                result['initial_positions'], result['final_positions'], 
                openmm_app.PDBFile(io.StringIO(fix_pdb(open('data_tmp/pdb/pdbstr.pdb').read()))), 'B'
            )
            
            print(f"{name:<25} {initial:<15.2f} {final:<15.2f} {change:<15.2f} {max_ligand*10:<15.4f}")
        except Exception as e:
            print(f"{name:<25} ERROR: {str(e)[:50]}")
    
    print("\n" + "="*70)
    print("分析结论")
    print("="*70)
    print("""
可能的原因分析：
1. 系统处于能量最小值 - 配体已经是优化后的构象
2. 氢键约束(HBonds)锁定了系统 - 尝试移除约束(测试4)
3. 力场对非标准残基处理 - 需要检查力场参数
4. 约束力过强 - 尝试降低刚度(测试5)
5. 最小化算法问题 - 尝试不同算法(测试6)
6. 口袋构象刚性 - 尝试允许受体主链移动(测试7)
7. 与CB1R对比 - 相同参数下测试(测试8)
""")
