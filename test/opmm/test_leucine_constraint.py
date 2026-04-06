#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import tempfile
import os
import pdbfixer
from openmm import app as openmm_app
from openmm import openmm, unit
import unittest

import pyrosetta
from mbapy.file import opts_file
from lazydock.pyrt.pose_utils import _Pose
from lazydock.pyrt.relax import RelaxPDBChain


pyrosetta.init('-mute all')


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


def test_ca_constraint():
    """Test 1: Basic Cα constraint test with simple peptide"""
    print("\nTest 1: Basic Cα constraint test with simple peptide...")
    
    # Create a simple peptide for basic testing
    pose_obj = _Pose(seq='AAAAAAAAAA')
    relaxer = RelaxPDBChain(max_iter=10)
    pose = relaxer(pose_obj.pose, ['A'])
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=True) as temp_pdb:
        pose.dump_pdb(temp_pdb.name)
        pdb_str = opts_file(temp_pdb.name)
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    # Create force field
    force_field = openmm_app.ForceField("charmm36.xml")
    constraints = openmm_app.HBonds
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    
    # Add constraint only to Cα atom
    force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    k = 10**10 * unit.kilojoules_per_mole / (unit.nanometer**2)  # Very high stiffness
    force.addGlobalParameter("k", k)
    for p in ["x0", "y0", "z0"]:
        force.addPerParticleParameter(p)
    
    # Find Cα atom index
    ca_index = None
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.name == "CA":
            ca_index = i
            force.addParticle(i, pdb.positions[i])
            print(f"  Found Cα atom at index {i}")
    
    system.addForce(force)
    
    # Set up simulation
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    platform = openmm.Platform.getPlatformByName("CUDA")  # Use CUDA for GPU acceleration
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    # Get initial positions and energy
    state_initial = simulation.context.getState(getEnergy=True, getPositions=True)
    pos_initial = state_initial.getPositions()
    energy_initial = state_initial.getPotentialEnergy()
    
    print(f"  Initial energy: {energy_initial.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    # Perform energy minimization
    print("  Performing energy minimization...")
    tolerance = 1 * unit.kilojoules_per_mole / unit.nanometer  # Lower tolerance for better optimization
    simulation.minimizeEnergy(maxIterations=2000, tolerance=tolerance)
    
    # Get final positions and energy
    state_final = simulation.context.getState(getEnergy=True, getPositions=True)
    pos_final = state_final.getPositions()
    energy_final = state_final.getPotentialEnergy()
    
    print(f"  Final energy: {energy_final.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    print(f"  Energy change: {(energy_final - energy_initial).value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    # Check positions
    ca_moved = False
    other_moved = False
    max_ca_move = 0.0
    max_other_move = 0.0
    
    for i, atom in enumerate(pdb.topology.atoms()):
        # Calculate distance between initial and final positions
        dx = pos_final[i][0] - pos_initial[i][0]
        dy = pos_final[i][1] - pos_initial[i][1]
        dz = pos_final[i][2] - pos_initial[i][2]
        distance = (dx**2 + dy**2 + dz**2)**0.5
        distance_angstrom = distance.value_in_unit(unit.angstrom)
        
        if atom.name == "CA":
            max_ca_move = distance_angstrom
            if distance_angstrom > 1e-4:  # More than 0.0001 Å
                ca_moved = True
            print(f"  Cα atom moved: {distance_angstrom:.6f} Å")
        else:
            if distance_angstrom > max_other_move:
                max_other_move = distance_angstrom
            if distance_angstrom > 0.01:  # More than 0.01 Å to consider significant movement
                other_moved = True
            if i < 10:  # Only print first 10 non-CA atoms for brevity
                print(f"  Atom {atom.name} moved: {distance_angstrom:.6f} Å")
    
    # Verify results
    print(f"\n  Results:")
    print(f"  - Cα atom maximum movement: {max_ca_move:.6f} Å")
    print(f"  - Other atoms maximum movement: {max_other_move:.6f} Å")
    
    if not ca_moved:
        print("  ✓ Cα atom remained fixed (as expected)")
    else:
        print("  ✗ Cα atom moved (unexpected!)")
        
    if other_moved:
        print("  ✓ Other atoms moved (as expected)")
    else:
        print("  ✗ Other atoms did not move (unexpected!)")
    
    return not ca_moved and other_moved


def test_cb1r_backbone_constraint():
    """Test 2: Backbone constraint test with CB1R_0.pdb"""
    print("\nTest 2: Backbone constraint test with CB1R_0.pdb...")
    
    # Read CB1R_0.pdb
    cb1r_path = 'data_tmp/docking/CB1R_0.pdb'
    if not os.path.exists(cb1r_path):
        print(f"  ❌ CB1R_0.pdb not found at {cb1r_path}")
        return False
    
    with open(cb1r_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    # Create force field
    force_field = openmm_app.ForceField("charmm36.xml")
    constraints = openmm_app.HBonds
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    
    # Add constraints to backbone atoms of chain B
    force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    k = 10**8 * unit.kilojoules_per_mole / (unit.nanometer**2)  # Higher stiffness to better restrain backbone
    force.addGlobalParameter("k", k)
    for p in ["x0", "y0", "z0"]:
        force.addPerParticleParameter(p)
    
    # Find backbone atoms in chain B
    backbone_atoms = ['N', 'CA', 'C', 'O']
    backbone_indices = []
    non_backbone_indices = []
    
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'B':
            if atom.name in backbone_atoms and atom.element.name != 'hydrogen':
                backbone_indices.append(i)
                force.addParticle(i, pdb.positions[i])
            elif atom.element.name != 'hydrogen':
                non_backbone_indices.append(i)
    
    print(f"  Found {len(backbone_indices)} backbone atoms in chain B")
    print(f"  Found {len(non_backbone_indices)} non-backbone atoms in chain B")
    
    if len(backbone_indices) == 0 or len(non_backbone_indices) == 0:
        print("  ❌ Not enough atoms found for meaningful test")
        return False
    
    system.addForce(force)
    
    # Set up simulation
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    platform = openmm.Platform.getPlatformByName("CUDA")
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    # Get initial positions and energy
    state_initial = simulation.context.getState(getEnergy=True, getPositions=True)
    pos_initial = state_initial.getPositions()
    energy_initial = state_initial.getPotentialEnergy()
    
    print(f"  Initial energy: {energy_initial.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    # Perform energy minimization
    print("  Performing energy minimization...")
    tolerance = 10 * unit.kilojoules_per_mole / unit.nanometer
    simulation.minimizeEnergy(maxIterations=1000, tolerance=tolerance)
    
    # Get final positions and energy
    state_final = simulation.context.getState(getEnergy=True, getPositions=True)
    pos_final = state_final.getPositions()
    energy_final = state_final.getPotentialEnergy()
    
    print(f"  Final energy: {energy_final.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    print(f"  Energy change: {(energy_final - energy_initial).value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    # Check positions
    backbone_moved = False
    non_backbone_moved = False
    max_backbone_move = 0.0
    max_non_backbone_move = 0.0
    
    # Check backbone atoms
    for i in backbone_indices:
        dx = pos_final[i][0] - pos_initial[i][0]
        dy = pos_final[i][1] - pos_initial[i][1]
        dz = pos_final[i][2] - pos_initial[i][2]
        distance = (dx**2 + dy**2 + dz**2)**0.5
        distance_angstrom = distance.value_in_unit(unit.angstrom)
        
        if distance_angstrom > max_backbone_move:
            max_backbone_move = distance_angstrom
        if distance_angstrom > 0.001:  # More than 0.001 Å to consider significant movement
            backbone_moved = True
    
    # Check non-backbone atoms
    for i in non_backbone_indices:
        dx = pos_final[i][0] - pos_initial[i][0]
        dy = pos_final[i][1] - pos_initial[i][1]
        dz = pos_final[i][2] - pos_initial[i][2]
        distance = (dx**2 + dy**2 + dz**2)**0.5
        distance_angstrom = distance.value_in_unit(unit.angstrom)
        
        if distance_angstrom > max_non_backbone_move:
            max_non_backbone_move = distance_angstrom
        if distance_angstrom > 0.001:  # More than 0.001 Å to consider significant movement
            non_backbone_moved = True
    
    # Verify results
    print(f"\n  Results:")
    print(f"  - Backbone atoms maximum movement: {max_backbone_move:.6f} Å")
    print(f"  - Non-backbone atoms maximum movement: {max_non_backbone_move:.6f} Å")
    
    if not backbone_moved:
        print("  ✓ Backbone atoms remained fixed (as expected)")
    else:
        print("  ✗ Backbone atoms moved (unexpected!)")
        
    if non_backbone_moved:
        print("  ✓ Non-backbone atoms moved (as expected)")
    else:
        print("  ✗ Non-backbone atoms did not move (unexpected!)")
    
    return not backbone_moved and non_backbone_moved


def test_cb1r_sidechain_optimization():
    """Test 3: Sidechain optimization test with CB1R_0.pdb"""
    print("\nTest 3: Sidechain optimization test with CB1R_0.pdb...")
    
    # Read CB1R_0.pdb
    cb1r_path = 'data_tmp/docking/CB1R_0.pdb'
    if not os.path.exists(cb1r_path):
        print(f"  ❌ CB1R_0.pdb not found at {cb1r_path}")
        return False
    
    with open(cb1r_path, 'r') as f:
        pdb_str = f.read()
    
    fixed_pdb_str = fix_pdb(pdb_str)
    pdb = openmm_app.PDBFile(io.StringIO(fixed_pdb_str))
    
    # Create force field
    force_field = openmm_app.ForceField("charmm36.xml")
    constraints = openmm_app.HBonds
    system = force_field.createSystem(pdb.topology, constraints=constraints)
    
    # Add constraints to backbone atoms of chain A (receptor)
    force = openmm.CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    k = 10**5 * unit.kilojoules_per_mole / (unit.nanometer**2)  # Lower stiffness to allow more movement
    force.addGlobalParameter("k", k)
    for p in ["x0", "y0", "z0"]:
        force.addPerParticleParameter(p)
    
    # Find backbone atoms in chain A
    backbone_atoms = ['N', 'CA', 'C', 'O']
    backbone_indices = []
    sidechain_indices = []
    
    for i, atom in enumerate(pdb.topology.atoms()):
        if atom.residue.chain.id == 'A':  # Receptor chain
            if atom.name in backbone_atoms and atom.element.name != 'hydrogen':
                backbone_indices.append(i)
                force.addParticle(i, pdb.positions[i])
            elif atom.element.name != 'hydrogen' and atom.name not in backbone_atoms:
                sidechain_indices.append(i)
    
    print(f"  Found {len(backbone_indices)} backbone atoms in chain A")
    print(f"  Found {len(sidechain_indices)} sidechain atoms in chain A")
    
    if len(backbone_indices) == 0 or len(sidechain_indices) == 0:
        print("  ❌ Not enough atoms found for meaningful test")
        return False
    
    system.addForce(force)
    
    # Set up simulation
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    platform = openmm.Platform.getPlatformByName("CUDA")
    simulation = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    
    # Get initial positions and energy
    state_initial = simulation.context.getState(getEnergy=True, getPositions=True)
    pos_initial = state_initial.getPositions()
    energy_initial = state_initial.getPotentialEnergy()
    
    print(f"  Initial energy: {energy_initial.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    # Perform energy minimization
    print("  Performing energy minimization...")
    tolerance = 1 * unit.kilojoules_per_mole / unit.nanometer  # Lower tolerance for better optimization
    simulation.minimizeEnergy(maxIterations=2000, tolerance=tolerance)
    
    # Get final positions and energy
    state_final = simulation.context.getState(getEnergy=True, getPositions=True)
    pos_final = state_final.getPositions()
    energy_final = state_final.getPotentialEnergy()
    
    print(f"  Final energy: {energy_final.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    print(f"  Energy change: {(energy_final - energy_initial).value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    # Check positions and energy improvement
    backbone_moved = False
    sidechain_moved = False
    max_backbone_move = 0.0
    max_sidechain_move = 0.0
    
    # Check backbone atoms
    for i in backbone_indices:
        dx = pos_final[i][0] - pos_initial[i][0]
        dy = pos_final[i][1] - pos_initial[i][1]
        dz = pos_final[i][2] - pos_initial[i][2]
        distance = (dx**2 + dy**2 + dz**2)**0.5
        distance_angstrom = distance.value_in_unit(unit.angstrom)
        
        if distance_angstrom > max_backbone_move:
            max_backbone_move = distance_angstrom
        if distance_angstrom > 0.1:  # Allow some movement for backbone
            backbone_moved = True
    
    # Check sidechain atoms
    for i in sidechain_indices:
        dx = pos_final[i][0] - pos_initial[i][0]
        dy = pos_final[i][1] - pos_initial[i][1]
        dz = pos_final[i][2] - pos_initial[i][2]
        distance = (dx**2 + dy**2 + dz**2)**0.5
        distance_angstrom = distance.value_in_unit(unit.angstrom)
        
        if distance_angstrom > max_sidechain_move:
            max_sidechain_move = distance_angstrom
        if distance_angstrom > 0.05:  # Significant movement for sidechains
            sidechain_moved = True
    
    # Verify results
    print(f"\n  Results:")
    print(f"  - Backbone atoms maximum movement: {max_backbone_move:.6f} Å")
    print(f"  - Sidechain atoms maximum movement: {max_sidechain_move:.6f} Å")
    
    energy_improvement = energy_initial - energy_final
    print(f"  - Energy improvement: {energy_improvement.value_in_unit(unit.kilocalories_per_mole):.4f} kcal/mol")
    
    if max_backbone_move < 0.5:  # Backbone should move less than sidechains
        print("  ✓ Backbone atoms remained relatively fixed")
    else:
        print("  ⚠ Backbone atoms moved more than expected")
        
    if sidechain_moved and energy_improvement.value_in_unit(unit.kilocalories_per_mole) > 0.1:
        print("  ✓ Sidechain atoms moved and energy improved (optimization successful)")
        return True
    else:
        print("  ✗ Sidechain optimization did not achieve expected results")
        return False


class TestLeucineConstraint(unittest.TestCase):
    """Test suite for leucine constraint functionality"""
    
    def test_basic_ca_constraint(self):
        """Test basic Cα constraint functionality"""
        success = test_ca_constraint()
        self.assertTrue(success, "Basic Cα constraint test failed")
    
    def test_cb1r_backbone_constraint(self):
        """Test backbone constraint with CB1R_0.pdb"""
        success = test_cb1r_backbone_constraint()
        self.assertTrue(success, "CB1R backbone constraint test failed")
    
    def test_cb1r_sidechain_optimization(self):
        """Test sidechain optimization with CB1R_0.pdb"""
        success = test_cb1r_sidechain_optimization()
        self.assertTrue(success, "CB1R sidechain optimization test failed")


if __name__ == "__main__":
    print("Starting Leucine Constraint test suite...")
    
    # Run individual tests first for debugging
    print("\n=== Running individual tests ===")
    
    test1_success = test_ca_constraint()
    test2_success = test_cb1r_backbone_constraint()
    test3_success = test_cb1r_sidechain_optimization()
    
    print(f"\n=== Test Results ===")
    print(f"Test 1 (Basic Cα constraint): {'✅ PASSED' if test1_success else '❌ FAILED'}")
    print(f"Test 2 (CB1R backbone constraint): {'✅ PASSED' if test2_success else '❌ FAILED'}")
    print(f"Test 3 (CB1R sidechain optimization): {'✅ PASSED' if test3_success else '❌ FAILED'}")
    
    # Run unittest suite
    print("\n=== Running unittest suite ===")
    unittest.main(verbosity=2)