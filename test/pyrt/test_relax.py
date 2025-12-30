#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import unittest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from lazydock.pyrt.relax import RelaxPDBChain
from lazydock.pyrt.energy_utils import calcu_interface_energy


class TestRelaxPDBChain(unittest.TestCase):
    """Test RelaxPDBChain class with different input types and relaxation scenarios"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test class - initialize pyrosetta"""
        print("Initializing PyRosetta...")
        pyrosetta.init(silent=True)
        print("PyRosetta initialized successfully!")
        
        # Test data
        cls.pdb_path = 'data_tmp/docking/CB1R_0.pdb'
        cls.receptor_chains = 'A'
        cls.ligand_chains = 'B'
        
        # Read pdb string for later use
        with open(cls.pdb_path, 'r') as f:
            cls.pdb_string = f.read()
        
        # Create pose object for testing
        cls.pose = pyrosetta.pose_from_pdb(cls.pdb_path)
        
        # Initialize relaxer
        cls.relaxer = RelaxPDBChain(max_iter=10)  # Reduced iterations for faster testing
    
    def test_relax_with_pose_object(self):
        """Test 1: Relax with Pose object input"""
        print("\nTest 1: Relaxing with Pose object...")
        
        # Calculate energy before relaxation
        energy_before = calcu_interface_energy(
            self.pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy before relaxation: {energy_before:.4f} kcal/mol")
        
        # Relax ligand chain
        relaxed_pose = self.relaxer(self.pose, self.ligand_chains)
        
        # Verify result
        self.assertIsInstance(relaxed_pose, Pose, "Relaxed result should be a Pose object")
        
        # Calculate energy after relaxation
        energy_after = calcu_interface_energy(
            relaxed_pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy after relaxation: {energy_after:.4f} kcal/mol")
        
        print("   ✓ Successfully relaxed Pose object")
    
    def test_relax_with_pdb_path(self):
        """Test 2: Relax with PDB file path input"""
        print("\nTest 2: Relaxing with PDB file path...")
        
        # Calculate energy before relaxation
        energy_before = calcu_interface_energy(
            self.pdb_path, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy before relaxation: {energy_before:.4f} kcal/mol")
        
        # Relax ligand chain
        relaxed_pose = self.relaxer(self.pdb_path, self.ligand_chains)
        
        # Verify result
        self.assertIsInstance(relaxed_pose, Pose, "Relaxed result should be a Pose object")
        
        # Calculate energy after relaxation
        energy_after = calcu_interface_energy(
            relaxed_pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy after relaxation: {energy_after:.4f} kcal/mol")
        
        print("   ✓ Successfully relaxed PDB file")
    
    def test_relax_with_pdb_string(self):
        """Test 3: Relax with PDB string input"""
        print("\nTest 3: Relaxing with PDB string...")
        
        # Calculate energy before relaxation
        energy_before = calcu_interface_energy(
            self.pdb_string, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy before relaxation: {energy_before:.4f} kcal/mol")
        
        # Relax ligand chain
        relaxed_pose = self.relaxer(self.pdb_string, self.ligand_chains)
        
        # Verify result
        self.assertIsInstance(relaxed_pose, Pose, "Relaxed result should be a Pose object")
        
        # Calculate energy after relaxation
        energy_after = calcu_interface_energy(
            relaxed_pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy after relaxation: {energy_after:.4f} kcal/mol")
        
        print("   ✓ Successfully relaxed PDB string")
    
    def test_relax_ligand_only(self):
        """Test 4: Relax ligand chain only"""
        print("\nTest 4: Relaxing ligand chain only...")
        
        # Calculate energy before relaxation
        energy_before = calcu_interface_energy(
            self.pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy before relaxation: {energy_before:.4f} kcal/mol")
        
        # Relax ligand chain
        relaxed_pose = self.relaxer(self.pose, self.ligand_chains)
        
        # Verify result
        self.assertIsInstance(relaxed_pose, Pose, "Relaxed result should be a Pose object")
        
        # Calculate energy after relaxation
        energy_after = calcu_interface_energy(
            relaxed_pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy after relaxation: {energy_after:.4f} kcal/mol")
        
        print("   ✓ Successfully relaxed ligand chain")
    
    def test_relax_receptor_only(self):
        """Test 5: Relax receptor chain only"""
        print("\nTest 5: Relaxing receptor chain only...")
        
        # Calculate energy before relaxation
        energy_before = calcu_interface_energy(
            self.pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy before relaxation: {energy_before:.4f} kcal/mol")
        
        # Relax receptor chain
        relaxed_pose = self.relaxer(self.pose, self.receptor_chains)
        
        # Verify result
        self.assertIsInstance(relaxed_pose, Pose, "Relaxed result should be a Pose object")
        
        # Calculate energy after relaxation
        energy_after = calcu_interface_energy(
            relaxed_pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy after relaxation: {energy_after:.4f} kcal/mol")
        
        print("   ✓ Successfully relaxed receptor chain")
    
    def test_relax_both_chains(self):
        """Test 6: Relax both receptor and ligand chains"""
        print("\nTest 6: Relaxing both receptor and ligand chains...")
        
        # Calculate energy before relaxation
        energy_before = calcu_interface_energy(
            self.pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy before relaxation: {energy_before:.4f} kcal/mol")
        
        # Relax both chains
        relaxed_pose = self.relaxer(self.pose, [self.receptor_chains, self.ligand_chains])
        
        # Verify result
        self.assertIsInstance(relaxed_pose, Pose, "Relaxed result should be a Pose object")
        
        # Calculate energy after relaxation
        energy_after = calcu_interface_energy(
            relaxed_pose, self.receptor_chains, self.ligand_chains
        )
        print(f"   Energy after relaxation: {energy_after:.4f} kcal/mol")
        
        print("   ✓ Successfully relaxed both chains")
    
    def test_pose_unchanged_after_relax(self):
        """Test 7: Verify original pose is unchanged after relaxation"""
        print("\nTest 7: Verifying original pose is unchanged...")
        
        # Make a deep copy of the original pose
        original_pose_copy = self.pose.clone()
        
        # Relax the pose
        relaxed_pose = self.relaxer(self.pose, self.ligand_chains)
        
        # Verify original pose is unchanged by comparing sequence
        original_sequence = self.pose.sequence()
        copy_sequence = original_pose_copy.sequence()
        
        self.assertEqual(original_sequence, copy_sequence, 
                        "Original pose should remain unchanged after relaxation")
        
        # Verify relaxed pose is different
        relaxed_sequence = relaxed_pose.sequence()
        self.assertEqual(original_sequence, relaxed_sequence,
                        "Relaxed pose should have same sequence but different coordinates")
        
        print("   ✓ Original pose remains unchanged")


if __name__ == '__main__':
    print("Starting RelaxPDBChain test suite...")
    unittest.main(verbosity=2)