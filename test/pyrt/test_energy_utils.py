#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import unittest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import pyrosetta
from lazydock.pyrt.energy_utils import calcu_interface_energy


class TestCalcuInterfaceEnergy(unittest.TestCase):
    """Test calcu_interface_energy function with different parameters"""
    
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
    
    def test_with_path_default_score(self):
        """Test 1: Calculate interface energy with file path and default score function"""
        print("\nTest 1: Calculating with file path (default score function)...")
        energy = calcu_interface_energy(
            self.pdb_path, 
            self.receptor_chains, 
            self.ligand_chains
        )
        self.assertIsInstance(energy, float, "Energy should be a float")
        print(f"   ✓ Interface energy: {energy:.4f} kcal/mol")
    
    def test_with_pdb_string(self):
        """Test 2: Calculate interface energy with pdb string"""
        print("\nTest 2: Calculating with pdb string...")
        energy = calcu_interface_energy(
            self.pdb_string, 
            self.receptor_chains, 
            self.ligand_chains
        )
        self.assertIsInstance(energy, float, "Energy should be a float")
        print(f"   ✓ Interface energy: {energy:.4f} kcal/mol")
    
    def test_with_different_score_function(self):
        """Test 3: Calculate interface energy with different score function"""
        print("\nTest 3: Calculating with score12 score function...")
        energy = calcu_interface_energy(
            self.pdb_path, 
            self.receptor_chains, 
            self.ligand_chains,
            scorefxn_name='score12'
        )
        self.assertIsInstance(energy, float, "Energy should be a float")
        print(f"   ✓ Interface energy (score12): {energy:.4f} kcal/mol")
    
    def test_with_non_existent_score_function(self):
        """Test 4: Test error handling with non-existent score function"""
        print("\nTest 4: Testing with non-existent score function (expected to fail)...")
        with self.assertRaises(Exception):
            calcu_interface_energy(
                self.pdb_path, 
                self.receptor_chains, 
                self.ligand_chains,
                scorefxn_name='non_existent_scorefxn'
            )
        print("   ✓ Expected error occurred")


if __name__ == '__main__':
    print("Starting test suite...")
    unittest.main(verbosity=2)
