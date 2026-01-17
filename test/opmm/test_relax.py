#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import os
import sys
import tempfile
import unittest
from io import StringIO
from openmm import app as openmm_app
# Add the project root to Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from lazydock.opmm.relax import ForceFieldMinimizer


class TestForceFieldMinimizer(unittest.TestCase):
    """Test ForceFieldMinimizer class with different input types and relaxation scenarios"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test class - prepare test data"""
        print("Preparing test data for ForceFieldMinimizer...")
        
        # Test data
        cls.pdb_path = 'data_tmp/docking/CB1R_0.pdb'
        
        # Read pdb string for later use
        with open(cls.pdb_path, 'r') as f:
            cls.pdb_string = f.read()
        
        # Initialize minimizer
        cls.minimizer = ForceFieldMinimizer(max_iterations=10, platform='CUDA')  # Use CUDA for testing, reduced iterations
        
        print("Test data prepared successfully!")
    
    def test_relax_with_pdb_path(self):
        """Test 1: Relax with PDB file path input"""
        print("\nTest 1: Relaxing with PDB file path...")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Relax using file path
            result = self.minimizer(self.pdb_path, output_path)
            
            # Verify result
            self.assertIsInstance(result, tuple, "Result should be a tuple")
            self.assertEqual(len(result), 2, "Result should contain 2 elements")
            
            relaxed_pdb_str, info = result
            self.assertIsInstance(relaxed_pdb_str, str, "Relaxed PDB should be a string")
            self.assertIsInstance(info, dict, "Info should be a dictionary")
            
            # Verify energy change (allow small increase due to numerical precision)
            energy_diff = info['efinal'] - info['einit']
            self.assertLessEqual(energy_diff, 1.0, "Energy should not increase by more than 1 kcal/mol")
            
            # Verify output file exists and contains data
            self.assertTrue(os.path.exists(output_path), "Output PDB file should exist")
            self.assertTrue(os.path.getsize(output_path) > 0, "Output PDB file should not be empty")
            
            print(f"   Initial energy: {info['einit']:.4f} kcal/mol")
            print(f"   Final energy: {info['efinal']:.4f} kcal/mol")
            print("   ✓ Successfully relaxed using PDB file path")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_relax_with_pdb_string(self):
        """Test 2: Relax with PDB string input"""
        print("\nTest 2: Relaxing with PDB string...")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Relax using PDB string
            result = self.minimizer(self.pdb_string, output_path)
            
            # Verify result
            self.assertIsInstance(result, tuple, "Result should be a tuple")
            self.assertEqual(len(result), 2, "Result should contain 2 elements")
            
            relaxed_pdb_str, info = result
            self.assertIsInstance(relaxed_pdb_str, str, "Relaxed PDB should be a string")
            self.assertIsInstance(info, dict, "Info should be a dictionary")
            
            # Verify energy change (allow small increase due to numerical precision)
            energy_diff = info['efinal'] - info['einit']
            self.assertLessEqual(energy_diff, 1.0, "Energy should not increase by more than 1 kcal/mol")
            
            print(f"   Initial energy: {info['einit']:.4f} kcal/mol")
            print(f"   Final energy: {info['efinal']:.4f} kcal/mol")
            print("   ✓ Successfully relaxed using PDB string")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_relax_with_different_platforms(self):
        """Test 3: Relax with different platforms"""
        print("\nTest 3: Testing different platforms...")
        
        # Test with CPU platform
        cpu_minimizer = ForceFieldMinimizer(max_iterations=5, platform='CPU')
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Relax using CPU
            result_cpu = cpu_minimizer(self.pdb_string, output_path)
            self.assertIsInstance(result_cpu, tuple, "CPU relaxation should return a tuple")
            
            print("   ✓ CPU platform works correctly")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_relax_with_different_stiffness(self):
        """Test 4: Relax with different stiffness values"""
        print("\nTest 4: Testing different stiffness values...")
        
        # Test with low stiffness
        low_stiffness_minimizer = ForceFieldMinimizer(stiffness=1.0, max_iterations=5, platform='CUDA')
        
        # Test with high stiffness
        high_stiffness_minimizer = ForceFieldMinimizer(stiffness=100.0, max_iterations=5, platform='CUDA')
        
        # Create temporary output files
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_low:
            output_path_low = tmp_low.name
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp_high:
            output_path_high = tmp_high.name
        
        try:
            # Relax with low stiffness
            result_low, info_low = low_stiffness_minimizer(self.pdb_string, output_path_low)
            
            # Relax with high stiffness
            result_high, info_high = high_stiffness_minimizer(self.pdb_string, output_path_high)
            
            # Verify both worked
            self.assertIsInstance(result_low, str, "Low stiffness relaxation should return a string")
            self.assertIsInstance(result_high, str, "High stiffness relaxation should return a string")
            
            print(f"   Low stiffness energy change: {info_low['einit'] - info_low['efinal']:.4f} kcal/mol")
            print(f"   High stiffness energy change: {info_high['einit'] - info_high['efinal']:.4f} kcal/mol")
            print("   ✓ Different stiffness values work correctly")
            
        finally:
            # Clean up
            if os.path.exists(output_path_low):
                os.unlink(output_path_low)
            if os.path.exists(output_path_high):
                os.unlink(output_path_high)
    
    def test_relax_without_return_info(self):
        """Test 5: Relax without returning info"""
        print("\nTest 5: Relaxing without returning info...")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Relax without returning info
            result = self.minimizer(self.pdb_string, output_path, return_info=False)
            
            # Verify result is just the PDB string
            self.assertIsInstance(result, str, "Result should be a string when return_info=False")
            
            print("   ✓ Successfully relaxed without returning info")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_fix_pdb(self):
        """Test 6: Test PDB fixing functionality"""
        print("\nTest 6: Testing PDB fixing functionality...")
        
        # Test the _fix method (this is a protected method, but we'll test it for completeness)
        fixed_pdb_str = self.minimizer._fix(self.pdb_string)
        
        # Verify the fixed PDB string is different from the original
        self.assertIsInstance(fixed_pdb_str, str, "Fixed PDB should be a string")
        self.assertNotEqual(fixed_pdb_str, self.pdb_string, "Fixed PDB should be different from original")
        
        print("   ✓ PDB fixing functionality works correctly")
    
    def test_energy_remarks_in_output(self):
        """Test 7: Verify energy remarks are added to output PDB"""
        print("\nTest 7: Verifying energy remarks in output...")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Relax the PDB
            relaxed_pdb_str, info = self.minimizer(self.pdb_string, output_path)
            
            # Check if energy remarks are present
            self.assertIn("REMARK   1  INITIAL ENERGY:", relaxed_pdb_str, "Initial energy remark should be present")
            self.assertIn("REMARK   1  FINAL ENERGY:", relaxed_pdb_str, "Final energy remark should be present")
            
            print("   ✓ Energy remarks are correctly added to output PDB")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)


    def test_restrain_chain_position(self):
        """Test 8: Restrain a chain and compare positions before and after relaxation"""
        print("\nTest 8: Restraining chain and comparing positions...")
        
        # Create a minimizer with very high stiffness to effectively freeze the chain
        high_stiffness_minimizer = ForceFieldMinimizer(stiffness=10**10, max_iterations=10, platform='CUDA')
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Restrain chain A (receptor) and relax
            relaxed_pdb_str, info = high_stiffness_minimizer(
                self.pdb_string, output_path, restrain_chain=['A']
            )
            
            # Compare positions of chain A atoms before and after relaxation
            pos_before = info['posinit']
            pos_after = info['pos']
            
            # Read the original PDB to get atom chain information
            pdb = openmm_app.PDBFile(io.StringIO(self.pdb_string))
            
            # Check positions for chain A atoms
            max_diff = 0.0
            for i, atom in enumerate(pdb.topology.atoms()):
                if atom.residue.chain.id == 'A' and atom.element.name != 'hydrogen':
                    # Calculate position difference
                    diff = ((pos_after[i][0] - pos_before[i][0])**2 + 
                            (pos_after[i][1] - pos_before[i][1])**2 + 
                            (pos_after[i][2] - pos_before[i][2])**2)**0.5
                    max_diff = max(max_diff, diff)
                    
                    # Position should be almost unchanged (within 1e-4 angstrom)
                    self.assertLessEqual(diff, 1e-4, 
                        f"Chain A atom {atom.name} in residue {atom.residue.name} {atom.residue.id} "
                        f"moved by {diff:.6f} angstrom, which exceeds the tolerance of 1e-4 angstrom"
                    )
            
            print(f"   ✓ Chain A positions are preserved (max difference: {max_diff:.6f} angstrom)")
            
            # Check positions for chain B atoms (should have moved)
            b_atom_moved = False
            max_b_diff = 0.0
            avg_b_diff = 0.0
            b_atom_count = 0
            
            for i, atom in enumerate(pdb.topology.atoms()):
                if atom.residue.chain.id == 'B' and atom.element.name != 'hydrogen':
                    # Calculate position difference
                    diff = ((pos_after[i][0] - pos_before[i][0])**2 + 
                            (pos_after[i][1] - pos_before[i][1])**2 + 
                            (pos_after[i][2] - pos_before[i][2])**2)**0.5
                    max_b_diff = max(max_b_diff, diff)
                    avg_b_diff += diff
                    b_atom_count += 1
                    
                    if diff > 0.000001:  # If any atom moved by more than 0.000001 angstrom (1e-6 Å)
                        b_atom_moved = True
            
            if b_atom_count > 0:
                avg_b_diff /= b_atom_count
                print(f"   ✓ Chain B atoms moved (max difference: {max_b_diff:.6f} angstrom, average difference: {avg_b_diff:.6f} angstrom)")
            
            # At least some chain B atoms should have moved
            self.assertTrue(b_atom_moved, "Chain B atoms should have moved during relaxation")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)
    
    def test_restrain_backbone_position(self):
        """Test 9: Restrain a chain's backbone and compare positions before and after relaxation"""
        print("\nTest 9: Restraining backbone and comparing positions...")
        
        # Create a minimizer with high stiffness to effectively freeze the backbone
        high_stiffness_minimizer = ForceFieldMinimizer(stiffness=10**10, max_iterations=10,
                                                       platform='CUDA', constraints=None)
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            output_path = tmp.name
        
        try:
            # Restrain chain B's backbone and relax
            relaxed_pdb_str, info = high_stiffness_minimizer(
                self.pdb_string, output_path, restrain_chain=['B'], restrain_backbone=True
            )
            
            # Compare positions before and after relaxation
            pos_before = info['posinit']
            pos_after = info['pos']
            
            # Read the original PDB to get atom chain information
            pdb = openmm_app.PDBFile(io.StringIO(self.pdb_string))
            
            # Check positions for chain B backbone atoms
            backbone_atoms = ['N', 'CA', 'C', 'O']
            max_backbone_diff = 0.0
            non_backbone_atoms_moved = False
            max_non_backbone_diff = 0.0
            avg_non_backbone_diff = 0.0
            non_backbone_atom_count = 0
            
            for i, atom in enumerate(pdb.topology.atoms()):
                if atom.residue.chain.id == 'B':
                    # Calculate position difference
                    diff = ((pos_after[i][0] - pos_before[i][0])**2 + 
                            (pos_after[i][1] - pos_before[i][1])**2 + 
                            (pos_after[i][2] - pos_before[i][2])**2)**0.5
                    
                    if atom.name in backbone_atoms and atom.element.name != 'hydrogen':
                        # Backbone atoms should be almost unchanged
                        max_backbone_diff = max(max_backbone_diff, diff)
                        
                        # Position should be almost unchanged (within 1e-4 angstrom)
                        self.assertLessEqual(diff, 1e-4, 
                            f"Chain B backbone atom {atom.name} in residue {atom.residue.name} {atom.residue.id} "
                            f"moved by {diff:.6f} angstrom, which exceeds the tolerance of 1e-4 angstrom"
                        )
                    elif atom.element.name != 'hydrogen':
                        # Non-backbone atoms should have moved
                        max_non_backbone_diff = max(max_non_backbone_diff, diff)
                        avg_non_backbone_diff += diff
                        non_backbone_atom_count += 1
                        
                        if diff > 0.01:  # If any non-backbone atom moved by more than 0.01 angstrom
                            non_backbone_atoms_moved = True
            
            print(f"   ✓ Chain B backbone positions are preserved (max difference: {max_backbone_diff:.6f} angstrom)")
            
            if non_backbone_atom_count > 0:
                avg_non_backbone_diff /= non_backbone_atom_count
                print(f"   ✓ Chain B non-backbone atoms moved (max difference: {max_non_backbone_diff:.6f} angstrom, average difference: {avg_non_backbone_diff:.6f} angstrom)")
            
            # At least some non-backbone atoms should have moved
            self.assertTrue(non_backbone_atoms_moved, "Chain B non-backbone atoms should have moved during relaxation")
            
        finally:
            # Clean up
            if os.path.exists(output_path):
                os.unlink(output_path)


if __name__ == '__main__':
    print("Starting ForceFieldMinimizer test suite...")
    unittest.main(verbosity=2)