import os
from typing import List, Union

import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

from lazydock.pyrt.pose_utils import load_pose


def calcu_interface_energy(pdb: str|Pose, receptor_chains: Union[str, List[str]],
                           ligand_chains: Union[str, List[str]], scorefxn_name: str = 'ref2015') -> float:
    """
    Calculate the interface energy between receptor and ligand using PyRosetta.
    
    Args:
        pdb: PDB file path or PDB string content or Pose object
        receptor_chains: Receptor chains (str or list of str)
        ligand_chains: Ligand chains (str or list of str)
        scorefxn_name: Energy function name. Available options:
            - 'ref2015': Default all-atom energy function (most recommended)
            - 'score12': Legacy all-atom energy function
            - 'ref2015_cart': Cartesian version of ref2015
            - 'beta_nov16': Requires -corrections::beta_nov16 flag
            - 'talaris2014': Requires -corrections::restore_talaris_behavior flag
    
    Returns:
        float: Interface energy (dG_separated) in kcal/mol
        
    Notes:
        The function populates pose.scores with additional valuable metrics:
        
        **Energy Metrics:**
        - dG_separated: Interface energy calculated from separated structures
        - dG_cross: Cross energy term
        - dG_cross/dSASAx100: Normalized cross energy per 100 Å² SASA
        - dG_separated/dSASAx100: Normalized interface energy per 100 Å² SASA
        - per_residue_energy_int: Average energy per interface residue
        
        **Surface Area Metrics:**
        - dSASA_int: Total solvent accessible surface area change at interface (Å²)
        - dSASA_hphobic: Hydrophobic SASA change at interface (Å²)
        - dSASA_polar: Polar SASA change at interface (Å²)
        
        **Hydrogen Bond Metrics:**
        - hbonds_int: Number of hydrogen bonds at interface
        - delta_unsatHbonds: Change in unsatisfied hydrogen bonds
        - hbond_E_fraction: Fraction of energy from hydrogen bonds
        
        **Structure Quality Metrics:**
        - sc_value: Shape complementarity score (0-1, higher is better)
        - packstat: Packing statistics (requires additional setup)
        
        **Residue Count Metrics:**
        - nres_all: Total number of residues in complex
        - nres_int: Number of interface residues
    """
    # Convert chains to list if they are strings
    if isinstance(receptor_chains, str):
        receptor_chains = [receptor_chains]
    if isinstance(ligand_chains, str):
        ligand_chains = [ligand_chains]
    
    # Create pose from pdb path or string
    pose = load_pose(pdb)
    
    # Create score function
    if scorefxn_name is None:
        scorefxn_name = 'ref2015'
    scorefxn = pyrosetta.create_score_function(scorefxn_name)
    
    # Create interface analyzer
    interface = ''.join(ligand_chains) + '_' + ''.join(receptor_chains)
    
    # Create InterfaceAnalyzerMover
    mover = InterfaceAnalyzerMover(interface, True, scorefxn)
    mover.set_pack_separated(True)
    mover.apply(pose)
    
    return pose.scores['dG_separated']


def calcu_single_energy(pdb: str|Pose, scorefxn_name: str = 'ref2015') -> float:
    """
    Calculate the energy of a single chain in a PDB file.
    
    Args:
        pdb: PDB file path or PDB string content or Pose object
        chain: Chain ID to calculate energy for
        scorefxn_name: Energy function name. Available options:
            - 'ref2015': Default all-atom energy function (most recommended)
            - 'score12': Legacy all-atom energy function
            - 'ref2015_cart': Cartesian version of ref2015
            - 'beta_nov16': Requires -corrections::beta_nov16 flag
            - 'talaris2014': Requires -corrections::restore_talaris_behavior flag
    
    Returns:
        float: Energy of the specified chain in kcal/mol
    """
    # Create score function
    scorefxn_name = scorefxn_name or 'ref2015'
    scorefxn = pyrosetta.create_score_function(scorefxn_name)
    # Create pose from pdb path or string
    pose = load_pose(pdb)
    # Calculate energy
    scorefxn(pose)
    return pose.energies().total_energy()


if __name__ == '__main__':
    # Initialize PyRosetta
    pyrosetta.init(silent=True)
    # Test the function with a PDB file
    pdb_path = 'data_tmp/docking/CB1R_0.pdb'
    receptor_chains = 'A'
    ligand_chains = 'B'
    
    # Calculate interface energy with score12 score function
    energy = calcu_interface_energy(pdb_path, receptor_chains, ligand_chains, scorefxn_name='score12')
    print(f"Interface energy (score12): {energy:.4f} kcal/mol")
    
    # Calculate single chain energy with ref2015 score function
    energy = calcu_single_energy(pdb_path, scorefxn_name='ref2015')
    print(f"Receptor energy (ref2015): {energy:.4f} kcal/mol")