'''
Date: 2025-02-20 10:00:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-20 10:00:00
Description: Prepare ligand files for docking
'''

import argparse
import os
from pathlib import Path
from typing import List

from mbapy_lite.base import put_err

from lazydock.scripts._script_utils_ import Command, excute_command


class smiles2pdb(Command):
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-s', '--smiles', type=str, required=True,
                          help='SMILES string to convert to PDB format. Required.')
        args.add_argument('-o', '--output', type=str, required=True,
                          help='Output PDB file path. Required.')
        return args
    
    def process_args(self):
        # Check if RDKit is available
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            put_err('RDKit is required for SMILES to PDB conversion. Please install it with: pip install rdkit', _exit=True)
        
        # Process output path
        self.args.output = Path(self.args.output).resolve()
        if not self.args.output.parent.exists():
            self.args.output.parent.mkdir(parents=True, exist_ok=True)
        
        # Store RDKit modules for later use
        self.Chem = Chem
        self.AllChem = AllChem
    
    def main_process(self):
        try:
            # Convert SMILES to molecule
            mol = self.Chem.MolFromSmiles(self.args.smiles)
            if mol is None:
                put_err(f'Failed to create molecule from SMILES: {self.args.smiles}', _exit=True)
            
            # Add hydrogens
            mol = self.Chem.AddHs(mol)
            
            # Generate 3D coordinates
            self.AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # Minimize the structure
            self.AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            
            # Write to PDB file
            pdb_block = self.Chem.MolToPDBBlock(mol)
            with open(self.args.output, 'w') as f:
                f.write(pdb_block)
            
            self.printf(f'Successfully converted SMILES to PDB: {self.args.output}')
        except Exception as e:
            put_err(f'Error during SMILES to PDB conversion: {str(e)}', _exit=True)


_str2func = {
    'smiles2pdb': smiles2pdb,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description='Prepare ligand files for docking')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))
    
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main()