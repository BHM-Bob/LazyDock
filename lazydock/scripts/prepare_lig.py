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

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web_utils.task import TaskPool
from pymol import cmd
from tqdm import tqdm

from lazydock.scripts._script_utils_ import Command, excute_command, process_batch_dir_lst


class smiles2pdb(Command):
    SOURCE_REPR = 'SMILES'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        # Check if RDKit is available
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError:
            put_err('RDKit is required for SMILES to PDB conversion. Please install it with: pip install rdkit', _exit=True)
        # Store RDKit modules for later use
        self.Chem = Chem
        self.AllChem = AllChem
        self.transfer_method = self.Chem.MolFromSmiles
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-i', '--input', type=str, required=True,
                          help=f'{smiles2pdb.SOURCE_REPR} string to convert to PDB format. Required.')
        args.add_argument('-o', '--output', type=str, required=True,
                          help=f'Output PDB file path. Required.')
        return args
    
    def process_args(self):
        
        # Process output path
        self.args.output = Path(self.args.output).resolve()
        if not self.args.output.parent.exists():
            self.args.output.parent.mkdir(parents=True, exist_ok=True)
    
    def main_process(self):
        try:
            # Convert SMILES to molecule
            mol = self.transfer_method(self.args.input)
            if mol is None:
                put_err(f'Failed to create molecule from {self.SOURCE_REPR}: {self.args.input}', _exit=True)
            
            # Add hydrogens
            mol = self.Chem.AddHs(mol, addResidueInfo=True)
            
            # Generate 3D coordinates
            self.AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # Minimize the structure
            self.AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            
            # Write to PDB file
            pdb_block = self.Chem.MolToPDBBlock(mol)
            with open(self.args.output, 'w') as f:
                f.write(pdb_block)
            
            self.printf(f'Successfully converted {self.SOURCE_REPR} to PDB: {self.args.output}')
        except Exception as e:
            put_err(f'Error during {self.SOURCE_REPR} to PDB conversion: {str(e)}', _exit=True)


class seq2pdb(smiles2pdb):
    SOURCE_REPR = 'Amino acid sequence'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.transfer_method = self.Chem.MolFromFASTA
        
        
class cif2pdb(Command):
    HELP = """"""
    def __init__(self, args, printf=print):
        super().__init__(args, printf, ['batch_dir'])
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--batch-dir', type=str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--main-name', type=str, required=True,
                          help='file in each sub-directory, such as model.cif.')
        args.add_argument('--new-name', type=str, default=None,
                          help='new name of the pdb, such as complex.pdb, default is %(default)s.')
        args.add_argument('--suffix', type=str, default=None,
                          help='suffix of the output pdb, such as _transfer, default is %(default)s.')

    def process_args(self):
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        if self.args.new_name:
            self.new_name_fn = self.new_name
        elif self.args.suffix:
            self.new_name_fn = self.add_suffix
        else:
            self.new_name_fn = self.only_pdb
        
    @staticmethod
    def only_pdb(cif_path: Path, *args, **kwargs):
        return cif_path.with_suffix('.pdb')
    
    @staticmethod
    def add_suffix(cif_path: Path, suffix: str, *args, **kwargs):
        return cif_path.with_suffix(f'{suffix}.pdb')
    
    @staticmethod
    def new_name(cif_path: Path, new_name: str, *args, **kwargs):
        return cif_path.with_name(f'{new_name}').with_suffix('.pdb')
        
    def main_process(self):
        # get complex paths
        cif_paths = get_paths_with_extension(self.args.batch_dir, [self.args.main_name], name_substr=self.args.main_name)
        put_log(f'get {len(cif_paths)} task(s)')
        # process each
        for cif_path in tqdm(cif_paths, total=len(cif_paths)):
            cif_path = Path(cif_path).resolve()
            cmd.reinitialize()
            cmd.set('connect_mode', 4)
            cmd.set('pdb_conect_all', 'on')
            cmd.load(str(cif_path))
            cmd.save(str(self.new_name_fn(cif_path, suffix=self.args.suffix, new_name=self.args.new_name)))


_str2func = {
    'smiles2pdb': smiles2pdb,
    'seq2pdb': seq2pdb,
    'cif2pdb': cif2pdb,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description='Prepare ligand files for docking')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))
    
    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    main()