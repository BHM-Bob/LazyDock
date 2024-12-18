'''
Date: 2024-12-13 20:18:59
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-18 19:03:01
Description: 
'''

import argparse
import os
import shutil
import time
from functools import wraps
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import Configs, put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.web import TaskPool, random_sleep
from pymol import cmd
from tqdm import tqdm

from lazydock.gmx.run import Gromacs
from lazydock.gmx.thirdparty.cgenff_charmm2gmx import run_transform
from lazydock.gmx.thirdparty.sort_mol2_bonds import sort_bonds
from lazydock.pml.align_to_axis import align_pose_to_axis
from lazydock.pml.utils import get_seq
from lazydock.scripts._script_utils_ import Command, clean_path, excute_command
from lazydock.web.cgenff import get_login_browser, get_result_from_CGenFF


class prepare_complex(Command):
    HELP = """
    prepare complex for GROMACS MDS.
    - input complex.pdb should have two chains, one for receptor and one for ligand.
    - complex.pdb should already add hydrogens by Avogadro or other software.
    - complex.pdb is supposed to be aligned with the axes to save space when MDS.
    
    STEPS:
    0. center complex.pdb by obabel, align with xyz axes by lazydock.
    1. extract receptor and ligand from complex.pdb.
    2. transfer ligand.pdb to mol2 by obabel.
    3. fix ligand name in mol2 file.
    4. sort mol2 bonds by lazydock.
    5. retrive str file from CGenFF.
    6. transfer str file to top and gro file by cgenff_charmm2gmx.py
    7. make receptor.top
    8. build complex.top and complex.gro
    9. run solve, ion
    10. make select index file and add restraints to complex.top and complex.gro
    11. run last MDS
    """
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                          help='complex directory. Default is %(default)s.')
        args.add_argument('-n', '--complex-name', type = str,
                          help='complex name in each sub-directory.')
        args.add_argument('--receptor-chain-name', type = str,
                          help='receptor chain name.')
        args.add_argument('--ligand-chain-name', type = str,
                          help='ligand chain name.')
        args.add_argument('--ff-dir', type = str,
                          help='force field files directory.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        
    def main_process(self):
        if os.path.isdir(self.args.dir):
            complexs_path = get_paths_with_extension(self.args.dir, [], name_substr=self.args.complex_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.config}, exit.', _exit=True)
        put_log(f'get {len(complexs_path)} complex(s)')
        for complex_path in tqdm(complexs_path, total=len(complexs_path)):
            complex_path = Path(complex_path).resolve()
            result_dir = complex_path.parent / f'_prepare_complex_{complex_path.stem}'
            os.makedirs(result_dir, exist_ok=True)
            gmx = Gromacs(working_dir=str(complex_path.parent))
            cmd.reinitialize()
            # STEP 0.1: center complex.pdb by obabel.
            ipath, opath = str(complex_path), str(complex_path.parent / f'0.1_{complex_path.stem}_center.pdb')
            if not os.path.exists(opath):
                os.system(f'obabel -ipdb {ipath} -opdb -O {opath} -c')
            # step 0.2: align complex.pdb with xyz axes by lazydock.
            ipath, opath = opath, str(complex_path.parent / f'0.2_{complex_path.stem}_center_align_axis.pdb')
            if not os.path.exists(opath):
                cmd.load(ipath, 'complex')
                align_pose_to_axis('complex')
                cmd.h_add('complex')
                cmd.save(opath, 'complex')
                cmd.reinitialize()
            # STEP 1: extract receptor and ligand from complex.pdb.
            success = True
            ipath, opath_r, opath_l = opath, str(complex_path.parent / f'1_{complex_path.stem}_receptor.pdb'), str(complex_path.parent / f'1_{complex_path.stem}_ligand.pdb')
            if not os.path.exists(opath_r) or not os.path.exists(opath_l):
                cmd.load(ipath, 'complex')
                for mol, opath in zip(['receptor', 'ligand'], [opath_r, opath_l]):
                    chain = getattr(self.args, f'{mol}_chain_name')
                    if cmd.select(mol, f'complex and chain {chain}') == 0:
                        put_err(f'{mol} chain {chain} has zero atom in {complex_path}, skip this complex.')
                        success = False
                    else:
                        cmd.save(opath, mol)
            if not success:
                continue
            # STEP 2: transfer ligand.pdb to mol2 by obabel.
            ipath, opath = opath_l, str(complex_path.parent / f'2_{complex_path.stem}_ligand.mol2')
            if not os.path.exists(opath):
                os.system(f'obabel -ipdb {ipath} -omol2 -O {opath}')
            # STEP 3: fix ligand name and residue name in mol2 file.
            ipath, opath = opath, str(complex_path.parent / f'3_{complex_path.stem}_ligand_named.mol2')
            if not os.path.exists(opath):
                lines = opts_file(ipath, way='lines')
                lines[1] = 'LIG\n'
                get_idx_fn = lambda content: lines.index(list(filter(lambda x: x.startswith(content), lines))[0])
                atom_st, atom_ed = get_idx_fn('@<TRIPOS>ATOM')+1, get_idx_fn('@<TRIPOS>BOND')
                for i in range(atom_st, atom_ed):
                    resn = lines[i].split()[7].strip()
                    resn_idx = lines[i].index(resn)
                    lines[i] = lines[i][:resn_idx] + 'LIG' + ' '*(min(0, len(resn) - 3)) + lines[i][resn_idx+len(resn):]
                opts_file(opath, 'w', way='lines', data=lines)
            # STEP 4: sort mol2 bonds by lazydock.
            ipath, opath = opath, str(complex_path.parent / f'4_{complex_path.stem}_ligand_sorted.mol2')
            if not os.path.exists(opath):
                opts_file(opath, 'w', data=sort_bonds(ipath))
            # STEP 5: retrive str file from CGenFF.
            ipath, opath_str, opath_mol2 = opath, str(complex_path.parent / f'5_{complex_path.stem}_ligand_sorted.str'), str(complex_path.parent / f'5_{complex_path.stem}_ligand_sorted.mol2')
            cgenff_path = Path(ipath).with_suffix('.zip')
            if not os.path.exists(cgenff_path):
                pass
            if not os.path.exists(opath_str) or not os.path.exists(opath_mol2):
                for file_name, content in opts_file(cgenff_path, 'r', way='zip').items():
                    opts_file(cgenff_path.parent / file_name.replace('4_', '5_'), 'wb', data=content)
            # STEP 6: transfer str file to top and gro file by cgenff_charmm2gmx.py
            ipath_str, ipath_mol2, opath_itp = opath_str, opath_mol2, str(complex_path.parent / f'lig.itp')
            ff_dir = complex_path.parent / Path(self.args.ff_dir).name
            if not ff_dir.exists():
                shutil.copytree(os.path.abspath(self.args.ff_dir), ff_dir, dirs_exist_ok=True)
            if not os.path.exists(opath_itp):
                run_transform('LIG', ipath_mol2, ipath_str, r"Z:\USERS\BHM\PROGS\MPN12\pdb\ligand_MDS\MPN12\DOR\LazyDock\charmm36.ff")
            # STEP 7: Prepare the Protein Topology
            ipath, opath_gro = opath_r, str(complex_path.parent / f'{complex_path.stem}_receptor.gro')
            if not os.path.exists(opath_gro):
                receptor_n_term = '1' if get_seq(opath_r, fasta=False)[0] == 'P' else '0'
                if receptor_n_term == '1':
                    put_log(f'using NH2 as N ternimal because the first residue of receptor is PRO.')
                gmx.run_command_with_expect(f'gmx pdb2gmx -f {Path(ipath).name} -o {Path(opath_gro).name} -ter', [{')': '1\r'}, {'None': '1\r'}, {'None': f'{receptor_n_term}\r'}, {'None': '1\r'}])
            # STEP 8: Prepare the Complex Topology
            ipath_r, ipath_l, opath_t, opath_g = opath_r, opath, str(complex_path.parent / f'{complex_path.stem}_complex.top'), str(complex_path.parent / f'{complex_path.stem}_complex.gro')
            if not os.path.exists(opath_t) or not os.path.exists(opath_g):
                pass
            # STEP 9: run solve, ion
            ipath_t, ipath_g, opath_t, opath_g = opath_t, opath_g, str(complex_path.parent / f'{complex_path.stem}_complex_solv.top'), str(complex_path.parent / f'{complex_path.stem}_complex_solv.gro')
            if not os.path.exists(opath_t) or not os.path.exists(opath_g):
                pass
            # STEP 10: make select index file and add restraints to complex.top and complex.gro
            ipath_t, ipath_g, opath_t, opath_g = opath_t, opath_g, str(complex_path.parent / f'{complex_path.stem}_complex_solv_index.ndx'), str(complex_path.parent / f'{complex_path.stem}_complex_solv_restr.top'), str(complex_path.parent / f'{complex_path.stem}_complex_solv_restr.gro')
            if not os.path.exists(opath_t) or not os.path.exists(opath_g):
                pass
            # STEP 11: run last MDS
            ipath_t, ipath_g, opath_t, opath_g = opath_t, opath_g, str(complex_path.parent / f'{complex_path.stem}_complex_solv_restr_index.ndx'), str(complex_path.parent / f'{complex_path.stem}_complex_solv_restr_final.top'), str(complex_path.parent / f'{complex_path.stem}_complex_solv_restr_final.gro')
            if not os.path.exists(opath_t) or not os.path.exists(opath_g):
                pass
            

_str2func = {
    'prepare-complex': prepare_complex,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    prepare_complex_args = prepare_complex.make_args(subparsers.add_parser('prepare-complex', description=prepare_complex.HELP))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    # pass
    main(r'prepare-complex -d data_tmp/gmx/complex -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir data_tmp/gmx/charmm36-jul2022.ff'.split())
    
    main()