'''
Date: 2024-12-13 20:18:59
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-23 16:50:59
Description: 
'''

import argparse
import os
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import Configs, put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from mbapy_lite.web import Browser, TaskPool, random_sleep
from pymol import cmd
from tqdm import tqdm

from lazydock.config import CONFIG_FILE_PATH, GlobalConfig
from lazydock.gmx.run import Gromacs
from lazydock.gmx.thirdparty.cgenff_charmm2gmx import run_transform
from lazydock.gmx.thirdparty.sort_mol2_bonds import sort_bonds
from lazydock.pml.align_to_axis import align_pose_to_axis
from lazydock.pml.utils import get_seq
from lazydock.scripts._script_utils_ import Command, clean_path
from lazydock.web.cgenff import get_login_browser, get_result_from_CGenFF


class protein(Command):
    HELP = """
    prepare single protein for GROMACS MDS.
    
    STEPS:
    0. center complex.pdb by obabel, align with xyz axes by lazydock.
    1. prepare protein topology.
    """
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                          help='protein directory. Default is %(default)s.')
        args.add_argument('-n', '--protein-name', type = str,
                          help='protein file name in each sub-directory.')
        args.add_argument('--ff-dir', type = str,
                          help='force field files directory.')
        args.add_argument('--n-term', type = str, default='auto',
                          help='N-Term type for gmx pdb2gmx. Default is %(default)s.')
        args.add_argument('--pdb2gmx-args', type = str, default="-ter -ignh",
                          help='args pass to pdb2gmx command, default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        
    def main_process(self):
        # get protein paths
        if os.path.isdir(self.args.dir):
            proteins_path = get_paths_with_extension(self.args.dir, [], name_substr=self.args.protein_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.config}, exit.', _exit=True)
        put_log(f'get {len(proteins_path)} protein(s)')
        # process each complex
        for protein_path in tqdm(proteins_path, total=len(proteins_path)):
            protein_path = Path(protein_path).resolve()
            gmx = Gromacs(working_dir=str(protein_path.parent))
            # STEP 0.1: center complex.pdb by obabel.
            ipath, opath = str(protein_path), str(protein_path.parent / f'0.1_{protein_path.stem}_center.pdb')
            if not os.path.exists(opath):
                os.system(f'obabel -ipdb {ipath} -opdb -O {opath} -c')
            # step 0.2: align complex.pdb with xyz axes by lazydock.
            ipath, opath = opath, str(protein_path.parent / f'0.2_{protein_path.stem}_center_align_axis.pdb')
            if not os.path.exists(opath):
                cmd.load(ipath, 'protein')
                align_pose_to_axis('protein')
                cmd.save(opath, 'protein')
                cmd.reinitialize()
            # STEP 1: Prepare the Protein Topology
            ipath, opath_rgro = opath, str(protein_path.parent / f'{protein_path.stem}.gro')
            if not os.path.exists(opath_rgro):
                if self.args.n_term == 'auto':
                    receptor_n_term = '1' if get_seq(ipath, fasta=False)[0] == 'P' else '0'
                    if receptor_n_term == '1':
                        put_log(f'using NH2 as N ternimal because the first residue of receptor is PRO.')
                else:
                    receptor_n_term = self.args.n_term
                gmx.run_command_with_expect(f'pdb2gmx -f {Path(ipath).name} -o {Path(opath_rgro).name} {self.args.pdb2gmx_args}',
                                            [{'dihedrals)': '1\r'}, {'None': '1\r'}, {'None': f'{receptor_n_term}\r'}, {'None': '0\r'}])


class ligand(Command):
    HELP = """
    prepare single ligand for GROMACS MDS.
    
    STEPS:
    0. center ligand.pdb by obabel
    1. alter chain code to "Z", align with xyz axes by lazydock.
    2. transfer ligand.pdb to mol2 by obabel.
    3. fix ligand name in mol2 file.
    4. sort mol2 bonds by lazydock.
    5. retrive str file from CGenFF.
    6. transfer str file to top and gro file by cgenff_charmm2gmx.py
    7. prepare ligand topology.
    """
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                          help='protein directory. Default is %(default)s.')
        args.add_argument('-n', '--ligand-name', type = str,
                          help='protein file name in each sub-directory.')
        args.add_argument('--ff-dir', type = str,
                          help='force field files directory, if before CGenFF step, \
can be a asbpath to charmmFF dir; \
if after CGenFF step, should be charmmFF dir name in each sub-directory, \
the program will use the ff-dir in sub-directory.')
        args.add_argument('--max-step', type = int, default=7,
                          help='max step to do. Default is %(default)s.')
        args.add_argument('--disable-browser', action='store_true',
                          help='whether to disable browser for CGenFF.')
        return args

    def process_args(self):
        self.args.dir = clean_path(self.args.dir)

    @staticmethod
    def insert_content(content: str, before: str, new_content: str):
        is_file, path = False, None
        if os.path.isfile(content):
            is_file, path = True, content
            content = opts_file(content)
        idx1 = content.find(before)
        if idx1 == -1:
            return put_err(f'{before} not found, skip.')
        content = content[:idx1+len(before)] + new_content + content[idx1+len(before):]
        if is_file:
            opts_file(path, 'w', data=content)
        return content

    @staticmethod
    def fix_name_in_mol2(ipath: str, opath: str):
        lines = opts_file(ipath, way='lines')
        lines[1] = 'LIG\n'
        get_idx_fn = lambda content: lines.index(list(filter(lambda x: x.startswith(content), lines))[0])
        atom_st, atom_ed = get_idx_fn('@<TRIPOS>ATOM')+1, get_idx_fn('@<TRIPOS>BOND')
        for i in range(atom_st, atom_ed):
            resn = lines[i].split()[7].strip()
            resn_idx = lines[i].index(resn)
            lines[i] = lines[i][:resn_idx] + 'LIG' + ' '*(min(0, len(resn) - 3)) + lines[i][resn_idx+len(resn):]
        opts_file(opath, 'w', way='lines', data=lines)

    @staticmethod
    def get_login_browser(download_dir: str):
        put_log(f'getting CGenFF account from {CONFIG_FILE_PATH}')
        email, password = GlobalConfig.named_accounts['CGenFF']['email'], GlobalConfig.named_accounts['CGenFF']['password']
        if email is None or password is None:
            return put_err('CGenFF email or password not found in config file, skip.')
        return get_login_browser(email, password, download_dir=download_dir)

    @staticmethod
    def get_str_from_CGenFF(mol2_path: str, zip_path: str, browser: Browser) -> Union[str, None]:
        put_log(f'getting str file from CGenFF for {mol2_path}')
        get_result_from_CGenFF(mol2_path, b=browser)
        download_path = Path(browser.download_path) / Path(mol2_path).with_suffix('.zip').name
        shutil.move(str(download_path), zip_path)
        
    def prepare_ligand(self, ligand_path: str, main_path: Path):
        # STEP 2: transfer ligand.pdb to mol2 by obabel.
        ipath, opath = ligand_path, str(main_path.parent / f'2_{main_path.stem}_ligand.mol2')
        if self.args.max_step >= 2 and (not os.path.exists(opath)):
            os.system(f'obabel -ipdb {ipath} -omol2 -O {opath}')
        # STEP 3: fix ligand name and residue name in mol2 file.
        ipath, opath = opath, str(main_path.parent / f'3_{main_path.stem}_ligand_named.mol2')
        if self.args.max_step >= 3 and (not os.path.exists(opath)):
            self.fix_name_in_mol2(ipath, opath)
        # STEP 4: sort mol2 bonds by lazydock.
        ipath, opath = opath, str(main_path.parent / f'4_{main_path.stem}_ligand_sorted.mol2')
        if self.args.max_step >= 4 and (not os.path.exists(opath)):
            opts_file(opath, 'w', data=sort_bonds(ipath))
        # STEP 5: retrive str file from CGenFF.
        ipath, opath_str, opath_mol2 = opath, str(main_path.parent / f'5_{main_path.stem}_ligand_sorted.str'), str(main_path.parent / f'5_{main_path.stem}_ligand_sorted.mol2')
        cgenff_path = Path(ipath).with_suffix('.zip')
        if self.args.max_step >= 5 and (not os.path.exists(cgenff_path)):
            self.get_str_from_CGenFF(ipath, cgenff_path, browser=self.browser)
        if self.args.max_step >= 5 and (not os.path.exists(opath_str) or not os.path.exists(opath_mol2)):
            for file_name, content in opts_file(cgenff_path, 'r', way='zip').items():
                opts_file(cgenff_path.parent / file_name.replace('4_', '5_'), 'wb', data=content)
        # STEP 6: transfer str file to top and gro file by cgenff_charmm2gmx.py
        ipath_str, ipath_mol2, opath_itp = opath_str, opath_mol2, str(main_path.parent / f'lig.itp')
        # check and copy ff-dir
        ## if ff-dir already in each sub-directory, do not overwrite it.
        if os.path.exists(main_path.parent / self.args.ff_dir):
            ff_dir = main_path.parent / self.args.ff_dir
            put_log(f'ff-dir already exists in {ff_dir}, skip.')
        ## if ff-dir is a asbpath to charmmFF dir, copy it to sub-directory.
        elif os.path.exists(self.args.ff_dir):
            ff_dir = main_path.parent / Path(self.args.ff_dir).name
            if self.args.max_step >= 6 and (not ff_dir.exists()):
                shutil.copytree(os.path.abspath(self.args.ff_dir), ff_dir, dirs_exist_ok=True)
        else:
            ff_dir = None
            put_err(f'cannot find ff_dir: {self.args.ff_dir} in {main_path.parent}, set ff_dir to None, skip.')
        if self.args.max_step >= 6 and (not os.path.exists(opath_itp)):
            run_transform('LIG', ipath_mol2, ipath_str, ff_dir)

    def main_process(self):
        # allocate browser for CGenFF
        if not self.args.disable_browser:
            self.browser = self.get_login_browser(str(self.args.dir))
        # get ligand paths
        if os.path.isdir(self.args.dir):
            ligands_path = get_paths_with_extension(self.args.dir, [], name_substr=self.args.ligand_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.config}, exit.', _exit=True)
        put_log(f'get {len(ligands_path)} ligand(s)')
        # process each ligand
        for ligand_path in tqdm(ligands_path, total=len(ligands_path)):
            ligand_path = Path(ligand_path).resolve()
            gmx = Gromacs(working_dir=str(ligand_path.parent))
            cmd.reinitialize()
            # STEP 0: center complex.pdb by obabel.
            ipath, opath = str(ligand_path), str(ligand_path.parent / f'0_{ligand_path.stem}_center.pdb')
            if not os.path.exists(opath):
                os.system(f'obabel -ipdb {ipath} -opdb -O {opath} -c')
            # step 1: align complex.pdb with xyz axes by lazydock.
            ipath, opath_l = opath, str(ligand_path.parent / f'1_{ligand_path.stem}_center_align_axis.pdb')
            if not os.path.exists(opath_l):
                cmd.load(ipath, 'ligand')
                cmd.alter('ligand', 'chain="Z"')
                align_pose_to_axis('ligand')
                cmd.save(opath_l, 'ligand')
                cmd.reinitialize()
            # STEP 2 ~ 6: prepare ligand topology.
            self.prepare_ligand(opath_l, ligand_path)
            # STEP 7: Prepare the ligand Topology
            ipath, opath_gro = str(ligand_path.parent / f'lig_ini.pdb'), str(ligand_path.parent / f'lig.gro')
            if self.args.max_step >= 7 and (not os.path.exists(opath_gro)):
                cmd.load(ipath, 'ligand')
                cmd.alter('ligand', 'chain="Z"')
                cmd.save(ipath, 'ligand')
                # because may do not have Gromacs installed, so just try
                try:
                    gmx.run_command_with_expect(f'editconf -f lig_ini.pdb -o lig.gro')
                except Exception as e:
                    put_err(f'pdb2gmx failed: {e}, skip.')


class complex(ligand):
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
    7. prepare protein topology.
    8. prepare ligand topology.
    9. merge receptor and ligand gro into complex.gro, prepare topol.top.
    """
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.browser = None
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type = str, default='.',
                          help='complex directory. Default is %(default)s.')
        args.add_argument('-n', '--complex-name', type = str,
                          help='complex name in each sub-directory.')
        args.add_argument('--max-step', type = int, default=9,
                          help='max step to do. Default is %(default)s.')
        args.add_argument('--receptor-chain-name', type = str,
                          help='receptor chain name.')
        args.add_argument('--ligand-chain-name', type = str,
                          help='ligand chain name.')
        args.add_argument('--ff-dir', type = str,
                          help='force field files directory.')
        args.add_argument('--disable-browser', action='store_true',
                          help='whether to disable browser for CGenFF.')
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        if self.args.max_step < 1:
            return put_err('max step should be greater or equal to 1.', _exit=True)
        
    @staticmethod
    def extract_receptor_ligand(ipath: str, receptor_chain_name: str, ligand_chain_name: str, opath_r: str, opath_l: str):
        cmd.load(ipath, 'complex')
        for mol, opath, chain in zip(['receptor', 'ligand'], [opath_r, opath_l], [receptor_chain_name, ligand_chain_name]):
            if cmd.select(mol, f'complex and chain {chain}') == 0:
                put_err(f'{mol} chain {chain} has zero atom in {ipath}, skip this complex.')
                return False
            else:
                cmd.save(opath, mol)
        return True

    @staticmethod
    def prepare_complex_topol(ipath_rgro: str, ipath_lgro: str, ipath_top: str, opath_cgro: str, opath_top: str):
        # merge receptor and ligand gro into complex.gro
        receptor_gro_lines = list(filter(lambda x: len(x.strip()), opts_file(ipath_rgro, 'r', way='lines')))
        lig_gro_lines = list(filter(lambda x: len(x.strip()), opts_file(ipath_lgro, 'r', way='lines')))
        complex_gro_lines = receptor_gro_lines[:-1] + lig_gro_lines[2:-1] + receptor_gro_lines[-1:]
        complex_gro_lines[1] = f'{int(receptor_gro_lines[1]) + int(lig_gro_lines[1])}\n'
        opts_file(opath_cgro, 'w', way='lines', data=complex_gro_lines)
        # inset ligand paramters in topol.top
        topol = opts_file(ipath_top)
        topol = complex.insert_content(topol, '#include "posre.itp"\n#endif\n',
                                    '\n; Include ligand topology\n#include "lig.itp"\n')
        topol = complex.insert_content(topol, '#include "./charmm36-jul2022.ff/forcefield.itp"\n',
                                    '\n; Include ligand parameters\n#include "lig.prm"\n')
        topol += 'LIG                 1'
        opts_file(opath_top, 'w', data=topol)
        
    def main_process(self):
        # allocate browser for CGenFF
        if not self.args.disable_browser:
            self.browser = self.get_login_browser(str(self.args.dir))
        # get complex paths
        if os.path.isdir(self.args.dir):
            complexs_path = get_paths_with_extension(self.args.dir, [], name_substr=self.args.complex_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.config}, exit.', _exit=True)
        put_log(f'get {len(complexs_path)} complex(s)')
        # process each complex
        for complex_path in tqdm(complexs_path, total=len(complexs_path)):
            complex_path = Path(complex_path).resolve()
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
                cmd.save(opath, 'complex')
                cmd.reinitialize()
            # STEP 1: extract receptor and ligand from complex.pdb.
            ipath, opath_r, opath_l = opath, str(complex_path.parent / f'1_{complex_path.stem}_receptor.pdb'), str(complex_path.parent / f'1_{complex_path.stem}_ligand.pdb')
            if self.args.max_step >= 1 and (not os.path.exists(opath_r) or not os.path.exists(opath_l)):
                if not self.extract_receptor_ligand(ipath, self.args.receptor_chain_name, self.args.ligand_chain_name, opath_r, opath_l):
                    continue
            # STEP 2 ~ 6: prepare ligand topology.
            self.prepare_ligand(opath_l, complex_path)
            # STEP 7: Prepare the Protein Topology
            ipath, opath_rgro = opath_r, str(complex_path.parent / f'{complex_path.stem}_receptor.gro')
            if self.args.max_step >= 7 and (not os.path.exists(opath_rgro)):
                receptor_n_term = '1' if get_seq(opath_r, fasta=False)[0] == 'P' else '0'
                if receptor_n_term == '1':
                    put_log(f'using NH2 as N ternimal because the first residue of receptor is PRO.')
                gmx.run_command_with_expect(f'pdb2gmx -f {Path(ipath).name} -o {Path(opath_rgro).name} -ter',
                                            [{'dihedrals)': '1\r'}, {'None': '1\r'}, {'None': f'{receptor_n_term}\r'}, {'None': '0\r'}])
            # STEP 8: Prepare the Ligand Topology
            opath_lgro = str(complex_path.parent / 'lig.gro')
            if self.args.max_step >= 8 and (not os.path.exists(opath_lgro)):
                gmx.run_command_with_expect('editconf -f lig_ini.pdb -o lig.gro')
            # STEP 9: Prepare the Complex Topology
            opath_cgro, opath_top = str(complex_path.parent / 'complex.gro'), str(complex_path.parent / 'topol.top')
            if self.args.max_step >= 9 and (not os.path.exists(opath_cgro)):
                self.prepare_complex_topol(opath_rgro, opath_lgro, opath_top, opath_cgro, opath_top)


_str2func = {
    'protein': protein,
    'ligand': ligand,
    'complex': complex,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    prepare_protein_args = protein.make_args(subparsers.add_parser('protein', description=protein.HELP))
    prepare_ligand_args = ligand.make_args(subparsers.add_parser('ligand', description=ligand.HELP))
    prepare_complex_args = complex.make_args(subparsers.add_parser('complex', description=complex.HELP))

    args = args_paser.parse_args(sys_args)
    if args.sub_command in _str2func:
        _str2func[args.sub_command](args).excute()


if __name__ == "__main__":
    # pass
    # main(r'complex -d data_tmp/gmx/complex -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir data_tmp/gmx/charmm36-jul2022.ff'.split())
    
    main()