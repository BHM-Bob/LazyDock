'''
Date: 2026-02-18 22:36:48
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2026-02-26 19:37:12
Description: steps most from http://www.mdtutorials.com/gmx/umbrella
'''
import argparse
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from pymol import cmd
from tqdm import tqdm

from lazydock.gmx.prepare_ff import insert_content
from lazydock.gmx.run import Gromacs
from lazydock.scripts._script_utils_ import make_args_and_excute
from lazydock.scripts.run_gmx import simple_complex as _simple_complex
from lazydock.scripts.run_gmx import simple_protein as _simple_protein


class protein(_simple_protein):
    HELP = """
    run simple protein GROMACS simulation
    1. gmx editconf -f protein.gro -o protein_newbox.gro -c -d 1.0 -bt cubic
    2. gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top
    3. gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr
    4. gmx genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
    5. gmx grompp -f minim.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr
    6. gmx mdrun -v -deffnm em
    7. gmx energy -f em.edr -o potential.xvg # At the prompt, type "10 0" to select Potential (10); zero (0) terminates input.
    8. gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    9. gmx mdrun -deffnm nvt
    10. gmx energy -f nvt.edr -o temperature.xvg # Type "16 0" at the prompt to select the temperature of the system and exit.
    11. gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    12. gmx mdrun -deffnm npt
    13. gmx energy -f npt.edr -o pressure.xvg # Type "18 0" at the prompt to select the pressure of the system and exit.
    14. gmx energy -f npt.edr -o density.xvg # using energy and entering "24 0" at the prompt.
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        self.indexs = {}
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        _simple_protein.make_args(args)
        for action in args._actions[:]:
            if action.dest == 'md_mdp':
                action.required = False
        return args

    def main_process(self):
        # get protein paths
        if os.path.isdir(self.args.batch_dir):
            proteins_path = get_paths_with_extension(self.args.batch_dir, [], name_substr=self.args.protein_name)
        else:
            put_err(f'dir argument should be a directory: {self.args.batch_dir}, exit.', _exit=True)
        put_log(f'get {len(proteins_path)} protein(s)')
        # check mdp files
        mdp_names = ['ion', 'em', 'nvt', 'npt']
        mdp_exist = list(map(lambda x: os.path.isfile(getattr(self.args, f'{x}_mdp')), mdp_names))
        if not all(mdp_exist):
            missing_names = [n for n, e in zip(mdp_names, mdp_exist) if not e]
            put_log(f'Warning: can not find mdp files in abspath: {", ".join(missing_names)}, skip.')
        # sleep until start time
        self.sleep_until_start_time()
        # process each complex
        for protein_path in tqdm(proteins_path, total=len(proteins_path)):
            protein_path = Path(protein_path).resolve()
            main_name = protein_path.stem
            # check if md.tpr exists, if yes, skip
            if os.path.exists(protein_path.parent / 'md.tpr'):
                put_log(f'{protein_path} already done with md.tpr, skip.')
                continue
            # prepare gmx env and mdp files
            gmx = Gromacs(working_dir=str(protein_path.parent))
            mdps = self.get_mdp(protein_path.parent, mdp_names)
            # STEP 1 ~ 4: make box, solvate, ions
            self.make_box(protein_path, main_name, gmx, mdps)
            # STEP 5 ~ 7: energy minimization
            self.energy_minimization(protein_path, main_name, gmx, mdps)
            # STEP 8 ~ 14: equilibration
            self.equilibration(protein_path, main_name, gmx, mdps)
            # STEP 15 ~ 16: production md if assigned md-mdp
            if self.args.md_mdp:
                md_mdp = self.get_mdp(protein_path.parent, ['md'])
                self.production_md(protein_path, main_name, gmx, md_mdp)


class complex(_simple_complex):
    # inherit from simple_complex
    # include special equilibration method
    HELP = _simple_complex.HELP.replace('protein', 'complex')
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        _simple_complex.make_args(args)
        for action in args._actions[:]:
            if action.dest == 'md_mdp':
                action.required = False
        return args
        
    def main_process(self):
        # pass self, so will call methods of simple_complex
        # which include special equilibration process
        return protein.main_process(self)



_str2func = {
    'protein': protein,
    'complex': complex,
}


def main(sys_args: List[str] = None):
    make_args_and_excute('tools for GROMACS.', _str2func, sys_args)


if __name__ == "__main__":
    # pass
    # main(r'complex -d data_tmp/gmx/complex -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir data_tmp/gmx/charmm36-jul2022.ff'.split())
    
    main()