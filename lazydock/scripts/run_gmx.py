

import argparse
import os
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import Configs, put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file
from pymol import cmd
from tqdm import tqdm

from lazydock.gmx.run import Gromacs
from lazydock.scripts._script_utils_ import Command, clean_path


class simple_protein(Command):
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
    15. gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
    16. gmx mdrun -v -ntomp 4 -deffnm md -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str,
                          help='directory to store the prepared files')
        args.add_argument('-n', '--protein-name', type = str,
                          help='protein name in each sub-directory, such as protein.gro.')
        args.add_argument('--ion-mdp', type = str,
                          help='energy minimization mdp file, if is a file-path, will copy to working directory; if is a file-name, will search in working directory.')
        args.add_argument('--em-mdp', type = str,
                          help='energy minimization mdp file, if is a file-path, will copy to working directory; if is a file-name, will search in working directory.')
        args.add_argument('--nvt-mdp', type = str,
                          help='nvt mdp file, if is a file-path, will copy to working directory; if is a file-name, will search in working directory.')
        args.add_argument('--npt-mdp', type = str,
                          help='npt mdp file, if is a file-path, will copy to working directory; if is a file-name, will search in working directory.')
        args.add_argument('--md-mdp', type = str,
                          help='production md mdp file, if is a file-path, will copy to working directory; if is a file-name, will search in working directory.')
        args.add_argument('--editconf-args', type = str, default="-c -d 1.2 -bt dodecahedron",
                          help='args pass to editconf command, default is %(default)s.')
        args.add_argument('--solvate-args', type = str, default="-cs spc216.gro",
                          help='args pass to solvate command, default is %(default)s.')
        args.add_argument('--genion-args', type = str, default="-pname NA -nname CL -neutral",
                          help='args pass to genion command, default is %(default)s.')
        args.add_argument('--mdrun-args', type = str, default="-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu",
                          help='args pass to mdrun command, default is %(default)s.')                          
        return args
    
    def process_args(self):
        self.args.dir = clean_path(self.args.dir)
        
    def get_mdp(self, working_dir: Path):
        mdps = {}
        for name in ['ion', 'em', 'nvt', 'npt','md']:
            mdp_file = getattr(self.args, f'{name}_mdp')
            if os.path.isfile(mdp_file):
                shutil.copy(mdp_file, working_dir)
                mdps[name] = Path(mdp_file).name
            elif os.path.isfile(working_dir / mdp_file):
                mdps[name] = (working_dir / mdp_file).name
            else:
                put_err(f'can not find {name} mdp file: {mdp_file} in {mdp_file} or {working_dir}, exit.', _exit=True)
        return mdps
        
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
            main_name = protein_path.stem
            gmx = Gromacs(working_dir=str(protein_path.parent))
            mdps = self.get_mdp(protein_path.parent)
            # STEP 1: editconf -f protein.gro -o protein_newbox.gro -c -d 1.0 -bt cubic
            gmx.run_command_with_expect(f'editconf {self.args.editconf_args}', f=f'{main_name}.gro', o=f'{main_name}_newbox.gro')
            # STEP 2: solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top
            gmx.run_command_with_expect(f'solvate {self.args.solvate_args}', cp=f'{main_name}_newbox.gro', o=f'{main_name}_solv.gro', p='topol.top')
            # STEP 3: grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr
            gmx.run_command_with_expect('grompp', f=mdps['ion'], c=f'{main_name}_solv.gro', p='topol.top', o='ions.tpr')
            # STEP 4: genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
            gmx.run_command_with_expect(f'genion {self.args.genion_args}', s='ions.tpr', o=f'{main_name}_solv_ions.gro', p='topol.top',
                                        expect_actions=[{'Select a group:': '13\r', 'No ions to add': '', '\\timeout': ''}],
                                        expect_settings={'start_timeout': 600})
            # STEP 5: grompp -f minim.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr
            gmx.run_command_with_expect('grompp', f=mdps['em'], c=f'{main_name}_solv_ions.gro', p='topol.top', o='em.tpr')
            # STEP 6: mdrun -v -deffnm em
            gmx.run_command_with_expect(f'mdrun', deffnm='em')
            # STEP 7: energy -f em.edr -o potential.xvg
            gmx.run_command_with_expect('energy', f='em.edr', o='potential.xvg',
                                        expect_actions=[{'T-rest': '11 0\r'}])
            os.system(f'cd "{protein_path.parent}" && dit xvg_show -f potential.xvg -o potential.png -smv')
            # STEP 8: grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
            gmx.run_command_with_expect('grompp', f=mdps['nvt'], c='em.gro', r='em.gro', p='topol.top', o='nvt.tpr')
            # STEP 9: mdrun -deffnm nvt
            gmx.run_command_with_expect('mdrun', deffnm='nvt')
            # STEP 10: energy -f nvt.edr -o temperature.xvg
            gmx.run_command_with_expect('energy', f='nvt.edr', o='temperature.xvg',
                                        expect_actions=[{'Lamb-non-Protein': '16 0\r'}])
            os.system(f'cd "{protein_path.parent}" && dit xvg_show -f temperature.xvg -o temperature.png -smv')
            # STEP 11: grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
            gmx.run_command_with_expect('grompp', f=mdps['npt'], c='nvt.gro', r='nvt.gro', t='nvt.cpt', p='topol.top', o='npt.tpr')
            # STEP 12: mdrun -deffnm npt
            gmx.run_command_with_expect('mdrun', deffnm='npt')
            # STEP 13: energy -f npt.edr -o pressure.xvg
            gmx.run_command_with_expect('energy', f='npt.edr', o='pressure.xvg',
                                        expect_actions=[{'Lamb-non-Protein': '17 0\r'}])
            os.system(f'cd "{protein_path.parent}" && dit xvg_show -f pressure.xvg -o pressure.png -smv')
            # STEP 14: energy -f npt.edr -o density.xvg
            gmx.run_command_with_expect('energy', f='npt.edr', o='density.xvg',
                                        expect_actions=[{'Lamb-non-Protein': '23 0\r'}])
            os.system(f'cd "{protein_path.parent}" && dit xvg_show -f density.xvg -o density.png -smv')
            # STEP 15: grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
            gmx.run_command_with_expect('grompp', f=mdps['md'], c='npt.gro', t='npt.cpt', p='topol.top', o='md.tpr')
            # STEP 16: mdrun -v -ntomp 4 -deffnm md -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu
            gmx.run_command_with_expect(f'mdrun {self.args.mdrun_args}', deffnm='md')


class simple_ligand(simple_protein):
    HELP = 'run ligand GROMACS simulation'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
            
    
class simple_complex(simple_ligand):
    HELP = 'run complex GROMACS simulation'
    def __init__(self, args, printf=print):
        super().__init__(args, printf)
        



_str2func = {
    'simple-protein': simple_protein,
    'simple-ligand': simple_ligand,
    'simple-complex': simple_complex,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'tools for GROMACS.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')

    simple_protein_args = simple_protein.make_args(subparsers.add_parser('simple-protein', description=simple_protein.HELP))
    simple_ligand_args = simple_ligand.make_args(subparsers.add_parser('simple-ligand', description=simple_ligand.HELP))
    simple_complex_args = simple_complex.make_args(subparsers.add_parser('simple-complex', description=simple_complex.HELP))

    args = args_paser.parse_args(sys_args)
    if args.sub_command in _str2func:
        _str2func[args.sub_command](args).excute()


if __name__ == "__main__":
    # pass
    # main(r'complex -d data_tmp/gmx/complex -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir data_tmp/gmx/charmm36-jul2022.ff'.split())
    
    main()