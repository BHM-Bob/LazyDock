'''
Date: 2024-11-27 17:24:03
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-31 10:01:03
Description: 
'''
import argparse
import os
from pathlib import Path
from typing import Callable, Dict, List, Set, Tuple

import pandas as pd
from mbapy_lite.base import put_err
from mbapy_lite.file import get_paths_with_extension, opts_file
from pymol import cmd
from tqdm import tqdm

if __name__ == '__main__':
    from lazydock.pml.autodock_utils import DlgFile
    from lazydock.pml.interaction_utils import SUPPORTED_MODE as pml_mode
    from lazydock.pml.interaction_utils import \
        calcu_receptor_poses_interaction as calc_fn_pml
    from lazydock.pml.ligplus_interaction import SUPPORTED_MODE as ligplus_mode
    from lazydock.pml.ligplus_interaction import \
        calcu_receptor_poses_interaction as calc_fn_ligplus
    from lazydock.pml.plip_interaction import SUPPORTED_MODE as plip_mode
    from lazydock.pml.plip_interaction import \
        calcu_receptor_poses_interaction as calc_fn_plip
    from lazydock.scripts._script_utils_ import Command, clean_path
else:
    from ..pml.autodock_utils import DlgFile
    from ..pml.interaction_utils import SUPPORTED_MODE as pml_mode
    from ..pml.interaction_utils import \
        calcu_receptor_poses_interaction as calc_fn_pml
    from ..pml.ligplus_interaction import SUPPORTED_MODE as ligplus_mode
    from ..pml.ligplus_interaction import \
        calcu_receptor_poses_interaction as calc_fn_ligplus
    from ..pml.plip_interaction import SUPPORTED_MODE as plip_mode
    from ..pml.plip_interaction import \
        calcu_receptor_poses_interaction as calc_fn_plip
    from ._script_utils_ import Command, clean_path
    
    
class simple_analysis(Command):
    METHODS: Dict[str, Tuple[Callable, List[str]]] = {'pymol': (calc_fn_pml, pml_mode),
                                                      'ligplus': (calc_fn_ligplus, ligplus_mode),
                                                      'plip': (calc_fn_plip, plip_mode)}
    def __init__(self, args: argparse.Namespace, printf=print) -> None:
        super().__init__(args, printf)
        self.tasks = []
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-r', '--receptor', type = str,
                          help="receptor pdb file path or file name, will be loaded by pymol.")
        args.add_argument('-l', '--ligand', type = str,
                          help=f"docking result path or file name, support Vina(pdbqt) and AutoDock(dlg) format.")
        args.add_argument('-bd', '--batch-dir', type = str, default=None,
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('--method', type = str, default='pymol', choices=simple_analysis.METHODS.keys(),
                          help='search input directory recursively, default is %(default)s.')
        args.add_argument('--mode', type = str, default='all',
                          help=f'interaction mode, multple modes can be separated by comma, all method support `\'all\'` model.\npymol: {",".join(pml_mode)}\nligplus: {",".join(ligplus_mode)}\nplip: {",".join(plip_mode)}')
        args.add_argument('--cutoff', type = float, default=4,
                          help='distance cutoff for interaction calculation, default is %(default)s.')
        args.add_argument('--hydrogen-atom-only', default=False, action='store_true',
                          help='only consider hydrogen bond acceptor and donor atoms, this only works when method is pymol, default is %(default)s.')
        args.add_argument('--output-style', type = str, default='receptor', choices=['receptor'],
                          help='output style\n receptor: resn resi distance')
        args.add_argument('--ref-res', type = str, default='',
                          help='reference residue name, input string shuld be like GLY300,ASP330, also support a text file contains this format string as a line.')
        return args

    def process_args(self):
        # process IO
        if self.args.batch_dir:
            self.args.batch_dir = clean_path(self.args.batch_dir)
        else:
            def check_file(path: str, name: str):
                path = clean_path(path)
                if not os.path.isfile(path):
                    put_err(f"{name} file is not a file path: {path}, exit.")
                    exit(1)
                if not os.path.exists(path):
                    put_err(f"{name} file is not exists: {path}, exit.")
                    exit(1)
                return path
            self.args.receptor = check_file(self.args.receptor, "receptor")
            self.args.ligand = check_file(self.args.ligand, "ligand")
        # check method and mode
        if ',' in self.args.mode:
            self.args.mode = [m.strip() for m in self.args.mode.split(',')]
        all_modes = set(simple_analysis.METHODS[self.args.method][1] + ['all'])
        if isinstance(self.args.mode, str) and self.args.mode not in all_modes:
            put_err(f"Unsupported mode: {self.args.mode}, supported mode: {simple_analysis.METHODS[self.args.method][1]}, exit.", _exit=True)
        elif isinstance(self.args.mode, list) and any(m not in all_modes for m in self.args.mode):
            unsuuported_mode = [m for m in self.args.mode if m not in all_modes]
            put_err(f"the mode you input has unsupported item(s): {unsuuported_mode}, supported mode: {simple_analysis.METHODS[self.args.method][1]}, exit.", _exit=True)
        # check ref_res
        split_fn = lambda x: set(map(lambda y: y.strip(), x.split(',')))
        if os.path.isfile(self.args.ref_res):
            self.args.ref_res = split_fn(opts_file(self.args.ref_res))
        elif self.args.ref_res:
            self.args.ref_res = split_fn(self.args.ref_res)
        else:
            self.args.ref_res = set()
            
        
    @staticmethod
    def output_fromater_receptor(inter_value: Dict[str, float], method: str):
        # pymol: [('receptor', '', 'GLY', '300', 'O', 2817), ('LIGAND_0', 'Z', 'UNK', '1', 'N', 74), 2.8066137153155943]
        if method == 'pymol':
            inter_value = sorted(inter_value, key=lambda x: int(x[0][3]))
            return '; '.join(f'{v[0][2]}{v[0][3]}-{v[2]:.2f}' for v in inter_value)
        elif method in {'ligplus', 'plip'}:
            inter_value = sorted(inter_value, key=lambda x: int(x[0][0]))
            return '; '.join(f'{v[0][1]}{v[0][0]}-{v[2]:.2f}' for v in inter_value)
        else:
            put_err(f"Unsupported method: {method}, exit.")
            exit(1)
            
    @staticmethod
    def calc_interaction_from_dlg(receptor_path: str, dlg_path: str, method: str, mode: List[str], cutoff: float,
                                  output_formater: Callable, hydrogen_atom_only: bool = True, ref_res: Set[str] = None) -> None:
        ref_res = ref_res or set()
        # set path
        bar = tqdm(desc='Calculating interaction', leave=False)
        root = os.path.abspath(os.path.dirname(dlg_path))
        w_dir = os.path.join(root, 'ligplus') if method == 'ligplus' else root
        # load receptor
        receptor_name = Path(receptor_path).stem
        cmd.load(receptor_path, receptor_name)
        cmd.alter(receptor_name, 'chain="A"')
        bar.set_description(f'receptor loaded from {receptor_path}')
        # load poses from dlg and perform analysis
        # load poses
        dlg = DlgFile(dlg_path, None, True, True)
        bar.set_description(f'dlg loaded from {dlg_path}')
        dlg.sort_pose() # sort by docking energy
        pose_names = []
        for i, pose in enumerate(dlg.pose_lst):
            pose_names.append(f'LIGAND_{i}')
            cmd.read_pdbstr(pose.as_pdb_string(), pose_names[-1])
            cmd.alter(pose_names[-1], 'chain="Z"')
            cmd.alter(pose_names[-1], 'type="HETATM"')
        bar.set_description(f'{len(pose_names)} pose loaded')
        # calcu interactions
        fn, _ = simple_analysis.METHODS[method]
        bar.set_description(f'performing {method} calculation')
        interactions, _ = fn(receptor_name, pose_names, mode=mode, cutoff=cutoff, verbose=True, force_cwd=True, w_dir=w_dir, hydrogen_atom_only=hydrogen_atom_only)
        if interactions is None:
            cmd.reinitialize()
            return put_err(f"No interactions found in {dlg_path}")
        if method == 'pymol':
            interactions = {k:v[-1] for k,v in interactions.items()}
        bar.set_description(f'{method} interactions calculated')
        # save interactions
        df = pd.DataFrame()
        for i, (pose, interaction) in enumerate(zip(dlg.pose_lst, interactions.values())):
            df.loc[i, 'energy'] = pose.energy
            df.loc[i, 'ref_res'] = ''
            for inter_mode, inter_value in interaction.items():
                fmt_string = output_formater(inter_value, method)
                df.loc[i, inter_mode] = fmt_string
                for r in ref_res:
                    if r in fmt_string and not r in df.loc[i,'ref_res']:
                        df.loc[i,'ref_res'] += r
        df.to_excel(os.path.join(root, f'{Path(dlg_path).stem}_{method}_interactions.xlsx'))
        bar.set_description(f'{method} interactions saved')
        # release all in pymol
        cmd.reinitialize()

    def main_process(self):
        # load origin dfs from data file
        if self.args.batch_dir:
            r_paths = get_paths_with_extension(self.args.batch_dir, ['.pdb', '.pdbqt'], name_substr=self.args.receptor)
            l_paths = get_paths_with_extension(self.args.batch_dir, ['.pdbqt', '.dlg'], name_substr=self.args.ligand)
            if len(r_paths) != len(l_paths):
                r_roots = [os.path.dirname(p) for p in r_paths]
                l_roots = [os.path.dirname(p) for p in l_paths]
                roots_count = {root: r_roots.count(root)+l_roots.count(root) for root in (set(r_roots) | set(l_roots))}
                invalid_roots = '\n'.join([root for root, count in roots_count.items() if count != 2])
                return put_err(f"The number of receptor and ligand files is not equal, please check the input files.\ninvalid roots:{invalid_roots}")
            for r_path, l_path in zip(r_paths, l_paths):
                self.tasks.append((r_path, l_path, self.args.method, self.args.mode))
        else:
            self.tasks.append((self.args.receptor, self.args.ligand, self.args.method, self.args.mode))
        # run tasks
        print(f'found {len(self.tasks)} tasks.')
        bar = tqdm(total=len(self.tasks), desc='Calculating interaction')
        for r_path, l_path, method, mode in self.tasks:
            wdir = os.path.dirname(l_path)
            bar.set_description(f"{wdir}: {os.path.basename(r_path)} and {os.path.basename(l_path)}")
            self.calc_interaction_from_dlg(r_path, l_path, method, mode, self.args.cutoff,
                                           getattr(self, f'output_fromater_{self.args.output_style}'),
                                           self.args.hydrogen_atom_only)
            bar.update(1)

_str2func = {
    'simple-analysis': simple_analysis,
}


def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser()
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    simple_analysis_args = simple_analysis.make_args(subparsers.add_parser('simple-analysis', description='perform simple analysis on docking result'))

    args = args_paser.parse_args(sys_args)
    if args.sub_command in _str2func:
        _str2func[args.sub_command](args).excute()


if __name__ == "__main__":
    # main(r'simple-analysis -r receptor.pdbqt -l dock.pdbqt -bd data_tmp/docking --method pymol --mode '.split() + ['polar contact'])
    
    main()