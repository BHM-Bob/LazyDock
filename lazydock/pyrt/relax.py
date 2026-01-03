# edited from https://github.com/luost26/diffab/blob/main/diffab/tools/relax/pyrosetta_relaxer.py and Kong X, Jia Y, Huang W, et al. Full-Atom Peptide Design with Geometric Latent Diffusion[C],2024.

import os
from typing import List, Union

import pyrosetta
from pyrosetta.rosetta import protocols
from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.select import residue_selector as selections
from pyrosetta.rosetta.core.select.movemap import (MoveMapFactory,
                                                   move_map_action)
from pyrosetta.rosetta.protocols.relax import FastRelax

from lazydock.pyrt.pose_utils import load_pose

init_options = ' '.join([
    '-mute', 'all',
    '-use_input_sc',
    '-ignore_unrecognized_res',
    '-ignore_zero_occupancy', 'false',
    '-load_PDB_components', 'false',
    '-relax:default_repeats', '2',
    '-no_fconfig',
    '-use_terminal_residues', 'true',
    '-in:file:silent_struct_type', 'binary',
])


class RelaxPDBChain:
    def __init__(self, max_iter: int = 3) :
        self.scorefxn = pyrosetta.create_score_function('ref2015')
        self.fast_relax = FastRelax()
        self.fast_relax.set_scorefxn(self.scorefxn)
        self.fast_relax.max_iter(max_iter)
        
    def __call__(self, pdb: Union[str, Pose], chain: Union[str, List[str]]):
        """Relax a PDB file or string with a given chain.
        
        Args:
            pdb (Union[str, Pose]): Path to the PDB file or a PDB string or a Pose object.
            chain (Union[str, List[str]]): Chain ID or a list of chain IDs to relax.
        
        Returns:
            pyrosetta.Pose: The relaxed pose.
        """
        # Create pose from pdb path or string
        original_pose = load_pose(pdb)
        # Create task factory
        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(operation.RestrictToRepacking()) # Only allow residues to repack No design at any position.
        # make chain selector
        chain_selectors = [selections.ChainSelector(chain) for chain in chain]
        if len(chain_selectors) == 1:
            gen_selector = chain_selectors[0]
        else:
            gen_selector = selections.OrResidueSelector(chain_selectors[0], chain_selectors[1])
            for selector in chain_selectors[2:]:
                gen_selector = selections.OrResidueSelector(gen_selector, selector)
        # make neighborhood selector
        nbr_selector = selections.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(gen_selector)
        nbr_selector.set_include_focus_in_subset(True)
        subset_selector = nbr_selector
        prevent_repacking_rlt = operation.PreventRepackingRLT()
        prevent_subset_repacking = operation.OperateOnResidueSubset(prevent_repacking_rlt, subset_selector, flip_subset =True,)
        tf.push_back(prevent_subset_repacking)
        fr = self.fast_relax
        pose = original_pose.clone()
        mmf = MoveMapFactory()
        mmf.add_bb_action(move_map_action.mm_enable, gen_selector)
        mmf.add_chi_action(move_map_action.mm_enable, subset_selector)
        mm = mmf.create_movemap_from_pose(pose)
        fr.set_movemap(mm)
        fr.set_task_factory(tf)
        fr.apply(pose)
        return pose
    
    
def relax_pdb(pdb_path: str, output_path: str, chain: Union[str, List[str]], max_iter: int = 1000):
    """Relax a PDB file with a given chain.
    
    Args:
        pdb_path (str): Path to the PDB file.
        output_path (str): Path to the output PDB file.
        chain (Union[str, List[str]]): Chain ID or a list of chain IDs to relax.
    """
    relaxer = RelaxPDBChain(max_iter)
    pose = relaxer(pdb_path, chain)
    pose.dump_pdb(output_path)
    return pose