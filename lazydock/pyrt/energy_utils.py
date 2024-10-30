from typing import List, Union

import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover


def calcu_interface_energy(pdb_path: str, receptor_chains: Union[str, List[str]],
                           ligand_chains: Union[str, List[str]]) -> float:
    """
    use pyrosetta to calculate the interface energy between receptor and ligand
    """
    pose = pyrosetta.pose_from_pdb(pdb_path)
    interface = ''.join(ligand_chains) + '_' + ''.join(receptor_chains)
    mover = InterfaceAnalyzerMover(interface)
    mover.set_pack_separated(True)
    mover.apply(pose)
    return pose.scores['dG_separated']