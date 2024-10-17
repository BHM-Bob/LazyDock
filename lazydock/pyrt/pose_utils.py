'''
Date: 2024-05-13 15:40:33
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-05-14 16:52:08
Description: 
'''
from typing import Union

import pyrosetta
from mbapy_lite.base import get_default_call_for_None, parameter_checker
from pyrosetta.rosetta.core.chemical import ChemicalManager
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.pose import (Pose, append_pose_to_pose,
                                         get_chain_from_chain_id,
                                         get_chain_id_from_chain,
                                         renumber_pdbinfo_based_on_conf_chains)
from pyrosetta.rosetta.protocols.grafting import delete_region
from pyrosetta.rosetta.protocols.simple_moves import (
    AddPDBInfoMover, DeleteChainMover, SwitchChainOrderMover,
    SwitchResidueTypeSetMover)
from pyrosetta.toolbox import mutate_residue


class _Pose(object):
    
    add_pdbinfo_mover = AddPDBInfoMover()
    delete_chain_mover = DeleteChainMover()
    swap_chain_mover = None
    chm = ChemicalManager.get_instance()
    rts = None
    to_centroid_mover = SwitchResidueTypeSetMover('centroid')
    to_fullatom_mover = SwitchResidueTypeSetMover('full_atom')
    
    def __init__(self, seq: str = None, pdb_path: str = None, rcsb: str = None) -> None:
        # init public attributes
        if self.swap_chain_mover is None:
            self.swap_chain_mover = SwitchChainOrderMover()
        if self.rts is None:
            self.rts = self.chm.residue_type_set("fa_standard")
        # init pose from seq or pdb
        if seq is not None:
            self.pose = pyrosetta.pose_from_sequence(seq)
        elif pdb_path is not None:
            self.pose = pyrosetta.pose_from_pdb(pdb_path)
        elif rcsb is not None:
            self.pose = pyrosetta.pose_from_rcsb(rcsb)
        else:
            self.pose = Pose()
        # init pdb_info from pose
        self.pdb_info = self.pose.pdb_info() or self.add_pdbinfo_mover.apply(self.pose)
            
    @parameter_checker(pos = lambda x: x in {'before', 'after'})
    def add_resi(self, resi_name3: str, insert_pos: int, pos: str = 'after'):
        new_rsd = ResidueFactory.create_residue(self.rts.name_map(resi_name3))
        if pos == 'before':
            self.pose.prepend_polymer_residue_before_seqpos(new_rsd, insert_pos, True)
        else:
            self.pose.append_polymer_residue_after_seqpos(new_rsd, insert_pos, True)
        renumber_pdbinfo_based_on_conf_chains(self.pose)
        return self
            
    def delete_resi(self, resi_id: int):
        self.pose.delete_polymer_residue(resi_id)
        renumber_pdbinfo_based_on_conf_chains(self.pose)
        return self
    
    def mutate_resi(self, resi_id: int, new_name1: str, pack_radius: float = 0.0):
        mutate_residue(self.pose, resi_id, new_name1, pack_radius)
            
    def add_chain(self, chain: Union['_Pose', Pose], jump_resi_id: int = None) -> '_Pose':
        jump_resi_id = get_default_call_for_None(jump_resi_id, self.pose.total_residue)
        chain = chain.pose if isinstance(chain, _Pose) else chain
        self.pose.append_pose_by_jump(chain, jump_anchor_residue=jump_resi_id)
        self.pose.update_residue_neighbors()
        self.pose.update_pose_chains_from_pdb_chains()
        self.pose.conformation().detect_disulfides()
        renumber_pdbinfo_based_on_conf_chains(self.pose)
        return self
    
    def del_region(self, start_resi_id: int, stop_resi_id: int):
        delete_region(self.pose, start_resi_id, stop_resi_id)
        renumber_pdbinfo_based_on_conf_chains(self.pose)
    
    def merge_chain(self, chain: Union['_Pose', Pose], new_chain: bool = False) -> '_Pose':
        if isinstance(chain, _Pose):
            chain = chain.pose
        append_pose_to_pose(self.pose, chain, new_chain)
        return self
    
    def delete_chain(self, chain_id: int = None, chain_name: str = None) -> '_Pose':
        if isinstance(chain_name, str):
            chain_id = self.get_chain_id(chain_name)
        elif not isinstance(chain_id, int):
            raise ValueError('Either chain_id or chain_name must be specified.')
        self.delete_chain_mover.chain_num(chain_id)
        self.delete_chain_mover.apply(self.pose)
        return self
        
    def swap_chain(self, order: str | list[str]) -> '_Pose':
        if isinstance(order, list) and not isinstance(order, str):
            order = ''.join(order)
        else:
            raise ValueError('order must be a str or list of chain names')
        self.swap_chain_mover.chain_order(order)
        self.swap_chain_mover.apply(self.pose)
        return self
    
    def split_chain(self, return_Pose: bool = False) -> list['_Pose'] | list[Pose]:
        chains = self.pose.split_by_chain()
        return [_Pose(chain) for chain in chains] if return_Pose else chains
    
    def get_resi_id_in_pose_via_pdb(self, chain_name: str, resi_id: int) -> int:
        return self.pdb_info.pdb2pose(chain_name, resi_id)
    
    def get_resi_id_in_pdb_via_pose(self, resi_id: int) -> tuple[str, int]:
        resi_info = self.pdb_info.pose2pdb(resi_id)
        resi_id, chain_name, _ = resi_info.split(' ') # '2 A ' -> ['2', 'A', '']
        return chain_name, int(resi_id)
    
    def get_chain_id(self, chain_name: str) -> int:
        return get_chain_id_from_chain(chain_name, self.pose)
    
    def get_chain_name(self, chain_id: int = None, resi_id: int = None) -> str:
        if chain_id is not None:
            return get_chain_from_chain_id(chain_id, self.pose)
        elif resi_id is not None:
            return self.pdb_info.chain(resi_id)
        else:
            raise ValueError('Either chain_id or resi_id must be specified.')
        
    def get_chain_seq(self, chain_id: int = None, chain_name: str = None) -> str:
        if isinstance(chain_id, int):
            return self.pose.chain_sequence(chain_id)
        elif isinstance(chain_name, str):
            return self.pose.chain_sequence(self.get_chain_id(chain_name))
        else:
            raise ValueError('Either chain_id or chain_name must be specified.')
        

__all__ = [
    '_Pose',
    'renumber_pdbinfo_based_on_conf_chains',
]


if __name__ == '__main__':
    pyrosetta.init()
    pose = _Pose(seq='AAAA')
    pose.get_resi_id_in_pdb_via_pose(2)