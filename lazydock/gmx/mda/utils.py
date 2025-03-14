'''
Date: 2025-03-14 16:12:21
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-14 16:14:01
Description: 
'''
from typing import List, Union
from MDAnalysis import Universe, AtomGroup


def filter_atoms_by_chains(atoms: AtomGroup, chains: Union[str, List[str]]) -> AtomGroup:
    '''
    过滤出指定链的AtomGroup

    Args:
        atoms (AtomGroup): 原始AtomGroup
        chains (Union[str, List[str]]): 指定链

    Returns:
        AtomGroup: 过滤后的AtomGroup
    '''
    if isinstance(chains, str):
        chains = [chains]
    chain_mask = atoms.chainIDs == chains[0]
    for chain_i in chains[1:]:
        chain_mask = chain_mask | (atoms.chainIDs == chain_i)
    return atoms[chain_mask]