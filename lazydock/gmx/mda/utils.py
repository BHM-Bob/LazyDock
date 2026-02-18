'''
Date: 2025-03-14 16:12:21
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-09-13 21:38:22
Description: 
'''
from typing import List, Union
from MDAnalysis import Universe, AtomGroup


def filter_atoms_by_chains(atoms: AtomGroup, chains: Union[str, List[str]],
                           return_mask: bool = False) -> AtomGroup:
    '''
    过滤出指定链的AtomGroup

    Args:
        atoms (AtomGroup): 原始AtomGroup
        chains (Union[str, List[str]]): 指定链
        return_mask (bool, optional): 是否返回mask. Defaults to False.

    Returns:
        AtomGroup: 过滤后的AtomGroup
    '''
    if isinstance(chains, str):
        chains = [chains]
    chain_mask = atoms.chainIDs == chains[0]
    for chain_i in chains[1:]:
        chain_mask = chain_mask | (atoms.chainIDs == chain_i)
    if return_mask:
        return chain_mask
    return atoms[chain_mask]


def filter_atoms_by_resns(atoms: AtomGroup, resns: Union[str, List[str]],
                          return_mask: bool = False) -> AtomGroup:
    '''
    过滤出指定resn的AtomGroup

    Args:
        atoms (AtomGroup): 原始AtomGroup
        resns (Union[str, List[str]]): 指定resn
        return_mask (bool, optional): 是否返回mask. Defaults to False.

    Returns:
        AtomGroup: 过滤后的AtomGroup
    '''
    if isinstance(resns, str):
        resns = [resns]
    resn_mask = atoms.resnames == resns[0]
    for resn_i in resns[1:]:
        resn_mask = resn_mask | (atoms.resnames == resn_i)
    if return_mask:
        return resn_mask
    return atoms[resn_mask]
