'''
Date: 2024-08-15 19:54:22
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-15 20:05:09
Description: print selected residue names and numbers as autodock flex receptor residue format
'''
from pymol import cmd

def __init_plugin__(app):
    app.menuBar.addmenuitem('Plugin', 'command',
        label='print sele res',
        command=lambda: print_sele_res(app.root))

def print_sele_res(parent):
    """output example: CB1:D:PHE108_PHE177_HIS178_LEU193_PHE379_SER383""" 
    attrs = ['model', 'chain', 'resn', 'resi']
    atoms = {k:[] for k in attrs}
    chains = {}
    
    cmd.iterate("sele", lambda atom: atoms.setdefault('model', []).append(atom.model))
    cmd.iterate("sele", lambda atom: atoms.setdefault('chain', []).append(atom.chain))
    cmd.iterate("sele", lambda atom: atoms.setdefault('resn', []).append(atom.resn))
    cmd.iterate("sele", lambda atom: atoms.setdefault('resi', []).append(atom.resi))
    
    for m, c, r, i in zip(atoms['model'], atoms['chain'], atoms['resn'], atoms['resi']):
        chains.setdefault(f'{m}:{c}', {})[f'{r}{i}'] = True

    final_output = ",".join(f"{k}:{'_'.join(v.keys())}" for k,v in chains.items())
    print(final_output)