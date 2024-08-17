'''
Date: 2024-08-15 19:54:22
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-17 10:20:50
Description: print selected residue names and numbers as autodock flex receptor residue format
'''
import os
import tkinter as tk
import uuid
from typing import List

import Pmw
from mbapy.base import put_err
from pymol import api, cmd

if __name__ == '__main__':
    import lazydock_pymol_plugin
    from lazydock_pymol_plugin._autodock_utils import MyFileDialog
else:
    from ._autodock_utils import MyFileDialog


def _get_res_info_from_sele(select: str, attrs: List[str] = None):
    """
    Parameters
    ----------
        - select : str, pymol selection
        - attrs : List[str], optional, attributes to get, by default ['model', 'chain', 'resn', 'resi']

    Returns
    -------
        - chains : dict
            - key: model_str:chain_str
            - value: dict
                - key: res_str:res_id_str
                - value: True
        - chains_data : dict
            - key: model_str
            - value: dict
                - key: chain_str
                - value: List[List[res_str, res_id]]
    """
    if attrs is None:
        attrs = ['model', 'chain', 'resn', 'resi']
    atoms = {k:[] for k in attrs}
    chains = {}
    chains_data = {}
    
    cmd.iterate(select, lambda atom: atoms.setdefault('model', []).append(atom.model))
    cmd.iterate(select, lambda atom: atoms.setdefault('chain', []).append(atom.chain))
    cmd.iterate(select, lambda atom: atoms.setdefault('resn', []).append(atom.resn))
    cmd.iterate(select, lambda atom: atoms.setdefault('resi', []).append(atom.resi))
    
    for m, c, r, i in zip(atoms['model'], atoms['chain'], atoms['resn'], atoms['resi']):
        chains.setdefault(f'{m}:{c}', {})[f'{r}{i}'] = True
        chains_data.setdefault(m, {c:[]})[c].append([r, i])
        
    return chains, chains_data
    

class LazyPocket:
    def __init__(self, app, _dev_mode: bool = False):
        if _dev_mode:
            parent = app
        else:
            parent = app.root
            
        self.sele = None
        self.sele_chains = None

        self.dialog = Pmw.Dialog(parent,
                                 buttons=('Exit',),
                                 title = 'LazyDock Pymol Plugin - Lazy Pocket',
                                 command = self._gui_withdraw)
        self.dialog.withdraw() # ???
        self.dialog.geometry('650x780')
        self.dialog.bind('<Return>', self._gui_withdraw)
        
        # the title
        self.title_label = tk.Label(self.dialog.interior(),
                                    text='LazyDock Pymol Plugin - Lazy Pocket\nBHM-Bob G\n<https://github.com/BHM-Bob/LazyDock/>',
                                    background='blue', foreground='white')
        self.title_label.pack(expand=0, fill='both', padx=4, pady=4)
        
        # main notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=3, pady=3)
        
        # build pages
        self.lazy_pocket_page = self.notebook.add('Lazy Pocket')

        ## Lazy-Pocket PAGE
        ### print_sele_around_res GROUP
        self.print_sele_around_res_group = Pmw.Group(self.lazy_pocket_page,
                                                     tag_text='Print sele (around pocket) in flex residue format')
        self.print_sele_around_res_group.pack(fill='both', padx=10, pady=5)
        self.ui_print_sele_res_sele_label = Pmw.LabeledWidget(self.print_sele_around_res_group.interior(),
                                                              labelpos = 'w', label_text = 'when radius is blank, only print sele')
        self.ui_print_sele_res_sele_label.pack(fill='both', expand=1, padx=10, pady=5)
        self.ui_print_sele_res_sele_input = Pmw.EntryField(self.print_sele_around_res_group.interior(),
                                                           labelpos='w', label_text='sele=', value='sele')
        self.ui_print_sele_res_sele_input.pack(fill='both', expand=1, padx=10, pady=5)
        self.ui_print_sele_res_radius_input = Pmw.EntryField(self.print_sele_around_res_group.interior(),
                                                             labelpos='w', label_text='radius=', value='')
        self.ui_print_sele_res_radius_input.pack(fill='both', expand=1, padx=10, pady=5)
        self.ui_print_sele_res_button = Pmw.ButtonBox(self.print_sele_around_res_group.interior(),
                                                     labelpos='nw', label_text = 'Option:')
        self.ui_print_sele_res_button.add('print', command = self.print_sele_around_res)
        self.ui_print_sele_res_button.add('save rigid', command = self.save_rigid_receptor)
        self.ui_print_sele_res_button.pack(fill='both', expand=1, padx=10, pady=5)
        
        # GUI 
        self.dialog.show() # ???
        
    
    def _gui_withdraw(self, result):
        self.dialog.deactivate(result)
        self.dialog.withdraw()
        
    def print_sele_around_res(self):
        """output example: CB1:D:PHE108_PHE177_HIS178_LEU193_PHE379_SER383"""
        sele = self.ui_print_sele_res_sele_input.getvalue()
        radius = self.ui_print_sele_res_radius_input.getvalue()
        if radius:
            try:
                radius = float(radius)
            except Exception as e:
                print(f'radius should be a float number, but got {radius}.\nException: {e}\n')
                return
            pocket = f'Pocket_{str(uuid.uuid4())[1:5]}'
            cmd.select(pocket, f'byres {sele} around {radius}')
            sele = pocket
        self.sele = sele
        chains, self.sele_chains = _get_res_info_from_sele(sele)
        final_output = ",".join(f"{k}:{'_'.join(v.keys())}" for k,v in chains.items())
        print('\nResidue names and numbers as autodock flex receptor residue format:')
        print(final_output)
        
    def save_rigid_receptor(self):
        if self.sele_chains is None:
            return put_err(f'Please run "Print sele (around pocket) in flex residue format" first.')
        pdb_path = MyFileDialog(types = [('PDB file', '*.pdb')],
                                initialdir=os.getcwd()).getsavefile()
        if pdb_path is None:
            return put_err('Please choose a file to save.')
        if not pdb_path.endswith('.pdb'):
            pdb_path += '.pdb'
        is_first = True
        for model in self.sele_chains:
            tmp_model_name = f'{model}_{str(uuid.uuid4())[:4]}'
            api.copy(tmp_model_name, model)
            for chain in self.sele_chains[model]:
                for _, resi in self.sele_chains[model][chain]:
                    tmp_sele_name = f'{tmp_model_name}_{chain}_{resi}'
                    api.select(tmp_sele_name, f'model {tmp_model_name} and (chain {chain} and resi {resi})')
                    cmd.remove(tmp_sele_name)
                    cmd.delete(tmp_sele_name)
            api.multisave(pdb_path, tmp_model_name, append = 0 if is_first else 1)
            is_first = False
        print(f'Rigid receptor saved to {pdb_path}.')
            
        
# dev mode
if __name__ == '__main__':
    root = tk.Tk()
    Pmw.initialise(root)
    root.title('LazyDock Pymol Plugin - Lazy Pocket - Dev Mode')

    exitButton = tk.Button(root, text = 'Exit', command = root.destroy)
    exitButton.pack(side = 'bottom')
    widget = LazyPocket(root, _dev_mode=True)
    root.mainloop()
    