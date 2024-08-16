'''
Date: 2024-08-15 19:54:22
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-16 15:25:49
Description: print selected residue names and numbers as autodock flex receptor residue format
'''
import tkinter as tk
import uuid
from typing import List

import Pmw
from pymol import cmd


def _get_res_info_from_sele(select: str, attrs: List[str] = None):
    if attrs is None:
        attrs = ['model', 'chain', 'resn', 'resi']
    atoms = {k:[] for k in attrs}
    chains = {}
    
    cmd.iterate(select, lambda atom: atoms.setdefault('model', []).append(atom.model))
    cmd.iterate(select, lambda atom: atoms.setdefault('chain', []).append(atom.chain))
    cmd.iterate(select, lambda atom: atoms.setdefault('resn', []).append(atom.resn))
    cmd.iterate(select, lambda atom: atoms.setdefault('resi', []).append(atom.resi))
    
    for m, c, r, i in zip(atoms['model'], atoms['chain'], atoms['resn'], atoms['resi']):
        chains.setdefault(f'{m}:{c}', {})[f'{r}{i}'] = True
        
    return chains
    

class LazyPocket:
    def __init__(self, app, _dev_mode: bool = False):
        if _dev_mode:
            parent = app
        else:
            parent = app.root

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
        chains = _get_res_info_from_sele(sele)
        final_output = ",".join(f"{k}:{'_'.join(v.keys())}" for k,v in chains.items())
        print('\nResidue names and numbers as autodock flex receptor residue format:')
        print(final_output)
        
        
# dev mode
if __name__ == '__main__':
    root = tk.Tk()
    Pmw.initialise(root)
    root.title('LazyDock Pymol Plugin - Lazy Pocket - Dev Mode')

    exitButton = tk.Button(root, text = 'Exit', command = root.destroy)
    exitButton.pack(side = 'bottom')
    widget = LazyPocket(root, _dev_mode=True)
    root.mainloop()
    