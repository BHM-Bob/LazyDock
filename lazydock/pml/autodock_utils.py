# copied and fixed from Autodock Plugin: http://www.pymolwiki.org/index.php/autodock_plugin


from __future__ import absolute_import, print_function

import re
import tkinter as tk
import tkinter.filedialog as tkFileDialog

# atom-type, atom-number, atom-name, residue-name, chain-name, residue-number, x, y, z, occupancy, temperature-factor
# ATOM      1  CA  LYS     7     136.747 133.408 135.880 -0.06 +0.10
PDB_PATTERN = r"(ATOM|HETATM) +(\d+) +(\w+) +(\w+) +(\w+)? +(\d+) +([\d\-\.]+) +([\d\-\.]+) +([\d\-\.]+) +([\+\-][\d\-\.]+) +([\+\-][\d\-\.]+)"


class ADModel:
    """ STORAGE CLASS FOR DOCKED LIGANDS """
    def __init__(self, lst=None, sort_atom_by_res: bool = False):
        self.atomlines = []
        self.energy = 0.
        self.name = ''
        self.poseN = 0
        self.info = ''
        self.lst = []
        self.num = 0
        self.as_string = ''
        self.pdb_lines = []
        if lst is not None:
            self.from_list(lst, sort_atom_by_res)
            
    def sort_atom_by_res(self, pdb_string: str = None):
        pdb_string = pdb_string or self.as_string
        pdb_string_lines = pdb_string.strip().split('\n')
        for idx, line in enumerate(pdb_string_lines):
            matches = re.findall(PDB_PATTERN, line)
            items = list(matches[0])
            items[1], items[5] = int(items[1]), int(items[5])
            self.pdb_lines.append(items+[idx]) # line items + line idx
        # sort by chain name first, then by residue number, then by atom number
        self.pdb_lines.sort(key=lambda x: (x[4], x[5], x[1]))
        # apply sorted line idx to as_string
        self.as_string = '\n'.join([pdb_string_lines[i[-1]] for i in self.pdb_lines])
        return self.as_string
            
    def from_list(self, lst, sort_atom_by_res: bool = False):
        self.info = '\n'.join(lst)
        self.lst = lst
        for line in lst:
            if line.startswith('ATOM') or \
                    line.startswith('HETATM'):
                self.atomlines.append(line)
                self.as_string += 'ATOM  ' + line[6:67] + '\n'
            elif line.startswith('USER'):
                if 'Free Energy of Binding' in line:
                    entr = line.split('=')[1]
                    self.energy = float(entr.split()[0])
            elif line.startswith('REMARK'):
                if 'VINA RESULT' in line:
                    entr = line.split(':')[1]
                    self.energy = float(entr.split()[0])
        if sort_atom_by_res:
            self.sort_atom_by_res()

    def as_pdb_string(self):
        return self.as_string
        
    def info_string(self):
        return self.info


def tk_file_dialog_wrapper(*args, **kwargs):
    def ret_wrapper(tk_file_dialog_func):
        def core_wrapper(*args, **kwargs):
            parent = tk.Tk()
            result = tk_file_dialog_func(*args, parent = parent, **kwargs)
            parent.withdraw()
            if result == "":
                return None
            else:
                return result
        return core_wrapper
    return ret_wrapper
            

class MyFileDialog:
    def __init__(self, types=[("Executable", "*")], initialdir: str = None):
        self.initialdir = initialdir
        self.types = types

    @tk_file_dialog_wrapper()
    def get_open_file(self, parent):
        return tkFileDialog.askopenfilename(parent=parent, initialdir = self.initialdir, filetypes=self.types)


    @tk_file_dialog_wrapper()
    def get_save_file(self, parent):
        return tkFileDialog.asksaveasfilename(parent=parent, initialdir = self.initialdir, filetypes=self.types)
        
    @tk_file_dialog_wrapper()
    def get_ask_dir(self, parent):
        return tkFileDialog.askdirectory(parent=parent, initialdir = self.initialdir)
