# copied and fixed from Autodock Plugin: http://www.pymolwiki.org/index.php/autodock_plugin


from __future__ import absolute_import, print_function

import re
import tkinter as tk
import tkinter.filedialog as tkFileDialog
from typing import Any, Callable, List
from mbapy.file import opts_file, decode_bits_to_str
from mbapy.game import BaseInfo

# atom-type, atom-number, atom-name, residue-name, chain-name, residue-number, x, y, z, occupancy, temperature-factor
# ATOM      1  CA  LYS     7     136.747 133.408 135.880 -0.06 +0.10
PDB_PATTERN = r"(ATOM|HETATM) +(\d+) +(\w+) +(\w+) +(\w+)? +(\d+) +([\d\-\.]+) +([\d\-\.]+) +([\d\-\.]+) +([\+\-][\d\-\.]+) +([\+\-][\d\-\.]+)"


class ADModel(BaseInfo):
    """STORAGE CLASS FOR DOCKED LIGAND"""
    def __init__(self, content: str = None, _sort_atom_by_res: bool = False,
                 _parse2std: bool = False):
        self.energy = 0.
        self.name = ''
        self.poseN = 0
        self.info = ''
        self.pdb_string = ''
        self.pdb_lines = []
        if content is not None:
            self.parse_content(content, _sort_atom_by_res, _parse2std)

    def parse_content(self, content: str, _sort_atom_by_res: bool = False,
                      _parse2std: bool = False):
        # parse pdb lines
        self.info = content
        self.pdb_lines = list(map(lambda x: x[0], re.findall(r'((ATOM|HETATM).+?\n)', self.info)))
        self.pdb_atoms = re.findall(PDB_PATTERN, self.info)
        if _sort_atom_by_res or _parse2std:
            pack = sorted(zip(self.pdb_lines, self.pdb_atoms), key = lambda x : (x[1][4], int(x[1][5]), int(x[1][1])))
            self.pdb_lines, self.pdb_atoms = zip(*pack)
        if _parse2std:
            self.pdb_string = '\n'.join([' '.join(line) for line in self.pdb_atoms])
        else:
            self.pdb_string = ''.join(self.pdb_lines)
        # parse energy could be str.find?
        energy_line = re.findall(r'USER.+?Free Energy of Binding.+?\n', self.info)
        if energy_line:
            entr = energy_line[0].split('=')[1]
            self.energy = float(entr.split()[0])
        else:
            energy_line = re.findall(r'REMARK.+?VINA RESULT.+?\n', self.info)
            if energy_line:
                entr = energy_line[0].split(':')[1]
                self.energy = float(entr.split()[0])
            else:
                self.energy = None

    def as_pdb_string(self):
        return self.pdb_string

    def info_string(self):
        return self.info


class DlgFile(BaseInfo):
    def __init__(self, path: str = None, content: str = None,
                 sort_pdb_line_by_res: bool = False,
                 parse2std: bool = False):
        # load content from file if path is provided
        self.path = path
        self.sort_pdb_line_by_res = sort_pdb_line_by_res
        self.parse2std = parse2std
        if path is not None:
            self.content = decode_bits_to_str(opts_file(path, 'rb'))
        elif content is not None:
            self.content = content
        # decode content to pose_lst
        if self.content is not None:
            self.pose_lst: List[ADModel] = self.decode_content()
            self.sort_pose()
        else:
            self.pose_lst: List[ADModel] = []
        self.n2i = {}
        
    def __len__(self):
        return len(self.pose_lst)
        
    def sort_pose(self, key: Callable[[ADModel], Any] = None,
                  inplace: bool = True, reverse: bool = False) -> None:
        if key is None:
            key = lambda x : x.energy
        if inplace:
            self.pose_lst.sort(key=key, reverse=reverse)
            ret = self.pose_lst
        else:
            ret = sorted(self.pose_lst, key=key, reverse=reverse)
        return ret
        
    def decode_content(self):
        dlg_pose_lst = []
        for model in re.findall('MODEL.+?ENDMDL', self.content, re.DOTALL):
            model = model.replace('\nDOCKED: ', '\n')
            dlg_pose_lst.append(ADModel(model, self.sort_pdb_line_by_res, self.parse2std))
        return dlg_pose_lst
    
    def asign_pose_name(self, pose_names: List[str]):
        if len(pose_names) != len(self.pose_lst):
            raise ValueError("Number of pose names must match number of poses")
        for i, name in enumerate(pose_names):
            self.n2i[name] = i
        
    def asign_prop(self, prop: str, value: List[Any]):
        if value is not None and len(value) == len(self.pose_lst):
            setattr(self, prop, value)
                
    def set_pose_prop(self, prop: str, value: Any, pose_name: str = None, pose_idx: int = None):
        if pose_name is None and pose_idx is None:
            raise ValueError("Either pose_name or pose_idx must be provided")
        if pose_name is not None:
            pose_idx = self.n2i[pose_name]
        if not hasattr(self, prop):
            setattr(self, prop, [None]*len(self.pose_lst))
        getattr(self, prop)[pose_idx] = value
            
    def get_pose(self, pose_name: str = None, pose_idx: int = None):
        if pose_name is None and pose_idx is None:
            raise ValueError("Either pose_name or pose_idx must be provided")
        if pose_name is not None:
            pose_idx = self.n2i[pose_name]
        return self.pose_lst[pose_idx]
    
    def get_pose_prop(self, prop: str, pose_name: str = None, pose_idx: int = None, default: Any = None):
        if pose_name is None and pose_idx is None:
            raise ValueError("Either pose_name or pose_idx must be provided")
        if pose_name is not None:
            pose_idx = self.n2i[pose_name]
        if hasattr(self, prop):
            return getattr(self, prop)[pose_idx]
        else:
            return default



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


if __name__ == '__main__':
    # from mbapy.base import TimeCosts
    # @TimeCosts(10)
    # def load_test(idx):
    #     DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=True)
    # load_test()
    normal = DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=True, parse2std=False)
    std = DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=True, parse2std=True)
    pass