'''
Date: 2024-08-31 21:40:56
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-09-03 16:52:25
Description: 
'''
from dataclasses import dataclass
from threading import Lock
from typing import Dict, List, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Colormap, TwoSlopeNorm
from mbapy_lite.plot import rgb2hex
from pymol import cmd

from lazydock.utils import uuid4


@dataclass
class ShaderAtom:
    model: str
    chain: str
    resi: int
    elem: str
    index: int
    c_value: float = None
    alpha: float = 1.0


@dataclass
class ShaderRes:
    model: str
    chain: str
    resi: int
    atoms: List[ShaderAtom] = None
    c_value: float = None
    alpha: float = 1.0

    def get_c_value(self):
        """
        Returns:
            - c_value, if it is not None.
            - sum of c_values of atoms, if atoms is not None and c_value is None.
        Notes: 
            - If atoms is updated after c_value calculation, c_value will not be updated.
        """
        if self.atoms is not None and self.c_value is None:
            self.c_value = sum(map(lambda x: x.c_value, filter(lambda x: x.c_value, self.atoms)))
        return self.c_value


class ShaderValues:
    def __init__(self, chains: Dict[str, List[ShaderRes]] = None) -> None:
        self.chains: Dict[str, List[ShaderRes]] = chains or {}

    def get_c_value(self, chain: str, resi: int, atom_index: int = None):
        if chain not in self.chains:
            raise ValueError(f'Chain {chain} not found in shader values.')
        res = [res for res in self.chains[chain] if res.resi == resi]
        if not res:
            raise ValueError(f'Residue {resi} not found in chain {chain}.')
        if atom_index is not None:
            atom = [atom for atom in res[0].atoms if atom.index == atom_index]
            if not atom:
                raise ValueError(f'Atom {atom_index} not found in residue {resi}.')
            return atom[0].c_value
        return res[0].get_c_value()
    
    def get_all_c_values(self, level: str = 'res'):
        if level not in {'res', 'atom'}:
            raise ValueError(f'Level must be "res" or "atom", got {level}.')
        if level == 'res':
            return [(res.model, res.chain, res.resi, res.alpha, res.get_c_value()) for chain in self.chains for res in self.chains[chain]]
        else:
            return [(atom.model, atom.chain, atom.resi, atom.index, res.alpha, atom.c_value) for chain in self.chains for res in self.chains[chain] for atom in res.atoms]
          
    def from_interaction_df(self, df: Union[str, pd.DataFrame], model: str, sum_axis: int = 0):
        """
        load values from interaction_df calculated by lazudock.interaction_utils.calcu_receptor_poses_interaction
        
        Parameters:
            df (Union[str, pd.DataFrame]): interaction_df or path to interaction_df.
            model (str): model name.
            sum_axis (int): 0 or 1, pass to df.sum, 0 means sum each rows and get a colum, 1 means sum each colums and get a row.
            
        Returns:
            ShaderValues: self.
        """
        if isinstance(df, str):
            df = pd.read_excel(df, index_col = 0)
        df = df.sum(axis=sum_axis)[1:] # first row is index and colum names
        for res, v in zip(df.index, df.values):
            chain, resi, _ = res.split(':')
            self.chains.setdefault(chain, []).append(ShaderRes(model, chain, int(resi), c_value=v))
        return self


class Shader:
    global_locker = Lock()
    global_name2col = {}
    def __init__(self, cmap: Union[str, Colormap] = 'coolwarm',
                 norm: TwoSlopeNorm = None,
                 col_name_prefix: str = 'COL_') -> None:
        self.cmap = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
        self.norm = norm or TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
        self.COL_NAME_PREFIX = col_name_prefix
        
    def get_col_from_c_value(self, c_value: float):
        """return rgba and color name with prefix for given c_value."""
        rgba = self.cmap(self.norm(c_value))
        name = rgb2hex(*[int(rgba[i]*255) for i in range(3)])
        return rgba, f'{self.COL_NAME_PREFIX}{name[1:]}'
                
    def _get_rgba_col_name(self, c_value: float, _col_name: str = None, _cmd = None):
        """get rgba and color name with prefix for given c_value, 
        store name in self.global_name2col if not exist."""
        _cmd = _cmd or cmd
        rgba, col_name = self.get_col_from_c_value(c_value)
        col_name = _col_name or col_name
        if col_name not in self.global_name2col:
            _cmd.set_color(col_name, list(rgba[:-1]))
            with self.global_locker: # NOTE: only lock-unlock when adding new color to avoid ReLock
                self.global_name2col[col_name] = rgba
        return rgba, col_name
        
    def create_colors_in_pml(self, values: ShaderValues, level: str = 'res',
                             names: List[str] = None, _cmd = None):
        """
        set color and it's name in pymol for each c_value in values.
        
        Parameters:
            values (ShaderValues): values to create colors for.
            level (str): 'atom' or'res', level to create colors for.
            names (List[str]): list of color names to use. If None, will generate names from c_values.
            cmd: pymol cmd module, can be pymol.cmd or lazudock.pml.server.PymolAPI
        Notes:
            If name is given, the color name prefix will not be added.
        """
        c_values = np.unique([v[-1] for v in values.get_all_c_values(level)])
        if isinstance(names, list) and len(names) != len(c_values):
            raise ValueError(f'If names is given, names must be equal to length of c_values, got {len(names)} and {len(c_values)}.')
        for i, c in enumerate(c_values):
            name = names[i] if names is not None else None
            if name not in self.global_name2col:
                self._get_rgba_col_name(c, name, _cmd)
                
    def auto_scale_norm(self, c_values):
        self.norm.autoscale(list(map(lambda x: x[-1], c_values)))
        
                
    def apply_shader_values(self, values: ShaderValues, level: str = 'res',
                            auto_detect_vlim: bool = True, alpha_mode: str = None,
                            _cmd = None):
        """
        apply shader values to pymol.
        
        Parameters:
            values (ShaderValues): values to apply.
            level (str): 'atom' or'res', level to apply values for.
            auto_detect_vlim (bool): if True, will set vmin and vmax for colormap based on min and max c_values.
            alpha_mode (str): pymol transparency mode, such as `cartoon_transparency`.
            _cmd: pymol cmd module, can be pymol.cmd or lazudock.pml.server.PymolAPI
        """
        _cmd = _cmd or cmd
        if level not in {'res', 'atom'}:
            raise ValueError(f'Level must be "res" or "atom", got {level}.')
        c_values = values.get_all_c_values(level)
        if auto_detect_vlim:
            self.auto_scale_norm(c_values)
        # loop through residues or atoms and apply color
        for c_value in c_values:
            if level =='res':
                model, chain, resi, alpha, c = c_value
                _, col_name = self._get_rgba_col_name(c, _cmd=_cmd)
                sele_exp = f'model {model} and (chain {chain} and resi {resi})'
                _cmd.color(col_name, sele_exp)
                if alpha_mode is not None:
                    _cmd.set(alpha_mode, alpha, sele_exp)
            else: # level == 'atom' and res.atoms is not None
                model, chain, resi, atom_index, alpha, c = c_value
                _, col_name = self._get_rgba_col_name(c, _cmd=_cmd)
                sele_exp = f'(model {model} and chain {chain}) and (resi {resi} and index {atom_index})'
                _cmd.color(col_name, sele_exp)
                if alpha_mode is not None:
                    _cmd.set(alpha_mode, alpha, sele_exp)
                    
    def apply_shader_values_to_selection(self, selection, c_value: float = None,
                                         alpha: float = 1.0, alpha_mode: str = None,
                                         _cmd = None):
        _cmd = _cmd or cmd
        if c_value is not None:
            _, col_name = self._get_rgba_col_name(c_value, _cmd=_cmd) # NOTE: do not use alpha from cmap
            _cmd.color(col_name, selection)
        if alpha_mode is not None:
            _cmd.set(alpha_mode, alpha, selection)
                    
    def apply_shader_values_to_sele(self, select_expression: str, c_value: float = None,
                                    alpha: float = 1.0, alpha_mode: str = None,
                                    _cmd = None):
        _cmd = _cmd or cmd
        selection = uuid4()
        _cmd.select(selection, select_expression)
        self.apply_shader_values_to_selection(selection, c_value, alpha, alpha_mode, _cmd)
        _cmd.delete(selection)
                    
    def __repr__(self):
        return f'{self.cmap.name}({self.norm.vmin:.2f},{self.norm.vcenter:.2f},{self.norm.vmax:.2f}){self.COL_NAME_PREFIX}'
    
    
__all__ = [
    'ShaderAtom',
    'ShaderRes',
    'ShaderValues',
    'Shader',
]
    
    
if __name__ == '__main__':
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    values = ShaderValues().from_interaction_df('data_tmp/interaction_df.xlsx', 'receptor')
    shader = Shader()
    shader.create_colors_in_pml(values)
    shader.apply_shader_values(values, alpha_mode='cartoon_transparency')
    print(shader)
    cmd.clip('near', -10)
    cmd.clip('slab', 100)
    cmd.zoom('model receptor')
    cmd.png('data_tmp/shader.png', 1200, 1200, 300)