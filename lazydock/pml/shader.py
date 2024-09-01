'''
Date: 2024-08-31 21:40:56
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-09-01 10:20:05
Description: 
'''
from dataclasses import dataclass, field
from typing import Dict, List, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Colormap
from mbapy.plot import rgb2hex, rgbs2hexs
from pymol import cmd


@dataclass
class ShaderAtom:
    model: str
    chain: str
    resi: int
    elem: str
    index: int
    c_value: float = None
    
@dataclass
class ShaderRes:
    model: str
    chain: str
    resi: int
    atoms: List[ShaderAtom] = None
    c_value: float = None
    
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
            return [(res.model, res.chain, res.resi, res.get_c_value()) for chain in self.chains for res in self.chains[chain]]
        else:
            return [(atom.model, atom.chain, atom.resi, atom.index, atom.c_value) for chain in self.chains for res in self.chains[chain] for atom in res.atoms]
          
    def from_interaction_df(self, df: pd.DataFrame, sum_axis: int = 1):
        pass


class Shader:
    def __init__(self, cmap: Union[str, Colormap] = 'coolwarm', col_name_prefix: str = 'COL_') -> None:
        self.cmap = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
        self.COL_NAME_PREFIX = col_name_prefix
        self.name2col = {}
        
    def get_col_from_c_value(self, c_value: float):
        """return rgba and color name with prefix for given c_value."""
        rgba = self.cmap(c_value)
        name = rgb2hex(*[int(rgba[i]*255) for i in range(3)])
        return rgba, f'{self.COL_NAME_PREFIX}{name[1:]}'
                
    def _get_rgba_col_name(self, c_value: float, _col_name: str = None):
        """get rgba and color name with prefix for given c_value, 
        store name in self.name2col if not exist."""
        rgba, col_name = self.get_col_from_c_value(c_value)
        col_name = _col_name or col_name
        if col_name not in self.name2col:
            cmd.set_color(col_name, list(rgba[:-1]))
            self.name2col[col_name] = rgba
        return rgba, col_name
        
    def create_colors_in_pml(self, values: ShaderValues, level: str = 'res', names: List[str] = None):
        """
        set color and it's name in pymol for each c_value in values.
        
        Parameters:
            values (ShaderValues): values to create colors for.
            level (str): 'atom' or'res', level to create colors for.
            names (List[str]): list of color names to use. If None, will generate names from c_values.
            
        Notes:
            If name is given, the color name prefix will not be added.
        """
        c_values = np.unique([v[-1] for v in values.get_all_c_values(level)])
        if isinstance(names, list) and len(names) != len(c_values):
            raise ValueError(f'If names is given, names must be equal to length of c_values, got {len(names)} and {len(c_values)}.')
        for i, c in enumerate(c_values):
            name = names[i] if names is not None else None
            if name not in self.name2col:
                self._get_rgba_col_name(c, name)
                
    def apply_shader_values(self, values: ShaderValues, level: str = 'res'):
        if level not in {'res', 'atom'}:
            raise ValueError(f'Level must be "res" or "atom", got {level}.')
        # loop through residues or atoms and apply color
        for res in [res for chain in values.chains.values() for res in chain]:
            if level =='res' or (level == 'atom' and res.atoms is None):
                _, col_name = self._get_rgba_col_name(res.get_c_value())
                cmd.color(col_name, f'model {res.model} and (chain {res.chain} and resi {res.resi})')
            else: # level == 'atom' and res.atoms is not None
                for atom in res.atoms:
                    _, col_name = self._get_rgba_col_name(atom.c_value)
                    cmd.color(col_name, f'(model {atom.model} and chain {atom.chain}) and (resi {atom.resi} and index {atom.elem})')
    
    
__all__ = [
    'ShaderAtom',
    'ShaderRes',
    'ShaderValues',
    'Shader',
]
    
    
if __name__ == '__main__':
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    res1 = ShaderRes('receptor', 'A', 46, None, 0.5)
    value1 = ShaderValues({'A': [res1]})
    shader = Shader()
    shader.create_colors_in_pml(value1)
    shader.apply_shader_values(value1)
    shader.apply_shader_values(value1, 'atom')