# Overview
This module provides a set of classes and functions designed to coloring molecular in PyMOL. It includes classes to represent atoms and residues with associated data, and a `Shader` class to apply color gradients based on calculated interaction values.

# Functions
## `get_c_value()`
### Description
Retrieves the `c_value` for a `ShaderRes` instance. If `c_value` is not set, it calculates the sum of `c_values` of its atoms.

### Parameters
- None

### Returns
- `float`: The `c_value` of the residue.

### Notes
- The `c_value` is only recalculated if `atoms` is updated after the initial calculation.

# Classes
## `ShaderAtom`
### Description
A data class representing an atom in a molecular model with properties such as model, chain, residue index, element, and an optional `c_value` for coloring.

### Members
- `model`: The model name.
- `chain`: The chain identifier.
- `resi`: The residue index.
- `elem`: The element type.
- `index`: The atom index.
- `c_value`: An optional float value associated with the atom.
- `alpha`: An optional float representing the opacity.

## `ShaderRes`
### Description
A data class representing a residue in a molecular model, containing a list of `ShaderAtom` instances and optional `c_value` and `alpha` properties.

### Methods
#### `get_c_value()`
##### Description
Retrieves the `c_value` for the residue. If not set, it calculates the sum of `c_values` of its atoms.

##### Parameters
- None

##### Returns
- `float`: The `c_value` of the residue.

### Members
- `atoms`: A list of `ShaderAtom` instances.
- `c_value`: An optional float value associated with the residue.
- `alpha`: An optional float representing the opacity.

## `ShaderValues`
### Description
A class to manage and retrieve `c_values` for chains, residues, and atoms in a molecular model.

### Methods
#### `get_c_value(chain: str, resi: int, atom_index: int = None) -> float`
##### Description
Retrieves the `c_value` for a specified chain, residue, and optionally an atom.

##### Parameters
- `chain`: The chain identifier.
- `resi`: The residue index.
- `atom_index`: The optional atom index.

##### Returns
- `float`: The `c_value` for the specified atom or residue.

#### `get_all_c_values(level: str = 'res') -> list`
##### Description
Retrieves all `c_values` for the specified level (residue or atom).

##### Parameters
- `level`: The level of detail ('res' for residue, 'atom' for atom).

##### Returns
- `list`: A list of tuples containing model, chain, residue index, atom index (if applicable), alpha, and `c_value`.

#### `from_interaction_df(df: Union[str, pd.DataFrame], model: str, sum_axis: int = 0) -> ShaderValues`
##### Description
Loads `c_values` from an interaction data frame.

##### Parameters
- `df`: The data frame or path to the data frame.
- `model`: The model name.
- `sum_axis`: The axis along which to sum the data frame.

##### Returns
- `ShaderValues`: The instance itself after loading the data.

## `Shader`
### Description
A class to apply color gradients to molecular models in PyMOL based on `c_values`.

### Methods
#### `get_col_from_c_value(c_value: float) -> (tuple, str)`
##### Description
Returns the RGBA color and a color name for a given `c_value`.

##### Parameters
- `c_value`: The `c_value` for which to generate the color.

##### Returns
- `tuple`: The RGBA color.
- `str`: The color name.

#### `_get_rgba_col_name(c_value: float, _col_name: str = None, _cmd = None) -> (tuple, str)`
##### Description
Generates the RGBA color and a color name for a given `c_value`, storing it globally if it doesn't exist.

##### Parameters
- `c_value`: The `c_value` for which to generate the color.
- `_col_name`: An optional predefined color name.
- `_cmd`: The PyMOL command module.

##### Returns
- `tuple`: The RGBA color.
- `str`: The color name.

#### `create_colors_in_pml(values: ShaderValues, level: str = 'res', names: List[str] = None, _cmd = None)`
##### Description
Creates colors in PyMOL for each `c_value` in the provided `ShaderValues`.

##### Parameters
- `values`: The `ShaderValues` instance.
- `level`: The level of detail ('res' for residue, 'atom' for atom).
- `names`: An optional list of predefined color names.
- `_cmd`: The PyMOL command module.

#### `auto_scale_norm(c_values)`
##### Description
Automatically scales the normalization based on the provided `c_values`.

##### Parameters
- `c_values`: A list of `c_values` to use for scaling.

#### `apply_shader_values(values: ShaderValues, level: str = 'res', auto_detect_vlim: bool = True, alpha_mode: str = None, _cmd = None)`
##### Description
Applies the shader values to a molecular model in PyMOL.

##### Parameters
- `values`: The `ShaderValues` instance.
- `level`: The level of detail ('res' for residue, 'atom' for atom).
- `auto_detect_vlim`: Whether to automatically detect the color limits.
- `alpha_mode`: An optional transparency mode for PyMOL.
- `_cmd`: The PyMOL command module.

#### `apply_shader_values_to_selection(selection, c_value: float = None, alpha: float = 1.0, alpha_mode: str = None, _cmd = None)`
##### Description
Applies the shader values to a specific selection in PyMOL.

##### Parameters
- `selection`: The selection to which the shader values will be applied.
- `c_value`: An optional `c_value` for the color.
- `alpha`: The opacity level.
- `alpha_mode`: An optional transparency mode for PyMOL.
- `_cmd`: The PyMOL command module.

#### `apply_shader_values_to_sele(select_expression: str, c_value: float = None, alpha: float = 1.0, alpha_mode: str = None, _cmd = None)`
##### Description
Applies the shader values to a selection expression in PyMOL.

##### Parameters
- `select_expression`: The selection expression to which the shader values will be applied.
- `c_value`: An optional `c_value` for the color.
- `alpha`: The opacity level.
- `alpha_mode`: An optional transparency mode for PyMOL.
- `_cmd`: The PyMOL command module.

#### `__repr__(self)`
##### Description
Returns a string representation of the `Shader` instance.

##### Returns
- `str`: The string representation.
