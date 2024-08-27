# Module Description

This Python module provides a suite of functions for calculating and managing molecular interactions, specifically hydrogen bonds, within molecular structures, likely used in the context of computational chemistry or structural biology. It includes functions for setting and getting hydrogen bond cutoff distances, retrieving distance information from PyMOL distance objects, calculating atom-level interactions, sorting these interactions based on molecular models, merging interaction data into a DataFrame, and filtering this DataFrame based on distance cutoffs.

# Functions

## get_h_bond_cutoff()-> (float, float)
### Function Introduction
Retrieves the current distance cutoff values for identifying hydrogen bonds, which include both the optimal (center) and the minimal acceptable (edge) distances.

### Parameters
None

### Return Value
A tuple containing two floats representing the center and edge distance cutoffs for hydrogen bonds.

### Notes
- Used to understand the current parameters for hydrogen bond identification.

### Example
```python
center, edge = get_h_bond_cutoff()
```

## set_h_bond_cutoff(center: float, edge: float)
### Function Introduction
Sets the distance cutoff values used to identify hydrogen bonds.

### Parameters
- `center`: The ideal distance for a hydrogen bond.
- `edge`: The maximum acceptable distance for a hydrogen bond.

### Return Value
None

### Notes
- Allows customization of hydrogen bond identification criteria.

### Example
```python
set_h_bond_cutoff(center=3.5, edge=4.0)
```

## get_distance_info(dist_name: str, state=1, selection='all', xyz2atom=None)-> list of tuples
### Function Introduction
Extracts detailed information about distances between atoms from a PyMOL distance object.

### Parameters
- `dist_name`: The name of the distance object.
- `state`: The state of the object to consider.
- `selection`: The selection of atoms to focus on.
- `xyz2atom`: A dictionary mapping coordinates to atom information.

### Return Value
A list of tuples, each containing two sets of atom information and a distance value.

### Notes
- The order of atoms in the tuples is not guaranteed to match the input selection order.

### Example
```python
distance_info = get_distance_info('distance_object_name')
```

## calcu_atom_level_interactions(sele1: str, sele2: str, mode=0, cutoff=3.6, xyz2atom: dict = None)
### Function Introduction
Calculates distances between two sets of atoms using PyMOL's `distance` command and retrieves the atom information and distances.

### Parameters
- `sele1`: The first atom selection.
- `sele2`: The second atom selection.
- `mode`: The mode for the `distance` command.
- `cutoff`: The maximum distance to consider.
- `xyz2atom`: A dictionary for mapping coordinates to atom information.

### Return Value
A dictionary containing atom information, mappings, residue information, and the calculated interactions.

### Notes
- Utilizes temporary selections and distance objects which are deleted after calculation.

### Example
```python
atoms, xyz2atom, residues, interactions = calcu_atom_level_interactions('sele1', 'sele2')
```

## sort_atom_level_interactions(interactions: dict, model1: str, model2: str)
### Function Introduction
Sorts interactions such that the specified models appear in a consistent order within the interaction tuples.

### Parameters
- `interactions`: The interactions dictionary to sort.
- `model1`: The first model name to prioritize in sorting.
- `model2`: The second model name to prioritize in sorting.

### Return Value
The sorted interactions dictionary.

### Notes
- Ensures that the order of models in the interaction tuples is predictable.

### Example
```python
interactions = sort_atom_level_interactions(interactions, 'modelA', 'modelB')
```

## merge_interaction_df(interaction: dict, interaction_df: pd.DataFrame, distance_cutoff: float, nagetive_factor: float)
### Function Introduction
Merges atom-level interactions into a DataFrame, which represents the interactions as a matrix of distances between residues.

### Parameters
- `interaction`: The dictionary of interactions to merge.
- `interaction_df`: The DataFrame to update.
- `distance_cutoff`: The cutoff distance used to scale the interaction values.
- `nagetive_factor`: A factor to adjust the interaction values for certain atom types.

### Return Value
The updated DataFrame with merged interaction data.

### Notes
- The function modifies the DataFrame in place to include interaction values.

### Example
```python
interaction_df = merge_interaction_df(interactions, interaction_df, 4.0, -1.0)
```

## calcu_receptor_poses_interaction(receptor: str, poses: list, mode: int = 0, cutoff: float = 4., nagetive_factor: float = -1.)
### Function Introduction
Calculates interactions between a receptor and multiple poses of a ligand, and compiles the results into a DataFrame.

### Parameters
- `receptor`: The receptor's PyMOL object name.
- `poses`: A list of ligand PyMOL object names.
- `mode`: The calculation mode for the distances.
- `cutoff`: The maximum distance for considering an interaction.
- `nagetive_factor`: A factor to adjust the interaction values.

### Return Value
A tuple containing a dictionary of interactions and the compiled DataFrame.

### Notes
- The function performs calculations for each pose and merges the results into a single DataFrame.

### Example
```python
interactions, interaction_df = calcu_receptor_poses_interaction('receptor', ['pose1', 'pose2'])
```

## filter_interaction_df(interaction_df: pd.DataFrame, colum_axis_min: float = None, row_axis_min: float = None, inplace: bool = False)
### Function Introduction
Filters a DataFrame of interactions based on minimum distance values for rows and columns.

### Parameters
- `interaction_df`: The DataFrame to filter.
- `colum_axis_min`: The minimum value for column filtering.
- `row_axis_min`: The minimum value for row filtering.
- `inplace`: Whether to filter in place or return a new DataFrame.

### Return Value
The filtered DataFrame.

### Notes
- Rows and/or columns with maximum values below the specified cutoffs are removed.

### Example
```python
interaction_df = filter_interaction_df(interaction_df, colum_axis_min=3.5, row_axis_min=3.5, inplace=True)
```

# Classes

There are no classes defined in the provided code snippet. The code consists only of functions and their associated logic.
