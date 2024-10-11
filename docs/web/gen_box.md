# Brief Description
This module provides functions to interact with the ProteinPlus service for downloading and parsing pocket information from protein-ligand complexes. It includes functionality to automate the process of downloading zip files containing pocket data, parsing these files to extract box information, and visualizing the pockets in PyMOL.

# Functions

## `get_pocket_box_from_ProteinPlus(receptor_path: str, ligand_path: str, browser: Browser = None)`-> `None`
### Function Overview
Downloads the pocket result zip file from ProteinPlus using the provided receptor and ligand file paths.

### Parameters
- `receptor_path`: `str`, the path to the receptor PDB file.
- `ligand_path`: `str`, the path to the ligand SDF file.
- `browser`: `lazydock.web.Browser`, an optional browser object to use. If not provided, a new browser will be created.

### Returns
`None`, the pocket result zip file will be downloaded to the download path specified in the browser.

### Notes
- The function assumes that the provided paths are valid and the files exist.
- If a browser object is not provided, one will be created and used for the download process.

### Example
```python
get_pocket_box_from_ProteinPlus('path/to/receptor.pdb', 'path/to/ligand.sdf')
```

## `parse_pocket_box_from_ProteinPlus(result_path: str, k: Union[int, List[int]] = None, reinitialize: bool = False, draw_box: bool = True, _cmd = None, method: str = 'extend')`-> `box_df: pandas.DataFrame, box_info: dict`
### Function Overview
Parses the zip result file from ProteinPlus and returns box information.

### Parameters
- `result_path`: `str`, the path to the zip file downloaded from ProteinPlus.
- `k`: `int` or `List[int]`, the index of the pocket to be parsed. If `None`, all pockets will be parsed.
- `reinitialize`: `bool`, whether to reinitialize PyMOL or not.
- `draw_box`: `bool`, whether to draw the box or not.
- `_cmd`: `pymol.cmd`, an optional PyMOL command object.
- `method`: `str`, the method to calculate box center and size. 'extend' uses PyMOL's get_extent function, 'mean' calculates mean coordinates of all atoms in the pocket.

### Returns
- `box_df`: `pandas.DataFrame`, a dataframe of pocket information. Only returned in the 'mean' method.
  - 'resn_resi_index': `str`, residue name, residue index, and atom index.
  - 'X', 'Y', 'Z': `float`, coordinates of the atom.
- `box_info`: `dict`, box information.
  - 'box_center': `list` of `float`, center of the box.
  - 'box_size': `list` of `float`, size of the box.
  - 'box_vertices': `list` of `list` of `float`, vertices of the box.

### Notes
- The function requires the result zip file to be in the specified path.
- If `reinitialize` is `True`, PyMOL will be reinitialized before parsing.
- The `method` parameter determines the approach to calculating the box center and size.

### Example
```python
df, box = parse_pocket_box_from_ProteinPlus('path/to/POCKET.zip', [1], True, method='mean')
print(df)
print(box)
```

# Classes

This module does not contain any classes.
