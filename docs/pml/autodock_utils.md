<!--
 * @Date: 2024-08-27 10:15:51
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-10-09 20:58:22
 * @Description: 
-->

# Module Description
This module provides a collection of classes and functions for handling and parsing AutoDock docing result file(dlg) and PDB (Protein Data Bank) formatted strings, which are commonly used to represent the three-dimensional structures of molecules.

# Description

This module is designed to facilitate the process of setting up docking runs with Autodock and viewing the results using PyMOL. It provides functionalities to parse and manipulate PDB formatted strings, manage docking poses, and interact with the file system through a graphical user interface.

# Functions

## tk_file_dialog_wrapper(*args, **kwargs)
### Function Description
This function is a wrapper that creates a dialog box for file operations. It initializes a Tkinter window, invokes the specified file dialog function, and ensures that the dialog is closed properly.

### Parameters
- `tk_file_dialog_func`: The file dialog function to be wrapped.

### Returns
A wrapped function that handles the file dialog operation.

### Notes
- This function is used to create a graphical file dialog that interacts with the user.

### Example
```python
@tk_file_dialog_wrapper()
def get_open_file(parent):
    return tkFileDialog.askopenfilename(parent=parent, filetypes=self.types)
```

# Classes

## ADModel
### Class Definition
`ADModel(content: str = None, _sort_atom_by_res: bool = False, _parse2std: bool = False)`

### Initialization Method
Initializes an ADModel instance with optional PDB content. Parses the content if provided.

### Members
- `energy`: Stores the energy of the docked ligand.
- `name`: Stores the name of the ligand.
- `poseN`: Stores the pose number.
- `info`: Stores additional information about the ligand.
- `pdb_string`: Stores the PDB formatted string.
- `pdb_lines`: Stores the lines of the PDB file.

### Methods

#### parse_content(content: str, _sort_atom_by_res: bool = False, _parse2std: bool = False)
##### Method Description
Parses the PDB content and organizes it into structured data.

##### Parameters
- `content`: The PDB content to be parsed.
- `_sort_atom_by_res`: Whether to sort atoms by residue.
- `_parse2std`: Whether to parse the content into standard PDB format.

##### Returns
- None

##### Notes
- This method is crucial for processing the PDB data and making it accessible.

##### Example
```python
ad_model = ADModel()
ad_model.parse_content(some_pdb_content)
```

#### as_pdb_string()
##### Method Description
Returns the PDB content as a formatted string.

##### Parameters
- None

##### Returns
A string representing the PDB content.

##### Notes
- This method is useful for exporting or displaying the PDB data.

##### Example
```python
pdb_str = ad_model.as_pdb_string()
```

#### info_string()
##### Method Description
Returns the information string associated with the model.

##### Parameters
- None

##### Returns
The information string.

##### Notes
- This method provides additional context about the model.

##### Example
```python
info = ad_model.info_string()
```

## DlgFile
### Class Definition
`DlgFile(path: str = None, content: str = None, sort_pdb_line_by_res: bool = False, parse2std: bool = False)`

### Initialization Method
Initializes a DlgFile instance with optional file path or content. Decodes the content into a list of poses.

### Members
- `path`: The file path.
- `sort_pdb_line_by_res`: Whether to sort PDB lines by residue.
- `parse2std`: Whether to parse to standard format.
- `pose_lst`: A list of ADModel instances.
- `n2i`: A dictionary mapping pose names to indices.

### Methods

#### __len__()
##### Method Description
Returns the number of poses.

##### Parameters
- None

##### Returns
The number of poses.

##### Notes
- This method allows the use of len() on DlgFile instances.

##### Example
```python
num_poses = len(dlg_file)
```

#### sort_pose(self, key: Callable[[ADModel], Any] = None, inplace: bool = True, reverse: bool = False)
##### Method Description
Sorts the poses based on a given key.

##### Parameters
- `key`: The function to extract the sorting key from an ADModel instance.
- `inplace`: Whether to sort the poses in place.
- `reverse`: Whether to reverse the sort order.

##### Returns
The sorted list of poses.

##### Notes
- This method is useful for organizing poses based on energy or other criteria.

##### Example
```python
sorted_poses = dlg_file.sort_pose(key=lambda x: x.energy)
```

#### decode_content()
##### Method Description
Decodes the content into a list of ADModel instances.

##### Parameters
- None

##### Returns
A list of ADModel instances.

##### Notes
- This method is essential for converting the raw content into a usable format.

##### Example
```python
poses = dlg_file.decode_content()
```

#### asign_pose_name(pose_names: List[str])
##### Method Description
Assigns names to the poses.

##### Parameters
- `pose_names`: A list of names to be assigned to the poses.

##### Returns
- None

##### Notes
- This method helps in identifying poses by their names.

##### Example
```python
dlg_file.asign_pose_name(['pose1', 'pose2'])
```

#### asign_prop(prop: str, value: List[Any])
##### Method Description
Assigns a property to all poses.

##### Parameters
- `prop`: The property name.
- `value`: A list of values for the property.

##### Returns
- None

##### Notes
- This method is used to add custom properties to poses.

##### Example
```python
dlg_file.asign_prop('custom_prop', [value1, value2])
```

#### set_pose_prop(prop: str, value: Any, pose_name: str = None, pose_idx: int = None)
##### Method Description
Sets a property for a specific pose.

##### Parameters
- `prop`: The property name.
- `value`: The value to be set.
- `pose_name`: The name of the pose.
- `pose_idx`: The index of the pose.

##### Returns
- None

##### Notes
- This method allows setting properties for individual poses.

##### Example
```python
dlg_file.set_pose_prop('energy', 5.0, pose_name='pose1')
```

#### get_pose(pose_name: str = None, pose_idx: int = None)
##### Method Description
Retrieves a specific pose.

##### Parameters
- `pose_name`: The name of the pose.
- `pose_idx`: The index of the pose.

##### Returns
The requested ADModel instance.

##### Notes
- This method provides access to individual poses.

##### Example
```python
pose = dlg_file.get_pose(pose_name='pose1')
```

#### get_pose_prop(prop: str, pose_name: str = None, pose_idx: int = None, default: Any = None)
##### Method Description
Gets a property of a specific pose.

##### Parameters
- `prop`: The property name.
- `pose_name`: The name of the pose.
- `pose_idx`: The index of the pose.
- `default`: The default value if the property is not found.

##### Returns
The value of the property.

##### Notes
- This method retrieves properties of individual poses.

##### Example
```python
energy = dlg_file.get_pose_prop('energy', pose_name='pose1')
```

## MyFileDialog
### Class Definition
`MyFileDialog(types=[("Executable", "*")], initialdir: str = None)`

### Initialization Method
Initializes a MyFileDialog instance with optional file types and initial directory.

### Members
- `initialdir`: The initial directory for the file dialog.
- `types`: A list of file types accepted by the dialog.

### Methods

#### get_open_file(parent)
##### Method Description
Opens a dialog to select an existing file.

##### Parameters
- `parent`: The parent window of the dialog.

##### Returns
The path of the selected file.

##### Notes
- This method is a convenient way to get file paths from the user.

##### Example
```python
file_path = my_file_dialog.get_open_file(parent_window)
```

#### get_save_file(parent)
##### Method Description
Opens a dialog to select a location and name for saving a file.

##### Parameters
- `parent`: The parent window of the dialog.

##### Returns
The path where the file will be saved.

##### Notes
- This method helps in saving files with user interaction.

##### Example
```python
save_path = my_file_dialog.get_save_file(parent_window)
```

#### get_ask_dir(parent)
##### Method Description
Opens a dialog to select a directory.

##### Parameters
- `parent`: The parent window of the dialog.

##### Returns
The path of the selected directory.

##### Notes
- This method allows the user to choose a directory interactively.

##### Example
```python
directory_path = my_file_dialog.get_ask_dir(parent_window)
```


