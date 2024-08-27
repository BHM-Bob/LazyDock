<!--
 * @Date: 2024-08-27 10:15:51
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-08-27 10:21:03
 * @Description: 
-->

# Module Description
This module provides a collection of classes and functions for handling and parsing AutoDock docing result file(dlg) and PDB (Protein Data Bank) formatted strings, which are commonly used to represent the three-dimensional structures of molecules.

# Functions

## tk_file_dialog_wrapper(*args, **kwargs)
### Function Introduction
A decorator function that wraps around the tkinter file dialog functions to provide a consistent interface for opening files, saving files, and asking for directories.

### Parameters
- `args`: Variable length argument list for the tkinter file dialog function.
- `kwargs`: Arbitrary keyword arguments for the tkinter file dialog function.

### Return Value
A wrapped function that creates a Tkinter window, executes the file dialog function, and properly closes the window after the operation.

### Notes
- Ensures that the file dialog is properly managed without leaving any open Tkinter windows.

### Example
```python
@tk_file_dialog_wrapper()
def get_file():
    return tkFileDialog.askopenfilename()
```

# Classes

## ADModel(BaseInfo)
### Class Introduction
A storage class designed to hold information about a docked ligand, including its PDB string representation, energy, and other relevant data.

### Constructor
- `__init__(content: str = None, _sort_atom_by_res: bool = False)`: Initializes a new instance of the ADModel class.

#### Parameters
- `content`: A string containing the PDB content.
- `_sort_atom_by_res`: A boolean indicating whether to sort atoms by residue.

### Members
- `energy`: The free energy of binding.
- `name`: The name of the ligand.
- `poseN`: The pose number.
- `info`: Additional information about the ligand.
- `pdb_string`: The PDB string representation of the ligand.
- `pdb_lines`: A list of lines from the PDB file.

### Methods

#### as_pdb_string()-> str
##### Method Introduction
Returns the PDB string representation of the ligand.

##### Parameters
None

##### Return Value
The PDB string as a single concatenated string.

##### Notes
- Useful for exporting or further processing the PDB data.

##### Example
```python
pdb_str = ad_model_instance.as_pdb_string()
```

#### info_string()-> str
##### Method Introduction
Returns the additional information string associated with the ligand.

##### Parameters
None

##### Return Value
The info string.

##### Notes
- Provides context or metadata about the ligand.

##### Example
```python
# 假设我们有一个包含 PDB 文件路径的列表
pdb_file_paths = ['path/to/pdb1.pdb', 'path/to/pdb2.pdb', 'path/to/pdb3.pdb']

# 创建一个空列表来存储 ADModel 实例
ad_models = []

# 遍历所有文件路径
for file_path in pdb_file_paths:
    # 读取 PDB 文件内容
    with open(file_path, 'r') as file:
        pdb_content = file.read()
    
    # 创建 ADModel 实例并解析内容
    ad_model = ADModel(pdb_content)
    ad_model.parse_content(pdb_content)
    
    # 将实例添加到列表中
    ad_models.append(ad_model)

# 现在 ad_models 列表包含了所有 PDB 文件的数据
# 你可以进行进一步的操作，例如按能量排序
ad_models.sort(key=lambda x: x.energy)

# 打印每个模型的能量值
for model in ad_models:
    print(f"Model name: {model.name}, Energy: {model.energy}")
```

## DlgFile(BaseInfo)
### Class Introduction
A class for handling file operations related to PDB files, including loading, decoding, and sorting the content.

### Constructor
- `__init__(path: str = None, content: str = None, sort_pdb_line_by_res: bool = False)`: Initializes a new instance of the DlgFile class with optional file path or content.

#### Parameters
- `path`: The file path to load the PDB content from.
- `content`: A string containing the PDB content.
- `sort_pdb_line_by_res`: A boolean indicating whether to sort PDB lines by residue.

### Members
- `path`: The file path of the PDB file.
- `sort_pdb_line_by_res`: A flag for sorting PDB lines by residue.
- `pose_lst`: A list of ADModel instances representing the poses.
- `n2i`: A dictionary for mapping pose names to indices.

### Methods

#### __len__()-> int
##### Method Introduction
Returns the number of poses in the pose list.

##### Parameters
None

##### Return Value
The count of poses.

##### Notes
- Provides a quick way to determine the number of docked ligands.

##### Example
```python
num_poses = len(dlg_file_instance)
```

#### sort_pose(key: Callable[[ADModel], Any] = None, inplace: bool = True, reverse: bool = False) -> List[ADModel]
##### Method Introduction
Sorts the pose list based on a given key function.

##### Parameters
- `key`: A function used to extract a comparison key from each pose.
- `inplace`: A boolean indicating whether to sort in place.
- `reverse`: A boolean indicating the sort order.

##### Return Value
The sorted pose list.

##### Notes
- Allows for flexible sorting based on various criteria such as energy.

##### Example
```python
sorted_poses = dlg_file_instance.sort_pose(key=lambda x: x.energy, reverse=True)
```

#### decode_content()-> List[ADModel]
##### Method Introduction
Decodes the content into a list of poses.

##### Parameters
None

##### Return Value
A list of ADModel instances.

##### Notes
- Parses the PDB content to create a structured list of poses.

##### Example
```python
poses = dlg_file_instance.decode_content()
```

#### asign_pose_name(pose_names: List[str])
##### Method Introduction
Assigns names to the poses based on the provided list.

##### Parameters
- `pose_names`: A list of names for the poses.

##### Return Value
None

##### Notes
- Ensures that the number of names matches the number of poses.

##### Example
```python
dlg_file_instance.asign_pose_name(['Pose1', 'Pose2'])
```

#### asign_prop(prop: str, value: List[Any])
##### Method Introduction
Assigns a property to all poses with the given value.

##### Parameters
- `prop`: The property name to assign.
- `value`: A list of values for the property.

##### Return Value
None

##### Notes
- Used to set properties such as energy or other metrics for all poses.

##### Example
```python
dlg_file_instance.asign_prop('energy', [0.5, -0.3])
```

#### set_pose_prop(prop: str, value: Any, pose_name: str = None, pose_idx: int = None)
##### Method Introduction
Sets a property for a specific pose by name or index.

##### Parameters
- `prop`: The property name to set.
- `value`: The value to set for the property.
- `pose_name`: The name of the pose.
- `pose_idx`: The index of the pose.

##### Return Value
None

##### Notes
- Allows for setting properties on individual poses.

##### Example
```python
dlg_file_instance.set_pose_prop('energy', -0.2, pose_name='Pose1')
```

#### get_pose(pose_name: str = None, pose_idx: int = None)
##### Method Introduction
Retrieves a pose by name or index.

##### Parameters
- `pose_name`: The name of the pose.
- `pose_idx`: The index of the pose.

##### Return Value
The ADModel instance representing the pose.

##### Notes
- Provides access to a specific pose within the list.

##### Example
```python
pose = dlg_file_instance.get_pose(pose_name='Pose1')
```

#### get_pose_prop(prop: str, pose_name: str = None, pose_idx: int = None, default: Any = None)
##### Method Introduction
Gets a property value for a specific pose by name or index.

##### Parameters
- `prop`: The property name to retrieve.
- `pose_name`: The name of the pose.
- `pose_idx`: The index of the pose.
- `default`: The default value if the property is not set.

##### Return Value
The property value for the specified pose.

##### Notes
- Retrieves individual property values for poses.

##### Example
```python
energy = dlg_file_instance.get_pose_prop('energy', pose_name='Pose1')
```

## MyFileDialog
### Class Introduction
A class that encapsulates file dialog operations, providing methods to open files, save files, and ask for directories with a consistent interface.

### Constructor
- `__init__(types=[("Executable", "*")], initialdir: str = None)`: Initializes a new instance of the MyFileDialog class with optional file types and initial directory.

#### Parameters
- `types`: A list of file types to filter by.
- `initialdir`: The initial directory to start the file dialog in.

### Methods

#### get_open_file(parent)
##### Method Introduction
Opens a dialog to select an existing file.

##### Parameters
- `parent`: The parent Tkinter window for the dialog.

##### Return Value
The path of the selected file or None if the dialog was cancelled.

##### Notes
- Utilizes the wrapped tkinter file dialog function.

##### Example
```python
file_path = my_file_dialog_instance.get_open_file()
```

#### get_save_file(parent)
##### Method Introduction
Opens a dialog to select a location to save a file.

##### Parameters
- `parent`: The parent Tkinter window for the dialog.

##### Return Value
The path to save the file or None if the dialog was cancelled.

##### Notes
- Utilizes the wrapped tkinter file dialog
