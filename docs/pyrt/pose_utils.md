# Module Description

This module provides a class `_Pose` for handling and manipulating protein structures using the PyRosetta API. It includes methods for initializing a pose from a sequence, PDB file, or RCSB identifier, adding or deleting residues, mutating residues, managing chains, and splitting or merging chains. The class encapsulates various PyRosetta functionalities to facilitate common tasks in protein structure analysis and modification.

# Functions

## renumber_pdbinfo_based_on_conf_chains(pose)
### Function Introduction
Renumbers the PDB information of a pose based on the actual chains present in the conformation.

### Parameters
- `pose`: The PyRosetta `Pose` object to renumber.

### Return Value
None

### Notes
- Ensures that the PDB information is consistent with the pose's chain configuration.

### Example
```python
renumber_pdbinfo_based_on_conf_chains(pose)
```

# Classes

## _Pose
### Class Introduction
A class that encapsulates PyRosetta pose operations, providing a high-level interface for manipulating protein structures.

### Constructor
- `__init__(seq: str = None, pdb_path: str = None, rcsb: str = None)`: Initializes a new `_Pose` instance from a sequence string, a PDB file path, or an RCSB identifier.

#### Parameters
- `seq`: A sequence string to create the pose.
- `pdb_path`: A file path to a PDB file to load the pose.
- `rcsb`: An RCSB identifier to fetch and create the pose.

#### Return Value
None

### Members
- `pose`: The PyRosetta `Pose` object representing the protein structure.
- `pdb_info`: The PDB information associated with the pose.

### Methods

#### add_resi(resi_name3: str, insert_pos: int, pos: str = 'after') -> '_Pose'
##### Method Introduction
Adds a residue to the pose at a specified position.

##### Method Parameters
- `resi_name3`: The three-letter name of the residue to add.
- `insert_pos`: The position in the sequence to insert the residue.
- `pos`: The position relative to the `insert_pos` ('before' or 'after').

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Modifies the pose in place.

##### Example
```python
pose_instance.add_resi('ALA', 5, 'after')
```

#### delete_resi(resi_id: int) -> '_Pose'
##### Method Introduction
Deletes a residue from the pose by its identifier.

##### Method Parameters
- `resi_id`: The identifier of the residue to delete.

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Modifies the pose in place.

##### Example
```python
pose_instance.delete_resi(5)
```

#### mutate_resi(resi_id: int, new_name1: str, pack_radius: float = 0.0)
##### Method Introduction
Mutates a residue to a different type.

##### Method Parameters
- `resi_id`: The identifier of the residue to mutate.
- `new_name1`: The new residue type.
- `pack_radius`: The packing radius for the mutation.

##### Method Return Value
None

##### Notes
- Mutates the residue and triggers a packing operation.

##### Example
```python
pose_instance.mutate_resi(5, 'GLY', 1.0)
```

#### add_chain(chain: Union['_Pose', Pose], jump_resi_id: int = None) -> '_Pose'
##### Method Introduction
Adds another chain to the pose.

##### Method Parameters
- `chain`: The chain to add, which can be a `_Pose` or `Pose` instance.
- `jump_resi_id`: The residue identifier where the new chain will be attached.

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Connects the new chain at the specified residue.

##### Example
```python
pose_instance.add_chain(other_pose_instance)
```

#### del_region(start_resi_id: int, stop_resi_id: int) -> '_Pose'
##### Method Introduction
Deletes a region of the pose between two residue identifiers.

##### Method Parameters
- `start_resi_id`: The starting residue identifier for the deletion.
- `stop_resi_id`: The ending residue identifier for the deletion.

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Removes the specified region from the pose.

##### Example
```python
pose_instance.del_region(5, 10)
```

#### merge_chain(chain: Union['_Pose', Pose], new_chain: bool = False) -> '_Pose'
##### Method Introduction
Merges another pose into the current one.

##### Method Parameters
- `chain`: The pose to merge, which can be a `_Pose` or `Pose` instance.
- `new_chain`: Whether to start a new chain for the merged pose.

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Combines the poses into a single structure.

##### Example
```python
pose_instance.merge_chain(other_pose_instance)
```

#### delete_chain(chain_id: int = None, chain_name: str = None) -> '_Pose'
##### Method Introduction
Deletes a chain from the pose.

##### Method Parameters
- `chain_id`: The identifier of the chain to delete.
- `chain_name`: The name of the chain to delete.

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Removes the specified chain from the pose.

##### Example
```python
pose_instance.delete_chain(chain_id=1)
```

#### swap_chain(order: str | list[str]) -> '_Pose'
##### Method Introduction
Swaps the order of chains in the pose.

##### Method Parameters
- `order`: A string or list specifying the new order of chains.

##### Method Return Value
The current instance of `_Pose`.

##### Notes
- Rearranges the chains according to the provided order.

##### Example
```python
pose_instance.swap_chain(['A', 'B'])
```

#### split_chain(return_Pose: bool = False) -> list['_Pose'] | list[Pose]
##### Method Introduction
Splits the pose into individual chains.

##### Method Parameters
- `return_Pose`: Whether to return a list of `_Pose` instances or PyRosetta `Pose` objects.

##### Method Return Value
A list of `_Pose` instances or PyRosetta `Pose` objects, depending on `return_Pose`.

##### Notes
- Divides the pose into separate chains.

##### Example
```python
pose_list = pose_instance.split_chain(return_Pose=True)
```

#### get_resi_id_in_pose_via_pdb(chain_name: str, resi_id: int) -> int
##### Method Introduction
Retrieves the internal pose residue identifier based on the PDB chain name and identifier.

##### Method Parameters
- `chain_name`: The PDB chain name.
- `resi_id`: The PDB residue identifier.

##### Method Return Value
The internal pose residue identifier.

##### Notes
- Maps PDB identifiers to internal pose identifiers.

##### Example
```python
internal_id = pose_instance.get_resi_id_in_pose_via_pdb('A', 5)
```

#### get_resi_id_in_pdb_via_pose(resi_id: int) -> tuple[str, int]
##### Method Introduction
Retrieves the PDB chain name and identifier based on the internal pose residue identifier.

##### Method Parameters
- `resi_id`: The internal pose residue identifier.

##### Method Return Value
A tuple containing the PDB chain name and identifier.

##### Notes
- Maps internal pose identifiers to PDB identifiers.

##### Example
```python
chain_name, pdb_id = pose_instance.get_resi_id_in_pdb_via_pose(5)
```

#### get_chain_id(chain_name: str) -> int
##### Method Introduction
Retrieves the internal pose chain identifier based on the chain name.

##### Method Parameters
- `chain_name`: The chain name.

##### Method Return Value
The internal pose chain identifier.

##### Notes
- Maps chain names to internal pose identifiers.

##### Example
```python
chain_id = pose_instance.get_chain_id('A')
```

#### get_chain_name(chain_id: int = None, resi_id: int = None) -> str
##### Method Introduction
Retrieves the chain name from the pose based on either a chain identifier or a residue identifier.

##### Method Parameters
- `chain_id`: The internal pose chain identifier.
- `resi_id`: The internal pose residue identifier.

##### Method Return Value
The chain name.

##### Notes
- Retrieves the chain name associated with the given identifier.

##### Example
```python
chain_name = pose_instance.get_chain_name(chain_id=1)
```

#### get_chain_seq(chain_id: int = None, chain_name: str = None) -> str
##### Method Introduction
Retrieves the sequence of a chain in the pose.

##### Method Parameters
- `chain_id`: The internal pose chain identifier.
- `chain_name`: The chain name.

##### Method Return Value
The sequence string of the specified chain.

##### Notes
- Provides the amino acid sequence of the chain.

##### Example
```python
chain_sequence = pose_instance.get_chain_seq(chain_name='A')
```

The class `_Pose` and the function `renumber_pdbinfo_based_on_conf_chains` are the main components of this module, offering a comprehensive set of tools for working with protein structures in PyRosetta.

