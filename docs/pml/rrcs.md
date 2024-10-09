<!--
 * @Date: 2024-10-09 21:02:25
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-10-09 21:03:47
 * @Description: 
-->
# Description

This module is designed to calculate the Relative Residue Contact Score (RRCS) in PyMOL, which is a measure used in the study of protein structures, particularly in the analysis of GPCR activation mechanisms. The RRCS score quantifies the degree of contact between residue pairs in a protein structure.

# Functions

## _test_close(dict_coord: Dict[str, Dict[int, Tuple[float, float, float, float]]], ires: str, jres: str) -> bool
### Function Description
Tests if any atom in residue `ires` is within 4.63 Angstroms of any atom in residue `jres`.

### Parameters
- `dict_coord`: A dictionary containing the coordinates of atoms.
- `ires`: The identifier for the first residue.
- `jres`: The identifier for the second residue.

### Returns
A boolean indicating whether any close contact is found.

### Notes
- This function is used to preliminarily filter residue pairs that are close enough to potentially interact.

## _calcu_score(dict_coord: Dict[str, Dict[int, Tuple[float, float, float, float]]], atomnum2name: Dict[int, str], ires: str, jres: str, check_name: bool, score_count: int) -> float
### Function Description
Calculates the RRCS score for a pair of residues based on the distance between their atoms.

### Parameters
- `dict_coord`: A dictionary containing the coordinates of atoms.
- `atomnum2name`: A dictionary mapping atom numbers to atom names.
- `ires`: The identifier for the first residue.
- `jres`: The identifier for the second residue.
- `check_name`: A boolean indicating whether to skip certain atom types.
- `score_count`: A reference to a counter for the number of scored pairs.

### Returns
The calculated RRCS score for the residue pair.

### Notes
- This function implements the scoring logic as described in the RRCS methodology.

## calcu_RRCS(model: str) -> pd.DataFrame
### Function Description
Calculates the RRCS for all residue pairs in a given model and returns the results in a DataFrame.

### Parameters
- `model`: The PyMOL model for which to calculate RRCS.

### Returns
A pandas DataFrame containing the RRCS scores for each residue pair.

### Notes
- This function coordinates the overall RRCS calculation process.

### Example
```python
from pymol import cmd
cmd.reinitialize()
cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
rrcs_scores = calcu_RRCS('receptor')
```