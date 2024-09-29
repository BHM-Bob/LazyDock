# Description
This module contains a collection of functions designed to interface with various online protein structure evaluation servers. It allows users to submit Protein Data Bank (PDB) files to these servers and retrieve evaluation scores and other information.

# Functions

## `get_score_from_proq(pdb_path: str, **kwargs) -> Dict[str, float]`
### Description
Returns a dictionary containing LGscore and MaxSub scores from the ProQ server.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.

### Returns
A dictionary with keys 'LGscore' and 'MaxSub', each associated with a float value.

### Notes
- Requires an active internet connection to submit to the ProQ server.

### Example
```python
lgscore, maxsub = get_score_from_proq('path/to/your/pdbfile.pdb')['LGscore'], get_score_from_proq('path/to/your/pdbfile.pdb')['MaxSub']
```

## `get_score_from_VoroMQA(pdb_path: str, browser: Browser = None, timeout: int = 1200, **kwargs) -> Dict[str, float]`
### Description
Returns a dictionary containing the VoroMQA score.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with the key 'Score' and its associated float value.

### Notes
- Requires an active internet connection to submit to the VoroMQA server.

### Example
```python
score = get_score_from_VoroMQA('path/to/your/pdbfile.pdb')['Score']
```

## `get_score_from_ProSA(pdb_path: str, browser: Browser = None, timeout: int = 1200, **kwargs) -> Dict[str, Union[float, bytes]]`
### Description
Returns a dictionary containing the Z-Score, model_quality_img_res, and res_score_img_res from the ProSA server.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with keys 'Z-Score', 'model_quality_img_res', and 'res_score_img_res'.

### Notes
- Requires an active internet connection to submit to the ProSA server.

### Example
```python
z_score, model_quality_img, res_score_img = get_score_from_ProSA('path/to/your/pdbfile.pdb')['Z-Score'], get_score_from_ProSA('path/to/your/pdbfile.pdb')['model_quality_img_res'], get_score_from_ProSA('path/to/your/pdbfile.pdb')['res_score_img_res']
```

## `get_score_from_ModEval(modkey: str, pdb_path: str, align_file_path: str, browser: Browser = None, timeout: int = 1800, **kwargs) -> Dict[str, Union[float, str]]`
### Description
Returns a dictionary containing RMSD, overlap, identity, and Z-DOPE scores from the ModEval server.

### Parameters
- `modkey`: The ModEval key for authentication.
- `pdb_path`: The path to the PDB file to be evaluated.
- `align_file_path`: The path to the alignment file.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with keys 'RMSD', 'overlap', 'identity', and 'Z-DOPE'.

### Notes
- Requires an active internet connection to submit to the ModEval server.

### Example
```python
rmsd, overlap, identity, z_dope = get_score_from_ModEval('your_modkey', 'path/to/your/pdbfile.pdb', 'path/to/your/alignfile.aln')['RMSD'], get_score_from_ModEval('your_modkey', 'path/to/your/pdbfile.pdb', 'path/to/your/alignfile.aln')['overlap'], get_score_from_ModEval('your_modkey', 'path/to/your/pdbfile.pdb', 'path/to/your/alignfile.aln')['identity'], get_score_from_ModEval('your_modkey', 'path/to/your/pdbfile.pdb', 'path/to/your/alignfile.aln')['Z-DOPE']
```

## `get_score_from_MolProbity(pdb_path: str, browser: Browser = None, timeout: int = 1800, **kwargs) -> Dict[str, Union[float, str]]`
### Description
Returns a dictionary containing the MolProbity Z-Score, Ramachandran Favorability, and Ramachandran Outerness.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with keys 'Z-Score', 'Ramachandran Favorability', and 'Ramachandran Outerness'.

### Notes
- Requires an active internet connection to submit to the MolProbity server.

### Example
```python
z_score, rama_favor, rama_outer = get_score_from_MolProbity('path/to/your/pdbfile.pdb')['Z-Score'], get_score_from_MolProbity('path/to/your/pdbfile.pdb')['Ramachandran Favorability'], get_score_from_MolProbity('path/to/your/pdbfile.pdb')['Ramachandran Outerness']
```

## `get_score_from_ProQ3(pdb_path: str, browser: Browser = None, timeout: int = 1200, **kwargs) -> Dict[str, Union[float, str]]`
### Description
Returns a dictionary containing ProQ2D, ProQRosCenD, ProQRosCenD, and ProQ3D scores from the ProQ3 server.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with keys 'ProQ2D', 'ProQRosCenD', 'ProQRosCenD', and 'ProQ3D'.

### Notes
- Requires an active internet connection to submit to the ProQ3 server.

### Example
```python
proq2d, proqroscend, proqroscend2, proq3d = get_score_from_ProQ3('path/to/your/pdbfile.pdb')['ProQ2D'], get_score_from_ProQ3('path/to/your/pdbfile.pdb')['ProQRosCenD'], get_score_from_ProQ3('path/to/your/pdbfile.pdb')['ProQRosCenD'], get_score_from_ProQ3('path/to/your/pdbfile.pdb')['ProQ3D']
```

## `get_score_from_SAVES(pdb_path: str, browser: Browser = None, timeout: int = 1200, **kwargs) -> Dict[str, Union[float, str]]`
### Description
Returns a dictionary containing SAVES evaluation scores.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with keys 'errat', 'whatcheck_bad', 'whatcheck_mid', 'whatcheck_good', 'procheck_errors', 'procheck_warnings', and 'procheck_pass'.

### Notes
- Requires an active internet connection to submit to the SAVES server.

### Example
```python
results = get_score_from_SAVES('path/to/your/pdbfile.pdb')
```

## `get_score_from_QMEANDisCo(pdb_path: str, browser: Browser = None, timeout: int = 1200, **kwargs)`
### Description
Returns a dictionary containing the QMEANDisCo score from the SwissModel server.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with the key 'QMEANDisCo Global'.

### Notes
- Requires an active internet connection to submit to the SwissModel server.

### Example
```python
qmean_disco_score = get_score_from_QMEANDisCo('path/to/your/pdbfile.pdb')['QMEANDisCo Global']
```

## `get_score_from_QMEAN(pdb_path: str, browser: Browser = None, timeout: int = 1200, **kwargs)`
### Description
Returns a dictionary containing the QMEAN score from the SwissModel server.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `browser`: An optional `Browser` object for automation.
- `timeout`: The timeout in seconds for waiting for results.

### Returns
A dictionary with the key 'QMEAN4'.

### Notes
- Requires an active internet connection to submit to the SwissModel server.

### Example
```python
qmean_score = get_score_from_QMEAN('path/to/your/pdbfile.pdb')['QMEAN4']
```

## `get_eval_info_from_servers(pdb_path: str, align_file_path: str = None, modkey: str = None, browser: Browser = None, servers: List[str] = None) -> Dict[str, Dict[str, Union[float, str, bytes]]]`
### Description
Consolidates evaluation results from multiple servers into a single dictionary.

### Parameters
- `pdb_path`: The path to the PDB file to be evaluated.
- `align_file_path`: Optional path to the alignment file (used by some servers).
- `modkey`: Optional ModEval key for authentication.
- `browser`: An optional `Browser` object for automation.
- `servers`: A list of server names to retrieve scores from.

### Returns
A dictionary where each key is a server name and each value is a dictionary of scores returned by that server.

### Notes
- If no servers are specified, it returns an error message.
- Requires an active internet connection to submit to the specified servers.

### Example
```python
servers = ['ProQ', 'VoroMQA', 'ProSA']
eval_results = get_eval_info_from_servers('path/to/your/pdbfile.pdb', servers=servers)
```

# Classes

## `_name2server`
### Description
A dictionary mapping server names to their respective scoring functions.

### Members
- `ProQ`: Maps to `get_score_from_proq`.
- `VoroMQA`: Maps to `get_score_from_VoroMQA`.
- `ProSA`: Maps to `get_score_from_ProSA`.
- `ModEval`: Maps to `get_score_from_ModEval`.
- `MolProbity`: Maps to `get_score_from_MolProbity`.
- `ProQ3`: Maps to `get_score_from_ProQ3`.
- `SAVES`: Maps to `get_score_from_SAVES`.
- `QMEANDisCo`: Maps to `get_score_from_QMEANDisCo`.
- `QMEAN`: Maps to `get_score_from_QMEAN`.

### Example
```python
# Example of using the _name2server dictionary to get the ProQ function
proq_function = _name2server['ProQ']
result = proq_function('path/to/your/pdbfile.pdb')
```

This module is a comprehensive tool for protein structure evaluation, automating the process of submitting PDB files to various online servers and collecting the results for further analysis.
