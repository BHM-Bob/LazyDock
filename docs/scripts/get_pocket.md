<!--
 * @Date: 2024-12-02 17:02:55
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-12-02 17:06:00
 * @Description: 
-->
## Command: Get Pocket Residues from ProteinPlus

### Overview
This script is designed to retrieve pocket residues from the ProteinPlus web server. It processes input receptor files (in mol2 format) or directories containing such files, and optionally includes ligand information (in sdf format). The script outputs the pocket residues in various formats, including JSON and PSE files, and can use different methods for pocket calculation.

### Parameters

- `-r`, `--receptor`: The input mol2 file or directory path. Defaults to the current directory `.`.
- `-l`, `--ligand`: The ligand sdf file path. Defaults to `None`.
- `-o`, `--output`: The output directory. Defaults to `None`.
- `-m`, `--method`: The pocket calculation method. Can be either 'extend' or 'mean'. Defaults to 'extend'.

### Behavior
The script performs the following actions:
1. Processes the input receptor path, which can be a single file or a directory containing multiple files.
2. If a directory is provided, it searches for `.pdb` files within that directory.
3. Displays the parsed arguments for verification.
4. For each receptor file, it creates an output directory and copies the receptor file into it.
5. Retrieves the pocket box from the ProteinPlus web server using the provided ligand information.
6. Parses the pocket box and saves the pocket residues in JSON format.
7. Saves the pocket box in PSE format using PyMOL.

### Notes
- The script assumes that the input receptor path is either a valid mol2 file or a directory containing mol2 files.
- If the output directory is not specified, the script will create a new directory for each receptor file with the suffix `_pocket_box`.
- The script requires internet access to retrieve data from the ProteinPlus web server.
- The `method` parameter affects how the pocket residues are calculated and parsed.

### Examples

To retrieve pocket residues for a single receptor file and use the 'mean' method:
```bash
lazydock-cli get-pocket --receptor path/to/receptor.pdb --method mean
```

To retrieve pocket residues for a directory of receptor files:
```bash
lazydock-cli get-pocket --receptor path/to/receptors/
```

To retrieve pocket residues with ligand information:
```bash
lazydock-cli get-pocket --receptor path/to/receptor.pdb --ligand path/to/ligand.sdf
```
