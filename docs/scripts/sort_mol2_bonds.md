# Mol2 Bond Sorter CLI Tool

## Overview
This Python script is a command-line interface (CLI) tool designed to sort the bonds in Mol2 files based on atom indices. It can process individual Mol2 files or directories containing multiple Mol2 files.

## Command: sort_mol2_bonds
### Brief
The `sort_mol2_bonds` command sorts the bonds in Mol2 files by atom index, ensuring that the bonds are listed in the correct order.

### Arguments
- `-i`, `--input`: The input Mol2 file or directory path. If a directory is provided, the script will process all Mol2 files within it. Defaults to the current working directory.
- `-s`, `--suffix`: The suffix to append to the sorted Mol2 file names. Defaults to "_sorted".
- `-r`, `--recursive`: A flag indicating whether to search the input directory recursively. Defaults to False.

### Behavior
1. The script processes the input path to determine if it is a single Mol2 file or a directory containing Mol2 files.
2. If the input is a directory, it retrieves a list of Mol2 files within the directory (and its subdirectories if recursive search is enabled).
3. For each Mol2 file, it sorts the bonds based on atom indices using the `sort_bonds` function.
4. The sorted Mol2 content is then saved to a new file with the specified suffix appended to the original file name.

### Notes
- Ensure that the input path is valid and contains the desired Mol2 files.
- The script assumes that the necessary dependencies, including the `mbapy_lite` and `lazydock` libraries, are installed and properly configured.
- The `sort_bonds` function is responsible for the actual bond sorting logic, which is not detailed here.

### Example
```bash
lazydock-cli sort_mol2_bonds -i /path/to/mol2_files -s _sorted_bonds -r
```