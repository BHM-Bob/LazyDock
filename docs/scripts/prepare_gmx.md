<!--
 * @Date: 2024-12-28 21:20:15
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-12-28 21:20:30
 * @Description: 
-->
## Command: protein

### Introduction
The `protein` command is designed to prepare a single protein for GROMACS Molecular Dynamics Simulations (MDS). It involves a series of steps to ensure the protein is properly centered, aligned, and has its topology prepared.

### Parameters
- `-d`, `--dir`: The directory containing the protein files. Default is the current directory.
- `-n`, `--protein-name`: The filename of the protein within each sub-directory.
- `--ff-dir`: The directory containing force field files.
- `--n-term`: The N-Term type for gmx pdb2gmx. Default is 'auto'.
- `--pdb2gmx-args`: Arguments passed to the pdb2gmx command. Default is "-ter -ignh".

### Behavior
1. The command searches for protein files within the specified directory.
2. It centers the protein structure using Open Babel.
3. The protein structure is then aligned with the xyz axes using LazyDock.
4. The protein topology is prepared using GROMACS pdb2gmx.

### Notes
- The directory argument must be a valid directory.
- The protein file name must be specified to identify the correct files.
- Force field files are required for topology preparation and must be in the specified directory.

### Example
```bash
lazydock-cli prepare-gmx protein -d /path/to/protein/directory -n protein.pdb --ff-dir /path/to/ff
```

## Command: ligand

### Introduction
The `ligand` command prepares a single ligand for GROMACS MDS. It includes steps for centering, aligning, converting to mol2 format, fixing names, sorting bonds, retrieving force field parameters, and preparing the ligand topology.

### Parameters
- `-d`, `--dir`: The directory containing the ligand files. Default is the current directory.
- `-n`, `--ligand-name`: The filename of the ligand within each sub-directory.
- `--ff-dir`: The directory containing force field files.
- `--max-step`: The maximum step to perform. Default is 7.
- `--disable-browser`: A flag to disable the browser for CGenFF.

### Behavior
1. The command searches for ligand files within the specified directory.
2. It centers the ligand structure using Open Babel.
3. The ligand structure is aligned with the xyz axes using LazyDock.
4. The ligand is converted to mol2 format.
5. The ligand name and residue name in the mol2 file are fixed.
6. Bonds in the mol2 file are sorted using LazyDock.
7. Force field parameters are retrieved from CGenFF.
8. The ligand topology is prepared using cgenff_charmm2gmx.py.

### Notes
- The directory argument must be a valid directory.
- The ligand file name must be specified to identify the correct files.
- Force field files are required for topology preparation and must be in the specified directory.
- The command may require a browser instance for CGenFF, which can be disabled with the `--disable-browser` flag.

### Example
```bash
lazydock-cli prepare-gmx ligand -d /path/to/ligand/directory -n ligand.pdb --ff-dir /path/to/ff
```

## Command: complex

### Introduction
The `complex` command prepares a protein-ligand complex for GROMACS MDS. It involves steps for centering, aligning, extracting receptor and ligand, preparing topologies, and merging files.

### Parameters
- `-d`, `--dir`: The directory containing the complex files. Default is the current directory.
- `-n`, `--complex-name`: The filename of the complex within each sub-directory.
- `--max-step`: The maximum step to perform. Default is 9.
- `--receptor-chain-name`: The name of the receptor chain.
- `--ligand-chain-name`: The name of the ligand chain.
- `--ff-dir`: The directory containing force field files.
- `--disable-browser`: A flag to disable the browser for CGenFF.

### Behavior
1. The command searches for complex files within the specified directory.
2. It centers the complex structure using Open Babel.
3. The complex structure is aligned with the xyz axes using LazyDock.
4. The receptor and ligand are extracted from the complex.
5. The ligand is prepared using the `ligand` command.
6. The protein topology is prepared using GROMACS pdb2gmx.
7. The ligand topology is prepared.
8. The receptor and ligand gro files are merged into a complex.gro file.
9. The topol.top file is prepared by merging receptor and ligand parameters.

### Notes
- The directory argument must be a valid directory.
- The complex file name must be specified to identify the correct files.
- Force field files are required for topology preparation and must be in the specified directory.
- The receptor and ligand chain names must be provided to extract the respective components from the complex.
- The command may require a browser instance for CGenFF, which can be disabled with the `--disable-browser` flag.

### Example
```bash
lazydock-cli prepare-gmx complex -d /path/to/complex/directory -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir /path/to/ff
```
