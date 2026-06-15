<!--
 * @Date: 2024-12-28 21:20:15
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2025-06-12
 * @Description: prepare_gmx 脚本文档
-->

## Command: protein

### Introduction
The `protein` command is designed to prepare a single protein for GROMACS Molecular Dynamics Simulations (MDS). It centers the protein, aligns it with axes, and prepares the topology using pdb2gmx.

### Parameters
| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | The directory containing the protein files. |
| `-n` | `--protein-name` | str | None | The filename of the protein within each sub-directory. |
| | `--ff-dir` | str | None | The directory containing force field files. |
| | `--n-term` | str (+) | `['0']` | N-Term type for gmx pdb2gmx. Use `auto` to automatically detect (1 if N-term is MET, else 0). Supports multiple values for multi-chain proteins. |
| | `--c-term` | str (+) | `['0']` | C-Term type for gmx pdb2gmx. Use `auto` to automatically detect (1 if C-term is MET, else 0). Supports multiple values for multi-chain proteins. |
| | `--chain-num` | int | `1` | Number of chains in the protein. |
| | `--pdb2gmx-args` | str | `"-ter -ignh"` | Arguments passed to the pdb2gmx command. |

### Behavior
1. **STEP 0.1**: Center the protein structure using Open Babel (`obabel -c`).
2. **STEP 0.2**: Align the protein structure with the xyz axes using LazyDock.
3. **STEP 1**: Prepare the protein topology using GROMACS pdb2gmx with terminal handling.

### Notes
- The directory argument must be a valid directory.
- The protein file name must be specified to identify the correct files.
- Force field files are required for topology preparation.
- For multi-chain proteins, use `--chain-num` to specify the number of chains, and provide `--n-term` and `--c-term` values for each chain.
- When using `auto` for terminals, the program checks if the terminal residue is Methionine (M) to decide the protonation state.

### Example
```bash
# Single chain protein
lazydock-cli prepare-gmx protein -d /path/to/protein/directory -n protein.pdb --ff-dir /path/to/ff

# Multi-chain protein with auto terminal detection
lazydock-cli prepare-gmx protein -d /path/to/protein/directory -n protein.pdb --ff-dir /path/to/ff --chain-num 2 --n-term auto auto --c-term auto auto

# Specify terminal types explicitly
lazydock-cli prepare-gmx protein -d /path/to/protein/directory -n protein.pdb --ff-dir /path/to/ff --n-term 1 --c-term 0
```

---

## Command: ligand

### Introduction
The `ligand` command prepares a single ligand for GROMACS MDS. It includes steps for centering, aligning, converting to mol2 format, fixing names, sorting bonds, retrieving force field parameters from CGenFF, and preparing the ligand topology.

### Parameters
| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | The directory containing the ligand files. |
| `-n` | `--ligand-name` | str | None | The filename of the ligand within each sub-directory. |
| | `--ff-dir` | str | None | Force field files directory. If before CGenFF step, can be an absolute path to charmmFF dir; if after CGenFF step, should be charmmFF dir name in each sub-directory. |
| | `--max-step` | int | `8` | The maximum step to perform. |
| | `--disable-browser` | flag | False | Whether to disable the browser for CGenFF (for offline mode or when CGenFF files already exist). |

### Behavior
1. **STEP 0**: Center the ligand structure using Open Babel.
2. **STEP 1**: Alter chain code to "Z" and align with xyz axes using LazyDock.
3. **STEP 2**: Transfer ligand.pdb to mol2 format by Open Babel.
4. **STEP 3**: Fix ligand name and residue name in the mol2 file.
5. **STEP 4**: Sort mol2 bonds by LazyDock.
6. **STEP 5**: Retrieve str file from CGenFF (requires browser login).
7. **STEP 6**: Transfer str file to top and gro file by cgenff_charmm2gmx.py.
8. **STEP 7**: Prepare the ligand gro file using editconf.
9. **STEP 8**: Copy ligand topology to system topology (topol.top).

### Notes
- The directory argument must be a valid directory.
- The ligand file name must be specified.
- CGenFF account credentials must be configured in `~/.lazydock/config.json`.
- Use `--disable-browser` if you already have CGenFF files or want to skip the download step.

### Example
```bash
# With CGenFF download
lazydock-cli prepare-gmx ligand -d /path/to/ligand/directory -n ligand.pdb --ff-dir /path/to/charmm36-jul2022.ff

# Skip browser (offline mode)
lazydock-cli prepare-gmx ligand -d /path/to/ligand/directory -n ligand.pdb --ff-dir /path/to/charmm36-jul2022.ff --disable-browser
```

---

## Command: complex

### Introduction
The `complex` command prepares a protein-ligand complex for GROMACS MDS. It involves steps for centering, aligning, extracting receptor and ligand, preparing topologies using CGenFF for the ligand, and merging files.

### Parameters
| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | The directory containing the complex files. |
| `-n` | `--complex-name` | str | None | The filename of the complex within each sub-directory. |
| | `--max-step` | int | `9` | The maximum step to perform. |
| `-rc` | `--receptor-chain-name` | str | None | The name of the receptor chain. |
| `-lc` | `--ligand-chain-name` | str | None | The name of the ligand chain. |
| | `--ff-dir` | str | None | Force field files directory. |
| | `--disable-browser` | flag | False | Whether to disable the browser for CGenFF. |
| | `--pdb2gmx-args` | str | `"-ter -ignh"` | Arguments passed to pdb2gmx command. |
| | `--n-term` | str | `"0"` | N-Term type for gmx pdb2gmx. Use `auto` for automatic detection. |
| | `--c-term` | str | `"0"` | C-Term type for gmx pdb2gmx. Use `auto` for automatic detection. |

### Behavior
1. **STEP 0.1**: Center complex.pdb by Open Babel.
2. **STEP 0.2**: Align complex.pdb with xyz axes by LazyDock.
3. **STEP 1**: Extract receptor and ligand from complex.pdb using chain names.
4. **STEP 2-6**: Prepare ligand topology (same as `ligand` command).
5. **STEP 7**: Prepare the protein topology using pdb2gmx with terminal handling.
6. **STEP 8**: Prepare the ligand gro file.
7. **STEP 9**: Merge receptor and ligand gro files into complex.gro and prepare topol.top.

### Notes
- The input complex.pdb should have two chains, one for receptor and one for ligand.
- The complex.pdb should already have hydrogens added by Avogadro or other software.
- The receptor and ligand chain names must be provided to extract components.
- Supports `auto` mode for terminal detection (assumes single protein receptor).

### Example
```bash
# Basic usage
lazydock-cli prepare-gmx complex -d /path/to/complex/directory -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir /path/to/ff

# With auto terminal detection
lazydock-cli prepare-gmx complex -d /path/to/complex/directory -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir /path/to/ff --n-term auto --c-term auto
```

---

## Command: complex-sobtop

### Introduction
The `complex-sobtop` command prepares a protein-ligand complex for GROMACS MDS using **sobtop** instead of CGenFF to prepare ligand topology. This is useful when CGenFF is not available or when you prefer to use sobtop for ligand parameter generation.

### Parameters
| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | The directory containing the complex files. |
| `-n` | `--complex-name` | str | None | The filename of the complex within each sub-directory. |
| `-rc` | `--receptor-chain-name` | str | None | The receptor chain name. |
| `-lc` | `--ligand-chain-name` | str | None | The ligand chain name. |
| | `--ff-dir` | str | None | Force field files directory. |
| | `--pdb2gmx-args` | str | `"-ter -ignh"` | Arguments passed to pdb2gmx command. |
| | `--n-term` | str | `"0"` | N-Term type for gmx pdb2gmx. Use `auto` for automatic detection. |
| | `--c-term` | str | `"0"` | C-Term type for gmx pdb2gmx. Use `auto` for automatic detection. |

### Behavior
1. **STEP 0.1**: Center complex.pdb by Open Babel.
2. **STEP 0.2**: Align complex.pdb with xyz axes by LazyDock.
3. **STEP 1**: Extract receptor and ligand from complex.pdb.
4. **STEP 2**: Transfer ligand.pdb to mol2 by Open Babel.
5. **STEP 3**: Fix ligand name in mol2 file.
6. **STEP 4**: Delete UNITY_ATOM_ATTR TRIPOS from mol2 file.
7. **STEP 5**: Prepare ligand topology by sobtop (GAFF atom type).
8. **STEP 6**: Fix ligand name in lig.itp.
9. **STEP 7**: Prepare the forcefield directory.
10. **STEP 8**: Prepare the protein topology using pdb2gmx with terminal handling.
11. **STEP 9**: (Optional) Add receptor position restraints.
12. **STEP 10**: Prepare the complex topology by merging receptor and ligand.

### Notes
- Input complex.pdb should have two chains, one for receptor and one for ligand.
- Complex.pdb should already have hydrogens added.
- Complex.pdb should be aligned with axes to save space during MDS.
- **sobtop** must be configured in `~/.lazydock/config.json` under `named_paths.sobtop_dir`.
- You should cite sobtop if you use it in your work.
- Supports `auto` mode for terminal detection (assumes single protein receptor).

### Example
```bash
# Basic usage
lazydock-cli prepare-gmx complex-sobtop -d /path/to/complex/directory -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir /path/to/ff

# With auto terminal detection
lazydock-cli prepare-gmx complex-sobtop -d /path/to/complex/directory -n complex.pdb --receptor-chain-name A --ligand-chain-name Z --ff-dir /path/to/ff --n-term auto --c-term auto
```
