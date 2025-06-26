### Command: trjconv  
#### Introduction  
Processes trajectory files using GROMACS `trjconv` to center the system and handle periodic boundary conditions (PBC). This command automates trajectory alignment for structural analysis across multiple directories.  

#### Parameters  
| Short | Long             | Type      | Default | Description |
|-------|------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`    | str (+)  | `['.']` | Directories containing trajectory files to process. |
| `-n`  | `--main-name`    | str      | `'md.tpr'` | Filename pattern for topology/input files (e.g., `md.tpr`). |
| `-g`  | `--groups`       | str (+)  | `['1','0']` | Atom groups used by `trjconv` for alignment. |
|       | `-pbc`           | str      | `'mol'` | PBC treatment method (`mol`, `atom`, `res`, etc.). |
|       | `-ur`            | str      | `'compact'` | Unit-cell representation (`rect`, `tric`, `compact`). |
| `-nw` | `--n-workers`    | int      | `1` | Number of parallel workers for processing. |
| `-F`  | `--force`        | flag     | `False` | Reprocess trajectories even if output exists. |
| `-D`  | `--delete`       | flag     | `False` | Delete existing output files before reprocessing. |

#### Behavior  
1. **Locates Files**: Searches for `main-name.tpr` in each directory.  
2. **Skips Existing**: Ignores directories where `*_center.xtc` exists (unless `--force`/`--delete` is used).  
3. **Runs `trjconv`**: Executes with parameters `-center`, `-pbc`, and `-ur`. Users are prompted interactively for group selections.  
4. **Parallel Processing**: Runs tasks concurrently using `n-workers`.  

#### Notes  
- Output file: `<main-name>_center.xtc`.  
- Requires GROMACS.  
- Works with `.tpr` topology files and `.xtc` trajectories.  

#### Example  
```bash
lazydock-cli trjconv -d /path/to/trajectories -n md -g 1 0 -pbc mol -ur compact -nw 4 -F
```

---

### Command: make_ndx  
#### Introduction  
Creates or updates index files for trajectory analysis using GROMACS `make_ndx`. Specifies atom groups (e.g., protein/ligand) for focused analysis.  

#### Parameters  
| Short | Long             | Type      | Default | Description |
|-------|------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`    | str (+)  | `['.']` | Directories containing topology files. |
| `-f`  | `--main-name`    | str      | `'md.tpr'` | Filename pattern for topology files. |
| `-g`  | `--groups`       | str (+)  | `['1','0']` | Groups to include in the index file. |
| `-o`  | `--output`       | str      | `'ana_index.ndx'` | Name for the output index file. |
| `-F`  | `--force`        | flag     | `False` | Overwrite existing index files. |
| `-D`  | `--delete`       | flag     | `False` | Delete existing index files first. |

#### Behavior  
1. **Locates Files**: Finds `main-name.tpr` in each directory.  
2. **Skipping**: Skips processing if output exists (unless `--force`/`--delete` used).  
3. **Runs `make_ndx`**: Adds specified groups to the index file. Automatically exits after group selection.  

#### Notes  
- Requires topology files (`.tpr`).  
- Output: `<output>.ndx`.  

#### Example  
```bash
lazydock-cli make_ndx -d /path/to/data -f md -g "Protein" "Ligand" -o complex.ndx
```

---

### Command: simple  
#### Introduction  
Performs standard MD analyses (RMSD, RMSF, Rg, Hbonds, SASA, PCA/DSSP) and plots results. Automates common GROMACS tools and visualization.  

#### Parameters  
| Short | Long                | Type      | Default               | Description |
|-------|---------------------|-----------|-----------------------|-------------|
| `-d`  | `--batch-dir`       | str (+)  | `['.']` | Directories to analyze. |
| `-n`  | `--main-name`       | str      | `'md.tpr'` | Topology filename pattern. |
| `-ndx`| `--index`           | str      | `None` | Input index file (e.g., `ana_index.ndx`). |
|       | `--methods`         | str (+)  | `[all]` | Analyses to run: `rms`, `rmsf`, `gyrate`, `hbond`, `sasa`, `covar`, `dssp`, `FEL`, `PDF`. |
|       | `--dit-style`       | str      | `None` | Path to `.mplstyle` file for plot customization. |
| `-rg` | `--rms-group`       | str      | `'4'` | Atom group for RMSD/RMSF/Rg calculations. |
| `-hg` | `--hbond-group`     | int (+)  | `[1,1]` | Groups for Hbond analysis. |
| `-sg` | `--sasa-group`      | str      | `'4'` | Group for SASA calculation. |
| `-eg` | `--eigenval-group`  | str      | `'4'` | Group for PCA. |
| `-dg` | `--dssp-group`      | str      | `'1'` | Group for DSSP secondary structure. |
|       | `--dssp-num`        | flag     | `False` | Calculate DSSP residue types. |
|       | `--dssp-clear`      | flag     | `False` | Clear DSSP selections. |
| `-F`  | `--force`           | flag     | `False` | Reprocess existing data. |
| `-D`  | `--delete`          | flag     | `False` | Delete existing results. |

#### Behavior  
1. **Runs Analysis**: Executes tools like `rms`, `rmsf`, `hbond`, etc., per `methods` list.  
2. **Generates Plots**: Uses `dit` for visualization. Plots include RMSD, FEL, and PDFs.  
3. **Handles Skips**: Ignores tasks where output exists unless `-F`/`-D` is specified.  

#### Notes  
- Requires `*_center.xtc` (preprocessed by `trjconv`).  
- `FEL`: Free-energy landscape via MD-DaVis.  
- `PDF`: Probability density function for RMSD/Rg.  

#### Example  
```bash
lazydock-cli simple -d /md_runs --methods rms gyrate FEL --hbond-group 1 1 -F
```

---

### Command: mmpbsa  
#### Introduction  
Calculates binding free energies using the MMPBSA method with `gmx_MMPBSA`. Defines receptor/ligand groups and runs calculations.  

#### Parameters  
| Short | Long                | Type      | Default  | Description |
|-------|---------------------|-----------|----------|-------------|
| `-d`  | `--batch-dir`       | str (+)  | `['.']`  | Directories containing simulation files. |
| `-i`  | `--input`           | str      | Required | MMPBSA config file (e.g., `mmpbsa.in`). |
| `-o`  | `--output`          | str      | `'MMPBSA_FINAL_RESULTS'` | Output prefix for results. |
| `-np` | `--np`              | int      | Required | Number of MPI processes for `gmx_MMPBSA`. |
| `-top`| `--top-name`        | str      | `'md.tpr'` | Topology filename pattern. |
| `-traj`| `--traj-name`      | str      | `'md_center.xtc'` | Trajectory filename pattern. |
|       | `--receptor-chain-name` | str | Required | Receptor chain ID (e.g., `"A"`). |
|       | `--ligand-chain-name` | str | Required | Ligand chain ID (e.g., `"LIG"`). |
| `-F`  | `--force`            | flag     | `False` | Reprocess even if output exists. |

#### Behavior  
1. **Defines Groups**: Automatically extracts receptor/ligand atom ranges from topology.  
2. **Creates Index**: Generates `mmpbsa.ndx` with receptor/ligand groups.  
3. **Runs MMPBSA**: Executes `gmx_MMPBSA` with MPI (`-np`) on trajectories.  

#### Notes  
- Requires `gmx_MMPBSA` installation and `*_center.xtc`.  
- Output: `<output>.dat` (energies) and `<output>.csv`.  

#### Example  
```bash
lazydock-cli mmpbsa -d /binding_sims -i mmpbsa.in -np 8 --receptor-chain-name A --ligand-chain-name LIG
```

---

### Command: interaction  
#### Introduction  
Identifies non-covalent interactions (e.g., H-bonds, hydrophobic) between receptor and ligand across frames. Supports PyMOL and PLIP.  

#### Parameters  
| Short | Long                | Type      | Default | Description |
|-------|---------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`       | str (+)  | `['.']` | Directories to analyze. |
| `-gro`| `--gro-name`        | str      | `'md.gro'` | Structure file (.gro). |
|       | `--alter-receptor-chain` | str | `None` | Override receptor chain ID. |
|       | `--alter-ligand-chain` | str | `None` | Override ligand chain ID. |
|       | `--alter-ligand-res` | str | `None` | Override ligand residue name. |
|       | `--alter-ligand-atm` | str | `None` | Override ligand atom type. |
|       | `--method`         | str      | `'pymol'` | Tool: `pymol` or `plip`. |
|       | `--mode`           | str      | `'all'` | Interaction types (comma-separated). |
|       | `--cutoff`         | float    | `4.0` | Distance threshold (Å). |
|       | `--hydrogen-atom-only` | flag | `False` | Only analyze polar contacts. |
|       | `--output-style`   | str      | `'receptor'` | Output table format. |
| `-nw` | `--n-workers`       | int      | `4` | Number of parallel processes. |
| `-b`  | `--begin-frame`     | int      | `1` | First frame to analyze. |
| `-e`  | `--end-frame`       | int      | `None` | Last frame to analyze. |
|       | `--plot-time-unit`  | int      | `100` | Frames per heatmap time unit. |

#### Behavior  
1. **Maps Interactions**: Per-frame detection using PyMOL/PLIP and stores results.  
2. **Filters/Plots**: Generates residue-level heatmaps showing interaction frequencies.  
3. **Parallelizes**: Frame processing distributed across workers.  

#### Notes  
- Output: CSV tables and heatmap plots (`*_interactions.png`).  
- Use `--alter-*` if topology chain IDs mismatch expectations.  

#### Example  
```bash
lazydock-cli interaction -d /complexes --method plip --mode all --cutoff 5.0 --nw 8
```

---

### Command: rrcs  
#### Introduction  
Computes Residue-Residue Contact Scores (RRCS) to quantify dynamic interactions within chains across frames. Generates heatmaps and line plots.  

#### Parameters  
| Short | Long              | Type      | Default | Description |
|-------|-------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`     | str (+)  | `['.']` | Directories to analyze. |
| `-top`| `--top-name`      | str      | `'md.tpr'` | Topology filename pattern. |
| `-gro`| `--gro-name`      | str      | `'md.gro'` | Structure file (.gro). |
| `-traj`| `--traj-name`    | str      | `'md_center.xtc'` | Trajectory filename pattern. |
| `-c`  | `--chains`        | str (+)  | `None` | Chains to analyze (e.g., `A`, `B`). |
| `-np` | `--n-workers`     | int      | `4` | Parallel worker count. |
| `-b`  | `--begin-frame`   | int      | `0` | First frame to process. |
| `-e`  | `--end-frame`     | int      | `None` | Last frame to process. |
|       | `--backend`       | str      | `'numpy'` | Compute backend: `numpy`, `torch`, or `cuda`. |

#### Behavior  
1. **Calculates RRCS**: Per-frame residue pair scoring.  
2. **Aggregates & Plots**: Generates three plots:  
   - Average heatmap of RRCS.  
   - Diagonal RRCS values vs. frame.  
   - Line plots for top-scoring residue pairs.  
3. **Fast Backends**: GPU acceleration via PyTorch (`--backend cuda`).  

#### Notes  
- Output: `.npz` file with scores and PNG plots.  

#### Example  
```bash
lazydock-cli rrcs -d /dynamics --chains A B --nw 4 --backend cuda
```

---

### Command: porcupine  
#### Introduction  
Generates PyMOL "porcupine plots" using `modevectors` to visualize structural dynamics. Shows displacement vectors between start/end states.  

#### Parameters  
| Short | Long              | Type      | Default | Description |
|-------|-------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`     | str (+)  | `['.']` | Directories to analyze. |
| `-top`| `--top-name`      | str      | `'md.tpr'` | Topology filename pattern. |
| `-traj`| `--traj-name`    | str      | `'md_center.xtc'` | Trajectory filename pattern. |
| `-b`  | `--begin-frame`   | int      | `0` | First frame. |
| `-e`  | `--end-frame`     | int      | `None` | Last frame. |
|       | `-D`             | flag      | `False` | Delete existing output. |

#### Behavior  
1. **Loads Trajectory**: Splits into start/end states.  
2. **Calculates Vectors**: Runs `modevectors` to generate displacements.  
3. **Saves Session**: Exports `.pse` file for PyMOL.  

#### Notes  
- Requires PyMOL and MDAnalysis.  
- Output: `<traj>_porcupine.pse`.  

#### Example  
```bash
lazydock-cli porcupine -d /motions -b 0 -e 1000
```