## Command: simple-protein

### Introduction
The `simple-protein` command automates a complete GROMACS simulation workflow for proteins. It performs a series of steps including box creation, solvation, ionization, energy minimization, equilibration, and production molecular dynamics.

### Parameters
| Short | Long                | Type      | Default | Description |
|-------|---------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`       | str (+)  | `['.']` | Directories containing protein files to process. |
| `-n`  | `--protein-name`    | str      | Required | Protein filename in each sub-directory (e.g., `protein.gro`). |
| `-st` | `--start-time`      | str      | `None` | Start time for execution (format: "2025-01-15 17:30:00"). |
|       | `--auto-box`        | flag     | `False` | Automatically generate rectangular bounding box using PyMOL. |
|       | `--auto-box-padding`| float    | `1.2` | Padding distance for automatic box generation. |
|       | `--ion-mdp`         | str      | Required | Ion addition MDP file. |
|       | `--em-mdp`          | str      | Required | Energy minimization MDP file. |
|       | `--nvt-mdp`         | str      | Required | NVT equilibration MDP file. |
|       | `--npt-mdp`         | str      | Required | NPT equilibration MDP file. |
|       | `--md-mdp`          | str      | Required | Production MD MDP file. |
|       | `--editconf-args`   | str      | `"-c -d 1.2 -bt dodecahedron"` | Arguments for editconf command. |
|       | `--solvate-args`    | str      | `"-cs spc216.gro"` | Arguments for solvate command. |
|       | `--genion-args`     | str      | `"-pname NA -nname CL -neutral"` | Arguments for genion command. |
|       | `--em-args`         | str      | `""` | Arguments for mdrun during energy minimization. |
|       | `--mdrun-args`      | str      | `"-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu"` | Arguments for production mdrun. |
|       | `--genion-groups`   | str      | `"SOL"` | Group for ion replacement. |
|       | `--potential-groups`| str      | `"11 0"` | Groups to plot potential energy. |
|       | `--temperature-groups`| str     | `"16 0"` | Groups to plot temperature. |
|       | `--pressure-groups` | str      | `"17 0"` | Groups to plot pressure. |
|       | `--density-groups`  | str      | `"23 0"` | Groups to plot density. |
|       | `--maxwarn`         | int      | `0` | Maximum warnings for grompp commands. |

### Behavior
1. **Box Creation**: Uses `editconf` to create a simulation box around the protein
2. **Solvation**: Adds water molecules using `solvate`
3. **Ionization**: Replaces water molecules with ions using `genion`
4. **Energy Minimization**: Minimizes system energy using steepest descent
5. **NVT Equilibration**: Equilibrates system at constant volume and temperature
6. **NPT Equilibration**: Equilibrates system at constant pressure and temperature
7. **Production MD**: Runs the production molecular dynamics simulation

### Notes
- The command automatically generates plots for potential, temperature, pressure, and density
- If `md.tpr` already exists in the directory, the simulation is skipped
- MDP files are required for each simulation step
- Supports GPU acceleration through mdrun arguments

### Example
```bash
lazydock-cli run-gmx simple-protein -d /path/to/protein/directory -n protein.gro --ion-mdp ions.mdp --em-mdp em.mdp --nvt-mdp nvt.mdp --npt-mdp npt.mdp --md-mdp md.mdp
```

---

## Command: simple-ligand

### Introduction
The `simple-ligand` command automates a complete GROMACS simulation workflow for ligands. It follows the same steps as `simple-protein` but is optimized for small molecule simulations with appropriate default parameters.

### Parameters
| Short | Long                | Type      | Default | Description |
|-------|---------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`       | str (+)  | `['.']` | Directories containing ligand files to process. |
| `-n`  | `--protein-name`    | str      | Required | Ligand filename in each sub-directory (e.g., `ligand.gro`). |
| `-ln` | `--ligand-name`     | str      | `'lig.gro'` | Ligand filename (used for compatibility with complex workflow). |
| `-st` | `--start-time`      | str      | `None` | Start time for execution (format: "2025-01-15 17:30:00"). |
|       | `--auto-box`        | flag     | `False` | Automatically generate rectangular bounding box using PyMOL. |
|       | `--auto-box-padding`| float    | `1.2` | Padding distance for automatic box generation. |
|       | `--ion-mdp`         | str      | Required | Ion addition MDP file. |
|       | `--em-mdp`          | str      | Required | Energy minimization MDP file. |
|       | `--nvt-mdp`         | str      | Required | NVT equilibration MDP file. |
|       | `--npt-mdp`         | str      | Required | NPT equilibration MDP file. |
|       | `--md-mdp`          | str      | Required | Production MD MDP file. |
|       | `--editconf-args`   | str      | `"-c -d 1.2 -bt dodecahedron"` | Arguments for editconf command. |
|       | `--solvate-args`    | str      | `"-cs spc216.gro"` | Arguments for solvate command. |
|       | `--genion-args`     | str      | `"-pname NA -nname CL -neutral"` | Arguments for genion command. |
|       | `--em-args`         | str      | `""` | Arguments for mdrun during energy minimization. |
|       | `--mdrun-args`      | str      | `"-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu"` | Arguments for production mdrun. |
|       | `--genion-groups`   | str      | `"SOL"` | Group for ion replacement. |
|       | `--potential-groups`| str      | `"10 0"` | Groups to plot potential energy (optimized for ligands). |
|       | `--temperature-groups`| str     | `"15 0"` | Groups to plot temperature (optimized for ligands). |
|       | `--pressure-groups` | str      | `"16 0"` | Groups to plot pressure (optimized for ligands). |
|       | `--density-groups`  | str      | `"22 0"` | Groups to plot density (optimized for ligands). |
|       | `--tc-groups`       | str      | `"2"` | Temperature coupling groups for ligand simulations. |
|       | `--maxwarn`         | int      | `0` | Maximum warnings for grompp commands. |

### Behavior
1. **Box Creation**: Uses `editconf` to create a simulation box around the ligand
2. **Solvation**: Adds water molecules using `solvate`
3. **Ionization**: Replaces water molecules with ions using `genion`
4. **Energy Minimization**: Minimizes system energy using steepest descent
5. **NVT Equilibration**: Equilibrates system at constant volume and temperature
6. **NPT Equilibration**: Equilibrates system at constant pressure and temperature
7. **Production MD**: Runs the production molecular dynamics simulation

### Notes
- Uses ligand-optimized default parameters for energy analysis groups
- Temperature coupling groups are set to "2" for ligand simulations
- If `md.tpr` already exists in the directory, the simulation is skipped
- MDP files are required for each simulation step
- Supports GPU acceleration through mdrun arguments

### Example
```bash
lazydock-cli run-gmx simple-ligand -d /path/to/ligand/directory -n ligand.gro --ion-mdp ions.mdp --em-mdp em.mdp --nvt-mdp nvt.mdp --npt-mdp npt.mdp --md-mdp md.mdp
```

---

## Command: simple-complex

### Introduction
The `simple-complex` command automates a complete GROMACS simulation workflow for protein-ligand complexes. It extends the protein workflow with additional steps for handling ligand position restraints and proper temperature coupling groups.

### Parameters
| Short | Long                | Type      | Default | Description |
|-------|---------------------|-----------|---------|-------------|
| `-d`  | `--batch-dir`       | str (+)  | `['.']` | Directories containing complex files to process. |
| `-n`  | `--protein-name`    | str      | Required | Complex filename in each sub-directory (e.g., `complex.gro`). |
| `-ln` | `--ligand-name`     | str      | `'lig.gro'` | Ligand filename in each sub-directory. |
| `-st` | `--start-time`      | str      | `None` | Start time for execution (format: "2025-01-15 17:30:00"). |
|       | `--auto-box`        | flag     | `False` | Automatically generate rectangular bounding box using PyMOL. |
|       | `--auto-box-padding`| float    | `1.2` | Padding distance for automatic box generation. |
|       | `--ion-mdp`         | str      | Required | Ion addition MDP file. |
|       | `--em-mdp`          | str      | Required | Energy minimization MDP file. |
|       | `--nvt-mdp`         | str      | Required | NVT equilibration MDP file. |
|       | `--npt-mdp`         | str      | Required | NPT equilibration MDP file. |
|       | `--md-mdp`          | str      | Required | Production MD MDP file. |
|       | `--editconf-args`   | str      | `"-c -d 1.2 -bt dodecahedron"` | Arguments for editconf command. |
|       | `--solvate-args`    | str      | `"-cs spc216.gro"` | Arguments for solvate command. |
|       | `--genion-args`     | str      | `"-pname NA -nname CL -neutral"` | Arguments for genion command. |
|       | `--em-args`         | str      | `""` | Arguments for mdrun during energy minimization. |
|       | `--mdrun-args`      | str      | `"-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu"` | Arguments for production mdrun. |
|       | `--genion-groups`   | str      | `"SOL"` | Group for ion replacement. |
|       | `--potential-groups`| str      | `"11 0"` | Groups to plot potential energy. |
|       | `--temperature-groups`| str     | `"16 0"` | Groups to plot temperature. |
|       | `--pressure-groups` | str      | `"17 0"` | Groups to plot pressure. |
|       | `--density-groups`  | str      | `"23 0"` | Groups to plot density. |
|       | `--lig-posres`      | str      | `'POSRES'` | Symbol for ligand position restraints. |
|       | `--tc-groups`       | str      | `'auto'` | Temperature coupling groups (auto-detects Protein and LIG). |
|       | `--maxwarn`         | int      | `0` | Maximum warnings for grompp commands. |

### Behavior
1. **Box Creation**: Uses `editconf` to create a simulation box around the complex
2. **Solvation**: Adds water molecules using `solvate`
3. **Ionization**: Replaces water molecules with ions using `genion`
4. **Energy Minimization**: Minimizes system energy using steepest descent
5. **Ligand Restraints**: Creates position restraints for the ligand:
   - Generates index file for ligand heavy atoms
   - Creates position restraint file with force constants of 1000 kJ/mol/nm²
   - Inserts restraint information into topology file
6. **Temperature Coupling Groups**: Creates custom temperature coupling groups:
   - Automatically detects Protein and LIG groups if set to 'auto'
   - Creates index file for proper temperature coupling
7. **NVT Equilibration**: Equilibrates system at constant volume and temperature
8. **NPT Equilibration**: Equilibrates system at constant pressure and temperature
9. **Production MD**: Runs the production molecular dynamics simulation

### Notes
- Automatically creates ligand position restraints to prevent unrealistic movements
- Temperature coupling groups are set up to allow different temperatures for protein and ligand
- If `md.tpr` already exists in the directory, the simulation is skipped
- MDP files are required for each simulation step
- Supports GPU acceleration through mdrun arguments

### Example
```bash
lazydock-cli run-gmx simple-complex -d /path/to/complex/directory -n complex.gro -ln ligand.gro --ion-mdp ions.mdp --em-mdp em.mdp --nvt-mdp nvt.mdp --npt-mdp npt.mdp --md-mdp md.mdp
```
