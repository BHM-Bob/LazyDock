## Command: simple-protein

### Introduction
The `simple-protein` command is used to run a simple GROMACS simulation for a protein. It automates a series of steps including editing the configuration, solvating, ionizing, energy minimization, and molecular dynamics runs.

### Parameters
- `-d`, `--dir`: The directory to store the prepared files.
- `-n`, `--protein-name`: The protein name in each sub-directory, such as `protein.gro`.
- `--auto-box`: A flag to automatically generate a rectangular bounding box using PyMOL.
- `--auto-box-padding`: The distance to padding the box size, default is `1.2`.
- `--ion-mdp`: The energy minimization mdp file.
- `--em-mdp`: The energy minimization mdp file.
- `--nvt-mdp`: The NVT mdp file.
- `--npt-mdp`: The NPT mdp file.
- `--md-mdp`: The production MD mdp file.
- `--editconf-args`: Arguments passed to the `editconf` command, default is `"-c -d 1.2 -bt dodecahedron"`.
- `--solvate-args`: Arguments passed to the `solvate` command, default is `"-cs spc216.gro"`.
- `--genion-args`: Arguments passed to the `genion` command, default is `"-pname NA -nname CL -neutral"`.
- `--mdrun-args`: Arguments passed to the `mdrun` command, default is `"-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu"`.

### Behavior
1. Edit the protein configuration to create a new box.
2. Solvate the protein within the box.
3. Add ions to the solvated system.
4. Perform energy minimization.
5. Run equilibration simulations (NVT and NPT).
6. Conduct production molecular dynamics simulation.

### Notes
- The directory argument must be a valid directory.
- The protein file name must be specified to identify the correct files.
- MDP files are required for each step and must be provided.

### Example
```bash
lazydock-cli run-gmx simple-protein -d /path/to/protein/directory -n protein.gro
```

## Command: simple-ligand

### Introduction
The `simple-ligand` command is similar to `simple-protein` but is designed for a ligand. It follows the same steps but is tailored for ligand-specific simulations.

### Parameters
- `-d`, `--dir`: The directory to store the prepared files.
- `-n`, `--protein-name`: The protein name in each sub-directory, such as `protein.gro`.
- `--auto-box`: A flag to automatically generate a rectangular bounding box using PyMOL.
- `--auto-box-padding`: The distance to padding the box size, default is `1.2`.
- `--ion-mdp`: The energy minimization mdp file.
- `--em-mdp`: The energy minimization mdp file.
- `--nvt-mdp`: The NVT mdp file.
- `--npt-mdp`: The NPT mdp file.
- `--md-mdp`: The production MD mdp file.
- `--editconf-args`: Arguments passed to the `editconf` command, default is `"-c -d 1.2 -bt dodecahedron"`.
- `--solvate-args`: Arguments passed to the `solvate` command, default is `"-cs spc216.gro"`.
- `--genion-args`: Arguments passed to the `genion` command, default is `"-pname NA -nname CL -neutral"`.
- `--mdrun-args`: Arguments passed to the `mdrun` command, default is `"-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu"`.

### Behavior
1. Edit the ligand configuration to create a new box.
2. Solvate the ligand within the box.
3. Add ions to the solvated system.
4. Perform energy minimization.
5. Run equilibration simulations (NVT and NPT).
6. Conduct production molecular dynamics simulation.

### Notes
- The directory argument must be a valid directory.
- The ligand file name must be specified to identify the correct files.
- MDP files are required for each step and must be provided.

### Example
```bash
lazydock-cli run-gmx simple-ligand -d /path/to/ligand/directory -n ligand.gro
```

## Command: simple-complex

### Introduction
The `simple-complex` command is used to run a simple GROMACS simulation for a protein-ligand complex. It extends the `simple-protein` command with additional steps for handling the ligand within the complex.

### Parameters
- `-d`, `--dir`: The directory to store the prepared files.
- `-n`, `--protein-name`: The protein name in each sub-directory, such as `protein.gro`.
- `--ln`, `--ligand-name`: The ligand name in each sub-directory, such as `lig.gro`.
- `--lig-posres`: The ligand position restraint symbol, default is `POSRES`.
- `--tc-groups`: The tc-grps to select, default is `"1 | 13"`.
- `--auto-box`: A flag to automatically generate a rectangular bounding box using PyMOL.
- `--auto-box-padding`: The distance to padding the box size, default is `1.2`.
- `--ion-mdp`: The energy minimization mdp file.
- `--em-mdp`: The energy minimization mdp file.
- `--nvt-mdp`: The NVT mdp file.
- `--npt-mdp`: The NPT mdp file.
- `--md-mdp`: The production MD mdp file.
- `--editconf-args`: Arguments passed to the `editconf` command, default is `"-c -d 1.2 -bt dodecahedron"`.
- `--solvate-args`: Arguments passed to the `solvate` command, default is `"-cs spc216.gro"`.
- `--genion-args`: Arguments passed to the `genion` command, default is `"-pname NA -nname CL -neutral"`.
- `--mdrun-args`: Arguments passed to the `mdrun` command, default is `"-v -ntomp 4 -update gpu -nb gpu -pme gpu -bonded gpu -pmefft gpu"`.

### Behavior
1. Edit the complex configuration to create a new box.
2. Solvate the complex within the box.
3. Add ions to the solvated system.
4. Perform energy minimization.
5. Run equilibration simulations (NVT and NPT).
6. Conduct production molecular dynamics simulation with additional steps for ligand restraints.

### Notes
- The directory argument must be a valid directory.
- The protein and ligand file names must be specified to identify the correct files.
- MDP files are required for each step and must be provided.

### Example
```bash
lazydock-cli run-gmx simple-complex -d /path/to/complex/directory -n protein.gro -ln lig.gro
```
