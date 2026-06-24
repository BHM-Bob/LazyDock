<!--
 * @Date: 2025-06-15
 * @LastEditors: BHM-Bob
 * @Description: run_amber script documentation
-->

# run_amber

The `run_amber` script automates AMBER molecular dynamics simulations using the "Golden Standard" 6-step equilibration protocol. It supports running simulations for proteins, ligands, and protein-ligand complexes using pmemd (CPU or GPU).

## Batch Directory Design

The `--batch-dir` parameter is designed to point to a **single parent directory** that contains **multiple task subdirectories**. Each subdirectory should contain one simulation system (prmtop + inpcrd pair).

### Directory Structure Example

```
/path/to/batch_dir/              <-- batch-dir
├── system_1/                    <-- task subdirectory
│   ├── complex.prmtop
│   └── complex.inpcrd
├── system_2/                    <-- task subdirectory
│   ├── complex.prmtop
│   └── complex.inpcrd
└── system_3/                    <-- task subdirectory
    ├── complex.prmtop
    └── complex.inpcrd
```

### Key Design Principles

1. **Single Batch Directory**: The `-d` parameter should point to one parent directory
2. **Multiple Task Subdirectories**: Each subdirectory is treated as an independent simulation task
3. **Automatic File Pairing**: The script uses `check_prmtop_inpcrd()` to automatically pair prmtop and inpcrd files based on their directory location
4. **Parallel Task Discovery**: All tasks are discovered upfront, validated, then executed sequentially

### Usage Example

```bash
# Run simulations for all systems in /data/simulations/
# Each subdirectory contains one system
lazydock-cli run-amber simple-md -d /data/simulations

# The script will:
# 1. Find all *.prmtop files in /data/simulations and its subdirectories
# 2. Pair each prmtop with its corresponding *.inpcrd in the same directory
# 3. Validate that each directory has exactly one prmtop and one inpcrd
# 4. Run the 6-step pipeline for each system sequentially
```

## Pipeline Overview

The script implements the recommended AMBER equilibration workflow consisting of 6 sequential steps:

```
Input: prmtop + inpcrd
    ↓
Step 1: Minimization with restraints (solute heavy atoms)
    ↓
Step 2: Minimization without restraints
    ↓
Step 3: Heating (0K → 300K) with restraints (protein backbone)
    ↓
Step 4: NPT equilibration with weak restraints (protein backbone)
    ↓
Step 5: NPT equilibration without restraints
    ↓
Step 6: Production MD
    ↓
Output: Production trajectory (md.nc) and restart files
```

Each step uses the output of the previous step as input, creating a continuous simulation workflow.

## Commands

### simple-md

Runs the complete 6-step AMBER MD simulation pipeline. This is the main command for running production simulations.

#### Parameters

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--batch-dir` | str (+) | `['.']` | Directories containing simulation files |
| `-n` | `--prmtop-name` | str | `*.prmtop` | Topology file name pattern |
| `-c` | `--inpcrd-name` | str | `*.inpcrd` | Coordinate file name pattern |
| | `--min1-maxcyc` | int | `5000` | Max cycles for minimization step 1 |
| | `--min1-ncyc` | int | `2500` | Steepest descent steps for min1 |
| | `--min1-restraint-wt` | float | `10.0` | Restraint weight for min1 (kcal/mol/Å²) |
| | `--min1-restraint-mask` | str | `!@H=` | Restraint mask for min1 |
| | `--min2-maxcyc` | int | `5000` | Max cycles for minimization step 2 |
| | `--min2-ncyc` | int | `2500` | Steepest descent steps for min2 |
| | `--heat-nstlim` | int | `50000` | Steps for heating (100ps @ 2fs) |
| | `--heat-temp0` | float | `300.0` | Target temperature for heating |
| | `--heat-restraint-wt` | float | `5.0` | Restraint weight for heating |
| | `--heat-restraint-mask` | str | `@CA,C,N` | Restraint mask for heating |
| | `--npt1-nstlim` | int | `50000` | Steps for NPT step 1 (100ps @ 2fs) |
| | `--npt1-restraint-wt` | float | `1.0` | Restraint weight for NPT step 1 |
| | `--npt1-restraint-mask` | str | `@CA,C,N` | Restraint mask for NPT step 1 |
| | `--npt2-nstlim` | int | `100000` | Steps for NPT step 2 (200ps @ 2fs) |
| | `--md-nstlim` | int | `50000000` | Steps for production MD (100ns @ 2fs) |
| | `--md-temp0` | float | `300.0` | Target temperature for production |
| | `--md-ntpr` | int | `5000` | Energy output frequency |
| | `--md-ntwx` | int | `5000` | Coordinate output frequency |
| | `--md-ntwr` | int | `50000` | Restart file output frequency |
| | `--cut` | float | `10.0` | Nonbonded cutoff (Å) |
| | `--dt` | float | `0.002` | Time step (ps) |
| | `--barostat` | int | `2` | Barostat: 1=Berendsen, 2=Monte Carlo |
| | `--start-step` | int | `1` | Start from step (1-6) |
| | `--end-step` | int | `6` | End at step (1-6) |

#### Pipeline Steps

##### Step 1: Minimization with Restraints

**Purpose**: Relax solvent and ions while keeping solute fixed

**Key Parameters**:
- `imin = 1`: Energy minimization
- `ntb = 1`: Constant volume
- `ntr = 1`: Enable restraints
- `restraint_wt = 10.0`: Strong restraints on heavy atoms
- `restraintmask = '!@H='`: Restrain all non-hydrogen atoms

**Input**: Initial coordinates from `prepare_amber`
**Output**: `step01_min1.rst` (restart file)

##### Step 2: Minimization without Restraints

**Purpose**: Relax entire system

**Key Parameters**:
- `imin = 1`: Energy minimization
- `ntb = 1`: Constant volume
- `ntr = 0`: No restraints

**Input**: `step01_min1.rst`
**Output**: `step02_min2.rst`

##### Step 3: Heating (0K → 300K)

**Purpose**: Gradually heat system to target temperature with backbone restraints

**Key Parameters**:
- `imin = 0`: Molecular dynamics
- `ntb = 1`: Constant volume (NVT)
- `tempi = 0.0`: Initial temperature
- `temp0 = 300.0`: Target temperature
- `ntt = 3`: Langevin thermostat
- `ntr = 1`: Enable restraints
- `restraintmask = '@CA,C,N'`: Restrain protein backbone

**Input**: `step02_min2.rst`
**Output**: `step03_heat.rst`, `step03_heat.nc` (trajectory)

##### Step 4: NPT Equilibration with Weak Restraints

**Purpose**: Equilibrate density with weak backbone restraints

**Key Parameters**:
- `ntb = 2`: Constant pressure (NPT)
- `ntp = 1`: Isotropic pressure scaling
- `pres0 = 1.0`: 1 atm pressure
- `barostat = 2`: Monte Carlo barostat (recommended for GPU)
- `restraint_wt = 1.0`: Weak restraints

**Input**: `step03_heat.rst`
**Output**: `step04_npt1.rst`, `step04_npt1.nc`

##### Step 5: NPT Equilibration without Restraints

**Purpose**: Final equilibration without restraints

**Key Parameters**:
- `ntr = 0`: No restraints
- Longer simulation (200ps) for complete equilibration

**Input**: `step04_npt1.rst`
**Output**: `step05_npt2.rst`, `step05_npt2.nc`

##### Step 6: Production MD

**Purpose**: Production simulation for data collection

**Key Parameters**:
- `nstlim = 50000000`: 100ns at 2fs timestep
- `ntpr = 5000`: Energy output every 10ps
- `ntwx = 5000`: Coordinates every 10ps
- `ntwr = 50000`: Restart file every 100ps

**Input**: `step05_npt2.rst`
**Output**: `step06_md.rst`, `step06_md.nc` (production trajectory)

#### Progress Monitoring

The script provides real-time progress monitoring with tqdm progress bars:

```
heat:  45%|████▌     | 22500/50000 [02:15<02:45, 166.67steps/s] T=285.3K | 45.2ns/day
```

Information displayed:
- Current step and progress percentage
- Elapsed and estimated remaining time
- Current temperature (for MD steps)
- Current pressure (for NPT steps)
- Performance in ns/day

#### Checkpoint and Restart

The script automatically:
- Checks for existing output files to avoid re-running completed steps
- Generates restart files (`.rst`) at each step for recovery
- Can resume from any step using `--start-step`

#### Example

```bash
# Complete simulation from step 1 to 6
lazydock-cli run-amber simple-md -d /path/to/simulation \
    -n complex.prmtop -c complex.inpcrd

# Restart from step 4 (NPT equilibration)
lazydock-cli run-amber simple-md -d /path/to/simulation \
    -n complex.prmtop -c complex.inpcrd \
    --start-step 4

# Run only minimization steps
lazydock-cli run-amber simple-md -d /path/to/simulation \
    -n complex.prmtop -c complex.inpcrd \
    --start-step 1 --end-step 2

# Custom production length (50ns instead of 100ns)
lazydock-cli run-amber simple-md -d /path/to/simulation \
    -n complex.prmtop -c complex.inpcrd \
    --md-nstlim 25000000

# Custom restraints for ligand simulations
lazydock-cli run-amber simple-md -d /path/to/ligand \
    -n ligand.prmtop -c ligand.inpcrd \
    --heat-restraint-mask '!@H=' \
    --npt1-restraint-mask '!@H='
```

---

## Simulation Protocol Details

### Temperature Control

Amber uses Langevin dynamics (`ntt = 3`) for temperature control:
- `gamma_ln = 2.0`: Collision frequency (ps⁻¹)
- No need for temperature coupling groups (unlike GROMACS)
- Same thermostat applied to entire system

### Pressure Control

Two barostat options available:

| Barostat | Value | Description | Recommended For |
|----------|-------|-------------|---------------|
| Berendsen | 1 | Weak coupling | CPU simulations |
| Monte Carlo | 2 | MC barostat | GPU simulations (default) |

### Constraints

- `ntc = 2`: SHAKE algorithm for hydrogen bonds
- `ntf = 2`: Bond interactions involving H omitted
- Allows 2fs timestep

### Restraint Masks

AMBER mask syntax for restraints:

| Mask | Description |
|------|-------------|
| `!@H=` | All heavy atoms (non-hydrogen) |
| `@CA,C,N` | Protein backbone atoms |
| `:LIG` | Ligand residue |
| `@*` | All atoms |

### Output Files

Each step generates:

| Suffix | Description |
|--------|-------------|
| `.out` | Log file with energies and properties |
| `.rst` | Restart file (coordinates and velocities) |
| `.nc` | NetCDF trajectory (MD steps only) |
| `.mdinfo` | Progress information file |

## GPU Acceleration

The script automatically detects and uses GPU-accelerated pmemd:

1. `pmemd.cuda_SPFP` - Single precision (fastest)
2. `pmemd.cuda_DPFP` - Double precision
3. `pmemd.cuda` - Default CUDA version
4. `pmemd` - CPU version (fallback)

### Performance Tips

- Use Monte Carlo barostat (`--barostat 2`) for GPU
- Typical performance: 50-150 ns/day for ~50k atoms on modern GPUs
- Monitor `ns/day` in progress bar

## Continuing Simulations

If a simulation is interrupted, it can be continued:

```bash
# Continue from last restart file
pmemd.cuda -O -i md.in -o md.out -p complex.prmtop \
    -c step06_md.rst -r md_continue.rst -x md_continue.nc

# Or use the script to restart from a specific step
lazydock-cli run-amber simple-md -d /path/to/simulation \
    -n complex.prmtop -c complex.inpcrd \
    --start-step 6
```

## GaMD Command

A new command `gamd_md` is available for running Amber GaMD in a two-stage workflow.

### What it does
- Requires a restart file containing velocities from standard Amber equilibration.
- Uses `simple_md` NPT equilibration output as the GaMD starting point.
- Runs GaMD equilibration to collect boost statistics and ramp up the bias.
- Runs GaMD production with the fixed bias parameters computed in equilibration.
- Supports batch directories with one `prmtop` / restart pair per task subdirectory.

### Required upstream workflow
1. `prepare_amber` to build `prmtop` and `inpcrd`
2. `run_amber simple_md` through NPT equilibration
   - the recommended starting restart is `step05_npt2.rst`
3. `run_amber gamd_md` using this restart as input

### Usage example

```bash
lazydock-cli run-amber simple_md -d /path/to/batch \
    -n complex.prmtop -c complex.inpcrd \
    --start-step 1 --end-step 5

lazydock-cli run-amber gamd_md -d /path/to/batch \
    -n complex.prmtop -c step05_npt2.rst \
    --eq-nstlim 3000000 --prod-nstlim 100000000 \
    --ntcmd 1000000 --nteb 1000000 \
    --igamd 3 --sigma0P 6.0 --sigma0D 6.0
```

### Phase control
- `--start-phase equil` / `--end-phase prod`: run both phases (default).
- `--start-phase equil --end-phase equil`: run only GaMD equilibration.
- `--start-phase prod --end-phase prod`: run only GaMD production, requiring an existing `gamd_equil.rst` from a prior equilibration.

### Key GaMD options
- `--igamd`: GaMD boost mode (see below). Default `3` is recommended for protein systems.
- `--ntcmd`: steps for statistic collection during equilibration.
- `--nteb`: steps for bias ramp-up during equilibration.
- `--ntave`: running average window size for boost parameter updates.
- `--sigma0P`, `--sigma0D`: standard deviation thresholds for total and dihedral potential boosts.

### igamd boost modes

| igamd | Mode | Description | When to use |
|-------|------|-------------|-------------|
| 1 | Total potential boost | Adds boost to total potential energy V only | Small systems (< 30k atoms); simple systems |
| 2 | Dihedral boost | Adds boost to dihedral energy Vd only | Systems where dihedral transitions are the bottleneck |
| 3 | Dual boost | Boosts both V and Vd independently | **Recommended default** for most protein/complex systems |
| 4 | Total potential (lower bound) | Variant of igamd=1 with lower-bound threshold | Similar to igamd=1 but with different threshold strategy |

**Why igamd=3 is recommended**: In dual-boost mode, the total potential and dihedral potential are boosted separately with independent parameters. Each boost component is smaller and more likely to follow a Gaussian distribution, which is essential for accurate reweighting. With igamd=1, the total potential boost can be very large in big systems, causing the boost distribution to deviate significantly from Gaussian (high anharmonicity).

### GaMD phase workflow

GaMD simulation consists of two sequential phases:

1. **Equilibration** (`gamd_equil.in`): Collects statistics on the potential energy distribution and gradually ramps up the boost potential. The boost parameters (kP, kD, EP, ED) are computed during this phase.
2. **Production** (`gamd_prod.in`): Runs with fixed boost parameters from equilibration. The boost potential is applied throughout, and `gamd.log` records ΔV_P and ΔV_D for each step.

Phase timing:
- `ntcmdprep`: Pre-collection steps (default 250000, 0.5 ns)
- `ntcmd`: Statistic collection steps (default 1000000, 2 ns)
- `ntebprep`: Pre-ramp steps (default 250000, 0.5 ns)
- `nteb`: Bias ramp-up steps (default 1000000, 2 ns)
- Total equilibration = ntcmdprep + ntcmd + ntebprep + nteb

### Anharmonicity and reweighting reliability

The accuracy of GaMD reweighting depends on the boost potential ΔV following a Gaussian distribution. The anharmonicity index γ quantifies this deviation:

$$\gamma = S_{max} - S_{\Delta V}$$

where $S_{max}$ is the entropy of a Gaussian with the same variance, and $S_{\Delta V}$ is the actual distribution entropy. This is the **entropy-based** definition from Miao et al., JCTC 2015.

| γ value | Quality | Typical example |
|---------|---------|-----------------|
| γ < 0.01 | Excellent | Alanine dipeptide (~0.002) |
| 0.01 ≤ γ < 0.05 | Acceptable | M3 receptor (~0.007) |
| γ ≥ 0.05 | Poor | Large systems with igamd=1 |

When γ is high, the cumulant expansion to 2nd order ("Gaussian Approximation") becomes unreliable, and reweighted free energies may be inaccurate.

**Common causes of high anharmonicity**:
1. **igamd=1 on large systems**: Total potential energy has very broad fluctuations, making ΔV_P non-Gaussian
2. **sigma0P too high**: Larger sigma0P allows more aggressive boosting, which can distort the distribution
3. **Insufficient equilibration**: Boost parameters may not have converged

**How to improve anharmonicity**:
1. Use `igamd=3` (dual boost) instead of `igamd=1` — each component is smaller and more Gaussian
2. Reduce `sigma0P` (e.g., from 6.0 to 3.0-4.0) — less aggressive boosting improves Gaussianity at the cost of reduced acceleration
3. Extend equilibration (`ntcmd`, `nteb`) to ensure boost parameters converge
4. For ligand-binding studies, consider LiGaMD (selective boost on ligand)

### sigma0P/sigma0D tuning

These parameters control the upper limit of the boost standard deviation. Lower values produce more conservative (smaller) boosts that are more likely to be Gaussian, while higher values allow more aggressive acceleration.

| sigma0P/sigma0D | Effect | Recommended |
|-----------------|--------|-------------|
| 3.0 | Conservative boost, better Gaussianity | Large systems (> 50k atoms) |
| 4.0-6.0 | Moderate boost | Medium systems (20k-50k atoms) |
| 6.0+ | Aggressive boost, may be non-Gaussian | Small systems only (< 20k atoms) |

Note: The Amber GaMD tutorial uses sigma0P=6.0 and sigma0D=6.0 for alanine dipeptide (a very small system). For protein-ligand complexes, lower values are typically needed.

### Troubleshooting GaMD

| Problem | Likely cause | Solution |
|---------|-------------|----------|
| kP ≈ 0 after equilibration | sigma0P too low or system too large | Increase sigma0P or use igamd=3 |
| SIGFPE crash during production | kP too small, GPU SPFP precision issue | Reduce sigma0P, or use igamd=2/3 |
| High anharmonicity (γ ≥ 0.05) | igamd=1 on large system | Switch to igamd=3, reduce sigma0P |
| Boost = 0 during production | equil/prod igamd mismatch | Ensure same igamd in both phases |
| Slow convergence of boost params | Insufficient ntcmd/nteb | Increase ntcmd and nteb |

### Current implementation boundaries
- The command assumes Amber system preparation is already complete; it does not build `prmtop` / `inpcrd` files.
- It does not perform automatic reweighting or free-energy analysis.
- It does not currently manage membrane insertion or special membrane equilibration protocols; use appropriate Amber preparation tools for membrane systems.
- Input files must be compatible and present in the same task directory.

## Notes

### Input Requirements

- **prmtop**: System topology from `prepare_amber`
- **inpcrd**: Initial coordinates from `prepare_amber`
- Files must be compatible (same atom count and order)

### System Compatibility

The protocol works for:
- **Proteins**: Standard protein systems
- **Ligands**: Small molecule simulations
- **Complexes**: Protein-ligand complexes
- **Membranes**: With appropriate restraints

### Common Issues

1. **High initial energy**: Check for clashes in input structure
2. **Heating failure**: Reduce heating rate or increase restraints
3. **NPT density issues**: Extend NPT equilibration time
4. **GPU memory error**: Reduce system size or use CPU

### Best Practices

1. Always visualize trajectories after equilibration
2. Check density convergence in NPT steps
3. Monitor temperature and pressure stability
4. Use appropriate restraints for your system type
5. Extend production MD based on property convergence
