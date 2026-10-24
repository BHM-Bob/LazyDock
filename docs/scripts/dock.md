
## Command: dock

### sub commmand: vina
#### Introduction
The vina command is used to perform molecular docking using the Vina software. It allows users to run Vina in parallel for multiple tasks.
#### Parameters  
- `-d` or `--dir`: The directory containing the config files for docking tasks. Each sub-directory is considered a separate task. Default is the current directory.
- `-c` or `--config-name`: The sub-string of name of the config file for each task. Default is config.txt.
- `-v` or `--vina-name`: The name of the Vina executable to call. Default is vina.
- `--vina-args`: Additional arguments to pass to the Vina executable. Default is `--log ./log.txt`.
- `-n` or `--n-workers`: The number of tasks to run in parallel. Default is 1.
#### Behavior
- The command first checks if the provided directory exists and is a valid directory.
- It then finds all config files with the specified name in the directory and its subdirectories.
- For each config file, it checks if the output file already exists. If it does, the task is skipped.
- The command runs the Vina executable with the specified arguments for each task in parallel using a thread pool with the specified number of workers.
- It waits for all tasks to complete before exiting.
#### Precautions
- Ensure that the Vina executable is properly installed and accessible.
- The config files should be in the correct format expected by Vina.
- The number of workers should be adjusted based on the system's capabilities to avoid overloading.
#### Example
```bash
lazydock-cli dock vina -d ./docking_tasks --config-name config.txt -v vina --vina-args "--log ./log.txt --exhaustiveness 8" --n-workers 4
```
This command will run Vina docking for all tasks in the docking_tasks directory using 4 workers in parallel, with the specified Vina arguments.

### sub commmand: hdock
#### Introduction
The hdock command is used to perform molecular docking using the HDOCK web server. It supports running docking tasks in batch mode.
#### Parameters
- `-d` or `--dir`: The directory containing the config files or receptor/ligand files for docking tasks. Default is the current directory.
- `-r` or `--receptor`: The sub-string of name of the receptor PDB file. If provided, it will override the config file.
- `-l` or `--ligand`: The sub-string of name of the ligand PDB file. If provided, it will override the config file.
- `--config-name`: The sub-string of name of the config file for each task. Default is config.txt.
- `-m` or `--method`: The docking method to use. Currently, only web is supported. Default is web.
- `--email`: The email address to use for the HDOCK web server.
- `-gui` or `--gui`: A flag to show the browser GUI. Default is False.
#### Behavior
- The command checks if the provided directory exists and is a valid directory.
- It finds all config files or receptor/ligand file pairs in the directory and its subdirectories.
- For each task, it checks if the output file already exists. If it does, the task is skipped.
- The command runs the HDOCK docking using the web server for each task. If the GUI flag is set, it will open a browser window.
- It waits for a random time between 3 to 5 minutes between tasks to avoid overloading the server.
#### Precautions
- Ensure that the HDOCK web server is accessible.
#### Example
```bash
lazydock-cli dock hdock -d ./docking_tasks --config-name config.txt --email user@example.com
```
This command will run HDOCK docking for all tasks in the docking_tasks directory using the config files, with the specified email address.

### sub commmand: hpepdock
#### Introduction
The hpepdock command is similar to hdock but specifically uses the HPEPDOCK web server for molecular docking.
#### Parameters
- `-d` or `--dir`: The directory containing the config files or receptor/ligand files for docking tasks. Default is the current directory.
- `-r` or `--receptor`: The sub-string of name of the receptor PDB file. If provided, it will override the config file.
- `-l` or `--ligand`: The sub-string of name of the ligand PDB file. If provided, it will override the config file.
- `--config-name`: The sub-string of name of the config file for each task. Default is config.txt.
- `-m` or `--method`: The docking method to use. Currently, only web is supported. Default is web.
- `--email`: The email address to use for the HPEPDOCK web server.
- `-gui` or `--gui`: A flag to show the browser GUI. Default is False.
#### Behavior
- The command checks if the provided directory exists and is a valid directory.
- It finds all config files or receptor/ligand file pairs in the directory and its subdirectories.
- For each task, it checks if the output file already exists. If it does, the task is skipped.
- The command runs the HPEPDOCK docking using the web server for each task. If the GUI flag is set, it will open a browser window.
- It waits for a random time between 3 to 5 minutes between tasks to avoid overloading the server.
#### Precautions
- Ensure that the HPEPDOCK web server is accessible.
#### Example
```bash
lazydock-cli dock hpepdock -d ./docking_tasks --config-name config.txt --email user@example.com
```
This command will run HPEPDOCK docking for all tasks in the docking_tasks directory using the config files, with the specified email address.

### sub commmand: dinc-ensemble
#### Introduction
The dinc-ensemble command is used to perform molecular docking using the DINC-Ensemble web server. It supports running docking tasks in batch mode and can use grid box parameters from config files.
#### Parameters
- `-d` or `--dir`: The directory containing the config files or receptor/ligand files for docking tasks. Default is the current directory.
- `-r` or `--receptor`: The sub-string of name of the receptor PDB file. If provided, it will override the config file.
- `-l` or `--ligand`: The sub-string of name of the ligand PDB file. If provided, it will override the config file.
- `--config-name`: The sub-string of name of the config file for each task. Default is config.txt.
- `-m` or `--method`: The docking method to use. Currently, only web is supported. Default is web.
- `--email`: The email address to use for the DINC-Ensemble web server. Required.
- `-gui` or `--gui`: A flag to show the browser GUI. Default is False.
- `--use-config-box`: A flag to use the grid box parameters from the config file. Default is False.
#### Behavior
- The command checks if the provided directory exists and is a valid directory.
- It finds all config files or receptor/ligand file pairs in the directory and its subdirectories.
- For each task, it checks if the output file already exists. If it does, the task is skipped.
- The command runs the DINC-Ensemble docking using the web server for each task. If the GUI flag is set, it will open a browser window.
- It waits for a random time between 3 to 5 minutes between tasks to avoid overloading the server.
- If the --use-config-box flag is set, it will use the grid box parameters from the config file for docking.
#### Precautions
- Ensure that the DINC-Ensemble web server is accessible.
- The email address is required and should be valid.
- The config file should contain the correct grid box parameters if using the --use-config-box flag.
#### Example
```bash
lazydock-cli dock dinc-ensemble -d ./docking_tasks --config-name config.txt --email user@example.com --use-config-box
```
This command will run DINC-Ensemble docking for all tasks in the docking_tasks directory using the config files and grid box parameters, with the specified email address.

### sub commmand: redock-gpu
#### Introduction
The redock-gpu command performs redocking with AutoDock-GPU. Prepared pdbqt files are converted
to a gpf with ADFRsuite's `prepare_gpf.py`, `autogrid4` builds the affinity maps, and the docking
itself runs on the GPU.

Each input complex is docked independently; poses come from `--n-poses` LGA runs.
#### Parameters
- `-d` or `--batch-dir`: Directories searched for complex PDB files. Default is `.`.
- `-n` or `--name`: Substring that the complex file name must contain.
- `-lc` or `--ligand-chain`: Chain ID of the ligand. Required.
- `-rc` or `--receptor-chain`: Chain IDs of the receptor. Default is every chain except the ligand.
- `-gb` or `--grid-buffer`: **Angstroms** added around the ligand to build the grid box. Default is 4.
- `--spacing`: Grid spacing in angstroms. **Only 0.375 is accepted**; the value is written by
  `prepare_gpf.py`, and non-default values are rejected up front.
- `-np` or `--n-poses`: Number of LGA runs, each yielding one pose. Default is 20.
- `--seed`: Random seed. Default is 0.
- `--nev`: Score evaluations (max) per LGA run. Default is 2500000.
- `--stopstd`: AutoStop energy standard deviation tolerance. Default is 0.15.
- `--heurmax`: Asymptotic heuristics evaluations limit. Default is 12000000.
- `--heuristics`: `1` (engine default) lets AutoDock-GPU derive the search method and the
  evaluation budget from the ligand torsion count; `0` makes `--nev` the exact budget.
  AutoStop cannot stop before 100% of the budget while this is `1`, so `0` is what makes
  `--stopstd` effective. Default is 1.
- `--ngen`: Maximum generations per LGA run. Default is 42000.
- `--autostop`: `1` stops a run once the best energies converge. Default is 1.
- `--gpus`: GPU ids to use. Default is `[0]`.
- `--n-task-per-gpu`: Concurrent tasks per GPU. Default is 1.
- `--prepare-backend`: `meeko` or `adfr`. Default is `meeko`.
- `--macrocycle` / `--max-ring-size`: Flexible ligand ring, see
  [macrocycle](#macrocycle-ligand-preparation---macrocycle).
- `--suffix`: Default is `_redock`.
- `--n-pdb`: How many poses to write as PDB. Default is 1.
- `--summary`: Default is `redock_summary.csv`.
- `-nw` or `--n-workers`: Parallel tasks. Default is 1.
#### Behavior
- The grid box is computed from the ligand plus `--grid-buffer` and converted to grid point
  counts (`npts = int(size / spacing)`) for `prepare_gpf.py`. See
  [grid box units](#grid-box-units-and-npts) for what that does to the realized box size.
- `--gpus` is expanded into a slot queue (`g0` x `--n-task-per-gpu`, `g1` x ...). Each worker
  takes one slot, so a saturated GPU never receives a new task; the slot is returned either way.
- Poses are written per input file as `<name>_redock.dlg`, plus `<name>_redock_<i>.pdb`.
- `--heuristics 0` is what you want when a fixed evaluation budget or `--stopstd` must be
  honoured; with the default `1` the heuristic estimate usually exceeds `--nev` anyway, so
  `--nev` still bounds the run but AutoStop stays pinned at 100%.
#### Precautions
- **The binary must support enough torsions.** AutoDock-GPU is compiled with a fixed
  `MAX_NUM_OF_ROTBONDS` and a matching internal array stride. A binary whose stride is too
  small for the ligand's torsion count can stall silently at generation 0 instead of raising,
  and there is no timeout: the task then occupies its GPU slot indefinitely. Verify the binary
  with a small test ligand before running a many-torsion one.
- `--spacing` cannot be changed here. The value is a gpf line, and this path does not forward a
  different value, so a non-default request would size the box by one value and describe it
  with another.
- Requires ADFRsuite (`prepare_gpf.py`, `autogrid4`, `pythonsh`). The bin directory is
  resolved from the located `prepare_gpf.py` and prepended to `$PATH` before each call, so it
  does not have to be on your `PATH` beforehand.
#### Example
```bash
# default: GPU 0, one task at a time
lazydock-cli dock redock-gpu -d ./docking_tasks -n complex.pdb -lc P
# two GPUs, two concurrent tasks each, fixed evaluation budget
lazydock-cli dock redock-gpu -d ./docking_tasks -n complex.pdb -lc P \
    --gpus 0 1 --n-task-per-gpu 2 --heuristics 0 --nev 5000000
# cyclic peptide with a flexible ring
lazydock-cli dock redock-gpu -d ./docking_tasks -n complex.pdb -lc P --macrocycle --gpus 0
```

### sub commmand: redock-flex-gpu
#### Introduction
The redock-flex-gpu command performs redocking with AutoDock-GPU while treating selected
receptor residues as flexible. The receptor is split by
[meeko](https://meeko.readthedocs.io/) into a rigid part (used for the affinity maps) and
flexible sidechains (passed to AutoDock-GPU via `--flexres`), so both the ligand and the
selected sidechains are sampled during docking.

The output PDB of each pose contains the rigid receptor plus that pose's own flexible
sidechain conformation and the ligand.
#### Parameters
Inherits all parameters of `redock-gpu` (`-d`, `-n`, `-lc`, `-rc`, `-gb`, `--seed`, `--nev`,
`--stopstd`, `--heurmax`, `--heuristics`, `--ngen`, `--autostop`, `-np`, `--n-pdb`, `-nw`,
`--gpus`, `--n-task-per-gpu`, `--macrocycle`), with these additional ones:
- `-fr` or `--flex-residues`: Explicit flexible residues. Accepts `<chain>:<resnum>(,<resnum>)`
  such as `A:68,50`, or `<chain>:<RESN><resnum>` such as `A:LYS68,LYS50`. Multiple values allowed.
- `-fd` or `--flex-distance`: Select receptor residues within this distance (angstrom) of any
  ligand heavy atom.
- `-fi` or `--flex-interaction`: Select receptor residues that have interactions with the ligand,
  detected by fplip.
- `--flex-interaction-mode`: Interaction types for fplip, `all` to include all modes. Default is all.
- `--flex-interaction-cutoff`: Interaction distance cutoff for fplip. Default is 4.0.
- `--prepare-backend`: Which toolkit prepares the ligand pdbqt, `meeko` or `adfr`.
  Default is `meeko`. The receptor is always split by meeko here, because ADT's
  `prepare_receptor` can not produce a rigid part + flexible sidechains pair at all
  (see [prepare-backend](#prepare-backend-receptor-and-ligand-preparation---prepare-backend)).
- `--macrocycle` / `--max-ring-size`: Flexible ligand ring, see
  [macrocycle](#macrocycle-ligand-preparation---macrocycle).
- `--suffix`: Default is `_redock_flex`.
- `--summary`: Default is `redock_flex_summary.csv`.
#### Behavior
- The three selection methods could be combined; duplicated residues are merged.
- The grid box is first computed from the ligand plus `--grid-buffer`, then automatically expanded
  to enclose all selected flexible residues.
- Affinity maps are computed only from the rigid part of the receptor.
- `--macrocycle` is supported: the glue term is part of AutoDock-GPU itself and needs no
  derived grid type, so a cyclic peptide keeps a movable ring here too.
#### Precautions
- AutoDock-GPU caps the number of atoms (512) and rotatable bonds (256); too many flexible
  residues may exceed these limits.
- Flexible residues must not overlap the ligand chain.
#### Example
```bash
# explicit residues
lazydock-cli dock redock-flex-gpu -d ./docking_tasks -n complex.pdb -lc P -rc A -fr "A:68,50"
# auto-select by distance to ligand
lazydock-cli dock redock-flex-gpu -d ./docking_tasks -n complex.pdb -lc P -rc A -fd 4
# auto-select by interaction with ligand
lazydock-cli dock redock-flex-gpu -d ./docking_tasks -n complex.pdb -lc P -rc A -fi
```

### sub commmand: redock
#### Introduction
The redock command performs redocking with AutoDock Vina. It has **two modes**, selected by
whether any flexible-receptor source is given:

| Mode | Trigger | Receptor handed to Vina |
|---|---|---|
| rigid (default) | none of `-fr` / `-fd` / `-fi` | one whole pdbqt |
| flexible | any of `-fr` / `-fd` / `-fi` | split by [meeko](https://meeko.readthedocs.io/) into a rigid part + flexible sidechains |

In flexible mode the two parts go to `Vina.set_receptor(rigid, flex)`, so both the ligand and
the selected sidechains are sampled. Vina always computes its own affinity maps internally —
there is no gpf/autogrid step in either mode.

> **Performance notice**: Vina searches on CPU, and flexible sidechains add many degrees of
> freedom. Measured on the same complex (`--exhaustiveness 2`, `--max-evals 20000`, 1 flexible
> residue), Vina took ~18 s while AutoDock-GPU finished in ~6 s with a much larger search
> budget. Expect the gap to widen with more flexible residues. Prefer `redock-flex-gpu` when
> the Vina scoring function is not specifically required.
#### Parameters
- `-d` or `--batch-dir`: Directories searched for complex PDB files. Default is `.`.
- `-n` or `--name`: Substring that the complex file name must contain.
- `-lc` or `--ligand-chain`: Chain ID of the ligand. Required.
- `-rc` or `--receptor-chain`: Chain IDs of the receptor. Default is every chain except the ligand.
- `-gb` or `--grid-buffer`: **Angstroms** added around the ligand to build the grid box. Default is 4.
- `--spacing`: Grid spacing in **angstroms**. Default is 0.375. Unlike the AutoDock-GPU paths,
  Vina honours any value here.
- `-sf` or `--score-name`: `vina`, `vinardo` or `ad4`. Default is `vina`.
- `-e` or `--exhaustiveness`: Search exhaustiveness. Default is 32.
- `-np` or `--n-poses`: Number of poses Vina writes. Default is 20.
- `--seed`: Random seed. Default is 0.
- `--max-evals`: Maximum evaluations; 0 means use Vina heuristics. Default is 0.
- `--min-rmsd`: Minimum RMSD between poses. Default is 1.
- `--cpu`: Threads for Vina. Default is unset, which means `--exhaustiveness` threads in rigid
  mode and **all cores** in flexible mode. `0` also means all cores.
- `--prepare-backend`: `meeko` or `adfr`. Default is `meeko`. See
  [prepare-backend](#prepare-backend-receptor-and-ligand-preparation---prepare-backend).
- `--macrocycle` / `--max-ring-size`: Flexible ligand ring, see
  [macrocycle](#macrocycle-ligand-preparation---macrocycle).
- `--suffix`: Default is `_redock`.
- `--n-pdb`: How many poses to write as PDB. Default is 1.
- `--summary`: Default is `redock_summary.csv`.
- `-nw` or `--n-workers`: Parallel tasks. Default is 1.

Flexible-receptor parameters (giving any of the first three switches the command into
flexible mode):
- `-fr` or `--flex-residues`: Explicit flexible residues. Accepts `<chain>:<resnum>(,<resnum>)`
  such as `A:68,50`, or `<chain>:<RESN><resnum>` such as `A:LYS68,LYS50`. Multiple values allowed.
- `-fd` or `--flex-distance`: Select receptor residues within this distance (angstrom) of any
  ligand heavy atom.
- `-fi` or `--flex-interaction`: Select receptor residues that have interactions with the ligand,
  detected by fplip.
- `--flex-interaction-mode`: Interaction types for fplip, `all` to include all modes. Default is all.
- `--flex-interaction-cutoff`: Interaction distance cutoff for fplip. Default is 4.0.
#### Behavior
- **Rigid mode**: the receptor is prepared as one pdbqt and the grid box encloses only the
  ligand plus `--grid-buffer`.
- **Flexible mode**: the three selection methods are combined and duplicated residues merged.
  The grid box is first computed from the ligand plus `--grid-buffer`, then automatically
  expanded to enclose all selected flexible residues. Affinity maps are computed only from the
  rigid part of the receptor.
- In flexible mode the receptor is **always** split by meeko, because ADT's `prepare_receptor`
  can not produce a rigid part + flexible sidechains pair at all. There `--prepare-backend`
  only selects the **ligand** backend, so `--prepare-backend adfr` with flex residues yields an
  ADT ligand with a meeko receptor.
- Output is a Vina-format PDBQT. In flexible mode each MODEL block carries the ligand **and**
  that pose's own flexible sidechain conformation, and the written complex PDB contains both.
- Glue pseudo-atoms from `--macrocycle` are stripped from the written PDB. Each row of the
  summary also reports a `closure` value (largest glue-pair distance in angstroms) when the
  ligand was prepared as a macrocycle; values above 3 Å are flagged as a warning.
#### Precautions
- **The two modes share the default output names.** Running a rigid redock and then a flexible
  one in the same directory overwrites the first run's `*_redock.pdbqt`, `*_redock_<i>.pdb` and
  summary. Pass a distinct `--suffix` / `--summary` if you need both.
- Flexible residues must not overlap the ligand chain.
- Keep the number of flexible residues small; each one adds rotatable bonds on the CPU search.
#### Example
```bash
# rigid receptor
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P
# explicit flexible residues
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P -rc A -fr "A:68,50"
# auto-select by distance to ligand, with the vinardo scoring function
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P -rc A -fd 4 -sf vinardo
# keep the two modes apart in the same directory
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P \
    -fr "A:68,50" --suffix _redock_flex --summary redock_flex_summary.csv
```

### prepare-backend: receptor and ligand preparation (`--prepare-backend`)
#### Introduction
`--prepare-backend` selects **which toolkit writes the pdbqt files** that are handed to the
docking engine. It is available on every redock subcommand (`redock`, `redock-gpu`,
`redock-flex-gpu`) and accepts:

- `meeko` — [meeko](https://meeko.readthedocs.io/) prepares both receptor and ligand.
  **This is the default.**
- `adfr` — ADFRsuite's `prepare_receptor` / `prepare_ligand` (AutoDockTools), i.e. the
  historical behaviour.

#### Behavior
- On the flexible paths (`redock` with `-fr`/`-fd`/`-fi`, and `redock-flex-gpu`) the
  **receptor** is always split by meeko, because ADT's `prepare_receptor` can not produce a
  rigid part + flexible sidechains pair at all. There the switch only selects the **ligand**
  backend.
- meeko adds the polar hydrogens itself in both cases, so a hydrogen-free receptor pdb is fine.
- Ring breaking is **disabled** on the non-`--macrocycle` path. meeko breaks rings by default
  over a 7–33 atom window, which also hits plain medium rings (measured: 7-, 10-, 12- and
  20-membered rings all get glue atoms), while ADT delivers them rigid. Leaving it on would both
  change the torsion tree and introduce glue atom types that have no AutoDock grid.
#### Precautions
- **The two backends are not score-comparable.** Their atom typing differs: on a measured
  receptor meeko classified 2 nitrogen acceptors and 4 polar hydrogens differently from ADT
  (meeko `NA=4, HD=221, N=164` vs ADT `NA=2, HD=225, N=166`; heavy-atom totals and all other
  types are identical). Vina's H-bond and N/NA terms read exactly these types, so scores shift.
  Since `meeko` is the default, **existing `redock` / `redock-gpu` results prepared by ADT are
  no longer reproducible unless you pass `--prepare-backend adfr`.**
- Pass the **same** backend when rescoring a pose that you docked with, otherwise the ligand is
  re-prepared with a different torsion tree and the rescored `total` can not be compared with
  the dock score.
- meeko must be importable by the interpreter that runs lazydock: all preparation runs
  in-process now, so there is no `--meeko-python` switch any more.
#### Example
```bash
# default: meeko prepares both receptor and ligand
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P
# reproduce a historical ADT-prepared run
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P --prepare-backend adfr
# GPU backend with a meeko ligand
lazydock-cli dock redock-gpu -d ./docking_tasks -n complex.pdb -lc P --gpus 0 1 --prepare-backend meeko
```

### grid box units and npts
#### Introduction
`--grid-buffer` and `--spacing` are both in **angstroms** on every subcommand. The two engines
however want different things: Vina takes a box size in angstrom, while AutoDock's
`prepare_gpf.py` wants a **count of grid points**. The conversion is `npts = size / spacing`.

Because `npts` must be an integer, the box actually built is quantized to the spacing: the grid
covers `npts * spacing`, which is generally **not** the size you asked for.
#### Behavior
The three writers round differently:

| Writer | Formula | Direction |
|---|---|---|
| Vina (`compute_vina_maps`) | `ceil(size / spacing)` | round **up** |
| lazydock (`prepare_gpf_and_map`) | `int(size / spacing)` | round **down** |
| meeko (`gridbox.get_gpf_string`) | `2 * int(size / (2 * spacing))` | round down, **even only** |

Measured at spacing 0.375. The Vina column is read from the live `Vina._voxels` after
`compute_vina_maps`; the AD4 column is the arithmetic `prepare_gpf_and_map` performs, which
was cross-checked against real gpf files on disk (`npts 46 46 76` at `spacing 0.375` → a
17.250 x 17.250 x 28.500 A box):

| requested | size/spacing | Vina npts | Vina box | err | AD4 npts | AD4 box | err |
|---|---|---|---|---|---|---|---|
| 20.492 | 54.645 | 55 | 20.625 | +0.133 | 54 | 20.250 | -0.242 |
| 20.692 | 55.179 | 56 | 21.000 | +0.308 | 55 | 20.625 | -0.067 |
| 20.892 | 55.712 | 56 | 21.000 | +0.108 | 55 | 20.625 | -0.267 |
| 22.492 | 59.979 | 60 | 22.500 | +0.008 | 59 | 22.125 | -0.367 |
| 22.500 | 60.000 | 60 | 22.500 | 0 | 60 | 22.500 | 0 |
| 22.875 | 61.000 | 61 | 22.875 | 0 | 61 | 22.875 | 0 |
| 23.000 | 61.333 | 62 | 23.250 | +0.250 | 61 | 22.875 | -0.125 |
| 33.300 | 88.800 | 89 | 33.375 | +0.075 | 88 | 33.000 | -0.300 |

The first four rows come from a real 143-atom ligand (`extent 12.492 A` on x) at
`--grid-buffer` 4.0 / 4.1 / 4.2 / 5.0; the rest are the arithmetic edge cases.

A full three-axis run on the same ligand with the default `--grid-buffer 4`:

```
requested   : 20.492 x 21.835 x 31.992 A
vina  npts  : (55, 59, 86) -> 20.625 x 22.125 x 32.250 A   (measured)
ad4   npts  : (54, 58, 85) -> 20.250 x 21.750 x 31.875 A   (arithmetic)
difference  : 1 point = 0.375 A per axis
```
#### Precautions
- **`A = npts * spacing` is an identity, not a way to recover your request.** It holds exactly
  for the realized box; the realized box equals the requested size only when that size is
  already a multiple of the spacing.
- **The two engines err in opposite directions.** For the same `--grid-buffer`, Vina's box is
  up to 0.75 A **larger** per axis than requested while AutoDock's is up to 0.375 A **smaller**.
- The practical consequence: vina and AD4 runs with the same `--grid-buffer` do **not** search
  the same volume. If a ligand sits right at the edge of the box, the AD4 run may exclude part
  of it. Raise `--grid-buffer` slightly when switching engines if that matters.
- meeko's even-only rule costs up to one extra point: an odd exact multiple (e.g. 22.875 A = 61
  points) becomes 60 points = 22.500 A.

### macrocycle ligand preparation (`--macrocycle`)
#### Introduction
`--macrocycle` is available on every redock subcommand (`redock`, `redock-gpu`,
`redock-flex-gpu`).
It makes meeko **cut a large ring open** in the ligand, so the ring conformation can adapt to
the pocket during docking.

This matters for cyclic peptides and other macrocycles: ADT's `prepare_ligand` has **no
ring/macrocycle option at all**, so every ring bond ends up frozen — the ring conformation
is then literally unable to move, no matter how the docking is configured. meeko instead
breaks one ring bond and attaches glue pseudo-atoms. Both engines pull those back together
during the search: Vina (>= 1.2) through its linear attraction potential, AutoDock-GPU through
its own glue term (`E = G·r`, no cutoff).

> **glue atoms are not scored by Vina**: the closure potential is only evaluated on the search
> gradient, never in `Vina.score()`. So the reported `total` stays a normal vina score.

> **This is not the same switch as `--prepare-backend`.** Since `meeko` is the default
> backend, the ligand is already prepared by meeko without `--macrocycle` — only with every
> ring kept rigid. `--macrocycle` flips ring breaking on. `--macrocycle` always wins for the
> **ligand** (it goes through meeko), while the **receptor** still follows `--prepare-backend`.
> Combining the two, `--prepare-backend adfr`, is rejected: ADT has no ring-breaking option, so
> it would hand back a ligand whose ring is still frozen.
#### Parameters
- `--macrocycle`: Enable ring breaking for the ligand. Default is False, so a default run
  keeps every ring rigid and existing results stay reproducible.
- `--max-ring-size`: Largest ring meeko may break. Default is 60.
#### Behavior
- The ligand is read as a PDB, converted with RDKit (explicit hydrogens are added), and
  prepared by meeko's Python API.
- `--max-ring-size` defaults to 60 because cyclic peptides reach 56-atom rings and meeko's
  own default of 33 misses most of them.
- Glue pseudo-atoms are removed from the output PDB, so downstream relaxation
  (PDBFixer/OpenMM) never sees a non-physical element.
- In flexible mode, each summary row also reports a `closure` column: the largest distance
  between a `CGx` atom and its paired `Gx` glue atom, in angstroms.
#### Precautions
- **meeko silently ignores rings larger than `--max-ring-size`.** There is no error and no
  warning: you get a normal-looking PDBQT, the docking runs, and you get scores — but the
  ring is frozen. Raise `--max-ring-size` until it covers the real ring size.
- Ligands prepared with `--macrocycle` have a different torsion tree than rigid-ring ones.
  vina's `0.05846 x N_rot` term therefore differs, so **scores are not comparable** across
  the two modes. To run an A/B comparison, use the same mode on both sides.
- Raising `--max-ring-size` adds ~30 rotatable bonds, which makes the search much more
  expensive. Set an explicit `--max-evals` rather than relying on vina's heuristics, which
  can ask for orders of magnitude more evaluations on high-rotor ligands.
- **A closed ring is not guaranteed just because the run finished.** AD4 does not by itself
  prefer a closed ring (measured: after removing the glue atoms, 19 of 20 poses were open),
  and a long linear peptide is a legitimate output. Read the `closure` column instead of
  assuming.
- On the AutoDock-GPU paths the binary must support enough torsions. The glue term needs no
  extra grid (`CGx` falls back to the carbon map), but a binary built with the upstream
  `MAX_NUM_OF_ROTBONDS`/stride pair can stall silently at generation 0 on a large-ring ligand
  instead of reporting an error.
- Requires `vina >= 1.2` on the vina paths: the ring closure relies on its linear attraction
  potential.
#### Example
```bash
# dock a cyclic peptide with a flexible ring
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P --macrocycle --max-evals 1000000
# flexible receptor residues + flexible ligand ring
lazydock-cli dock redock -d ./docking_tasks -n complex.pdb -lc P -fr "A:68,50" --macrocycle
# AutoDock-GPU with a flexible ring
lazydock-cli dock redock-gpu -d ./docking_tasks -n complex.pdb -lc P --macrocycle --gpus 0
```

### sub commmand: convert-result
Introduction
The convert-result command is used to convert docking result files from various formats to PDB format.
#### Parameters
- `-d` or `--dir`: The directory containing the input files to convert. Default is the current directory.
- `-n` or `--name`: The sub-string of name of the input file to convert. Default is an empty string, which means all files in the directory will be considered.
- `-i` or `--input-type`: The type of the input files. Multiple types can be specified separated by commas. Default is pdbqt,dlg.
- `-o` or `--output-type`: The type of the output files. Currently, only pdb is supported. Default is pdb.
- `-s` or `--suffix`: The suffix to add to the output file names. Default is an empty string.
- `-m` or `--method`: The conversion method to use. Currently, lazydock and obabel are supported. Default is lazydock.
- `--n-workers`: The number of tasks to run in parallel. Default is 4.
#### Behavior
- The command checks if the provided directory exists and is a valid directory.
- It finds all input files with the specified types in the directory and its subdirectories.
- For each input file, it converts the file to the specified output type using the chosen method.
- The command runs the conversion tasks in parallel using a process pool with the specified number of workers.
- It waits for all tasks to complete before exiting.
#### Precautions
Ensure that the conversion tools (e.g., lazydock, obabel) are properly installed and accessible.
The input files should be in the correct format expected by the conversion method.
#### Example
```bash
lazydock-cli dock convert-result -d ./docking_results -m lazydock --n-workers 2
```
This command will convert all pdbqt and dlg files in the docking_results directory to PDB format using the lazydock method with 2 workers in parallel.