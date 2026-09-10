<!--
 * @Date: 2025-02-06 16:18:51
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2025-03-13 21:34:26
 * @Description: 
-->
## install lazydock python package itself 
- pypi
```bash
pip install lazydock -U
```
- source
```bash
pip install git+https://github.com/BHM-Bob/BA_PY.git
```
```bash
pip install git+https://gitee.com/BHM-Bob/BA_PY.git
```
- dependencies
  1. mbapy: **coding base**.
  2. pygame: needed by `mbapy.BaseInfo`, part of coding base.
  3. nicegui: needed by `lazydock.pml_plugin` to run web ui.
  4. openbabel-wheel: needed by `lazydock.scripts.prepare_gmx.py` to center molecule, can be commented out if installed openbabel software.
  5. networkx: needed by `lazydock.gmx.thirdparty.cgenff_charmm2gmx.py`, can be changed to other version if compatible.
  6. DuIvyTools: needed by `lazydock.scripts.run_gmx` & `ana_gmx` to visualize gmx output.
  7. lazydock_md_task: needed by `lazydock.scripts.ana_gmx` to run MD-TASK.
  8. compas: needed by `lazydock.pml.align_to_axis` to calculate bonding box.
  9. md-davis: needed by `lazydock.scripts.ana_gmx` to run MD-DaVis.
  10. expect: a software for running shell command, needed by `lazydock.gmx.run.Gromacs`.
  11. autodocktools-py3: a python package as AutoDockTooks.

## install lazydock-pymol-plugin
In pymol plugin installation GUI, install `path-to-site-packages/lazydock_pymol_plugin/__init__.py`.

## install lazydock dependencies
#### install plip for lazydock
```bash
conda install -c conda-forge openbabel
pip install plip -U --no-deps
```

#### install MD-TASK for lazydock
Now all needed functionality of MD-TASK for lazydock is included in lazydock_md_task package. Already installed by requirements.

#### install MD-DaVis for lazydock
```bash
pip install biopandas h5py
pip install md-davis -U --no-deps
```

#### install gmx_MMPBSA for lazydock
**Warning**: gmx_MMPBSA is better to be installed in a separate conda environment to avoid lib conflicts.
```bash
conda install -c conda-forge "mpi4py==4.0.1" "ambertools<=23.3" "parmed==4.2.2" pocl
pip install gmx_MMPBSA -U --no-deps
```
Version are from gmx_MMPBSA's docs.


#### install expect for lazydock
```bash
sudo apt install expect
```

#### install autodocktools-py3
```bash
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3"
```

## install ABFE (FEP) dependencies
ABFE (prepare-abfe + run-abfe) needs a python stack for small-molecule parameterization
(OpenFF), Boresch restraints and alchemlyb analysis. All extras live in the
[requirements.json](../requirements.json).

### 1. conda channel packages (OpenFF + OpenMM stack)
`openff-toolkit` (and `OpenMM`) are **conda packages** - install them from
conda-forge *before* the pip extras, so the heavy compiled dependencies
(amberff, rdkit, numpy, scipy) are resolved by conda rather than pip:

```bash
conda install -c conda-forge openff-toolkit openmm  # versions from https://docs.openforcefield.org/projects/toolkit/en/stable/installation.html
```

If you already run in a purely pip-managed environment, `openff-toolkit` is
also available on PyPI (`pip install openff-toolkit`), but the conda route is
the officially recommended one and avoids toolchain conflicts.

### 2. pip extras (from the `fep` requirements group)
```bash
# manually:
pip install 'MDRestraintsGenerator>=0.2.1' --no-deps
# bellow are already installed by requirements.json
pip install 'alchemlyb>=2.0.0' 'pymbar>=4.0.1' 'parmed>=4.1.0' 'toff==0.2.0' rdkit openmmforcefields pyyaml
```

> **Known issue**: `MDRestraintsGenerator==0.2.1` declares `scipy < 1.8` in its
> metadata (its older internal implementation), but the code only uses
> `scipy.stats` circular statistics and works fine with modern scipy
> (verified with scipy 1.12 + python 3.10/3.12). 
> `toff==0.2.0` also has incomplete metadata: it only declares pyyaml/parmed/
> rdkit but **hard-imports** `openff.toolkit` and `openmmforcefields` at module
> top-level - install those explicitly (step 1/2 above).

## python env compatibility
### matplotlib
- matplotlib==3.7.5
- contourpy==1.1.0