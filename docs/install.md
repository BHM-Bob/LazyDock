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

## python env compatibility
### matplotlib
- matplotlib==3.7.5
- contourpy==1.1.0