<!--
 * @Date: 2024-05-06 17:19:11
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2025-02-05 10:46:13
 * @Description: 
-->


<h1 style="text-align:center;">LazyDock: Lazy Dock, Easy Life</h1>

<p style="text-align:center;">
<img src="https://static.pepy.tech/badge/lazydock" alt="Downloads" style="display:inline-block; margin-left:auto; margin-right:auto;" />
<img src="https://img.shields.io/pypi/dm/lazydock" alt="Downloads" style="display:inline-block; margin-left:auto; margin-right:auto;" />
<img src="https://img.shields.io/github/downloads/BHM-Bob/LazyDock/total?label=GitHub%20all%20releases%20downloads" alt="Downloads" style="display:inline-block; margin-left:auto; margin-right:auto;" />
</p>

<p style="text-align:center;">
<a href="https://github.com/BHM-Bob/LazyDock/"><img src="https://img.shields.io/github/repo-size/BHM-Bob/LazyDock" alt="repo size" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://github.com/BHM-Bob/LazyDock/"><img src="https://img.shields.io/github/languages/code-size/BHM-Bob/LazyDock" alt="code size" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://github.com/BHM-Bob/LazyDock/releases"><img src="https://img.shields.io/github/v/release/BHM-Bob/LazyDock?label=GitHub%20Release" alt="GitHub release (latest by date)" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://github.com/BHM-Bob/LazyDock/releases"><img src="https://img.shields.io/github/commit-activity/m/BHM-Bob/LazyDock" alt="GitHub Commit Activity" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a><img src="https://img.shields.io/github/last-commit/BHM-Bob/LazyDock?label=GitHub%20Last%20Commit" alt="GitHub last commit" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
</p>

<p style="text-align:center;">
<a href="https://pypi.org/project/lazydock/"><img src="https://img.shields.io/pypi/status/lazydock?label=PyPI%20Status" alt="PyPI Status" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://pypi.org/project/lazydock/"><img src="https://img.shields.io/pypi/v/lazydock?label=PyPI%20Release" alt="PyPI" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://pypi.org/project/lazydock/"><img src="https://img.shields.io/pypi/pyversions/lazydock" alt="python versions" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
</p>

<p style="text-align:center;">
<img alt="GitHub language count" src="https://img.shields.io/github/languages/count/BHM-Bob/LazyDock">
<a href="https://github.com/BHM-Bob/LazyDock/"><img src="https://img.shields.io/readthedocs/ba-py" alt="docs" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://github.com/BHM-Bob/LazyDock/"><img src="https://img.shields.io/github/license/BHM-Bob/LazyDock" alt="license" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
</p>

<p style="text-align:center;">
<a href="https://github.com/BHM-Bob/BA_PY/"><img src="https://img.shields.io/badge/Windows-support-<COLOR>.svg" alt="windows" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://github.com/BHM-Bob/BA_PY/"><img src="https://img.shields.io/badge/Linux-support-<COLOR>.svg" alt="linux" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
</p>


<h2 style="text-align:center;">Get start</h2>

## install 
#### install lazydock python package itself 
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

#### install lazydock-pymol-plugin
In pymol plugin installation GUI, install `path-to-site-packages/lazydock_pymol_plugin/__init__.py`.

#### install plip for lazydock
```bash
pip install plip -U --no-deps
```

#### install MD-TASK for lazydock
Now all needed functionality of MD-TASK for lazydock is included in lazydock_md_task package. Already installed by requirements.

#### install MD-DaVis for lazydock
```bash
pip install md-davis -U --no-deps
```

#### install gmx_MMPBSA for lazydock
```bash
conda install -c conda-forge "mpi4py==4.0.1" "ambertools<=23.3"
conda install -c conda-forge "gromacs<=xxxx.x" pocl 
pip install gmx_MMPBSA -U --no-deps
```
gromacs version is needed to match the version for your Gromacs installation, others' version are from gmx_MMPBSA's docs.

## Docs
- [read the docs: lazydock.rtfd.io](https://lazydock.rtfd.io)
