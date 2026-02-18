<!--
 * @Date: 2024-05-06 17:19:11
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2025-08-06 10:12:54
 * @Description: 
-->


<h1 style="text-align:center;">LazyDock: Faster, Easier, and All-In-One Experence in Docking & MD</h1>

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
#### Details
- [How to install](docs/install.md)

## Docs
- [read the docs: lazydock.rtfd.io](https://lazydock.rtfd.io)

## Quick View
**LazyDock is a python package that provides a fast, easy, and all-in-one experience in docking and MD.**  
- For molecular docking, it has tools help pocket detect, parameters file distrbute, batchlize docking perform (via Vina, HDOCK...), docking result convertng, batchlize interaction analysis (support PyMOL, PLIP, LigPlus).  
- For molecular dynamic simulation (GROMACS), it has tools help box settings, ligand, protein(s), ligand-protein system preparation, energy minimization, equilibration, product MD run in batch mode.  
- For molecular dynamic simulation analysis, it has tools help basic analysis (GROMACS), Elastic Network Analysis, PCA, interaction analysis (support PyMOL, PLIP), RRCS, and so on, many of which has CUDA accelerate.
