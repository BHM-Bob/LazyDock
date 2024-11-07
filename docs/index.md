<!--
 * @Date: 2024-08-24 22:24:36
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-11-07 22:38:21
 * @Description: 
-->


<h1 style="text-align:center;">Lazy Dock</h1>

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
<a href="https://github.com/BHM-Bob/LazyDock/"><img src="https://camo.githubusercontent.com/c292429e232884db22e86c2ea2ea7695bc49dc4ae13344003a95879eeb7425d8/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f57696e646f77732d3030373844363f7374796c653d666f722d7468652d6261646765266c6f676f3d77696e646f7773266c6f676f436f6c6f723d7768697465" alt="windows" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
<a href="https://github.com/BHM-Bob/LazyDock/"><img src="https://camo.githubusercontent.com/7eefb2ba052806d8a9ce69863c2eeb3b03cd5935ead7bd2e9245ae2e705a1adf/68747470733a2f2f696d672e736869656c64732e696f2f62616467652f4c696e75782d4643433632343f7374796c653d666f722d7468652d6261646765266c6f676f3d6c696e7578266c6f676f436f6c6f723d626c61636b" alt="linux" style="display:inline-block; margin-left:auto; margin-right:auto;" /></a>
</p>


<h2 style="text-align:center;">Get start</h2>

## install 
Now, lazydock only support pypi install:  
```
pip install lazydock
```


## Contains

### pymol utils - pml
#### autodock
- [autodock](pml/autodock_utils.md): parse autodock results to pdb pose. 
#### interaction 
- [pymol](pml/interaction_utils.md): calculate interactions between receptor and poses.  
- [plip](pml/plip_interaction.md): calculate interactions between receptor and poses.
- [ligplus](pml/ligplus_interaction.md): calculate interactions between receptor and poses.
#### pymol utils
- [server](pml/server.md): provide a server-client interface to communicate with pymol.  
- [shader](pml/shader.md): provide a shader for pymol. 
- [RRCS](pml/rrcs.md)  

### pyroseeta utils - pyrt
- [pose](pyrt/pose.md): provide a wrapper pose class to opt pose. 

### web - web-server-tools
- [model_eval](web/model_eval.md): evaluate model build result.  
- [gen_box](web/gen_box.md): generate pocket/gridbox info.

### lazydock pymol plugin
#### command line user
- a script entry named `lazydock-pml` is provided in the python/Scripts folder.
#### pymol user
- install the plugin in pymol manually.




## release note
- [0.6.1](release_note/0.6.1.md)
- [0.6.0](release_note/0.6.0.md)
- [0.5.0](release_note/0.5.0.md)
- [0.4.0](release_note/0.4.0.md)
- [0.3.0](release_note/0.3.0.md)
- [0.2.0](release_note/0.2.0.md)