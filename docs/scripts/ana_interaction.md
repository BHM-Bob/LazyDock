<!--
 * @Date: 2024-11-30 20:44:04
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-12-02 22:51:59
 * @Description: 
-->
## Command: simple-analysis

### Introduction
The `simple-analysis` command is designed to perform a simple analysis on docking results. It calculates interactions between a receptor and a ligand using various methods and outputs the results in a specified format.

### Parameters

- `-r`, `--receptor` (str): The path or file name of the receptor PDB file, which will be loaded by PyMOL.
- `-l`, `--ligand` (str): The path or file name of the docking result, supporting Vina (pdbqt) and AutoDock (dlg) formats.
- `-bd`, `--batch-dir` (str): The directory containing multiple sub-folders, each with docking result files.
- `--method` (str): The method to use for interaction analysis. Defaults to 'pymol'. Choices are 'pymol', 'ligplus', and 'plip'.
- `--mode` (str): The interaction mode. Multiple modes can be separated by commas. All methods support the 'all' mode. Default is 'all'.
  - pymol: 'all', 'bond distances', 'polar contact', 'all distance_exclusion', 'centroids', 'pi-pi and pi-cation', 'pi-pi interactions', 'pi-cation interactions', 'ratio distance_exclusion'  
  - ligplus: 'all', 'Hydrogen Bonds', 'Non-bonded Interactions'  
  - plip: 'all', 'Hydrophobic Interactions', 'Hydrogen Bonds', 'Water Bridges', 'Salt Bridges', 'pi-Stacking', 'pi-Cation Interactions', 'Halogen Bonds', 'Metal Complexes'  
- `--cutoff` (float): The distance cutoff for interaction calculation. Default is 4.
- `--output-style` (str): The output style. Currently, only 'receptor' is supported, which outputs in the format 'resn resi distance'.
- `--hydrogen-atom-only` (flag): only consider hydrogen bond acceptor and donor atoms, this only works when method is pymol. Default is False.

### Behavior
The `simple-analysis` command processes the input receptor and ligand files or directories, calculates the interactions based on the specified method and mode, and outputs the results in the specified format. If a batch directory is provided, it will process all sub-folders containing docking result files.

### Notes
- Ensure that the receptor and ligand files exist and are in the correct format.
- If using the batch directory option, make sure each sub-folder contains both receptor and ligand files.
- The interaction modes should be one of the supported modes for the chosen method.
- The distance cutoff should be a positive float value.

### Examples

1. Analyze a single receptor-ligand pair using the PyMOL method with default settings:
   ```
   lazydock-cli ana-interaction simple-analysis -r receptor.pdbqt -l dock.pdbqt
   ```

2. Analyze a batch of docking results using the LigPlus method with a specific mode and distance cutoff:
   ```
   lazydock-cli ana-interaction simple-analysis -bd data_tmp\docking --method ligplus --mode 'all' --cutoff 4
   ```

3. Analyze a single receptor-ligand pair using the PLIP method with a custom interaction mode:
   ```
   lazydock-cli ana-interaction simple-analysis -r receptor.pdbqt -l dock.pdbqt --method plip --mode "Hydrogen Bonds,Salt Bridges"
   ```
