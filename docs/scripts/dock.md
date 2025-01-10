
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