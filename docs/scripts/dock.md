<!--
 * @Date: 2024-12-10 11:22:31
 * @LastEditors: BHM-Bob 2262029386@qq.com
 * @LastEditTime: 2024-12-10 11:24:05
 * @Description: 
-->
## Command: vina
### Brief
The `vina` command is used to perform molecular docking using the Vina software. It processes configuration files in specified directories, running docking tasks in parallel.

### Parameters
- `--dir` (`-d`): The directory containing Vina configuration files (named "config.txt"). Each sub-directory represents a task. Default is the current directory.
- `--vina-name` (`-v`): The name of the Vina executable to call. Default is "vina".
- `--n-workers` (`-n`): The number of tasks to run in parallel. Default is 1.

### Behavior
The `vina` command scans the specified directory for configuration files, then runs Vina docking tasks in parallel using the specified number of workers. If a task has already been completed (indicated by the presence of a 'dock.pdbqt' file), it will be skipped.

### Notes
- Ensure that the Vina executable is accessible and correctly named.
- The number of workers should be a positive integer.

### Example
```bash
lazydock-cli dock vina -d /dir/to/configs -v vina -n 4
```

## Command: hdock
### Brief
The `hdock` command performs molecular docking using the HDOCK web service. It can process configuration files or directly use receptor and ligand file paths.

### Parameters
- `--dir` (`-d`): The directory containing HDOCK configuration files (named "config.txt"). Each sub-directory represents a task. Default is the current directory.
- `--receptor` (`-r`): The receptor PDB file name. If provided, the command will ignore "config.txt".
- `--ligand` (`-l`): The ligand PDB file name. If provided, the command will ignore "config.txt".
- `--method` (`-m`): The docking method. Currently, only "web" is supported. Default is "web".
- `--email`: The email address for the HDOCK web server. Default is None.

### Behavior
The `hdock` command processes the specified directory for configuration files or directly uses the provided receptor and ligand file paths to perform docking using the HDOCK web service. If a task has already been completed (indicated by the presence of an 'HDOCK_all_results.tar.gz' file), it will be skipped.

### Notes
- Ensure that the receptor and ligand file names are correctly specified if using the `--receptor` and `--ligand` options.
- The email address is optional but recommended for the HDOCK web service.

### Example
if you only have a receptor and a ligand file in each sub-directory, you can use the following command:
```bash
lazydock-cli dock hdock -r receptor.pdb -l ligand.pdb -m web --email user@example.com
```
if you have a vina configuration file, you can use the following command, the tool will use the receptor and ligand file names in the config.txt file:
```bash
lazydock-cli dock hpepdock -d dir/to/configs
```

## Command: hpepdock
### Brief
The `hpepdock` command performs molecular docking using the HPEPDOCK web service. It processes configuration files or directly uses receptor and ligand file paths.

### Parameters
- `--dir` (`-d`): The directory containing HPEPDOCK configuration files (named "config.txt"). Each sub-directory represents a task. Default is the current directory.
- `--receptor` (`-r`): The receptor PDB file name. If provided, the command will ignore "config.txt".
- `--ligand` (`-l`): The ligand PDB file name. If provided, the command will ignore "config.txt".
- `--method` (`-m`): The docking method. Currently, only "web" is supported. Default is "web".

### Behavior
The `hpepdock` command processes the specified directory for configuration files or directly uses the provided receptor and ligand file paths to perform docking using the HPEPDOCK web service. If a task has already been completed (indicated by the presence of an 'HPEPDOCK_all_results.tar.gz' file), it will be skipped.

### Notes
- Ensure that the receptor and ligand file names are correctly specified if using the `--receptor` and `--ligand` options.
- The `hpepdock` command currently only supports the "web" method.

### Example
if you only have a receptor and a ligand file in each sub-directory, you can use the following command:
```bash
lazydock-cli dock hpepdock -r receptor.pdb -l ligand.pdb
```
if you have a vina configuration file, you can use the following command, the tool will use the receptor and ligand file names in the config.txt file:
```bash
lazydock-cli dock hpepdock -d dir/to/configs
```

## Command: convert-result
### Brief
The `convert-result` command converts docking result files to PDB format using specified conversion tools.

### Parameters
- `--dir` (`-d`): The input directory. Default is the current directory.
- `--name` (`-n`): The input file name. Default is an empty string.
- `--input-type` (`-i`): The input file type. Default is "pdbqt,dlg".
- `--output-type` (`-o`): The output file type. Default is "pdb".
- `--suffix` (`-s`): The output file suffix. Default is an empty string.
- `--method` (`-m`): The conversion tool to use. Currently supports "lazydock" and "obabel". Default is "lazydock".
- `--n-workers`: The number of tasks to run in parallel. Default is 4.

### Behavior
The `convert-result` command scans the specified directory for input files of the specified types, then converts them to PDB format using the specified conversion tool and number of workers.

### Notes
- Ensure that the input directory and file types are correctly specified.
- The number of workers should be a positive integer.

### Example
```bash
lazydock-cli dock convert-result -d /path/to/results -n result -i pdbqt,dlg -o pdb -s _converted -m lazydock --n-workers 4
```
