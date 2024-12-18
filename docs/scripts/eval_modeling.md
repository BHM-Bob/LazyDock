# Model Evaluation CLI Tool

## Overview
This Python script serves as a command-line interface (CLI) tool for evaluating protein models using various web servers. It automates the process of submitting models, retrieving evaluation results, and compiling them into a user-friendly format.

## Command: model_eval
### Brief
The `model_eval` command evaluates protein models in a specified directory using multiple web servers.

### Arguments
- `-d`, `--dir`: The input directory containing the models to be evaluated. Defaults to the current working directory.
- `-n`, `--name`: A substring to filter model filenames. Only files containing this substring will be processed. Defaults to an empty string, meaning all files are considered.
- `-t`, `--type`: The file type of the models to search for in the input directory. Defaults to "pdb".
- `--servers`: A comma-separated list of web servers to use for evaluation. Available servers include ProQ, VoroMQA, ProSA, MolProbity, ProQ3, SAVES, QMEANDisCo, and QMEAN. Defaults to all supported servers.
- `--n-workers`: The number of worker threads for parallel processing. Defaults to 4.
- `-r`, `--recursive`: A flag indicating whether to search the input directory recursively. Defaults to False.
- `--disable-cache`: A flag to disable caching of results. If enabled, the tool will re-evaluate models even if results already exist. Defaults to False.

### Behavior
1. The script processes the input directory and filters files based on the specified type and name substring.
2. It initializes a task pool with the specified number of worker threads.
3. For each model file, it checks if cached results exist (unless caching is disabled). If cached results are found, it skips evaluation and adds the result path to the task list.
4. If no cached results are found, it alters the chain code of the model using PyMOL and saves the modified file.
5. The script then submits a task to the task pool for evaluation using the specified web servers.
6. Once all tasks are completed, it gathers the results, compiles them into a pandas DataFrame, and exports the DataFrame to an Excel file named `eval_<timestamp>.xlsx` in the input directory.

### Notes
- Ensure that the input directory exists and contains the desired model files.
- The script assumes that the necessary dependencies, including PyMOL and the `mbapy_lite` and `lazydock` libraries, are installed and properly configured.
- Be aware of the terms of use and rate limits of the web servers to avoid any potential issues.

### Example
```bash
lazydock-cli eval-modeling -d /path/to/models -n my_model -t pdb --servers ProQ,ProSA --n-workers 8 -r