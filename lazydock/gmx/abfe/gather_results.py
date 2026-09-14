import glob
import json
import os
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd

from lazydock.gmx.abfe.tools import PathLike, sum_uncertainty_propagation


def get_fep_stats(replica_paths: List[PathLike]) -> dict:
    """Takes all the replica path and extract free energy statistics

    Parameters
    ----------
    replica_paths : List[PathLike]
        A list with all the replica paths

    Returns
    -------
    dict
        A dictionary with keywords:

            #. <estimator>: the average value of the estimator
            #. <estimator>_sem: Standard error of the mean
            #. <estimator>_uncertainty_propagation: Propagate the uncertainties after the average.
            This use the estimated uncertainties from alchemlyb (Check :func:`lazydock.gmx.abfe.fep_analysis.run_alchemlyb`)
            #. <estimator>_num_replicas: The number of replicas employed.

    """
    # Get for each estimator the corresponded values and standard deviations
    estimator_result = dict()
    for replica_path in replica_paths:
        df = pd.read_csv(replica_path, index_col=0)
        for estimator in df.columns:
            if estimator in estimator_result:
                estimator_result[estimator].append((df.loc['value', estimator], df.loc['std_dev', estimator]))
            else:
                estimator_result[estimator] = [(df.loc['value', estimator], df.loc['std_dev', estimator])]

    # Build the final result
    final_result = dict()
    for estimator in estimator_result:
        mean_value = np.mean([value_error[0] for value_error in estimator_result[estimator]])
        mean_std_dev = sum_uncertainty_propagation(
            errors=[value_error[1] for value_error in estimator_result[estimator]],
            coefficients=len(estimator_result[estimator]) * [1/len(estimator_result[estimator])]
        )
        # Save results
        final_result[estimator] = mean_value
        final_result[f"{estimator}_sem"] = pd.Series([value_error[0] for value_error in estimator_result[estimator]]).sem(ddof=1)
        final_result[f"{estimator}_uncertainty_propagation"] = mean_std_dev
        final_result[f"{estimator}_num_replicas"] = len(estimator_result[estimator])

    return final_result


def get_all_fep_dgs(root_folder_path: PathLike, out_csv: PathLike = None) -> pd.DataFrame:
    """Get the independent FEP results and gather them.
    Average and standard error of the mean are reported.

    Parameters
    ----------
    root_folder_path : PathLike
        Where the simulation run. Inside it should be the files:
        root_folder_path + "/replica_*/dG_results.csv" (LazyDock layout
        2026-09-14: no {ligand_name}/ layer).
    out_csv : PathLike, optional
        If given a pandas.DataFrame will be written as csv file, by default None

    Returns
    -------
    pd.DataFrame
        All gather results. If there are not dG_results.csv; It will return an empty DataFrame
    """
    # Get all dG_results.csv files (one per replica)
    root_folder_path = str(root_folder_path)
    dg_files = []
    for dg_file in glob.glob(root_folder_path + "/replica_*/dG_results.csv"):
        dg_files.append(os.path.normpath(dg_file))

    gathered_results = []
    if dg_files:
        statistics = get_fep_stats(dg_files)
        # single-ligand run: the ligand column is the working directory name
        statistics['ligand'] = os.path.basename(os.path.normpath(root_folder_path))
        gathered_results.append(statistics)
        gathered_results = pd.DataFrame(gathered_results)
        # Put the column 'ligand' at the beginning
        columns = ['ligand'] + [col for col in gathered_results.columns if col != 'ligand']
        gathered_results = gathered_results[columns]
        # Safe data on request
        if out_csv:
            gathered_results.to_csv(out_csv)
        return gathered_results
    else:
        print(f"There is not dG_results.csv yet on {root_folder_path}/replica_*")
        return pd.DataFrame()


def get_raw_fep_data(root_folder_path: PathLike, out_csv: PathLike = None) -> pd.DataFrame:
    """Generate raw dat for an FEP calculation.

    Parameters
    ----------
    root_folder_path : PathLike
        Where the simulation run. Inside it should be the files:
        root_folder_path + "/replica_*/complex/fep/ana/dg_complex_contributions.json"
        (LazyDock layout 2026-09-14: no {ligand_name}/ layer).
    out_csv : PathLike, optional
        If given a pandas.DataFrame will be written as csv file, by default None

    Returns
    -------
     pd.DataFrame
        Raw FEP data, all the contributions for all ligand/replicas
    """
    sample_data = []
    root_folder_path = Path(root_folder_path).resolve()
    for item1 in root_folder_path.iterdir():
        if item1.is_dir() and item1.name.startswith('replica_'):
            replica = item1.stem
            complex_json = item1/"complex/fep/ana/dg_complex_contributions.json"
            ligand_json = item1/"ligand/fep/ana/dg_ligand_contributions.json"
            if complex_json.is_file() and ligand_json.is_file():
                with open(complex_json, 'r') as cj:
                    complex_data = json.load(cj)
                with open(ligand_json, 'r') as lj:
                    ligand_data = json.load(lj)

                sample_data.append(
                    [
                        os.path.basename(root_folder_path),  # ligand = workdir name
                        replica,
                        complex_data['vdw']['MBAR']['value'],
                        complex_data['coul']['MBAR']['value'],
                        complex_data['bonded']['MBAR']['value'],
                        ligand_data['vdw']['MBAR']['value'],
                        ligand_data['coul']['MBAR']['value'],

                        complex_data['vdw']['TI']['value'],
                        complex_data['coul']['TI']['value'],
                        complex_data['bonded']['TI']['value'],
                        ligand_data['vdw']['TI']['value'],
                        ligand_data['coul']['TI']['value'],

                        complex_data['boresch'],

                        complex_data['vdw']['MBAR']['error'],
                        complex_data['coul']['MBAR']['error'],
                        complex_data['bonded']['MBAR']['error'],
                        ligand_data['vdw']['MBAR']['error'],
                        ligand_data['coul']['MBAR']['error'],

                        complex_data['vdw']['TI']['error'],
                        complex_data['coul']['TI']['error'],
                        complex_data['bonded']['TI']['error'],
                        ligand_data['vdw']['TI']['error'],
                        ligand_data['coul']['TI']['error'],
                    ]
                )
    df = pd.DataFrame(
        sample_data,
        columns=[
            'ligand', 'replica',
            'mbar_complex_vdw_value', 'mbar_complex_coul_value', 'mbar_complex_bonded_value', 'mbar_ligand_vdw_value', 'mbar_ligand_coul_value',
            'ti_complex_vdw_value', 'ti_complex_coul_value', 'ti_complex_bonded_value', 'ti_ligand_vdw_value', 'ti_ligand_coul_value',
            'boresch',
            'mbar_complex_vdw_error', 'mbar_complex_coul_error', 'mbar_complex_bonded_error', 'mbar_ligand_vdw_error', 'mbar_ligand_coul_error',
            'ti_complex_vdw_error', 'ti_complex_coul_error', 'ti_complex_bonded_error', 'ti_ligand_vdw_error', 'ti_ligand_coul_error',
        ]
        )
    if out_csv:
        df.to_csv(out_csv)
    return df


if __name__ == '__main__':
    pass