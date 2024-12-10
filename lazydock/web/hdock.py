'''
Date: 2024-12-07 20:24:12
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-10 10:09:30
Description: 
'''
import os
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
import requests
from mbapy_lite.web import Browser, download_streamly, random_sleep


def gridbox2residues(gridbox: str) -> List[str]:
    pass


def run_dock_on_HDOCK(receptor_path: str, ligand_path: str, w_dir: str = None,
                      email: str = None, browser: Browser = None) -> Dict[str, Union[str, float]]:
    """"""
    receptor_path, ligand_path = str(Path(receptor_path).resolve()), str(Path(ligand_path).resolve())
    w_dir = w_dir or Path(receptor_path).resolve().parent
    b = browser or Browser()
    b.browser.set_page_load_timeout(10000)
    b.get('http://hdock.phys.hust.edu.cn/')
    # upload receptor and ligand pdb files
    for path, xpath in zip([receptor_path, ligand_path], ['//*[@id="pdbfile1"]', '//*[@id="pdbfile2"]']):
        btn = b.find_elements(xpath)[0]
        btn.send_keys(path)
    # set email
    if email is not None:
        b.send_keys('//*[@id="emailaddress"]', email)
    # submit and wait for result
    b.click(element='//*[@id="form1"]/table/tbody/tr[9]/td/input[1]')
    xpath_result_url = '/html/body/center/table[1]/tbody/tr/td/div/a[25]'
    while not b.find_elements(xpath_result_url):
        random_sleep(100)
        b.browser.refresh()
    # download result
    result_url = b.browser.current_url + 'all_results.tar.gz'
    result_path = os.path.join(w_dir, 'HDOCK_all_results.tar.gz')
    download_streamly(result_url, result_path, requests.Session())
    # parse score, if 'too many atoms in the lig residue UNK', result file will not contain score
    tabel = b.find_elements('/html/body/center/table[3]/tbody')[0].text
    scores = list(map(lambda x: x.split(), tabel.split('\n')))
    df_data = list(map(list, zip(*[s[-10:] for s in scores])))
    df = pd.DataFrame(df_data, columns=[' '.join(s[:-10]) for s in scores])
    df.to_excel(os.path.join(w_dir, 'HDOCK_score.xlsx'), index=False)
    return {'Score': df}


def run_dock_on_HPEPDOCK(receptor_path: str, ligand_path: str, w_dir: str = None,
                         email: str = None, browser: Browser = None) -> Dict[str, Union[str, float]]:
    """"""
    receptor_path, ligand_path = str(Path(receptor_path).resolve()), str(Path(ligand_path).resolve())
    w_dir = w_dir or Path(receptor_path).resolve().parent
    b = browser or Browser()
    b.browser.set_page_load_timeout(10000)
    b.get('http://huanglab.phys.hust.edu.cn/hpepdock/')
    # upload receptor and ligand pdb files
    for path, xpath in zip([receptor_path, ligand_path], ['//*[@id="pdbfile1"]', '//*[@id="pdbfile2"]']):
        btn = b.find_elements(xpath)[0]
        btn.send_keys(path)
    # set email
    if email is not None:
        b.send_keys('//*[@id="emailaddress"]', email)
    # submit and wait for result
    b.click(element='//*[@id="form1"]/table/tbody/tr[8]/td/input[1]')
    xpath_result_url = '/html/body/center/table[1]/tbody/tr/td/div/a[24]'
    while not b.find_elements(xpath_result_url):
        random_sleep(100)
        b.browser.refresh()
    # download result
    result_url = b.browser.current_url + 'all_results.tar.gz'
    result_path = os.path.join(w_dir, 'HPEPDOCK_all_results.tar.gz')
    download_streamly(result_url, result_path, requests.Session())
    # parse score, if 'too many atoms in the lig residue UNK', result file will not contain score
    tabel = b.find_elements('/html/body/center/table[3]/tbody')[0].text
    scores = list(map(lambda x: x.split(), tabel.split('\n')))
    df_data = list(map(list, zip(*[s[-10:] for s in scores])))
    df = pd.DataFrame(df_data, columns=[' '.join(s[:-10]) for s in scores])
    df.to_excel(os.path.join(w_dir, 'HPEPDOCK_score.xlsx'), index=False)
    return {'Score': df}


if __name__ == '__main__':
    receptor_path = 'data_tmp/docking/ligand1/receptor.pdb'
    ligand_path = 'data_tmp/docking/ligand1/ligand.pdb'
    from mbapy_lite.base import Configs
    b = Browser(options=[f"--user-agent={Configs.web.chrome_driver_path}"])
    run_dock_on_HPEPDOCK(receptor_path, ligand_path, browser=b)