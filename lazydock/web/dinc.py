'''
Date: 2025-01-08 14:37:20
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-01-09 18:38:12
Description: 
'''
import os
import shutil
from pathlib import Path
from typing import Dict, List, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file
from mbapy_lite.web import Browser, random_sleep


def gridbox2residues(gridbox: str) -> List[str]:
    pass


def run_dock_on_DINC_ensemble(receptor_path: str, ligand_path: str, email: str,
                              box_center: Union[str, Dict[str, float]] = 'receptor', box_size: Union[str, Dict[str, float]] = 'ligand',
                              w_dir: str = None,
                              browser: Browser = None) -> Dict[str, Union[str, float]]:
    """
    Parameters:
        - receptor_path (str), path to receptor pdb file
        - ligand_path (str), path to ligand pdb file
        - box_center (Union[str, Dict[str, float]]), center of docking box, can be 'ligand' or'receptor' or a dict with keys 'x', 'y', 'z'
        - box_size (Union[str, Dict[str, float]]), size of docking box, can be 'ligand' as ligand size or a dict with keys 'x', 'y', 'z'
        - email (str), email address for DINC ensemble
        - w_dir (str), working directory, default is the parent directory of receptor_path
        - browser (Browser), a Browser object, default is None, if None, a new Browser object will be created with default options
    """
    receptor_path, ligand_path = str(Path(receptor_path).resolve()), str(Path(ligand_path).resolve())
    w_dir = w_dir or Path(receptor_path).resolve().parent
    b = browser if (browser and browser.download_path) is not None else Browser(download_path=str(w_dir))
    b.browser.set_page_load_timeout(10000)
    b.get('https://dinc-ensemble.kavrakilab.rice.edu/')
    # upload receptor and ligand pdb files
    for path, xpath in zip([receptor_path, ligand_path], ['//*[@id=":r3:"]', '//*[@id=":r0:"]']):
        btn = b.find_elements(xpath)[0]
        btn.send_keys(path)
    # set box center
    if isinstance(box_center, str):
        if box_center == 'ligand':
            b.click(element='//*[@id="root"]/div/div/div[2]/section/div/form/div[4]/div[2]/div/label[1]/span[1]', executor='element')
        elif box_center == 'receptor':
            b.click(element='//*[@id="root"]/div/div/div[2]/section/div/form/div[4]/div[2]/div/label[2]/span[1]', executor='element')
        else:
            return put_err(f'Invalid box_center: {box_center} as string, should be "ligand" or "receptor", return None')
    elif isinstance(box_center, dict):
        b.click('//*[@id="root"]/div/div/div[2]/section/div/form/div[4]/div[2]/div/label[3]/span[1]', executor='element')
        for xpath, k in zip(['//*[@id=":rs:"]', '//*[@id=":rt:"]', '//*[@id=":ru:"]'], ['x', 'y', 'z']):
            b.send_key(str(box_center[k]), element=xpath)
    else:
        return put_err(f'Invalid box_center: {box_center}, should be str or dict, return None')
    # set box size
    if isinstance(box_size, str) and box_size == 'ligand':
        b.click(element='//*[@id="root"]/div/div/div[2]/section/div/form/div[5]/div[2]/div/label[1]/span[1]', executor='element')
    elif isinstance(box_size, dict):
        b.click(element='//*[@id="root"]/div/div/div[2]/section/div/form/div[5]/div[2]/div/label[2]/span[1]', executor='element')
        for xpath, k in zip(['//*[@id=":rs:"]', '//*[@id=":rt:"]', '//*[@id=":ru:"]'], ['x', 'y', 'z']):
            b.send_key(str(box_size[k]), element=xpath)
    else:
        return put_err(f'Invalid box_size: {box_size}, should be str or dict, return None')
    # set email
    b.send_key(email, element='//*[@id=":rh:"]', executor='element')
    # submit and wait for result
    b.click(element='//*[@id="root"]/div/div/div[2]/section/div/form/button[2]', executor='element')
    # wait and click accept alert, wait for result page url
    random_sleep(15, 10)
    b.browser.switch_to.alert.accept()
    random_sleep(15, 10)
    result_url_xpath = '//*[@id="root"]/div/section/div/div[2]/a'
    while not b.wait_element(result_url_xpath):
        random_sleep(100)
        b.browser.refresh()
    b.get(b.find_elements(result_url_xpath)[0].get_attribute('href'))
    # download result
    b.click(element='//*[@id="root"]/div/section/div[1]/button', executor='element')
    result_path = os.path.join(b.download_path, 'file.zip')
    dist_path = os.path.join(w_dir, 'DINC-Ensemble_result.zip')
    shutil.move(result_path, dist_path)
    # extract result
    path_key = os.path.join('dince_res', 'analysis', 'all_info_results.csv')
    df = opts_file(dist_path, 'r:', way='tar')[path_key] # actually a tar file
    df.to_excel(os.path.join(w_dir, 'DINC_ensemble_score.xlsx'), index=False)
    return {'Score': df}


if __name__ == '__main__':
    receptor_path = 'data_tmp/docking/ligand1/receptor.pdb'
    ligand_path = 'data_tmp/docking/ligand1/ligand.pdb'
    from mbapy_lite.base import Configs
    b = Browser(options=[f"--user-agent={Configs.web.chrome_driver_path}"], use_undetected=False)
    run_dock_on_DINC_ensemble(receptor_path, ligand_path, browser=b, email='2262029386@qq.com')