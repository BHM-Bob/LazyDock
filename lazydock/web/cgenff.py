'''
Date: 2024-12-13 16:17:09
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-01-15 21:00:58
Description: 
'''

import os
from pathlib import Path
from typing import Dict, List, Union

from mbapy.base import Configs
from mbapy_lite.web import Browser, random_sleep


def get_login_browser(username: str, password: str, b: Browser = None, download_dir: str = None, timeout: float = 300.0) -> Browser:
    b = b or Browser(download_path=download_dir, use_undetected=True)
    # TODO: set popup blocker
    b.get('https://app.cgenff.com/login')
    b.find_elements('//*[@id="email"]')[0].send_keys(username)
    b.find_elements('//*[@id="password"]')[0].send_keys(password)
    b.find_elements('//*[@id="root"]/main/div/div/div/div/div[5]/button')[0].click()
    b.wait_element(element='//*[@id="root"]/div/main/div/div/div[2]/div[1]/div/p/div/div[3]/div[1]/div/button/img', timeout=timeout)
    return b


def get_result_from_CGenFF(mol2_path: str, username: str = None, password: str = None,
                           b: Browser = None, timeout: float = 30.0) -> Dict[str, Union[str, List[str]]]:
    mol2_path = str(Path(mol2_path).absolute())
    b = b or get_login_browser(username, password, download_dir=os.path.dirname(mol2_path), timeout=timeout)
    b.get('https://app.cgenff.com/homepage')
    # upload mol2 file and click "Submit" button
    b.click(element='//*[@id="root"]/div/main/div/div/div[2]/div[1]/div/p/div/div[3]/div[1]/div/button')
    upload = b.find_elements('/html/body/div[2]/div[3]/div/div[1]/div/input')[0]
    upload.send_keys(mol2_path)
    b.click(element='/html/body/div[2]/div[3]/div/div[2]/button[2]')
    # check and run
    b.click(element='//*[@id="root"]/div/main/div/div[1]/div/div[3]/div[2]/button[2]')
    b.click(element='//*[@id="root"]/div/main/div/div[2]/button[1]')
    # download result
    b.click(element='//*[@id="root"]/div/main/div/div[2]/button[2]')
    # run generate gmx result and wait, click "Download" button
    b.click(element='//*[@id="root"]/div/main/div/div[2]/button[3]')
    while b.find_elements('//*[@id="root"]/div/main/div/div[2]/button[3]')[0].text != 'DOWNLOAD GROMACS FORMAT':
        random_sleep(2)
    b.click(element='//*[@id="root"]/div/main/div/div[2]/button[3]')


if __name__ == '__main__':
    from mbapy.file import opts_file
    config = opts_file('data_tmp/config.json', way='json')
    mol2_path = 'data_tmp/pdb/ligand.mol2'
    mol2_path = str(Path(mol2_path).absolute())
    b = Browser(options = [f"--user-agent={Configs.web.chrome_driver_path}", "--disable-popup-blocking"], download_path=os.path.dirname(mol2_path), use_undetected=True)
    b = get_login_browser(config['CGENFF']['username'], config['CGENFF']['password'], b=b, download_dir=os.path.dirname(mol2_path))
    get_result_from_CGenFF(mol2_path, b=b)