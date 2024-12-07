'''
Date: 2024-12-07 21:04:01
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-07 22:18:50
Description: 
'''

import os
from pathlib import Path
from typing import Dict

import pandas as pd
import requests
from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import opts_file
from mbapy_lite.web import Browser, random_sleep


def get_score_from_SwissADME(lig_SMILES: str, result_dir: str, browser: Browser = None, timeout: int = 1200, **kwargs) -> Dict[str, str]:
    """return Score in dict"""
    b = browser or Browser()
    b.browser.fullscreen_window()
    b.get('http://www.swissadme.ch/index.php')
    b.send_key(key=lig_SMILES, element='//*[@id="smiles"]')
    b.execute_script("arguments[0].removeAttribute('disabled')", b.find_elements(element='//*[@id="submitButton"]')[0])
    b.click(element='//*[@id="submitButton"]')
    xpath_result_img1 = '//*[@id="content"]/div[11]/div[1]/div[3]/div[2]/img'
    if not b.wait_element(xpath_result_img1, timeout=timeout):
        return put_err('Timeout', {})
    # download img1
    img1_url = b.find_elements(xpath_result_img1)[0].get_attribute('src')
    img1_path = os.path.join(result_dir, 'img1.png')
    opts_file(img1_path, 'wb', data=requests.get(img1_url).content)
    # cap img2
    b.click(element='//*[@id="content"]/div[8]/span[1]/button')
    img2 = b.find_elements('//*[@id="placeholder"]/canvas[2]')[0].screenshot_as_png
    img2_path = os.path.join(result_dir, 'img2.png')
    opts_file(img2_path, 'wb', data=img2)
    # download scv file
    csv_link = b.find_elements('//*[@id="content"]/div[7]/a[1]')[0].get_attribute('href')
    csv_path = os.path.join(result_dir,'result.csv')
    opts_file(csv_path, 'wb', data=requests.get(csv_link).content)
    # transpose csv file
    csvT_path = os.path.join(result_dir,'resultT.csv')
    pd.read_csv(csv_path).T.to_csv(csvT_path, index=True)
    return {'img1_path': img1_path, 'img2_path': img2_path, 'csv_path': csv_path, 'csvT_path': csvT_path}


def get_score_from_SwissTargetPrediction(lig_SMILES: str, result_dir: str, browser: Browser = None, timeout: int = 1200, **kwargs) -> Dict[str, str]:
    """return Score in dict"""
    result_dir = str(Path(result_dir).resolve())
    if not browser:
        b = Browser(download_path=result_dir)
        put_log(f'using download path: {result_dir}')
    else:
        b = browser
    b.browser.fullscreen_window()
    b.get('http://www.swisstargetprediction.ch/')
    b.send_key(key=lig_SMILES, element='//*[@id="smilesBox"]')
    b.execute_script("arguments[0].removeAttribute('disabled')", b.find_elements(element='//*[@id="submitButton"]')[0])
    b.click(element='//*[@id="submitButton"]')
    xpath_result_img1 = '//*[@id="content"]/div[3]/img'
    if not b.wait_element(xpath_result_img1, timeout=timeout):
        return put_err('Timeout', {})
    # download img1
    img1_url = b.find_elements(xpath_result_img1)[0].get_attribute('src')
    img1_path = os.path.join(result_dir, 'img1.png')
    opts_file(img1_path, 'wb', data=requests.get(img1_url).content)
    # cap img2
    img2_url = b.find_elements('//*[@id="placePieChart"]/a/img')[0].get_attribute('src')
    img2_path = os.path.join(result_dir, 'img2.png')
    opts_file(img2_path, 'wb', data=requests.get(img2_url).content)
    # download scv file
    b.click(element='//*[@id="exportButtons"]/div/button[2]')
    b.click(element='//*[@id="exportButtons"]/div/button[3]')
    return {'img1_path': img1_path, 'img2_path': img2_path}


if __name__ == '__main__':
    lig_SMILES = 'CC1=C(C2=C3N1[C@@H](COC3=CC=C2)CN4CCOCC4)C(=O)C5=CC=CC6=CC=CC=C65'
    result_dir = 'data_tmp/web'
    get_score_from_SwissTargetPrediction(lig_SMILES, result_dir)