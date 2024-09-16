'''
Date: 2024-09-15 22:05:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-09-16 11:16:32
Description: 
'''
import re
from pathlib import Path
from typing import Dict, List, Union

import requests
from lxml import etree
from mbapy.base import put_err
from mbapy.file import opts_file
from mbapy.web import BROWSER_HEAD, Browser, random_sleep


def get_score_from_proq(pdb_path: str, **kwargs) -> Dict[str, float]:
    """return LGscore and MaxSub in dict"""
    session = requests.Session()
    session.headers.update({'User-Agent': BROWSER_HEAD})
    response = session.post('https://proq.bioinfo.se/cgi-bin/ProQ/ProQ.cgi',
                            files = {'pdbfile': opts_file(pdb_path, 'rb')})
    xml = etree.HTML(response.content.decode())
    LGscore = float(xml.xpath('/html/body/b[1]/text()')[0])
    MaxSub = float(xml.xpath('/html/body/b[2]/text()')[0])
    return {'LGscore': LGscore, 'MaxSub': MaxSub}


def get_score_from_VoroMQA(pdb_path: str, browser: Browser = None, **kwargs) -> Dict[str, float]:
    """return Score in dict"""
    pdb_path = str(Path(pdb_path).resolve())
    b = browser or Browser()
    b.get('https://bioinformatics.lt/wtsam/voromqa/submit')
    btn = b.find_elements('//input[@type="file"]')[0]
    btn.send_keys(pdb_path)
    b.click(element='//button[text()="Submit"]')
    xpath_result_line = '//div[text()="Results:"]/..//div[2]'
    while not b.find_elements(xpath_result_line):
        random_sleep(3)
    score_text = b.find_elements(xpath_result_line)[0].text
    score = float(re.findall(r'Score: ([-\d.]+)', score_text)[0])
    return {'Score': score}


def get_score_from_ProSA(pdb_path: str, browser: Browser = None, **kwargs) -> Dict[str, Union[float, bytes]]:
    """return Z-Score, model_quality_img_res, res_score_img_res in dict"""
    pdb_path = str(Path(pdb_path).resolve())
    b = browser or Browser()
    b.get('https://prosa.services.came.sbg.ac.at/prosa.php')
    btn = b.find_elements('//input[@type="file"]')[0]
    btn.send_keys(pdb_path)
    b.click(element='//input[@type="submit"]')
    xpath_result_line = '/html/body/div[2]/span[2]'
    while not b.find_elements(xpath_result_line):
        random_sleep(3)
    score = float(b.find_elements(xpath_result_line)[0].text)
    model_quality_img_link = b.find_elements('/html/body/div[2]/p/a')[0].get_attribute('href')
    model_quality_img_res = requests.get(model_quality_img_link)
    res_score_img_link = b.find_elements('/html/body/div[3]/p/span/a')[0].get_attribute('href')
    res_score_img_res = requests.get(res_score_img_link)
    return {'Z-Score': score,
            'model_quality_img_res': model_quality_img_res,
            'res_score_img_res': res_score_img_res}


def get_score_from_ModEval(modkey: str, pdb_path: str, align_file_path: str,
                           browser: Browser = None, **kwargs) -> Dict[str, Union[float, str]]:
    """return RMSD, overlap, identity, Z-DOPE in dict"""
    pdb_path = str(Path(pdb_path).resolve())
    align_file_path = str(Path(align_file_path).resolve())
    b = browser or Browser()
    b.get('https://modbase.compbio.ucsf.edu/evaluation/')
    b.send_key(key=modkey, element='//input[@name="modkey"]')
    btn = b.find_elements('//input[@type="file" and @name="pdb_file"]')[0]
    btn.send_keys(pdb_path)
    btn = b.find_elements('//input[@type="file" and @name="alignment_file"]')[0]
    btn.send_keys(align_file_path)
    b.click(element='//input[@type="submit"]')
    while not b.find_elements('//*[@id="resulttable"]/table/tbody/tr[7]/td[2]'):
        random_sleep(3)
    rmsd_text = b.find_elements('//*[@id="resulttable"]/table/tbody/tr[7]/td[2]')[0].text.strip()
    rmsd = float(re.findall(r'[0-9.]+', rmsd_text)[0])
    overlap_text = b.find_elements('//*[@id="resulttable"]/table/tbody/tr[8]/td[2]')[0].text.strip()
    overlap = float(re.findall(r'[0-9.]+', overlap_text)[0])
    identity_text = b.find_elements('//*[@id="resulttable"]/table/tbody/tr[11]/td[2]')[0].text.strip()
    identity = float(re.findall(r'[0-9.]+', identity_text)[0])
    z_dope_text = b.find_elements('//*[@id="resulttable"]/table/tbody/tr[12]/td[2]')[0].text.strip()
    z_dope = float(re.findall(r'[0-9.]+', z_dope_text)[0])
    return {'RMSD': rmsd, 'overlap': overlap,
            'identity': identity, 'Z-DOPE': z_dope}


def get_score_from_MolProbity(pdb_path: str, browser: Browser = None, **kwargs) -> Dict[str, Union[float, str]]:
    """return Z-Score, Ramachandran Favorability, Ramachandran Outerness in dict"""
    pdb_path = str(Path(pdb_path).resolve())
    b = browser or Browser()
    b.get('http://molprobity.biochem.duke.edu/index.php')
    btn = b.find_elements('//input[@type="file" and @name="uploadFile"]')[0]
    btn.send_keys(pdb_path)
    b.click(element='//input[@type="submit" and @value="Upload >"]')
    while not b.find_elements('//input[@type="submit" and @value="Continue >"]'):
        random_sleep(3)
    b.click(element='//input[@type="submit" and @value="Continue >"]')
    result_link = b.find_elements('/html/body/table/tbody/tr[2]/td[2]/div[1]/div/table/tbody/tr/td[2]/table/tbody/tr[4]/td[2]/a')[0].get_attribute('href')
    b.get(result_link)
    b.click(element='//input[@type="submit" and @value="Run programs to perform these analyses >"]')
    xpath_z_score = '/html/body/table/tbody/tr[2]/td/div/p[1]/table/tbody/tr[5]/td[2]'
    while not b.find_elements(xpath_z_score):
        random_sleep(3)
    z_score_text = b.find_elements(xpath_z_score)[0].text
    rama_favor_text = b.find_elements('/html/body/table/tbody/tr[2]/td/div/p[1]/table/tbody/tr[4]/td[3]')[0].text
    rama_outer_text = b.find_elements('/html/body/table/tbody/tr[2]/td/div/p[1]/table/tbody/tr[3]/td[3]')[0].text
    return {'Z-Score': z_score_text,
            'Ramachandran Favorability': rama_favor_text,
            'Ramachandran Outerness': rama_outer_text}


def get_score_from_ProQ3(pdb_path: str, browser: Browser = None, **kwargs) -> Dict[str, Union[float, str]]:
    """return ProQ2D, ProQRosCenD, ProQRosCenD, ProQ3D in dict"""
    pdb_path = str(Path(pdb_path).resolve())
    b = browser or Browser()
    b.get('https://proq3.bioinfo.se/pred/')
    btn = b.find_elements('//input[@type="file" and @name="modelfile"]')[0]
    btn.send_keys(pdb_path)
    b.click(element='//input[@type="submit" and @value="Submit"]')
    xpath_result_line = '//*[@id="jobtable"]/tbody'
    while not b.find_elements(xpath_result_line):
        random_sleep(3)
    z_score_text = b.find_elements(xpath_result_line)[0].text
    lines = z_score_text.split('\n')
    return {k:float(v) if re.findall('[0-9.]+', v) else v \
            for k,v in zip(lines[0].split(' '), lines[1].split(' '))}


def get_score_from_SAVES(pdb_path: str, browser: Browser = None, **kwargs) -> Dict[str, Union[float, str]]:
    """return ProQ2D, ProQRosCenD, ProQRosCenD, ProQ3D in dict"""
    pdb_path = str(Path(pdb_path).resolve())
    b = browser or Browser()
    b.get('https://saves.mbi.ucla.edu/')
    btn = b.find_elements('//input[@type="file" and @name="pdbfile"]')[0]
    btn.send_keys(pdb_path)
    b.click(element='//input[@type="submit" and @value="Run programs"]')
    while not b.find_elements('/html/body/h2[1]/a[1]'):
        random_sleep(3)
    # start check progress
    b.click(element='//*[@id="errat"]')
    b.click(element='//*[@id="whatcheck"]')
    b.click(element='//*[@id="procheck"]')
    # wait for result
    while (not b.find_elements('//*[@id="werrat"]/div/center/center/h1')) or (not b.find_elements('//*[@id="wwhatcheck"]/div/center')) or (not b.find_elements('//*[@id="wprocheck"]/div/center')):
        random_sleep(3)
    errat_text = float(b.find_elements('//*[@id="werrat"]/div/center/center/h1')[0].text.strip())
    whatcheck_bad = len(b.find_elements('//*[@id="wwhatcheck"]/div/center//span[@class="wcitem bad"]'))
    whatcheck_mid = len(b.find_elements('//*[@id="wwhatcheck"]/div/center//span[@class="wcitem mid"]'))
    whatcheck_good = len(b.find_elements('//*[@id="wwhatcheck"]/div/center//span[@class="wcitem good"]'))
    procheck_errors = int(b.find_elements('//*[@id="wprocheck"]/div/center/ul/li[1]/span')[0].text.split(':')[-1])
    procheck_warnings = int(b.find_elements('//*[@id="wprocheck"]/div/center/ul/li[2]/span')[0].text.split(':')[-1])
    procheck_pass = int(b.find_elements('//*[@id="wprocheck"]/div/center/ul/li[3]/span')[0].text.split(':')[-1])
    # parse whatcheck result
    return {'errat': errat_text, 'whatcheck_bad': whatcheck_bad, 'whatcheck_mid': whatcheck_mid, 'whatcheck_good': whatcheck_good,
            'procheck_errors': procheck_errors, 'procheck_warnings': procheck_warnings, 'procheck_pass': procheck_pass}
    
    
_name2server = {
    'ProQ': get_score_from_proq,
    'VoroMQA': get_score_from_VoroMQA,
    'ProSA': get_score_from_ProSA,
    'ModEval': get_score_from_ModEval,
    'MolProbity': get_score_from_MolProbity,
    'ProQ3': get_score_from_ProQ3,
    'SAVES': get_score_from_SAVES,
}
    
    
def get_eval_info_from_servers(pdb_path: str, align_file_path: str = None,
                                modkey: str = None, browser: Browser = None,
                                servers: List[str] = None) -> Dict[str, Dict[str, Union[float, str, bytes]]]:
    """return a dict of evaluation results"""
    servers = servers or _name2server.keys()
    if not servers:
        return put_err('No server specified, return empty result', {})
    pdb_path = str(Path(pdb_path).resolve())
    b = browser or Browser()
    results = {}
    for server in servers:
        fn = _name2server[server]
        try:
            result = fn(pdb_path=pdb_path, align_file_path=align_file_path, modkey=modkey, browser=b)
        except Exception as e:
            result = {'Error': str(e)}
        results[server] = result
    return results


__all__ = [
    'get_score_from_proq',
    'get_score_from_VoroMQA',
    'get_score_from_ProSA',
    'get_score_from_ModEval',
    'get_score_from_MolProbity',
    'get_score_from_ProQ3',
    'get_score_from_SAVES',
    '_name2server',
    'get_eval_info_from_web',
]


if __name__ == '__main__':
    get_score_from_SAVES('data_tmp/pdb/RECEPTOR.pdb')