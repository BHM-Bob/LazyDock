
from pathlib import Path
from typing import Dict, List, Union

import pandas as pd
from mbapy_lite.base import put_log
from mbapy_lite.file import opts_file
from mbapy_lite.web import Browser, random_sleep
from pymol import cmd

if __name__ == '__main__':
    from lazydock.pml.thirdparty.draw_bounding_box import draw_bounding_box
    from lazydock.utils import uuid4
else:
    from ..pml.thirdparty.draw_bounding_box import draw_bounding_box
    from ..utils import uuid4


def get_pocket_box_from_ProteinPlus(receptor_path: str, ligand_path: str = None,
                                    browser: Browser = None):
    """
    download pocket result zip file from ProteinPlus.
    
    Parameters
        - receptor_path: str, path to the receptor pdb file.
        - ligand_path: str, path to the ligand sdf file.
        - browser: lazydock.web.Browser, browser object to use. If None, a new browser will be created.
        
    Returns
        None, the pocket result zip file will be downloaded to the download_path specified in the browser.
    """
    receptor_path = str(Path(receptor_path).resolve())
    ligand_path = str(Path(ligand_path).resolve()) if ligand_path else None
    if not browser:
        download_path = str(Path(receptor_path).resolve().parent)
        b = Browser(download_path=str(Path(receptor_path).resolve().parent))
        put_log(f'using download path: {download_path}')
    else:
        b = browser
    b.get('https://proteins.plus/')
    btn = b.find_elements('//*[@id="pdb_file_pathvar"]')[0]
    btn.send_keys(receptor_path)
    if ligand_path:
        btn = b.find_elements('//*[@id="pdb_file_userligand"]')[0]
        btn.send_keys(ligand_path)
    b.click(element='//*[@id="new_pdb_file"]/input[6]')
    # now in the result page
    b.click(element='//*[@id="headingTwo"]/h4/a/span')
    b.click(element='/html/body/div[4]/div[3]/div/div[1]/div/div/div[2]/div[2]/div/a')
    b.click(element='/html/body/div[4]/div[3]/div/div[3]/div[3]/form/div[8]/div/input')
    while not b.wait_element(element='/html/body/div[4]/div[3]/div/div[3]/div[3]/div[3]/div[1]/form/input[1]'):
        random_sleep(3)
    b.click(element='/html/body/div[4]/div[3]/div/div[3]/div[3]/div[3]/div[1]/form/input[1]')
    # distroy browser
    if not browser:
        b.browser.quit()
    return None

def parse_pocket_box_from_ProteinPlus(result_path: str, k: Union[int, List[int]] = None,
                                      reinitialize: bool = False, draw_box: bool = True,
                                      _cmd = None, method: str = 'extend'):
    """
    parse zip result file from ProteinPlus and return box information.
    
    Parameters
        - result_path: str, path to the zip file downloaded from ProteinPlus.
        - k: int or list of int, index of the pocket to be parsed. If None, all pockets will be parsed.
        - reinitialize: bool, whether to reinitialize pymol or not.
        - draw_box: bool, whether to draw box or not.
        - _cmd: pymol.cmd, pymol command object, can be pml.server.PymolAPI
        - method: str, method to calculate box center and size. 'extend' to use pymol's get_extent function, 'mean' to calculate mean coordinates of all atoms in the pocket.
        
    Returns
        box_df: pandas.DataFrame, dataframe of pocket information. Only return in the method of 'mean'.
            - 'resn_resi_index': str, residue name, residue index and atom index.
            -  'X', 'Y', 'Z': float, coordinates of the atom.
        box_info: dict, box information.
            - 'box_center': list of float, center of the box.
            - 'box_size': list of float, size of the box.
            - 'box_vertices': list of list of float, vertices of the box.
        
    """
    _cmd = _cmd or cmd
    if reinitialize:
        _cmd.reinitialize()
    # load result file
    files = opts_file(result_path, way='zip')
    pdbs = list(filter(lambda x: '_P_' in x and x.endswith('.pdb'), files))
    # check parameter k
    if k is None:
        k = [0]
    elif isinstance(k, int):
        k = [k]
    elif not isinstance(k, list):
        raise TypeError(f'k must be int or list, got {type(k)}.')
    if isinstance(k, list) and any((not any(f'_P_{i}_' in p for p in pdbs)) for i in k):
        raise ValueError(f'k contains invalid index: {k}.')
    # load pocket into pymol
    pocket_name = f'pocket_{uuid4()}'
    for idx, i in enumerate(k):
        pdb_i = list(filter(lambda x: f'_P_{i}_' in x, pdbs))[0]
        _cmd.read_pdbstr(files[pdb_i], f'{pocket_name}_{i}')
        if idx == 0:
            _cmd.select(pocket_name, f'{pocket_name}_{i}')
        else:
            _cmd.select(pocket_name, f'{pocket_name} or {pocket_name}_{i}')
    # draw box
    if draw_box:
        draw_bounding_box(pocket_name, quiet=1, _cmd = _cmd)
    # calcu box
    if method == 'extend':
        ([minX, minY, minZ], [maxX, maxY, maxZ]) = _cmd.get_extent(pocket_name)
        _cmd.delete(pocket_name)
        return {
            'box_center': [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2],
            'box_size': [maxX - minX, maxY - minY, maxZ - minZ],
            'box_vertices': [[minX, minY, minZ], [maxX, maxY, maxZ]],
        }
    elif method == 'mean':
        df = pd.DataFrame(columns=['resn_resi_index', 'X', 'Y', 'Z'])
        df.set_index(['resn_resi_index'])
        cmd.iterate_state(1, pocket_name, 'df.loc[f"{resn}_{resi}_{index}", ("X", "Y", "Z")] = [x, y, z]',
                          space=locals())
        _cmd.delete(pocket_name)
        return df, {
            'box_center': df[['X', 'Y', 'Z']].mean().tolist(),
            'box_size': (df[['X', 'Y', 'Z']].max() - df[['X', 'Y', 'Z']].min()).tolist(),
            'box_vertices': [df[['X', 'Y', 'Z']].min().tolist(), df[['X', 'Y', 'Z']].max().tolist()],
        }
    
    
__all__ = [
    'get_pocket_box_from_ProteinPlus',
    'parse_pocket_box_from_ProteinPlus',
]

if __name__ == '__main__':
    # from mbapy_lite.base import Configs
    # b = Browser(options=['--no-sandbox', f"--user-agent={Configs.web.chrome_driver_path}"],
    #             download_path=os.path.abspath('data_tmp/pdb'))
    # get_pocket_box_from_ProteinPlus('data_tmp/pdb/RECEPTOR.pdb', 'data_tmp/pdb/LIGAND.sdf', browser=b)
    df, box = parse_pocket_box_from_ProteinPlus('data_tmp/pdb/POCKET.zip', [1], True, method='mean')
    print(parse_pocket_box_from_ProteinPlus('data_tmp/pdb/POCKET.zip', [1, 2, 3], True))
    print(parse_pocket_box_from_ProteinPlus('data_tmp/pdb/POCKET.zip', [1, 2, 3, 4, 5], True))