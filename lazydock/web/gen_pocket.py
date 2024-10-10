
import os
from pathlib import Path
from typing import Dict, List, Union

from pymol import cmd
from mbapy.file import opts_file
from mbapy.web import BROWSER_HEAD, Browser, random_sleep

if __name__ == '__main__':
    from lazydock.utils import uuid4
    from lazydock.pml.thirdparty.draw_bounding_box import draw_bounding_box
else:
    from ..utils import uuid4
    from ..pml.thirdparty.draw_bounding_box import draw_bounding_box


def get_pocket_box_from_ProteinPlus(receptor_path: str, ligand_path: str,
                                    browser: Browser = None):
    receptor_path = str(Path(receptor_path).resolve())
    ligand_path = str(Path(ligand_path).resolve())
    b = browser or Browser({'profile.default_content_settings.popups': 0,
                            'download.default_directory': os.path.dirname(receptor_path)})
    b.get('https://proteins.plus/')
    btn = b.find_elements('//*[@id="pdb_file_pathvar"]')[0]
    btn.send_keys(receptor_path)
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
    return None

def parse_pocket_box_from_ProteinPlus(result_path: str, k: Union[int, List[int]] = None,
                                      reinitialize: bool = False, draw_box: bool = True, _cmd = None):
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
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = _cmd.get_extent(pocket_name)
    _cmd.delete(pocket_name)
    return {
        'box_center': [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2],
        'box_size': [maxX - minX, maxY - minY, maxZ - minZ],
        'box_vertices': [[minX, minY, minZ], [maxX, maxY, maxZ]],
    }
    
    
__all__ = [
    'get_pocket_box_from_ProteinPlus',
    'parse_pocket_box_from_ProteinPlus',
]

if __name__ == '__main__':
    # from mbapy.base import Configs
    # b = Browser(options=['--no-sandbox', f"--user-agent={Configs.web.chrome_driver_path}"])
    # get_pocket_box_from_ProteinPlus('data_tmp/pdb/RECEPTOR.pdb', 'data_tmp/pdb/LIGAND.sdf', browser=b)
    print(parse_pocket_box_from_ProteinPlus('data_tmp/pdb/POCKET.zip', [1], True))
    print(parse_pocket_box_from_ProteinPlus('data_tmp/pdb/POCKET.zip', [1, 2, 3], True))
    print(parse_pocket_box_from_ProteinPlus('data_tmp/pdb/POCKET.zip', [1, 2, 3, 4, 5], True))