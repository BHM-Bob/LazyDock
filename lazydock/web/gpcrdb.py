'''
Date: 2024-11-13 22:01:55
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-14 16:33:22
Description: 
'''
import os
from typing import Dict, List, Union

import requests
from mbapy_lite.base import put_err
from mbapy_lite.game import BaseInfo


class Residue(BaseInfo):
    def __init__(self, sequence_number: int = 0, amino_acid: str = '', protein_segment: str = '',
                 display_generic_number: str = '', alternative_generic_numbers: List[Dict[str, Union[int, str]]] = None, **kwargs):
        super().__init__()
        self.idx = sequence_number
        self.name = amino_acid
        self.seg = protein_segment
        self.generic_number = display_generic_number
        self.alternative_generic_numbers = alternative_generic_numbers or []


class Protein(BaseInfo):
    def __init__(self, entry_name: str = '', from_cache: bool = True):
        super().__init__()
        # try to load from cache first if from_cache is True and cache exists
        self.cache_path = os.path.expanduser(f'~/.lazydock/cache/gpcrdb/{entry_name}.json')
        loaded_from_cache = False
        if from_cache and os.path.exists(self.cache_path):
            try:
                self.from_json(self.cache_path)
                loaded_from_cache = True
            except Exception as e:
                put_err(f"Error loading cache: {e}")
        # if not loaded_from_cache, load from web
        if not loaded_from_cache:
            req = requests.get(f"https://gpcrdb.org/services/residues/extended/{entry_name}/")
            self.req = req.json()
            self.residuies = list(map(lambda x: Residue(**x), self.req))
        # dump to cache if not loaded from cache
        if not loaded_from_cache:
            self.to_json(self.cache_path) # makedirs inner this func


if __name__ == '__main__':
    p = Protein('cnr1_human')
    print(p.residuies)