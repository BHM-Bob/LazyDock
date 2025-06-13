'''
Date: 2024-11-05 16:27:41
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-06-12 16:34:54
Description: 
'''
import os
from pathlib import Path

from mbapy_lite.game import BaseInfo

CONFIG_FILE_PATH = Path('~/.lazydock/lazydock_config.json').expanduser()

class Config(BaseInfo):
    def __init__(self, ):
        super().__init__()
        self.named_paths = {
            'ligplus_dir': None,
            'sobtop_dir': None,
        }
        self.named_accounts = {
            'CGenFF': {
                'email': None,
                'password': None,},
        }
        
    def load(self, path: str = CONFIG_FILE_PATH):
        self.from_json(path)
        
    def dump(self, path: str = CONFIG_FILE_PATH):
        self.to_json(path)
        

GlobalConfig = Config()

if not os.path.exists(CONFIG_FILE_PATH.parent):
    os.makedirs(CONFIG_FILE_PATH.parent)
    
if not os.path.exists(CONFIG_FILE_PATH):
    GlobalConfig.to_json(CONFIG_FILE_PATH)
    
GlobalConfig.from_json(CONFIG_FILE_PATH)