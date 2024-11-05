'''
Date: 2024-11-05 16:27:41
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-05 16:38:24
Description: 
'''
import os
from pathlib import Path

from mbapy_lite.game import BaseInfo


class Config(BaseInfo):
    def __init__(self, ):
        super().__init__()
        self.named_paths = {
            'ligplus_dir': None,
        }
        
configs = Config()
CONFIG_FILE_PATH = Path('~/.lazydock/lazydock_config.json').expanduser()

if not os.path.exists(CONFIG_FILE_PATH.parent):
    os.makedirs(CONFIG_FILE_PATH.parent)
    
if not os.path.exists(CONFIG_FILE_PATH):
    configs.to_json(CONFIG_FILE_PATH)
    
configs.from_json(CONFIG_FILE_PATH)