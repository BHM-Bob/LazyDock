'''
Date: 2024-08-16 09:36:38
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-19 23:30:48
Description: LazyDock Pymol Plugin
'''

import os

os.environ['MBAPY_FAST_LOAD'] = 'True'
os.environ['MBAPY_AUTO_IMPORT_TORCH'] = 'False'

import sys

sys.path.append(os.path.dirname(__file__))

import _autodock_utils, _interaction_utils, _utils
from lazy_dlg import LazyDLG
from lazy_pocket import LazyPocket
from main import GUILauncher


def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('LazyDock', GUILauncher)
        return
    except Exception as e:
        print(e)
    self.menuBar.addmenuitem('Plugin', 'command', 'lazydock',
                             label = 'LazyDock', command = lambda s=self : GUILauncher(s)) 