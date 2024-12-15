'''
Date: 2024-08-16 09:36:38
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-15 19:34:35
Description: LazyDock Pymol Plugin
'''

import os

os.environ['MBAPY_FAST_LOAD'] = 'True'

import sys

sys.path.append(os.path.dirname(__file__))

from pymol import cmd

# do all init in frist level code block
import lazydock.pml_plugin._init_ as lazydock_pml_plugin_init
# GUILauncher is the main function of the plugin, so just init in this file.
from lazydock.pml_plugin.main import GUILauncher


def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('LazyDock', GUILauncher)
        return
    except Exception as e:
        print(e)
    self.menuBar.addmenuitem('Plugin', 'command', 'lazydock',
                             label = 'LazyDock', command = lambda s=self : GUILauncher(s)) 
    
