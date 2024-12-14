'''
Date: 2024-08-16 09:36:38
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-14 21:45:28
Description: LazyDock Pymol Plugin
'''

import os

os.environ['MBAPY_FAST_LOAD'] = 'True'

import sys

sys.path.append(os.path.dirname(__file__))

from main import GUILauncher
from pymol import cmd

from lazydock.pml.server import VServer
from lazydock.pml.align_to_axis import align_pose_to_axis


def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('LazyDock', GUILauncher)
        return
    except Exception as e:
        print(e)
    self.menuBar.addmenuitem('Plugin', 'command', 'lazydock',
                             label = 'LazyDock', command = lambda s=self : GUILauncher(s)) 
    


def start_lazydock_server(host: str = 'localhost', port: int = 8085, quiet: int = 1):
    print(f'Starting LazyDock server on {host}:{port}, quiet={quiet}')
    VServer(host, port, not bool(quiet))
    
    
cmd.extend('start_lazydock_server', start_lazydock_server)
cmd.extend('align_pose_to_axis', align_pose_to_axis)