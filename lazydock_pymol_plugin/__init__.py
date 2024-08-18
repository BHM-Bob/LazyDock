'''
Date: 2024-08-16 09:36:38
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-18 19:02:02
Description: LazyDock Pymol Plugin
'''

import os

os.environ['MBAPY_FAST_LOAD'] = 'True'
os.environ['MBAPY_AUTO_IMPORT_TORCH'] = 'False'


from . import _utils, _autodock_utils, _interaction_utils
from .lazy_dlg import LazyDLG
from .lazy_pocket import LazyPocket


def __init__(self):
	self.menuBar.addcascademenu('Plugin', 'LazyDockPymolPlugin', 'LazyDockPymolPlugin', label = 'LazyDockPymolPlugin')
	self.menuBar.addmenuitem('LazyDockPymolPlugin', 'command', 'LazyPocket',
                          label = 'LazyPocket', command = lambda s=self : LazyPocket(s))
	self.menuBar.addmenuitem('LazyDockPymolPlugin', 'command', 'LazyDLG',
                          label = 'LazyDLG', command = lambda s=self : LazyDLG(s))