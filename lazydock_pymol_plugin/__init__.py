'''
Date: 2024-08-16 09:36:38
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-16 11:05:48
Description: 
'''

from .lazy_pocket import LazyPocket


def __init__(self):
	self.menuBar.addcascademenu('Plugin', 'LazyDockPymolPlugin', 'LazyDockPymolPlugin', label = 'LazyDockPymolPlugin')
	self.menuBar.addmenuitem('LazyDockPymolPlugin', 'command', 'LazyPocket',
                          label = 'LazyPocket', command = lambda s=self : LazyPocket(s))