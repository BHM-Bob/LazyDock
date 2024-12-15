'''
Date: 2024-11-02 22:09:17
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-03 16:17:45
Description: 
'''
import io
import sys
from typing import List

from nicegui import ui
from pymol import CmdException, api, cmd


class LazyPml:
    def __init__(self, app, _dev_mode: bool = False):
        
        self._app = app
        
        sys.stdout = self.captured_output = io.StringIO()
        
        self.cmd_string = ''
        self.ui_log = None
        
    def build_cmd_gui(self):
        with ui.row().classes('w-full'):
            ui.input(label='pymol python api call').bind_value_to(self, 'cmd_string').classes('flex flex-grow')
            ui.button('execute', on_click=self.execute_cmd).classes('flex flex-grow')
        with ui.row():
            self.ui_log = ui.log()
            
    def execute_cmd(self):
        try:
            exec(self.cmd_string, globals(), locals())
            self.ui_log.push(self.captured_output.getvalue().split('\n')[-2])
        except CmdException as e:
            self.ui_log.push(str(e))
        
    def build_gui(self):
        with ui.tabs().classes('w-full').props('align=left active-bg-color=blue') as tabs:
            self.ui_cmd_tab = ui.tab('command').props('no-caps')
        with ui.tab_panels(tabs, value=self.ui_cmd_tab).classes('w-full'):
            with ui.tab_panel(self.ui_cmd_tab):
                self.build_cmd_gui()
        # return self
        return self
            
        
# dev mode
if __name__ in {"__main__", "__mp_main__"}:    
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/LIGAND.pdb', 'ligand')
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    
    from main import GUILauncher
    GUILauncher()