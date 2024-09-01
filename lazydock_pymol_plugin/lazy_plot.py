'''
Date: 2024-08-31 20:08:07
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-09-01 17:23:54
Description: 
'''

from typing import Dict, List

from nicegui import ui
from pymol import CmdException, api, cmd

from lazydock.pml.shader import Shader, ShaderValues

NULL_CHAIN = "''"
    
    
class ShaderPage:
    def __init__(self, app):
        self._app = app
        self._app.ui_update_func.append(self.ui_update_ui)
        # select receptor
        self.ui_cmap_name = None
        self.ui_shader = None
        self.cmap_name = 'coolwarm'
        self.shaders: Dict[str, Shader] = dict()
        self.now_shader: Shader = None
        # vitulazation control
        self.fig = None
        
    def ui_update_ui(self):
        # update shaders' repr and self.shaders
        for name in list(self.shaders.keys()):
            new_name = repr(self.shaders[name])
            if new_name!= name:
                self.shaders[new_name] = self.shaders.pop(name)
                del self.shaders[name]
        # update ui_shader
        self.ui_shader.set_options(list(self.shaders.keys()))
        
    def build_gui(self):
        with ui.splitter(value = 20).classes('w-full h-full') as splitter:
            with splitter.before:
                # shader
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('Shader')
                        with ui.row().classes('w-full'):
                            self.ui_cmap_name = ui.input(label = 'cmap', value = self.cmap_name).bind_value_to(self, 'cmap_name').classes('w-2/3')
                            ui.button('create', on_click=self.create_shader).classes('w-1/3')
                        self.ui_shader = ui.select(list(self.shaders.keys()), label = 'select a shader').bind_value_to(self, 'now_shader').classes('w-full')
                # load shader values
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('Shader values').classes('w-full')
            # interaction vitualization
            with splitter.after:
                with ui.row():
                    ui.label('Interaction: ')
                    
    def create_shader(self):
        shader = Shader(self.cmap_name)
        self.shaders[repr(shader)] = shader


class LazyPlot:
    def __init__(self, app, _dev_mode: bool = False):
        self._app = app
        
        self.shader_page = ShaderPage(self._app)
        
    def build_gui(self):
        with ui.tabs().classes('w-full').props('align=left active-bg-color=blue') as tabs:
            self.ui_shader_tab = ui.tab('Pose Shader').props('no-caps')
        with ui.tab_panels(tabs, value=self.ui_shader_tab).classes('w-full'):
            with ui.tab_panel(self.ui_shader_tab):
                self.shader_page.build_gui()
        # return self
        return self
            
        
# dev mode
if __name__ in {"__main__", "__mp_main__"}:    
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/LIGAND.pdb', 'ligand')
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    
    from main import GUILauncher
    GUILauncher()
    