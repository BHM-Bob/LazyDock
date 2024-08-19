'''
Date: 2024-08-19 10:41:56
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-08-19 11:07:47
Description: 
'''
from nicegui import app, ui

from .lazy_dlg import LazyDLG
from .lazy_pocket import LazyPocket


class GUILauncher:
    def __init__(self, app = None):
        self.lazy_pocket = LazyPocket(self)
        self.lazy_dlg = LazyDLG(self)
        
        self.build_gui()
        ui.run(title='LazyDock', url = 'http://localhost', port=8090)
        
    def run_gui(self, ):
        with ui.header(elevated=True).style('background-color: #3874c8'):
            ui.label('LazyDock | Pymol Plugin').classes('text-h4')
            ui.space()
            ui.button('Exit', on_click=app.shutdown, icon='power')
        with ui.splitter(value=10).classes('w-full h-56') as splitter:
            with splitter.before:
                with ui.tabs().props('vertical').classes('w-full') as tabs:
                    lazy_pocket_tab = ui.tab('Pocket')
                    lazy_dlg_tab = ui.tab('DLG')
            with splitter.after:
                with ui.tab_panels(tabs, value=lazy_pocket_tab) \
                        .props('vertical').classes('w-full h-full'):
                    with ui.tab_panel(lazy_pocket_tab):
                        self.lazy_pocket.build_gui()
                    with ui.tab_panel(lazy_dlg_tab):
                        self.lazy_pocket.build_gui()
                
    

if __name__ in {"__main__", "__mp_main__"}:
    # dev code
    GUILauncher()