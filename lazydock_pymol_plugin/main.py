'''
Date: 2024-08-19 10:41:56
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-17 22:53:08
Description: 
'''
import os
from functools import partial
from pathlib import Path
from threading import Thread

from _nicegui.local_file_picker import local_file_picker
from lazy_dlg import LazyDLG
from lazy_plot import LazyPlot
from lazy_pml import LazyPml
from lazy_pocket import LazyPocket
from mbapy_lite.file import decode_bits_to_str, opts_file
from mbapy_lite.web import TaskPool
from nicegui import app, ui
from pymol import cmd

from lazydock.pml.thirdparty.draw_bounding_box import draw_box


class GUILauncher:
    def __init__(self, app = None, n_threads=4):
        
        self._now_molecule = cmd.get_names_of_type('object:molecule') or []
        self._now_selection = cmd.get_names_of_type('selection') + ['sele']
        self.ui_update_func = []
        ui.timer(1, self.ui_update_content_from_pymol)
        
        self.lazy_pml = LazyPml(self)
        self.lazy_pocket = LazyPocket(self)
        self.lazy_dlg = LazyDLG(self)
        self.lazy_plot = LazyPlot(self)
        
        self.build_gui()
        
        self.taskpool = TaskPool('threads', n_threads).start()
        self.host = Thread(name = 'GUI-Host', target=ui.run,
                            kwargs = dict(title='LazyDock', host = 'localhost', port=8090, reload = False),
                            daemon=False)
        self.host.start()
        
    def ui_update_content_from_pymol(self, ):
        self._now_molecule = cmd.get_names_of_type('object:molecule')
        self._now_selection = cmd.get_names_of_type('selection') + ['sele']
        
        for fn in self.ui_update_func:
            fn()
            
    def load_gpf(self, content: str):
        lines = content.split('\n')
        grid_size_line = list(filter(lambda x: x.startswith('npts'), lines))[0]
        grid_center_line = list(filter(lambda x: x.startswith('gridcenter'), lines))[0]
        getter_fn = lambda line: line.split(' ')[1:4]
        grid_center = list(map(float, getter_fn(grid_center_line)))
        grid_size = list(map(int, getter_fn(grid_size_line)))
        min_x, min_y, min_z = [grid_center[i] - grid_size[i] / 2 for i in range(3)]
        max_x, max_y, max_z = [grid_center[i] + grid_size[i] / 2 for i in range(3)]
        draw_box(min_x, min_y, min_z, max_x, max_y, max_z)

    async def uni_load(self, event):
        if event is None:
            paths = await local_file_picker(os.path.abspath('.'), upper_limit=None,
                                            multiple=True, show_hidden_files=True,
                                            file_extensions=['.pdb', '.pdbqt', '.dlg'])
            if paths is None:
                return ui.notify('No file selected')
            for path in paths:
                if path.endswith('.pdb') or path.endswith('.pdbqt'):
                    cmd.load(path)
                elif path.endswith('.dlg'):
                    await self.lazy_dlg.pose_page.load_dlg_file(None, path=path)
                elif path.endswith('.gpf'):
                    self.load_gpf(opts_file(path))
        else:
            if event.name.endswith('.pdb') or event.name.endswith('.pdbqt'):
                pdbstr = decode_bits_to_str(event.content.read())
                cmd.read_pdbstr(pdbstr, Path(event.name).stem)
            elif event.name.endswith('.dlg'):
                await self.lazy_dlg.pose_page.load_dlg_file(event, path=None)
            elif event.name.endswith('.gpf'):
                self.load_gpf(decode_bits_to_str(event.content.read()))
                    
    def build_gui(self):
        with ui.header(elevated=True).style('background-color: #3874c8'):
            ui.label('LazyDock | Pymol Plugin').classes('text-h4')
            ui.space()
            ui.upload(label='upload', auto_upload=True, on_upload=self.uni_load).style('max-width: 300px;max-height: 35px')
            ui.button('Open', on_click=partial(self.uni_load, event=None), icon='file_open')
            ui.button('Exit', on_click=app.shutdown, icon='power')
        with ui.splitter(value=10).classes('w-full h-full') as splitter:
            with splitter.before:
                with ui.tabs().props('vertical align=left active-bg-color=blue').classes('w-full') as tabs:
                    lazy_pml_tab = ui.tab('Pymol')
                    lazy_pocket_tab = ui.tab('Pocket')
                    lazy_dlg_tab = ui.tab('DLG')
                    lazy_plot_tab = ui.tab('Plot')
            with splitter.after:
                with ui.tab_panels(tabs, value=lazy_pml_tab) \
                        .props('vertical').classes('w-full h-full'):
                    with ui.tab_panel(lazy_pml_tab):
                        self.lazy_pml.build_gui()
                    with ui.tab_panel(lazy_pocket_tab):
                        self.lazy_pocket.build_gui()
                    with ui.tab_panel(lazy_dlg_tab):
                        self.lazy_dlg.build_gui()
                    with ui.tab_panel(lazy_plot_tab):
                        self.lazy_plot.build_gui()


def exec_from_cli():
    app = GUILauncher()

    
if __name__ in {"__main__", "__mp_main__"}:
    # dev code
    app = GUILauncher()