
import os
import re
import tkinter as tk
import uuid
from pathlib import Path
from typing import Dict, List, Union

from mbapy.base import put_err
from mbapy.file import opts_file, decode_bits_to_str

from nicegui import ui

from pymol import cmd


from _autodock_utils import ADModel
from _utils import uuid4
    

class LazyPose:
    """load dlg to pose, show pose, save pose"""
    def __init__(self, app):
        self._app = app
        
        self.dlg_pose: Dict[str, Dict[str, Union[List[ADModel]]]] = {}
        self.now_dlg_name = None
        self.last_pose_name = None
        self.now_pose_name = None
        
        self.sort_pdb_by_res = False
        self.show_best_per = 0
        
    def build_gui(self):
        """called by LazyDLG.__init__"""
        with ui.row().classes('w-full'):
            ui.upload(label = 'Load DLG', multiple=True, auto_upload=True,
                    on_upload=self.load_dlg_file, on_multi_upload=self.load_dlg_file).props('no-caps')
            ui.checkbox('sort pdb by res', value=self.sort_pdb_by_res).bind_value_to(self, 'sort_pdb_by_res')
        self.ui_make_dlg_file_list()
        
    @ui.refreshable
    def ui_make_dlg_file_list(self):
        with ui.row().classes('w-full h-full'):
            # file list
            with ui.tabs(value = self.now_dlg_name).props('vertical active-bg-color=blue').classes('w-1/8 h-2/3') as file_tabs:
                for n in self.dlg_pose:
                    ui.tab(n).props('no-caps').classes('w-full').tooltip(n)
            file_tabs.on('click', self.handle_dlg_tab_click)
            if not self.now_dlg_name:
                return file_tabs
            # pose list
            self.now_pose_name = self.dlg_pose[self.now_dlg_name]['name_lst'][0]
            with ui.tabs(value = self.now_pose_name).props('vertical').classes('w-1/8 h-2/3') as pose_tabs:
                for n in self.dlg_pose[self.now_dlg_name]['pose']:
                    text_col = 'text-blue' if self.dlg_pose[self.now_dlg_name]['is_show'][n] else 'text-black'
                    self.dlg_pose[self.now_dlg_name]['ui_tab'][n] = \
                        ui.tab(n).classes(f'w-full {text_col} no-caps').tooltip(n).on('click', self.handle_pose_tab_click)
            pose_tabs.bind_value_to(self, 'now_pose_name')
            # pose info
            with ui.column().classes('w-2/5'):
                ui.label().props('no-cap').bind_text_from(self, 'now_dlg_name', lambda x: f'DLG: {x}')
                # pose control
                with ui.row().classes('w-full'):
                    ui.button('show best 10', on_click=self.show_best_10_poses).props('no-caps').classes('flex flex-grow')
                    ui.number('best percentage', min=0, max=100, step=0.1, value=self.show_best_per).bind_value_to(self, 'show_best_per').classes('flex flex-grow')
                    ui.button('show percentage', on_click=self.show_percentage_poses).props('no-caps').classes('flex flex-grow')
                    ui.button('show all', on_click=self.show_all_poses).props('no-caps').classes('flex flex-grow')
                    ui.button('hide all', on_click=self.hide_all_poses).props('no-caps').classes('flex flex-grow')
                # pose info
                with ui.card().classes('w-full h-full'):
                    self.build_pose_info_gui()
                    
    def handle_dlg_tab_click(self, event):
        click_dlg = event.sender._props['model-value']
        self.now_dlg_name = click_dlg
        self.now_pose_name = self.dlg_pose[click_dlg]['name_lst'][0]
        self.ui_make_dlg_file_list.refresh()
                    
    def handle_pose_tab_click(self, event):
        click_pose = event.sender._props['name']
        if self.dlg_pose[self.now_dlg_name]['is_show'][click_pose]:
            self.now_pose_name = None
            self.hide_pose(self.now_dlg_name, click_pose)
        else:
            self.now_pose_name = click_pose
            self.show_pose(self.now_dlg_name, click_pose)
        self.build_pose_info_gui.refresh()
                    
    @ui.refreshable
    def build_pose_info_gui(self):
        ui.label().props('no-cap').bind_text_from(self, 'now_pose_name', lambda x: f'Pose: {x}')
        if self.now_dlg_name and self.now_pose_name:
            pml_name = self.dlg_pose[self.now_dlg_name]['pml_name'][self.now_pose_name]
            energy = self.dlg_pose[self.now_dlg_name]['pose'][self.now_pose_name].energy
            ui.label(f'{pml_name} docked Energy: {energy:10.4f} kcal/mol')
            ui.textarea().bind_value_from(self, 'now_dlg_name', lambda x: self.dlg_pose[x]['pose'][self.now_pose_name].info_string()).classes('w-full h-full flex flex-grow')
        
    async def load_dlg_file(self, event):
        if hasattr(event, 'name'):
            names, contents = [event.name], [event.content]
        else:
            names, contents = event.names, event.contents
        for file_name, content in zip(names, contents):
            dlg_name = Path(file_name).stem
            if dlg_name in self.dlg_pose:
                continue # TODO: don't kown why, will load again on same file, cause read null at the second time
            dlg_content = decode_bits_to_str(content.read())
            dlg_pose_lst = []
            for model in re.findall('MODEL.+?ENDMDL', dlg_content, re.DOTALL):
                model = model.replace('\nDOCKED: ', '\n')
                dlg_pose_lst.append(ADModel(model.split('\n'), self.sort_pdb_by_res))
            dlg_pose_lst.sort(key = lambda x: x.energy)
            self.dlg_pose[dlg_name] = {'pose':{}, 'is_show':{}, 'name_lst': [], 'pml_name': {}, 'ui_tab': {}}
            for i in range(len(dlg_pose_lst)):
                dlg_pose_lst[i].poseN = i + 1
                pose_name = dlg_name + '::%d' % (i + 1)
                self.dlg_pose[dlg_name]['pose'][pose_name] = dlg_pose_lst[i]
                self.dlg_pose[dlg_name]['is_show'][pose_name] = False
                self.dlg_pose[dlg_name]['name_lst'].append(pose_name)
                self.dlg_pose[dlg_name]['pml_name'][pose_name] = pose_name.replace(' ', '_').replace('::', '_')
            self.now_dlg_name = dlg_name
        self.ui_make_dlg_file_list.refresh()
        
    def show_pose(self, dlg_name: str, pose_name: str):
        """
            - set self.ui_dlg_pose[dlg_name]['is_pose_show'][pose_name] to True
            - create pymol obj neamed pose_name.replace(' ', '_').replace('::', '_')
            - show the obj in pymol in sticks representation
        """
        self.dlg_pose[dlg_name]['is_show'][pose_name] = True
        ui_tab = self.dlg_pose[dlg_name]['ui_tab'][pose_name]
        ui_tab.classes(add='text-blue', remove='text-black')
        view = cmd.get_view()
        pose = self.dlg_pose[dlg_name]['pose'][pose_name]
        cmd.read_pdbstr(pose.as_pdb_string(), pose_name)
        pml_name = self.dlg_pose[dlg_name]['pml_name'][pose_name]
        cmd.show('sticks', pml_name)
        cmd.set_view(view)
        cmd.zoom(pml_name)

    def show_all_poses(self):
        for pose_name in self.dlg_pose[self.now_dlg_name]['name_lst']:
            self.show_pose(self.now_dlg_name, pose_name)

    def show_best_10_poses(self, n = 10):
        for pose_name in self.dlg_pose[self.now_dlg_name]['name_lst'][:n]:
            self.show_pose(self.now_dlg_name, pose_name)
            
    def show_percentage_poses(self):
        n = int(len(self.dlg_pose[self.now_dlg_name]['name_lst']) * self.show_best_per / 100)
        self.show_best_10_poses(n)
        
    def hide_pose(self, dlg_name: str, pose_name: str):
        if self.dlg_pose[dlg_name]['is_show'][pose_name]:
            ui_tab = self.dlg_pose[dlg_name]['ui_tab'][pose_name]
            ui_tab.classes(add='text-black', remove='text-blue')
            self.dlg_pose[dlg_name]['is_show'][pose_name] = False
            cmd.delete(self.dlg_pose[dlg_name]['pml_name'][pose_name])

    def hide_all_poses(self):
        for pose_name in self.dlg_pose[self.now_dlg_name]['is_show']:
            self.hide_pose(self.now_dlg_name, pose_name)

    def delete_dlg(self):
        cmd.delete(self.now_dlg_name + '_*')
        del self.dlg_pose[self.now_dlg_name]
        self.now_dlg_name = self.now_dlg_name[list(self.dlg_pose.keys())[0]] if self.now_dlg_name else None
        
        
class InteractionPage:
    def __init__(self, app):
        self._app = app
        
    def build_gui(self):
        pass
        
    def ui_select_receptor(self, receptor_name: str):
        pass
    
    def ui_select_ligand(self, ligand_name: str):
        pass
    

class LazyDLG:
    def __init__(self, app, _dev_mode: bool = False):
        self._app = app
        
        self.pose_page = LazyPose(self._app)
        self.analysis_page = InteractionPage(self._app)
        
    def build_gui(self):
        with ui.tabs().classes('w-full').props('align=left active-bg-color=blue') as tabs:
            self.ui_loader_tab = ui.tab('DLG Pose Loader').props('no-caps')
            self.ui_analysis_tab = ui.tab('DLG Pose Analysis').props('no-caps')
        with ui.tab_panels(tabs, value=self.ui_loader_tab).classes('w-full'):
            with ui.tab_panel(self.ui_loader_tab):
                self.pose_page.build_gui()
            with ui.tab_panel(self.ui_analysis_tab):
                self.analysis_page.build_gui()
        # return self
        return self
        
        
# dev mode
if __name__ in {"__main__", "__mp_main__"}:
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    
    from main import GUILauncher
    GUILauncher()
    