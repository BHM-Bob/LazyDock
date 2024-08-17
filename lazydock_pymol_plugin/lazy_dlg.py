
import os
import re
import tkinter as tk
import uuid
from pathlib import Path
from typing import Dict, List, Union

import Pmw
from mbapy.base import put_err
from mbapy.file import opts_file

try:
    from Pmw.Pmw_2_0_1.lib.PmwComboBox import ComboBox as PmwComboBox
except:
    from Pmw import ComboBox as PmwComboBox

from pymol import cmd

if __name__ == '__main__':
    import lazydock_pymol_plugin
    from lazydock_pymol_plugin._autodock_utils import (
        ADModel, FileDialogButtonClassFactory, quickFileValidation)
else:
    from ._autodock_utils import (ADModel, FileDialogButtonClassFactory,
                                  quickFileValidation)
    

class PosePage:
    """load dlg to pose, show pose, save pose"""
    def __init__(self, notebook):
        self.config = self._get_default_config()
        self.dlg_pose: Dict[str, List[ADModel]] = {}
        self.ui_dlg_pose: Dict[str, Dict[str, Union[List[ADModel], PmwComboBox, List[str]]]] = {}
        self.now_dlg_name = None
        
        self.page = notebook.add('DLG Pose')
        self.page.pack(fill='both', expand=1, padx=3, pady=3)

        # dlg file
        self.ui_file_group = Pmw.Group(self.page, tag_text='DLG File')
        self.ui_file_group.pack(fill='both', expand=1, padx=3, pady=3)
        self.dlg_file_path = Pmw.EntryField(self.ui_file_group.interior(),
                                            labelpos='w',
                                            label_pyclass=FileDialogButtonClassFactory.get(self.set_pose_filename, filter=[("DLG File", "*.dlg")]),
                                            validate={'validator': quickFileValidation, },
                                            value=self.config['dlg_input_file'],
                                            label_text='Browse:')
        self.dlg_file_path.pack(side='left', fill='x', expand=1, padx=1, pady=5)
        self.load_dlg_file_buttonbox = Pmw.ButtonBox(self.ui_file_group.interior(), padx=0)
        self.load_dlg_file_buttonbox.pack(side='bottom', expand=1, padx=10, pady=5)
        self.load_dlg_file_buttonbox.add('Load', command=self.load_dlg_file)

        # pose viewer
        self.ui_pose_viewer_group = Pmw.Group(self.page, tag_text='Poses')
        self.ui_pose_viewer_group.pack(fill='both', expand=1, padx=10, pady=0)
        self.ui_pose_viewer_notebook = Pmw.NoteBook(self.ui_pose_viewer_group.interior())
        self.ui_pose_viewer_notebook.pack(fill='both', expand=1, padx=3, pady=3)
        
    def _get_default_config(self):
        return dict(dlg_input_file = '')
    
    def set_pose_filename(self, path: str):
        self.dlg_file_path.setvalue(path)
        
    def update_GUI_from_new_dlg(self):
        try:
            self.ui_pose_viewer_notebook.delete(self.now_dlg_name)
        except:
            pass
        self.ui_dlg_pose[self.now_dlg_name] = {'name': self.ui_pose_viewer_notebook.add(self.now_dlg_name)}
        pose_names = list(self.dlg_pose[self.now_dlg_name].keys())
        pose_names.sort(key = lambda x: int(x.split('::')[1]))
        self.ui_dlg_pose[self.now_dlg_name].update({'pose_names': pose_names})
        pose_pml_name = list(map(lambda x: x.replace(' ', '_').replace('::', '_'), pose_names))
        self.ui_dlg_pose[self.now_dlg_name].update({'pose_pml_name': pose_pml_name})
        self.ui_dlg_pose[self.now_dlg_name].update({'is_pose_show': {n:False for n in pose_names}})
        # control button
        self.pose_viewer_buttonbox = Pmw.ButtonBox(self.ui_dlg_pose[self.now_dlg_name]['name'], padx=3)
        self.pose_viewer_buttonbox.add('Show best 10', command=self.show_best_10_poses)
        self.pose_viewer_buttonbox.add('Show all', command=self.show_all_poses)
        self.pose_viewer_buttonbox.add('Hide all', command=self.hide_all_poses)
        self.pose_viewer_buttonbox.add('Delete', command=self.delete_dlg)
        self.pose_viewer_buttonbox.pack(fill='x', side='top')
        self.pose_viewer_buttonbox.alignbuttons()
        # pose list for now dlg
        self.ui_dlg_pose[self.now_dlg_name]['combo'] = Pmw.ComboBox(self.ui_dlg_pose[self.now_dlg_name]['name'],
                                                                   label_text='Docked', labelpos='nw',
                                                                   scrolledlist_items=pose_names,
                                                                   selectioncommand=self.pose_combo_box_selected,
                                                                   listbox_height=10, listbox_width=15, dropdown=False)
        self.ui_dlg_pose[self.now_dlg_name]['combo'].pack(side='left', padx=3, anchor='n')
        # now selected text info
        self.ui_dlg_pose[self.now_dlg_name]['text'] = Pmw.ScrolledText(self.ui_dlg_pose[self.now_dlg_name]['name'],
                                                                       borderframe=5,
                                                                       vscrollmode='dynamic',
                                                                       hscrollmode='dynamic',
                                                                       labelpos='n',
                                                                       label_text=self.now_dlg_name,
                                                                       text_width=150, text_height=15,
                                                                       text_wrap='none',
                                                                       text_background='white',
                                                                       text_foreground='black'
                                                                       )
        self.ui_dlg_pose[self.now_dlg_name]['text'].pack()
        # update ui_pose_viewer_notebook now page
        self.ui_pose_viewer_notebook.selectpage(self.now_dlg_name)
        
    def load_dlg_file(self):
        dlg_path = self.dlg_file_path.getvalue()
        if not os.path.exists(dlg_path):
            return put_err(f'dlg file not found: {dlg_path}')
        dlg_name = Path(dlg_path).stem
        dlg_content = opts_file(dlg_path)
        dlg_pose_lst = []
        for model in re.findall('MODEL.+?ENDMDL', dlg_content, re.DOTALL):
            model = model.replace('\nDOCKED: ', '\n')
            dlg_pose_lst.append(ADModel(lst=model.split('\n')))
        dlg_pose_lst.sort(key = lambda x: x.energy)
        self.dlg_pose[dlg_name] = {}
        for i in range(len(dlg_pose_lst)):
            dlg_pose_lst[i].poseN = i + 1
            self.dlg_pose[dlg_name][dlg_name + '::%d' % (i + 1)] = dlg_pose_lst[i]
        self.now_dlg_name = dlg_name
        self.update_GUI_from_new_dlg()
        
    def show_pose(self, dlg_name: str, pose_name: str):
        """
            - set self.ui_dlg_pose[dlg_name]['is_pose_show'][pose_name] to True
            - create pymol obj neamed pose_name.replace(' ', '_').replace('::', '_')
            - show the obj in pymol in sticks representation
        """
        self.ui_dlg_pose[dlg_name]['is_pose_show'][pose_name] = True
        view = cmd.get_view()
        pose = self.dlg_pose[dlg_name][pose_name]
        cmd.read_pdbstr(pose.as_pdb_string(), pose_name)
        pml_name = pose_name.replace(' ', '_').replace('::', '_')
        cmd.show('sticks', pml_name)
        cmd.set_view(view)
        cmd.zoom(pml_name)
        # show pose test
        self.ui_dlg_pose[self.now_dlg_name]['text'].clear()
        self.ui_dlg_pose[self.now_dlg_name]['text'].insert('end', pose.info_string())
        print(f'{pml_name} docked Energy: {pose.energy:8.2f} kcal/mol')
        
    def pose_combo_box_selected(self, value: str):
        self.now_dlg_name = self.ui_pose_viewer_notebook.getcurselection()
        self.show_pose(self.now_dlg_name, value)

    def show_all_poses(self):
        self.now_dlg_name = self.ui_pose_viewer_notebook.getcurselection()
        for pose_name in self.ui_dlg_pose[self.now_dlg_name]['pose_names']:
            self.show_pose(self.now_dlg_name, pose_name)

    def show_best_10_poses(self):
        self.now_dlg_name = self.ui_pose_viewer_notebook.getcurselection()
        for pose_name in self.ui_dlg_pose[self.now_dlg_name]['pose_names'][:10]:
            self.show_pose(self.now_dlg_name, pose_name)

    def hide_all_poses(self):
        self.now_dlg_name = self.ui_pose_viewer_notebook.getcurselection()
        for i, (pose_name, statue) in enumerate(self.ui_dlg_pose[self.now_dlg_name]['is_pose_show'].items()):
            if statue:
                cmd.delete(self.ui_dlg_pose[self.now_dlg_name]['pose_pml_name'][i])
                self.ui_dlg_pose[self.now_dlg_name]['is_pose_show'][pose_name] = False

    def delete_dlg(self):
        name = self.ui_pose_viewer_notebook.getcurselection()
        cmd.delete(name + '_*')
        self.ui_pose_viewer_notebook.delete(name)
        del self.pose_viewer_ligand_pages[name]
        del self.pose_viewer_ligand_dic[name]
        self.status_line.configure(text='Deleted %s' % name)
        how = self.score_table_radiobuttons.getvalue()
        self.update_score_table(how)
        
    

class LazyDLG:
    def __init__(self, app, _dev_mode: bool = False):
        if _dev_mode:
            parent = app
        else:
            parent = app.root

        self.dialog = Pmw.Dialog(parent,
                                 buttons=('Exit',),
                                 title = 'LazyDock Pymol Plugin - Lazy DLG',
                                 command = self._gui_withdraw)
        self.dialog.withdraw() # ???
        self.dialog.geometry('650x780')
        self.dialog.bind('<Return>', self._gui_withdraw)
        
        # the title
        self.title_label = tk.Label(self.dialog.interior(),
                                    text='LazyDock Pymol Plugin - Lazy DLG\nBHM-Bob G\n<https://github.com/BHM-Bob/LazyDock/>',
                                    background='orange', foreground='white')
        self.title_label.pack(expand=0, fill='both', padx=4, pady=4)
        
        # main notebook
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=3, pady=3)
        
        # build pages
        self.pose_page = PosePage(self.notebook)
        
        # GUI 
        self.dialog.show() # ???
        
    
    def _gui_withdraw(self, result):
        self.dialog.deactivate(result)
        self.dialog.withdraw()
        
        
# dev mode
if __name__ == '__main__':
    root = tk.Tk()
    Pmw.initialise(root)
    root.title('LazyDock Pymol Plugin - Lazy DLG - Dev Mode')

    exitButton = tk.Button(root, text = 'Exit', command = root.destroy)
    exitButton.pack(side = 'bottom')
    widget = LazyDLG(root, _dev_mode=True)
    root.mainloop()
    