
import gzip
import os
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Union

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from mbapy_lite.base import put_err
from mbapy_lite.file import decode_bits_to_str, opts_file
from mbapy_lite.plot import get_palette, save_show
from mbapy_lite.stats import pca
from nicegui import ui
from pymol import api, cmd

from lazydock.pml.autodock_utils import DlgFile, MyFileDialog
from lazydock.pml.interaction_utils import SUPPORTED_MODE as MODE_PML
from lazydock.pml.interaction_utils import \
    calcu_receptor_poses_interaction as calc_pml_interaction
from lazydock.pml.interaction_utils import filter_interaction_df
from lazydock.pml.ligplus_interaction import SUPPORTED_MODE as MODE_LIGPLUS
from lazydock.pml.ligplus_interaction import \
    calcu_receptor_poses_interaction as calc_ligplus_interaction

try:
    from lazydock.pml.plip_interaction import SUPPORTED_MODE as MODE_PLIP
    from lazydock.pml.plip_interaction import \
        calcu_receptor_poses_interaction as calc_plip_interaction
except ImportError:
    # if plip is not installed, use interaction_utils instead
    put_err('plip is not installed, just use pymol as plip instead')
    from lazydock.pml.interaction_utils import \
        calcu_receptor_poses_interaction as calc_plip_interaction
    MODE_PLIP = MODE_PML

from lazydock.utils import uuid4


class LazyPose:
    """load dlg to pose, show pose, save pose"""
    def __init__(self, app):
        self._app = app
        self._app.ui_update_func.append(self.ui_update_ui)
        # pose storage
        self.dlg_pose: Dict[str, DlgFile] = {}
        self.now_dlg_name = None
        self.last_pose_name = None
        self.now_pose_name = None
        # show pose control
        self.sort_pdb_by_res = True
        self.parse2std = True
        self.show_best_per = 10
        self.save_lig_with_HETATM = True
        self.save_with_each_header = False
        self.save_original_pdbstring = False
        # select receptor
        self.ui_molecule = None
        self.ui_sele = None
        self.sele_molecule = None
        self.sele_selection = None
        
    def ui_update_ui(self):
        if self.ui_sele is not None:
            self.ui_sele.set_options(self._app._now_selection)
        if self._app.lazy_dlg.pose_page.dlg_pose:
            dlg_pose: Dict[str, DlgFile] = self._app.lazy_dlg.pose_page.dlg_pose
            all_ligands = [dlg_pose[dlg].get_pose_prop('pml_name', lig) for dlg in dlg_pose for lig in dlg_pose[dlg].n2i if dlg_pose[dlg].get_pose_prop('is_show', lig, default = False)]
            self.ui_molecule.set_options(list(set(self._app._now_molecule) - set(all_ligands)))
        else:
            self.ui_molecule.set_options(self._app._now_molecule)
        
    def build_gui(self):
        """called by LazyDLG.__init__"""
        with ui.row().classes('w-full'):
            ui.upload(label = 'Load DLG', multiple=True, auto_upload=True,
                    on_upload=self.load_dlg_file, on_multi_upload=self.load_dlg_file).props('no-caps')
            ui.button('refresh GUI', on_click=self.ui_make_dlg_file_list.refresh).props('no-caps')
            # load and save control
            with ui.column():
                ui.checkbox('sort pdb by res', value=self.sort_pdb_by_res).bind_value_to(self, 'sort_pdb_by_res')
                ui.checkbox('parse2std', value=self.parse2std).bind_value_to(self, 'parse2std')
                ui.checkbox('save lig with HETATM', value=self.save_lig_with_HETATM).bind_value_to(self,'save_lig_with_HETATM')
                ui.checkbox('save with each header', value=self.save_with_each_header).bind_value_to(self, 'save_with_each_header')
                ui.checkbox('save orignal pdbstring', value=self.save_original_pdbstring).bind_value_to(self,'save_original_pdbstring')
            # sele receptor
            with ui.card().classes('w-1/2'):
                with ui.row().classes('w-full'):
                    ui.label('select a receptor')
                    ui.button('save showed complex', on_click=self.save_showed_complex).props('no-caps').bind_enabled_from(self, 'now_dlg_name')
                with ui.row().classes('w-full'):
                    self.ui_molecule = ui.select(self._app._now_molecule,
                                                label = 'select a molecule').bind_value_to(self, 'sele_molecule').classes('w-2/5').props('use-chips')
                    # ui.label('OR').classes('w-1/5')
                    self.ui_sele = ui.select(self._app._now_selection,
                                            label = 'select a selection').bind_value_to(self, 'sele_selection').classes('w-2/5').props('use-chips')
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
            self.now_pose_name = list(self.dlg_pose[self.now_dlg_name].n2i.keys())[0]
            with ui.tabs(value = self.now_pose_name).props('vertical').classes('w-1/8 h-2/3') as pose_tabs:
                for n in self.dlg_pose[self.now_dlg_name].n2i:
                    text_col = 'text-blue' if self.dlg_pose[self.now_dlg_name].get_pose_prop('is_show', n) else 'text-black'
                    ui_tab = ui.tab(n).classes(f'w-full {text_col} no-caps').tooltip(n).on('click', self.handle_pose_tab_click)
                    self.dlg_pose[self.now_dlg_name].set_pose_prop('ui_tab', ui_tab, n)
                        
            pose_tabs.bind_value_to(self, 'now_pose_name')
            # pose info
            with ui.column().classes('w-3/5'):
                ui.label().props('no-cap').bind_text_from(self, 'now_dlg_name', lambda x: f'DLG: {x}')
                # pose control
                with ui.row().classes('w-full'):
                    ui.button('show best 10', on_click=self.show_best_10_poses).props('no-caps').classes('flex flex-grow')
                    ui.number('best percentage', min=0, max=100, value=self.show_best_per).bind_value_to(self, 'show_best_per').classes('flex flex-grow')
                    ui.button('show percentage', on_click=self.show_percentage_poses).props('no-caps').classes('flex flex-grow')
                    ui.button('show all', on_click=self.show_all_poses).props('no-caps').classes('flex flex-grow')
                    ui.button('hide all', on_click=self.hide_all_poses).props('no-caps').classes('flex flex-grow')
                    ui.button('save showed', on_click=self.save_showed_pose).props('no-caps').classes('flex flex-grow')
                # pose info
                with ui.card().classes('w-full h-full'):
                    self.build_pose_info_gui()
                    
    def handle_dlg_tab_click(self, event):
        click_dlg = event.sender._props['model-value']
        self.now_dlg_name = click_dlg
        self.now_pose_name = list(self.dlg_pose[click_dlg].n2i.keys())[0]
        self.ui_make_dlg_file_list.refresh()
                    
    def handle_pose_tab_click(self, event):
        click_pose = event.sender._props['name']
        if self.dlg_pose[self.now_dlg_name].get_pose_prop('is_show', click_pose):
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
            pml_name = self.dlg_pose[self.now_dlg_name].get_pose_prop('pml_name', self.now_pose_name)
            energy = self.dlg_pose[self.now_dlg_name].get_pose(self.now_pose_name).energy
            ui.label(f'{pml_name} docked Energy: {energy:10.4f} kcal/mol')
            ui.textarea().bind_value_from(self, 'now_dlg_name', lambda x: self.dlg_pose[x].get_pose(self.now_pose_name).info_string()).classes('w-full h-full flex flex-grow')
        
    async def load_dlg_file(self, event, path: str = None):
        if path is not None:
            names, contents = [os.path.basename(path)], [open(path, 'rb')]
        elif hasattr(event, 'name'):
            names, contents = [event.name], [event.content]
        else:
            names, contents = event.names, event.contents
        for file_name, content in zip(names, contents):
            dlg_name = Path(file_name).stem
            if dlg_name in self.dlg_pose:
                continue # TODO: don't kown why, will load again on same file, cause read null at the second time
            dlg_content = decode_bits_to_str(content.read())
            content.close()
            self.dlg_pose[dlg_name] = DlgFile(content=dlg_content,
                                              sort_pdb_line_by_res=self.sort_pdb_by_res,
                                              parse2std=self.parse2std)
            self.dlg_pose[dlg_name].sort_pose()
            self.dlg_pose[dlg_name].asign_pose_name(list(map(lambda x: f'{dlg_name}::{x+1}',
                                                             range(len(self.dlg_pose[dlg_name].pose_lst)))))
            self.dlg_pose[dlg_name].asign_prop('is_show', [False] * len(self.dlg_pose[dlg_name].pose_lst))
            self.dlg_pose[dlg_name].asign_prop('pml_name', list(map(lambda x: x.replace(' ', '_').replace('::', '_'),
                                                                       self.dlg_pose[dlg_name].n2i.keys())))
            self.now_dlg_name = dlg_name
        self.ui_make_dlg_file_list.refresh()
        
    def show_pose(self, dlg_name: str, pose_name: str):
        """
            - set self.ui_dlg_pose[dlg_name]['is_pose_show'][pose_name] to True
            - create pymol obj neamed pose_name.replace(' ', '_').replace('::', '_')
            - show the obj in pymol in sticks representation
        """
        self.dlg_pose[dlg_name].set_pose_prop('is_show', True, pose_name)
        ui_tab = self.dlg_pose[dlg_name].get_pose_prop('ui_tab', pose_name)
        ui_tab.classes(add='text-blue', remove='text-black')
        view = cmd.get_view()
        pose = self.dlg_pose[dlg_name].get_pose(pose_name)
        cmd.read_pdbstr(pose.as_pdb_string(), pose_name)
        pml_name = self.dlg_pose[dlg_name].get_pose_prop('pml_name', pose_name)
        cmd.show('sticks', pml_name)
        cmd.set_view(view)
        cmd.zoom(pml_name)

    def show_all_poses(self):
        for pose_name in list(self.dlg_pose[self.now_dlg_name].n2i.keys()):
            self.show_pose(self.now_dlg_name, pose_name)

    def show_best_10_poses(self, n = 10):
        for pose_name in list(self.dlg_pose[self.now_dlg_name].n2i.keys())[:n]:
            self.show_pose(self.now_dlg_name, pose_name)
            
    def show_percentage_poses(self):
        n = int(len(self.dlg_pose[self.now_dlg_name].n2i) * self.show_best_per / 100)
        self.show_best_10_poses(n)
        
    def hide_pose(self, dlg_name: str, pose_name: str):
        if self.dlg_pose[dlg_name].get_pose_prop('is_show', pose_name):
            ui_tab = self.dlg_pose[dlg_name].get_pose_prop('ui_tab', pose_name)
            ui_tab.classes(add='text-black', remove='text-blue')
            self.dlg_pose[dlg_name].set_pose_prop('is_show', False, pose_name)
            cmd.delete(self.dlg_pose[dlg_name].get_pose_prop('pml_name', pose_name))

    def hide_all_poses(self):
        for pose_name in self.dlg_pose[self.now_dlg_name].n2i:
            self.hide_pose(self.now_dlg_name, pose_name)
            
    def save_showed_pose(self):
        if any(self.dlg_pose[self.now_dlg_name].get_pose_prop('is_show', pose_name)\
                for pose_name in self.dlg_pose[self.now_dlg_name].n2i):
            save_dir = MyFileDialog().get_ask_dir()
            if not save_dir:
                ui.notify('Please select a directory to save')
                return None
        else:
            ui.notify('Please show at least one pose')
            return None
        for pose_name in self.dlg_pose[self.now_dlg_name].n2i:
            if self.dlg_pose[self.now_dlg_name].get_pose_prop('is_show', pose_name):
                pml_name = self.dlg_pose[self.now_dlg_name].get_pose_prop('pml_name', pose_name)
                if self.save_lig_with_HETATM:
                    cmd.alter(pml_name, 'type="HETATM"')
                if self.save_original_pdbstring:
                    pose = self.dlg_pose[self.now_dlg_name].get_pose(pose_name)
                    opts_file(os.path.join(save_dir, f'{pml_name}.pdb'), 'w', data = pose.as_pdb_string())
                else:
                    cmd.save(os.path.join(save_dir, f'{pml_name}.pdb'), pml_name)
                print(f'{pml_name} saved to {save_dir}')

    def delete_dlg(self):
        cmd.delete(self.now_dlg_name + '_*')
        del self.dlg_pose[self.now_dlg_name]
        self.now_dlg_name = self.now_dlg_name[list(self.dlg_pose.keys())[0]] if self.now_dlg_name else None
        
    def save_showed_complex(self):
        if self.sele_molecule is None and self.sele_selection is None:
            ui.notify('Please select a receptor and a selection first')
            return None
        receptor = self.sele_molecule or self.sele_selection
        save_dir = MyFileDialog().get_ask_dir()
        if not save_dir:
            ui.notify('Please select a directory to save')
            return None
        for name in self.dlg_pose[self.now_dlg_name].n2i:
            if self.dlg_pose[self.now_dlg_name].get_pose_prop('is_show', name):
                pml_name = self.dlg_pose[self.now_dlg_name].get_pose_prop('pml_name', name)
                if self.save_lig_with_HETATM:
                    cmd.alter(pml_name, 'type="HETATM"')
                pdb_path = os.path.join(save_dir, f'{receptor}_{pml_name}.pdb')
                if self.save_with_each_header:
                    api.multisave(pdb_path, receptor, append = 0)
                    api.multisave(pdb_path, pml_name, append = 1)
                else:
                    tmp_sele = uuid4()
                    cmd.select(tmp_sele, f'{receptor} or {pml_name}')
                    cmd.save(pdb_path, tmp_sele)
                    cmd.delete(tmp_sele)
                print(f'{receptor} and {pml_name} saved to {pdb_path}')
                
                
class ResultPage(LazyPose):
    """plot energy curve and RMSD heatmap for pose in dlg"""
    def __init__(self, app):
        self._app = app
        self._app.ui_update_func.append(self.ui_update_ui)
        self.energy_fig = None
        self.rmsd_heatmap_fig = None
        self.cluster_fig = None
        self.cluster_scatter_palette = None
        self.energy_data = None
        self.rmsd_mat = None
        # select receptor
        self.ui_dlg = None
        self.sele_dlg = None
        # cluster
        self.cluster_n = 2
        self.max_cluster_n = 7
        
    def ui_update_ui(self):
        self.ui_dlg.set_options(list(self._app.lazy_dlg.pose_page.dlg_pose.keys()))
        
    @ui.refreshable
    def plot_energy_curve(self):
        # calculate energy curve
        if self.sele_dlg is None:
            return 
        dlg_pose: DlgFile = self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg]
        idx_energy = [[pose.run_idx, pose.energy] for pose in dlg_pose.pose_lst]
        sorted_energies_by_idx = list(map(lambda x: x[1], sorted(idx_energy, key=lambda x: x[0])))
        sorted_energies = sorted(sorted_energies_by_idx, reverse=True)
        pooled_energies_by_idx = [min(sorted_energies_by_idx[:i]) for i in range(1, len(sorted_energies_by_idx)+1)]
        self.energy_data = (sorted_energies_by_idx, sorted_energies, pooled_energies_by_idx)
        ui.notify(f'Energy curve calculated for {len(dlg_pose.pose_lst)} poses, [{min(sorted_energies_by_idx)} ~ {max(sorted_energies_by_idx)}]')
        # plot energy curve
        plt.close(self.energy_fig)
        self.energy_fig = None
        if self.energy_data is not None:
            sorted_energies_by_idx, sorted_energies, pooled_energies_by_idx = self.energy_data
            mean_energy = np.mean(sorted_energies)
            with ui.pyplot(figsize=(12, 8), close=False) as fig:
                self.energy_fig = fig.fig
                gs = fig.fig.add_gridspec(1, 2,  width_ratios=(4, 1), left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
                # Create the Axes.
                ax = self.energy_fig.add_subplot(gs[0, 0])
                ax_histy = self.energy_fig.add_subplot(gs[0, 1], sharey=ax)
                # Draw the scatter plot and marginals.
                ax.plot(pooled_energies_by_idx, linewidth=2, label='Minimum Energy in Run Order')
                ax.plot([mean_energy] * len(pooled_energies_by_idx), c='black', linewidth=2, label=f'Mean Energy: {mean_energy:.4f} kcal/mol')
                ax.scatter(list(range(len(sorted_energies_by_idx))), sorted_energies_by_idx,
                           alpha=0.4, c='green', s=50, label='Individual Energy in Run Order')
                ax.scatter(list(range(len(sorted_energies))), sorted_energies,
                           alpha=0.2, c='red', s=30, label='Individual Energy in Descending Order')
                ax.set_xlabel('Pose Index', fontdict={'size': 14})
                ax.set_ylabel('Energy (kcal/mol)', fontdict={'size': 14})
                ax.set_title(f'Energy Curve for {self.sele_dlg}', fontdict={'size': 16})
                # plot hist
                ax_histy.tick_params(axis="y", labelleft=False)
                ax_histy.hist(sorted_energies, bins=int(max(sorted_energies)-min(sorted_energies)), orientation='horizontal')
                # minor works
                ax.tick_params(axis='both', which='major', labelsize=12)
                ax_histy.tick_params(axis='both', which='minor', labelsize=12)
                ax.legend(fontsize=12)
                plt.tight_layout()
 
    async def calcu_RMSD(self):
        if self.sele_dlg is None:
            return ui.notify('Please select a receptor and a selection first')
        dlg_pose: DlgFile = self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg]
        self.rmsd_mat = dlg_pose.rmsd('numpy')
        ui.notify(f'RMSD mat calculated for {len(dlg_pose.pose_lst)} poses')

    @ui.refreshable
    def plot_cluster_curve(self):
        if self.sele_dlg is None or self.rmsd_mat is None:
            return 
        dlg_pose: DlgFile = self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg]
        rs = []
        for k in range(2, self.max_cluster_n+1):
            groups_idx = dlg_pose.rmsd_cluster(k)
            sse, ssr, r = dlg_pose.calcu_SSE_SSR(self.rmsd_mat, groups_idx)
            rs.append(r)
        plt.close(self.cluster_fig)
        with ui.pyplot(close=False, figsize=(8, 6)) as fig:
            self.cluster_fig = fig.fig
            fig.fig.gca().plot(list(range(2, self.max_cluster_n+1)), rs)
        fig.update()

    @ui.refreshable
    def plot_cluster_2d_scatter(self):
        if self.sele_dlg is None or self.rmsd_mat is None:
            return 
        dlg_pose: DlgFile = self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg]
        groups_idx = dlg_pose.rmsd_cluster(int(self.cluster_n))
        dots = pca(self.rmsd_mat, 2)
        plt.close(self.cluster_fig)
        with ui.pyplot(close=False, figsize=(10, 8)) as fig:
            ax = fig.fig.gca()
            self.cluster_fig = fig.fig
            for i, c in zip(np.unique(groups_idx), get_palette(len(np.unique(groups_idx)), self.cluster_scatter_palette)):
                ax.scatter(dots[groups_idx==i, 0], dots[groups_idx==i, 1], c=c, alpha=0.3, label=f'Cluster {i}')
            ax.legend()
            ax.set_xlabel('PC1', fontdict={'size': 14})
            ax.set_ylabel('PC2', fontdict={'size': 14})
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_title(f'PCA of RMSD for {self.sele_dlg}', fontdict={'size': 16})
        fig.update()
        
    @ui.refreshable
    def plot_cluster_3d_scatter(self):
        if self.sele_dlg is None or self.rmsd_mat is None:
            return 
        dlg_pose: DlgFile = self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg]
        groups_idx = dlg_pose.rmsd_cluster(int(self.cluster_n))
        dots = pca(self.rmsd_mat, 3)
        plt.close(self.cluster_fig)
        with ui.pyplot(close=False, figsize=(10, 8)) as fig:
            ax = fig.fig.add_subplot(projection='3d')
            self.cluster_fig = fig.fig
            for i, c in zip(np.unique(groups_idx), get_palette(len(np.unique(groups_idx)), self.cluster_scatter_palette)):
                ax.scatter(dots[groups_idx==i, 0], dots[groups_idx==i, 1], dots[groups_idx==i, 2], c=c, alpha=0.3, label=f'Cluster {i}')
            ax.legend()
            ax.set_xlabel('PC1', fontdict={'size': 14})
            ax.set_ylabel('PC2', fontdict={'size': 14})
            ax.set_zlabel('PC3', fontdict={'size': 14})
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.set_title(f'PCA of RMSD for {self.sele_dlg}', fontdict={'size': 16})
        plt.show()
            
    def build_gui(self):
        with ui.splitter(value = 20).classes('w-full h-full') as splitter:
            with splitter.before:
                # choose ligand
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('select a ligand').classes('w-full')
                        self.ui_dlg = ui.select(self._app.lazy_dlg.pose_page.now_dlg_name or [],
                                                label = 'select a dlg').bind_value_to(self, 'sele_dlg').classes('w-full').props('use-chips')
            # vitualization
            with splitter.after:
                with ui.tabs().classes('w-full').props('align=left active-bg-color=blue') as tabs:
                    energy_tab = ui.tab('Energy Curve').props('no-caps')
                    cluster_tab = ui.tab('RMSD Cluster').props('no-caps')
                with ui.tab_panels(tabs, value=energy_tab).classes('w-full'):
                    with ui.tab_panel(energy_tab):
                        ui.button('plot energy curve', on_click=self.plot_energy_curve.refresh).props('no-caps')
                        self.plot_energy_curve()
                    with ui.tab_panel(cluster_tab):
                        with ui.row().classes('w-full flex flex-grow'):
                            ui.button('calculate RMSD', on_click=self.calcu_RMSD).props('no-caps')
                            ui.number(label='cluster number', value=self.cluster_n, min=2, max=10, step=1).classes('flex flex-grow').bind_value_to(self, 'cluster_n')
                            ui.number(label='max cluster number', value=self.max_cluster_n, min=2, max=10, step=1).classes('flex flex-grow').bind_value_to(self,'max_cluster_n')
                            avaliabel_palette = ['hls', 'Set1', 'Set2', 'Set3', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'tab10', 'tab20', 'tab20b', 'tab20c']
                            ui.select(label='palette', value='hls', options=avaliabel_palette).classes('flex flex-grow').bind_value_to(self, 'cluster_scatter_palette')
                        with ui.row().classes('w-full flex flex-grow'):
                            ui.button('plot cluster SSE/SST curve', on_click=self.plot_cluster_curve.refresh).props('no-caps')
                            ui.button('plot cluster 2D scatter', on_click=self.plot_cluster_2d_scatter.refresh).props('no-caps')
                            ui.button('plot cluster 3D scatter', on_click=self.plot_cluster_3d_scatter.refresh).props('no-caps')
                        self.plot_cluster_curve()
                        self.plot_cluster_2d_scatter()
                        self.plot_cluster_3d_scatter()


class InteractionPage(LazyPose):
    def __init__(self, app):
        self._app = app
        self._app.ui_update_func.append(self.ui_update_ui)
        # select receptor
        self.ui_molecule = None
        self.ui_sele = None
        self.ui_dlg = None
        self.sele_molecule = None
        self.sele_selection = None
        self.sele_dlg = None
        # interaction
        self.interaction_method = 'pymol'
        self.support_mode: Dict[str, List[str]] = {
            'pymol': MODE_PML,
            'PLIP': MODE_PLIP,
            'LigPlus': MODE_LIGPLUS,
        }
        self.ui_mode = None
        self.interactions = None
        self.interaction_df = None
        self.interaction_mode = 0
        self.distance_cutoff = 4
        self.nagetive_factor = -0.5
        # vitulazation control
        self.fig = None
        self.align_vlim = False
        self.plot_cluster = False
        self.use_one_letter = True
        self.min_ligand_interaction = 0 # abs
        self.min_receptor_interaction = 0 # abs
        self.fig_w = 12
        self.fig_h = 7
        
    def ui_update_ui(self):
        super().ui_update_ui()
        self.ui_dlg.set_options(list(self._app.lazy_dlg.pose_page.dlg_pose.keys()))
            
    def calculate_interaction(self):
        # get receptor
        if self.sele_molecule is None and self.sele_selection is None:
            ui.notify('Please select a receptor')
            return None
        # get ligand
        if self.sele_dlg is None:
            ui.notify('Please select a DLG')
        ligands = []
        for name in self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg].n2i:
            if self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg].get_pose_prop('is_show', name):
                ligands.append(self._app.lazy_dlg.pose_page.dlg_pose[self.sele_dlg].get_pose_prop('pml_name', name))
        if not ligands:
            ui.notify('Please select at least one ligand')
            return None
        # calculate interactions
        if self.interaction_method == 'pymol':
            calc_fn = calc_pml_interaction
        elif self.interaction_method == 'PLIP':
            calc_fn = calc_plip_interaction
        elif self.interaction_method == 'LigPlus':
            calc_fn = calc_ligplus_interaction
        else:
            assert False, f'Unknown interaction method: {self.interaction_method}'
        self.interactions, self.interaction_df = calc_fn(self.sele_molecule or self.sele_selection, ligands,
                                                         mode=self.interaction_mode, cutoff=self.distance_cutoff,
                                                         nagetive_factor=self.nagetive_factor)
        if self.interactions is not None:
            ui.notify(f'Interactions calculated for {len(ligands)} ligands')
        else:
            ui.notify(f'No interactions found between {self.sele_molecule or self.sele_selection} and {len(ligands)} ligands')
        return self.interactions
            
    @ui.refreshable
    def plot_interaction(self):
        plt.close(self.fig)
        self.fig = None
        if (self.interaction_df is not None) and (not self.interaction_df.empty):
            
            tmp_interaction_df = filter_interaction_df(self.interaction_df, self.min_receptor_interaction,
                                                       self.min_ligand_interaction, inplace=False)
            # plot
            with ui.pyplot(figsize=(self.fig_w, self.fig_h), close=False) as fig:
                vmax, vmin = tmp_interaction_df.max().max(), tmp_interaction_df.min().min()
                if self.align_vlim:
                    vlim = max(abs(vmax), abs(vmin))
                    vmax, vmin = vlim, -vlim
                if self.plot_cluster:
                    grid = sns.clustermap(tmp_interaction_df, cmap='coolwarm', figsize=(self.fig_w, self.fig_h),
                                        annot=True, fmt='.2f', vmax=vmax, vmin=vmin,
                                        linewidths=0.5, linecolor='black', cbar_kws={"shrink": 0.5})
                    self.fig = fig.fig = grid._figure
                else:
                    self.fig = fig.fig
                    grid = sns.heatmap(tmp_interaction_df, cmap='coolwarm', ax = fig.fig.gca(),
                                    annot=True, fmt='.2f', vmax=vmax, vmin=vmin,
                                    linewidths=0.5, linecolor='black', cbar_kws={"shrink": 0.5})
                plt.tight_layout()
                
    def _handle_method_change(self, e):
        self.ui_mode.set_options(self.support_mode[self.interaction_method])
        # set multi-select
        if self.interaction_method in {'LigPlus', 'PLIP'}:
            self.ui_mode.multiple = True
            self.ui_mode.props['multiple'] = True
            self.ui_mode.props['use-chips'] = True
        else:
            self.ui_mode.multiple = False
            self.ui_mode.props['multiple'] = False
            self.ui_mode.props['use-chips'] = False
        
    def build_gui(self):
        with ui.splitter(value = 20).classes('w-full h-full') as splitter:
            with splitter.before:
                # choose receptor
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('select a receptor')
                        self.ui_molecule = ui.select(self._app._now_molecule,
                                                    label='select a molecule').bind_value_to(self, 'sele_molecule').classes('w-full').props('use-chips')
                        self.ui_sele = ui.select(self._app._now_selection,
                                                label='select a selection').bind_value_to(self, 'sele_selection').classes('w-full').props('use-chips')
                # choose ligand
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('select a ligand').classes('w-full')
                        self.ui_dlg = ui.select(self._app.lazy_dlg.pose_page.now_dlg_name or [],
                                                label='select a dlg').bind_value_to(self, 'sele_dlg').classes('w-full').props('use-chips')
                # interaction calculation config
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('interaction calculation config').classes('w-full')
                        ui.select(['pymol', 'PLIP', 'LigPlus'], label='method', value=self.interaction_method).classes('w-full').bind_value_to(self, 'interaction_method').on_value_change(self._handle_method_change)
                        ui.number('distance cutoff', value=self.distance_cutoff).classes('w-full').bind_value_to(self, 'distance_cutoff')
                        self.ui_mode = ui.select(options=self.support_mode[self.interaction_method], label='mode').classes('w-full').bind_value_to(self, 'interaction_mode')
                        ui.number('nagetive factor', value=self.nagetive_factor).classes('w-full').bind_value_to(self, 'nagetive_factor')
                # plot config
                with ui.card().classes('w-full'):
                    with ui.column().classes('w-full'):
                        ui.label('plot config').classes('w-full')
                        ui.checkbox('plot cluster', value=self.plot_cluster).bind_value_to(self, 'plot_cluster')
                        ui.checkbox('align vlim', value=self.align_vlim).bind_value_to(self, 'align_vlim')
                        ui.checkbox('use one letter', value=self.use_one_letter).bind_value_to(self, 'use_one_letter')
                        ui.number('min ligand interaction (abs)', value=self.min_ligand_interaction).bind_value_to(self,'min_ligand_interaction')
                        ui.number('min receptor interaction (abs)', value=self.min_receptor_interaction).bind_value_to(self,'min_receptor_interaction')
                        ui.number('fig width', value=self.fig_w).bind_value_to(self, 'fig_w')
                        ui.number('fig height', value=self.fig_h).bind_value_to(self, 'fig_h')
            # interaction vitualization
            with splitter.after:
                with ui.row():
                    ui.label('Interaction: ')
                    ui.button('calculate', on_click=self.calculate_interaction).classes('flex flex-grow').props('no-caps')
                    ui.button('save calcu result', on_click=self.save_calcu_result).classes('flex flex-grow').props('no-caps')
                    ui.button('plot', on_click=self.plot_interaction.refresh).classes('flex flex-grow').props('no-caps')
                    ui.button('save Plot', on_click=self.save_plot).classes('flex flex-grow').props('no-caps')
                self.plot_interaction()
        
    def save_calcu_result(self):
        if self.interaction_df is None:
            ui.notify('Please calculate the interaction first')
            return None
        data_path = MyFileDialog(types = [('LazyInteraction file', '*.pkl')],
                                 initialdir=os.getcwd()).get_save_file()
        if data_path is None:
            ui.notify('Please select a file to save')
            return None
        if not data_path.endswith('.pkl'):
            data_path += '.pkl'
        with open(data_path, 'wb') as f:
            f.write(gzip.compress(pickle.dumps(self.interactions)))            
            self.interaction_df.to_pickle(f, compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1})
        ui.notify(f'Interactions saved to {data_path}')
    
    def save_plot(self):
        if self.fig is None:
            ui.notify('Please plot the interaction first')
            return None
        fig_path = MyFileDialog(initialdir=os.getcwd()).get_save_file()
        if fig_path is None:
            ui.notify('Please select a file to save')
            return None
        save_show(fig_path, dpi = 600, show = False)
    

class LazyDLG:
    def __init__(self, app, _dev_mode: bool = False):
        self._app = app
        
        self.pose_page = LazyPose(self._app)
        self.result_page = ResultPage(self._app)
        self.analysis_page = InteractionPage(self._app)
        
    def build_gui(self):
        with ui.tabs().classes('w-full').props('align=left active-bg-color=blue') as tabs:
            self.ui_loader_tab = ui.tab('DLG Pose Loader').props('no-caps')
            self.ui_result_tab = ui.tab('DLG Pose Result').props('no-caps')
            self.ui_analysis_tab = ui.tab('DLG Pose Analysis').props('no-caps')
        with ui.tab_panels(tabs, value=self.ui_loader_tab).classes('w-full'):
            with ui.tab_panel(self.ui_loader_tab):
                self.pose_page.build_gui()
            with ui.tab_panel(self.ui_result_tab):
                self.result_page.build_gui()
            with ui.tab_panel(self.ui_analysis_tab):
                self.analysis_page.build_gui()
        # return self
        return self
        
        
# dev mode
if __name__ in {"__main__", "__mp_main__"}:
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    
    from main import GUILauncher
    GUILauncher(None, 'process', 4)
    