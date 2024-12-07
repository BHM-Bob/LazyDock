from __future__ import absolute_import, print_function

import re
import time
import tkinter as tk
import tkinter.filedialog as tkFileDialog
from pathlib import Path
from typing import Any, Callable, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from mbapy_lite.base import put_err
from mbapy_lite.file import decode_bits_to_str, opts_file
from mbapy_lite.game import BaseInfo
from mbapy_lite.stats.cluster import KMeans, KMeansBackend
from mbapy_lite.web_utils.task import TaskPool
from pymol import cmd
from tqdm import tqdm

from lazydock.utils import uuid4

# atom-type, atom-number, atom-name, residue-name, chain-name, residue-number, x, y, z, occupancy, temperature-factor
# ATOM      1  CA  LYS     7     136.747 133.408 135.880 -0.06 +0.10 OA
# HETATM   79  CB  PRO     9     132.763 138.092 170.987+13.19 -0.02    +0.037 C  \n
# HETATM   40  H   UNK     1     -25.811   4.632  -7.417  1.00 20.00     0.163 HD
PDB_PATTERN = r"(ATOM|HETATM) +(\d+) +(\w+) +(\w+) +(\w+)? +(\d+) +([\d\-\.]+) +([\d\-\.]+) +([\d\-\.]+)([ \+\-]+[\d\-\.]+)([ \+\-]+[\d\-\.]+) [ \.\-\+\d]+ ([A-Z]+)"
PDB_FORMAT = "{:6s}{:>5s}  {:<3s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}{:>6s}{:>6s}          {:>2s}  "
PDB_FORMAT2= "{:6s}{:>5s} {:<4s} {:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}{:>6s}{:>6s}          {:>2s}  "


def calcu_RMSD(pose1: Union[str, 'ADModel'], pose2: Union[str, 'ADModel'], _cmd = None):
    _cmd = _cmd or cmd
    def _get_pose(pose):
        if isinstance(pose, ADModel):
            name = uuid4()
            _cmd.read_pdbstr(pose.as_pdb_string(), name)
            return name, True
        return pose, False
    pml_pose1, is_loaded1 = _get_pose(pose1)
    pml_pose2, is_loaded2 = _get_pose(pose2)
    rmsd = _cmd.rms(pml_pose1, pml_pose2, cycles=0) # cutoff = float: outlier rejection cutoff (only if cycles>0) {default: 2.0}
    if is_loaded1:
        _cmd.delete(pml_pose1)
    if is_loaded2:
        _cmd.delete(pml_pose2)
    return rmsd


class ADModel(BaseInfo):
    """STORAGE CLASS FOR DOCKED LIGAND"""
    def __init__(self, content: str = None, _sort_atom_by_res: bool = False,
                 _parse2std: bool = False, default_chain: str = 'Z'):
        self.energy = 0.
        self.name = ''
        self.run_idx = 0
        self.info = ''
        self.pdb_string = ''
        self.default_chain = default_chain
        self.pdb_lines = []
        self.coords = None
        self.coords_only_bb = None
        if content is not None:
            self.parse_content(content, _sort_atom_by_res, _parse2std)

    def parse_content(self, content: str, _sort_atom_by_res: bool = False,
                      _parse2std: bool = False):
        # parse pdb lines
        self.info = content
        self.pdb_lines = list(map(lambda x: x[0], re.findall(r'((ATOM|HETATM).+?\n)', self.info)))
        self.pdb_atoms = list(map(list, re.findall(PDB_PATTERN, self.info)))
        if len(self.pdb_atoms) != len(self.pdb_lines):
            put_err(f'pdb_atoms and pdb_lines length not match: {len(self.pdb_atoms)} vs {len(self.pdb_lines)}, just use pdb_lines')
            _sort_atom_by_res = _parse2std = False
        if _sort_atom_by_res:
            pack = sorted(zip(self.pdb_lines, self.pdb_atoms), key = lambda x : (x[1][4], int(x[1][5]), int(x[1][1])))
            self.pdb_lines, self.pdb_atoms = zip(*pack)
        if _parse2std:
            for atom in self.pdb_atoms:
                if not atom[4]:
                    atom[4] = self.default_chain
                if len(atom[2]) == 4:
                    atom[-1] = atom[2][1] # such as 1HH1
                elif len(atom[2]) >= 1:
                    if atom[2][0].isdigit():
                        atom[-1] = atom[2][1] # such as 1HH
                    else:
                        atom[-1] = atom[2][0] # such as H, CA, NH1
            self.pdb_lines = [PDB_FORMAT2.format(*line) if len(line[2])==4 else PDB_FORMAT.format(*line) for line in self.pdb_atoms]
            self.pdb_string = '\n'.join(self.pdb_lines)
        else:
            self.pdb_string = ''.join(self.pdb_lines)
        # parse energy could be str.find?
        energy_line = re.findall(r'USER.+?Free Energy of Binding.+?\n', self.info)
        if energy_line:
            entr = energy_line[0].split('=')[1]
            self.energy = float(entr.split()[0])
        else:
            energy_line = re.findall(r'REMARK.+?VINA RESULT.+?\n', self.info)
            if energy_line:
                entr = energy_line[0].split(':')[1]
                self.energy = float(entr.split()[0])
            else:
                self.energy = None
        # parse run idx
        self.run_idx = int(re.findall(r'MODEL +(\d+)', self.info)[0])

    def as_pdb_string(self):
        return self.pdb_string

    def info_string(self):
        return self.info
    
    def parse_coords(self, only_bb: bool = True):
        if self.coords is None or only_bb != self.coords_only_bb:
            self.coords_only_bb = only_bb
            self.coords = np.array([list(map(float, line[6:9])) for line in self.pdb_atoms if line[2] in ['C', 'N', 'O', 'CA'] or not only_bb])
        if self.coords.shape[0] == 0:
            # considerate that atom name may not be 'C', 'N', 'O', 'CA', try to use last column
            self.coords = np.array([list(map(float, line[6:9])) for line in self.pdb_atoms if line[-1] in ['C', 'N', 'O', 'CA'] or not only_bb])
        return self.coords
    
    def calcu_rmsd_by_coords(self, other: 'ADModel', only_bb: bool = True):
        self_coords = self.parse_coords(only_bb)
        other_coords = other.parse_coords(only_bb)
        d2 = np.sum((self_coords - other_coords)**2, axis=1)
        return np.sqrt(np.mean(d2))
    
    def rmsd(self, other: 'ADModel', backend: str = 'pymol', only_bb: bool = True):
        if backend == 'pymol':
            return calcu_RMSD(self, other)
        elif backend == 'rdkit':
            from rdkit import Chem
            from rdkit.Chem import AllChem
            pose1 = Chem.MolFromPDBBlock(self.as_pdb_string(), removeHs=True, sanitize=False)
            pose2 = Chem.MolFromPDBBlock(other.as_pdb_string(), removeHs=True, sanitize=False)
            return AllChem.GetBestRMS(pose1, pose2)
        elif backend == 'numpy':
            return self.calcu_rmsd_by_coords(other, only_bb)
        else:
            raise ValueError(f'Unsupported backend: {backend}')


class RmsdClusterBackend(KMeansBackend):
    def __init__(self) -> None:
        super().__init__('scipy')
    def cdist(self, data, centers):
        return np.sqrt(np.mean(np.sum((data[:, None] - centers[None, :])**2, axis=-1), axis=-1))


class DlgFile(BaseInfo):
    def __init__(self, path: str = None, content: str = None,
                 sort_pdb_line_by_res: bool = False,
                 parse2std: bool = False):
        # load content from file if path is provided
        self.path = path
        self.sort_pdb_line_by_res = sort_pdb_line_by_res
        self.parse2std = parse2std
        if path is not None:
            self.content = decode_bits_to_str(opts_file(path, 'rb'))
        elif content is not None:
            self.content = content
        else:
            self.content = None
        # decode content to pose_lst
        if self.content is not None:
            self.pose_lst: List[ADModel] = self.decode_content()
            self.sort_pose()
        else:
            self.pose_lst: List[ADModel] = []
        self.n2i = {}
        
    def __len__(self):
        return len(self.pose_lst)
    
    def merge(self, other: 'DlgFile', inplace: bool = True):
        """only merge poses, not content or other attributes"""
        if not inplace:
            new_dlg = DlgFile()
            new_dlg.pose_lst = self.pose_lst + other.pose_lst
            return new_dlg
        else:
            self.pose_lst += other.pose_lst
            return self
        
    def sort_pose(self, key: Callable[[ADModel], Any] = None,
                  inplace: bool = True, reverse: bool = False) -> List[ADModel]:
        if key is None:
            key = lambda x : x.energy
        if inplace:
            self.pose_lst.sort(key=key, reverse=reverse)
            ret = self.pose_lst
        else:
            ret = sorted(self.pose_lst, key=key, reverse=reverse)
        return ret
        
    def decode_content(self):
        dlg_pose_lst = []
        for model in re.findall('MODEL.+?ENDMDL', self.content, re.DOTALL):
            model = model.replace('\nDOCKED: ', '\n')
            dlg_pose_lst.append(ADModel(model, self.sort_pdb_line_by_res, self.parse2std))
        return dlg_pose_lst
    
    def asign_pose_name(self, pose_names: List[str]):
        if len(pose_names) != len(self.pose_lst):
            raise ValueError("Number of pose names must match number of poses")
        for i, name in enumerate(pose_names):
            self.n2i[name] = i
        
    def asign_prop(self, prop: str, value: List[Any]):
        if value is not None and len(value) == len(self.pose_lst):
            setattr(self, prop, value)
                
    def set_pose_prop(self, prop: str, value: Any, pose_name: str = None, pose_idx: int = None):
        if pose_name is None and pose_idx is None:
            raise ValueError("Either pose_name or pose_idx must be provided")
        if pose_name is not None:
            pose_idx = self.n2i[pose_name]
        if not hasattr(self, prop):
            setattr(self, prop, [None]*len(self.pose_lst))
        getattr(self, prop)[pose_idx] = value
            
    def get_pose(self, pose_name: str = None, pose_idx: int = None):
        if pose_name is None and pose_idx is None:
            raise ValueError("Either pose_name or pose_idx must be provided")
        if pose_name is not None:
            pose_idx = self.n2i[pose_name]
        return self.pose_lst[pose_idx]
    
    def get_pose_prop(self, prop: str, pose_name: str = None, pose_idx: int = None, default: Any = None):
        if pose_name is None and pose_idx is None:
            raise ValueError("Either pose_name or pose_idx must be provided")
        if pose_name is not None:
            pose_idx = self.n2i[pose_name]
        if hasattr(self, prop):
            return getattr(self, prop)[pose_idx]
        else:
            return default
    
    def parse_energies(self):
        """
        calculate energies of all poses, and sort them by idx and energy
        
        - self.sorted_energies_by_idx: [N_pose], sorted energies by pose run index
        - self.sorted_energies: [N_pose], sorted energies in descending order, first is the biggest
        - self.pooled_energies_by_idx: [N_pose], pooled energies by pose index
        - self.mean_energy: float, mean energy of all poses
        """
        idx_energy = [[pose.run_idx, pose.energy] for pose in self.pose_lst]
        self.sorted_energies_by_idx = list(map(lambda x: x[1], sorted(idx_energy, key=lambda x: x[0])))
        self.sorted_energies = sorted(self.sorted_energies_by_idx, reverse=True)
        self.pooled_energies_by_idx = [min(self.sorted_energies_by_idx[:i]) for i in range(1, len(self.sorted_energies_by_idx)+1)]
        self.mean_energy = np.mean(self.sorted_energies)
        
    def plot_energy_curve(self, figsize: Tuple[int, int] = (8, 6)):
        """
        plot energy curve of all pose
        - plot pooled energy curve
        - plot individual energy curve
        - plot hist of individual energy
        - plot mean energy
        - plot best energy
        
        Returns:
            - fig, ax, ax_histy
        """
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(1, 2,  width_ratios=(4, 1), left=0.1, right=0.9, bottom=0.1, top=0.9,
                wspace=0.05, hspace=0.05)
        # Create the Axes.
        ax = fig.add_subplot(gs[0, 0])
        ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
        # Draw the scatter plot and marginals.
        ax.plot(self.pooled_energies_by_idx, linewidth=2, label='Minimum Energy in Run Order')
        ax.plot([self.mean_energy] * len(self.pooled_energies_by_idx), c='black', linewidth=2, label=f'Mean Energy: {self.mean_energy:.4f} kcal/mol')
        ax.plot([self.sorted_energies[-1]] * len(self.pooled_energies_by_idx), c='gray', linestyle='--',
                linewidth=2, label=f'Best Energy: {self.sorted_energies[-1]:.4f} kcal/mol')
        ax.scatter(list(range(len(self.sorted_energies_by_idx))), self.sorted_energies_by_idx,
                    alpha=0.4, c='green', s=50, label='Individual Energy in Run Order')
        ax.scatter(list(range(len(self.sorted_energies))), self.sorted_energies,
                    alpha=0.2, c='red', s=30, label='Individual Energy in Descending Order')
        ax.set_xlabel('Pose Index', fontdict={'size': 14})
        ax.set_ylabel('Energy (kcal/mol)', fontdict={'size': 14})
        ax.set_title(f'Energy Curve for {Path(self.path).stem}', fontdict={'size': 16})
        # plot hist
        ax_histy.tick_params(axis="y", labelleft=False)
        ax_histy.hist(self.sorted_energies, bins=int(max(self.sorted_energies)-min(self.sorted_energies)), orientation='horizontal')
        # minor works
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax_histy.tick_params(axis='both', which='minor', labelsize=12)
        ax.legend(fontsize=12, draggable=True, framealpha=0.8)
        plt.tight_layout()
        return fig, ax, ax_histy
        
    def rmsd(self, backend: str = 'pymol', taskpool: TaskPool = None, verbose: bool = False):
        # if backend is numpy, use matrix calculation to calculate rmsd for acceleration
        if backend == 'numpy':
            # [N_pose, M_atom, D_3]
            coords = np.concatenate([[pose.parse_coords()] for pose in self.pose_lst], axis=0)
            # d2 = [N_pose, N_pose, M_atom, D_3] => [N_pose, N_pose, M_atom]: np.sum((coords[None, :, :] - coords[:, None, :])**2, axis=-1)
            # sum+devide => mean => [N_pose, N_pose, M_atom] => [N_pose, N_pose]
            return np.sqrt(np.mean(np.sum((coords[None, :, :] - coords[:, None, :])**2, axis=-1), axis=-1))
        # if is other backend, loop over all pairs and calculate rmsd
        rmsd_mat = np.zeros((len(self.pose_lst), len(self.pose_lst))).tolist()
        for i, pose_i in tqdm(enumerate(self.pose_lst), desc='Calculating RMSD', total=len(self.pose_lst), disable=not verbose):
            for j, pose_j in enumerate(self.pose_lst):
                if i > j:
                    continue
                if taskpool is None:
                    rmsd_mat[i][j] = pose_i.rmsd(pose_j, backend=backend)
                else:
                    while taskpool.count_waiting_tasks() > 0:
                        time.sleep(0.1)
                    rmsd_mat[i][j] = taskpool.add_task(f'{i}-{j}', pose_i.rmsd, pose_j, backend=backend)
        for i in range(len(rmsd_mat)):
            for j in range(len(rmsd_mat)):
                if i > j:
                    rmsd_mat[i][j] = rmsd_mat[j][i]
                elif taskpool is not None:
                    rmsd_mat[i][j] = taskpool.query_task(f'{i}-{j}', block=True, timeout=999)
        return rmsd_mat
    
    def rmsd_cluster(self, n_clusters: int):
        coords = np.concatenate([[pose.parse_coords()] for pose in self.pose_lst], axis=0)
        cluster = KMeans(n_clusters, backend=RmsdClusterBackend())
        return cluster.fit_predict(coords)
    
    def calcu_SSE_SSR(self, rmsd_mat: np.ndarray, groups_idx: np.ndarray,
                      centers_idx: Optional[np.ndarray] = None, center_idx: Optional[int] = None):
        """
        Parameters
            - rmsd_mat: [N_pose, N_pose], RMSD matrix return by DlgFile.rmsd
            - groups_idx: [N_pose], clustering result return by DlgFile.rmsd_cluster or other clustering algorithm
            - centers_idx: [N_cluster], clustering center index, **must in the rmsd matrix**, if not provided, use the minimum of each group
            - center_idx: int, center of the all poses, **must in the rmsd matrix**, if not provided, use the minimum of all poses
            
        Returns
            - SSE: 组内误差平方和
            - SSR: 组间误差平方和
            - SSR/SST: 
        """
        if centers_idx is None:
            centers_idx = np.array([np.argmin(rmsd_mat[groups_idx==i]) for i in np.unique(groups_idx)])
        if center_idx is None:
            center_idx = np.argmin(rmsd_mat)
        sse, ssr = 0, 0
        for i in np.unique(groups_idx):
            sse += np.sum(rmsd_mat[groups_idx==i, centers_idx[i]]**2)
            ssr += np.sum(groups_idx==i) * rmsd_mat[centers_idx[i], center_idx]**2
        return sse, ssr, ssr/(sse+ssr)
    
    def get_top_k_pose(self, method: str = 'energy', k: int = 1, rmsd_mat: np.ndarray = None,
                       groups_idx: np.ndarray = None) -> List[ADModel]:
        if method == 'energy':
            return self.sort_pose(key=lambda x: x.energy, inplace=False)[:k]
        elif method == 'rmsd':
            if rmsd_mat is None:
                rmsd_mat = self.rmsd('numpy')
            if groups_idx is None:
                raise ValueError("groups_idx must be provided if method is 'rmsd'")
            rmsd_vec = rmsd_mat.sum(axis=0)
            idx_vec = np.arange(rmsd_vec.shape[0])
            top_poses = {}
            for group_idx in np.unique(groups_idx):
                sub_idx, sub_rmsd_vec = idx_vec[groups_idx==group_idx], rmsd_vec[groups_idx==group_idx]
                idx = sub_idx[np.argsort(sub_rmsd_vec)[:k]] # increasing order
                top_poses[group_idx] = [self.pose_lst[i] for i in idx]
            return top_poses
        else:
            raise ValueError(f'Unsupported method: {method}')


def merge_dlg(dlg_lst: List[DlgFile]) -> DlgFile:
    if not dlg_lst:
        return put_err('Empty input list, return None')
    new_dlg = DlgFile()
    for dlg in dlg_lst:
        new_dlg.merge(dlg, inplace=True)
    return new_dlg


def tk_file_dialog_wrapper(*args, **kwargs):
    def ret_wrapper(tk_file_dialog_func):
        def core_wrapper(*args, **kwargs):
            parent = tk.Tk()
            result = tk_file_dialog_func(*args, parent = parent, **kwargs)
            parent.withdraw()
            if result == "":
                return None
            else:
                return result
        return core_wrapper
    return ret_wrapper
            

class MyFileDialog:
    def __init__(self, types=[("Executable", "*")], initialdir: str = None):
        self.initialdir = initialdir
        self.types = types

    @tk_file_dialog_wrapper()
    def get_open_file(self, parent):
        return tkFileDialog.askopenfilename(parent=parent, initialdir = self.initialdir, filetypes=self.types)


    @tk_file_dialog_wrapper()
    def get_save_file(self, parent):
        return tkFileDialog.asksaveasfilename(parent=parent, initialdir = self.initialdir, filetypes=self.types)
        
    @tk_file_dialog_wrapper()
    def get_ask_dir(self, parent):
        return tkFileDialog.askdirectory(parent=parent, initialdir = self.initialdir)


if __name__ == '__main__':
    # from mbapy_lite.base import TimeCosts
    # @TimeCosts(10)
    # def load_test(idx):
    #     DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=True)
    # load_test()
    normal = DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=False, parse2std=False)
    std = DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=False, parse2std=True)
    std.parse_energies()
    std.plot_energy_curve()
    plt.show()
    from pymol import cmd
    from rdkit import Chem
    cmd.reinitialize()
    cmd.read_pdbstr(std.pose_lst[0].as_pdb_string(), 'std')
    print(std.pose_lst[0].as_pdb_string(), '\n\n')
    print(cmd.get_pdbstr('std'))
    ligand = Chem.MolFromPDBBlock(cmd.get_pdbstr('std'), removeHs=True, sanitize=False)
    
    import seaborn as sns
    from mbapy_lite.stats import pca
    dlg = DlgFile(path='data_tmp/dlg/1000run.dlg', sort_pdb_line_by_res=True, parse2std=True)
    dlg.sort_pose(inplace=True)
    energies = np.array([pose.energy for pose in dlg.pose_lst])
    rmsd = dlg.rmsd('numpy')
    
    sns.clustermap(rmsd)
    plt.show()
    
    np.random.seed(3407)
    rs = []
    for k in range(2, 7):
        groups_idx = dlg.rmsd_cluster(k)
        sse, ssr, r = dlg.calcu_SSE_SSR(rmsd, groups_idx)
        rs.append(r)
    plt.plot(list(range(2, 7)), rs)
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    group_cmaps = ['Reds', 'Blues', 'Greens']
    groups_idx = dlg.rmsd_cluster(len(group_cmaps))
    top_poses = dlg.get_top_k_pose(method='rmsd', k=1, rmsd_mat=rmsd, groups_idx=groups_idx)
    dots = pca(rmsd, 3)
    for i, group_cmap in enumerate(group_cmaps):
        ax.scatter(dots[groups_idx==i, 0], dots[groups_idx==i, 1], dots[groups_idx==i, 2],
                   c=energies[groups_idx==i], cmap=group_cmap, alpha=0.4, s=40)
    plt.show()
    