import os
import tempfile
import warnings
from typing import List, Union

import dgl
import networkx as nx
import numpy as np
import pymol
import torch
from rdkit import Chem, RDLogger
from scipy.spatial import distance_matrix

RDLogger.DisableLog('rdApp.*')
np.set_printoptions(threshold=np.inf)
warnings.filterwarnings('ignore')

if __name__ == '__main__':
    from lazydock.nn.EHIGN.utils import angle, area_triangle, cal_dist
else:
    from .utils import angle, area_triangle, cal_dist


def _mols_to_graph(ligand: Chem.Mol, pocket: Chem.Mol, dis_threshold: float = 5.0) -> dgl.DGLHeteroGraph:
    """
    从配体和口袋分子对象构建DGL异构图。

    Parameters
    ----------
    ligand : Chem.Mol
        配体分子对象
    pocket : Chem.Mol
        口袋分子对象
    dis_threshold : float, default 5.0
        距离阈值（埃）

    Returns
    -------
    dgl.DGLHeteroGraph
        DGL异构图
    """
    atom_num_l = ligand.GetNumAtoms()
    atom_num_p = pocket.GetNumAtoms()

    x_l, edge_index_l, edge_attr_l = mol2graph(ligand)
    x_p, edge_index_p, edge_attr_p = mol2graph(pocket)
    (edge_index_l2p, edge_attr_l2p), (edge_index_p2l, edge_attr_p2l) = inter_graph(
        ligand, pocket, dis_threshold=dis_threshold)

    graph_data = {
        ('ligand', 'intra_l', 'ligand'): (edge_index_l[0], edge_index_l[1]),
        ('pocket', 'intra_p', 'pocket'): (edge_index_p[0], edge_index_p[1]),
        ('ligand', 'inter_l2p', 'pocket'): (edge_index_l2p[0], edge_index_l2p[1]),
        ('pocket', 'inter_p2l', 'ligand'): (edge_index_p2l[0], edge_index_p2l[1])
    }
    g = dgl.heterograph(graph_data, num_nodes_dict={"ligand": atom_num_l, "pocket": atom_num_p})
    g.nodes['ligand'].data['h'] = x_l
    g.nodes['pocket'].data['h'] = x_p
    g.edges['intra_l'].data['e'] = edge_attr_l
    g.edges['intra_p'].data['e'] = edge_attr_p
    g.edges['inter_l2p'].data['e'] = edge_attr_l2p
    g.edges['inter_p2l'].data['e'] = edge_attr_p2l

    return g


def _generate_pocket_from_receptor(receptor_pdbstr: str, ligand_pdbstr: str,
                                    distance: float = 5.0, remove_hs: bool = True) -> Union[str, None]:
    """
    使用pymol从受体PDB字符串生成口袋PDB字符串。

    Parameters
    ----------
    receptor_pdbstr : str
        受体的PDB格式字符串
    ligand_pdbstr : str
        配体的PDB格式字符串
    distance : float, default 5.0
        距离阈值（埃）
    remove_hs : bool, default True
        是否移除氢原子

    Returns
    -------
    str or None
        口袋的PDB字符串，失败时返回None
    """
    pymol.cmd.reinitialize()
    pymol.cmd.read_pdbstr(receptor_pdbstr, 'rec')
    pymol.cmd.remove('resn HOH')
    pymol.cmd.remove('inorganic')
    pymol.cmd.read_pdbstr(ligand_pdbstr, 'lig')
    if remove_hs:
        pymol.cmd.remove('hydrogens')
    pymol.cmd.select('Pocket', f'byres lig around {distance}')
    pocket_pdbstr = pymol.cmd.get_pdbstr('Pocket')
    pymol.cmd.delete('all')

    if pocket_pdbstr.strip() == '':
        return None
    return pocket_pdbstr


def pdb_path_to_graph(ligand_path: str, pocket_path: str = None,
                       receptor_path: str = None, dis_threshold: float = 5.0) -> Union[dgl.DGLHeteroGraph, None]:
    """
    从PDB文件路径构建单个DGL异构图。

    Parameters
    ----------
    ligand_path : str
        配体PDB文件路径
    pocket_path : str, optional
        口袋PDB文件路径。如果提供，则忽略receptor参数
    receptor_path : str, optional
        受体PDB文件路径（用于生成口袋）
    dis_threshold : float, default 5.0
        距离阈值（埃）

    Returns
    -------
    dgl.DGLHeteroGraph or None
        DGL异构图，失败时返回None
    """
    ligand = Chem.MolFromPDBFile(ligand_path, removeHs=True)
    if ligand is None:
        print(f"Unable to process ligand: {ligand_path}")
        return None

    if pocket_path is not None:
        pocket = Chem.MolFromPDBFile(pocket_path, removeHs=True)
        if pocket is None:
            print(f"Unable to process pocket: {pocket_path}")
            return None
    elif receptor_path is not None:
        with open(receptor_path, 'r') as f:
            receptor_pdbstr = f.read()
        with open(ligand_path, 'r') as f:
            ligand_pdbstr = f.read()

        pocket_pdbstr = _generate_pocket_from_receptor(
            receptor_pdbstr, ligand_pdbstr, distance=dis_threshold
        )
        if pocket_pdbstr is None:
            print(f"Unable to generate pocket from receptor: {receptor_path}")
            return None

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pocket_pdbstr)
            pocket_path_temp = f.name

        try:
            pocket = Chem.MolFromPDBFile(pocket_path_temp, removeHs=True)
        finally:
            os.unlink(pocket_path_temp)

        if pocket is None:
            print("Unable to process generated pocket")
            return None
    else:
        raise ValueError("Either pocket_path or receptor_path must be provided")

    return _mols_to_graph(ligand, pocket, dis_threshold)


def pdb_path_list_to_batch(ligand_path_list: List[str],
                            pocket_path_list: List[str] = None,
                            receptor_path_list: List[str] = None,
                            dis_threshold: float = 5.0):
    """
    从PDB文件路径列表构建DGL异构batch图。

    Parameters
    ----------
    ligand_path_list : list of str
        配体PDB文件路径列表
    pocket_path_list : list of str, optional
        口袋PDB文件路径列表（与ligand_path_list一一对应）
    receptor_path_list : list of str, optional
        受体PDB文件路径列表（与ligand_path_list一一对应）
    dis_threshold : float, default 5.0
        距离阈值（埃）

    Returns
    -------
    dgl.DGLHeteroGraph
        DGL异构batch图
    list of str
        成功处理的样本ID列表（使用文件名作为ID）
    """
    graphs = []
    success_ids = []

    if pocket_path_list is not None and len(pocket_path_list) == len(ligand_path_list):
        for ligand_path, pocket_path in zip(ligand_path_list, pocket_path_list):
            try:
                g = pdb_path_to_graph(ligand_path, pocket_path=pocket_path, dis_threshold=dis_threshold)
                if g is not None:
                    graphs.append(g)
                    sample_id = os.path.splitext(os.path.basename(ligand_path))[0]
                    success_ids.append(sample_id)
            except Exception as e:
                print(f"Error processing {ligand_path}: {e}")
                continue
    elif receptor_path_list is not None and len(receptor_path_list) == len(ligand_path_list):
        for ligand_path, receptor_path in zip(ligand_path_list, receptor_path_list):
            try:
                g = pdb_path_to_graph(
                    ligand_path,
                    receptor_path=receptor_path,
                    dis_threshold=dis_threshold
                )
                if g is not None:
                    graphs.append(g)
                    sample_id = os.path.splitext(os.path.basename(ligand_path))[0]
                    success_ids.append(sample_id)
            except Exception as e:
                print(f"Error processing {ligand_path}: {e}")
                continue
    else:
        raise ValueError("pocket_path_list or receptor_path_list must be provided and match ligand_path_list length")

    if len(graphs) == 0:
        return None, []

    return dgl.batch(graphs), success_ids


def pdbstr_to_graph(ligand_pdbstr: str, pocket_pdbstr: str = None,
                     receptor_pdbstr: str = None, dis_threshold: float = 5.0):
    """
    从PDB字符串构建单个DGL异构图。

    Parameters
    ----------
    ligand_pdbstr : str
        配体的PDB格式字符串
    pocket_pdbstr : str, optional
        口袋的PDB格式字符串。如果提供，则忽略receptor参数
    receptor_pdbstr : str, optional
        受体的PDB格式字符串（用于生成口袋）
    dis_threshold : float, default 5.0
        配体-蛋白质相互作用的距离阈值（埃）

    Returns
    -------
    dgl.DGLHeteroGraph or None
        DGL异构图，失败时返回None
    """
    ligand = Chem.MolFromPDBBlock(ligand_pdbstr, removeHs=True)

    if ligand is None:
        print("Unable to process ligand pdbstr")
        return None

    if pocket_pdbstr is not None:
        pocket = Chem.MolFromPDBBlock(pocket_pdbstr, removeHs=True)

        if pocket is None:
            print("Unable to process pocket pdbstr")
            return None
    elif receptor_pdbstr is not None:
        pocket_pdbstr = _generate_pocket_from_receptor(
            receptor_pdbstr, ligand_pdbstr, distance=dis_threshold
        )
        if pocket_pdbstr is None:
            print("Unable to generate pocket from receptor")
            return None

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pocket_pdbstr)
            pocket_path_temp = f.name

        try:
            pocket = Chem.MolFromPDBFile(pocket_path_temp, removeHs=True)
        finally:
            os.unlink(pocket_path_temp)

        if pocket is None:
            print("Unable to process generated pocket")
            return None
    else:
        raise ValueError("Either pocket_pdbstr or receptor_pdbstr must be provided")

    return _mols_to_graph(ligand, pocket, dis_threshold)


def pdbstr_list_to_batch(ligand_pdbstr_list: List[str],
                          pocket_pdbstr_list: List[str] = None,
                          receptor_pdbstr_list: List[str] = None,
                          dis_threshold: float = 5.0):
    """
    从PDB字符串列表构建DGL异构batch图。

    Parameters
    ----------
    ligand_pdbstr_list : list of str
        配体PDB字符串列表
    pocket_pdbstr_list : list of str, optional
        口袋PDB字符串列表（与ligand_pdbstr_list一一对应）
    receptor_pdbstr_list : list of str, optional
        受体PDB字符串列表（与ligand_pdbstr_list一一对应）
    dis_threshold : float, default 5.0
        距离阈值（埃）

    Returns
    -------
    dgl.DGLHeteroGraph
        DGL异构batch图
    list of str
        成功处理的样本ID列表（使用索引作为ID）
    """
    graphs = []
    success_ids = []

    if pocket_pdbstr_list is not None and len(pocket_pdbstr_list) == len(ligand_pdbstr_list):
        for i, (ligand_str, pocket_str) in enumerate(zip(ligand_pdbstr_list, pocket_pdbstr_list)):
            try:
                g = pdbstr_to_graph(ligand_str, pocket_pdbstr=pocket_str, dis_threshold=dis_threshold)
                if g is not None:
                    graphs.append(g)
                    success_ids.append(str(i))
            except Exception as e:
                print(f"Error processing sample {i}: {e}")
                continue
    elif receptor_pdbstr_list is not None and len(receptor_pdbstr_list) == len(ligand_pdbstr_list):
        for i, (ligand_str, receptor_str) in enumerate(zip(ligand_pdbstr_list, receptor_pdbstr_list)):
            try:
                g = pdbstr_to_graph(
                    ligand_str,
                    receptor_pdbstr=receptor_str,
                    dis_threshold=dis_threshold
                )
                if g is not None:
                    graphs.append(g)
                    success_ids.append(str(i))
            except Exception as e:
                print(f"Error processing sample {i}: {e}")
                continue
    else:
        raise ValueError("pocket_pdbstr_list or receptor_pdbstr_list must be provided and match ligand_pdbstr_list length")

    if len(graphs) == 0:
        return None, []

    return dgl.batch(graphs), success_ids


def complex_to_graph(complex: str, ligand_chain: str, receptor_chain: str, dis_threshold: float = 5.0):
    """
    从复杂PDB文件路径构建单个DGL异构图。

    Parameters
    ----------
    complex : str
        PDB格式字符串（包含配体和受体） OR PDB文件路径
    ligand_chain : str
        配体链ID
    receptor_chain : str
        受体链ID
    dis_threshold : float, default 5.0
        配体-蛋白质相互作用的距离阈值（埃）

    Returns
    -------
    dgl.DGLHeteroGraph or None
        DGL异构图，失败时返回None
    """
    pymol.cmd.reinitialize()
    if os.path.isfile(complex):
        pymol.cmd.load(complex, 'complex')
    else:
        pymol.cmd.read_pdbstr(complex, 'complex')
    pymol.cmd.select('Ligand', f'chain {ligand_chain}')
    pymol.cmd.select('Receptor', f'chain {receptor_chain}')
    ligand_str = pymol.cmd.get_pdbstr('Ligand')
    receptor_str = pymol.cmd.get_pdbstr('Receptor')
    
    pocket_pdbstr = _generate_pocket_from_receptor(
        receptor_str,
        ligand_str,
        distance=dis_threshold
    )
    return pdbstr_to_graph(
        ligand_str,
        pocket_pdbstr,
        dis_threshold=dis_threshold
    )


def graphs_to_batch(graphs: List[dgl.DGLHeteroGraph]):
    """
    将图列表打包为batch图。

    Parameters
    ----------
    graphs : list of dgl.DGLHeteroGraph
        DGL异构图列表

    Returns
    -------
    dgl.DGLHeteroGraph
        DGL异构batch图
    """
    if len(graphs) == 0:
        return None
    return dgl.batch(graphs)


def one_of_k_encoding(k, possible_values):
    if k not in possible_values:
        raise ValueError(f"{k} is not a valid value in {possible_values}")
    return [k == e for e in possible_values]


def one_of_k_encoding_unk(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))


def atom_features(mol: Chem.Mol, graph, atom_symbols=['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I'], explicit_H=True):
    for atom in mol.GetAtoms():
        results = one_of_k_encoding_unk(atom.GetSymbol(), atom_symbols + ['Unknown']) + \
                one_of_k_encoding_unk(atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6]) + \
                one_of_k_encoding_unk(atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5, 6]) + \
                one_of_k_encoding_unk(atom.GetHybridization(), [
                    Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
                    Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
                    Chem.rdchem.HybridizationType.SP3D2
                    ]) + [atom.GetIsAromatic()]
        if explicit_H:
            results = results + one_of_k_encoding_unk(atom.GetTotalNumHs(), [0, 1, 2, 3, 4])

        atom_feats = np.array(results).astype(np.float32)
        graph.add_node(atom.GetIdx(), feats=torch.from_numpy(atom_feats))


def edge_features(mol, graph):
    geom = mol.GetConformers()[0].GetPositions()
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        angles_ijk = []
        areas_ijk = []
        dists_ik = []
        for neighbor in mol.GetAtomWithIdx(j).GetNeighbors():
            k = neighbor.GetIdx()
            if mol.GetBondBetweenAtoms(j, k) is not None and i != k:
                vector1 = geom[j] - geom[i]
                vector2 = geom[k] - geom[i]

                angles_ijk.append(angle(vector1, vector2))
                areas_ijk.append(area_triangle(vector1, vector2))
                dists_ik.append(cal_dist(geom[i], geom[k]))

        angles_ijk = np.array(angles_ijk) if angles_ijk != [] else np.array([0.])
        areas_ijk = np.array(areas_ijk) if areas_ijk != [] else np.array([0.])
        dists_ik = np.array(dists_ik) if dists_ik != [] else np.array([0.])
        dist_ij1 = cal_dist(geom[i], geom[j], ord=1)
        dist_ij2 = cal_dist(geom[i], geom[j], ord=2)
        geom_feats = [
            angles_ijk.max() * 0.1,
            angles_ijk.sum() * 0.01,
            angles_ijk.mean() * 0.1,
            areas_ijk.max() * 0.1,
            areas_ijk.sum() * 0.01,
            areas_ijk.mean() * 0.1,
            dists_ik.max() * 0.1,
            dists_ik.sum() * 0.01,
            dists_ik.mean() * 0.1,
            dist_ij1 * 0.1,
            dist_ij2 * 0.1,
        ]

        bond_type = bond.GetBondType()
        basic_feats = [
            bond_type == Chem.rdchem.BondType.SINGLE,
            bond_type == Chem.rdchem.BondType.DOUBLE,
            bond_type == Chem.rdchem.BondType.TRIPLE,
            bond_type == Chem.rdchem.BondType.AROMATIC,
            bond.GetIsConjugated(),
            bond.IsInRing()
        ]

        graph.add_edge(i, j, feats=torch.tensor(basic_feats + geom_feats).float())


def mol2graph(mol):
    graph = nx.Graph()
    atom_features(mol, graph)
    edge_features(mol, graph)

    graph = graph.to_directed()
    x = torch.stack([feats['feats'] for _, feats in graph.nodes(data=True)])
    edge_index = torch.stack([torch.LongTensor((u, v)) for u, v, _ in graph.edges(data=True)]).T
    edge_attr = torch.stack([feats['feats'] for _, _, feats in graph.edges(data=True)])

    return x, edge_index, edge_attr


def geom_feat(pos_i, pos_j, pos_k, angles_ijk, areas_ijk, dists_ik):
    vector1 = pos_j - pos_i
    vector2 = pos_k - pos_i
    angles_ijk.append(angle(vector1, vector2))
    areas_ijk.append(area_triangle(vector1, vector2))
    dists_ik.append(cal_dist(pos_i, pos_k))


def geom_feats(pos_i, pos_j, angles_ijk, areas_ijk, dists_ik):
    angles_ijk = np.array(angles_ijk) if angles_ijk != [] else np.array([0.])
    areas_ijk = np.array(areas_ijk) if areas_ijk != [] else np.array([0.])
    dists_ik = np.array(dists_ik) if dists_ik != [] else np.array([0.])
    dist_ij1 = cal_dist(pos_i, pos_j, ord=1)
    dist_ij2 = cal_dist(pos_i, pos_j, ord=2)
    geom = [
        angles_ijk.max() * 0.1,
        angles_ijk.sum() * 0.01,
        angles_ijk.mean() * 0.1,
        areas_ijk.max() * 0.1,
        areas_ijk.sum() * 0.01,
        areas_ijk.mean() * 0.1,
        dists_ik.max() * 0.1,
        dists_ik.sum() * 0.01,
        dists_ik.mean() * 0.1,
        dist_ij1 * 0.1,
        dist_ij2 * 0.1,
    ]
    return geom


def inter_graph(ligand, pocket, dis_threshold=5.):
    graph_l2p = nx.DiGraph()
    graph_p2l = nx.DiGraph()
    pos_l = ligand.GetConformers()[0].GetPositions()
    pos_p = pocket.GetConformers()[0].GetPositions()
    dis_matrix = distance_matrix(pos_l, pos_p)
    node_idx = np.where(dis_matrix < dis_threshold)

    for i, j in zip(node_idx[0], node_idx[1]):
        ks = node_idx[0][node_idx[1] == j]
        angles_ijk = []
        areas_ijk = []
        dists_ik = []
        for k in ks:
            if k != i:
                geom_feat(pos_l[i], pos_p[j], pos_l[k], angles_ijk, areas_ijk, dists_ik)
        geom = geom_feats(pos_l[i], pos_p[j], angles_ijk, areas_ijk, dists_ik)
        bond_feats = torch.FloatTensor(geom)
        graph_l2p.add_edge(i, j, feats=bond_feats)

        ks = node_idx[1][node_idx[0] == i]
        angles_ijk = []
        areas_ijk = []
        dists_ik = []
        for k in ks:
            if k != j:
                geom_feat(pos_p[j], pos_l[i], pos_p[k], angles_ijk, areas_ijk, dists_ik)
        geom = geom_feats(pos_p[j], pos_l[i], angles_ijk, areas_ijk, dists_ik)
        bond_feats = torch.FloatTensor(geom)
        graph_p2l.add_edge(j, i, feats=bond_feats)

    edge_index_l2p = torch.stack([torch.LongTensor((u, v)) for u, v, _ in graph_l2p.edges(data=True)]).T
    edge_attr_l2p = torch.stack([feats['feats'] for _, _, feats in graph_l2p.edges(data=True)])

    edge_index_p2l = torch.stack([torch.LongTensor((u, v)) for u, v, _ in graph_p2l.edges(data=True)]).T
    edge_attr_p2l = torch.stack([feats['feats'] for _, _, feats in graph_p2l.edges(data=True)])

    return (edge_index_l2p, edge_attr_l2p), (edge_index_p2l, edge_attr_p2l)


if __name__ == '__main__':
    print(complex_to_graph('data_tmp/pdb/pdbstr.pdb', 'B', 'A'))