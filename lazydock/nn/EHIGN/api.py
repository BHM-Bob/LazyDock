from typing import Dict, List, Optional, Tuple, Union

import dgl
import torch
import torch.nn as nn

if __name__ == '__main__':
    from lazydock.nn.EHIGN.data import (complex_to_graph, graphs_to_batch,
                                        pdb_path_list_to_batch,
                                        pdb_path_to_graph,
                                        pdbstr_list_to_batch, pdbstr_to_graph)
    from lazydock.nn.EHIGN.model import DTIPredictor
else:
    from .data import (complex_to_graph, graphs_to_batch,
                       pdb_path_list_to_batch, pdb_path_to_graph,
                       pdbstr_list_to_batch, pdbstr_to_graph)
    from .model import DTIPredictor

NODE_FEAT_SIZE = 35
EDGE_FEAT_SIZE = 17
INTER_EDGE_FEAT_SIZE = 11


def build_model(
    node_feat_size: int = NODE_FEAT_SIZE,
    edge_feat_size: int = EDGE_FEAT_SIZE,
    hidden_feat_size: int = 256,
    layer_num: int = 3
) -> DTIPredictor:
    """
    构建EHIGN模型。

    Parameters
    ----------
    node_feat_size : int, default 35
        节点特征维度（原子特征维度）
    edge_feat_size : int, default 17
        边特征维度（分子内边特征维度）
    hidden_feat_size : int, default 256
        隐藏层维度
    layer_num : int, default 3
        图卷积层数

    Returns
    -------
    DTIPredictor
        EHIGN模型实例
    """
    model = DTIPredictor(
        node_feat_size=node_feat_size,
        edge_feat_size=edge_feat_size,
        hidden_feat_size=hidden_feat_size,
        layer_num=layer_num
    )
    return model


def load_checkpoint(
    model: nn.Module,
    checkpoint_path: str,
    device: Optional[Union[str, torch.device]] = None
) -> None:
    """
    加载模型checkpoint。

    Parameters
    ----------
    model : nn.Module
        模型实例
    checkpoint_path : str
        checkpoint文件路径
    device : str or torch.device, optional
        设备（'cpu', 'cuda', 或具体的torch.device）
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'] if 'model_state_dict' in checkpoint else checkpoint)
    model.to(device)
    model.eval()


def evaluate(
    model: nn.Module,
    batch_graph: dgl.DGLHeteroGraph,
    device: Optional[Union[str, torch.device]] = None
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    使用模型对batch图进行推理。

    Parameters
    ----------
    model : nn.Module
        模型实例
    batch_graph : dgl.DGLHeteroGraph
        DGL异构batch图
    device : str or torch.device, optional
        设备（'cpu', 'cuda', 或具体的torch.device）

    Returns
    -------
    torch.Tensor
        ligand→pocket方向的亲和力分数
    torch.Tensor
        pocket→ligand方向的亲和力分数
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    model.eval()
    with torch.no_grad():
        batch_graph = batch_graph.to(device)
        score_lp, score_pl = model(batch_graph)
        score_lp = score_lp.cpu()
        score_pl = score_pl.cpu()

    return score_lp, score_pl


class EHIGN:
    """
    EHIGN模型推理类，提供统一的API接口。

    Examples
    --------
    >>> from lazydock.nn.EHIGN import EHIGN
    >>>
    >>> # 方式1: 使用文件路径
    >>> ehign = EHIGN(checkpoint_path='path/to/ckp.pt', device='cuda')
    >>> batch, ids = ehign.build_batch_from_paths(
    ...     ligand_path_list=['ligand1.pdb', 'ligand2.pdb'],
    ...     pocket_path_list=['pocket1.pdb', 'pocket2.pdb']
    ... )
    >>> scores = ehign.evaluate(batch)
    >>>
    >>> # 方式2: 使用receptor自动生成pocket
    >>> batch, ids = ehign.build_batch_from_paths(
    ...     ligand_path_list=['ligand1.pdb', 'ligand2.pdb'],
    ...     receptor_path_list=['receptor1.pdb', 'receptor2.pdb']
    ... )
    >>>
    >>> # 方式3: 使用PDB字符串
    >>> batch, ids = ehign.build_batch_from_strs(
    ...     ligand_pdbstr_list=[ligand_str1, ligand_str2],
    ...     pocket_pdbstr_list=[pocket_str1, pocket_str2]
    ... )
    """

    def __init__(
        self,
        model_config: Optional[Dict] = None,
        checkpoint_path: Optional[str] = None,
        device: Optional[Union[str, torch.device]] = None
    ):
        """
        初始化EHIGN模型。

        Parameters
        ----------
        model_config : dict, optional
            模型配置字典，包含：
            - node_feat_size: int, 节点特征维度（默认35）
            - edge_feat_size: int, 边特征维度（默认17）
            - hidden_feat_size: int, 隐藏层维度（默认128）
            - layer_num: int, 图卷积层数（默认3）
        checkpoint_path : str, optional
            checkpoint文件路径，如果提供则自动加载模型权重
        device : str or torch.device, optional
            设备（'cpu', 'cuda', 或具体的torch.device）
        """
        if device is None:
            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.device = device

        if model_config is None:
            model_config = {}
        self.model_config = {
            'node_feat_size': model_config.get('node_feat_size', NODE_FEAT_SIZE),
            'edge_feat_size': model_config.get('edge_feat_size', EDGE_FEAT_SIZE),
            'hidden_feat_size': model_config.get('hidden_feat_size', 256),
            'layer_num': model_config.get('layer_num', 3)
        }

        self.model = build_model(**self.model_config)
        self.model.to(self.device)

        self.checkpoint_path = checkpoint_path
        if checkpoint_path is not None:
            self.load_model(checkpoint_path)

    def load_model(self, checkpoint_path: Optional[str] = None) -> None:
        """
        加载模型checkpoint。

        Parameters
        ----------
        checkpoint_path : str, optional
            checkpoint文件路径，如果为None则使用初始化时提供的路径
        """
        if checkpoint_path is None:
            checkpoint_path = self.checkpoint_path

        if checkpoint_path is None:
            raise ValueError("checkpoint_path must be provided")

        load_checkpoint(self.model, checkpoint_path, self.device)

    def build_graph_from_path(
        self,
        ligand_path: str,
        pocket_path: str = None,
        receptor_path: str = None,
        dis_threshold: float = 5.0
    ) -> Optional[dgl.DGLHeteroGraph]:
        """
        从PDB文件路径构建单个DGL异构图。

        Parameters
        ----------
        ligand_path : str
            配体PDB文件路径
        pocket_path : str, optional
            口袋PDB文件路径
        receptor_path : str, optional
            受体PDB文件路径（用于生成口袋）
        dis_threshold : float, default 5.0
            距离阈值（埃）

        Returns
        -------
        dgl.DGLHeteroGraph or None
            DGL异构图
        """
        return pdb_path_to_graph(
            ligand_path=ligand_path,
            pocket_path=pocket_path,
            receptor_path=receptor_path,
            dis_threshold=dis_threshold
        )

    def build_batch_from_paths(
        self,
        ligand_path_list: List[str],
        pocket_path_list: List[str] = None,
        receptor_path_list: List[str] = None,
        dis_threshold: float = 5.0
    ):
        """
        从PDB文件路径列表构建DGL异构batch图。

        Parameters
        ----------
        ligand_path_list : list of str
            配体PDB文件路径列表
        pocket_path_list : list of str, optional
            口袋PDB文件路径列表
        receptor_path_list : list of str, optional
            受体PDB文件路径列表
        dis_threshold : float, default 5.0
            距离阈值（埃）

        Returns
        -------
        dgl.DGLHeteroGraph
            DGL异构batch图
        list of str
            成功处理的样本ID列表
        """
        return pdb_path_list_to_batch(
            ligand_path_list=ligand_path_list,
            pocket_path_list=pocket_path_list,
            receptor_path_list=receptor_path_list,
            dis_threshold=dis_threshold
        )

    def build_graph_from_str(
        self,
        ligand_pdbstr: str,
        pocket_pdbstr: str = None,
        receptor_pdbstr: str = None,
        dis_threshold: float = 5.0
    ) -> Union[dgl.DGLHeteroGraph, None]:
        """
        从PDB字符串构建单个DGL异构图。

        Parameters
        ----------
        ligand_pdbstr : str
            配体PDB字符串
        pocket_pdbstr : str, optional
            口袋PDB字符串
        receptor_pdbstr : str, optional
            受体PDB字符串（用于生成口袋）
        dis_threshold : float, default 5.0
            距离阈值（埃）

        Returns
        -------
        dgl.DGLHeteroGraph or None
            DGL异构图
        """
        return pdbstr_to_graph(
            ligand_pdbstr=ligand_pdbstr,
            pocket_pdbstr=pocket_pdbstr,
            receptor_pdbstr=receptor_pdbstr,
            dis_threshold=dis_threshold
        )

    def build_batch_from_strs(
        self,
        ligand_pdbstr_list: List[str],
        pocket_pdbstr_list: List[str] = None,
        receptor_pdbstr_list: List[str] = None,
        dis_threshold: float = 5.0
    ):
        """
        从PDB字符串列表构建DGL异构batch图。

        Parameters
        ----------
        ligand_pdbstr_list : list of str
            配体PDB字符串列表
        pocket_pdbstr_list : list of str, optional
            口袋PDB字符串列表
        receptor_pdbstr_list : list of str, optional
            受体PDB字符串列表
        dis_threshold : float, default 5.0
            距离阈值（埃）

        Returns
        -------
        dgl.DGLHeteroGraph
            DGL异构batch图
        list of str
            成功处理的样本ID列表
        """
        return pdbstr_list_to_batch(
            ligand_pdbstr_list=ligand_pdbstr_list,
            pocket_pdbstr_list=pocket_pdbstr_list,
            receptor_pdbstr_list=receptor_pdbstr_list,
            dis_threshold=dis_threshold
        )

    def batch_graphs(self, graphs: List[dgl.DGLHeteroGraph]) -> Union[dgl.DGLHeteroGraph, None]:
        """
        将图列表打包为batch图。

        Parameters
        ----------
        graphs : list of dgl.DGLHeteroGraph
            DGL异构图列表

        Returns
        -------
        dgl.DGLHeteroGraph or None
            DGL异构batch图
        """
        return graphs_to_batch(graphs)

    def evaluate(self, batch_graph: dgl.DGLHeteroGraph) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        对batch图进行推理。

        Parameters
        ----------
        batch_graph : dgl.DGLHeteroGraph
            DGL异构batch图

        Returns
        -------
        torch.Tensor
            ligand→pocket方向的亲和力分数
        torch.Tensor
            pocket→ligand方向的亲和力分数
        """
        return evaluate(self.model, batch_graph, self.device)

    def predict_from_paths(
        self,
        ligand_path_list: List[str],
        pocket_path_list: List[str] = None,
        receptor_path_list: List[str] = None,
        dis_threshold: float = 5.0
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        对PDB文件路径列表进行预测。

        Parameters
        ----------
        ligand_path_list : list of str
            配体PDB文件路径列表
        pocket_path_list : list of str, optional
            口袋PDB文件路径列表
        receptor_path_list : list of str, optional
            受体PDB文件路径列表
        dis_threshold : float, default 5.0
            距离阈值（埃）

        Returns
        -------
        torch.Tensor
            平均亲和力分数_lp
        torch.Tensor
            平均亲和力分数_pl
        """
        batch, _ = self.build_batch_from_paths(
            ligand_path_list=ligand_path_list,
            pocket_path_list=pocket_path_list,
            receptor_path_list=receptor_path_list,
            dis_threshold=dis_threshold
        )

        if batch is None:
            return torch.tensor([]), torch.tensor([])

        score_lp, score_pl = self.evaluate(batch)
        return score_lp.mean(), score_pl.mean()

    def predict_from_strs(
        self,
        ligand_pdbstr_list: List[str],
        pocket_pdbstr_list: List[str] = None,
        receptor_pdbstr_list: List[str] = None,
        dis_threshold: float = 5.0
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        对PDB字符串列表进行预测。

        Parameters
        ----------
        ligand_pdbstr_list : list of str
            配体PDB字符串列表
        pocket_pdbstr_list : list of str, optional
            口袋PDB字符串列表
        receptor_pdbstr_list : list of str, optional
            受体PDB字符串列表
        dis_threshold : float, default 5.0
            距离阈值（埃）

        Returns
        -------
        torch.Tensor
            平均亲和力分数_lp
        torch.Tensor
            平均亲和力分数_pl
        """
        batch, _ = self.build_batch_from_strs(
            ligand_pdbstr_list=ligand_pdbstr_list,
            pocket_pdbstr_list=pocket_pdbstr_list,
            receptor_pdbstr_list=receptor_pdbstr_list,
            dis_threshold=dis_threshold
        )

        if batch is None:
            return torch.tensor([]), torch.tensor([])

        score_lp, score_pl = self.evaluate(batch)
        return score_lp.mean(), score_pl.mean()
    
    def predict_from_complexes(
        self,
        complex_list: List[str],
        ligand_chain: str,
        receptor_chain: str,
        dis_threshold: float = 5.0
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        对complex列表进行预测。

        Parameters
        ----------
        complex_list : list of str
            complex字符串列表（包含配体和受体）OR complex文件路径列表
        ligand_chain : str
            配体链ID
        receptor_chain : str
            受体链ID
        dis_threshold : float, default 5.0
            距离阈值（埃）
        
        Returns
        -------
        torch.Tensor
            平均亲和力分数_lp
        torch.Tensor
            平均亲和力分数_pl
        """
        graphs = []
        for complex in complex_list:
            graph = complex_to_graph(complex, ligand_chain, receptor_chain, dis_threshold)
            if graph is not None:
                graphs.append(graph)
        
        if len(graphs) == 0:
            return torch.tensor([]), torch.tensor([])

        batch = graphs_to_batch(graphs)
        score_lp, score_pl = self.evaluate(batch)
        return score_lp.mean(), score_pl.mean()

if __name__ == '__main__':
    predictor = EHIGN(checkpoint_path='data_tmp/nn/EHIGN_ckp.pt', device='cuda')
    print(predictor.predict_from_complexes(
        complex_list=['data_tmp/pdb/pdbstr.pdb'],
        ligand_chain='B',
        receptor_chain='A'
    ))
