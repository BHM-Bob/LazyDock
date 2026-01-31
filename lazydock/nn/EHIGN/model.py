from functools import partial

import dgl
import dgl.function as fn
import torch
import torch as th
import torch.nn as nn
import torch.nn.functional as F
from dgl.nn.pytorch import edge_softmax
from dgl.utils import expand_as_pair
from torch import nn
from torch.nn import functional as F


class CIGConv(nn.Module):

    def __init__(self, input_dim, output_dim, drop=0.1):
        super(CIGConv, self).__init__()

        self.mlp = nn.Sequential(
            nn.Linear(input_dim, output_dim),
            nn.Dropout(drop),
            nn.LeakyReLU(),
            nn.BatchNorm1d(output_dim)
        )

    def message(self, edges):
        return {'m': F.relu(edges.src['hn'] + edges.data['he'])}

    def forward(self, graph, node_feat, edge_feat):
        with graph.local_scope():
            feat_src, feat_dst = expand_as_pair(node_feat, graph)
            graph.srcdata['hn'] = feat_src
            graph.edata['he'] = edge_feat
            graph.update_all(self.message, fn.sum('m', 'neigh'))
            rst = feat_dst + graph.dstdata['neigh']
            rst = self.mlp(rst)

            return rst



class NIGConv(nn.Module):
    def __init__(self,
                 in_feats,
                 out_feats,
                 feat_drop=0.,
                 bias=True):
        super(NIGConv, self).__init__()

        self._in_src_feats, self._in_dst_feats = expand_as_pair(in_feats)
        self._out_feats = out_feats
        self.feat_drop = nn.Dropout(feat_drop)
        self.fc_neigh = nn.Linear(self._in_src_feats, out_feats, bias=False)
        self.fc_self = nn.Linear(self._in_dst_feats, out_feats, bias=False)
        if bias:
            self.bias = nn.parameter.Parameter(torch.zeros(self._out_feats))
        else:
            self.register_buffer('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        gain = nn.init.calculate_gain('relu')
        nn.init.xavier_uniform_(self.fc_self.weight, gain=gain)
        nn.init.xavier_uniform_(self.fc_neigh.weight, gain=gain)
        
    def message_func(edges):
        return {'m': edges.src['h']}

    def forward(self, graph, feat, edge_weight=None):
        with graph.local_scope():
            if isinstance(feat, tuple):
                feat_src = self.feat_drop(feat[0])
                feat_dst = self.feat_drop(feat[1])
            else:
                feat_src = feat_dst = self.feat_drop(feat)
                if graph.is_block:
                    feat_dst = feat_src[:graph.number_of_dst_nodes()]
            msg_fn = self.message_func
            if edge_weight is not None:
                assert edge_weight.shape[0] == graph.number_of_edges()
                graph.edata['_edge_weight'] = edge_weight
                msg_fn = fn.u_mul_e('h', '_edge_weight', 'm')

            h_self = feat_dst

            # Handle the case of graphs without edges
            if graph.number_of_edges() == 0:
                graph.dstdata['neigh'] = torch.zeros(
                    feat_dst.shape[0], self._in_src_feats).to(feat_dst)

            # Determine whether to apply linear transformation before message passing A(XW)
            lin_before_mp = self._in_src_feats > self._out_feats

            # Message Passing
            graph.srcdata['h'] = self.fc_neigh(feat_src) if lin_before_mp else feat_src
            graph.update_all(msg_fn, fn.mean('m', 'neigh'))
            h_neigh = graph.dstdata['neigh']
            if not lin_before_mp:
                h_neigh = self.fc_neigh(h_neigh)
  
            rst = self.fc_self(h_self) + h_neigh

            # bias term
            if self.bias is not None:
                rst = rst + self.bias

            return rst

class HeteroGraphConv(nn.Module):
    r"""A generic module for computing convolution on heterogeneous graphs.

    The heterograph convolution applies sub-modules on their associating
    relation graphs, which reads the features from source nodes and writes the
    updated ones to destination nodes. If multiple relations have the same
    destination node types, their results are aggregated by the specified method.
    If the relation graph has no edge, the corresponding module will not be called.

    Pseudo-code:

    .. code::

        outputs = {nty : [] for nty in g.dsttypes}
        # Apply sub-modules on their associating relation graphs in parallel
        for relation in g.canonical_etypes:
            stype, etype, dtype = relation
            dstdata = relation_submodule(g[relation], ...)
            outputs[dtype].append(dstdata)

        # Aggregate the results for each destination node type
        rsts = {}
        for ntype, ntype_outputs in outputs.items():
            if len(ntype_outputs) != 0:
                rsts[ntype] = aggregate(ntype_outputs)
        return rsts

    Examples
    --------

    Create a heterograph with three types of relations and nodes.

    >>> import dgl
    >>> g = dgl.heterograph({
    ...     ('user', 'follows', 'user') : edges1,
    ...     ('user', 'plays', 'game') : edges2,
    ...     ('store', 'sells', 'game')  : edges3})

    Create a ``HeteroGraphConv`` that applies different convolution modules to
    different relations. Note that the modules for ``'follows'`` and ``'plays'``
    do not share weights.

    >>> import dgl.nn.pytorch as dglnn
    >>> conv = dglnn.HeteroGraphConv({
    ...     'follows' : dglnn.GraphConv(...),
    ...     'plays' : dglnn.GraphConv(...),
    ...     'sells' : dglnn.SAGEConv(...)},
    ...     aggregate='sum')

    Call forward with some ``'user'`` features. This computes new features for both
    ``'user'`` and ``'game'`` nodes.

    >>> import torch as th
    >>> h1 = {'user' : th.randn((g.number_of_nodes('user'), 5))}
    >>> h2 = conv(g, h1)
    >>> print(h2.keys())
    dict_keys(['user', 'game'])

    Call forward with both ``'user'`` and ``'store'`` features. Because both the
    ``'plays'`` and ``'sells'`` relations will update the ``'game'`` features,
    their results are aggregated by the specified method (i.e., summation here).

    >>> f1 = {'user' : ..., 'store' : ...}
    >>> f2 = conv(g, f1)
    >>> print(f2.keys())
    dict_keys(['user', 'game'])

    Call forward with some ``'store'`` features. This only computes new features
    for ``'game'`` nodes.

    >>> g1 = {'store' : ...}
    >>> g2 = conv(g, g1)
    >>> print(g2.keys())
    dict_keys(['game'])

    Call forward with a pair of inputs is allowed and each submodule will also
    be invoked with a pair of inputs.

    >>> x_src = {'user' : ..., 'store' : ...}
    >>> x_dst = {'user' : ..., 'game' : ...}
    >>> y_dst = conv(g, (x_src, x_dst))
    >>> print(y_dst.keys())
    dict_keys(['user', 'game'])

    Parameters
    ----------
    mods : dict[str, nn.Module]
        Modules associated with every edge types. The forward function of each
        module must have a `DGLHeteroGraph` object as the first argument, and
        its second argument is either a tensor object representing the node
        features or a pair of tensor object representing the source and destination
        node features.
    aggregate : str, callable, optional
        Method for aggregating node features generated by different relations.
        Allowed string values are 'sum', 'max', 'min', 'mean', 'stack'.
        The 'stack' aggregation is performed along the second dimension, whose order
        is deterministic.
        User can also customize the aggregator by providing a callable instance.
        For example, aggregation by summation is equivalent to the follows:

        .. code::

            def my_agg_func(tensors, dsttype):
                # tensors: is a list of tensors to aggregate
                # dsttype: string name of the destination node type for which the
                #          aggregation is performed
                stacked = torch.stack(tensors, dim=0)
                return torch.sum(stacked, dim=0)

    Attributes
    ----------
    mods : dict[str, nn.Module]
        Modules associated with every edge types.
    """
    def __init__(self, mods, aggregate='sum'):
        super(HeteroGraphConv, self).__init__()
        self.mods = nn.ModuleDict(mods)
        # Do not break if graph has 0-in-degree nodes.
        # Because there is no general rule to add self-loop for heterograph.
        for _, v in self.mods.items():
            set_allow_zero_in_degree_fn = getattr(v, 'set_allow_zero_in_degree', None)
            if callable(set_allow_zero_in_degree_fn):
                set_allow_zero_in_degree_fn(True)
        if isinstance(aggregate, str):
            self.agg_fn = get_aggregate_fn(aggregate)
        else:
            self.agg_fn = aggregate

    def forward(self, g, inputs, mod_args=None, mod_kwargs=None):
        """Forward computation

        Invoke the forward function with each module and aggregate their results.

        Parameters
        ----------
        g : DGLHeteroGraph
            Graph data.
        inputs : dict[str, Tensor] or pair of dict[str, Tensor]
            Input node features.
        mod_args : dict[str, tuple[any]], optional
            Extra positional arguments for the sub-modules.
        mod_kwargs : dict[str, dict[str, any]], optional
            Extra key-word arguments for the sub-modules.

        Returns
        -------
        dict[str, Tensor]
            Output representations for every types of nodes.
        """
        if mod_args is None:
            mod_args = {}
        if mod_kwargs is None:
            mod_kwargs = {}
        outputs = {nty : [] for nty in g.dsttypes}
        if isinstance(inputs, tuple) or g.is_block:
            if isinstance(inputs, tuple):
                src_inputs, dst_inputs = inputs
            else:
                src_inputs = inputs
                dst_inputs = {k: v[:g.number_of_dst_nodes(k)] for k, v in inputs.items()}

            for stype, etype, dtype in g.canonical_etypes:
                rel_graph = g[stype, etype, dtype]
                if stype not in src_inputs or dtype not in dst_inputs:
                    continue
                dstdata = self.mods[etype](
                    rel_graph,
                    (src_inputs[stype], dst_inputs[dtype]),
                    *mod_args.get(etype, ()),
                    **mod_kwargs.get(etype, {}))
                outputs[dtype].append(dstdata)
        else:
            for stype, etype, dtype in g.canonical_etypes:
                rel_graph = g[stype, etype, dtype]
                
                edge_feat = rel_graph.edata['e']
                if stype not in inputs:
                    continue
                # modify the function so that it can process edge features
                dstdata = self.mods[etype](
                    rel_graph,
                    (inputs[stype], inputs[dtype]),
                    edge_feat,
                    *mod_args.get(etype, ()),
                    **mod_kwargs.get(etype, {}))
                outputs[dtype].append(dstdata)
        rsts = {}
        for nty, alist in outputs.items():
            if len(alist) != 0:
                rsts[nty] = self.agg_fn(alist, nty)
        return rsts


def _max_reduce_func(inputs, dim):
    return th.max(inputs, dim=dim)[0]

def _min_reduce_func(inputs, dim):
    return th.min(inputs, dim=dim)[0]

def _sum_reduce_func(inputs, dim):
    return th.sum(inputs, dim=dim)

def _mean_reduce_func(inputs, dim):
    return th.mean(inputs, dim=dim)

def _stack_agg_func(inputs, dsttype): # pylint: disable=unused-argument
    if len(inputs) == 0:
        return None
    return th.stack(inputs, dim=1)

def _agg_func(inputs, dsttype, fn): # pylint: disable=unused-argument
    if len(inputs) == 0:
        return None
    stacked = th.stack(inputs, dim=0)
    return fn(stacked, dim=0)

def get_aggregate_fn(agg):
    """Internal function to get the aggregation function for node data
    generated from different relations.

    Parameters
    ----------
    agg : str
        Method for aggregating node features generated by different relations.
        Allowed values are 'sum', 'max', 'min', 'mean', 'stack'.

    Returns
    -------
    callable
        Aggregator function that takes a list of tensors to aggregate
        and returns one aggregated tensor.
    """
    if agg == 'sum':
        fn = _sum_reduce_func
    elif agg == 'max':
        fn = _max_reduce_func
    elif agg == 'min':
        fn = _min_reduce_func
    elif agg == 'mean':
        fn = _mean_reduce_func
    elif agg == 'stack':
        fn = None  # will not be called
    else:
        print('Invalid cross type aggregator. Must be one of '
                       '"sum", "max", "min", "mean" or "stack". But got "%s"' % agg)
    if agg == 'stack':
        return _stack_agg_func
    else:
        return partial(_agg_func, fn=fn)

class HeteroLinear(nn.Module):
    """Apply linear transformations on heterogeneous inputs.

    Parameters
    ----------
    in_size : dict[key, int]
        Input feature size for heterogeneous inputs. A key can be a string or a tuple of strings.
    out_size : int
        Output feature size.
    bias : bool, optional
        If True, learns a bias term. Defaults: ``True``.

    Examples
    --------

    >>> import dgl
    >>> import torch
    >>> from dgl.nn import HeteroLinear

    >>> layer = HeteroLinear({'user': 1, ('user', 'follows', 'user'): 2}, 3)
    >>> in_feats = {'user': torch.randn(2, 1), ('user', 'follows', 'user'): torch.randn(3, 2)}
    >>> out_feats = layer(in_feats)
    >>> print(out_feats['user'].shape)
    torch.Size([2, 3])
    >>> print(out_feats[('user', 'follows', 'user')].shape)
    torch.Size([3, 3])
    """
    def __init__(self, in_size, out_size, bias=True):
        super(HeteroLinear, self).__init__()

        self.linears = nn.ModuleDict()
        for typ, typ_in_size in in_size.items():
            self.linears[str(typ)] = nn.Linear(typ_in_size, out_size, bias=bias)

    def forward(self, feat):
        """Forward function

        Parameters
        ----------
        feat : dict[key, Tensor]
            Heterogeneous input features. It maps keys to features.

        Returns
        -------
        dict[key, Tensor]
            Transformed features.
        """
        out_feat = dict()
        for typ, typ_feat in feat.items():
            out_feat[typ] = self.linears[str(typ)](typ_feat)

        return out_feat


class HeteroEmbedding(nn.Module):
    """Create a heterogeneous embedding table.

    It internally contains multiple ``torch.nn.Embedding`` with different dictionary sizes.

    Parameters
    ----------
    num_embeddings : dict[key, int]
        Size of the dictionaries. A key can be a string or a tuple of strings.
    embedding_dim : int
        Size of each embedding vector.

    Examples
    --------

    >>> import dgl
    >>> import torch
    >>> from dgl.nn import HeteroEmbedding

    >>> layer = HeteroEmbedding({'user': 2, ('user', 'follows', 'user'): 3}, 4)
    >>> # Get the heterogeneous embedding table
    >>> embeds = layer.weight
    >>> print(embeds['user'].shape)
    torch.Size([2, 4])
    >>> print(embeds[('user', 'follows', 'user')].shape)
    torch.Size([3, 4])

    >>> # Get the embeddings for a subset
    >>> input_ids = {'user': torch.LongTensor([0]),
    ...              ('user', 'follows', 'user'): torch.LongTensor([0, 2])}
    >>> embeds = layer(input_ids)
    >>> print(embeds['user'].shape)
    torch.Size([1, 4])
    >>> print(embeds[('user', 'follows', 'user')].shape)
    torch.Size([2, 4])
    """
    def __init__(self, num_embeddings, embedding_dim):
        super(HeteroEmbedding, self).__init__()

        self.embeds = nn.ModuleDict()
        self.raw_keys = dict()
        for typ, typ_num_rows in num_embeddings.items():
            self.embeds[str(typ)] = nn.Embedding(typ_num_rows, embedding_dim)
            self.raw_keys[str(typ)] = typ

    @property
    def weight(self):
        """Get the heterogeneous embedding table

        Returns
        -------
        dict[key, Tensor]
            Heterogeneous embedding table
        """
        return {self.raw_keys[typ]: emb.weight for typ, emb in self.embeds.items()}

    def forward(self, input_ids):
        """Forward function

        Parameters
        ----------
        input_ids : dict[key, Tensor]
            The row IDs to retrieve embeddings. It maps a key to key-specific IDs.

        Returns
        -------
        dict[key, Tensor]
            The retrieved embeddings.
        """
        embeds = dict()
        for typ, typ_ids in input_ids.items():
            embeds[typ] = self.embeds[str(typ)](typ_ids)

        return embeds


class DTIPredictor(nn.Module):
    def __init__(self, node_feat_size, edge_feat_size, hidden_feat_size, layer_num=3):
        super(DTIPredictor, self).__init__()

        self.convs = nn.ModuleList()

        for _ in range(layer_num):
            convl = CIGConv(hidden_feat_size, hidden_feat_size)
            convp = CIGConv(hidden_feat_size, hidden_feat_size)
            convlp = NIGConv(hidden_feat_size, hidden_feat_size, feat_drop=0.1)
            convpl = NIGConv(hidden_feat_size, hidden_feat_size, feat_drop=0.1)
            conv = HeteroGraphConv(
                    {
                        'intra_l' : convl,
                        'intra_p' : convp,
                        'inter_l2p' : convlp,
                        'inter_p2l': convpl
                    }
                )
            self.convs.append(conv)

        self.lin_node_l = nn.Linear(node_feat_size, hidden_feat_size)
        self.lin_node_p = nn.Linear(node_feat_size, hidden_feat_size)
        self.lin_edge_ll = nn.Linear(edge_feat_size, hidden_feat_size)
        self.lin_edge_pp = nn.Linear(edge_feat_size, hidden_feat_size)

        self.lin_edge_lp = nn.Linear(11, hidden_feat_size)
        self.lin_edge_pl = nn.Linear(11, hidden_feat_size)

        # atom-atom affinities
        self.inter_atompairs = AtomAtomAffinities(hidden_feat_size, hidden_feat_size, hidden_feat_size)

        # bias correction
        self.bias_ligandpocket = BiasCorrectionLigandPocket(hidden_feat_size, hidden_feat_size, hidden_feat_size)
        self.bias_pocketligand = BiasCorrectionPocketLigand(hidden_feat_size, hidden_feat_size, hidden_feat_size)

    def forward(self, bg):
        atom_feats = bg.ndata['h']
        bond_feats = bg.edata['e']

        atom_feats = {
            'ligand':self.lin_node_l(atom_feats['ligand']),
            'pocket':self.lin_node_p(atom_feats['pocket'])
        }
        bond_feats = {
            ('ligand', 'intra_l', 'ligand'):self.lin_edge_ll(bond_feats[('ligand', 'intra_l', 'ligand')]),
            ('pocket', 'intra_p', 'pocket'):self.lin_edge_pp(bond_feats[('pocket', 'intra_p', 'pocket')]),       
            ('ligand', 'inter_l2p', 'pocket'):self.lin_edge_lp(bond_feats[('ligand', 'inter_l2p', 'pocket')]),    
            ('pocket', 'inter_p2l', 'ligand'):self.lin_edge_pl(bond_feats[('pocket', 'inter_p2l', 'ligand')]),        
        }

        bg.edata['e'] = bond_feats

        rsts = atom_feats
        for conv in self.convs:
            rsts = conv(bg, rsts)

        bg.nodes['ligand'].data['h'] = rsts['ligand']
        bg.nodes['pocket'].data['h'] = rsts['pocket']

        # atom-atom affinities
        atompairs_lp, atompairs_pl = self.inter_atompairs(bg)

        # bias correction
        bias_lp = self.bias_ligandpocket(bg)
        bias_pl = self.bias_pocketligand(bg)

        return (atompairs_lp - bias_lp).view(-1), (atompairs_pl - bias_pl).view(-1)

class AtomAtomAffinities(nn.Module):
    def __init__(self, node_feat_size, edge_feat_size, hidden_feat_size):
        super(AtomAtomAffinities, self).__init__()
        self.prj_lp_src = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_lp_dst = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_lp_edge = nn.Linear(edge_feat_size, hidden_feat_size)

        self.prj_pl_src = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_pl_dst = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_pl_edge = nn.Linear(edge_feat_size, hidden_feat_size)

        self.fc_lp = nn.Linear(hidden_feat_size, 1)
        self.fc_pl = nn.Linear(hidden_feat_size, 1)

    def apply_interactions(self, edges):
        return {'i' : edges.data['e'] * edges.src['h'] * edges.dst['h']}

    def forward(self, g):
        with g.local_scope():
            node_ligand_feats = g.nodes['ligand'].data['h']
            node_pocket_feats = g.nodes['pocket'].data['h']
            edge_lp_feat = g.edges['inter_l2p'].data['e']
            edge_pl_feat = g.edges['inter_p2l'].data['e']

            g.nodes['ligand'].data['h'] = self.prj_lp_src(node_ligand_feats)
            g.nodes['pocket'].data['h'] = self.prj_lp_dst(node_pocket_feats)
            g.edges['inter_l2p'].data['e'] = self.prj_lp_edge(edge_lp_feat)
            g.apply_edges(self.apply_interactions, etype='inter_l2p')
            logit_lp = self.fc_lp(g.edges['inter_l2p'].data['i'])
            g.edges['inter_l2p'].data['logit_lp'] = logit_lp
            logit_lp = dgl.sum_edges(g, 'logit_lp', etype='inter_l2p')

            g.nodes['ligand'].data['h'] = self.prj_pl_src(node_ligand_feats)
            g.nodes['pocket'].data['h'] = self.prj_pl_dst(node_pocket_feats)
            g.edges['inter_p2l'].data['e'] = self.prj_pl_edge(edge_pl_feat)
            g.apply_edges(self.apply_interactions, etype='inter_p2l')
            logit_pl = self.fc_pl(g.edges['inter_p2l'].data['i'])
            g.edges['inter_p2l'].data['logit_pl'] = logit_pl
            logit_pl = dgl.sum_edges(g, 'logit_pl', etype='inter_p2l')

            return logit_lp, logit_pl

class BiasCorrectionLigandPocket(nn.Module):
    def __init__(self, node_feat_size, edge_feat_size, hidden_feat_size):
        super(BiasCorrectionLigandPocket, self).__init__()
        self.prj_src = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_dst = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_edge = nn.Linear(edge_feat_size, hidden_feat_size)

        self.w_src = nn.Linear(node_feat_size, hidden_feat_size)
        self.w_dst = nn.Linear(node_feat_size, hidden_feat_size)
        self.w_edge = nn.Linear(edge_feat_size, hidden_feat_size)

        self.lin_att = nn.Sequential(
            nn.PReLU(),
            nn.Linear(hidden_feat_size, 1)
        )

        self.fc = FC(hidden_feat_size, 200, 2, 0.1, 1)

    def get_weight(self, edges):
        w = edges.src['h'] + edges.dst['h'] + edges.data['e']
        w = self.lin_att(w)

        return {'w' : w}
    
    def apply_scores(self, edges):
        return {'l' : edges.data['a'] * edges.data['e'] * edges.src['h'] * edges.dst['h']}

    def forward(self, g):
        with g.local_scope():
            node_ligand_feats = g.nodes['ligand'].data['h']
            node_pocket_feats = g.nodes['pocket'].data['h']
            edge_feat = g.edges['inter_l2p'].data['e']

            g.nodes['ligand'].data['h'] = self.prj_src(node_ligand_feats)
            g.nodes['pocket'].data['h'] = self.prj_dst(node_pocket_feats)
            g.edges['inter_l2p'].data['e'] = self.prj_edge(edge_feat)
            g.apply_edges(self.get_weight, etype='inter_l2p')
            scores = edge_softmax(g['inter_l2p'], g.edges['inter_l2p'].data['w'])

            g.edges['inter_l2p'].data['a'] = scores
            g.nodes['ligand'].data['h'] = self.w_src(node_ligand_feats)
            g.nodes['pocket'].data['h'] = self.w_dst(node_pocket_feats)
            g.edges['inter_l2p'].data['e'] = self.w_edge(edge_feat)
            g.apply_edges(self.apply_scores, etype='inter_l2p')

            bias = self.fc(dgl.sum_edges(g, 'l', etype='inter_l2p')) 
            
            return bias

class BiasCorrectionPocketLigand(nn.Module):
    def __init__(self, node_feat_size, edge_feat_size, hidden_feat_size):
        super(BiasCorrectionPocketLigand, self).__init__()
        self.prj_src = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_dst = nn.Linear(node_feat_size, hidden_feat_size)
        self.prj_edge = nn.Linear(edge_feat_size, hidden_feat_size)

        self.w_src = nn.Linear(node_feat_size, hidden_feat_size)
        self.w_dst = nn.Linear(node_feat_size, hidden_feat_size)
        self.w_edge = nn.Linear(edge_feat_size, hidden_feat_size)

        self.lin_att = nn.Sequential(
            nn.PReLU(),
            nn.Linear(hidden_feat_size, 1)
        )

        self.fc = FC(hidden_feat_size, 200, 2, 0.1, 1)

    def get_weight(self, edges):
        w = edges.src['h'] + edges.dst['h'] + edges.data['e']
        w = self.lin_att(w)

        return {'w' : w}
    
    def apply_scores(self, edges):
        return {'l' : edges.data['a'] * edges.data['e'] * edges.src['h'] * edges.dst['h']}

    def forward(self, g):
        with g.local_scope():
            node_ligand_feats = g.nodes['ligand'].data['h']
            node_pocket_feats = g.nodes['pocket'].data['h']
            edge_feat = g.edges['inter_p2l'].data['e']

            g.nodes['ligand'].data['h'] = self.prj_src(node_ligand_feats)
            g.nodes['pocket'].data['h'] = self.prj_dst(node_pocket_feats)
            g.edges['inter_p2l'].data['e'] = self.prj_edge(edge_feat)
            g.apply_edges(self.get_weight, etype='inter_p2l')
            scores = edge_softmax(g['inter_p2l'], g.edges['inter_p2l'].data['w'])

            g.edges['inter_p2l'].data['a'] = scores
            g.nodes['ligand'].data['h'] = self.w_src(node_ligand_feats)
            g.nodes['pocket'].data['h'] = self.w_dst(node_pocket_feats)
            g.edges['inter_p2l'].data['e'] = self.w_edge(edge_feat)
            g.apply_edges(self.apply_scores, etype='inter_p2l')

            bias = self.fc(dgl.sum_edges(g, 'l', etype='inter_p2l')) 
            
            return bias

class FC(nn.Module):
    def __init__(self, d_graph_layer, d_FC_layer, n_FC_layer, dropout, n_tasks):
        super(FC, self).__init__()
        self.d_graph_layer = d_graph_layer
        self.d_FC_layer = d_FC_layer
        self.n_FC_layer = n_FC_layer
        self.dropout = dropout
        self.predict = nn.ModuleList()
        for j in range(self.n_FC_layer):
            if j == 0:
                self.predict.append(nn.Linear(self.d_graph_layer, self.d_FC_layer))
                self.predict.append(nn.Dropout(self.dropout))
                self.predict.append(nn.LeakyReLU())
                self.predict.append(nn.BatchNorm1d(d_FC_layer))
            if j == self.n_FC_layer - 1:
                self.predict.append(nn.Linear(self.d_FC_layer, n_tasks))
            else:
                self.predict.append(nn.Linear(self.d_FC_layer, self.d_FC_layer))
                self.predict.append(nn.Dropout(self.dropout))
                self.predict.append(nn.LeakyReLU())
                self.predict.append(nn.BatchNorm1d(d_FC_layer))

    def forward(self, h):
        for layer in self.predict:
            h = layer(h)

        return h
