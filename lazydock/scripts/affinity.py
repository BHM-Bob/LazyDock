import argparse
import os
from typing import Callable, Dict, List, Set, Tuple

import pandas as pd
import torch
from mbapy_lite.base import put_err
from mbapy_lite.file import get_paths_with_extension
from mbapy_lite.web_utils.task import TaskPool
from pymol import cmd
from tqdm import tqdm

try:
    from dgl.multiprocessing import Queue
except:
    from multiprocessing import Queue

from lazydock.nn.EHIGN import (build_model, evaluate, graphs_to_batch,
                               load_checkpoint, pdbstr_to_graph)
from lazydock.scripts._script_utils_ import (Command, make_args_and_excute,
                                             process_batch_dir_lst)


def mp_pool_init_fn(q: Queue):
    global queue
    queue = q

def complex2graph(pdb_path, args):
    cmd.reinitialize()
    cmd.load(pdb_path, 'complex')
    if not args.receptor_chain:
        rec_chains = list(set(cmd.get_chains()) - set([args.ligand_chain]))
    else:
        rec_chains = args.receptor_chain
    if cmd.select('receptor', f'complex and (' + ' or '.join([f'chain {c}' for c in rec_chains]) + ')') == 0:
        put_err(f'{pdb_path}: can not find receptor with chain {rec_chains}')
        return None
    try:
        graph = pdbstr_to_graph(cmd.get_pdbstr(f'chain {args.ligand_chain}'),
                                receptor_pdbstr=cmd.get_pdbstr('receptor'), dis_threshold=args.distance)
    except Exception as e:
        put_err(f'{pdb_path}: error in pdbstr_to_graph: {e}')
        return None
    if not graph:
        put_err(f'{pdb_path}: graph is None')
        return None
    return queue.put([pdb_path, graph])


class EHIGN(Command):
    HELP = """perform affinity prediction on a complex.pdb file with EHIGN,
EHIGN is from 
Yang, Z., Zhong, W., Lv, Q., Dong, T., Chen, G., and Chen, C.Y.-C. (2024). Interaction-Based Inductive Bias in Graph Neural Networks: Enhancing Protein-Ligand Binding Affinity Predictions From 3D Structures. IEEE Transactions on Pattern Analysis and Machine Intelligence, 1-18. 10.1109/TPAMI.2024.3400515.
"""
    def __init__(self, args: argparse.Namespace, printf=print) -> None:
        super().__init__(args, printf, ['batch_dir'])
        self.tasks = []
        self.df = pd.DataFrame(columns=['path', 'rel_path', "affinity"])
        self.df.set_index('path', inplace=True)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('-n', '--name', type = str, required=True,
                          help="name for input complex file, such as `complex.pdb`.")
        args.add_argument('-lc', '--ligand-chain', type = str, required=True,
                          help=f"ligand chain.")
        args.add_argument('-rc', '--receptor-chain', type = str, nargs='+', default=None,
                          help="receptor chain, if not specified, will be others except ligand chain.")
        args.add_argument('-ckp', '--check-point', type = str, required=True,
                          help="checkpoint file, default is %(default)s.")
        args.add_argument('-bs', '--batch-size', type = int, default=128,
                          help='batch size, default is %(default)s. for EHIGN.')
        args.add_argument('-dist', '--distance', type = float, default=5,
                          help='distance threshold for pocket, default is %(default)s.')
        args.add_argument('-o', '--output', type = str, default='EHIGN_score.csv',
                          help="output file, default is %(default)s.")
        args.add_argument('-nw', '--n-workers', type = int, default=1,
                          help="number of workers for structure pre-process, default is %(default)s.")
        args.add_argument('-device', type = str, default='cuda:0',
                          help="device for model inference, default is %(default)s.")
        return args

    def process_args(self):
        # process IO
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        
    def predict_one_batch(self, model, graph_lst):
        paths, graphs = zip(*graph_lst)
        graphs = graphs_to_batch(graphs)
        with torch.no_grad():
            affinity1, affinity2 = evaluate(model, graphs, device=self.args.device)
            affinity = (affinity1 + affinity2) / 2
        affinity = affinity.cpu().numpy().tolist()
        for path, affi in zip(paths, affinity):
            self.df.loc[path] = [os.path.relpath(path, self.args.batch_dir), affi]  # type: ignore
        
    def main_process(self):
        # search for pdb files
        paths = get_paths_with_extension(self.args.batch_dir, ['.pdb'],
                                         name_substr=self.args.name, sort='natsort')
        if not paths:
            put_err(f'can not find any pdb file with name {self.args.name} in {self.args.batch_dir}')
            return
        # if no ckp exist, use default ckp
        if not os.path.exists(self.args.check_point):
            put_err(f'can not find checkpoint file {self.args.check_point}')
            return
        # load model
        model = build_model()
        load_checkpoint(model, self.args.check_point, device=self.args.device)
        # submit and run tasks parallel
        queue, graph_lst = Queue(), []
        pool = TaskPool('process', self.args.n_workers, report_error=True,
                        mp_pool_init_kwargs={'initializer': mp_pool_init_fn, 'initargs': (queue,)}).start()
        for path in tqdm(paths):
            pool.add_task(path, complex2graph, path, self.args)
            while not queue.empty():
                pdb_path, graph = queue.get()
                graph_lst.append([pdb_path, graph])
                if len(graph_lst) >= self.args.batch_size:
                    self.predict_one_batch(model, graph_lst)
                    graph_lst.clear()
            pool.wait_till_free()
        pool.wait_till_all_done()
        self.predict_one_batch(model, graph_lst)
        graph_lst.clear()
        self.df.to_csv(os.path.join(self.args.batch_dir, self.args.output), index=True)
        pool.close(1)


_str2func = {
    'ehign': EHIGN,
}


def main(sys_args: List[str] = None):
    torch.multiprocessing.set_sharing_strategy('file_system')
    make_args_and_excute('tools for predict affinity', _str2func, sys_args)


if __name__ == "__main__":    
    main()