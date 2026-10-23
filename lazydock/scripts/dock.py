'''
Date: 2024-12-04 20:58:39
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-18 16:01:18
Description: 
'''
import argparse
import multiprocessing
import os
import traceback
from pathlib import Path
from typing import Dict, List, Tuple, Union

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from mbapy_lite.base import Configs, put_err, put_log
from mbapy_lite.file import get_paths_with_extension, opts_file, write_sheets
from mbapy_lite.plot import save_show
from mbapy_lite.web import TaskPool, random_sleep
from tqdm import tqdm

from lazydock.dock.config import (MACROCYCLE_MAX_RING_SIZE,
                                  RING_CLOSURE_CLOSED_A)
from lazydock.dock.dock_adgpu import (perform_autodock_gpu_dock,
                                      perform_autodock_gpu_flex_dock)
from lazydock.dock.dock_vina import perform_vina_dock, perform_vina_flex_dock
from lazydock.dock.flex_utils import (build_flex_complex_pdbstr,
                                      parse_flex_residues)
from lazydock.dock.meeko_utils import measure_ring_closure, strip_glue_atoms
from lazydock.pml.autodock_utils import DlgFile
from lazydock.scripts._script_utils_ import (Command, check_memory_usage,
                                             clean_path, excute_command,
                                             process_batch_dir_lst)
from lazydock.web.dinc import run_dock_on_DINC_ensemble
from lazydock.web.hdock import run_dock_on_HDOCK, run_dock_on_HPEPDOCK

# 全局 GPU 槽位队列: 主进程写入可用 gpu_id, 子进程取走后用完归还。
# 某卡并发已满时队列自然不会向它放新的 gpu_id, 避免长任务卡拖满任务池。
_GPU_SLOT_QUEUE: 'multiprocessing.Queue' = None

def mp_pool_init_fn(slot_queue):
    """进程池子进程初始化: 继承主进程创建的 GPU 槽位队列"""
    global _GPU_SLOT_QUEUE
    _GPU_SLOT_QUEUE = slot_queue

def perform_autodock_gpu_dock_on_slot(*args, **kwargs):
    """从槽位队列领一个 gpu_id 执行 docking, 无论成败归还槽位"""
    return perform_dock_on_slot(perform_autodock_gpu_dock, *args, **kwargs)

def perform_dock_on_slot(dock_fn, *args, **kwargs):
    """从槽位队列领一个 gpu_id 执行 docking, 无论成败归还槽位"""
    global _GPU_SLOT_QUEUE
    if _GPU_SLOT_QUEUE is None:
        return dock_fn(*args, **kwargs)
    gpu_id = _GPU_SLOT_QUEUE.get()
    try:
        return dock_fn(*args, gpu_id=gpu_id, **kwargs)
    except Exception:
        traceback.print_exc()
    finally:
        _GPU_SLOT_QUEUE.put(gpu_id)


class vina(Command):
    def __init__(self, args, printf = print):
        super().__init__(args, printf, ['batch_dir'])
        self.taskpool = None
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-c', '--config-name', type = str, default='config.txt',
                                help='config file name, each config is a task. Default is %(default)s.')
        args.add_argument('-v', '--vina-name', type = str, default='vina',
                                help='vina executable name to call. Default is %(default)s.')
        args.add_argument('--vina-args', type = str, default='--log ./log.txt',
                                help='args for vina executable. Default is %(default)s.')
        args.add_argument('-n', '--n-workers', type=int, default=1,
                                help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.args.vina_args = '' if self.args.vina_args in {"''", '""', ''} else self.args.vina_args
        
    @staticmethod
    def run_vina(config_path: Path, vina_name: str, vina_args: str):
        print(f'current: {config_path}')
        config = opts_file(config_path, way='lines')
        try:
            out_name = list(filter(lambda x: x.startswith('out'), config))[0].split('=')[1].strip() # type: ignore
        except:
            put_err(f'config {config_path} has no out_name: {config}.')
            return 
        if (config_path.parent / out_name).exists():
            print(f'{config_path.parent} has done, skip')
            return 
        cmd_string = f'cd "{config_path.parent}" && {vina_name} --config ./{config_path.name} {vina_args}'
        print(f'running: {cmd_string}')
        os.system(cmd_string)
        
    def main_process(self):
        if not check_memory_usage(self.args.n_workers):
            return put_log('aborted by user.')
        configs_path = get_paths_with_extension(self.args.batch_dir, ['.txt'],
                                                name_substr=self.args.config_name, sort='natsort')
        print(f'get {len(configs_path)} config(s) for docking')
        self.taskpool = TaskPool('threads', self.args.n_workers, report_error=True).start()
        for config_path in tqdm(configs_path, total=len(configs_path)):
            self.taskpool.add_task(None, self.run_vina, Path(config_path), self.args.vina_name, self.args.vina_args)
            self.taskpool.wait_till_free()
        self.taskpool.wait_till_all_done()
        self.taskpool.close(1)


def check_macrocycle_args(args: argparse.Namespace):
    """reject --macrocycle combined with the ADT preparation backend.

    meeko is the only backend that can break a ring: ADT's `prepare_ligand` has no
    macrocycle option at all, so `--prepare-backend adfr --macrocycle` would hand
    back a ligand whose ring is still completely frozen while the rest of the
    pipeline behaved as though it were free.  Rejecting the combination is better
    than silently honouring only one of the two flags.

    Note what this check is NOT: it does not care whether the ENGINE can use glue.
    AutoDock-GPU can -- its AD4 glue term is a linear attraction with no cutoff, the
    retyped CGx atoms fall back on the carbon affinity map, and Gx is excluded from
    the intermolecular sum -- so the GPU subcommands accept --macrocycle too.
    """
    if not getattr(args, 'macrocycle', False):
        return
    if getattr(args, 'prepare_backend', 'meeko') == 'adfr':
        put_err('--macrocycle requires meeko: ADT has no ring-breaking option, so '
                '`--prepare-backend adfr` would deliver a ligand whose ring is still '
                'frozen. Drop --macrocycle, or drop `--prepare-backend adfr`. Exit.',
                _exit=True)


def check_gpu_spacing(args: argparse.Namespace, engine: str = 'redock-gpu'):
    """reject a non-default --spacing on every AutoDock-GPU subcommand.

    Both GPU paths read the grid spacing back from the gpf/fld, and neither can be
    handed another value:

    - `redock-gpu` builds its gpf with `prepare_gpf.py` + `autogrid4`.  This project
      does not forward `spacing` to the gpf, so the file keeps prepare_gpf.py's own
      0.375 default while `npts` was computed with the requested value -- the box
      would be sized by one spacing and described by another.
    - `redock-flex-gpu` takes its maps from the gpf that meeko writes, and meeko fixes
      the spacing at 0.375 A there (`meeko/gridbox.py`, `get_gpf_string`) with no CLI
      knob to change it.

    Failing loudly is better than accepting a value that has no effect.
    """
    spacing = getattr(args, 'spacing', None)
    if spacing is None or abs(spacing - 0.375) < 1e-9:
        return
    # `sub_command` is the name the user actually typed, which matters because
    # redock-flex-gpu inherits the call site from redock-gpu
    engine = getattr(args, 'sub_command', None) or engine
    put_err(f'--spacing {spacing} can not be honoured by {engine}: the AutoDock-GPU '
            f'paths take the grid spacing from the gpf, which is fixed at 0.375 A '
            f'(Angstrom). Use the vina-based subcommands if another spacing is '
            f'required. Exit.', _exit=True)


class redock(Command):
    HELP = """perform redock through vina or other method"""
    def __init__(self, args: argparse.Namespace, printf=print) -> None:
        super().__init__(args, printf, ['batch_dir'])
        self.tasks = []
        # closure is the ring-closure distance in ANGSTROM (see
        # `meeko_utils.measure_ring_closure`); None for a ligand without glue pairs
        self.summary_df = pd.DataFrame(columns=['path', 'rel_path', 'energy', 'closure'])
        self.summary_df.set_index('path', inplace=True)
        
    @staticmethod
    def make_common_args(args: argparse.ArgumentParser):
        """add every argument the redock subcommands share.

        The vina and the AutoDock-GPU subcommands need the same input selection,
        output naming, grid box, macrocycle switch and preparation backend; only the
        search knobs differ, and those are added by each subclass.  Declaring them
        once here (instead of a set of `add_*_args` helpers called from several
        `make_args`) is what keeps the two engine families from drifting apart: a
        subclass gets the shared set by delegating, not by re-listing it.
        """
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help=f"dir which contains many sub-folders, each sub-folder contains docking result files.")
        args.add_argument('-n', '--name', type = str, required=True,
                          help="name for input complex file, such as `complex.pdb`.")
        args.add_argument('-lc', '--ligand-chain', type = str, required=True,
                          help=f"ligand chain.")
        args.add_argument('-rc', '--receptor-chain', type = str, nargs='+', default=None,
                          help="receptor chain, if not specified, will be others except ligand chain.")
        args.add_argument('-gb', '--grid-buffer', type = float, default=4,
                          help='grid box padding in ANGSTROM, added on every side of the ligand '
                               '(and of every flexible residue), default is %(default)s.')
        args.add_argument('--spacing', type = float, default=0.375,
                          help='grid spacing in ANGSTROM, defines how many grid points a given '
                               'box size becomes (npts = size / spacing); it must match the '
                               'spacing the engine is run with, default is %(default)s.')
        args.add_argument('--macrocycle', action = 'store_true', default=False,
                          help='prepare the ligand with meeko so that one bond of a large ring is '
                               'broken open and the ring skeleton can adapt during the search. '
                               'Both engines pull it back shut: vina >= 1.2 and AutoDock-GPU both '
                               'score a linear glue attraction. Ligands prepared this way are not '
                               'score-comparable with rigid-ring ones, default is %(default)s.')
        args.add_argument('--max-ring-size', type = int, default=MACROCYCLE_MAX_RING_SIZE,
                          help='largest ring meeko may break. Meeko silently ignores bigger rings '
                               'and leaves them rigid, so this must cover the real ring size '
                               '(cyclic peptides reach 56 atoms), default is %(default)s.')
        args.add_argument('--prepare-backend', type = str, default='meeko',
                          choices=['meeko', 'adfr'],
                          help='toolkit that prepares the receptor and ligand pdbqt. meeko and '
                               'adfr type a few nitrogens/polar hydrogens differently, so scores '
                               'from the two are not comparable. Can not be combined with '
                               '--macrocycle, because ADT has no ring-breaking option. Note '
                               '`redock-flex-*` always splits the receptor with meeko and this '
                               'switch then only selects the ligand backend, default is %(default)s.')
        args.add_argument('-np', '--n-poses', type = int, default=20,
                          help='number of poses, default is %(default)s.')
        args.add_argument('--suffix', type = str, default='_redock',
                          help="suffix for the output file, default is %(default)s.")
        args.add_argument('--n-pdb', type = int, default=1,
                          help="number of pdb files to extract from the result, sorted by score, "
                               "such as a_redock_0.pdb, default is %(default)s.")
        args.add_argument('--summary', type = str, default='redock_summary.csv',
                          help='energy summary table, default is %(default)s.')
        args.add_argument('-nw', '--n-workers', type = int, default=1,
                          help="number of workers, default is %(default)s.")
        return args

    @staticmethod
    def make_flex_args(args: argparse.ArgumentParser):
        """add the flexible-residue arguments shared by both flex subcommands.

        Declared next to `make_common_args` and for the same reason: `redock-flex-gpu`
        and `redock-flex-vina` sit on different branches of the inheritance tree, so a
        single definition on `redock` is what keeps their flex flags identical.

        The three sources are combined, not exclusive: an explicit residue list,
        residues within a distance of the ligand, and residues that fplip finds
        interacting with the ligand.
        """
        args.add_argument('-fr', '--flex-residues', type = str, nargs='+', default=None,
                          help="flexible receptor residues, could be '<chain>:<resnum>(,<resnum>)' such as 'A:68,50', "
                               "or '<chain>:<RESN><resnum>' such as 'A:LYS68,LYS50'.")
        args.add_argument('-fd', '--flex-distance', type = float, default=None,
                          help='select receptor residues within this distance of any ligand heavy atom, default is %(default)s.')
        args.add_argument('-fi', '--flex-interaction', action = 'store_true', default=False,
                          help='select receptor residues which have interaction with ligand by fplip, default is %(default)s.')
        args.add_argument('--flex-interaction-mode', type = str, default='all',
                          help='interaction types for fplip, "all" to include all modes, default is %(default)s.')
        args.add_argument('--flex-interaction-cutoff', type = float, default=4.,
                          help='interaction distance cutoff for fplip, default is %(default)s.')
        return args

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = redock.make_common_args(args)
        args.add_argument('-sf', '--score-name', default='vina', choices=['vina', 'vinardo', 'ad4'],
                          help='Score name for vina, default is %(default)s.')
        args.add_argument('-e', '--exhaustiveness', type = int, default=32,
                          help='exhaustiveness, default is %(default)s.')
        return args

    def process_args(self):
        # process IO
        self.args.batch_dir = process_batch_dir_lst(self.args.batch_dir)
        check_macrocycle_args(self.args)
    
    def process_result(self, result_path: str, result):
        if not isinstance(result, tuple) or len(result) < 5:
            return put_err(f'invalid result: {result}', False)
        # redock/redock-gpu 返回 5 元组, redock-flex-gpu 返回 6 元组(多一个 rigid pdbqt)
        result_pdbstr, rec_pdbstr = result[0], result[1]
        rigid_pdbqtstr = result[2] if len(result) >= 6 else None
        opts_file(result_path, 'w', way='str', data = result_pdbstr)
        dlg = DlgFile(None, result_pdbstr, True, True)
        dlg.sort_pose() # sort by docking energy
        for i, pose in enumerate(dlg.pose_lst[:self.args.n_pdb]):
            path = Path(result_path)
            pdb_path = str(path.parent / (path.stem + f'_{i}.pdb'))
            if rigid_pdbqtstr is None:
                rec_pdbstr_ = rec_pdbstr
            else:
                # flexible docking: 柔性残基的 CA 与侧链已在 pose 中,
                # 受体需去掉这些原子(含其上的氢), 避免悬空氢
                rec_pdbstr_ = build_flex_complex_pdbstr(rec_pdbstr, rigid_pdbqtstr, pose.as_pdb_string())
            rec_pdbstr_ = rec_pdbstr_.replace('END', '')
            # 大环对接的 glue 伪原子会在 parse2std 下变成假元素 G, 下游 PDBFixer/OpenMM
            # 会因此失败, 必须在落盘前剔除
            pose_pdbstr = strip_glue_atoms(pose.as_pdb_string())
            # avoid null line, which will cause openbabel error
            opts_file(pdb_path, 'w', way='str', data = f'{rec_pdbstr_}\nTER\n' + pose_pdbstr.strip('TER\n').strip('\n'))
            # 环闭合度只能从原始 pdbqt 量: parse2std 会把 CG0/G0 的 type 改写成 C/G,
            # 配对信息在那之后就没有了
            closure = measure_ring_closure(pose.info_string())
            if self.args.macrocycle and closure is not None and closure > RING_CLOSURE_CLOSED_A:
                put_log(f'{pdb_path}: ring closure {closure:.2f} A > '
                        f'{RING_CLOSURE_CLOSED_A} A, the macrocycle did not close in this pose')
            self.summary_df.loc[pdb_path] = [os.path.relpath(pdb_path, self.args.batch_dir),
                                             pose.energy, closure]
        return True
        
    def main_process(self):
        if not check_memory_usage(self.args.n_workers):
            return put_log('aborted by user.')
        # search for pdb files
        paths = get_paths_with_extension(self.args.batch_dir, ['.pdb'],
                                         name_substr=self.args.name, sort='natsort')
        if not paths:
            put_err(f'can not find any pdb file with name {self.args.name} in {self.args.batch_dir}')
            return
        # submit and run tasks parallel
        pool = TaskPool('process', self.args.n_workers, report_error=True,
                        mp_pool_init_kwargs={'maxtasksperchild': 100}).start()
        for path in tqdm(paths):
            path = Path(path)
            result_path = str(path.parent / (path.stem + f'{self.args.suffix}.pdbqt'))
            pool.add_task(result_path, perform_vina_dock, self.args.score_name, self.args.grid_buffer, opts_file(path),
                          self.args.ligand_chain, self.args.receptor_chain, self.args.exhaustiveness, self.args.n_poses,
                          spacing=self.args.spacing,
                          macrocycle=self.args.macrocycle, max_ring_size=self.args.max_ring_size,
                          prepare_backend=self.args.prepare_backend)
            while True:
                result = pool.pull_task(return_name=True, block=False)
                if not result or result == pool.TASK_NOT_RETURNED or result == pool.TASK_NOT_SUCCEEDED:
                    break
                result_path, result = result
                self.process_result(result_path, result)
                self.summary_df.to_csv(os.path.join(self.args.batch_dir, self.args.summary), index=True)
            pool.wait_till_free()
        pool.wait_till_all_done()
        # retrieve results
        for result_path in list(pool.tasks.keys()):
            result = pool.query_task(result_path, True, 30) # type: ignore
            self.process_result(result_path, result)
        self.summary_df.to_csv(os.path.join(self.args.batch_dir, self.args.summary), index=True)


class redock_gpu(redock):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = redock.make_common_args(args)
        # 共用参数的语义在 GPU 上不同, 就地覆盖 help
        for act in args._actions:
            if act.dest == 'n_poses':
                act.help = ('number of LGA runs, each yielding one pose, '
                            'default is %(default)s.')
            elif act.dest == 'suffix':
                act.help = "suffix for output dlg file, default is %(default)s."
        args.add_argument('--seed', type = int, default=0,
                          help='seed for random number, default is %(default)s.')
        args.add_argument('--nev', type = int, default=2500000,
                          help='Score evaluations (max.) per LGA run, default is %(default)s.')
        args.add_argument('--stopstd', type = float, default=0.15,
                          help='AutoStop energy standard deviation tolerance, default is %(default)s.')
        args.add_argument('--heurmax', type = int, default=12000000,
                          help='Asymptotic heuristics evals limit (smooth limit), default is %(default)s.')
        args.add_argument('--heuristics', type = int, default=1, choices=[0, 1],
                          help='1 (the engine default) lets AutoDock-GPU derive the search method '
                               'and the evaluation budget from the ligand torsion count; 0 makes '
                               '--nev the exact budget. AutoStop can not stop before 100%% of the '
                               'budget while this is 1, so 0 is what makes --stopstd effective, '
                               'default is %(default)s.')
        args.add_argument('--ngen', type = int, default=42000,
                          help='maximum generations per LGA run, AutoDock-GPU default is '
                               '%(default)s; raise it when --nev rather than --ngen should be '
                               'what ends the run')
        args.add_argument('--autostop', type = int, default=1, choices=[0, 1],
                          help='stop a run once the best energies converge, default is %(default)s.')
        args.add_argument('--gpus', type = int, nargs='+', default=[0],
                          help="gpu ids, default is %(default)s.")
        args.add_argument('--n-task-per-gpu', type = int, default=1,
                          help="number of tasks per gpu, default is %(default)s.")
        return args

    def process_args(self):
        super().process_args()
        # --macrocycle 在本路径可用: AD4 的 glue 项是无 cutoff 的线性吸引, CGx 会
        # 退回碳图, Gx 不计入分子间求和, 所以这里既有的 map 不需要任何改动.
        # 仍然拒绝 spacing 的非默认值, 它由 gpf 决定, 兑现不了
        check_gpu_spacing(self.args, 'redock-gpu')

    def main_process(self):
        # search for pdb files
        paths = get_paths_with_extension(self.args.batch_dir, ['.pdb'], name_substr=self.args.name, sort='natsort')
        if not paths:
            put_err(f'can not find any pdb file with name {self.args.name} in {self.args.batch_dir}')
            return
        self._run_dock_pool(paths,
                            lambda p: (perform_autodock_gpu_dock_on_slot, opts_file(p),
                                       self.args.ligand_chain, self.args.receptor_chain, self.args.grid_buffer),
                            dict(seed=self.args.seed, nrun=self.args.n_poses, nev=self.args.nev,
                                 stopstd=self.args.stopstd, heurmax=self.args.heurmax,
                                 heuristics=self.args.heuristics, ngen=self.args.ngen,
                                 autostop=self.args.autostop,
                                 macrocycle=self.args.macrocycle,
                                 max_ring_size=self.args.max_ring_size,
                                 prepare_backend=self.args.prepare_backend))

    def _run_dock_pool(self, paths, task_args_fn, task_kwargs: dict):
        """通用 GPU 任务池: 槽位队列 + TaskPool 提交/回收 + 结果处理 + summary 落盘"""
        # GPU 槽位队列: 展开为 [g0 x n_task_per_gpu, g1 x n_task_per_gpu, ...],
        # 子进程取走一个 gpu_id 执行任务, 无论成败都归还, 天然形成背压负载均衡:
        # 某卡并发已满时队列没有它的空闲槽位, 长任务卡不会再被塞入新任务
        slot_queue = multiprocessing.Queue()
        for g in self.args.gpus:
            for _ in range(max(1, self.args.n_task_per_gpu)):
                slot_queue.put(g)
        # submit and run tasks parallel
        ## adfr use pymol.cmd to optrate pdb, so use process to isolate each task
        pool = TaskPool('process', self.args.n_workers, report_error=True,
                        mp_pool_init_kwargs={'maxtasksperchild': 100, 'initializer': mp_pool_init_fn, 'initargs': (slot_queue,)}).start()
        # launch tasks
        for path in tqdm(paths, total=len(paths)):
            path = Path(path)
            result_path = str(path.parent / (path.stem + f'{self.args.suffix}.dlg'))
            pool.add_task(result_path, *task_args_fn(path), **task_kwargs)
            while True:
                result = pool.pull_task(return_name=True, block=False)
                if not result or result == pool.TASK_NOT_RETURNED or result == pool.TASK_NOT_SUCCEEDED:
                    break
                result_path, result = result
                self.process_result(result_path, result)
                self.summary_df.to_csv(os.path.join(self.args.batch_dir, self.args.summary), index=True)
            pool.wait_till_free()
            pool.wait_till(lambda : not slot_queue.empty(), 0.01)
        pool.wait_till_all_done()
        # retrieve results
        for result_path in list(pool.tasks.keys()):
            result = pool.query_task(result_path, True, 30) # type: ignore
            self.process_result(result_path, result)
        self.summary_df.to_csv(os.path.join(self.args.batch_dir, self.args.summary), index=True)
def check_flex_args(args: argparse.Namespace):
    """validate flexible-residue arguments shared by all redock-flex-* subcommands"""
    # 至少需要一种柔性残基来源
    if args.flex_residues is None and args.flex_distance is None and not args.flex_interaction:
        put_err('at least one of --flex-residues, --flex-distance, --flex-interaction should be provided, exit.', _exit=True)
    # 提前校验格式, 避免每个任务重复报错
    if args.flex_residues is not None:
        if not parse_flex_residues(args.flex_residues):
            put_err(f'invalid --flex-residues: {args.flex_residues}, exit.', _exit=True)


class redock_flex_gpu(redock_gpu):
    """redock with AutoDock-GPU, receptor flexible residues are supported"""
    HELP = """redock with AutoDock-GPU and flexible receptor residues"""
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = redock_gpu.make_args(args)
        # flex 版本使用独立默认值, 避免与 redock/redock-gpu 的输出互相覆盖
        for act in args._actions:
            if act.dest == 'suffix':
                act.default = '_redock_flex'
                act.help = "suffix for output dlg file, default is %(default)s."
            elif act.dest == 'summary':
                act.default = 'redock_flex_summary.csv'
        return redock.make_flex_args(args)

    def process_args(self):
        super().process_args()
        check_flex_args(self.args)
        # --spacing 的拒绝由 redock_gpu.process_args 完成(此处继承), 不再重复检查
        
    def main_process(self):
        # search for pdb files
        paths = get_paths_with_extension(self.args.batch_dir, ['.pdb'], name_substr=self.args.name, sort='natsort')
        if not paths:
            put_err(f'can not find any pdb file with name {self.args.name} in {self.args.batch_dir}')
            return
        self._run_dock_pool(
            paths,
            lambda p: (perform_dock_on_slot, perform_autodock_gpu_flex_dock, opts_file(p),
                       self.args.ligand_chain, self.args.receptor_chain, self.args.grid_buffer),
            dict(flex_residues=self.args.flex_residues,
                 flex_distance=self.args.flex_distance,
                 flex_interaction=self.args.flex_interaction,
                 flex_interaction_mode=self.args.flex_interaction_mode,
                 flex_interaction_cutoff=self.args.flex_interaction_cutoff,
                 seed=self.args.seed, nrun=self.args.n_poses, nev=self.args.nev,
                 stopstd=self.args.stopstd, heurmax=self.args.heurmax,
                 heuristics=self.args.heuristics, ngen=self.args.ngen,
                 autostop=self.args.autostop,
                 macrocycle=self.args.macrocycle,
                 max_ring_size=self.args.max_ring_size,
                 prepare_backend=self.args.prepare_backend))


class redock_flex_vina(redock):
    """redock with vina, receptor flexible residues are supported"""
    HELP = """redock with vina and flexible receptor residues
    
NOTICE: vina searches on CPU and flexible sidechains add many degrees of freedom,
so this is much slower than `redock-flex-gpu`.  Prefer `redock-flex-gpu` unless
the vina scoring function is specifically required."""
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = redock.make_args(args)
        # flex 版本使用独立默认值, 避免与 redock 的输出互相覆盖
        for act in args._actions:
            if act.dest == 'suffix':
                act.default = '_redock_flex'
                act.help = "suffix for output pdbqt file, default is %(default)s."
            elif act.dest == 'summary':
                act.default = 'redock_flex_summary.csv'
        redock.make_flex_args(args)
        args.add_argument('--seed', type = int, default=0,
                          help='seed for random number, default is %(default)s.')
        args.add_argument('--max-evals', type = int, default=0,
                          help='max evaluations, 0 means use vina heuristics, default is %(default)s.')
        args.add_argument('--min-rmsd', type = float, default=1,
                          help='minimum RMSD difference between poses, default is %(default)s.')
        args.add_argument('--cpu', type = int, default=0,
                          help='cpu number for vina, 0 means use all cpus, default is %(default)s.')
        return args

    def process_args(self):
        super().process_args()
        check_flex_args(self.args)

    def main_process(self):
        # search for pdb files
        paths = get_paths_with_extension(self.args.batch_dir, ['.pdb'], name_substr=self.args.name, sort='natsort')
        if not paths:
            put_err(f'can not find any pdb file with name {self.args.name} in {self.args.batch_dir}')
            return
        # submit and run tasks parallel
        pool = TaskPool('process', self.args.n_workers, report_error=True,
                        mp_pool_init_kwargs={'maxtasksperchild': 100}).start()
        for path in tqdm(paths, total=len(paths)):
            path = Path(path)
            result_path = str(path.parent / (path.stem + f'{self.args.suffix}.pdbqt'))
            pool.add_task(result_path, perform_vina_flex_dock,
                          self.args.score_name, self.args.grid_buffer, opts_file(path),
                          self.args.ligand_chain, self.args.receptor_chain,
                          flex_residues=self.args.flex_residues,
                          flex_distance=self.args.flex_distance,
                          flex_interaction=self.args.flex_interaction,
                          flex_interaction_mode=self.args.flex_interaction_mode,
                          flex_interaction_cutoff=self.args.flex_interaction_cutoff,
                          exhaustiveness=self.args.exhaustiveness, n_poses=self.args.n_poses,
                          min_rmsd=self.args.min_rmsd, max_evals=self.args.max_evals,
                          seed=self.args.seed, cpu=self.args.cpu,
                          spacing=self.args.spacing,
                          macrocycle=self.args.macrocycle, max_ring_size=self.args.max_ring_size,
                          prepare_backend=self.args.prepare_backend)
            while True:
                result = pool.pull_task(return_name=True, block=False)
                if not result or result == pool.TASK_NOT_RETURNED or result == pool.TASK_NOT_SUCCEEDED:
                    break
                result_path, result = result
                self.process_result(result_path, result)
                self.summary_df.to_csv(os.path.join(self.args.batch_dir, self.args.summary), index=True)
            pool.wait_till_free()
        pool.wait_till_all_done()
        # retrieve results
        for result_path in list(pool.tasks.keys()):
            result = pool.query_task(result_path, True, 30) # type: ignore
            self.process_result(result_path, result)
        self.summary_df.to_csv(os.path.join(self.args.batch_dir, self.args.summary), index=True)


def hdock_run_fn_warpper(result_prefix: str = 'HDOCK', result_name: str = 'HDOCK_all_results.tar.gz'):
    def ret_warpper(func):
        def core_wrapper(*args, **kwargs):
            config_path = args[0] if len(args) > 0 else kwargs.get('config_path', None)
            if config_path is None:
                return put_err('config_path is required, skip.')
            # get parameters from config file
            if isinstance(config_path, Path):
                root = config_path.parent
                parameters = hdock.get_paramthers_from_config(config_path)
                parameters['receptor_path'] = config_path.parent / parameters['receptor']
                parameters['ligand_path'] = config_path.parent / parameters['ligand']
            # OR get parameters from receptor and ligand path in tuple
            elif isinstance(config_path, tuple):
                root = Path(config_path[0]).parent
                parameters = {'receptor_path': config_path[0], 'ligand_path': config_path[1]}
            # un-expected type
            else:
                return put_err(f'config_path type not support: {type(config_path)}, skip.')
            # check if done
            if (root / f'{result_prefix}_all_results.tar.gz').exists() or (root / result_name).exists():
                return print(f'{root} has done, skip')
            # perform docking
            print(f'current: {root}')
            ret =  func(*args, parameters=parameters, **kwargs)
            return ret
        return core_wrapper
    return ret_warpper


class hdock(vina):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
    
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-r', '--receptor', type = str, default=None,
                          help="receptor pdb file name, optional. If provided, will ignore config.txt.")
        args.add_argument('-l', '--ligand', type = str, default=None,
                          help="ligand pdb file name, optional. If provided, will ignore config.txt.")
        args.add_argument('--config-name', type=str, default=None,
                          help='Vina config for specifing parameters. Default is %(default)s.')
        args.add_argument('-m', '--method', type = str, default='web', choices=['web'],
                          help='docking method. Currently support "web". Default is %(default)s.')
        args.add_argument('--email', type = str, default=None,
                          help='email address for HDOCK web server. Default is %(default)s.')
        args.add_argument('-gui', '--gui', action='store_true', default=False,
                          help='show browser GUI. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
    
    @staticmethod
    def get_paramthers_from_config(config_path: Path) -> Dict:
        config_lines = config_path.read_text().split('\n')
        paramthers = {line.split('=')[0].strip():line.split('=')[1].strip() for line in config_lines if not line.startswith('#')}
        return {k:v[:v.find('#')] for k,v in paramthers.items()}

    @staticmethod
    @hdock_run_fn_warpper('HDOCK')
    def run_hdock_web(config_path: Union[Path, Tuple[str, str]], parameters: Dict[str, str] = None, email=None, browser=None,
                      args: argparse.ArgumentParser = None):
        w_dir = config_path.parent if isinstance(config_path, Path) else Path(config_path[0]).parent
        run_dock_on_HDOCK(receptor_path=parameters['receptor_path'], ligand_path=parameters['ligand_path'],
                          w_dir=w_dir, email=email, browser=browser)

    @staticmethod
    @hdock_run_fn_warpper('HDOCK')
    def run_hdock_local(config_path: Union[Path, Tuple[str, str]], parameters: Dict[str, str] = None, **kwargs):
        raise NotImplementedError('local docking not implemented yet.')

    def main_process(self):
        if not os.path.isdir(self.args.batch_dir):
            return put_err(f'dir argument should be a directory: {self.args.config}.')
        if self.args.receptor is not None and self.args.ligand is not None:
            r_paths = get_paths_with_extension(self.args.batch_dir, [],
                                               name_substr=self.args.receptor, sort='natsort')
            l_paths = get_paths_with_extension(self.args.batch_dir, [],
                                               name_substr=self.args.ligand, sort='natsort')
            if len(r_paths) != len(l_paths):
                r_roots = [os.path.dirname(p) for p in r_paths]
                l_roots = [os.path.dirname(p) for p in l_paths]
                roots_count = {root: r_roots.count(root)+l_roots.count(root) for root in (set(r_roots) | set(l_roots))}
                invalid_roots = '\n'.join([root for root, count in roots_count.items() if count != 2])
                return put_err(f"The number of receptor and ligand files is not equal, please check the input files.\ninvalid roots:\n{invalid_roots}")
            configs_path = [(r, l) for r, l in zip(r_paths, l_paths)]
        elif self.args.config_name is not None:
            configs_path = get_paths_with_extension(self.args.batch_dir, ['.txt'],
                                                    name_substr=self.args.config_name, sort='natsort')
        else:
            return put_err('config_name or receptor and ligand should be provided, skip.')
        print(f'get {len(configs_path)} config(s) for docking')
        # allow browser gui
        if self.args.gui:
            from mbapy_lite.web import Browser
            browser = Browser(options=[f"--user-agent={Configs.web.chrome_driver_path}"])
        else:
            browser = None
        # perform docking
        dock_fn = getattr(self, f'run_hdock_{self.args.method}')
        for config_path in tqdm(configs_path, total=len(configs_path)):
            dock_fn(Path(config_path).resolve() if isinstance(config_path, str) else config_path,
                    email=self.args.email, browser=browser, args=self.args)
            random_sleep(300, 180) # sleep 3~5 minutes to avoid overloading the server


class hpepdock(hdock):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = hdock.make_args(args)
        # make sure hpepdock only support web method
        method_arg_idx = args._actions.index(list(filter(lambda x: x.dest =='method', args._actions))[0])
        args._actions[method_arg_idx].choices = ['web']
        args._actions[method_arg_idx].default = 'web'
        args._actions[method_arg_idx].help = 'docking method. Currently support "web". Default is %(default)s.'
        return args
    
    @staticmethod
    @hdock_run_fn_warpper('HPEPDOCK')
    def run_hdock_web(config_path: Path, parameters: Dict[str, str] = None, email=None, browser=None,
                      args: argparse.ArgumentParser = None):
        w_dir = config_path.parent if isinstance(config_path, Path) else Path(config_path[0]).resolve().parent
        run_dock_on_HPEPDOCK(receptor_path=parameters['receptor_path'], ligand_path=parameters['ligand_path'],
                             w_dir=w_dir, email=email, browser=browser)


class dinc_ensemble(hdock):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        
    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args = hdock.make_args(args)
        args.add_argument('--use-config-box', action='store_true', default=False,
                          help='FLAG, whether to use grid box from config file. Default is %(default)s.')
        # make sure only support web method
        method_arg_idx = args._actions.index(list(filter(lambda x: x.dest =='method', args._actions))[0])
        args._actions[method_arg_idx].choices = ['web']
        args._actions[method_arg_idx].default = 'web'
        args._actions[method_arg_idx].help = 'docking method. Currently support "web". Default is %(default)s.'
        # make sure eamil is required
        email_arg_idx = args._actions.index(list(filter(lambda x: x.dest == 'email', args._actions))[0])
        args._actions[email_arg_idx].required = True
        args._actions[email_arg_idx].help = 'email address for HDOCK web server. Required.'
        return args
    
    def process_args(self):
        if self.args.receptor is not None and self.args.ligand is not None and self.args.config_name is not None:
            return put_err('receptor, ligand and config_name cannot be provided at the same time, exit.', _exit=True)
        return super().process_args()
    
    @staticmethod
    @hdock_run_fn_warpper('', 'DINC-Ensemble_result.zip')
    def run_hdock_web(config_path: Path, parameters: Dict[str, str] = None, email=None, browser=None,
                      args: argparse.ArgumentParser = None):
        w_dir = config_path.parent if isinstance(config_path, Path) else Path(config_path[0]).resolve().parent
        put_log(f'working in {w_dir}')
        if args.use_config_box: # type: ignore
            box_center = {k:parameters[f'center_{k}'] for k in ['x', 'y', 'z']}
            box_size = {k:parameters[f'size_{k}'] for k in ['x', 'y', 'z']}
        else:
            box_center, box_size = 'receptor', 'ligand'
        put_log(f'parameters: {parameters}, box_center: {box_center}, box_size: {box_size}')
        run_dock_on_DINC_ensemble(receptor_path=parameters['receptor_path'], ligand_path=parameters['ligand_path'],
                                  email=email, box_center=box_center, box_size=box_size, w_dir=str(w_dir), browser=browser)


def convert_result_run_convert(input_path: Path, output_path: Path, method: str):
    if method in {'lazydock', 'obabel'}:
        getattr(convert_result, f'run_convert_{method}')(input_path, output_path)
    else:
        put_err(f'unsupported convert method: {method}, exit.', _exit=True)


class convert_result(vina):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.taskpool = None

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--name', type = str, default='',
                          help='input file name. Default is %(default)s.')
        args.add_argument('-i', '--input-type', type = str, default='pdbqt,dlg',
                          help='input file type. Default is %(default)s.')
        args.add_argument('-o', '--output-type', type = str, default='pdb', choices=['pdb'],
                          help='input file type. Default is %(default)s.')
        args.add_argument('-s', '--suffix', type = str, default='',
                          help='output file suffix. Default is %(default)s.')
        args.add_argument('-m', '--method', type = str, default='lazydock', choices=['lazydock', 'obabel'],
                          help='convert tools to use. Currently support "lazydock, obabel". Default is %(default)s.')
        args.add_argument('--n-workers', type=int, default=4,
                          help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = list(map(clean_path, self.args.batch_dir))
        self.args.input_type = self.args.input_type.split(',')
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.taskpool = TaskPool('process', self.args.n_workers).start()
        
    @staticmethod
    def run_convert_lazydock(input_path: Path, output_path: Path):
        dlg = DlgFile(input_path)
        pdb_string = '\n'.join([f'MODEL {i+1}\n'+pose.as_pdb_string().replace('\r\n', '\n')+'ENDMDL' for i, pose in enumerate(dlg.pose_lst)])
        opts_file(output_path, 'w', data=pdb_string)

    @staticmethod
    def run_convert_obabel(input_path: Path, output_path: Path):
        ty1, ty2 = input_path.suffix.lower()[1:], output_path.suffix.lower()[1:]
        os.system(f'obabel -i{ty1} "{str(input_path)}" -o{ty2} -O "{str(output_path)}"')

    def main_process(self):
        input_paths = get_paths_with_extension(self.args.batch_dir, self.args.input_type,
                                               name_substr=self.args.name, sort='natsort')
        print(f'get {len(input_paths)} input(s) for convert:\n', '\n'.join([f'{i+1}. {x}' for i, x in enumerate(input_paths)]))
        if input('start convert? (y/n) ').lower() != 'y':
            return 
        for input_path in tqdm(input_paths, total=len(input_paths)):
            input_path = Path(input_path)
            output_path = input_path.parent / f'{input_path.stem}{self.args.suffix}.{self.args.output_type}'
            self.taskpool.add_task(None, convert_result_run_convert, input_path, output_path, self.args.method)
            self.taskpool.wait_till_free()
        self.taskpool.wait_till_all_done()
        self.taskpool.close(1)


class cluster_result(vina):
    def __init__(self, args, printf = print):
        super().__init__(args, printf)
        self.taskpool = None

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '-bd', '--batch-dir', type = str, nargs='+', default=['.'],
                          help="dir which contains many sub-folders, each sub-folder contains input files, default is %(default)s.")
        args.add_argument('-n', '--name', type = str, required=True,
                          help='input file name, such as "dock.pdbqt".')
        args.add_argument('-m', '--method', type = str, default='pose', choices=['pose', 'interaction'],
                          help='cluster method to use. Currently support "pose, interaction". Default is %(default)s.')
        args.add_argument('--range', type = str, default='2,7',
                          help='range of cluster groups, input as "min,max". Default is %(default)s.')
        args.add_argument('--repeat', type=int, default=3,
                          help='number of repeats for clustering. Default is %(default)s.')
        args.add_argument('--n-workers', type=int, default=4,
                          help='number of tasks to parallel docking. Default is %(default)s.')
        return args
    
    def process_args(self):
        self.args.batch_dir = self.process_batch_dir_lst(self.args.batch_dir)
        self.args.range = [int(x) for x in self.args.range.split(',')]
        if len(self.args.range) != 2 or self.args.range[0] >= self.args.range[1]:
            put_err(f'range should be "min,max", and min < max, got {self.args.range}, exit.', _exit=True)
        if self.args.n_workers <= 0:
            put_err(f'n_workers must be positive integer, got {self.args.n_workers}, exit.', _exit=True)
        self.taskpool = TaskPool('process', self.args.n_workers).start()
        
    @staticmethod
    def cluster_on_pose(input_path: Path, _range: Tuple[int, int], repeat: int):
        wdir = input_path.parent / 'cluster_pose'
        os.makedirs(str(wdir), exist_ok=True)
        dlg = DlgFile(input_path)
        dlg.sort_pose(inplace=True)
        rmsd = dlg.rmsd('numpy')
        sns.clustermap(rmsd)
        save_show(str(wdir / 'RMSD_cluster.png'), 600, show=False)
        rs_df = pd.DataFrame(columns=['repeat', 'k', 'sse','ssr', 'r'])
        dfs = {'rs_df': rs_df}
        for i in range(repeat):
            for k in range(_range[0], _range[1]+1):
                groups_idx = dlg.rmsd_cluster(k) # [N, ]
                groups_idx_lst = list(groups_idx)
                group_sizes = {idx: groups_idx_lst.count(idx) for idx in range(k)}
                sse, ssr, r = dlg.calcu_SSE_SSR(rmsd, groups_idx)
                rs_df.loc[len(rs_df)] = [i, k, sse, ssr, r]
                dfs[f'repeat_{i}_k_{k}'] = pd.DataFrame(data={'pose_idx': range(len(dlg.pose_lst)),
                                                              'groups_idx': groups_idx,
                                                              'group_size': [group_sizes[groups_idx[idx]] for idx in range(len(dlg))],
                                                              'energy': list(map(lambda x: x.energy, dlg.pose_lst))})
        fig = plt.figure()
        sns.lineplot(data=rs_df, x='k', y='r', hue='repeat', ax=fig.gca())
        save_show(str(wdir / 'SSR-SSE_lines.png'), 600, show=False)
        plt.close(fig)
        write_sheets(str(wdir / 'cluster.xlsx'), dfs)
        put_log(f'cluster result saved to {wdir}.')

    @staticmethod
    def cluster_on_interaction(input_path: Path):
        raise NotImplementedError('interaction clustering not implemented yet.')

    def main_process(self):
        input_paths = get_paths_with_extension(self.args.batch_dir, [],
                                               name_substr=self.args.name, sort='natsort')
        for input_path in tqdm(input_paths, total=len(input_paths)):
            input_path = Path(input_path)
            self.taskpool.add_task(None, getattr(self, f'cluster_on_{self.args.method}'),
                                                input_path, self.args.range, self.args.repeat)
            self.taskpool.wait_till_free()
        self.taskpool.wait_till_all_done()
        self.taskpool.close(1)


_str2func = {
    'vina': vina,
    'redock': redock,
    'redock-gpu': redock_gpu,
    'redock-flex-gpu': redock_flex_gpu,
    'redock-flex-vina': redock_flex_vina,
    'hdock': hdock,
    'hpepdock': hpepdock,
    'dinc-ensemble': dinc_ensemble,
    'convert-result': convert_result,
    'cluster-result': cluster_result,
}

def main(sys_args: List[str] = None):
    args_paser = argparse.ArgumentParser(description = 'perform docking with Vina or other docking software.')
    subparsers = args_paser.add_subparsers(title='subcommands', dest='sub_command')
    
    for k, v in _str2func.items():
        v.make_args(subparsers.add_parser(k))

    excute_command(args_paser, sys_args, _str2func)


if __name__ == "__main__":
    # main('convert-result -d data_tmp/docking/ligand1 -m lazydock --n-workers 1'.split())
    # main('dinc-ensemble -d data_tmp/docking/ligand1 -m web --email 2262029386@qq.com --config-name config.txt --use-config-box'.split())
    # main('cluster-result -d data_tmp/docking/ligand1 -n dock.pdbqt'.split())
    
    main()