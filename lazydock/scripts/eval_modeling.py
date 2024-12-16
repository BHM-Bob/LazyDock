'''
Date: 2024-12-16 15:32:10
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-12-16 17:51:35
Description: 
'''
import argparse
import os
import time
from pathlib import Path
from typing import Dict, List

import pandas as pd
from tqdm import tqdm
from mbapy_lite.base import get_fmt_time, put_err
from mbapy_lite.file import (get_paths_with_extension, is_jsonable, opts_file,
                             write_sheets)
from mbapy_lite.web import TaskPool, random_sleep
from pymol import cmd

from lazydock.scripts._script_utils_ import clean_path, show_args
from lazydock.web.model_eval import get_eval_info_from_servers


def get_info(result_path, *args, **kwargs):
    result = get_eval_info_from_servers(*args, **kwargs)
    opts_file(result_path, 'wb', way='pkl', data=result)
    return result


def asign_chain_code(pdb_path: str, save_path: str, name: str, code: str):
    cmd.load(pdb_path, name)
    cmd.alter(f'{name}', f'chain="{code}"')
    cmd.save(save_path, name)
    cmd.delete(name)


def main(sys_args: List[str] = None):
    SUPPORT_SERVERS = ['ProQ', 'VoroMQA', 'ProSA', 'MolProbity', 'ProQ3', 'SAVES', 'QMEANDisCo', 'QMEAN']
    args_paser = argparse.ArgumentParser(description = 'Run model evaluation on web servers.')
    args_paser.add_argument('-d', '--dir', type = str, default='.',
                            help='input directory. Default is %(default)s.')
    args_paser.add_argument('-n', '--name', type = str, default='',
                            help='sub string of file name. Default is %(default)s.')
    args_paser.add_argument('-t', '--type', type = str, default='pdb',
                            help='file type to search in input directory. Default is %(default)s.')
    args_paser.add_argument('--servers', type=str, default=','.join(SUPPORT_SERVERS),
                            help=f'web servers to evaluate, split by comma. Default is {SUPPORT_SERVERS}.')
    args_paser.add_argument('--n-workers', type=int, default=4,
                            help='number of workers to parallel. Default is %(default)s.')
    args_paser.add_argument('-r', '--recursive', action='store_true', default=False,
                            help='FLAG, recursive search. Default is %(default)s.')
    args_paser.add_argument('--disable-cache', action='store_true', default=False,
                            help='FLAG, whether to disable cache. Default is %(default)s.')
    args = args_paser.parse_args(sys_args)
    # process IO path
    args.dir = clean_path(args.dir)
    if os.path.isdir(args.dir):
        paths = get_paths_with_extension(args.dir, [args.type], name_substr=args.name, recursive=args.recursive)
        paths = list(map(lambda p: Path(p).resolve(), paths))
    else:
        return put_err(f'Input dir {args.dir} is not a dir, exit.', _exit=True)
    # process servers
    args.servers = args.servers.split(',')
    if not all(s in SUPPORT_SERVERS for s in args.servers):
        unsupports = list(set(args.servers) - set(SUPPORT_SERVERS))
        return put_err(f'Unsupported servers: {unsupports}.', _exit=True)
    # show args
    show_args(args, ['dir', 'name', 'n_workers','recursive', 'disable_cache'])
    # run tasks
    cmd.reinitialize()
    tasks = []
    taskpool = TaskPool('threads', args.n_workers).start()
    for path in tqdm(paths, total=len(paths), desc='Running tasks'):
        chain_alter_path = path.parent / f'{path.stem}_chain_alter{path.suffix}'
        result_path = os.path.join(args.dir, f'{chain_alter_path}.pkl')
        if not args.disable_cache and os.path.exists(result_path):
            print(f'Skip {path}, result file exists: {result_path}.')
            tasks.append(result_path)
            continue
        while taskpool.count_waiting_tasks() > 0:
            time.sleep(5)
        asign_chain_code(path, chain_alter_path, path.stem, 'A')
        task_name = taskpool.add_task(result_path, get_info, # task name is the result path
                                      result_path, chain_alter_path,
                                      servers=args.servers)
        tasks.append(task_name)
    taskpool.wait_till_tasks_done(tasks)
    # gather results
    df = pd.DataFrame(columns=['path'])
    df.set_index(['path'], inplace=True)
    for result_path, path in zip(tasks, paths):
        data = opts_file(result_path, 'rb', way='pkl')
        for server, result in data.items():
            for name, value in result.items():
                if is_jsonable(value):
                    df.loc[path, f'{server}:{name}'] = value
    df.sort_index(inplace=True)
    df.to_excel(os.path.join(args.dir, f'eval_{get_fmt_time()}.xlsx'))


if __name__ == "__main__":
    main()