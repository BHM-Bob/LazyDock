'''
Date: 2024-09-01 20:33:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-09-01 22:09:38
Description: 
'''
import pickle
import socket
import time
from threading import Lock, Thread
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from mbapy.base import get_fmt_time
from pymol import api as _api_
from pymol import cmd as _cmd_


class PymolAPI(object):
    """client module for pymol.api and pymol.cmd"""
    def __init__(self, client: 'VClient', api: str = 'cmd') -> None:
        self._client_ = client
        self._api_ = api
    def __getattribute__(self, name: str) -> Any:
        if name in {'_client_', '_api_'}:
            return super().__getattribute__(name)
        return lambda *args, **kwargs: self._client_.send_action(self._api_, name, *args, **kwargs)
        
        
class VServer:
    """run by pymol"""
    def __init__(self, ip: str = 'localhost', port: int = 8085, verbose: bool = False) -> None:
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.bind((ip, port))
        self.server = Thread(target=self._run_server, daemon=True)
        self.server.start()
        self.lock = Lock()
        self._logs: List[str] = []
        self._copied_log_len = 0
        self._verbose = verbose
        
    def copy_logs(self):
        with self.lock:
            return self._logs.copy()
        
    def copy_new_logs(self):
        with self.lock:
            new_logs = self._logs[self._copied_log_len:].copy()
            self._copied_log_len = len(self._logs)
            return new_logs
        
    def _add_log(self, log: str) -> None:
        with self.lock:
            self._logs.append(f'[{get_fmt_time()}] {log}')
            if self._verbose:
                print(f'[{get_fmt_time()}] {log}')
        
    def _run_server(self) -> None:
        """main loop run in a thread"""
        self.socket.listen(1)
        conn, addr = self.socket.accept()
        while True:
            # get and un-pickle data from client
            recv_data = conn.recv(4096)
            if not recv_data:
                continue
            api, fn, args, kwargs = pickle.loads(recv_data)
            # execute action
            if api not in {'cmd', 'api'}:
                self._add_log(f'Error in executing {api}.{fn}: api not found')
            try:
                if api == 'cmd':
                    getattr(_cmd_, fn)(*args, **kwargs)
                elif api == 'api':
                    getattr(_api_, fn)(*args, **kwargs)
                self._add_log(f'{api}.{fn} executed successfully')
            except Exception as e:
                self._add_log(f'Error in executing {api}.{fn}: {e}')
            # sleep to avoid busy loop
            time.sleep(0.05)
        
        
class VClient:
    """run by user in anthor python script"""
    def __init__(self, ip: str = 'localhost', port: int = 8085) -> None:
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect((ip, port))
    
    def send_action(self, api: str, fn: str, *args, **kwargs) -> Any:
        self.socket.sendall(pickle.dumps((api, fn, args, kwargs)))
        

def _test_server():
    _cmd_.reinitialize()
    server = VServer()
    while True:
        time.sleep(0.1)
        print(server.copy_new_logs())
        

def _test_client():
    client = VClient()
    vcmd = PymolAPI(client, 'cmd')
    vcmd.load('data_tmp/pdb/LIGAND.pdb', 'ligand')
    

if __name__ == '__main__':
    from multiprocessing import Process
    p1 = Process(target=_test_server)
    p2 = Process(target=_test_client)
    p1.start()
    time.sleep(3)
    p2.start()
    p1.join()
    p2.join()
    