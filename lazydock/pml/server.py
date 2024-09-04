'''
Date: 2024-09-01 20:33:00
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-09-04 21:26:46
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
from pymol import api as pml_api
from pymol import cmd as pml_cmd


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
        self.ip = ip
        self.port = int(port)
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 1 * 1024 * 1024) # set socket recv buffer size to 1MB
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
        
    def _add_log(self, log: str, verbose: bool = False) -> None:
        with self.lock:
            self._logs.append(f'[{get_fmt_time()}] {log}')
            if self._verbose or verbose:
                print(f'[{get_fmt_time()}] {log}')

    def _recvall(self, sock: socket.socket, count: int):
        buf = b''
        while count:
            newbuf = sock.recv(count)
            if not newbuf: return None
            buf += newbuf
            count -= len(newbuf)
        return buf
    
    def _sendall(self, scok: socket.socket, data: bytes) -> None:
        length = len(data)
        scok.sendall(length.to_bytes(4, 'big'))
        scok.sendall(data)
        
    def _run_server(self) -> None:
        """main loop run in a thread"""
        self.socket.bind((self.ip, self.port))
        self.socket.listen(1)
        while True:
            conn, addr = self.socket.accept()
            self._add_log(f'Connected by {addr}', verbose=True)
            while True:
                # get and un-pickle data from client
                recv_data = conn.recv(4)
                if not recv_data:
                    continue
                data_length = int.from_bytes(recv_data, 'big')
                recv_data = self._recvall(conn, data_length)
                if not recv_data:
                    self._add_log(f'VServer: Failed to receive data from client for length {data_length}')
                    continue
                if recv_data == b'quit':
                    self._add_log('VServer: Client quit')
                    break
                api, fn, args, kwargs = pickle.loads(recv_data)
                # execute action
                if api not in {'cmd', 'api'}:
                    self._add_log(f'VServer: Error in executing {api}.{fn}: api not found')
                try:
                    if api == 'cmd':
                        ret = getattr(pml_cmd, fn)(*args, **kwargs)
                    elif api == 'api':
                        ret = getattr(pml_api, fn)(*args, **kwargs)
                    self._sendall(conn, pickle.dumps(ret))
                    self._add_log(f'VServer: {api}.{fn} executed successfully')
                except Exception as e:
                    self._add_log(f'VServer: Error in executing {api}.{fn}: {e}')
                # sleep to avoid busy loop
                time.sleep(0.05)
            # close connection
            conn.close()
            self._add_log(f'Connect closed by {addr}', verbose=True)
        
        
class VClient(VServer):
    """run by user in anthor python script"""
    def __init__(self, ip: str = 'localhost', port: int = 8085) -> None:
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_SNDBUF, 1 * 1024 * 1024) # set socket send buffer size to 1MB
        self.socket.connect((ip, port))
        
    def _run_server(self) -> None:
        raise NotImplementedError('VClient cannot run as a server now')
    
    def send_action(self, api: str, fn: str, *args, **kwargs) -> Any:
        data = pickle.dumps((api, fn, args, kwargs))
        self._sendall(self.socket, data)
        try:
            ret_len = int.from_bytes(self.socket.recv(4), 'big')
            ret = self._recvall(self.socket, ret_len)
            return pickle.loads(ret)
        except Exception as e:
            print(f'Error in receiving data: {e}')
            return None
        
    def sen_quit(self) -> None:
        data = b'quit'
        self._sendall(self.socket, data)
        self.socket.close()
        
        
__all__ = [
    'PymolAPI',
    'VServer',
    'VClient',
]
        

def _test_server():
    pml_cmd.reinitialize()
    server = VServer()
    while True:
        time.sleep(0.1)
        print(server.copy_new_logs())
        

def _test_client():
    client = VClient()
    vcmd = PymolAPI(client, 'cmd')
    print(vcmd.load('data_tmp/pdb/LIGAND.pdb', 'ligand'))
    

if __name__ == '__main__':
    from multiprocessing import Process
    p1 = Process(target=_test_server)
    p2 = Process(target=_test_client)
    p1.start()
    time.sleep(3)
    p2.start()
    p1.join()
    p2.join()
    