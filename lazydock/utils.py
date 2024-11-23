'''
Date: 2024-08-18 13:26:20
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-23 19:55:05
Description: 
'''
import os
import uuid


def uuid4():
    """return a random 4-character string"""
    return uuid.uuid4().hex[:4]

def get_storage_path(sub_path:str):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'storage', sub_path)