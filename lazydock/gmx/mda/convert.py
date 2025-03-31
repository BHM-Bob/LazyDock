'''
Date: 2025-02-05 14:26:31
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-03-31 19:27:24
Description: 
'''
import warnings
from typing import Dict

import MDAnalysis
import numpy as np
from mbapy_lite.base import put_err
from mbapy_lite.web_utils.task import TaskPool
from MDAnalysis import AtomGroup, Universe
from MDAnalysis.coordinates.base import IOBase
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib import util


class FakeIOWriter:
    def __init__(self):
        self.str_lst = []

    def write(self, content: str):
        self.str_lst.append(content)


class FakeAtomGroup(PDBWriter):
    def __init__(self, ag: MDAnalysis.AtomGroup, reindex: bool = False):
        # 保留原有属性提取
        self.positions = self.convert_pos_to_native(ag.positions.copy(), inplace=False)
        
        # 预提取所有必要属性（固定长度字符串类型）
        self.resnames = self._get_attr(ag, 'resnames', 'UNK', 'U4')
        self.chainIDs = self._process_chainIDs(ag)
        self.resids = self._get_attr(ag, 'resids', 1, int)
        self.names = self._get_attr(ag, 'names', 'X', 'U4')
        self.elements = self._get_attr(ag, 'elements', '', 'U2')
        self.altLocs = self._get_attr(ag, 'altLocs', ' ', 'U1')
        self.occupancies = self._get_attr(ag, 'occupancies', 1.0, float)
        self.tempfactors = self._get_attr(ag, 'tempfactors', 0.0, float)
        self.segids = self._get_attr(ag, 'segids', ' ', 'U4')
        self.charges = self._get_attr(ag, 'formalcharges', 0, int)
        self.icodes = self._get_attr(ag, 'icodes', ' ', 'U1')
        
        # 处理序列号（向量化操作）
        self.ids = self._process_ids(ag, reindex)

    def _get_attr(self, ag, attr, default, dtype):
        """统一属性获取方法"""
        if hasattr(ag, attr):
            arr = getattr(ag, attr)
            return arr.astype(dtype) if dtype else arr.copy()
        return np.full(len(ag), default, dtype=dtype)

    def _process_chainIDs(self, ag):
        """完整保留原始chainID以便替换操作"""
        if hasattr(ag, 'chainIDs'):
            # 存储原始chainID（最多4字符）
            return ag.chainIDs.astype('U4')
        return np.full(len(ag), 'X', dtype='U4')

    def _process_ids(self, ag, reindex):
        """处理原子序列号"""
        if reindex or not hasattr(ag, 'ids'):
            return np.arange(1, len(ag)+1, dtype=np.int32)
        return ag.ids.astype(np.int32)


class PDBConverter(PDBWriter):
    def __init__(self, ag: MDAnalysis.AtomGroup, reindex: bool = False):
        # 增加类型检查和多进程支持
        if not isinstance(ag, FakeAtomGroup):
            ag = FakeAtomGroup(ag, reindex)
        self.fake_ag = ag
        
        # 保持父类初始化参数
        self.convert_units = False
        self._reindex = reindex
        self.pdbfile = FakeIOWriter()

        # 初始化计数器
        self.frames_written = 0

    def check_charges(self, fag: FakeAtomGroup):
        """检查电荷值是否合法并格式化为字符串"""
        charges = np.full(len(fag.charges), '', dtype='U3')
        pos_mask = fag.charges > 0
        neg_mask = fag.charges < 0
        # 异常检查（向量化）
        if np.any(fag.charges > 9) or np.any(fag.charges < -9):
            invalid_charges = np.unique(fag.charges[(fag.charges > 9) | (fag.charges < -9)])
            raise ValueError(f"Invalid formal charges: {invalid_charges}")
        # 向量化字符串生成
        charges[pos_mask] = np.char.add(
            fag.charges[pos_mask].astype('U1'), 
            np.full(pos_mask.sum(), '+', dtype='U1')
        )
        charges[neg_mask] = np.char.add(
            (-fag.charges[neg_mask]).astype('U1'), 
            np.full(neg_mask.sum(), '-', dtype='U1')
        )
        return charges

    def _write_single_timestep_fast(self, alter_chain=None, alter_res=None, alter_atm=None):
        # 批量替换函数
        def batch_replace(arr, mapping):
            if mapping:
                for k, v in mapping.items():
                    arr[arr == k] = v
            return arr

        # 获取预处理数据
        fag = self.fake_ag
        resnames = batch_replace(fag.resnames.copy(), alter_res or {})
        chainIDs = batch_replace(fag.chainIDs.copy(), alter_chain or {})
        # 在写入前截断为1个字符（符合PDB规范）
        chainIDs = np.array([c[:1] for c in chainIDs], dtype='U1')
        
        # 生成记录类型
        record_types = np.full_like(chainIDs, 'ATOM', dtype='U6')
        if alter_atm:
            for chain, rec_type in alter_atm.items():
                record_types[chainIDs == chain] = rec_type

        # 预生成格式化数据
        serials = np.vectorize(util.ltruncate_int)(fag.ids, 5)
        resSeqs = np.vectorize(util.ltruncate_int)(fag.resids, 4)
        elements = np.char.upper(fag.elements)
        charges = self.check_charges(fag)
        
        # 生成所有原子行
        for i in range(len(fag.ids)):
            vals = {
                'serial': serials[i],
                'name': self._deduce_PDB_atom_name(fag.names[i], resnames[i]),
                'altLoc': fag.altLocs[i],
                'resName': resnames[i][:4],
                'chainID': chainIDs[i],
                'resSeq': resSeqs[i],
                'iCode': fag.icodes[i],
                'pos': fag.positions[i],
                'occupancy': fag.occupancies[i],
                'tempFactor': fag.tempfactors[i],
                'segID': fag.segids[i][:4],
                'element': elements[i],
                'charge': charges[i]
            }
            
            try:
                self.pdbfile.write(self.fmt[record_types[i]].format(**vals))
            except KeyError:
                raise ValueError(f"Invalid record type: {record_types[i]}")

        self.frames_written += 1
    
    def fast_convert(self, alter_chain: Dict[str,str] = None,
                     alter_res: Dict[str,str] = None, alter_atm: Dict[str,str] = None):
        """
        Convert the AtomGroup to a PDB string.
        
        Parameters
            - alter_chain: Dict[str,str]: key is orignal chain name, value is target chain name
            - alter_res: Dict[str,str]: key is orignal res name, value is target res name
            - alter_atm: Dict[str,str]: key is target chain name, value is target atom record type
            
        Returns
        str
            The PDB string.
        """
        # self.ts = self.obj.universe.trajectory.ts
        # self.frames_written = 1
        self._write_single_timestep_fast(alter_chain, alter_res, alter_atm)
        return ''.join(self.pdbfile.str_lst)
    
    def close(self):
        pass
