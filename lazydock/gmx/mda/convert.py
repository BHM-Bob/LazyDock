'''
Date: 2025-02-05 14:26:31
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-06 15:51:37
Description: 
'''
import MDAnalysis
from MDAnalysis.coordinates.PDB import PDBWriter
from mbapy_lite.base import put_err


class FakeIOWriter:
    def __init__(self):
        self.string = ''

    def write(self, content: str):
        self.string += content


class PDBConverter(PDBWriter):
    """
    Convert an MDAnalysis AtomGroup to a PDB string.
    """
    def __init__(self, ag: MDAnalysis.AtomGroup, reindex: bool = False):
        """
        Parameters
        ----------
        ag : MDAnalysis.core.groups.AtomGroup
            The AtomGroup to convert.
        reindex : bool, optional
            Whether to reindex the AtomGroup, by default False
        """
        self.obj = ag
        self.convert_units = True
        self._multiframe = self.multiframe
        self.bonds = "conect"
        self._reindex = reindex
        
        self.start = self.frames_written = 0
        self.step = 1
        self.remarks = '"Created by MDAnalysis.coordinates.PDB.PDBWriter"'
        
        self.pdbfile = FakeIOWriter()
        self.has_END = False
        self.first_frame_done = False
        
    def convert(self):
        """
        Convert the AtomGroup to a PDB string.
        Returns
        -------
        str
            The PDB string.
        """
        self._update_frame(self.obj)
        self._write_pdb_header()
        try:
            ts = self.ts
        except AttributeError:
            return put_err("no coordinate data to write to trajectory file, return None")
        self._check_pdb_coordinates()
        self._write_timestep(ts)
        return self.pdbfile.string
    
    def fast_convert(self):
        """
        Convert the AtomGroup to a PDB string.
        Returns
        -------
        str
            The PDB string.
        """
        self.ts = self.obj.universe.trajectory.ts
        self.frames_written = 1
        self._write_timestep(self.ts, multiframe=False)
        return self.pdbfile.string
