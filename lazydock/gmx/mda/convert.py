'''
Date: 2025-02-05 14:26:31
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-19 16:32:45
Description: 
'''
from typing import Dict
import warnings
import MDAnalysis
from MDAnalysis.coordinates.PDB import PDBWriter
from MDAnalysis.lib import util
from MDAnalysis.exceptions import NoDataError
from mbapy_lite.base import put_err
import numpy as np


class FakeIOWriter:
    def __init__(self):
        self.str_lst = []

    def write(self, content: str):
        self.str_lst.append(content)


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
        return ''.join(self.pdbfile.str_lst)
    
    def _write_single_timestep_fast(self, alter_chain: Dict[str,str] = None,
                                    alter_res: Dict[str,str] = None, alter_atm: Dict[str,str] = None):
        alter_chain = alter_chain or {}
        alter_res = alter_res or {}
        alter_atm = alter_atm or {}
        atoms = self.obj.atoms
        pos = atoms.positions
        if self.convert_units:
            pos = self.convert_pos_to_native(pos, inplace=False)

        # Make zero assumptions on what information the AtomGroup has!
        # theoretically we could get passed only indices!
        def get_attr(attrname, default):
            """Try and pull info off atoms, else fake it

            attrname - the field to pull of AtomGroup (plural!)
            default - default value in case attrname not found
            """
            try:
                return getattr(atoms, attrname)
            except AttributeError:
                if self.frames_written == 0:
                    warnings.warn("Found no information for attr: '{}'"
                                  " Using default value of '{}'"
                                  "".format(attrname, default))
                return np.array([default] * len(atoms))
        altlocs = get_attr('altLocs', ' ')
        resnames = get_attr('resnames', 'UNK')
        icodes = get_attr('icodes', ' ')
        segids = get_attr('segids', ' ')
        chainids = get_attr('chainIDs', '')
        resids = get_attr('resids', 1)
        occupancies = get_attr('occupancies', 1.0)
        tempfactors = get_attr('tempfactors', 0.0)
        atomnames = get_attr('names', 'X')
        elements = get_attr('elements', ' ')
        record_types = list(get_attr('record_types', 'ATOM'))
        for alter_attr, alter_dict in zip([resnames, chainids], [alter_res, alter_chain]):
            for k, v in alter_dict.items():
                alter_attr[alter_attr == k] = v
        formal_charges = self._format_PDB_charges(get_attr('formalcharges', 0))

        def validate_chainids(chainids, default):
            """Validate each atom's chainID

            chainids - np array of chainIDs
            default - default value in case chainID is considered invalid
            """
            invalid_length_ids = False
            invalid_char_ids = False
            missing_ids = False

            for (i, chainid) in enumerate(chainids):
                if chainid == "":
                    missing_ids = True
                    chainids[i] = default
                elif len(chainid) > 1:
                    invalid_length_ids = True
                    chainids[i] = default
                elif not chainid.isalnum():
                    invalid_char_ids = True
                    chainids[i] = default

            if invalid_length_ids:
                warnings.warn("Found chainIDs with invalid length."
                              " Corresponding atoms will use value of '{}'"
                              "".format(default))
            if invalid_char_ids:
                warnings.warn("Found chainIDs using unnaccepted character."
                              " Corresponding atoms will use value of '{}'"
                              "".format(default))
            if missing_ids:
                warnings.warn("Found missing chainIDs."
                              " Corresponding atoms will use value of '{}'"
                              "".format(default))
            return chainids

        chainids = validate_chainids(chainids, "X")

        # If reindex == False, we use the atom ids for the serial. We do not
        # want to use a fallback here.
        if not self._reindex:
            try:
                atom_ids = atoms.ids
            except AttributeError:
                raise NoDataError(
                    'The "id" topology attribute is not set. '
                    'Either set the attribute or use reindex=True.'
                )
        else:
            atom_ids = np.arange(len(atoms)) + 1

        for i in range(len(atoms)):
            vals = {}
            vals['serial'] = util.ltruncate_int(atom_ids[i], 5)  # check for overflow here?
            vals['name'] = self._deduce_PDB_atom_name(atomnames[i], resnames[i])
            vals['altLoc'] = altlocs[i][:1]
            vals['resName'] = resnames[i][:4]
            vals['resSeq'] = util.ltruncate_int(resids[i], 4)
            vals['iCode'] = icodes[i][:1]
            vals['pos'] = pos[i]  # don't take off atom so conversion works
            vals['occupancy'] = occupancies[i]
            vals['tempFactor'] = tempfactors[i]
            vals['segID'] = segids[i][:4]
            vals['chainID'] = chainids[i]
            vals['element'] = elements[i][:2].upper()
            vals['charge'] = formal_charges[i]

            # record_type attribute, if exists, can be ATOM or HETATM
            if chainids[i] in alter_atm:
                record_types[i] = alter_atm[chainids[i]]
            try:
                self.pdbfile.write(self.fmt[record_types[i]].format(**vals))
            except KeyError:
                errmsg = (f"Found {record_types[i]} for the record type, but "
                          f"only allowed types are ATOM or HETATM")
                raise ValueError(errmsg) from None

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
        self.ts = self.obj.universe.trajectory.ts
        self.frames_written = 1
        self._write_single_timestep_fast(alter_chain, alter_res, alter_atm)
        return ''.join(self.pdbfile.str_lst)
