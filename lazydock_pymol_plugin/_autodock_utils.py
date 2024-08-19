# copied and fixed from Autodock Plugin: http://www.pymolwiki.org/index.php/autodock_plugin


from __future__ import absolute_import, print_function

import os
import re
import sys
import tkinter.filedialog as tkFileDialog
from tkinter import *


# atom-type, atom-number, atom-name, residue-name, chain-name, residue-number, x, y, z, occupancy, temperature-factor
# ATOM      1  CA  LYS     7     136.747 133.408 135.880 -0.06 +0.10
PDB_PATTERN = r"(ATOM|HETATM) +(\d+) +(\w+) +(\w+) +(\w+)? +(\d+) +([\d\-\.]+) +([\d\-\.]+) +([\d\-\.]+) +([\+\-][\d\-\.]+) +([\+\-][\d\-\.]+)"


class ADModel:
    """ STORAGE CLASS FOR DOCKED LIGANDS """
    def __init__(self, lst=None, sort_atom_by_res: bool = False):
        self.atomlines = []
        self.energy = 0.
        self.name = ''
        self.poseN = 0
        self.info = ''
        self.lst = []
        self.num = 0
        self.as_string = ''
        self.pdb_lines = []
        if lst is not None:
            self.from_list(lst, sort_atom_by_res)
            
    def sort_atom_by_res(self, pdb_string: str = None):
        pdb_string = pdb_string or self.as_string
        pdb_string_lines = pdb_string.strip().split('\n')
        for idx, line in enumerate(pdb_string_lines):
            matches = re.findall(PDB_PATTERN, line)
            items = list(matches[0])
            items[1], items[5] = int(items[1]), int(items[5])
            self.pdb_lines.append(items+[idx]) # line items + line idx
        # sort by chain name first, then by residue number, then by atom number
        self.pdb_lines.sort(key=lambda x: (x[4], x[5], x[1]))
        # apply sorted line idx to as_string
        self.as_string = '\n'.join([pdb_string_lines[i[-1]] for i in self.pdb_lines])
        return self.as_string
            

    def from_list(self, lst, sort_atom_by_res: bool = False):
        self.info = '\n'.join(lst)
        self.lst = lst
        for line in lst:
            if line.startswith('ATOM') or \
                    line.startswith('HETATM'):
                self.atomlines.append(line)
                self.as_string += 'ATOM  ' + line[6:67] + '\n'
            elif line.startswith('USER'):
                if 'Free Energy of Binding' in line:
                    entr = line.split('=')[1]
                    self.energy = float(entr.split()[0])
            elif line.startswith('REMARK'):
                if 'VINA RESULT' in line:
                    entr = line.split(':')[1]
                    self.energy = float(entr.split()[0])
        if sort_atom_by_res:
            self.sort_atom_by_res()

    def as_pdb_string(self):
        return self.as_string
        

    def info_string(self):
        return self.info


#---------------------------------------------------------------------------

class ADGridMap:

    """ CLASS FOR HANDLING AUTODOCK GRID MAP FILES"""

    def __init__(self, fp=None, name='map'):

        self.name = ''
        self.npts = [0, 0, 0]
        self.n = [0, 0, 0]
        self.center = [0, 0, 0]
        self.origin = [0, 0, 0]
        self.nelem = 0
        self.spacing = 0.
        self.values = []
        self.datafile = ''
        self.molecule = ''
        self.paramfile = ''
        self.precision = 0.0001
        if fp is not None:
            self.read(fp, name)

    def read(self, fp, name='map'):
        self.name = name
        for i in range(6):
            line = fp.readline()
            if i == 0:
                self.paramfile = line.split()[1]
            elif i == 1:
                self.datafile = line.split()[1]
            elif i == 2:
                self.molecule = line.split()[1]
            elif i == 3:
                self.spacing = float(line.split()[1])
            elif i == 4:
                self.npts = [int(x) for x in line.split()[1:]]
            elif i == 5:
                self.center = [float(x) for x in line.split()[1:]]
        for i in range(3):
            self.n[i] = self.npts[i] + 1
        self.nelem = self.n[0] * self.n[1] * self.n[2]
        i = 0
        while i < self.nelem:
            val = float(fp.readline())
            self.values.append(val)
            i += 1
        for i in range(3):
            self.origin[i] = self.center[i] - self.npts[i] / 2 * self.spacing

    def meta(self):
        s = 'GRID_PARAMETER_FILE %s\n' % self.paramfile + \
           'GRID_DATA_FILE %s\n' % self.datafile +\
           'MACROMOLECULE %s\n' % self.molecule +\
           'SPACING %4.3f\n' % self.spacing +\
           'NELEMENTS %d %d %d\n' % (self.npts[0], self.npts[1], self.npts[2]) +\
           'CENTER %5.3f %5.3f %5.3f\n' % (self.center[0], self.center[1], self.center[2])
        return s

    def write(self, fp):
        print('GRID_PARAMETER_FILE %s' % self.paramfile, file=fp)
        print('GRID_DATA_FILE %s' % self.datafile, file=fp)
        print('MACROMOLECULE %s' % self.molecule, file=fp)
        print('SPACING %4.3f' % self.spacing, file=fp)
        print('NELEMENTS %d %d %d' % (self.npts[0], self.npts[1], self.npts[2]), file=fp)
        print('CENTER %5.3f %5.3f %5.3f' % (self.center[0], self.center[1], self.center[2]), file=fp)
        for x in self.values:
            if abs(x) < self.precision:
                print('0.', file=fp)
            else:
                print('%.3f' % x, file=fp)

    def writeDX(self, fname):
        fp = open(fname, 'w')
        nx = self.n[0]
        ny = self.n[1]
        nz = self.n[2]
        ori = self.origin
        spacing = self.spacing
        vals = self.values

        print('#==================================', file=fp)
        print('# AutoGrid Map File: %s' % self.name, file=fp)
        print('# Receptor File Name: %s' % self.molecule, file=fp)
        print('#==================================', file=fp)
        print('object 1 class gridpositions counts %d %d %d' % (nx, ny, nz), file=fp)
        print('origin %12.5E %12.5E %12.5E' % (ori[0], ori[1], ori[2]), file=fp)
        print('delta %12.5E %12.5E %12.5E' % (spacing, 0, 0), file=fp)
        print('delta %12.5E %12.5E %12.5E' % (0, spacing, 0), file=fp)
        print('delta %12.5E %12.5E %12.5E' % (0, 0, spacing), file=fp)
        print('object 2 class gridconnections counts %d %d %d' % (nx, ny, nz), file=fp)
        print('object 3 class array type double rank 0 items %d data follows' % len(vals), file=fp)
        for k in range(nz):
            col = 0
            for j in range(ny):
                for i in range(nx):
                    fp.write(" %12.5E" % vals[i * ny * nz + j * nz + k])
                    col += 1
                    if col == 3:
                        print(file=fp)
                        col = 0
        print('attribute \"dep\" string \"positions\"', file=fp)
        print('object \"regular positions regular connections\" class field', file=fp)
        print('component \"positions\" value 1', file=fp)
        print('component \"connections\" value 2', file=fp)
        print('component \"data\" value 3', file=fp)
        fp.close()

#==========================================================================
#
#    CLASSES FOR INTERNAL HANDLING OF RECEPTORS AND LIGANDS


class Receptor:

    """CONTAINS ALL INFORMATION ABOUT A DEFINED RECEPTOR"""

    def __init__(self):
        self.selection = ''
        self.pdb_file = ''
        self.receptor_pdbqt = ''
        self.receptor_rigid = ''
        self.receptor_flexible = ''
        self.flexible_residues = {}
        self.resi_dic = {}

    def flex_res_string(self):
        flex_res_str = ''
        for key, val in self.flexible_residues.items():
            s = os.path.basename(self.receptor_pdbqt)[:-6]
            lst = []
            for idx, resname in val:
                lst.append(resname + str(idx))
            s += ':' + key + ':' + '_'.join(lst)
            flex_res_str += s
        return flex_res_str
##         lst = []
# for idx, resname in self.flexible_residues:
##             lst.append( resname+str(idx) )
# return '_'.join(lst)

    def info(self):
        s = '#===============================================\n'
        s += ' > Receptor                 : "%s"\n' % self.selection
        s += ' > Generated from pdb file  : "%s"\n' % self.pdb_file
        s += ' > Receptor file            : "%s"\n' % self.receptor_pdbqt
        s += ' > Receptor rigid           : "%s"\n' % self.receptor_rigid
        s += ' > Receptor flexible        : "%s"\n' % self.receptor_flexible
#        s+=' > Number of flex. residues : %d\n' % len(self.flexible_residues)
        s += '#===============================================\n'
        return s


class Ligand:

    """CONTAINS ALL INFORMATION ABOUT A LIGAND"""

    def __init__(self):
        self.name = ''
        self.selection = ''
        self.input_file = ''
        self.ligand_pdbqt = ''
        self.outfile_poses = ''

    def info(self):
        s = '#=============================================\n'
        s += ' > Ligand                  : %s\n' % self.name
        s += ' > Generated from file     : %s\n' % self.input_file
        s += ' > Ligand pdbqt file       : %s\n' % self.ligand_pdbqt
        s += ' > Poses output file       : %s\n' % self.outfile_poses
        s += '#=============================================\n'
        return s


class MyFileDialog:

    def __init__(self, types=[("Executable", "*")], initialdir: str = None):
        self.initialdir = initialdir
        self.types = types

    def getopenfile(self):
        result = tkFileDialog.askopenfilename(initialdir = self.initialdir, filetypes=self.types)
        if result == "":
            return None
        else:
            return result

    def getsavefile(self):
        result = tkFileDialog.asksaveasfilename(initialdir = self.initialdir, filetypes=self.types)
        if result == "":
            return None
        else:
            return result