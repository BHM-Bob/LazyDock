'''
Date: 2024-11-20 21:45:38
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-11-23 19:46:20
Description: transfer from http://www.mdtutorials.com/gmx/complex/Files/sort_mol2_bonds.pl
'''

import os
import sys

from mbapy_lite.file import opts_file


def sort_bonds(mol2: str) -> str:
    if os.path.exists(mol2):
        mol2_lines = opts_file(mol2).split('\n')
    else:
        mol2_lines = mol2.split('\n')
    # Check for header lines
    if not mol2_lines[0].startswith('@<TRIPOS>'):
        sys.exit(f"Nonstandard header found: {mol2_lines[0]}. Please delete header lines until the TRIPOS molecule definition.")
    # output
    result_lines = []
    # Get number of atoms and bonds from mol2 file
    natom, nbond = map(int, mol2_lines[2].strip().split()[:2])
    print(f"Found {natom} atoms in the molecule, with {nbond} bonds.")
    # Print out everything up until the bond section
    i = 0
    while not mol2_lines[i].startswith('@<TRIPOS>BOND'):
        result_lines.append(mol2_lines[i])
        i += 1
    # Print the bond section header line to output
    result_lines.append(mol2_lines[i])
    i += 1
    # Read in the bonds and sort them
    bondfmt = "%6d%6d%6d%5s"
    tmparray = []
    for j in range(nbond):
        parts = mol2_lines[i + j].split()
        ai, aj = int(parts[1]), int(parts[2])
        # Reorder if second atom number < first
        if aj < ai:
            ai, aj = aj, ai
        tmparray.append((ai, aj, int(parts[0]), parts[3]))
    # Sort the bonds
    tmparray.sort()
    # Loop over tmparray to find each atom number
    nbond = 0
    for x in range(1, natom + 1):
        bondarray = [t for t in tmparray if t[0] == x]
        if bondarray:
            for y, (ai, aj, bi, bj) in enumerate(sorted(bondarray, key=lambda x: x[1])):
                nbond += 1
                result_lines.append(bondfmt % (bi, ai, aj, bj))
    return '\n'.join(result_lines)


if __name__ == "__main__":
    # dev code
    sorted_mol2 = sort_bonds('data_tmp/pdb/ligand.mol2')
    opts_file('data_tmp/pdb/ligand_fix.mol2', 'w', data = sorted_mol2)
    print(sorted_mol2)