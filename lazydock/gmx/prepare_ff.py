'''
Date: 2025-06-12 16:42:14
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-06-29 13:10:55
Description: 
'''
import os
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Union

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import opts_file
from pymol import cmd

from lazydock.config import GlobalConfig
from lazydock.gmx.run import Gromacs
from lazydock.pml.align_to_axis import align_pose_to_axis


def insert_content(content: str, before: str, new_content: str):
    is_file, path = False, None
    if os.path.isfile(content):
        is_file, path = True, content
        content = opts_file(content)
    idx1 = content.find(before)
    if idx1 == -1:
        return put_err(f'{before} not found, skip.')
    content = content[:idx1+len(before)] + new_content + content[idx1+len(before):]
    if is_file:
        opts_file(path, 'w', data=content)
    return content


def fix_name_in_mol2(ipath: str, opath: str):
    lines = opts_file(ipath, way='lines')
    lines[1] = 'LIG\n'
    get_idx_fn = lambda content, offset: lines.index(list(filter(lambda x: x.startswith(content), lines[offset:]))[0])  # pyright: ignore[reportAttributeAccessIssue]
    atom_st = get_idx_fn('@<TRIPOS>ATOM', 0)+1
    atom_ed = get_idx_fn('@<TRIPOS>', atom_st+1)
    for i in range(atom_st, atom_ed):
        resn = lines[i].split()[7].strip()  # pyright: ignore[reportAttributeAccessIssue]
        resn_idx = lines[i].index(resn)  # pyright: ignore[reportAttributeAccessIssue]
        lines[i] = lines[i][:resn_idx] + 'LIG' + ' '*(min(0, len(resn) - 3)) + lines[i][resn_idx+len(resn):]
    opts_file(opath, 'w', way='lines', data=lines)
    

def prepare_ff_dir(w_dir: Path, ff_dir: str):
    if ff_dir is None:
        return put_err('ff-dir is None, skip transform.')
    if os.path.exists(ff_dir):
        _ff_dir = w_dir / Path(ff_dir).name
        if not _ff_dir.exists():
            shutil.copytree(os.path.abspath(ff_dir), _ff_dir, dirs_exist_ok=True)
            put_log(f'copy {ff_dir} to {_ff_dir}')
        else:
            put_log(f'ff-dir(repeat in sub-directory) already exists in {_ff_dir}, skip.')
    ## if ff-dir already in each sub-directory, do not overwrite it.
    elif os.path.exists(w_dir / ff_dir):
        _ff_dir = w_dir / ff_dir
        put_log(f'ff-dir already exists in (sub directory) {_ff_dir}, skip.')
    else:
        _ff_dir = None
        put_err(f'cannot find ff_dir: {ff_dir} in {w_dir}, set ff_dir to None, skip.')
    return _ff_dir


def prepare_complex_topol(ipath_rgro: str, ipath_lgro: str, ipath_top: str, opath_cgro: str, opath_top: str,
                            insert_itp: bool = True, insert_prm: bool = True):
    # merge receptor and ligand gro into complex.gro
    receptor_gro_lines = list(filter(lambda x: len(x.strip()), opts_file(ipath_rgro, 'r', way='lines')))  # pyright: ignore[reportAttributeAccessIssue]
    lig_gro_lines = list(filter(lambda x: len(x.strip()), opts_file(ipath_lgro, 'r', way='lines')))  # pyright: ignore[reportAttributeAccessIssue]
    complex_gro_lines = receptor_gro_lines[:-1] + lig_gro_lines[2:-1] + receptor_gro_lines[-1:]
    complex_gro_lines[1] = f'{int(receptor_gro_lines[1]) + int(lig_gro_lines[1])}\n'
    opts_file(opath_cgro, 'w', way='lines', data=complex_gro_lines)
    # inset ligand paramters in topol.top
    topol = opts_file(ipath_top)
    if (Path(ipath_top).parent / 'lig.itp').exists() and insert_itp:
        _topol = complex.insert_content(topol, '#include "posre.itp"\n#endif\n',  # pyright: ignore[reportAttributeAccessIssue]
                                    '\n; Include ligand topology\n#include "lig.itp"\n')
        if _topol is not None:
            topol = _topol
    if (Path(ipath_top).parent / 'lig.prm').exists() and insert_prm:
        _topol = complex.insert_content(topol, '#include "./charmm36-jul2022.ff/forcefield.itp"\n',  # pyright: ignore[reportAttributeAccessIssue]
                                    '\n; Include ligand parameters\n#include "lig.prm"\n')
        if _topol is not None:
            topol = _topol
    topol += 'LIG                 1\n'  # pyright: ignore[reportOperatorIssue]
    opts_file(opath_top, 'w', data=topol)


def abs_path(path: str):
    return str(Path(path).resolve().absolute())


def prepare_complex_with_sobtop(complex_path: str, ff_dir: str, rec_chain: str, lig_chain: str,
                                pdb2gmx: str = "-ter -ignh", n_term: str = '0', c_term: str = '0'):
    """root: the root contains sub-folders of ligand.mol and protein.pdb files"""
    complex_dir = os.path.dirname(complex_path)
    gmx = Gromacs(working_dir=abs_path(complex_dir))
    
    # check if the complex.gro file exists
    if os.path.exists(abs_path(os.path.join(complex_dir, 'complex.gro'))):
        return put_log(f'complex.gro already exists in {complex_dir}, skip.')
    
    # check if the sobtop exists
    if not os.path.exists(GlobalConfig.named_paths["sobtop_dir"]):
        return put_err(f'sobtop not found from ~/.lazydock/config.json, skip.')
    
    if os.path.exists(abs_path(os.path.join(complex_dir, 'FAILED'))):
        return put_log(f'FAILED already exists in {complex_dir}, skip.')
    
    cmd.reinitialize()
    
    # STEP 0.1: center complex.pdb by obabel.
    ipath, opath = os.path.join(complex_dir, 'complex.pdb'), os.path.join(complex_dir, f'0.1_complex_ori_centered.pdb')
    os.system(f'obabel -ipdb {ipath} -opdb -O {opath} -c')
        
    # step 0.2: align complex.pdb with xyz axes by lazydock.
    ipath, opath = opath, os.path.join(complex_dir, f'0.2_complex_center_align_axis.pdb')
    cmd.reinitialize()
    cmd.load(ipath, 'complex')
    align_pose_to_axis('complex')
    cmd.save(opath, 'complex')
        
    # STEP 1: extract receptor and ligand from complex.pdb.
    ipath, opath_r, opath_l = opath, os.path.join(complex_dir, f'1_complex_receptor.pdb'), os.path.join(complex_dir, f'1_complex_ligand.pdb')
    if cmd.select('rec', f'chain {rec_chain}') == 0:
        raise ValueError(f'No receptor chains found for {complex_path}.')
    if cmd.select('lig', f'chain {lig_chain}') == 0:
        raise ValueError(f'No ligand chains found for {complex_path}.')
    cmd.save(opath_r, 'rec')
    cmd.save(opath_l, 'lig')
        
    # STEP 2: transfer ligand.pdb to mol2 by obabel.
    ipath, opath = opath_l, os.path.join(complex_dir, f'2_complex_ligand.mol2')
    os.system(f'obabel -ipdb {ipath} -omol2 -O {opath}')
    # STEP 3: fix ligand name and residue name in mol2 file.
    ipath, opath = opath, os.path.join(complex_dir, f'3_complex_ligand_named.mol2')
    fix_name_in_mol2(ipath, opath)
    # STEP 4: delete UNITY_ATOM_ATTR TRIPOS.
    ipath, opath = opath, os.path.join(complex_dir, f'4_complex_ligand_sobtop.mol2')
    lig_mol2_lines = opts_file(ipath).split('\n')  # pyright: ignore[reportAttributeAccessIssue]
    try:
        attr_line_id = lig_mol2_lines.index('@<TRIPOS>UNITY_ATOM_ATTR')
        next_attr_id = [i for i, line in enumerate(lig_mol2_lines[attr_line_id+1:]) if line.startswith('@<TRIPOS>')][0]+attr_line_id+1
        opts_file(opath, 'w', data='\n'.join(lig_mol2_lines[:attr_line_id]+lig_mol2_lines[next_attr_id:]))
    except:
        shutil.copy(ipath, opath)

    # STEP 5: prepare the ligand topology.
    ipath, opath_lgro, opath_ltop, opath_litp = opath, os.path.join(complex_dir, f'lig.gro'), os.path.join(complex_dir, f'lig.top'), os.path.join(complex_dir, f'lig.itp')
    gmx.run_command_with_expect(f'cd {GlobalConfig.named_paths["sobtop_dir"]} && ./sobtop',
                                [{'.mol2 should be used': f'{abs_path(ipath)}\r'},
                                 {'10 Other functions': '2\r'},
                                 {'current folder': f'{abs_path(opath_lgro)}\r'},
                                 {'10 Other functions': '1\r'},
                                 {'9 Show all current atom types on screen': '2\r'},
                                 {'by mSeminario method': '4\r'},
                                 {'current folder': f'{abs_path(opath_ltop)}\r'},
                                 {'current folder': f'{abs_path(opath_litp)}\r'},
                                 {'10 Other functions': '0\r'},])
    
    # STEP 6: fix ligand name in lig.itp.
    lig_itp = opts_file(opath_litp).replace('4_complex_ligand_sobtop', 'LIG').replace('MOL', 'LIG')  # pyright: ignore[reportAttributeAccessIssue]
    opts_file(opath_litp, 'w', data=lig_itp)
    
    # STEP 7: Prepare the forcefiled directory.
    prepare_ff_dir(Path(complex_dir), ff_dir)
    
    # STEP 8: Prepare the Protein Topology
    ipath, opath_rgro = opath_r, os.path.join(complex_dir, f'complex_receptor.gro')
    ## fix pdb residues heavy atoms by pdbfixer
    gmx.run_cmd_with_expect(f'pdbfixer 1_complex_receptor.pdb --output 1_complex_receptor.pdb --add-atoms heavy --add-residues')
    ## convert pdb to gro
    if not os.path.exists(opath_rgro):
        expect_acts = [{'dihedrals)': '1\r'}, {'None': '1\r'}]
        # assign n-term and c-term to expect_acts.
        expect_acts.append({'None': f'{n_term}\r'})
        expect_acts.append({'None': f'{c_term}\r'})
        # run pdb2gmx
        gmx.run_gmx_with_expect(f'pdb2gmx -f {Path(ipath).name} -o {Path(opath_rgro).name} {pdb2gmx}', expect_acts)
    if not os.path.exists(opath_rgro):
        opts_file(os.path.join(complex_dir, f'FAILED'), 'w', data='pdb2gmx failed.')
        return put_err(f'pdb2gmx failed for {complex_path}.')
    # STEP 9: add receptor position restraints?
    
    # STEP 10: Prepare the Complex Topology
    opath_cgro, opath_top = os.path.join(complex_dir, 'complex.gro'), os.path.join(complex_dir, 'topol.top')
    prepare_complex_topol(opath_rgro, opath_lgro, opath_top, opath_cgro, opath_top, False, False)
    ## inset ligand paramters in topol.top
    ## lig.itp contains atomtypes, so put it before any [moleculetype] directive.
    ## https://manual.gromacs.org/current/user-guide/run-time-errors.html#invalid-order-for-directive-xxx
    topol = opts_file(opath_top)
    ff_fir_name = os.path.basename(ff_dir)
    topol = insert_content(topol, f'#include "./{ff_fir_name}/forcefield.itp"',
                                '\n\n; Include ligand topology\n#include "lig.itp"\n')
    opts_file(opath_top, 'w', data=topol)


__all__ = [
    'insert_content',
    'fix_name_in_mol2',
    'prepare_ff_dir',
    'prepare_complex_topol',
    'abs_path',
    'prepare_complex_with_sobtop',
]