'''
Date: 2026-10-09
LastEditors: BHM-Bob
Description: meeko based receptor/ligand preparation utilities.

This module is a LOW-LEVEL building block: it must not import other lazydock.pml
modules, so that `lazydock.pml.adfr` (and anything else) can depend on it without
creating an import cycle.
'''
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file

# flexible residue is represented as (chain, resname, resnum) in this module,
# resname could be empty string when user only provide residue number.
FlexResidue = Tuple[str, str, int]


# first meeko release where --flexres takes one comma-separated sequence
MEEKO_FLAT_FLEXRES_VERSION = (0, 7, 0)


def _get_resname_map(pdbstr: str) -> dict:
    """return {(chain, resnum): resname} from pdb string"""
    result = {}
    for line in pdbstr.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            try:
                result[(line[21], int(line[22:26]))] = line[17:20].strip()
            except ValueError:
                continue
    return result


def parse_flex_residues(flex_residues: Union[str, List[str]]) -> Optional[List[FlexResidue]]:
    """
    parse flexible residue strings into a list of (chain, resname, resnum).
    
    supported formats:
        - '<chain>:<resnum>,<resnum>'       e.g. 'A:68,50'  or 'A:68;A:50'
        - '<chain>:<resname><resnum>,...'   e.g. 'A:LYS68,LYS50' (agfr style)
        - '<chain>:<resname><resnum>_...'   e.g. 'D:PHE108_PHE177' (lazy_pocket style)
        - '<chain>_<resname><resnum>'       e.g. 'A_LYS68' or 'A_68'
        - '<chain>_<resname>-<resnum>'      e.g. 'A_LYS-68'
    
    ';' separates chain groups, ',' and '_' separate residues inside a group.
    
    :param flex_residues: str or list of str, flexible residues
    :return: list of (chain, resname, resnum), resname could be '' if not provided; None if failed
    """
    if isinstance(flex_residues, str):
        flex_residues = [flex_residues]
    result: List[FlexResidue] = []
    # ';' and ',' separate tokens; a token containing ':' always starts a new chain group
    for raw in flex_residues:
        chain = None
        for token in re.split(r'[;,]', raw):
            token = token.strip()
            if not token:
                continue
            if ':' in token:
                chain, rest = token.split(':', 1)
                chain = chain.strip()
                if not chain:
                    return put_err(f'invalid flexible residue token: {token}, chain is empty.')
            elif chain is not None:
                rest = token
            elif '_' in token:
                # '<chain>_<resname><resnum>' style
                chain, rest = token.split('_', 1)
                chain = chain.strip()
            else:
                return put_err(f'invalid flexible residue token: {token}, chain is missing '
                               f'(expected "<chain>:<residues>" or "<chain>_<residues>")')
            # residues inside a group could be separated by '_' or whitespace
            for res_token in re.split(r'[_\s]+', rest.strip()):
                res_token = res_token.strip()
                if not res_token:
                    continue
                # split resname (leading letters) and resnum (trailing digits)
                idx = 0
                while idx < len(res_token) and not res_token[idx].isdigit():
                    idx += 1
                resname, resnum_str = res_token[:idx].strip(' -'), res_token[idx:].strip()
                if not resnum_str or not resnum_str.isdigit():
                    return put_err(f'invalid flexible residue token: {res_token} '
                                   f'(from "{token}"), residue number should be integer, got "{resnum_str}".')
                result.append((chain, resname.upper(), int(resnum_str)))
    if not result:
        return put_err(f'can not parse any flexible residue from: {flex_residues}')
    return result


def get_meeko_version(caller: Optional[str] = None) -> Optional[Tuple[int, ...]]:
    """
    detect the installed meeko version, using the given python interpreter.

    meeko changed its `--flexres` syntax incompatibly: <=0.6 grouped residues per
    chain and joined the groups with ';' (e.g. 'A:78,79;T:22,24'), while >=0.7
    expects ONE comma-separated sequence with the chain prefix repeated for every
    group (e.g. 'A:78,79,T:22,24').  See `flex_residues_to_meeko` for details.

    :param caller: python interpreter to probe, None means the current one
    :return: version tuple such as (0, 8, 0), or None if meeko can not be probed
    """
    caller = caller or sys.executable
    code = ('import meeko; '
            'print(getattr(meeko, "__version__", "0.0.0"))')
    try:
        proc = subprocess.run([caller, '-c', code], stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True, timeout=120)
    except (subprocess.TimeoutExpired, OSError):
        return None
    if proc.returncode != 0:
        return None
    raw = (proc.stdout or '').strip().splitlines()
    if not raw:
        return None
    m = re.match(r'(\d+)\.(\d+)(?:\.(\d+))?', raw[-1].strip())
    if not m:
        return None
    return tuple(int(x) for x in m.groups(default='0'))


def flex_residues_to_meeko(flex_residues: List[FlexResidue],
                           meeko_version: Optional[Tuple[int, ...]] = None) -> str:
    """
    convert flexible residues to the meeko `--flexres` string.

    The syntax differs between meeko generations, so the caller should pass the
    detected version whenever several interpreters may be in play:

        - meeko < 0.7  : '<chain>:<resnum>,<resnum>;<chain2>:...'
                         e.g. [(A, '', 78), (A, '', 79), (T, '', 22)] -> 'A:78,79;T:22'
        - meeko >= 0.7 : '<chain>:<resnum>,<resnum>,<chain2>:<resnum>,...'
                         e.g. [(A, '', 78), (A, '', 79), (T, '', 22)] -> 'A:78,79,T:22'

    meeko only accepts residue numbers, resname is ignored.  Passing the older
    ';' form to meeko >= 0.7 makes it read the whole tail as a chain name and fail
    with `ValueError: residue_id='100;T:51' not in valid monomers`.

    :param meeko_version: detected meeko version, None -> treat as the new syntax
                          (see `get_meeko_version` for how it is obtained)
    """
    chains: dict = {}
    for chain, _, resnum in flex_residues:
        chains.setdefault(chain, []).append(resnum)
    groups = [f'{c}:{",".join(str(n) for n in sorted(set(v)))}'
              for c, v in chains.items()]
    if meeko_version is not None and meeko_version < MEEKO_FLAT_FLEXRES_VERSION:
        return ';'.join(groups)
    return ','.join(groups)


def flex_residues_to_agfr(flex_residues: List[FlexResidue], pdbstr: str = None) -> str:
    """
    convert flexible residues to agfr format: '<chain>:<RESN><resnum>,<RESN><resnum>;<chain2>:...'
    
    agfr requires resname, if resname is missing in flex_residues, it will be
    looked up from pdbstr.
    """
    resname_map = _get_resname_map(pdbstr) if pdbstr else {}
    chains: dict = {}
    for chain, resname, resnum in flex_residues:
        if not resname:
            resname = resname_map.get((chain, resnum), '')
        if not resname:
            put_err(f'can not determine resname for flexible residue {chain}:{resnum}, '
                    f'agfr requires resname.', _exit=True)
        chains.setdefault(chain, []).append(f'{resname}{resnum}')
    return ';'.join(f'{c}:{",".join(sorted(set(v)))}' for c, v in chains.items())


def prepare_receptor_with_flexres(tmpdir: str, receptor_pdbstr: str, flex_residues: List[FlexResidue],
                                  stress: float = None, bufsize: float = None,
                                  python_caller: str = None, timeout: int = 3600):
    """
    call meeko to prepare receptor pdbqt files: rigid part and flexible sidechains.
    
    meeko will add missing polar hydrogens automatically, so input receptor could
    be a hydrogen-free pdb string.
    
    :param tmpdir: working dir
    :param receptor_pdbstr: receptor pdb string
    :param flex_residues: flexible residues, list of (chain, resname, resnum)
    :param stress: grid box size (x, y, z), used to write gpf file by meeko
    :param bufsize: grid box center (x, y, z), used to write gpf file by meeko
    :param python_caller: python interpreter to call meeko, None means `sys.executable`
    :param timeout: subprocess timeout in seconds
    :return: (rigid_pdbqt_path, flex_pdbqt_path), None if failed
    """
    python_caller = python_caller or sys.executable
    rec_path = os.path.join(tmpdir, 'receptor.pdb')
    opts_file(rec_path, 'w', way='str', data=receptor_pdbstr)
    # meeko only accepts '<chain>:<resnum>' format, and the separator between
    # chain groups changed at 0.7 (see `flex_residues_to_meeko`): probe the
    # interpreter we are about to call, not the one running this code.
    meeko_version = get_meeko_version(python_caller)
    flex_str = flex_residues_to_meeko(flex_residues, meeko_version)
    cmd_list = [python_caller, '-m', 'meeko.cli.mk_prepare_receptor',
                '--read_pdb', rec_path, '-o', 'receptor', '-p',
                # the flex string may contain ';' (meeko < 0.7) or ',' (>= 0.7);
                # the command is joined with shell=True below, so it MUST be
                # quoted -- an unquoted ';' is a shell separator and the tail of
                # the list is then executed as a command (rc=127), while even the
                # quoted ';' still needs shlex to survive `' '.join`.
                '-f', shlex.quote(flex_str)]
    if stress and bufsize:
        # meeko writes gpf only when box is provided, but we compute maps by ourself
        cmd_list += ['-g', '-v',
                     '--box_size', f'{stress[0]:.3f}', f'{stress[1]:.3f}', f'{stress[2]:.3f}',
                     '--box_center', f'{bufsize[0]:.3f}', f'{bufsize[1]:.3f}', f'{bufsize[2]:.3f}']
    try:
        proc = subprocess.run(' '.join(cmd_list), shell=True, cwd=tmpdir,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return put_err(f'meeko prepare_receptor timeout after {timeout} seconds')
    output_str = proc.stdout + proc.stderr
    rigid_path = os.path.join(tmpdir, 'receptor_rigid.pdbqt')
    flex_path = os.path.join(tmpdir, 'receptor_flex.pdbqt')
    if proc.returncode != 0 or not Path(rigid_path).exists():
        return put_err(f'meeko prepare_receptor failed (returncode={proc.returncode}, '
                       f'meeko_version={meeko_version}, flexres="{flex_str}"):\n{output_str}')
    if not Path(flex_path).exists():
        return put_err(f'meeko did not generate flexible receptor file, '
                       f'please check flexible residues {flex_str}:\n{output_str}')
    return rigid_path, flex_path