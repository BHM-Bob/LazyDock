'''
Date: 2026-10-09
LastEditors: BHM-Bob
Description: meeko based receptor/ligand preparation utilities.

This module is a LOW-LEVEL building block: it must not import other lazydock.pml
modules, so that `lazydock.pml.adfr` (and anything else) can depend on it without
creating an import cycle.

Contents:
    - rigid receptor preparation by meeko Python API (`prepare_receptor_with_meeko`)
    - flexible receptor preparation by meeko CLI (`prepare_receptor_with_flexres`)
    - ligand preparation by meeko Python API (`prepare_ligand_with_meeko`)
    - macrocycle / glue bookkeeping (`inspect_macrocycle`, `strip_glue_atoms`)
    - flexible residue parsing and formatting helpers
'''
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file
from meeko import MoleculePreparation, PDBQTWriterLegacy
from meeko.polymer import Polymer, ResidueChemTemplates

# flexible residue is represented as (chain, resname, resnum) in this module,
# resname could be empty string when user only provide residue number.
FlexResidue = Tuple[str, str, int]


# first meeko release where --flexres takes one comma-separated sequence
MEEKO_FLAT_FLEXRES_VERSION = (0, 7, 0)


# meeko breaks rings only when min_ring_size <= size <= max_ring_size; larger rings
# are SILENTLY ignored (no error, no warning) and stay rigid.  The default 33 misses
# most cyclic peptides: a measured scan gives the largest main ring at 56 atoms, and
# a per-molecule probe shows the threshold reacts sharply (a 48-mer ring is still
# rigid at max_ring_size=45 and only breaks at 48).  Hence 60, with headroom.
MACROCYCLE_MAX_RING_SIZE = 60

# meeko's ring breaking is ON BY DEFAULT for the whole 7..33 atom window, which is
# easy to miss: it does NOT only affect big macrocycles but also plain medium rings
# (a measured probe opens 7-, 10-, 12- and 20-membered rings with stock settings,
# while benzene stays closed only because the lower bound is 7).  That silently
# changes a ligand that ADT's `prepare_ligand` would deliver as a rigid ring, and
# the glue pseudo-atoms it adds (G0/CG0) have no grid in an AutoDock map.  Passing
# `rigid_macrocycles=True` makes every ring rigid, which is the ADT-equivalent
# behaviour; meeko documents it as "every ring is rigid without macrocycle option".
RIGID_MACROCYCLES = True

# disulfide bonds are covalent cross-links, so they must not become freely rotating
# torsions once the ring opens.  meeko matches the SMARTS and freezes the bond
# between its atoms: the indices are POSITIONS INSIDE THE SMARTS (0-based here,
# 1-based on the CLI), NEVER atom serial numbers.
# A SMARTS with no match is harmless, so this can be passed unconditionally.
SS_RIGIDIFY_SMARTS = ['[SX2]-[SX2]']
SS_RIGIDIFY_INDICES = [(0, 1)]

# macrocycle glue pseudo-atoms are named "G" with zero radius/charge (atom types
# G0..G3).  "G" is not an atom name of any standard residue, so dropping it needs no
# opt-in; it MUST happen before a pose reaches a pdb consumer: `ADModel(parse2std=True)`
# re-derives the element from the atom NAME and writes it into a 2-char column, so a
# glue atom would reach the output pdb as a fake element "G" and break
# PDBFixer/OpenMM and any re-preparation of the pose.
GLUE_ATOM_NAMES = ('G',)


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


def prepare_receptor_with_meeko(receptor_pdbstr: str, flex_residues: List[FlexResidue] = None,
                                allow_bad_res: bool = True):
    """
    prepare a receptor pdbqt string with meeko, without touching the filesystem.

    This is the in-process Python API path, so meeko must be importable by the
    CURRENT interpreter; the CLI twin `prepare_receptor_with_flexres` exists for the
    cases where a separate interpreter has to be called.

    Measured byte-level equivalence with `python -m meeko.cli.mk_prepare_receptor`
    for both the rigid-only and the flexible-residue cases (97440 bytes and 1687
    bytes respectively, identical output), so the two paths are interchangeable.

    NOTICE: meeko's atom typing is NOT identical to ADT's -- on a real receptor the
    heavy atoms match one-by-one but 2 nitrogen acceptors are classified differently
    (meeko: NA=4, N=164 / ADT: NA=2, N=166) and 4 polar hydrogens are typed
    differently (meeko HD=221 / ADT HD=225).  vina's H-bond and N/NA terms depend on
    exactly these types, so scores from the two backends are NOT comparable and
    existing ADT-prepared results can not be reproduced with meeko.

    :param receptor_pdbstr: receptor pdb string, hydrogens are added by meeko
    :param flex_residues: flexible residues as list of (chain, resname, resnum), None
        or empty means a fully rigid receptor.  meeko only needs the residue numbers,
        but the chain grouping follows `flex_residues_to_meeko`.
    :param allow_bad_res: keep residues that meeko can not match to a template, so
        that one unusual residue does not abort the whole preparation.  This mirrors
        the `-a/--allow_bad_res` flag the CLI twin does not pass and that the probes
        used for the equivalence measurement (which is why the default is True here).
    :return: (pdbqt_string, flex_residues) where the second item is the flexible part
        (empty string when `flex_residues` is None/empty), None if failed.
    """
    try:
        templates = ResidueChemTemplates.create_from_defaults()
        mk_prep = MoleculePreparation()
        polymer = Polymer.from_pdb_string(receptor_pdbstr, templates, mk_prep,
                                          allow_bad_res=allow_bad_res)
    except Exception as e:
        return put_err(f'meeko failed to build a polymer from the receptor: {e}')
    if flex_residues:
        # meeko keys its monomers as '<chain>:<resnum>', NOT as the (chain, resname,
        # resnum) tuples this module uses everywhere else
        monomer_ids = {f'{chain}:{resnum}' for chain, _, resnum in flex_residues}
        unknown = sorted(rid for rid in monomer_ids if rid not in polymer.monomers)
        if unknown:
            return put_err(f'flexible residues {unknown} are not present in the '
                           f'receptor, can not make them flexible')
        for rid in sorted(monomer_ids):
            polymer.flexibilize_sidechain(rid, mk_prep)
    try:
        rigid_pdbqt, flex_dict = PDBQTWriterLegacy.write_from_polymer(polymer)
    except Exception as e:
        return put_err(f'meeko failed to write the receptor pdbqt: {e}')
    # one entry per flexible residue, each holding that sidechain's atoms; only the
    # non-empty ones are real (meeko also lists rotatable terminal groups here)
    flex_pdbqt = '\n'.join(v for v in flex_dict.values() if v)
    return rigid_pdbqt, flex_pdbqt


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


def pdbstr_to_rdkit_mol(pdbstr: str):
    """
    build an RDKit mol from a pdb string, ready for meeko ligand preparation.

    meeko's ligand readers take SDF/MOL2/MOL and require EXPLICIT hydrogens, while
    the pipeline hands us a pdb string, so this conversion is mandatory.  The bond
    orders are perceived from distances (`proximityBonding=True`), which is the
    only option for a pdb without CONECT records.

    :param pdbstr: ligand pdb string
    :return: rdkit.Chem.rdchem.Mol with explicit hydrogens, None if failed
    """
    try:
        from rdkit import Chem
    except ImportError:
        return put_err('rdkit is required for meeko ligand preparation, '
                       'please install it with `pip install rdkit`')
    try:
        mol = Chem.MolFromPDBBlock(pdbstr, removeHs=False,
                                   proximityBonding=True, sanitize=False)
    except Exception as e:
        return put_err(f'rdkit failed to parse ligand pdb: {e}')
    if mol is None:
        return put_err('rdkit failed to parse ligand pdb, it returned None')
    try:
        # SANITIZE_PROPERTIES is skipped because the perceived bond orders may not
        # satisfy the valence model; everything else (rings, aromaticity) is needed
        # by meeko to detect the macrocycle.
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_PROPERTIES)
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception as e:
        return put_err(f'rdkit failed to sanitize/add hydrogens to ligand: {e}')
    return mol


def inspect_macrocycle(pdbqtstr: str) -> Dict[str, Union[int, bool]]:
    """
    summarise the macrocycle state of a prepared ligand pdbqt string.

    meeko silently ignores rings larger than `max_ring_size` (no error, no warning)
    and then emits a ligand whose ring is frozen -- the dock still runs and still
    returns scores, so a silent miss is indistinguishable from success unless the
    product is inspected (see the module-level MACROCYCLE_MAX_RING_SIZE comment).

    :param pdbqtstr: ligand pdbqt string
    :return: dict with 'n_glue' (number of glue pseudo-atoms, 2 per broken bond),
        'n_glue_pairs' (n_glue // 2), 'n_branch' (BRANCH records), 'n_cg' (atoms
        retyped CG0..CG3), 'torsdof' (TORSDOF value) and 'ring_flexible'
        (True when any glue atom is present)
    """
    n_glue, n_cg, n_branch, torsdof = 0, 0, 0, 0
    for line in pdbqtstr.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            name = line[12:16].strip()
            atype = line[77:].strip() if len(line) > 77 else ''
            if name in GLUE_ATOM_NAMES:
                n_glue += 1
            if atype.startswith('CG'):
                n_cg += 1
        elif line.startswith('BRANCH'):
            n_branch += 1
        elif line.startswith('TORSDOF'):
            try:
                torsdof = int(line.split()[1])
            except (IndexError, ValueError):
                pass
    return dict(n_glue=n_glue, n_glue_pairs=n_glue // 2, n_branch=n_branch,
                n_cg=n_cg, torsdof=torsdof, ring_flexible=n_glue > 0)


def prepare_ligand_with_meeko(ligand_pdbstr: str,
                              max_ring_size: int = MACROCYCLE_MAX_RING_SIZE,
                              rigid_ss: bool = True, break_rings: bool = True):
    """
    prepare a ligand pdbqt with meeko, optionally breaking rings open.

    This is the in-process Python API path (the receptor path still shells out to
    the meeko CLI), so meeko must be importable by the CURRENT interpreter; there is
    no `python_caller` knob here on purpose.

    Ring breaking works by cutting one ring bond and attaching glue pseudo-atoms,
    which a vina >= 1.2 linear attraction potential pulls back together during the
    search -- the ring thus becomes a normal torsion tree and its conformation can
    adapt to the pocket.

    :param ligand_pdbstr: ligand pdb string (explicit hydrogens are added here)
    :param max_ring_size: largest ring meeko may break, default 60 (see
        MACROCYCLE_MAX_RING_SIZE), ignored when `break_rings` is False
    :param rigid_ss: freeze S-S bonds with `SS_RIGIDIFY_SMARTS`, otherwise an opened
        ring makes the disulfide a freely rotating torsion, which is unphysical
    :param break_rings: break rings at all.  PASS False FOR EVERY NON-MACROCYCLE RUN:
        meeko's ring breaking is ON BY DEFAULT and acts on the whole 7..33 atom
        window, so leaving it on silently opens plain medium rings (measured: 7-, 10-,
        12- and 20-membered rings all get glue atoms) that ADT's `prepare_ligand`
        would deliver rigid -- see the module-level RIGID_MACROCYCLES comment.
    :return: (pdbqt_string, macrocycle_info), None if failed.  The caller MUST check
        `macrocycle_info['ring_flexible']` when a flexible ring was intended.
    """
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError:
        return put_err('meeko is required for macrocycle ligand preparation, '
                       'please install it with `pip install meeko gemmi`')
    mol = pdbstr_to_rdkit_mol(ligand_pdbstr)
    if mol is None:
        return None
    if break_rings:
        kwargs = dict(max_ring_size=max_ring_size)
    else:
        # meeko documents this as "every ring is rigid without macrocycle option";
        # it short-circuits before min/max_ring_size are even looked at
        kwargs = dict(rigid_macrocycles=RIGID_MACROCYCLES)
    if rigid_ss:
        kwargs.update(rigidify_bonds_smarts=SS_RIGIDIFY_SMARTS,
                      rigidify_bonds_indices=SS_RIGIDIFY_INDICES)
    try:
        preparator = MoleculePreparation(**kwargs)
        molsetup_list = preparator.prepare(mol)
    except Exception as e:
        return put_err(f'meeko failed to prepare ligand (kwargs={kwargs}): {e}')
    if not molsetup_list:
        return put_err('meeko returned no molsetup for the ligand')
    if len(molsetup_list) > 1:
        # only reactive docking with several matches yields multiple setups, which
        # this path does not use; taking the first would silently drop the others
        return put_err(f'meeko returned {len(molsetup_list)} molsetups, expected 1')
    pdbqtstr, success, error_msg = PDBQTWriterLegacy.write_string(molsetup_list[0])
    if not success:
        return put_err(f'meeko failed to write ligand pdbqt: {error_msg}')
    return pdbqtstr, inspect_macrocycle(pdbqtstr)


def strip_glue_atoms(pdbstr: str) -> str:
    """
    drop macrocycle glue pseudo-atoms (element "G") from a pdb/pdbqt string.

    meeko names its glue pseudo-atoms "G", and `ADModel(parse2std=True)` re-derives
    the element from the atom NAME, so a glue atom survives into the output pdb as a
    fake element "G" -- which breaks PDBFixer/OpenMM and any later re-preparation of
    the pose from that pdb.  The atoms they were attached to (retyped CG0..CG3) are
    REAL atoms and are kept.

    "G" is not the atom name of any standard residue, so this is safe to run on any
    pdb and is idempotent.

    :param pdbstr: pdb or pdbqt string
    :return: the same string without glue pseudo-atom lines
    """
    kept = []
    for line in pdbstr.split('\n'):
        if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() in GLUE_ATOM_NAMES:
            continue
        kept.append(line)
    return '\n'.join(kept)