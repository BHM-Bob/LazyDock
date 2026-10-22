import re
from typing import List, NamedTuple, Optional, Tuple, Union

from mbapy_lite.base import put_err
from pymol import cmd

from lazydock.dock.config import MEEKO_FLAT_FLEXRES_VERSION
from lazydock.pml.utils import normalize_receptor_chain

# flexible residue is represented as (chain, resname, resnum) in this module,
# resname could be empty string when user only provide residue number.
FlexResidue = Tuple[str, str, int]


class FlexDockingSetup(NamedTuple):
    """result of `prepare_flex_docking_setup`"""
    flex_residues: List[FlexResidue]
    grid_center: List[float]
    grid_size: List[float]
    rigid_pdbqt_path: str
    flex_pdbqt_path: str
    rigid_gpf_path: str


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


def select_flex_residues_by_distance(complex_pdbstr: str, ligand_chain: str,
                                     receptor_chain: List[str], distance: float) -> Optional[List[FlexResidue]]:
    """
    select receptor residues within `distance` (angstrom) of any ligand heavy atom.
    
    both sides are restricted to heavy atoms (`not elem H`), otherwise hydrogens
    would make the effective cutoff shorter than requested.
    
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list
    :param distance: distance cutoff in angstrom
    :return: list of (chain, resname, resnum), None if failed
    """
    receptor_chain = normalize_receptor_chain(complex_pdbstr, ligand_chain, receptor_chain)
    cmd.reinitialize()
    cmd.read_pdbstr(complex_pdbstr, 'complex')
    rec_sel = ' or '.join([f'chain {c}' for c in receptor_chain])
    cmd.select('lig_heavy', f'complex and chain {ligand_chain} and not elem H')
    if cmd.count_atoms('lig_heavy') == 0:
        return put_err(f'can not find any heavy atom on ligand chain {ligand_chain}')
    # `within` is a top-level binary operator and can NOT be mixed with `and`,
    # so both operands must be parenthesized selection expressions
    cmd.select('rec_near', f'(complex and ({rec_sel}) and not elem H) within {distance} of lig_heavy')
    result = set()
    cmd.iterate('rec_near', 'result.add((chain, resn, int(resi)))', space={'result': result})
    if not result:
        return put_err(f'can not find any receptor residue within {distance} A of ligand chain {ligand_chain}')
    return sorted(result, key=lambda x: (x[0], x[2]))


def select_flex_residues_by_interaction(complex_pdbstr: str, ligand_chain: str,
                                        receptor_chain: List[str], cutoff: float = 4.,
                                        mode: Union[str, List[str]] = 'all') -> Optional[List[FlexResidue]]:
    """
    select receptor residues which have interaction with ligand, using fplip (PLIP port).
    
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list
    :param cutoff: interaction distance cutoff in angstrom
    :param mode: interaction types to be calculated, 'all' to include all modes.
    :return: list of (chain, resname, resnum), None if failed
    """
    from lazydock.pml.plip_interaction import (check_support_mode,
                                               run_plip_analysis)
    receptor_chain = normalize_receptor_chain(complex_pdbstr, ligand_chain, receptor_chain)
    mode = check_support_mode(mode)
    interactions = run_plip_analysis(complex_pdbstr, receptor_chain, ligand_chain, mode, cutoff)
    if isinstance(interactions, Exception):
        return put_err(f'fplip analysis failed: {interactions}')
    if not interactions or not any(interactions.values()):
        return put_err(f'fplip found no interaction between receptor {receptor_chain} and ligand {ligand_chain}')
    # fplip keeps the residue numbering and resname of the input pdb (verified:
    # 0 mismatch on a real complex), so its values are used as-is.
    receptor_chain_set = set(receptor_chain)
    result = set()
    for _type, values in interactions.items():
        for rec_res, _lig_res, _dist in values:
            resi, resn, chain = rec_res
            if chain in receptor_chain_set:
                result.add((chain, resn, int(resi)))
    if not result:
        return put_err(f'fplip found no interaction residue on receptor chain {receptor_chain}')
    return sorted(result, key=lambda x: (x[0], x[2]))


def collect_flex_residues(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                          flex_residues: Union[str, List[str], List[FlexResidue]] = None,
                          flex_distance: float = None, flex_interaction: bool = False,
                          flex_interaction_mode: Union[str, List[str]] = 'all',
                          flex_interaction_cutoff: float = 4.) -> Optional[List[FlexResidue]]:
    """
    collect flexible residues from three sources, duplicated residues are merged.

        - flex_residues: explicit residue list, such as 'A:LYS68,A:LYS50' or [(chain, resname, resnum), ...]
        - flex_distance: select receptor residues within this distance of any ligand heavy atom
        - flex_interaction: select receptor residues which have interaction with ligand by fplip
    
    :return: list of (chain, resname, resnum); None if nothing is found or a source failed
    """
    result: List[FlexResidue] = []
    seen = set()

    def _add(flex):
        if not flex:
            return
        for c, r, n in flex:
            if (c, n) not in seen:
                seen.add((c, n))
                result.append((c, r, n))

    # 1. explicit
    if flex_residues:
        if isinstance(flex_residues, str) or (isinstance(flex_residues, list) and isinstance(flex_residues[0], str)):
            _add(parse_flex_residues(flex_residues))
        else:
            _add(flex_residues)
    # 2. by distance
    if flex_distance is not None:
        _add(select_flex_residues_by_distance(complex_pdbstr, ligand_chain, receptor_chain, flex_distance))
    # 3. by interaction
    if flex_interaction:
        _add(select_flex_residues_by_interaction(complex_pdbstr, ligand_chain, receptor_chain,
                                                 flex_interaction_cutoff, flex_interaction_mode))
    if not result:
        return None
    return result


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


def get_grid_with_flexres(lig_pdbstr: str, flex_residues: List[FlexResidue],
                          complex_pdbstr: str, grid_buffer: float):
    """
    compute the grid box covering the ligand AND every flexible residue, in angstrom.

    The rigid-receptor paths use `adfr.get_grid`, which wraps the ligand only; this
    variant exists because a flexible sidechain outside the box would have no grid to
    score against.

    The box is axis-aligned: center and size are the bounding box of
    (ligand + flexible residues) expanded by `grid_buffer` on every side.  A flexible
    residue that can not be found in `complex_pdbstr` is reported and skipped, so the
    box then covers the ligand only.

    :param lig_pdbstr: ligand pdb string
    :param flex_residues: flexible residues, list of (chain, resname, resnum)
    :param complex_pdbstr: complex pdb string, the residues are looked up here
    :param grid_buffer: box padding in ANGSTROM, added on every side
    :return: (grid_center, grid_size), both in ANGSTROM -- the unit that vina's
        `compute_vina_maps` and meeko's `--box_size` both document.  The one consumer
        that wants grid point counts instead divides by the spacing itself (see
        `adfr.prepare_gpf_and_map`); meeko's own gpf writer converts internally.
    """
    cmd.reinitialize()
    cmd.read_pdbstr(lig_pdbstr, 'ligand')
    ([minX, minY, minZ], [maxX, maxY, maxZ]) = cmd.get_extent('ligand')
    if flex_residues:
        cmd.read_pdbstr(complex_pdbstr, 'complex')
        sele = ' or '.join([f'(chain {c} and resi {n})' for c, _, n in flex_residues])
        if cmd.select('flexres', f'complex and ({sele})') > 0:
            ([fminX, fminY, fminZ], [fmaxX, fmaxY, fmaxZ]) = cmd.get_extent('flexres')
            minX, minY, minZ = min(minX, fminX), min(minY, fminY), min(minZ, fminZ)
            maxX, maxY, maxZ = max(maxX, fmaxX), max(maxY, fmaxY), max(maxZ, fmaxZ)
        else:
            put_err(f'can not find flexible residues {flex_residues} in complex, '
                    f'grid box will only cover ligand.')
    grid_center = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    grid_size = [maxX - minX + grid_buffer * 2, maxY - minY + grid_buffer * 2, maxZ - minZ + grid_buffer * 2]
    return grid_center, grid_size


def get_pdb_atom_keys(pdbstr: str) -> set:
    """return {(chain, resnum, atomname)} of all atoms in a pdb/pdbqt string"""
    result = set()
    for line in pdbstr.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            try:
                result.add((line[21], int(line[22:26]), line[12:16].strip()))
            except ValueError:
                continue
    return result


def build_flex_complex_pdbstr(receptor_pdbstr: str, rigid_pdbqtstr: str, pose_pdbstr: str) -> str:
    """
    build the receptor part of a flexible-docking complex pdb string.
    
    AutoDock only keeps the backbone of flexible residues in the rigid receptor, and
    moves CA + sidechain into the dlg pose. So we must drop, from the original receptor
    pdb, both the moved atoms AND the atoms bonded to them (e.g. HA on the moved CA),
    otherwise those hydrogens would dangle in space.
    
    the rule: for a flexible residue, only keep the atoms that meeko kept in the rigid
    pdbqt; a residue is flexible when the pose contains atoms of it that are not rigid.
    then the caller should append the dlg pose (ligand + flexible sidechains) after this string.
    
    :param receptor_pdbstr: original receptor pdb string (all atoms, full hydrogen)
    :param rigid_pdbqtstr: rigid receptor pdbqt string generated by meeko
    :param pose_pdbstr: one pose pdb string from the dlg file
    :return: receptor pdb string, flexible residues only keep their rigid atoms
    """
    rigid_keys = get_pdb_atom_keys(rigid_pdbqtstr)
    pose_keys = get_pdb_atom_keys(pose_pdbstr)
    # residues which have atoms moved out of the rigid receptor by AutoDock;
    # ligand residues are also included here, but they do not exist in receptor pdb.
    flex_residue_keys = {(c, n) for c, n, a in pose_keys if (c, n, a) not in rigid_keys}
    kept = []
    for line in receptor_pdbstr.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            try:
                key = (line[21], int(line[22:26]), line[12:16].strip())
            except ValueError:
                kept.append(line)
                continue
            if key[:2] in flex_residue_keys and key not in rigid_keys:
                continue
        kept.append(line)
    return '\n'.join(kept)


def build_flex_pose_complex_pdbstr(receptor_pdbstr: str, rigid_pdbqtstr: str,
                                   pose_pdbstr: str, ligand_chain: str,
                                   flex_residues: Optional[List[FlexResidue]] = None) -> str:
    """
    join a flex-docking receptor and one dlg pose into a single complex pdb string.

    `build_flex_complex_pdbstr` returns the receptor part with the moved atoms of
    every flexible residue removed, and its docstring tells the caller to append
    the pose.  A plain concatenation is NOT enough: with `--flexres`, AutoDock
    keeps a flexible residue's backbone (N, C, O, H) in the rigid receptor and
    moves CA + sidechain into the pose, so concatenating leaves the two halves of
    ONE residue in different parts of the file.  PDBFixer (used by the relaxation
    stage) reads residues as *contiguous atom blocks*: it then sees two residues
    where there is one, rebuilds a complete duplicate of each flexible residue,
    and the duplicates sit ~0.02 A apart.  Measured on a real flex run: 256 blocks
    for 245 residues, 194 duplicated atom keys, einit == efinal ~ 1e26 kJ/mol
    (the minimisation never ran) and every pose was rejected as failed.
    Reordering the atoms this way fixes it: 245 blocks, 0 duplicates,
    efinal ~ -2.7e4 kJ/mol.

    Every receptor residue is emitted as one contiguous block, that residue's
    pose-only atoms follow immediately, and the ligand is written last.

    :param receptor_pdbstr: original receptor pdb string (all atoms, full hydrogen)
    :param rigid_pdbqtstr: rigid receptor pdbqt string generated by meeko
    :param pose_pdbstr: one pose pdb string from the dlg file
    :param ligand_chain: ligand chain id (separates the ligand from the flex sidechains)
    :return: complex pdb string, ready for relaxation / PLIP / drift metrics
    """
    keep = build_flex_complex_pdbstr(receptor_pdbstr, rigid_pdbqtstr, pose_pdbstr)

    def _rows(text: str) -> List[str]:
        return [l for l in text.splitlines() if l.startswith(('ATOM', 'HETATM'))]

    def _res_key(line: str):
        return (line[21], int(line[22:26]))

    rec_rows = _rows(keep)
    pose_rows = _rows(pose_pdbstr)
    ligand_rows = [l for l in pose_rows if l[21] == ligand_chain]
    # atoms the pose adds on top of the receptor: the moved sidechains.  The
    # ligand is excluded here and written at the end.
    have = {(l[21], int(l[22:26]), l[12:16].strip()) for l in rec_rows}
    side_by_res: dict = {}
    for l in pose_rows:
        if l[21] == ligand_chain:
            continue
        k = (l[21], int(l[22:26]), l[12:16].strip())
        if k not in have:
            side_by_res.setdefault(_res_key(l), []).append(l)

    out: List[str] = []
    emitted = set()
    i = 0
    while i < len(rec_rows):
        rk = _res_key(rec_rows[i])
        while i < len(rec_rows) and _res_key(rec_rows[i]) == rk:
            out.append(rec_rows[i])
            i += 1
        if rk in side_by_res and rk not in emitted:
            out.extend(side_by_res[rk])
            emitted.add(rk)
    # flexible residues that have no backbone left in the receptor at all
    rec_keys = {_res_key(l) for l in rec_rows}
    for rk, ls in side_by_res.items():
        if rk not in emitted and rk not in rec_keys:
            out.extend(ls)
    if out:
        out.append('TER')
    out.extend(ligand_rows)
    out.append('END')
    return '\n'.join(out) + '\n'