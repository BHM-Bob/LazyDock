import os
import re
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file
from pymol import cmd

from lazydock.config import GlobalConfig
from lazydock.pml.utils import new_pml_context

# flexible residue is represented as (chain, resname, resnum) in this module,
# resname could be empty string when user only provide residue number.
FlexResidue = Tuple[str, str, int]


def get_pdbqtstr_from_pdbstr(pdbstr: str, caller: str = 'prepare_ligand', entry_param: str = '-l',
                             options: List[Tuple['str', str]] = None):
    """
    call ADFR programe to prepare ligand pdbqt file or receptor pdbqt file
    
    :param pdbstr: pdb string
    :param caller: ADFR programe name, could be 'prepare_ligand' or 'prepare_receptor', or absolute path
    :param options: options for ADFR programe, if None, will be [('-A', 'hydrogens')]
    :return: FLAG, pdbqt string, FLAG is True if success
    """
    options = options or [('-A', 'hydrogens')]
    cmd_options = []
    [cmd_options.extend(list(i)) for i in options]
    
    with tempfile.TemporaryDirectory() as tmpdir:
        opts_file(os.path.join(tmpdir, 'mol.pdb'), 'w', way='str', data=pdbstr)
        cmd = f"cd {tmpdir} && {caller} {entry_param} mol.pdb -o mol.pdbqt {' '.join(cmd_options)}"
        proc = subprocess.Popen(cmd, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = proc.communicate()
        output_str = stdout + stderr
        if proc.returncode != 0 or not Path(f'{tmpdir}/mol.pdbqt').exists():
            return False, output_str
        return True, opts_file(f'{tmpdir}/mol.pdbqt', way='str')


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


# first meeko release where --flexres takes one comma-separated sequence
MEEKO_FLAT_FLEXRES_VERSION = (0, 7, 0)


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


def _normalize_receptor_chain(complex_pdbstr: str, ligand_chain: str,
                              receptor_chain: List[str] = None) -> List[str]:
    """if receptor_chain is None, use all chains except ligand chain"""
    if receptor_chain:
        return receptor_chain
    cmd.reinitialize()
    cmd.read_pdbstr(complex_pdbstr, 'complex')
    return sorted(set(cmd.get_chains()) - {ligand_chain})


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
    receptor_chain = _normalize_receptor_chain(complex_pdbstr, ligand_chain, receptor_chain)
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
    receptor_chain = _normalize_receptor_chain(complex_pdbstr, ligand_chain, receptor_chain)
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


def _prepare_ligand_receptor(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str]):
    cmd.reinitialize()
    cmd.read_pdbstr(complex_pdbstr, 'complex')
    # receptor
    if not receptor_chain:
        receptor_chain = list(set(cmd.get_chains()) - set([ligand_chain]))
    if not receptor_chain:
        put_err(f'there is no other chain exclude ligand chain {ligand_chain}.')
        return None
    # set receptor
    if cmd.select('receptor', f'complex and (' + ' or '.join([f'chain {c}' for c in receptor_chain]) + ')') == 0:
        put_err(f'can not find receptor with chain {receptor_chain}')
        return None
    receptor_pdbstr = cmd.get_pdbstr('receptor')
    success, receptor_pdbqtstr = get_pdbqtstr_from_pdbstr(receptor_pdbstr, 'prepare_receptor', '-r')
    if not success:
        put_err(f"Failed to prepare receptor_pdbqtstr: {receptor_pdbqtstr}")
        return None
    # ligand
    if cmd.select('ligand', f'complex and chain {ligand_chain}') == 0:
        put_err(f'can not find ligand with chain {ligand_chain}')
        return None
    ligand_pdbstr = cmd.get_pdbstr(f'chain {ligand_chain}')
    success, ligand_pdbqtstr = get_pdbqtstr_from_pdbstr(ligand_pdbstr, 'prepare_ligand', '-l')
    if not success:
        put_err(f"Failed to prepare ligand_pdbqtstr: {ligand_pdbqtstr}")
        return None
    return receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def _get_grid(lig_pdbstr: str, grid_buffer: float):
    cmd.reinitialize()
    cmd.read_pdbstr(lig_pdbstr, 'ligand')
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent('ligand')
    grid_center = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2]
    grid_size = [maxX - minX + grid_buffer*2, maxY - minY + grid_buffer*2, maxZ - minZ + grid_buffer*2]
    return grid_center, grid_size

def _prepare_ligand_receptor_grid(tmpdir: str, v, complex_pdbstr: str,
                                  ligand_chain: str, receptor_chain: List[str], grid_buffer: float):
    """
    prepare ligand and receptor pdbqt file and grid
    
    return: v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr
    """
    pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain)
    if not pdbstr_pack:
        return None
    receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
    # set receptor
    rec_path = os.path.join(tmpdir, 'receptor.pdbqt')
    opts_file(rec_path, 'w', way='str', data=receptor_pdbqtstr)
    v.set_receptor(rec_path)
    # set ligand
    v.set_ligand_from_string(ligand_pdbqtstr)
    # set grid; not bounding box, just box in axis
    grid_center, grid_size = _get_grid(ligand_pdbstr, grid_buffer)
    v.compute_vina_maps(center=grid_center, box_size=grid_size, spacing=0.375)
    return v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def calc_vina_score(score_name: str, grid_buffer: float,
                    complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str], cpu: int = 1):
    """"total", "lig_inter", "flex_inter", "other_inter", "flex_intra", "lig_intra", "torsions", "-lig_intra"
    calculate vina score of a complex pdb file.
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box extend, box will calculate with pymol.cmd.extend without align to axis.
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param cpu: cpu number
    
    :return: score list, if None, means failed
    """
    from vina import Vina
    v = Vina(sf_name=score_name, cpu=cpu, verbosity=False)    
    with tempfile.TemporaryDirectory() as tmpdir:
        v = _prepare_ligand_receptor_grid(tmpdir, v, complex_pdbstr, ligand_chain, receptor_chain, grid_buffer)
        if not v:
            return None
        # compute and retrun score
        return v[0].score().tolist()


def perform_vina_dock(score_name: str, grid_buffer: float,
                      complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                      exhaustiveness: int, n_poses: int, min_rmsd: float = 1, max_evals: int = 0,
                      seed: int = None, energy_range: float = None, cpu: Optional[int] = None):
    """
    perform docking with vina
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: box extend, box will calculate with pymol.cmd.extend without align to axis.
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param exhaustiveness: exhaustiveness, equal to cpu number
    :param n_poses: n_poses
    :param min_rmsd: min_rmsd
    :param max_evals: max_evals
    :param seed: random seed for vina, None means vina decides (non-reproducible).
        Added for the redock pipeline, which promises reproducibility (its
        pipeline.json stores a seed); without it the same input gave -8.62 and
        -7.36 kcal/mol on two consecutive runs.  Kept optional so that every
        existing caller keeps its old behaviour.
    :param energy_range: maximum energy difference (kcal/mol) a pose may have
        from the best one to still be written.  None means vina's own default
        (3.0), which SILENTLY DROPS poses: requesting n_poses=20 returned only
        14 in a measured run, quietly breaking the pipeline's 160-pose-per-round
        contract.  Pass a large value (e.g. 100) to actually get `n_poses`.
    :param cpu: threads for vina, None keeps the historical behaviour of tying
        the thread count to `exhaustiveness`.  Those are different things and
        tying them bites hard: a caller asking for a cheap search
        (exhaustiveness=2) silently got a 2-thread dock, measured at ~10 min per
        task where 8 threads needed ~2 min.  The redock pipeline passes its own
        `vina_cpu` here (see `redock/loop._vina_dock`); every existing caller is
        unaffected because the default reproduces the old value.
    
    :return: True if success, False otherwise
    """
    from vina import Vina

    # historically the thread count was tied to `exhaustiveness`; keep that as
    # the default so existing callers behave exactly as before, but let the
    # caller override (the redock pipeline sizes its own thread budget)
    if cpu is None:
        cpu = exhaustiveness
    v = Vina(sf_name=score_name, cpu=cpu, seed=seed, verbosity=False)    
    with tempfile.TemporaryDirectory() as tmpdir:
        v = _prepare_ligand_receptor_grid(tmpdir, v, complex_pdbstr, ligand_chain, receptor_chain, grid_buffer)
        if not v:
            return False
        v, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = v
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses, min_rmsd=min_rmsd, max_evals=max_evals)
        output_path = os.path.join(tmpdir, 'poses.pdbqt')
        write_kwargs = {} if energy_range is None else {'energy_range': float(energy_range)}
        v.write_poses(output_path, n_poses, overwrite=True, **write_kwargs)
        result_str = opts_file(output_path, 'r', way='str')
        os.remove(output_path)
        return result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


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


def perform_vina_flex_dock(score_name: str, grid_buffer: float,
                           complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],
                           flex_residues: Union[str, List[str], List[FlexResidue]] = None,
                           flex_distance: float = None, flex_interaction: bool = False,
                           flex_interaction_mode: Union[str, List[str]] = 'all', flex_interaction_cutoff: float = 4.,
                           meeko_python: Optional[str] = None, gpf_caller: Optional[str] = None,
                           exhaustiveness: int = 32, n_poses: int = 20,
                           min_rmsd: float = 1, max_evals: int = 0,
                           seed: int = 0, cpu: int = 0, energy_range: float = None):
    """
    perform docking with vina, receptor flexible residues are supported.
    
    flexible residues could be provided in three ways, and they could be combined:
        - flex_residues: explicit residue list, such as 'A:LYS68,A:LYS50' or [(chain, resname, resnum), ...]
        - flex_distance: select receptor residues within this distance of any ligand heavy atom
        - flex_interaction: select receptor residues which have interaction with ligand by fplip
    
    The receptor is split into rigid / flexible parts by meeko, then handed to
    `Vina.set_receptor(rigid, flex)`.  Unlike the AutoDock-GPU backend no affinity
    maps files are involved: vina computes its own maps internally from the rigid part.
    
    NOTICE: vina performs the search on CPU and flexible sidechains add many degrees
    of freedom, so this is much slower than `perform_autodock_gpu_flex_dock`.
    
    :param score_name: vina score name, could be vina, vinardo, ad4
    :param grid_buffer: grid size buffer for grid box, will expand the grid box from ligand by buffer size
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param flex_residues: explicit flexible residues
    :param flex_distance: distance cutoff, select residues within this distance of ligand heavy atoms
    :param flex_interaction: FLAG, whether to select flexible residues by fplip interaction
    :param flex_interaction_mode: interaction types for fplip, 'all' to include all modes
    :param flex_interaction_cutoff: interaction distance cutoff for fplip
    :param meeko_python: python interpreter to call meeko, None means current interpreter
    :param gpf_caller: prepare_gpf.py path, only used to locate ADFRsuite bin (prepare_ligand/prepare_receptor)
    :param exhaustiveness: exhaustiveness, default is 32
    :param n_poses: n_poses
    :param min_rmsd: min_rmsd
    :param max_evals: max_evals, 0 means use vina heuristics
    :param seed: random seed
    :param cpu: cpu number for vina, 0 means use all cpus
    :param energy_range: max energy difference (kcal/mol) from the best pose for
        a pose to be written; None keeps vina's default (3.0), which silently
        drops poses (see `perform_vina_dock`)
    
    :return: (result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr),
        None if failed. `result_str` is a vina-format pdbqt string (NOT a dlg),
        in which every MODEL contains the ligand AND the flexible sidechains.
        A pose string can be turned into a complex by `build_flex_pose_complex_pdbstr`,
        or parsed by `lazydock.pml.autodock_utils.DlgFile` (it understands the vina
        `REMARK VINA RESULT` score line).
    """
    # parse flexible residues from three sources
    flex_residue_lst = collect_flex_residues(complex_pdbstr, ligand_chain, receptor_chain,
                                             flex_residues, flex_distance,
                                             flex_interaction, flex_interaction_mode,
                                             flex_interaction_cutoff)
    if not flex_residue_lst:
        return put_err('no flexible residue is found or provided.')
    gpf_caller = _resolve_gpf_caller(gpf_caller)
    if not gpf_caller:
        return put_err('gpf_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $PREPARE_GPF')
    if not Path(gpf_caller).exists():
        return put_err(f'gpf_caller does not exist: {gpf_caller}')
    _prepend_adfr_bin_to_path(gpf_caller)
    from vina import Vina
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor pdb string, and ligand pdbqt
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, _, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        ligand_path = os.path.join(tmpdir, 'ligand.pdbqt')
        opts_file(ligand_path, 'w', way='str', data=ligand_pdbqtstr)
        # grid box covers ligand and all flexible residues
        grid_center, grid_size = _get_grid_with_flexres(ligand_pdbstr, flex_residue_lst,
                                                        complex_pdbstr, grid_buffer)
        # meeko splits receptor into rigid part and flexible sidechains
        pdbqt_pack = prepare_receptor_with_flexres(tmpdir, receptor_pdbstr, flex_residue_lst,
                                                   grid_size, grid_center, meeko_python)
        if not pdbqt_pack:
            return None
        rigid_pdbqt_path, flex_pdbqt_path = pdbqt_pack
        # vina scores the flexible sidechains itself, no gpf/autogrid needed
        v = Vina(sf_name=score_name, cpu=cpu, seed=seed, verbosity=False)
        v.set_receptor(rigid_pdbqt_path, flex_pdbqt_path)
        v.set_ligand_from_file(ligand_path)
        v.compute_vina_maps(center=grid_center, box_size=grid_size, spacing=0.375)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses,
               min_rmsd=min_rmsd, max_evals=max_evals)
        output_path = os.path.join(tmpdir, 'poses.pdbqt')
        # energy_range default (3.0 kcal/mol) silently drops poses, see
        # `perform_vina_dock`; the flex path must honour n_poses as well
        v.write_poses(output_path, n_poses, overwrite=True,
                      **({} if energy_range is None else {'energy_range': float(energy_range)}))
        result_str = opts_file(output_path, 'r', way='str')
        rigid_pdbqtstr = opts_file(rigid_pdbqt_path, 'r', way='str')
        flex_pdbqtstr = opts_file(flex_pdbqt_path, 'r', way='str')
        return result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def prepare_gpf_and_map(docking_dir: str, caller: str, grid_center: list[float], grid_size: list[int]):
    grid_center = ','.join(map(lambda x: f'{x:.3f}', grid_center))
    grid_size = ','.join(map(lambda x: f'{x:.0f}', grid_size))
    os.system(f'cd {docking_dir} && pythonsh {caller} -r receptor.pdbqt -o receptor.gpf -p gridcenter={grid_center} -p npts={grid_size} -p ligand_types="A,C,OA,N,NA,SA,HD"')
    os.system(f'cd {docking_dir} && autogrid4 -p receptor.gpf -l receptor.glg')


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


def _get_grid_with_flexres(lig_pdbstr: str, flex_residues: List[FlexResidue],
                           complex_pdbstr: str, grid_buffer: float):
    """
    calculate grid box from ligand, then expand it to cover all flexible residues.
    
    the box is axis-aligned: center and size are the bounding box of (ligand + flexres)
    expanded by grid_buffer.
    
    :return: grid_center, grid_size
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


def perform_autodock_gpu_flex_dock_core(docking_dir: str, caller: str, ffile: str, lfile: str, flexfile: str,
                                        seed: int = 0, nrun: int = 5, nev: int = 500000,
                                        stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    """
    perform docking with AutoDock-GPU with flexible receptor sidechains
    
    :param docking_dir: docking dir, must be unique for parallel docking
    :param caller: caller, such as autodock_gpu_64wi_R128
    :param ffile: ffile, such as receptor_rigid.maps.fld
    :param lfile: lfile, such as ligand.pdbqt
    :param flexfile: flexible residue pdbqt file, such as receptor_flex.pdbqt
    :param seed: seed
    :param nrun: nrun
    :param nev: nev
    :param stopstd: stopstd
    :param heurmax: heurmax
    :param gpu_id: gpu id
    
    :return: dlg file path, None if failed
    """
    cmd_list = [
        'cd', docking_dir, '&&',
        caller,
        '--ffile', ffile,
        '--lfile', lfile,
        '--flexres', flexfile,
        '--seed', str(seed),
        '--resnam', 'dock',
        '--nrun', str(nrun),
        '--nev', str(nev),
        '--stopstd', str(stopstd),
        '--heurmax', str(heurmax),
        '--devnum', str(gpu_id + 1),
    ]
    proc = subprocess.Popen(' '.join(cmd_list), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = proc.communicate()
    output_str = stdout + stderr
    dlg_path = os.path.join(docking_dir, 'dock.dlg')
    if proc.returncode != 0 or not Path(dlg_path).exists():
        put_err(f'docking failed, output_str: {output_str}')
        return None
    return dlg_path


def perform_autodock_gpu_flex_dock(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str], grid_buffer: float,
                                   flex_residues: Union[str, List[str], List[FlexResidue]] = None,
                                   flex_distance: float = None, flex_interaction: bool = False,
                                   flex_interaction_mode: Union[str, List[str]] = 'all', flex_interaction_cutoff: float = 4.,
                                   autodock_caller: Optional[str] = None, gpf_caller: Optional[str] = None,
                                   meeko_python: Optional[str] = None,
                                   seed: int = 0, nrun: int = 5, nev: int = 500000,
                                   stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    """
    perform docking with AutoDock-GPU with flexible receptor sidechains.
    
    flexible residues could be provided in three ways, and they could be combined:
        - flex_residues: explicit residue list, such as 'A:LYS68,A:LYS50' or [(chain, resname, resnum), ...]
        - flex_distance: select receptor residues within this distance of any ligand heavy atom
        - flex_interaction: select receptor residues which have interaction with ligand by fplip
    
    :param complex_pdbstr: complex pdb string
    :param ligand_chain: ligand chain
    :param receptor_chain: receptor chain list, if None, will be all chains except ligand chain
    :param grid_buffer: grid size buffer for grid box, will expand the grid box from ligand by buffer size
    :param flex_residues: explicit flexible residues
    :param flex_distance: distance cutoff, select residues within this distance of ligand heavy atoms
    :param flex_interaction: FLAG, whether to select flexible residues by fplip interaction
    :param flex_interaction_mode: interaction types for fplip, 'all' to include all modes
    :param flex_interaction_cutoff: interaction distance cutoff for fplip
    :param autodock_caller: autodock_gpu executable path, if None, will be read from config or $AUTODOCK_GPU
    :param gpf_caller: prepare_gpf.py path, if None, will be read from config or $PREPARE_GPF
    :param meeko_python: python interpreter to call meeko, None means current interpreter
    :param seed, nrun, nev, stopstd, heurmax, gpu_id: AutoDock-GPU options
    
    :return: (result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr),
        None if failed. The 6th item is the ligand chain pdb string.
        Use `build_flex_complex_pdbstr` to build the complex pdb from the first 3 items
        and one pose pdb string from the dlg.
    """
    # parse flexible residues from three sources
    flex_residue_lst = collect_flex_residues(complex_pdbstr, ligand_chain, receptor_chain,
                                             flex_residues, flex_distance,
                                             flex_interaction, flex_interaction_mode,
                                             flex_interaction_cutoff)
    if not flex_residue_lst:
        return put_err('no flexible residue is found or provided.')
    # resolve ADFRsuite tools first: prepare_ligand/prepare_receptor/autogrid4
    gpf_caller = _resolve_gpf_caller(gpf_caller)
    if not gpf_caller:
        return put_err('gpf_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $PREPARE_GPF')
    if not Path(gpf_caller).exists():
        return put_err(f'gpf_caller does not exist: {gpf_caller}')
    _prepend_adfr_bin_to_path(gpf_caller)
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor pdb string, and ligand pdbqt
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, _, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        opts_file(os.path.join(tmpdir, 'ligand.pdbqt'), 'w', way='str', data=ligand_pdbqtstr)
        # prepare rigid receptor and flexible sidechain pdbqt by meeko
        # grid box is calculated from ligand and flexible residues
        grid_center, grid_size = _get_grid_with_flexres(ligand_pdbstr, flex_residue_lst,
                                                        complex_pdbstr, grid_buffer)
        pdbqt_pack = prepare_receptor_with_flexres(tmpdir, receptor_pdbstr, flex_residue_lst,
                                                   grid_size, grid_center, meeko_python)
        if not pdbqt_pack:
            return None
        rigid_pdbqt_path, flex_pdbqt_path = pdbqt_pack
        # maps are computed only from the rigid part, meeko already wrote the gpf
        # for rigid receptor, so do not use prepare_gpf here.
        rigid_gpf = os.path.join(tmpdir, 'receptor_rigid.gpf')
        if not Path(rigid_gpf).exists():
            return put_err(f'meeko did not write gpf file {rigid_gpf}, '
                           f'box center and size must be provided to meeko.')
        _run_autogrid(rigid_gpf, tmpdir, gpf_caller)
        # perform docking with flexible residues
        if not autodock_caller:
            autodock_caller = GlobalConfig.named_paths["autodock_gpu"] or os.environ.get("AUTODOCK_GPU")
        if not autodock_caller:
            return put_err('autodock_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $AUTODOCK_GPU')
        docking_path = perform_autodock_gpu_flex_dock_core(
            tmpdir, autodock_caller, 'receptor_rigid.maps.fld', 'ligand.pdbqt',
            os.path.basename(flex_pdbqt_path), seed, nrun, nev, stopstd, heurmax, gpu_id)
        if not docking_path:
            return None
        result_str = opts_file(docking_path, 'r', way='str')
        rigid_pdbqtstr = opts_file(rigid_pdbqt_path, 'r', way='str')
        flex_pdbqtstr = opts_file(flex_pdbqt_path, 'r', way='str')
        return result_str, receptor_pdbstr, rigid_pdbqtstr, flex_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


def _resolve_gpf_caller(gpf_caller: Optional[str] = None) -> Optional[str]:
    """resolve prepare_gpf.py path from argument, GlobalConfig or environment"""
    if not gpf_caller:
        gpf_caller = GlobalConfig.named_paths["prepare_gpf"] or os.environ.get("PREPARE_GPF")
    if not gpf_caller:
        return None
    return os.path.expanduser(gpf_caller)


def _prepend_adfr_bin_to_path(gpf_caller: str):
    """
    prepend ADFRsuite bin dir to $PATH, so that `prepare_ligand` and
    `prepare_receptor` could be found by the shared receptor preparation code.
    """
    bin_dir = os.path.dirname(os.path.abspath(gpf_caller))
    if not Path(bin_dir).is_dir():
        return
    path = os.environ.get('PATH', '')
    if bin_dir not in path.split(os.pathsep):
        os.environ['PATH'] = f'{bin_dir}{os.pathsep}{path}'


def _run_autogrid(gpf_path: str, workdir: str, gpf_caller: str):
    """run autogrid4 with the gpf file written by meeko"""
    # autogrid4 is expected to be in the same dir as prepare_gpf.py (ADFRsuite bin)
    autogrid_caller = 'autogrid4'
    caller_dir = os.path.dirname(os.path.abspath(os.path.expanduser(gpf_caller)))
    candidate = os.path.join(caller_dir, 'autogrid4')
    if Path(candidate).exists():
        autogrid_caller = candidate
    os.system(f'cd {workdir} && {autogrid_caller} -p {os.path.basename(gpf_path)} -l receptor_rigid.glg')


def perform_autodock_gpu_dock_core(docking_dir: str, caller: str, ffile: str, lfile: str,
                                   seed: int = 0, nrun: int = 5, nev: int = 500000,
                                   stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    """
    peform docking with AutoDock-GPU
    :param docking_dir: docking dir, must be unique for parallel docking
    :param caller: caller, such as autodock_gpu_64wi_R128
    :param ffile: ffile, such as receptor.maps.fld
    :param lfile: lfile, such as ligand.pdbqt
    :param seed: seed
    :param nrun: nrun
    :param nev: nev
    :param stopstd: stopstd
    :param heurmax: heurmax
    :param gpu_id: gpu id
    
    :return: True if success, False otherwise
    """
    cmd = [
        'cd', docking_dir, '&&',
        caller,
        '--ffile', ffile,
        '--lfile', lfile,
        '--seed', str(seed),
        '--resnam', f'dock',
        '--nrun', str(nrun),
        '--nev', str(nev),
        '--stopstd', str(stopstd),
        '--heurmax', str(heurmax),
        '--devnum', str(gpu_id+1),
    ]
    proc = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = proc.communicate()
    output_str = stdout + stderr
    # return path
    if proc.returncode != 0 or not (docking_dir / Path(f'dock.dlg')).exists():
        print(f'Error: docking failed, output_str: {output_str}')
    return os.path.join(docking_dir, 'dock.dlg')


def perform_autodock_gpu_dock(complex_pdbstr: str, ligand_chain: str, receptor_chain: List[str],grid_buffer: float,
                              autodock_caller: Optional[str] = None, gpf_caller: Optional[str] = None,
                              seed: int = 0, nrun: int = 5, nev: int = 500000,
                              stopstd: float = 0.5, heurmax: int = 3000000, gpu_id: int = 0):
    with tempfile.TemporaryDirectory() as tmpdir:
        # prepare ligand and receptor
        pdbstr_pack = _prepare_ligand_receptor(complex_pdbstr, ligand_chain, receptor_chain)
        if not pdbstr_pack:
            return None
        receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr = pdbstr_pack
        opts_file(os.path.join(tmpdir, 'receptor.pdbqt'), 'w', way='str', data=receptor_pdbqtstr)
        opts_file(os.path.join(tmpdir, 'ligand.pdbqt'), 'w', way='str', data=ligand_pdbqtstr)        
        # get grid center and size
        grid_center, grid_size = _get_grid(ligand_pdbstr, grid_buffer)
        # setup docking files
        if not gpf_caller:
            gpf_caller = GlobalConfig.named_paths["prepare_gpf"] or os.environ.get("PREPARE_GPF")
        if not gpf_caller:
            return put_err('gpf_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $PREPARE_GPF')
        prepare_gpf_and_map(tmpdir, gpf_caller, grid_center, grid_size)
        # perform docking
        if not autodock_caller:
            autodock_caller = GlobalConfig.named_paths["autodock_gpu"] or os.environ.get("AUTODOCK_GPU")
        if not autodock_caller:
            return put_err('autodock_caller is None and can not be found with ~/.lazydock/lazydock_config.json or $AUTODOCK_GPU')
        docking_path = perform_autodock_gpu_dock_core(tmpdir, autodock_caller, 'receptor.maps.fld', 'ligand.pdbqt',
                                                      seed, nrun, nev, stopstd, heurmax, gpu_id)
        if not docking_path:
            return None
        # return docking path
        result_str = opts_file(docking_path, 'r', way='str')
        return result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr


if __name__ == '__main__':
    pdbqtstr = get_pdbqtstr_from_pdbstr("")
    complex_pdbstr = opts_file('data_tmp/pdb/complex_cycpep.pdb')
    result_str, receptor_pdbstr, receptor_pdbqtstr, ligand_pdbstr, ligand_pdbqtstr \
        = perform_autodock_gpu_dock(complex_pdbstr, 'P', ['A'], 4)
