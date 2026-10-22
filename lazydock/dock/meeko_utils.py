'''
Date: 2026-10-09
LastEditors: BHM-Bob
Description: meeko based receptor/ligand preparation utilities.

Every function here drives the meeko PYTHON API in process.  meeko is a standard
dependency of lazydock, so there is no subprocess route and no `python_caller` knob
anywhere: a molecule prepared by another interpreter can simply be passed in as a
pdb string.

This module must NOT import `lazydock.dock.adfr`: that module imports this one, so
the reverse would be a cycle.  Everything layer-related is documented in
`docs/dev/dock/dock_module_architecture.md`.

Contents:
    - rigid receptor preparation (`prepare_receptor_with_meeko`)
    - flexible receptor preparation: the rigid/flex split plus the rigid gpf
      (`prepare_receptor_with_flexres`)
    - ligand preparation, optionally breaking a macrocycle (`prepare_ligand_with_meeko`)
    - macrocycle / glue bookkeeping (`inspect_macrocycle`, `strip_glue_atoms`)
    - the shared preamble of the flexible-receptor docking paths
      (`prepare_flex_docking_setup`)
    - a meeko CLI `--flexres` string formatter (`get_meeko_version`, used only by
      `flex_utils.flex_residues_to_meeko`; nothing on the docking path needs it, since
      the API above takes residue ids directly)
'''
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Tuple, Union

from mbapy_lite.base import put_err
from mbapy_lite.file import opts_file
from meeko import MoleculePreparation, PDBQTWriterLegacy, gridbox
from meeko.polymer import Polymer, ResidueChemTemplates

from lazydock.dock.config import (GLUE_ATOM_NAMES, MACROCYCLE_MAX_RING_SIZE,
                                  RIGID_MACROCYCLES, SS_RIGIDIFY_INDICES,
                                  SS_RIGIDIFY_SMARTS)
from lazydock.dock.flex_utils import (FlexDockingSetup, FlexResidue,
                                      collect_flex_residues,
                                      get_grid_with_flexres)

# the two atom-type lists are copied verbatim from meeko's own CLI
# (`meeko/cli/mk_prepare_receptor.py`, the `rec_types` and `any_lig_base_types`
# literals) because `get_gpf_string` needs them and meeko exposes no constant for
# them.  They define which affinity maps autogrid4 computes, so a drift from meeko
# would silently change what the GPU engine can score.
_GPF_RECEPTOR_TYPES = ['HD', 'C', 'A', 'N', 'NA', 'OA', 'F', 'P', 'SA', 'S',
                       'Cl', 'Br', 'I', 'Mg', 'Ca', 'Mn', 'Fe', 'Zn']
_GPF_LIGAND_TYPES = ['HD', 'C', 'A', 'N', 'NA', 'OA', 'F', 'P', 'SA', 'S',
                     'Cl', 'Br', 'I', 'Si', 'B']


def get_meeko_version(caller: Optional[str] = None) -> Optional[Tuple[int, ...]]:
    """
    detect the installed meeko version, using the given python interpreter.

    ONLY needed to format a meeko CLI `--flexres` string, whose syntax changed
    incompatibly at 0.7: <=0.6 grouped residues per chain and joined the groups with
    ';' (e.g. 'A:78,79;T:22,24'), while >=0.7 expects ONE comma-separated sequence
    with the chain prefix repeated for every group (e.g. 'A:78,79,T:22,24').  See
    `flex_utils.flex_residues_to_meeko` for details.

    The Python API path used by this module does NOT go through that string at all:
    `prepare_receptor_with_meeko` hands meeko the residue ids as a set, so nothing
    under `lazydock/dock/` calls this function.  It is kept because
    `flex_utils.flex_residues_to_meeko` is a public formatter that other code (and the
    older `lazydock.pml.meeko_utils`) still uses together with a meeko CLI.

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


# ============================
# === receptor preparation ===
# ============================


def prepare_receptor_with_meeko(receptor_pdbstr: str, flex_residues: List[FlexResidue] = None,
                                allow_bad_res: bool = True):
    """
    prepare a receptor pdbqt string with meeko, without touching the filesystem.

    meeko must be importable by the CURRENT interpreter (it is a standard dependency,
    so this always holds).  The output is byte-identical to what
    `python -m meeko.cli.mk_prepare_receptor` writes: measured 97440 bytes for the
    rigid-only case here, and 96720 + 993 for the rigid and flex parts of a
    one-residue case.

    NOTICE: meeko's atom typing is NOT identical to ADT's -- on a real receptor the
    heavy atoms match one-by-one but 2 nitrogen acceptors are classified differently
    (meeko: NA=4, N=164 / ADT: NA=2, N=166) and 4 polar hydrogens are typed
    differently (meeko HD=221 / ADT HD=225).  vina's H-bond and N/NA terms depend on
    exactly these types, so scores from the two backends are NOT comparable and
    existing ADT-prepared results can not be reproduced with meeko.

    :param receptor_pdbstr: receptor pdb string, hydrogens are added by meeko
    :param flex_residues: flexible residues as list of (chain, resname, resnum), None
        or empty means a fully rigid receptor.  Only the chain and the residue number
        are used (meeko keys its monomers as '<chain>:<resnum>'); the resname is
        accepted for symmetry with the rest of this module but ignored.
    :param allow_bad_res: what to do with a residue meeko can not match to a template.
        True (the default) means IGNORE it: the residue is dropped from the polymer and
        a warning is logged, so one unusual residue does not abort a whole batch.
        False means raise.  The name is misleading -- True causes REMOVAL, not
        retention -- and meeko's CLI exposes the same switch as `-a/--allow_bad_res`.
        For a flexible residue, note that `flex_residues` below is checked separately:
        a residue that was dropped for being bad can not be made flexible either.
    :return: (pdbqt_string, flex_pdbqt_string) where the second item is the flexible
        part (empty string when `flex_residues` is None/empty), None if failed.
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
    # non-empty ones are real (meeko also lists rotatable terminal groups here).
    # Concatenated with NO separator, which is what meeko's own CLI does
    # (`all_flex_pdbqt += flexres_pdbqt`): every value already ends with a newline, so
    # joining with '\n' would insert a blank line between the residues and produce a
    # measurably different file (measured on a 4-residue probe: 3251 bytes joined with
    # '' against 3254 with '\n', the extra 3 being the inserted blank lines).
    flex_pdbqt = ''.join(v for v in flex_dict.values() if v)
    return rigid_pdbqt, flex_pdbqt


def prepare_receptor_with_flexres(tmpdir: str, receptor_pdbstr: str, flex_residues: List[FlexResidue],
                                  box_size: List[float] = None, box_center: List[float] = None,
                                  gpf_name: str = 'receptor_rigid.gpf'):
    """
    prepare a receptor pdbqt pair (rigid part + flexible sidechains) and its gpf.

    This is a filesystem wrapper around `prepare_receptor_with_meeko`: it writes what
    the AutoDock-GPU flexible path consumes, because that path hands the files to
    autogrid4 and to the GPU executable by name.

        - `<tmpdir>/receptor_rigid.pdbqt`   rigid part, the one the maps are built from
        - `<tmpdir>/receptor_flex.pdbqt`    the flexible sidechains
        - `<tmpdir>/<gpf_name>`             the gpf for the rigid part, fed to autogrid4
        - `<tmpdir>/boron-silicon-atom_par.dat`  the atom-parameter file the gpf references

    Measured byte-level equivalence with `python -m meeko.cli.mk_prepare_receptor` on
    the same box: rigid 96720, flex 993, gpf 749 and `.dat` 142 bytes, all identical.

    meeko adds the missing polar hydrogens itself, so a hydrogen-free pdb string is
    fine.

    NOTICE on the grid point count: `meeko.gridbox.get_gpf_string` rounds each axis
    down to an EVEN number of points (`2 * int(size / (2 * spacing))`), while
    `adfr.prepare_gpf_and_map` uses a plain floor (`int(size / spacing)`).  When a box
    size is not an exact multiple of the spacing the two grid writers therefore
    disagree by one point per axis.  That is a property of meeko's writer, not a bug
    here, and it is why this path must not reuse that helper.

    :param tmpdir: working dir
    :param receptor_pdbstr: receptor pdb string
    :param flex_residues: flexible residues, list of (chain, resname, resnum)
    :param box_size: grid box size (x, y, z) in ANGSTROM.  The gpf is written only
        when BOTH box arguments are given; meeko converts angstrom to npts itself
        with its own 0.375 spacing
    :param box_center: grid box center (x, y, z) in ANGSTROM
    :param gpf_name: file name for the gpf, inside `tmpdir`.  The GPU path feeds the
        file to autogrid4, which derives the map names from the gpf's `gridfld` line,
        so renaming it also renames the maps it will write
    :return: (rigid_pdbqt_path, flex_pdbqt_path, gpf_path), None if failed.
        `gpf_path` is None when no box was given, because no gpf is written then; the
        caller runs autogrid4 on it (see `adfr.run_autogrid4`)
    """
    prepared = prepare_receptor_with_meeko(receptor_pdbstr, flex_residues)
    if not prepared:
        return None
    rigid_pdbqtstr, flex_pdbqtstr = prepared
    if not flex_pdbqtstr:
        return put_err(f'meeko produced no flexible sidechain pdbqt for '
                       f'{flex_residues}, the receptor can not be split')
    rigid_path = os.path.join(tmpdir, 'receptor_rigid.pdbqt')
    flex_path = os.path.join(tmpdir, 'receptor_flex.pdbqt')
    opts_file(rigid_path, 'w', way='str', data=rigid_pdbqtstr)
    opts_file(flex_path, 'w', way='str', data=flex_pdbqtstr)
    if not (box_size and box_center):
        return rigid_path, flex_path, None
    # the gpf needs the atom-parameter file to sit next to it, because autogrid4
    # resolves `parameter_file` relative to its own working dir
    gpf_path = os.path.join(tmpdir, gpf_name)
    gpf_str, _ = gridbox.get_gpf_string(
        box_center, box_size, os.path.basename(rigid_path),
        _GPF_RECEPTOR_TYPES, _GPF_LIGAND_TYPES,
        ff_param_fname='boron-silicon-atom_par.dat')
    opts_file(os.path.join(tmpdir, 'boron-silicon-atom_par.dat'), 'w',
              way='str', data=gridbox.boron_silicon_atompar)
    opts_file(gpf_path, 'w', way='str', data=gpf_str)
    return rigid_path, flex_path, gpf_path


def prepare_flex_docking_setup(tmpdir: str, complex_pdbstr: str, ligand_chain: str,
                               receptor_chain: List[str], receptor_pdbstr: str,
                               ligand_pdbstr: str, grid_buffer: float,
                               flex_residues: Union[str, List[str], List[FlexResidue]] = None,
                               flex_distance: float = None, flex_interaction: bool = False,
                               flex_interaction_mode: Union[str, List[str]] = 'all',
                               flex_interaction_cutoff: float = 4.):
    """
    shared preamble of the two flexible-receptor docking paths (AutoDock-GPU and vina).

    Both engines run the same three steps in the same order, and each used to carry its
    own copy of them:

      1. collect the flexible residues from the three sources (explicit / by distance /
         by interaction)
      2. compute the grid box covering the ligand AND every flexible residue
      3. split the receptor into a rigid pdbqt and a flexible-sidechain pdbqt with
         meeko, which also writes the rigid gpf.  The GPU path feeds that gpf to
         autogrid4; vina builds its own maps from the rigid pdbqt and ignores the gpf

    ADFRsuite is deliberately NOT required here.  It matters only to
    `adfr.get_pdbqtstr_from_pdbstr` (for `prepare_ligand` / `prepare_receptor`) and to
    `adfr.run_autogrid4` / `adfr.prepare_gpf_and_map`, and each of those resolves the
    bin dir itself, so a meeko-only run never needs ADFRsuite configured.

    The caller runs `adfr._prepare_ligand_receptor` itself and passes the result in,
    because the two engines hand it different arguments (only the vina path takes
    `macrocycle`).

    :param tmpdir: working dir, shared with the caller
    :param complex_pdbstr: complex pdb string, used to select the flexible residues
    :param ligand_chain: ligand chain, used to select the flexible residues
    :param receptor_chain: receptor chain list, used to select the flexible residues
    :param receptor_pdbstr: receptor pdb string, split by meeko into rigid + flexible
    :param ligand_pdbstr: ligand pdb string, used to compute the grid box
    :param grid_buffer: grid padding in ANGSTROM (see `flex_utils.get_grid_with_flexres`);
        `grid_center`/`grid_size` come back in angstrom as well
    :param flex_residues / flex_distance / flex_interaction / flex_interaction_mode /
        flex_interaction_cutoff: the three sources of flexible residues, merged and
        deduplicated by `flex_utils.collect_flex_residues`
    :return: FlexDockingSetup (see `flex_utils.FlexDockingSetup`), None if failed
        (each step reports its own error)
    """
    flex_residue_lst = collect_flex_residues(complex_pdbstr, ligand_chain, receptor_chain,
                                             flex_residues, flex_distance,
                                             flex_interaction, flex_interaction_mode,
                                             flex_interaction_cutoff)
    if not flex_residue_lst:
        return put_err('no flexible residue is found or provided.')
    grid_center, grid_size = get_grid_with_flexres(ligand_pdbstr, flex_residue_lst,
                                                   complex_pdbstr, grid_buffer)
    # meeko writes the rigid gpf itself, from the same box; it takes angstrom
    pdbqt_pack = prepare_receptor_with_flexres(tmpdir, receptor_pdbstr, flex_residue_lst,
                                               grid_size, grid_center)
    if not pdbqt_pack:
        return None
    rigid_pdbqt_path, flex_pdbqt_path, rigid_gpf_path = pdbqt_pack
    if not rigid_gpf_path:
        return put_err('meeko did not write a rigid gpf; the grid box must be '
                       'provided for the GPU path to build its maps')
    return FlexDockingSetup(flex_residue_lst, grid_center, grid_size,
                            rigid_pdbqt_path, flex_pdbqt_path, rigid_gpf_path)


# ==========================
# === ligand preparation ===
# ==========================


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

    meeko must be importable by the CURRENT interpreter.  There is no `python_caller`
    knob anywhere in this module: every path uses the meeko API.

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