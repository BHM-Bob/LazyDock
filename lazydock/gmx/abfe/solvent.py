'''
Date: 2026-08-28
LastEditors: BHM-Bob 2262029386@qq.com
Description: solvation tools for ABFE, migrated from BindFlow bindflow.preparation.solvent.
             The gmx command layer is replaced with LazyDock's Gromacs class.
             Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
import logging
import os
import shutil
import tempfile
from pathlib import Path
from typing import Iterable, List, Tuple, Union

import yaml
from parmed import Structure

from lazydock.gmx.abfe._home import gmx_water_models_data_dir
from lazydock.gmx.run import Gromacs

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]


def get_atom_types(top: PathLike) -> dict:
    """Return the atomtypes section as a dict with key atom type and value the
    corresponding line. Include statements are not taken into account."""
    atom_types = {}
    with open(top, 'r') as f:
        lines = f.readlines()
    section_found = False
    for line in lines:
        if line.startswith('[ atomtypes ]'):
            section_found = True
            continue
        if section_found:
            if line.startswith(';'):
                continue
            elif (not line.strip() or line.startswith('[')) and not line.startswith('[ atomtypes ]'):
                section_found = False
                continue
            fields = line.split()
            if len(fields) >= 6:
                atom_type = fields[0]
                atom_types[atom_type] = line
    return atom_types


def get_molecule_names(input_topology: PathLike, section: str = 'molecules') -> list:
    """Get the molecule names specified in input_topology [molecules] or [moleculetype]."""
    if section not in ['molecules', 'moleculetype']:
        raise ValueError(f"section must be 'molecules', 'moleculetype'. {section} was provided.")
    with open(input_topology, 'r') as f:
        lines = f.readlines()
    molecules = []
    i = 0
    while i < len(lines):
        if section in lines[i]:
            i += 1
            while ("[" not in lines[i]):
                if not lines[i].startswith(';'):
                    split_line = lines[i].split()
                    if len(split_line) == 2:
                        molecules.append(split_line[0])
                i += 1
                if i >= len(lines):
                    break
        i += 1
    return molecules


def add_posres_section(input_topology: PathLike, molecules: Iterable[str], out_file: PathLike = None):
    """Add to the topology file the POSRES sections for the provided molecules."""
    with open(input_topology, "r") as f:
        top_lines = f.readlines()
    add_sol_internally = False
    if 'SOL' not in molecules:
        molecules = list(molecules) + ['SOL']
        add_sol_internally = True
    look_out_flag = False
    out_lines = []
    for line in top_lines:
        if not line.startswith("[ molecules ]"):
            for molecule in molecules:
                if molecule in line and " 3\n" in line:
                    look_out_flag = True
                    mol_name = line.split()[0]
                if look_out_flag and ('[ moleculetype ]' in line or '[ system ]' in line):
                    if mol_name == 'SOL' and add_sol_internally:
                        out_lines.append("\n#ifdef POSRES_WATER\n")
                        out_lines.append("; Position restraint for each water oxygen\n")
                        out_lines.append("[ position_restraints ]\n")
                        out_lines.append(";  i funct       fcx        fcy        fcz\n")
                        out_lines.append("   1    1       1000       1000       1000\n")
                        out_lines.append("#endif\n\n")
                    else:
                        out_lines.append("\n#ifdef POSRES\n")
                        out_lines.append(f'#include "posres_{mol_name}.itp"\n')
                        out_lines.append("#endif\n\n")
                    look_out_flag = False
        out_lines.append(line)
    if not out_file:
        out_file = input_topology
    with open(out_file, "w") as w:
        w.write("".join(out_lines))


def make_posres(input_topology: PathLike, molecules: Iterable[str], out_dir: PathLike,
                f_xyz: tuple = (2500, 2500, 2500), all_atoms: bool = False):
    """Make position restraint files out of input_topology for all the molecules
    specified on molecules.

    By default only heavy atoms (mass > 3) are restrained (BindFlow behaviour).
    With ``all_atoms=True`` every atom including hydrogens is restrained - this
    is needed for the equilibration of a fully-decoupled (lambda=0 vdw window)
    flexible/cyclic-peptide ligand: free H atoms drift far away from their
    restrained heavy atoms over ~10 ps (measured HA displacement 0.16 nm vs
    CA 0.06 nm at the 02_npt crash frame), eventually breaking the LINCS
    constraints.
    """
    for molecule in molecules:
        atom_flag = False
        with open(input_topology, "r") as f:
            top_lines = f.readlines()
        posres_filename = f"posres_{molecule}.itp"
        with open(Path(out_dir) / posres_filename, "w") as posres_file:
            posres_file.write("[ position_restraints ]\n")
            for i, _ in enumerate(top_lines):
                if f"{molecule}  " in (top_lines[i]) and " 3\n" in (top_lines[i]):
                    j = i + 1
                    while j < len(top_lines):
                        if '[ atoms ]' in top_lines[j]:
                            j += 1
                            atom_flag = True
                        if top_lines[j].startswith('['):
                            break
                        if atom_flag:
                            if not top_lines[j].startswith("\n") and not top_lines[j].startswith(";") and not top_lines[j].startswith("#"):
                                atoms_cols = top_lines[j].split()
                                if all_atoms or float(atoms_cols[7]) > 3:
                                    posres_str = f"{atoms_cols[0]} 1 {f_xyz[0]} {f_xyz[1]} {f_xyz[2]}\n"
                                    posres_file.write(posres_str)
                        j += 1
                    break
    add_posres_section(input_topology=input_topology, molecules=molecules, out_file=None)


def _tip3p_settles_to_constraints(top: PathLike, molecule: str, out_top: Union[PathLike, None] = None) -> None:
    """Change the settles section of `molecule` to tip3p-like triangular constraints."""
    constraints_section = "; https://gromacs.bioexcel.eu/t/how-to-treat-specific-water-molecules-as-ligand/3470/9\n"\
        "[ constraints ]\n"\
        "; ai aj funct length\n"\
        "1 2 1 0.09572\n"\
        "1 3 1 0.09572\n"\
        "2 3 1 0.15139\n\n"
    with open(top, 'r') as f:
        lines = f.readlines()
    idx_begins, idx_ends = None, None
    section_found = False
    i = 0
    while not lines[i].startswith('[ molecules ]') and i < len(lines):
        if molecule in lines[i] and " 3\n" in lines[i]:
            j = i
            while not lines[j].startswith('[ moleculetype ]') and j < len(lines):
                if lines[j].startswith('[ settles ]'):
                    section_found = True
                    idx_begins = j
                    j += 1
                if section_found and lines[j].startswith(('[', '#')):
                    idx_ends = j
                    break
                j += 1
            break
        i += 1
    if not out_top:
        out_top = top
    with open(out_top, 'w') as f:
        f.write("".join(lines[:idx_begins]) + constraints_section + "".join(lines[idx_ends:]))


class Solvate:
    """Solvate GMX systems. Migrated from BindFlow, using LazyDock Gromacs class for gmx calls.

    Available water models (data under abfe/data/gmx_water_models):
        amber: spc, spce, tip3p, tip4p, tip4pew, tip5p
        charmm: spc, spce, tip3p, tips3p, tip4p, tip5p
        oplsaa: spc, spce, tip3p, tip4p, tip4pew, tip5p, tip5pe
    """

    def __init__(self, water_model_code: str, builder_dir: PathLike = '.solvate',
                 load_dependencies: List[str] = None, gmx: Gromacs = None, cwd: str = None,
                 maxwarn: Union[int, None] = None):
        self.load_dependencies = load_dependencies
        self.gmx = gmx or Gromacs(working_dir=str(Path(builder_dir).resolve()))
        self.maxwarn = maxwarn
        self.builder_dir = Path(builder_dir).resolve()
        self.builder_dir.mkdir(exist_ok=True, parents=True)
        self.solvated_dir = self.builder_dir / 'solvated_sys'

        with open(Path(gmx_water_models_data_dir()) / 'water_models.yml', 'r') as f:
            self.water_models_data = yaml.safe_load(f)

        force_field_family, water_model = water_model_code.split('/')
        if force_field_family not in self.water_models_data:
            raise ValueError(f"Invalid force field family: {force_field_family}. Choose from {self.water_models_data.keys()}")
        elif water_model not in self.water_models_data[force_field_family]:
            raise ValueError(f"Invalid water model ({water_model}) for {force_field_family}."
                             f"Choose from {self.water_models_data[force_field_family].keys()}")

        self.force_field_family = force_field_family
        self.water_model = water_model
        self.water_itp, self.ions_itp, self.ffnonbonded_itp, self.water_gro = self._get_gmx_water_model()
        self.cwd = cwd or os.getcwd()

    def _get_gmx_water_model(self) -> Tuple[PathLike]:
        ff_dir = Path(gmx_water_models_data_dir())
        water_itp = (ff_dir / str(self.force_field_family) / f"{self.water_model}.itp").resolve()
        ions_itp = (ff_dir / str(self.force_field_family) / "ions.itp").resolve()
        ffnonbonded_itp = (ff_dir / str(self.force_field_family) / "ffnonbonded.itp").resolve()
        water_gro = (ff_dir / "configurations" / self.water_models_data[self.force_field_family][self.water_model]).resolve()
        return water_itp, ions_itp, ffnonbonded_itp, water_gro

    def _include_all_atom_types(self, top: PathLike) -> None:
        """Add all the atom types of the force field family to the first [ atomtypes ] section."""
        with open(top, 'r') as f:
            lines = f.readlines()
        idx_begins, idx_ends = None, None
        section_found = False
        for i, line in enumerate(lines):
            if line.startswith('[ atomtypes ]'):
                idx_begins = i + 1
                section_found = True
                continue
            if section_found:
                if line.startswith('['):
                    idx_ends = i
                    break
        if idx_begins is not None and idx_ends is not None:
            atom_types = get_atom_types(top)
            for atom_type_name, atom_type_info in get_atom_types(self.ffnonbonded_itp).items():
                if atom_type_name not in atom_types:
                    atom_types[atom_type_name] = atom_type_info
            with open(top, 'w') as f:
                f.write("".join(lines[:idx_begins] + list(atom_types.values()) + ["\n\n"] + lines[idx_ends:]))

    def _include_water_ions_params(self, top: PathLike) -> None:
        """Add include statements for water and ions itp files."""
        include_statements = [
            f"#include \"{self.water_itp}\"\n",
            f"#include \"{self.ions_itp}\"\n",
        ]
        with open(top, 'r') as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith('[ system ]'):
                break
        with open(top, 'w') as f:
            f.write("".join(lines[:idx] + include_statements + lines[idx:]))

    def _add_water_and_ions(self, gro: PathLike, top: PathLike, bt: str = "triclinic",
                            box: list = None, angles: list = None, d: float = None,
                            c: bool = False, pname: str = "NA", nname: str = "CL",
                            ion_conc: float = 150E-3, rmin: float = 1.0) -> None:
        """Make box, solvate and add ions to the system using LazyDock Gromacs class."""
        os.chdir(self.solvated_dir)
        self.gmx = Gromacs(working_dir=str(self.solvated_dir))

        editconf_kwargs = dict(f=gro, o=gro, bt=bt)
        if box:
            editconf_kwargs['box'] = ' '.join([str(i) for i in box])
        if angles:
            editconf_kwargs['angles'] = ' '.join([str(i) for i in angles])
        if d:
            editconf_kwargs['d'] = d
        if c:
            editconf_kwargs['c'] = True

        # First write an ions.mdp file
        with open("ions.mdp", "w") as file:
            file.write("; Neighbor searching\n"
                       "cutoff-scheme           = Verlet\n"
                       "rlist                   = 1.1\n"
                       "pbc                     = xyz\n"
                       "verlet-buffer-tolerance = -1\n"
                       "\n; Electrostatics\n"
                       "coulombtype             = cut-off\n"
                       "\n; VdW\n"
                       "rvdw                    = 1.0\n")

        # Execute the GMX functions via LazyDock's Gromacs class
        grompp_kwargs = dict(f="ions.mdp", c=gro, p=top, o="ions.tpr")
        if self.maxwarn is not None:
            grompp_kwargs['maxwarn'] = self.maxwarn
        self.gmx.run_gmx_with_expect('editconf', **editconf_kwargs)
        self.gmx.run_gmx_with_expect('solvate', cp=gro, p=top, cs=self.water_gro, o=gro)
        self.gmx.run_gmx_with_expect('grompp', **grompp_kwargs)
        # genion needs SOL group selection via expect
        self.gmx.run_gmx_with_expect('genion', s="ions.tpr", p=top, o=gro, neutral=True,
                                     pname=pname, nname=nname, rmin=rmin, conc=ion_conc,
                                     expect_actions=[{'Select a group:': 'SOL\r'}])

        # Just to clean the topology: build a monolithic topology
        struc = _read_parmed_molecule(top_file=top, gro_file=gro)
        struc.save(str(top), overwrite=True)
        struc.save(str(gro), overwrite=True)

        os.chdir(self.cwd)

    def clean(self, directory: Union[None, PathLike] = None) -> None:
        os.chdir(self.cwd)
        dir2delete = directory if directory else self.builder_dir
        try:
            shutil.rmtree(dir2delete)
        except FileNotFoundError:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.clean()

    def __call__(self, structure: Structure, bt: str = "triclinic", box: list = None,
                 angles: list = None, d: float = None, c: bool = False, pname: str = "NA",
                 nname: str = "CL", ion_conc: float = 150E-3, rmin: float = 1.0,
                 exclusion_list: list = None, out_dir: PathLike = '.', out_name: str = 'solvated',
                 f_xyz: tuple = (2500, 2500, 2500), settles_to_constraints_on: Union[PathLike, str] = None,
                 all_atoms: bool = False) -> None:
        if exclusion_list is None:
            exclusion_list = ["SOL", "NA", "CL"]
        out_dir = Path(out_dir)
        self.clean(self.solvated_dir)
        self.solvated_dir.mkdir(exist_ok=True, parents=True)
        out_dir.mkdir(exist_ok=True, parents=True)

        init_top = self.solvated_dir / 'init.top'
        init_gro = self.solvated_dir / 'init.gro'

        structure.save(str(init_top), overwrite=True)
        structure.save(str(init_gro), overwrite=True)

        self._include_water_ions_params(init_top)
        self._include_all_atom_types(init_top)
        self._add_water_and_ions(gro=init_gro, top=init_top, bt=bt, box=box, angles=angles,
                                 d=d, c=c, pname=pname, nname=nname, ion_conc=ion_conc, rmin=rmin)
        molecules = list(set(get_molecule_names(init_top)) - set(exclusion_list))
        make_posres(input_topology=init_top, molecules=molecules, out_dir=out_dir, f_xyz=f_xyz,
                    all_atoms=all_atoms)
        if settles_to_constraints_on:
            _tip3p_settles_to_constraints(top=init_top, molecule=settles_to_constraints_on, out_top=None)

        shutil.copy(init_top, out_dir / f'{out_name}.top')
        shutil.copy(init_gro, out_dir / f'{out_name}.gro')


def _read_parmed_molecule(top_file: PathLike, gro_file: PathLike) -> Structure:
    """Read a GROMACS top+gro into a single parmed Structure (monolithic)."""
    from parmed import load_file
    return load_file(str(top_file), xyz=str(gro_file))


def get_molecule_atom_ranges(top_file: PathLike) -> List[Tuple[str, int, int]]:
    """Parse a GROMACS top and return [(moleculetype, first_atom, n_atoms), ...]
    in the order of the [ molecules ] section.  Atom indices are 0-based (as in
    .gro / tpr).  This mirrors ana_gmx.mmpbsa: ligand/receptor groups are built
    from atom ranges (moleculetype blocks), never from residue names, so peptide
    ligands (whose residues legitimately keep protein names like GLU/PRO) are
    selected exactly like small-molecule ligands."""
    top_file = Path(top_file)
    lines = top_file.read_text().splitlines()

    moltype_atoms = {}      # moleculetype -> n_atoms (from [ atoms ] section)
    cur_mol = None
    cur_section = None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('['):
            cur_section = stripped[1:stripped.find(']')].strip()
            if cur_section == 'moleculetype':
                # next non-comment line is the name
                continue
        if cur_section == 'moleculetype' and not stripped.startswith((';', '#')) and stripped:
            cur_mol = stripped.split()[0]
            continue
        if cur_section == 'atoms' and not stripped.startswith((';', '#')) and stripped:
            parts = stripped.split()
            if len(parts) >= 7 and parts[0].isdigit():
                # skip comment lines; count atoms of current molecule
                moltype_atoms.setdefault(cur_mol, 0)
                moltype_atoms[cur_mol] += 1

    # order from [ molecules ]
    ranges = []
    offset = 0
    in_molecules = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('[ molecules ]'):
            in_molecules = True
            continue
        if in_molecules:
            if stripped.startswith('['):
                break
            if not stripped or stripped.startswith(';') or stripped.startswith('#'):
                continue
            parts = stripped.split()
            molname, count = parts[0], int(parts[1])
            n_atoms = moltype_atoms.get(molname, 0)
            ranges.append((molname, offset, n_atoms))
            offset += n_atoms * count
    return ranges


def write_atomic_index(complex_top: PathLike, complex_gro: PathLike,
                       ligand_moltype: str = 'LIG',
                       receptor_moltype: str = 'Protein',
                       ndxout: PathLike = 'index.ndx') -> None:
    """Write an index file with RECEPTOR / LIGAND groups selected by **atom
    ranges of moleculetypes**, exactly in the spirit of ana_gmx.mmpbsa (which
    builds groups from atom indices, not residue names).

    For ABFE the ligand is the moleculetype named ``ligand_moltype`` (the
    peptide is renamed to LIG at the topology level, see
    ``rename_peptide_moleculetype``); the receptor is everything else (excluding
    water/ions), i.e. the remaining protein moleculetype(s).

    Parameters
    ----------
    complex_top : topology with explicit [ molecules ] section
    complex_gro : matching structure (only used to cross-check atom count)
    ligand_moltype : name of the moleculetype to select as LIGAND
    receptor_moltype : name of the moleculetype(s) selected as RECEPTOR
                       (usually 'Protein' - in a soluble complex the pdb2gmx
                       protein gets moleculetype name 'Protein_chain_X' after
                       renaming, but by default the receptor is simply the
                       non-LIG non-water part)
    ndxout : output index file
    """
    from parmed import load_file
    top_file = Path(complex_top)
    gro_file = Path(complex_gro)

    ranges = get_molecule_atom_ranges(top_file)
    # receptor = all non-ligand, non-water/ion moleculetypes
    receptor_atoms, ligand_atoms = [], []
    for molname, first, n in ranges:
        molname = molname.strip()
        if molname == ligand_moltype:
            ligand_atoms.extend(range(first, first + n))
        elif molname not in ('SOL', 'NA', 'CL', 'W', 'ION'):
            receptor_atoms.extend(range(first, first + n))

    # sanity check against the gro atom count
    with open(gro_file) as f:
        for _ in range(2):
            line = f.readline()
        n_atoms_gro = int(line.split()[0])
    if len(receptor_atoms) + len(ligand_atoms) + 0 > n_atoms_gro:
        raise ValueError(f"index atom count {len(receptor_atoms)+len(ligand_atoms)} "
                         f"exceeds .gro atom count {n_atoms_gro}")
    if not ligand_atoms:
        raise ValueError(f"no atoms found for moleculetype {ligand_moltype!r} in {top_file}")

    def _fmt(indices):
        return ''.join(f'{i+1:8d}' for i in indices)  # 1-based for GROMACS

    txt = (f'[ RECEPTOR ]\n{_fmt(sorted(receptor_atoms))}\n\n'
           f'[ LIGAND ]\n{_fmt(sorted(ligand_atoms))}\n')
    Path(ndxout).write_text(txt)
    logger.info(f"index written to {ndxout}: RECEPTOR {len(receptor_atoms)} atoms, "
                f"LIGAND {len(ligand_atoms)} atoms")


def rename_peptide_moleculetype(top_file: PathLike, new_name: str = 'LIG',
                                n_ligand_atoms: int = None) -> Union[str, None]:
    """Rename the moleculetype of the peptide ligand in a topological file to
    ``new_name`` (default LIG), **in place**.  Returns the old moleculetype
    name, or None if nothing was renamed.

    This is the key step that makes a peptide (whose residues keep their
    protein names GLU/PRO/...) appear as the ligand for FEP: the mdp template
    uses ``couple-moltype = LIG``, so the moleculetype of the ligand must be
    named LIG.  It only rewrites:
      - the moleculetype name in every ``[ moleculetype ]`` header, and
      - the corresponding line in ``[ molecules ]``,
    leaving all atom/residue names, charges and interactions untouched.

    The peptide moleculetype is identified among the *non-water*, *non-ion*
    moleculetypes listed in [ molecules ] taking the following priority:
      1. the one whose atom count matches ``n_ligand_atoms`` (exact, preferred),
      2. the one with the fewest atoms (typical: receptor > peptide),
      3. the single remaining moleculetype.
    In the standard pipeline the final topology is written by ParmEd, which
    names any multi-residue molecule 'system{N}'; the receptor is 'system1',
    the peptide 'system2' (unless a multi-chain receptor splits into several
    moleculetypes - then the atom-count match is the robust selector).
    """
    top_file = Path(top_file)
    text = top_file.read_text()

    # collect moleculetype names from [ molecules ] section
    mol_lines = []
    in_mol = False
    for line in text.splitlines():
        if line.strip().startswith('[ molecules ]'):
            in_mol = True
            continue
        if in_mol:
            if line.strip().startswith('['):
                break
            if line.strip() and not line.strip().startswith(';'):
                mol_lines.append(line.strip().split()[0])

    # find a moleculetype that is not water/ion and not already LIG
    biotypes = {n for n in mol_lines if n not in ('SOL', 'NA', 'CL', 'W', 'ION')}
    target = None
    ranges = get_molecule_atom_ranges(top_file)
    by_atoms = {name: n for name, _, n in ranges if name in biotypes}
    if n_ligand_atoms is not None:
        # exact atom-count match (most robust for multi-chain receptors)
        for name, n_atoms in by_atoms.items():
            if n_atoms == n_ligand_atoms:
                target = name
                break
    if target is None:
        if len(biotypes) == 1:
            target = next(iter(biotypes))
        elif len(biotypes) >= 2 and by_atoms:
            # receptor is usually the one with most atoms (protein), peptide the smaller one
            target = min(by_atoms, key=by_atoms.get)
    if target is None or target == new_name:
        return None
    if target.lower() in ('water', 'ions'):
        raise ValueError(f"cannot rename water/ion moleculetype {target!r}")

    # 1) [ moleculetype ] headers: replace name line
    new_text = []
    in_mt = False
    for line in text.splitlines():
        if line.strip().startswith('[ moleculetype ]'):
            in_mt = True
            new_text.append(line)
            continue
        if in_mt:
            stripped = line.strip()
            if not stripped or stripped.startswith((';', '#')):
                # comment/blank line inside the header; keep it and keep waiting
                new_text.append(line)
                continue
            parts = line.split()
            if parts and parts[0] == target:
                parts[0] = new_name
                new_text.append('  '.join(parts))
            else:
                new_text.append(line)
            in_mt = False
            continue
        new_text.append(line)
    text = '\n'.join(new_text) + '\n'

    # 2) [ molecules ] section: rename references
    lines = text.splitlines()
    in_mol = False
    out = []
    for line in lines:
        if line.strip().startswith('[ molecules ]'):
            in_mol = True
            out.append(line)
            continue
        if in_mol:
            if line.strip().startswith('['):
                in_mol = False
            elif line.strip() and not line.strip().startswith(';'):
                parts = line.strip().split()
                if parts[0] == target:
                    # preserve original spacing
                    out.append(line.replace(parts[0], new_name, 1))
                    continue
        out.append(line)
    top_file.write_text('\n'.join(out) + '\n')
    logger.info(f"renamed moleculetype {target!r} -> {new_name!r} in {top_file}")
    return target


if __name__ == '__main__':
    pass


def index_for_soluble_system(configuration_file: PathLike, ndxout: PathLike = "index.ndx",
                             ligand_name: str = "LIG", host_name: str = "Protein",
                             gmx: Gromacs = None, cwd: str = None) -> None:
    """Make the index file for soluble systems (RECEPTOR + LIGAND groups)."""
    tmpopt = tempfile.NamedTemporaryFile(suffix='.opt')
    tmpndx = tempfile.NamedTemporaryFile(suffix='.ndx')

    sele_RECEPTOR = f"\"RECEPTOR\" group {host_name}"
    sele_LIGAND = f"\"LIGAND\" resname {ligand_name}"
    logger.info("Groups in the index.ndx file:")
    logger.info(f"\t{sele_RECEPTOR}")
    logger.info(f"\t{sele_LIGAND}")
    sele_RECEPTOR += ";\n"
    sele_LIGAND += ";\n"
    with open(tmpopt.name, "w") as opt:
        opt.write(sele_RECEPTOR + sele_LIGAND)

    wdir = cwd or str(Path(configuration_file).resolve().parent)
    _gmx = gmx or Gromacs(working_dir=wdir)
    _gmx.run_gmx_with_expect('make_ndx', f=str(configuration_file), o=tmpndx.name,
                             expect_actions=[{'>': 'q\r'}])
    _gmx.run_gmx_with_expect('select', s=str(configuration_file), sf=tmpopt.name,
                             n=tmpndx.name, on=str(ndxout))

    # deleting the _f0_t0.000 marker in the file
    with open(ndxout, "r") as index:
        data = index.read().replace("_f0_t0.000", "")
    with open(ndxout, "w") as index:
        index.write(data)

    tmpopt.close()
    tmpndx.close()


if __name__ == '__main__':
    pass
