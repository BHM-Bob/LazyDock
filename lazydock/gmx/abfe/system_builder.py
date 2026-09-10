'''
Date: 2026-08-28
LastEditors: BHM-Bob 2262029386@qq.com
Description: system builder for ABFE, migrated from BindFlow bindflow.preparation.system_builder.
             gmx command layer replaced with LazyDock's Gromacs class; supports peptide ligands.
             Core logic migrated from BindFlow (https://github.com/IFMIFMIF/BindFlow)
'''
import logging
import os
import shutil
import socket
import warnings
from itertools import chain
from pathlib import Path
from typing import List, Union

from parmed.structure import Structure
from parmed.tools.actions import HMassRepartition
from toff import Parameterize

from lazydock.gmx.abfe import solvent
from lazydock.gmx.abfe._home import gmx_ff_data_dir
from lazydock.gmx.abfe.solvent import Solvate
from lazydock.gmx.run import Gromacs
from mbapy_lite.base import put_log, put_err

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]


def _pdb2gmx_args_to_kwargs(args_list: List[str]) -> dict:
    """Convert a list of CLI-style pdb2gmx args ('-flag' or '-flag value')
    into kwargs for the Gromacs gen_command layer (bool -> bare flag)."""
    kwargs = {}
    i, n = 0, len(args_list)
    while i < n:
        tok = args_list[i]
        if tok.startswith('-'):
            key = tok.lstrip('-')
            if i + 1 < n and not args_list[i + 1].startswith('-'):
                kwargs[key] = args_list[i + 1]
                i += 2
            else:
                kwargs[key] = True
                i += 1
        else:
            i += 1
    return kwargs


def get_gmx_ff(ff_code: str, out_dir: PathLike = '.') -> PathLike:
    """Get the GROMACS force field tarball from abfe data dir and extract it."""
    import tarfile
    out_dir = Path(out_dir).resolve()
    supported_ff = ['Slipids_2020', 'amber99sb-star-ildn']
    if ff_code not in supported_ff:
        raise ValueError(f"ff_code = {ff_code} is not valid. Choose between: {supported_ff}")
    fname = Path(gmx_ff_data_dir()) / f'{ff_code}.ff.tar.gz'
    if not fname.exists():
        raise FileNotFoundError(f"Gromacs force field not bundled: {fname}")
    tar = tarfile.open(fname, "r:gz")
    tar.extractall(out_dir)
    tar.close()
    return out_dir / f'{ff_code}.ff'


def system_combiner(**md_elements) -> Structure:
    """Sum up all the elements provided as keyword arguments (ParmEd structures)."""
    md_system = None
    for _, element in md_elements.items():
        if element:
            element_copy = element[:]
            if md_system is None:
                md_system = element_copy
            else:
                md_system += element_copy
    if md_system is None:
        raise RuntimeError(f"system_combiner failed with the inputs: {md_elements}")
    logger.info(f"The system was constructed as follows: {' + '.join([k for k, v in md_elements.items() if v])}")
    return md_system


def make_bindflow_dir(out_dir: PathLike, ligand_dir: PathLike, sys_dir: PathLike):
    """Copy and paste function to create the structure of the BindFlow directory."""
    out_dir = Path(out_dir)
    ligand_dir = Path(ligand_dir)
    sys_dir = Path(sys_dir)

    complex_out = out_dir / "complex"
    ligand_out = out_dir / "ligand"
    complex_out.mkdir(exist_ok=True, parents=True)
    ligand_out.mkdir(exist_ok=True, parents=True)

    for itp_ndx_file in chain(ligand_dir.rglob("*.itp"), ligand_dir.rglob("*.ndx")):
        shutil.copy(src=itp_ndx_file, dst=ligand_out)
    shutil.copyfile(src=ligand_dir / "solvated.gro", dst=ligand_out / "ligand.gro")
    shutil.copyfile(src=ligand_dir / "solvated.top", dst=ligand_out / "ligand.top")

    for itp_ndx_file in chain(sys_dir.rglob("*.itp"), sys_dir.rglob("*.ndx")):
        shutil.copy(src=itp_ndx_file, dst=complex_out)
    shutil.copyfile(src=sys_dir / "solvated.gro", dst=complex_out / "complex.gro")
    shutil.copyfile(src=sys_dir / "solvated.top", dst=complex_out / "complex.top")


def _finalize_peptide_top(top_path: PathLike, posres_dir: PathLike = None,
                          n_ligand_atoms: int = None) -> None:
    """Rename the peptide moleculetype in ``top_path`` to LIG (in place) and
    sync the matching position-restraint file (posres_<old>.itp ->
    posres_LIG.itp), updating the ``#include`` references in the top.

    Only the moleculetype identity is touched: residue names (GLU/PRO/...),
    atom types, charges and all interactions stay protein, which is essential
    for stable GROMACS runs and downstream analysis.  The moleculetype name is
    what ``couple-moltype = LIG`` in the FEP mdps and what the index writer
    key on.  ``n_ligand_atoms`` pins the exact ligand moleculetype by atom
    count (robust even when the receptor splits into several moleculetypes).
    """
    top_path = Path(top_path)
    old = solvent.rename_peptide_moleculetype(top_path, new_name='LIG',
                                              n_ligand_atoms=n_ligand_atoms)
    if old is None:
        return
    posres_dir = Path(posres_dir) if posres_dir else top_path.parent
    old_posres = posres_dir / f'posres_{old}.itp'
    if old_posres.exists():
        old_posres.rename(posres_dir / 'posres_LIG.itp')
    text = top_path.read_text()
    if f'posres_{old}.itp' in text:
        top_path.write_text(text.replace(f'posres_{old}.itp', 'posres_LIG.itp'))
    logger.info(f"{top_path}: peptide moleculetype {old!r} -> LIG (posres synced)")


class CRYST1:
    """https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1 """
    def __init__(self, line: str = None):
        if line:
            self.a = float(line[6:15])
            self.b = float(line[15:24])
            self.c = float(line[24:33])
            self.alpha = float(line[33:40])
            self.beta = float(line[40:47])
            self.gamma = float(line[47:54])
            self.sGroup = line[55:66]
            try:
                self.z = int(line[66:70])
            except ValueError:
                self.z = ""
            self.__is_init = True
        else:
            self.__is_init = False

    def from_pdb(self, file: PathLike):
        with open(file, 'r') as f:
            for line in f.readlines():
                if line.startswith('CRYST1'):
                    self.__init__(line)
                    self.__is_init = True
                    break
        if not self.__is_init:
            warnings.warn('from_pdb was not able to initialize CRYST1')

    def __getitem__(self, key):
        return self.__dict__[key]

    def string(self):
        return "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-12s%4s\n" % \
            (self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.sGroup, self.z)

    def __repr__(self):
        return self.string()


class MakeInputs:
    """Build the systems for an ABFE calculation (complex + ligand topologies).

    Migrated from BindFlow's MakeInputs; gmx calls go through LazyDock's Gromacs
    class, and ``peptide_definition`` adds support for peptide ligands treated
    as a protein (pdb2gmx) instead of a small molecule.
    """

    def __init__(self, protein: dict = None, host_name: str = "Protein", membrane: dict = None,
                 cofactor: dict = None, cofactor_selection: str = "resname COF",
                 cofactor_on_protein: bool = True, water_model: str = 'amber/tip3p',
                 custom_ff_path: Union[None, PathLike] = None, hmr_factor: Union[float, None] = None,
                 fix_protein: bool = False, solv_d: float = 1.5, solv_bt: str = "dodecahedron",
                 solv_rmin: float = 1, solv_ion_conc: float = 150E-3,
                 builder_dir: PathLike = 'builder', gmx: Gromacs = None,
                 peptide_definition: dict = None,
                 pdb2gmx_args: Union[None, List[str]] = None,
                 n_term: int = 0, c_term: int = 0, maxwarn: Union[int, None] = None):
        self.pdb2gmx_args = pdb2gmx_args or []
        self.n_term = 0 if n_term is None else int(n_term)
        self.c_term = 0 if c_term is None else int(c_term)
        self.maxwarn = maxwarn
        self.protein = protein
        self.host_name = host_name
        self.membrane = membrane
        self.cofactor = cofactor
        self.cofactor_on_protein = cofactor_on_protein
        self.cofactor_selection = cofactor_selection
        self.hmr_factor = hmr_factor
        self.water_model = water_model
        self.fix_protein = fix_protein
        self.solv_d = solv_d
        self.solv_bt = solv_bt
        self.solv_rmin = solv_rmin
        self.solv_ion_conc = solv_ion_conc
        self.peptide_definition = peptide_definition
        self._peptide_n_atoms = None
        self.wd = Path(builder_dir).resolve()
        self.wd.mkdir(exist_ok=True, parents=True)
        self.__self_was_called = False
        self.gmx = gmx

        if custom_ff_path:
            self.custom_ff_path = Path(custom_ff_path).resolve()
            # GMXLIB must point to the directory that CONTAINS the .ff dirs
            # (pdb2gmx -ff <name> searches <GMXLIB>/<name>.ff). If the user
            # passed the .ff dir itself, use its parent instead.
            if self.custom_ff_path.name.endswith('.ff'):
                gmxlib = str(self.custom_ff_path.parent)
            else:
                gmxlib = str(self.custom_ff_path)
            os.environ["GMXLIB"] = gmxlib
            put_log(f'custom_ff_path={self.custom_ff_path} -> GMXLIB={gmxlib}', head='ABFE')
        else:
            self.custom_ff_path = None

        # Initialize vectors and angles from PDB only for membrane systems
        if self.membrane:
            cryst_info = CRYST1()
            cryst_info.from_pdb(self.membrane['conf'])
            self.vectors = (cryst_info.a / 10, cryst_info.b / 10, cryst_info.c / 10)
            self.angles = (cryst_info.alpha, cryst_info.beta, cryst_info.gamma)
            logger.info(f"This is a membrane system. Crystal information was taken from: {self.membrane}")
        else:
            self.vectors, self.angles = None, None

        self.cwd = os.getcwd()

    def small_mol_process(self, mol_definition: dict, name: str = "MOL", safe_naming_prefix: str = None) -> Structure:
        """Get parameters for small molecules: ligands, cofactors... (BindFlow logic)."""
        force_field_code_default = {
            'openff': 'openff_unconstrained-2.0.0.offxml',
            'gaff': 'gaff-2.11',
            'espaloma': 'espaloma-0.3.1'
        }
        dict_to_work = {'top': None, 'ff': {'type': 'openff', 'code': None}}
        if mol_definition:
            _recursive_update_dict(dict_to_work, mol_definition)
            if isinstance(dict_to_work['ff']['type'], str):
                dict_to_work['ff']['type'] = str(dict_to_work['ff']['type']).lower()
            if dict_to_work['ff']['type'] not in force_field_code_default and not dict_to_work['top']:
                raise ValueError(f"Molecule {dict_to_work} has non valid type for the force field. "
                                 f"Choose from {force_field_code_default.keys()}.")
            if not dict_to_work['ff']['code'] and dict_to_work['ff']['type']:
                dict_to_work['ff']['code'] = force_field_code_default[dict_to_work['ff']['type']]
        else:
            raise ValueError(f"Molecule {mol_definition} has a wrong definition.")
        if dict_to_work['conf']:
            if dict_to_work['top']:
                logger.info(f"Using supplied: {dict_to_work['top']} for {dict_to_work['conf']}")
            else:
                logger.info(f"Getting {dict_to_work['ff']['code']} (type = {dict_to_work['ff']['type']}) "
                            f"parameters for: {dict_to_work['conf']}")
        else:
            raise ValueError(f"Molecule {mol_definition} has a wrong configuration")

        provided_top_flag = False
        if dict_to_work['top']:
            top_file = Path(dict_to_work['top']).resolve()
            provided_top_flag = True
            if Path(dict_to_work['conf']).suffix == '.gro':
                gro_file = Path(dict_to_work['conf']).resolve()
            else:
                raise ValueError("For safety reasons, if top is provided for small molecule; "
                                 f"the gro file must be provided. Provided: {dict_to_work['conf']}.")
        else:
            parameterizer = Parameterize(
                force_field_code=dict_to_work['ff']['code'],
                force_field_type=dict_to_work['ff']['type'],
                ext_types=['top', 'gro'],
                hmr_factor=self.hmr_factor,
                overwrite=True,
                safe_naming_prefix=safe_naming_prefix,
                out_dir=self.wd,
            )
            parameterizer(input_mol=dict_to_work['conf'], mol_resi_name=name)
            top_file = self.wd / f"{name}.top"
            gro_file = self.wd / f"{name}.gro"

        parmed_system = _read_parmed_molecule(top_file=top_file, gro_file=gro_file)
        if provided_top_flag and self.hmr_factor:
            HMassRepartition(parmed_system, self.hmr_factor).execute()
        return parmed_system

    def gmx_process(self, mol_definition: dict, is_membrane: bool = False) -> Structure:
        """Process the compatible biomolecules (protein, peptide, membrane) using pdb2gmx via LazyDock Gromacs."""
        dict_to_work = {'top': None, 'ff': {'code': 'Slipids_2020' if is_membrane else 'amber99sb-ildn'}}
        if mol_definition:
            _recursive_update_dict(dict_to_work, mol_definition)
        else:
            return None
        if dict_to_work['conf']:
            dict_to_work['conf'] = Path(dict_to_work['conf']).resolve()
            name, ext = dict_to_work['conf'].stem, dict_to_work['conf'].suffix
            if dict_to_work['top']:
                dict_to_work['top'] = Path(dict_to_work['top']).resolve()
                logger.info(f"Using supplied: {dict_to_work['top']} for {dict_to_work['conf']}")
            else:
                logger.info(f"Getting {dict_to_work['ff']['code']} parameters for: {dict_to_work['conf']}")
        else:
            return None

        gro_out = self.wd / f'{name}.gro'
        top_out = self.wd / f'{name}.top'
        posre_out = self.wd / f'{name}_posre.itp'

        if is_membrane:
            if dict_to_work['ff']['code'] == 'Slipids_2020':
                if self.custom_ff_path:
                    if 'Slipids_2020' not in list(self.custom_ff_path.iterdir()):
                        get_gmx_ff('Slipids_2020', out_dir=self.wd)
                else:
                    get_gmx_ff('Slipids_2020', out_dir=self.wd)

        os.chdir(self.wd)
        gmx = self.gmx or Gromacs(working_dir=str(self.wd))
        check_box = False
        if dict_to_work['top']:
            shutil.copy(dict_to_work['top'], top_out)
            if ext == '.pdb':
                check_box = True
                gmx.run_gmx_with_expect('editconf', f=str(dict_to_work['conf']), o=str(gro_out))
            elif ext == '.gro':
                shutil.copy(dict_to_work['conf'], gro_out)
            else:
                raise ValueError(f"Extension of {dict_to_work['conf']} must be .gro or .pdb")
        else:
            if is_membrane:
                gmx.run_gmx_with_expect('pdb2gmx', f=str(dict_to_work['conf']), ff=dict_to_work['ff']['code'],
                                        water="none", o=str(gro_out), p=str(top_out), i=str(posre_out))
            else:
                _extra_pdb2gmx = _pdb2gmx_args_to_kwargs(self.pdb2gmx_args)
                if self.fix_protein:
                    logger.info(f"fix_protein = {self.fix_protein}; therefore pdbfixer will be used "
                                "with flags --add-atoms=all --replace-nonstandard and pdb2gmx with -ignh. "
                                "The protonation of your protein may have changed!")
                    env_prefix = os.environ.get("CONDA_PREFIX", "")
                    fixed_pdb = self.wd / f"{name}_fixed.pdb"
                    _run_cmd(f"{env_prefix}/bin/pdbfixer {dict_to_work['conf']} "
                             f"--output={fixed_pdb} --add-atoms=all --replace-nonstandard")
                    _strip_pdbfixer_oxt(dict_to_work['conf'], fixed_pdb)
                    # pdb2gmx: -ff is not used; -ignh and -merge all, with term selection via expect
                    # NOTE: -ff/-water are passed on the command line, so pdb2gmx will
                    # NOT prompt for force field / water model.  Only the termini
                    # prompts (with -ter) remain; matching patterns for non-existent
                    # prompts would hang forever.
                    _pdb2gmx_expect = [
                        {'Select start terminus type': f'{self.n_term}\r'},
                        {'Select end terminus type': f'{self.c_term}\r'},
                    ]
                    # NOTE: explicit flags below are defaults; user-supplied
                    # --pdb2gmx-args overwrite them (e.g. -noignh).
                    _pdb2gmx_kwargs = dict(f=str(fixed_pdb), merge="all",
                                           ff=dict_to_work['ff']['code'], water="none",
                                           o=str(gro_out), p=str(top_out), i=str(posre_out),
                                           ter=True, ignh=True)
                    _pdb2gmx_kwargs.update(_extra_pdb2gmx)
                    gmx.run_gmx_with_expect('pdb2gmx', expect_actions=_pdb2gmx_expect,
                                            expect_settings={'timeout': 180}, **_pdb2gmx_kwargs)
                else:
                    _pdb2gmx_expect = [
                        {'Select start terminus type': f'{self.n_term}\r'},
                        {'Select end terminus type': f'{self.c_term}\r'},
                    ]
                    # NOTE: explicit flags below are defaults; user-supplied
                    # --pdb2gmx-args overwrite them (e.g. -noignh).
                    _pdb2gmx_kwargs = dict(f=str(dict_to_work['conf']), merge="all",
                                           ff=dict_to_work['ff']['code'], water="none",
                                           o=str(gro_out), p=str(top_out), i=str(posre_out),
                                           ter=True, ignh=True)
                    _pdb2gmx_kwargs.update(_extra_pdb2gmx)
                    gmx.run_gmx_with_expect('pdb2gmx', expect_actions=_pdb2gmx_expect,
                                            expect_settings={'timeout': 180}, **_pdb2gmx_kwargs)
        os.chdir(self.cwd)

        system = solvent._read_parmed_molecule(top_file=top_out, gro_file=gro_out)
        _check_and_fix_box(system, check_box)
        if self.hmr_factor:
            HMassRepartition(system, self.hmr_factor).execute()
        system.write(str(self.wd / f'{name}_final.top'))
        return system

    def make_system(self, ligand_definition: dict):
        """Create self.sys_ligand, self.sys_protein, self.sys_membrane and self.md_system."""
        logger.info("Processing system components")
        # Peptide ligand: treated as a protein (pdb2gmx).  Its residue names
        # must stay protein (GLU/PRO/...) so that GROMACS/Boresch/analysis treat
        # it as a peptide; the ligand identity for FEP (couple-moltype LIG) is
        # enforced by renaming ONLY the moleculetype (see make_bindflow_dir).
        if self.peptide_definition:
            pep_def = dict(self.peptide_definition)
            pep_def['conf'] = Path(pep_def['conf']).resolve()
            logger.info(f"Processing peptide ligand as protein: {pep_def['conf']}")
            self.sys_ligand = self.gmx_process(mol_definition=pep_def)
            self._peptide_n_atoms = self.sys_ligand.n_atoms if hasattr(self.sys_ligand, 'n_atoms') \
                else len(self.sys_ligand.atoms)
        else:
            self.sys_ligand = self.small_mol_process(mol_definition=ligand_definition, name="LIG",
                                                     safe_naming_prefix='x')
            self._peptide_n_atoms = None

        if self.__self_was_called:
            logger.info("Reusing components from cache")
        else:
            if self.cofactor:
                self.sys_cofactor = self.small_mol_process(mol_definition=self.cofactor, name="COF",
                                                           safe_naming_prefix='z')
            else:
                self.sys_cofactor = None
            self.sys_protein = self.gmx_process(mol_definition=self.protein)
            self.sys_membrane = self.gmx_process(mol_definition=self.membrane, is_membrane=True)
        logger.info("Merging Components")
        self.md_system = system_combiner(protein=self.sys_protein, membrane=self.sys_membrane,
                                         ligand=self.sys_ligand, cofactor=self.sys_cofactor)

    def clean(self):
        """Small cleaner: the intermediate steps saved on builder_dir will be deleted."""
        os.chdir(self.cwd)
        try:
            shutil.rmtree(self.wd)
        except FileNotFoundError:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.clean()

    def __call__(self, ligand_definition: Union[dict, PathLike], out_dir: str = 'fep'):
        """The call implementation; builds all components and solvates complex and ligand."""
        logger.info(39 * "-")
        logger.info(f"Running on compute host: {socket.gethostname()}")
        if not isinstance(ligand_definition, dict):
            ligand_definition = {'conf': ligand_definition}
        logger.info(f"Processing ligand: {ligand_definition['conf']}")
        self.out_dir = Path(out_dir)
        self.out_dir.mkdir(exist_ok=True, parents=True)

        self.make_system(ligand_definition)
        system_dir = self.wd / 'system'
        ligand_dir = self.wd / 'ligand'

        logger.info(f"Solvating with {self.water_model}:")
        f_xyz_complex = 3 * [2500] if not self.membrane else 3 * ['POSRES_DYNAMIC']

        gmx = self.gmx or Gromacs(working_dir=str(self.wd))
        with Solvate(self.water_model, builder_dir=self.wd / '.solvating', gmx=gmx,
                     cwd=self.cwd, maxwarn=self.maxwarn) as SolObj:
            logger.info(f"Ligand in: {ligand_dir}")
            SolObj(structure=self.sys_ligand, bt=self.solv_bt, d=self.solv_d, rmin=self.solv_rmin,
                   ion_conc=self.solv_ion_conc, out_dir=ligand_dir, out_name='solvated',
                   f_xyz=3 * [2500], all_atoms=bool(self.peptide_definition))
            logger.info(f"Complex in: {system_dir}")
            settles_to_constraints_on = None
            if self.cofactor and self.cofactor.get('is_water'):
                warnings.warn(f"Provided cofactor {self.cofactor} was labeled as water (is_water = True). "
                              "Its settles section will be changed to tip3p-like triangular constraints.")
                settles_to_constraints_on = 'COF'
            if self.membrane:
                SolObj(structure=self.md_system, bt='triclinic', box=self.vectors, angles=self.angles,
                       rmin=self.solv_rmin, ion_conc=self.solv_ion_conc, out_dir=system_dir,
                       out_name='solvated', f_xyz=f_xyz_complex,
                       settles_to_constraints_on=settles_to_constraints_on)
            else:
                SolObj(structure=self.md_system, bt=self.solv_bt, d=self.solv_d, rmin=self.solv_rmin,
                       ion_conc=self.solv_ion_conc, out_dir=system_dir, out_name='solvated',
                       f_xyz=f_xyz_complex, settles_to_constraints_on=settles_to_constraints_on)

        # Make index file
        if self.membrane:
            solvent.index_for_membrane_system(
                configuration_file=system_dir / "solvated.gro", ndxout=system_dir / "index.ndx",
                ligand_name="LIG", host_name=self.host_name,
                cofactor_selection=self.cofactor_selection if self.cofactor else None,
                cofactor_on_protein=self.cofactor_on_protein, gmx=gmx, cwd=str(system_dir))
        else:
            (system_dir / "index.ndx").touch(exist_ok=True)
            if self.peptide_definition:
                # Peptide ligand: keep the peptide as a full protein in the
                # topology (residue names GLU/PRO/..., atom types, charges,
                # ring bonds) for stable GROMACS runs and downstream analysis.
                # For FEP the ligand moleculetype must be named 'LIG'
                # (couple-moltype = LIG in the mdps), so we rename ONLY the
                # moleculetype here (in place, before the index is written and
                # before make_bindflow_dir copies the files), and the index is
                # written by atom ranges of moleculetypes (mmpbsa-style)
                # instead of 'resname LIG'.
                _finalize_peptide_top(system_dir / "solvated.top", posres_dir=system_dir,
                                       n_ligand_atoms=self._peptide_n_atoms)
                solvent.write_atomic_index(complex_top=system_dir / "solvated.top",
                                           complex_gro=system_dir / "solvated.gro",
                                           ligand_moltype='LIG',
                                           receptor_moltype=self.host_name,
                                           ndxout=system_dir / "index.ndx")
            else:
                solvent.index_for_soluble_system(
                    configuration_file=system_dir / "solvated.gro", ndxout=system_dir / "index.ndx",
                    ligand_name="LIG", host_name=self.host_name, gmx=gmx, cwd=str(system_dir))

        # iso-solvent box for the ligand (peptide): its moleculetype also has to
        # be called LIG for FEP
        if self.peptide_definition:
            _finalize_peptide_top(ligand_dir / "solvated.top", posres_dir=ligand_dir,
                                  n_ligand_atoms=self._peptide_n_atoms)

        logger.info(f"Final build of BindFlow directory on: {self.out_dir}")
        make_bindflow_dir(out_dir=self.out_dir, ligand_dir=ligand_dir, sys_dir=system_dir)

        self.__self_was_called = True
        logger.info("--------- Building Completed ----------\n")


#############################################################################
# small helpers
#############################################################################

def _recursive_update_dict(target: dict, new_elements: dict):
    """Recursively update target with new_elements."""
    for key, value in new_elements.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _recursive_update_dict(target[key], value)
        else:
            target[key] = value


def _read_parmed_molecule(top_file: PathLike, gro_file: PathLike, check_box: bool = False) -> Structure:
    """Read a GROMACS top+gro into a single parmed Structure."""
    from parmed import load_file
    return load_file(str(top_file), xyz=str(gro_file))


def _check_and_fix_box(system: Structure, check_box: bool):
    """If check_box and the structure lacks box info, set a zero box (avoid solvation crashes)."""
    if check_box:
        try:
            if system.box is None or system.box[0] == 0:
                system.box = [0, 0, 0, 90, 90, 90]
        except Exception:
            pass


def _run_cmd(cmd: str):
    """Run a shell command (used for pdbfixer), raising on failure."""
    import subprocess
    put_log(cmd, head='ABFE')
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"command failed ({ret}): {cmd}")


def _strip_pdbfixer_oxt(input_conf: PathLike, fixed_pdb: PathLike):
    """Remove OXT atoms that pdbfixer added to chains without an OXT in the input.

    pdbfixer --add-atoms=all treats every chain as a linear peptide and adds a
    C-terminal OXT, which breaks cyclic peptides (pdb2gmx then fails with
    'Atom OXT ... not found in rtp entry').  Only chains whose input PDB
    actually has an OXT keep it; the ring-closure decision is left to pdb2gmx.
    """
    input_conf = Path(input_conf)
    fixed_pdb = Path(fixed_pdb)
    input_chains_with_oxt = set()
    for line in input_conf.read_text().splitlines():
        if line.startswith('ATOM') and line[12:16].strip() == 'OXT':
            input_chains_with_oxt.add(line[21])
    lines = fixed_pdb.read_text().splitlines()
    new_lines = []
    for line in lines:
        if line.startswith('ATOM') and line[12:16].strip() == 'OXT' and line[21] not in input_chains_with_oxt:
            continue
        new_lines.append(line)
    if len(new_lines) != len(lines):
        fixed_pdb.write_text('\n'.join(new_lines) + '\n')
        put_log(f'stripped pdbfixer-added OXT from chains '
                f'{sorted(set(l[21] for l in lines if l.startswith("ATOM") and l[12:16].strip()=="OXT") - input_chains_with_oxt)} '
                f'in {fixed_pdb}', head='ABFE')


if __name__ == '__main__':
    pass
