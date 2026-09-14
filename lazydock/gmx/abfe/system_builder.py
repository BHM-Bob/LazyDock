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
from itertools import chain
from pathlib import Path
from typing import Union

from parmed.structure import Structure
from parmed.tools.actions import HMassRepartition
from toff import Parameterize

from lazydock.gmx.abfe import solvent
from lazydock.gmx.abfe.solvent import Solvate
from lazydock.gmx.run import Gromacs
from mbapy_lite.base import put_log

logger = logging.getLogger(__name__)

PathLike = Union[str, Path]


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


class MakeInputs:
    """Build the systems for an ABFE calculation (complex + ligand topologies).

    Migrated from BindFlow's MakeInputs; gmx calls go through LazyDock's Gromacs
    class, and ``peptide_definition`` adds support for peptide ligands treated
    as a protein (pdb2gmx) instead of a small molecule.
    """

    def __init__(self, protein: dict = None, host_name: str = "Protein",
                 water_model: str = 'amber/tip3p',
                 custom_ff_path: Union[None, PathLike] = None, hmr_factor: Union[float, None] = None,
                 solv_d: float = 1.5, solv_bt: str = "dodecahedron",
                 solv_rmin: float = 1, solv_ion_conc: float = 150E-3,
                 builder_dir: PathLike = 'builder', gmx: Gromacs = None,
                 peptide_definition: dict = None,
                 maxwarn: Union[int, None] = None,
                 fep_rlist_ligand: Union[float, None] = None,
                 solv_d_ligand: Union[float, None] = None,
                 solv_d_complex: Union[float, None] = None,
                 auto_box: bool = False, auto_box_pad: Union[float, list] = 1.2,
                 mono_lock=None, ligand_chain: str = None,
                 whole_complex: bool = False):
        self.maxwarn = maxwarn
        self.protein = protein
        self.host_name = host_name
        self.hmr_factor = hmr_factor
        self.water_model = water_model
        self.solv_d = solv_d
        self.solv_bt = solv_bt
        self.solv_rmin = solv_rmin
        self.solv_ion_conc = solv_ion_conc
        self.solv_d_ligand = solv_d_ligand
        self.solv_d_complex = solv_d_complex
        self.fep_rlist_ligand = fep_rlist_ligand
        self.auto_box = auto_box
        self.auto_box_pad = auto_box_pad
        self.mono_lock = mono_lock
        self.ligand_chain = ligand_chain
        self.whole_complex = whole_complex
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

    def gmx_process(self, mol_definition: dict) -> Union[Structure, None]:
        """Read an already-prepared GROMACS system (top + gro) into a parmed
        Structure, apply HMR and rewrite a monolithic top.

        All pdb2gmx/prepare work happens in prepare-gmx (protein/complex);
        this method only *consumes* its outputs (top/gro pair).  Membrane and
        small-molecule paths were removed (user review 2026-09-12): membrane
        handling will not be part of ABFE; small molecules will go through
        prepare-gmx complex.
        """
        dict_to_work = {'top': None, 'ff': {'code': 'amber99sb-ildn'}}
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
                raise ValueError(f"gmx_process requires a pre-existing topology "
                                 f"('top' key) in {mol_definition}: all pdb2gmx is "
                                 f"done by prepare-gmx.")
        else:
            return None

        gro_out = self.wd / f'{name}.gro'
        top_out = self.wd / f'{name}.top'

        os.chdir(self.wd)
        if dict_to_work['top']:
            shutil.copy(dict_to_work['top'], top_out)
            if ext == '.pdb':
                gmx = self.gmx or Gromacs(working_dir=str(self.wd))
                gmx.run_gmx_with_expect('editconf', f=str(dict_to_work['conf']), o=str(gro_out))
            elif ext == '.gro':
                shutil.copy(dict_to_work['conf'], gro_out)
            else:
                raise ValueError(f"Extension of {dict_to_work['conf']} must be .gro or .pdb")
        os.chdir(self.cwd)

        system = solvent._read_parmed_molecule(top_file=top_out, gro_file=gro_out)
        if self.hmr_factor:
            HMassRepartition(system, self.hmr_factor).execute()
        system.write(str(self.wd / f'{name}_final.top'))
        return system

    def make_system(self, ligand_definition: dict):
        """Create self.sys_ligand, self.sys_protein and self.md_system.

        Only the whole-complex (peptide) mode is supported: complex leg is the
        whole complex.pdb prepared by prepare-gmx protein once (-merge no keeps
        the two moleculetypes, relative coordinates are naturally preserved);
        sys_ligand is the independent peptide protein topology prepared by
        prepare-gmx protein in ligand_pep/ (user decision 2026-09-11).  Small
        molecules will go through prepare-gmx complex (toff/CGenFF) in the
        future - for now raise NotImplementedError instead of silently falling
        into the old protein-based path (user review 2026-09-12).
        """
        logger.info("Processing system components")
        if self.whole_complex:
            self.sys_protein = self.gmx_process(mol_definition=self.protein)
            self.sys_membrane = None
            self.sys_cofactor = None
            self.md_system = self.sys_protein
            if isinstance(ligand_definition, dict) and ligand_definition.get('conf'):
                _lg = solvent._read_parmed_molecule(top_file=ligand_definition.get('top'),
                                                    gro_file=ligand_definition['conf'])
                self.sys_ligand = _lg
                self._peptide_n_atoms = len(_lg.atoms)
                self._ligand_parmed = _lg
                logger.info(f"whole-complex: sys_protein = whole "
                            f"({len(self.sys_protein.atoms)} atoms), sys_ligand = "
                            f"independent peptide protein ({self._peptide_n_atoms} atoms)")
            else:
                # 兜底: 按配体链原子数从整链切出配体
                n_lig = self._find_ligand_n_atoms_by_chain(self.sys_protein)
                n_total = self.sys_protein.n_atoms if hasattr(self.sys_protein, 'n_atoms') \
                    else len(self.sys_protein.atoms)
                self.sys_ligand = self.sys_protein[n_total - n_lig:n_total]
                self._peptide_n_atoms = n_lig
                logger.info(f"whole-complex (fallback slice): sys_ligand = "
                            f"{n_total - n_lig}:{n_total} ({n_lig} atoms)")
        else:
            raise NotImplementedError(
                'only whole_complex (peptide) mode is supported by make_system; '
                'small-molecule ABFE will go through prepare-gmx complex '
                '(toff/CGenFF) - see docs/dev/abfe/small_mol_future_prepare_gmx_complex.md')

    def _find_ligand_n_atoms_by_chain(self, whole) -> int:
        """Fallback: infer the ligand atom count by counting atoms of the
        ligand chain from the whole-complex parmed residue/atom names."""
        if self.ligand_chain is None:
            raise ValueError("whole-complex mode requires ligand_chain or a "
                             "ligand_definition with a known atom count")
        try:
            cnt = 0
            for a in whole.atoms:
                if str(a.chain) == self.ligand_chain:
                    cnt += 1
            if cnt > 0:
                return cnt
        except Exception as e:
            logger.warning(f"ligand-chain atom counting failed: {e}")
        raise ValueError(f"cannot determine ligand atom count from whole-complex "
                         f"(chain {self.ligand_chain})")

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
        f_xyz_complex = 3 * [2500]

        gmx = self.gmx or Gromacs(working_dir=str(self.wd))
        with Solvate(self.water_model, builder_dir=self.wd / '.solvating', gmx=gmx,
                     cwd=self.cwd, maxwarn=self.maxwarn) as SolObj:
            logger.info(f"Ligand in: {ligand_dir}")
            SolObj(structure=self.sys_ligand, bt=self.solv_bt, d=self.solv_d_ligand or self.solv_d,
                   rmin=self.solv_rmin, ion_conc=self.solv_ion_conc, out_dir=ligand_dir,
                   out_name='solvated', f_xyz=3 * [2500], all_atoms=bool(self.peptide_definition),
                   rlist_min=self.fep_rlist_ligand, auto_box=self.auto_box,
                   auto_box_pad=self.auto_box_pad, mono_lock=self.mono_lock)
            logger.info(f"Complex in: {system_dir}")
            SolObj(structure=self.md_system, bt=self.solv_bt, d=self.solv_d_complex or self.solv_d,
                   rmin=self.solv_rmin, ion_conc=self.solv_ion_conc, out_dir=system_dir,
                   out_name='solvated', f_xyz=f_xyz_complex,
                   rlist_min=self.fep_rlist_ligand, auto_box=self.auto_box,
                   auto_box_pad=self.auto_box_pad, mono_lock=self.mono_lock)

        # Make index file
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


def _read_parmed_molecule(top_file: PathLike, gro_file: PathLike) -> Structure:
    """Read a GROMACS top+gro into a single parmed Structure."""
    from parmed import load_file
    return load_file(str(top_file), xyz=str(gro_file))


if __name__ == '__main__':
    pass
