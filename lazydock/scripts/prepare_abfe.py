'''
Date: 2026-08-28
LastEditors: BHM-Bob 2262029386@qq.com
Description: prepare ABFE systems (complex + ligand topologies).
             Supports receptor multi-chain and peptide ligands (as protein).
'''
import argparse
import logging
import os
from pathlib import Path
from typing import List, Union

from mbapy_lite.base import put_err, put_log
from mbapy_lite.file import get_paths_with_extension

# system_builder imported lazily in _get_builder (avoids openff probe on CLI parse)
from lazydock.gmx.run import Gromacs
from lazydock.gmx.abfe import system_builder
from lazydock.scripts._script_utils_ import Command, make_args_and_excute, process_batch_dir_lst

logger = logging.getLogger(__name__)


def _fix_ligand_for_abfe(pdb_path: Path, ligand_resname: str = 'LIG') -> Path:
    """Return the pdb path with the ligand residue renamed (in-place edit of a copy)."""
    from mbapy_lite.file import opts_file
    lines = opts_file(pdb_path, way='lines')
    new_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            line = f'{line[:17]}{ligand_resname:<3s}{line[20:]}'
        new_lines.append(line)
    out = pdb_path.with_suffix('.ligrenamed.pdb')
    opts_file(out, 'w', way='lines', data=new_lines)
    return out


def _sanitize_pdb_for_pymol(pdb_path: Path) -> Path:
    """Copy a PDB, replacing NaN/Inf/huge numeric fields with sane values.

    Some third-party PDBs carry garbage in the numeric columns (coordinates,
    occupancy, b-factor). PyMOL's Map code (voxel maps built e.g. during
    load/save) prints "clamping Min/Max ..." and — on some PyMOL builds —
    segfaults on such data. Sanitizing before handing the file to PyMOL
    makes the pipeline robust anywhere. Returns the sanitized copy path.
    """
    import math

    def _clean_number(raw: str, default: float) -> str:
        try:
            v = float(raw)
            if not math.isfinite(v):
                v = default
        except ValueError:
            v = default
        return v

    out = pdb_path.with_name(pdb_path.stem + '_sanitized.pdb')
    new_lines = []
    with open(pdb_path, 'r', errors='replace') as fh:
        for line in fh:
            if line.startswith(('ATOM', 'HETATM')) and len(line) >= 66:
                xyz = [_clean_number(line[30:38], 0.0),
                       _clean_number(line[38:46], 0.0),
                       _clean_number(line[46:54], 0.0)]
                occ = _clean_number(line[54:60], 1.0)
                bf = _clean_number(line[60:66], 0.0)
                line = (f'{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}'
                        f'{occ:6.2f}{bf:6.2f}{line[66:]}')
            new_lines.append(line)
    out.write_text(''.join(new_lines))
    logger.info(f'sanitized PDB for PyMOL: {out}')
    return out


class Complex(Command):
    HELP = """
    prepare ABFE system: complex + ligand.
    
    INPUT:
        complex.pdb with two (or more) chains: receptor chain(s) + ligand chain.
        - receptor can have MULTIPLE chains (--receptor-chains A B ...)
        - ligand can be a PEPTIDE (--peptide, treated as protein via pdb2gmx)
          or a small molecule (espaloma/openff/gaff parameterization)

    OUTPUT (BindFlow-compatible dir structure):
        input/complex/  complex.gro complex.top + itps + index.ndx
        input/ligand/   ligand.gro ligand.top + itps
    """
    def __init__(self, args, printf=print):
        super().__init__(args, printf, iter_run_arg=['dir'])

    @staticmethod
    def make_args(args: argparse.ArgumentParser):
        args.add_argument('-d', '--dir', type=str, nargs='+', default=['.'],
                          help='directory containing complex.pdb. Default is %(default)s.')
        args.add_argument('-n', '--name', type=str, default='complex.pdb',
                          help='complex structure file name in each dir. Default is %(default)s.')
        args.add_argument('-rc', '--receptor-chains', type=str, nargs='+', required=True,
                          help='receptor chain name(s), e.g. A or A B C (multi-chain supported).')
        args.add_argument('-lc', '--ligand-chain', type=str, required=True,
                          help='ligand chain name, e.g. Z or P.')
        args.add_argument('--peptide', action='store_true',
                          help='treat ligand as a PEPTIDE (protein-protein form, pdb2gmx).')
        args.add_argument('--ligand-ff', type=str, default='espaloma', choices=['espaloma', 'openff', 'gaff'],
                          help='small molecule force field backend (ignored if --peptide). Default %(default)s.')
        args.add_argument('--ff-dir', type=str, default=None,
                          help='force field directory (custom_ff_path -> GMXLIB). '
                               'Can be a charmm36 dir (e.g. charmm36-jul2022.ff) or amber dir. '
                               'If None, pdb2gmx will use its installed force fields.')
        args.add_argument('--protein-ff', type=str, default='amber99sb-ildn',
                          help='force field code for pdb2gmx (protein/peptide), e.g. amber99sb-ildn, '
                               'charmm36-mar2019, charmm36-jul2022. '
                               'A path ending with .ff is also accepted (the dir name is used). '
                               'Default %(default)s.')
        args.add_argument('--pdb2gmx-args', type=str, default='-ter -ignh',
                          help='args passed to pdb2gmx. Default %(default)s.')
        args.add_argument('--n-term', type=int, default=0, help='N-term type for pdb2gmx (0 or 1). Default %(default)s.')
        args.add_argument('--c-term', type=int, default=0, help='C-term type for pdb2gmx (0 or 1). Default %(default)s.')
        args.add_argument('--water-model', type=str, default='amber/tip3p',
                          help='water model. Default %(default)s.')
        args.add_argument('--solv-d', type=float, default=1.5, help='editconf -d. Default %(default)s.')
        args.add_argument('--solv-bt', type=str, default='dodecahedron', help='box type. Default %(default)s.')
        args.add_argument('--ion-conc', type=float, default=150e-3, help='ion concentration (M). Default %(default)s.')
        args.add_argument('--hmr-factor', type=float, default=2.5,
                          help='HMR factor (Hydrogen Mass Repartition). Default %(default)s '
                               'matching BindFlow FEP driver (fep_full_run.py hard-codes 2.5). '
                               'Must be used with 4 fs FEP timestep (dt_max >= 0.004). '
                               'Pass 0 to disable HMR (then use dt_max <= 0.002).')
        args.add_argument('--maxwarn', type=int, default=0, help='maxwarn for grompp. Default %(default)s.')
        args.add_argument('--fix-protein', action='store_true', default=False,
                          help='run pdbfixer (--add-atoms=all --replace-nonstandard) before pdb2gmx. '
                               'Default off: pdbfixer adds a C-terminal OXT to every chain and breaks '
                               'cyclic peptides (see lazydock/opmm/relax.py).')
        args.add_argument('--builder-dir', type=str, default='builder', help='builder dir name. Default %(default)s.')
        return args

    def process_args(self):
        self.args.dir = process_batch_dir_lst(self.args.dir)
        if self.args.peptide and not self.args.ff_dir:
            put_err('--peptide requires --ff-dir (protein force field dir for pdb2gmx).', _exit=True)
        # espaloma/openff availability is checked lazily in main_process (keeps CLI parse fast).

    def main_process(self):
        # --protein-ff: accept either a force field name or a .ff dir path
        ff = self.args.protein_ff
        if isinstance(ff, str) and (ff.rstrip('/').endswith('.ff') or '/' in ff):
            self.args.protein_ff = Path(ff.rstrip('/')).stem
            put_log(f'--protein-ff looks like a path, using force field name: {self.args.protein_ff}')
        if not self.args.peptide and self.args.ligand_ff == 'espaloma':
            try:
                from espaloma import get_model
                get_model('0.3.1')
                put_log('espaloma model 0.3.1 OK.')
            except Exception as e:
                put_err(f'espaloma check failed: {e}', _exit=True)
        # find complex file
        if os.path.isdir(self.args.dir):
            complex_paths = get_paths_with_extension(self.args.dir, [], name_substr=self.args.name)
        else:
            put_err(f'dir argument should be a directory: {self.args.dir}', _exit=True)
        if len(complex_paths) == 0:
            put_err(f'no {self.args.name} found in {self.args.dir}', _exit=True)
        complex_path = Path(complex_paths[0]).resolve()
        wdir = complex_path.parent
        put_log(f'processing ABFE prepare for: {complex_path}')

        # sanitize numeric columns (NaN/Inf/huge occupancy or b-factor) before
        # handing the file to PyMOL: its Map code may segfault on such data
        complex_path = _sanitize_pdb_for_pymol(complex_path)

        # extract receptor and ligand chains via pymol
        from pymol import cmd
        import shutil
        cmd.reinitialize()
        cmd.load(str(complex_path), 'complex')
        rec_sel = '(' + ' or '.join([f'chain {c}' for c in self.args.receptor_chains]) + ')'
        if cmd.select('receptor', f'complex and {rec_sel}') == 0:
            put_err(f'receptor chains {self.args.receptor_chains} have zero atoms.', _exit=True)
        if cmd.select('ligand', f'complex and chain {self.args.ligand_chain}') == 0:
            put_err(f'ligand chain {self.args.ligand_chain} has zero atoms.', _exit=True)
        receptor_pdb = wdir / 'receptor.pdb'
        ligand_pdb = wdir / 'ligand.pdb'
        cmd.save(str(receptor_pdb), 'receptor')
        cmd.save(str(ligand_pdb), 'ligand')
        cmd.reinitialize()
        put_log(f'receptor saved: {receptor_pdb}, ligand saved: {ligand_pdb}')

        # ligand resname fix for small molecule
        if not self.args.peptide:
            ligand_pdb = _fix_ligand_for_abfe(ligand_pdb, 'LIG')

        gmx = Gromacs(working_dir=str(wdir))

        # build MakeInputs
        protein_def = {'conf': str(receptor_pdb), 'ff': {'code': self.args.protein_ff}}
        peptide_def = None
        ligand_def = None
        if self.args.peptide:
            peptide_def = {'conf': str(ligand_pdb), 'ff': {'code': self.args.protein_ff}}
        else:
            ligand_def = {'conf': str(ligand_pdb), 'ff': {'type': self.args.ligand_ff}}

        out_input = wdir / 'input'
        builder = system_builder.MakeInputs(
            protein=protein_def,
            host_name='Protein',
            water_model=self.args.water_model,
            custom_ff_path=self.args.ff_dir,
            hmr_factor=self.args.hmr_factor,
            fix_protein=self.args.fix_protein,
            solv_d=self.args.solv_d,
            solv_bt=self.args.solv_bt,
            solv_ion_conc=self.args.ion_conc,
            builder_dir=wdir / self.args.builder_dir,
            gmx=gmx,
            peptide_definition=peptide_def,
            pdb2gmx_args=self.args.pdb2gmx_args.split(),
            n_term=self.args.n_term,
            c_term=self.args.c_term,
            maxwarn=self.args.maxwarn,
        )
        with builder:
            builder(ligand_definition=ligand_def, out_dir=out_input)

        put_log(f'ABFE prepare completed. Output in: {out_input}')


_str2func = {
    'complex': Complex,
}


def main(sys_args: List[str] = None):
    make_args_and_excute('tools for ABFE (absolute binding free energy) preparation.', _str2func, sys_args)


if __name__ == '__main__':
    main()
