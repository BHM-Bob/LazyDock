<!--
 * @Date: 2025-06-15
 * @LastEditors: BHM-Bob
 * @Description: prepare_amber script documentation
-->

# prepare_amber

The `prepare_amber` script prepares molecular systems for AMBER molecular dynamics simulations. It follows the standard AMBER workflow for preparing protein, ligand, and complex systems with proper force field parameters.

## Pipeline Overview

The script implements a complete preparation pipeline that:

1. **Protein Preparation**: Cleans PDB files, adds missing atoms, and assigns force field parameters
2. **Ligand Parameterization**: Generates GAFF/GAFF2 force field parameters using antechamber
3. **Complex Assembly**: Combines receptor and ligand while preserving original coordinates
4. **Solvation**: Adds water box and neutralizes the system with ions

## Commands

### protein

Prepares a single protein for AMBER MD simulations.

#### Parameters

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | Directory containing protein files |
| `-n` | `--protein-name` | str | None | Protein filename (e.g., `protein.pdb`) |
| | `--protein-ff` | str | `ff19SB` | Protein force field (ff19SB, ff14SB) |
| | `--water-model` | str | `opc` | Water model (opc, tip3p, tip4pew) |
| | `--box-type` | str | `oct` | Box type (oct, box) |
| | `--box-size` | float | `12.0` | Box size in Angstroms |
| | `--neutralize` | flag | `True` | Add ions to neutralize system |
| | `--ion-conc` | float | `0.0` | Additional ion concentration (M) |

#### Pipeline Steps

1. **STEP 0**: Clean PDB file using `pdb4amber`
   - Removes water molecules (`--dry`)
   - Adds missing hydrogens (`--reduce`)
   - Fixes atom naming and residue numbering

2. **STEP 1**: Generate tleap input script
   - Loads protein force field
   - Loads water model and ion parameters
   - Creates solvated system

3. **STEP 2**: Run tleap to generate topology and coordinates
   - Output: `complex.prmtop`, `complex.inpcrd`

#### Example

```bash
# Basic usage
lazydock-cli prepare-amber protein -d /path/to/protein -n protein.pdb

# Custom force field and water model
lazydock-cli prepare-amber protein -d /path/to/protein -n protein.pdb \
    --protein-ff ff14SB --water-model tip3p

# Larger box with ion concentration
lazydock-cli prepare-amber protein -d /path/to/protein -n protein.pdb \
    --box-size 15.0 --ion-conc 0.15
```

---

### ligand

Prepares a single ligand for AMBER MD simulations by generating GAFF/GAFF2 force field parameters.

#### Parameters

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | Directory containing ligand files |
| `-n` | `--ligand-name` | str | None | Ligand filename (e.g., `ligand.pdb`) |
| | `--ligand-ff` | str | `gaff2` | Ligand force field (gaff, gaff2) |
| | `--ligand-charge-method` | str | `bcc` | Charge method (bcc, resp) |
| | `--ligand-net-charge` | int | `0` | Net charge of ligand |
| | `--ligand-residue-name` | str | `LIG` | Residue name for ligand |

#### Pipeline Steps

1. **STEP 0**: Convert ligand to MOL2 format using Open Babel

2. **STEP 1**: Generate GAFF parameters using antechamber
   - Generates GAFF-format MOL2 file (`*_gaff.mol2`)
   - Generates PREP file for compatibility
   - Assigns partial charges using specified method

3. **STEP 2**: Generate frcmod file using parmchk2
   - Identifies missing parameters
   - Generates force field modifications

#### Output Files

- `{ligand}_gaff.mol2`: GAFF-format MOL2 with correct atom types
- `{ligand}.prep`: PREP format file
- `{ligand}.frcmod`: Force field modification file

#### Example

```bash
# Basic usage
lazydock-cli prepare-amber ligand -d /path/to/ligand -n ligand.pdb

# Ligand with net charge
lazydock-cli prepare-amber ligand -d /path/to/ligand -n ligand.pdb \
    --ligand-net-charge 1

# Custom residue name
lazydock-cli prepare-amber ligand -d /path/to/ligand -n ligand.pdb \
    --ligand-residue-name MOL
```

---

### complex

Prepares a protein-ligand complex for AMBER MD simulations. This is the main command for preparing docked complexes.

#### Parameters

| Short | Long | Type | Default | Description |
|-------|------|------|---------|-------------|
| `-d` | `--dir` | str | `.` | Directory containing complex files |
| `-n` | `--complex-name` | str | None | Complex filename (e.g., `complex.pdb`) |
| `-rc` | `--receptor-chain-name` | str | None | Receptor chain identifier |
| `-lc` | `--ligand-chain-name` | str | None | Ligand chain identifier |
| | `--protein-ff` | str | `ff19SB` | Protein force field |
| | `--ligand-ff` | str | `gaff2` | Ligand force field |
| | `--water-model` | str | `opc` | Water model |
| | `--box-type` | str | `oct` | Box type |
| | `--box-size` | float | `12.0` | Box size in Angstroms |
| | `--neutralize` | flag | `True` | Neutralize system with ions |
| | `--ion-conc` | float | `0.0` | Additional ion concentration |
| | `--ligand-charge-method` | str | `bcc` | Charge calculation method |
| | `--ligand-net-charge` | int | `0` | Ligand net charge |
| | `--ligand-residue-name` | str | `LIG` | Ligand residue name |

#### Pipeline Steps

1. **STEP 0.1**: Center complex using Open Babel
2. **STEP 0.2**: Align complex with xyz axes using LazyDock

3. **STEP 1**: Extract receptor and ligand
   - Uses chain identifiers to separate components
   - Receptor: chain specified by `--receptor-chain-name`
   - Ligand: chain specified by `--ligand-chain-name`

4. **STEP 2**: Prepare receptor
   - Clean PDB using pdb4amber
   - Remove ligand atoms from receptor
   - Add missing atoms and hydrogens

5. **STEP 3**: Prepare ligand
   - Generate GAFF MOL2 file with correct atom types
   - Generate frcmod file for missing parameters

6. **STEP 4**: Build complex using tleap
   - Load protein force field
   - Load ligand force field (GAFF/GAFF2)
   - Load water model and ion parameters
   - Load receptor PDB
   - Load ligand MOL2 (preserves original coordinates)
   - Combine receptor and ligand
   - Save gas-phase complex

7. **STEP 5**: Solvate and add ions
   - Add water box
   - Neutralize with Na+/Cl-
   - Add additional ions if specified

#### Coordinate Preservation

The script uses the standard AMBER tutorial method to preserve ligand coordinates:

```bash
# Load receptor
receptor = loadpdb receptor.pdb

# Load ligand (using loadmol2 to preserve coordinates)
loadamberparams ligand.frcmod
lig = loadmol2 ligand_gaff.mol2

# Combine (maintains relative positions)
complex = combine {receptor lig}
```

This ensures the ligand maintains its original position relative to the receptor, which is critical for docked poses.

#### Output Files

- `complex.prmtop`: System topology
- `complex.inpcrd`: Initial coordinates
- `complex_solv.pdb`: Solvated system in PDB format
- `complex_gas.pdb`: Gas-phase complex (for verification)
- `tleap_complex.in`: tleap input script (for reference)

#### Example

```bash
# Basic usage with chain names
lazydock-cli prepare-amber complex -d /path/to/complex -n complex.pdb \
    --receptor-chain-name A --ligand-chain-name Z

# Custom force fields and box
lazydock-cli prepare-amber complex -d /path/to/complex -n complex.pdb \
    --receptor-chain-name A --ligand-chain-name Z \
    --protein-ff ff14SB --ligand-ff gaff \
    --water-model tip3p --box-size 15.0

# Charged ligand
lazydock-cli prepare-amber complex -d /path/to/complex -n complex.pdb \
    --receptor-chain-name A --ligand-chain-name Z \
    --ligand-net-charge 1 --ligand-residue-name MOL
```

---

## Force Field Compatibility

### Protein Force Fields

| Force Field | Description | Recommended For |
|-------------|-------------|---------------|
| `ff19SB` | Latest AMBER protein force field with CMAP | General use (default) |
| `ff14SB` | Previous generation protein force field | Legacy compatibility |

### Water Models

| Model | Description | Ion Parameters |
|-------|-------------|----------------|
| `opc` | Optimized Point Charge (recommended) | Li/Merz ions included |
| `tip3p` | Standard TIP3P | Joung-Cheatham ions |
| `tip4pew` | TIP4P-Ew | Joung-Cheatham ions |

### Ligand Force Fields

| Force Field | Description |
|-------------|-------------|
| `gaff2` | General AMBER Force Field 2.0 (recommended) |
| `gaff` | Original GAFF (legacy) |

## Notes

### Input Requirements

- **Protein PDB**: Should have standard residue names, hydrogens can be added by pdb4amber
- **Complex PDB**: Must have distinct chain identifiers for receptor and ligand
- **Ligand**: Should have all hydrogens, proper bond connectivity

### Coordinate Preservation

The script is designed to preserve the original ligand coordinates from the input complex:

1. Extracts ligand from complex PDB
2. Generates GAFF-format MOL2 with antechamber (preserves coordinates)
3. Loads MOL2 using `loadmol2` (not `loadamberprep` which rebuilds coordinates)
4. Combines with receptor using `combine`

This is critical for maintaining docked poses and relative binding orientations.

### Water Model Selection

- **OPC** (default): Recommended for new simulations, better bulk water properties
- **TIP3P**: Use for compatibility with older force fields
- **TIP4P-Ew**: Use for specific systems requiring 4-site water model

### Ion Parameters

Each water model includes compatible ion parameters:
- OPC: Li/Merz 12-6 parameters (`frcmod.ionslm_126_opc`)
- TIP3P: Joung-Cheatham parameters (`frcmod.ionsjc_tip3p`)

These are automatically loaded when sourcing the water model leaprc file.
