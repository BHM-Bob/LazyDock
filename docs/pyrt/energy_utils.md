# PyRT Energy Utils Documentation

## Overview

The `energy_utils` module provides utilities for calculating interface energies between receptor and ligand using PyRosetta. The core functionality is implemented in the `calcu_interface_energy` function, which calculates the binding energy at the interface of a protein complex.

## Function Reference

### `calcu_interface_energy`

```python
def calcu_interface_energy(pdb: str, receptor_chains: Union[str, List[str]],
                           ligand_chains: Union[str, List[str]], scorefxn_name: str = 'ref2015') -> float:
```

#### Description

Calculate the interface energy between receptor and ligand using PyRosetta's InterfaceAnalyzerMover. The function can accept either a PDB file path or a PDB string content, providing flexibility for different use cases.

#### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb` | `str` | PDB file path or PDB string content |
| `receptor_chains` | `Union[str, List[str]]` | Receptor chains (single chain as string or multiple chains as list) |
| `ligand_chains` | `Union[str, List[str]]` | Ligand chains (single chain as string or multiple chains as list) |
| `scorefxn_name` | `str` | Energy function name (default: 'ref2015') |

#### Returns

- `float`: Interface energy (dG_separated) in kcal/mol

## Available Score Functions

| Score Function | Description | Requirements |
|----------------|-------------|--------------|
| `ref2015` | Default all-atom energy function (recommended for most cases) | None |
| `score12` | Legacy all-atom energy function | None |
| `ref2015_cart` | Cartesian version of ref2015 | None |
| `beta_nov16` | Beta version of the 2016 energy function | Requires `-corrections::beta_nov16 true` flag |
| `talaris2014` | Talaris 2014 energy function | Requires `-corrections::restore_talaris_behavior true` flag |

## Scores Dictionary Metrics

The function populates the pose.scores dictionary with valuable metrics about the protein-protein interface. These metrics provide insights into the structural and energetic properties of the interface.

### Energy Metrics

| Metric | Description | Unit |
|--------|-------------|------|
| `dG_separated` | Interface energy calculated from separated structures | kcal/mol |
| `dG_cross` | Cross energy term between the two proteins | kcal/mol |
| `dG_cross/dSASAx100` | Normalized cross energy per 100 Å² SASA | kcal/mol / 100 Å² |
| `dG_separated/dSASAx100` | Normalized interface energy per 100 Å² SASA | kcal/mol / 100 Å² |
| `per_residue_energy_int` | Average energy per interface residue | kcal/mol/residue |

### Surface Area Metrics

| Metric | Description | Unit |
|--------|-------------|------|
| `dSASA_int` | Total solvent accessible surface area change at interface | Å² |
| `dSASA_hphobic` | Hydrophobic SASA change at interface | Å² |
| `dSASA_polar` | Polar SASA change at interface | Å² |

### Hydrogen Bond Metrics

| Metric | Description | Unit |
|--------|-------------|------|
| `hbonds_int` | Number of hydrogen bonds at interface | Count |
| `delta_unsatHbonds` | Change in unsatisfied hydrogen bonds upon complex formation | Count |
| `hbond_E_fraction` | Fraction of total energy from hydrogen bonds | Ratio |

### Structure Quality Metrics

| Metric | Description | Unit |
|--------|-------------|------|
| `sc_value` | Shape complementarity score (higher is better) | 0-1 |
| `packstat` | Packing statistics (requires additional setup) | 0-1 |

### Residue Count Metrics

| Metric | Description | Unit |
|--------|-------------|------|
| `nres_all` | Total number of residues in complex | Count |
| `nres_int` | Number of interface residues | Count |

## Usage Examples

### Basic Usage (from file path)

```python
from lazydock.pyrt.energy_utils import calcu_interface_energy

# Calculate interface energy using default score function
energy = calcu_interface_energy('complex.pdb', 'A', 'B')
print(f"Interface energy: {energy:.4f} kcal/mol")
```

### Using PDB String

```python
from lazydock.pyrt.energy_utils import calcu_interface_energy

# Get PDB content as string (e.g., from API or other source)
pdb_string = get_pdb_from_api('1abc')

# Calculate interface energy from string
energy = calcu_interface_energy(pdb_string, 'A', 'B')
print(f"Interface energy: {energy:.4f} kcal/mol")
```

### Using Different Score Function

```python
from lazydock.pyrt.energy_utils import calcu_interface_energy

# Calculate interface energy using score12
energy = calcu_interface_energy('complex.pdb', 'A', 'B', scorefxn_name='score12')
print(f"Interface energy (score12): {energy:.4f} kcal/mol")
```

### Using Multiple Chains

```python
from lazydock.pyrt.energy_utils import calcu_interface_energy

# Calculate interface energy with multiple receptor chains
energy = calcu_interface_energy('complex.pdb', ['A', 'B'], 'C')
print(f"Interface energy: {energy:.4f} kcal/mol")
```

## Advanced Usage: Accessing Additional Scores

While the function returns only the `dG_separated` value, you can access all metrics by directly using PyRosetta's InterfaceAnalyzerMover. Here's an example:

```python
import pyrosetta
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

# Initialize PyRosetta
pyrosetta.init(silent=True)

# Create pose
pose = pyrosetta.pose_from_pdb('complex.pdb')

# Create InterfaceAnalyzerMover
interface = 'B_A'  # ligand_chain_receptor_chain
mover = InterfaceAnalyzerMover(interface)
mover.set_pack_separated(True)
mover.apply(pose)

# Access all scores
for score_name, score_value in pose.scores.items():
    print(f"{score_name}: {score_value}")
```

## Best Practices

1. **Score Function Selection**:
   - Use `ref2015` for most cases (default)
   - Use `score12` for legacy compatibility
   - Use `ref2015_cart` for Cartesian space calculations

2. **Chain Specification**:
   - Use single character strings for individual chains
   - Use lists for multiple chains
   - Ensure chain identifiers match those in the PDB file

3. **Performance Considerations**:
   - For large complexes, the calculation may take several seconds
   - Reuse PyRosetta initialization across multiple calls
   - Consider using parallel processing for batch calculations

4. **Interpreting Results**:
   - Negative `dG_separated` values indicate favorable binding
   - Higher `sc_value` indicates better shape complementarity
   - Higher `hbonds_int` indicates more hydrogen bonds at the interface
   - Normalized energies (per 100 Å² SASA) allow comparison between different interfaces

## Error Handling

The function handles common errors gracefully:

- Invalid PDB paths will raise a `RuntimeError`
- Invalid chain identifiers will be handled by PyRosetta
- Unsupported energy functions will raise appropriate exceptions

## Dependencies

- `pyrosetta`: PyRosetta is required for the energy calculations
- `typing`: For type annotations

## Installation

```bash
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
```

## See Also

- [PyRosetta Documentation](https://www.pyrosetta.org/)
- [InterfaceAnalyzerMover Documentation](https://www.pyrosetta.org/documentation/html/classprotocols_1_1analysis_1_1InterfaceAnalyzerMover.html)
- [Rosetta Energy Functions](https://www.ncbi.nlm.nih.gov/pubmed/28430426)

## License

This module is part of the LazyDock package and follows its licensing terms.
