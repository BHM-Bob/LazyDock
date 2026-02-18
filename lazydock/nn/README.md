# LazyDock Neural Network Models

This sub-module contains deep learning models for AIDD (AI-driven Drug Discovery).

## EHIGN (Enhanced Heterogeneous Interaction Graph Network)

EHIGN is a deep learning model for protein-ligand interaction prediction using heterogeneous graph neural networks.

### Overview

EHIGN constructs a heterogeneous graph representation of protein-ligand complexes, where:
- **Ligand nodes** represent atoms in the ligand molecule
- **Pocket nodes** represent atoms in the protein binding pocket
- **Edges** capture both intramolecular (within ligand/pocket) and intermolecular (ligand-pocket) interactions

The model uses PyMOL for automatic pocket generation from receptor structures and RDKit for molecular feature extraction.

### Installation Requirements

```bash
# Required dependencies
pip install torch dgl rdkit-pypi pymol-open-source scipy networkx numpy
```

### Module Structure

```
lazydock/nn/EHIGN/
├── __init__.py          # Package exports and API interface
├── api.py               # High-level API and model management
├── data.py              # Graph construction and data processing
├── model.py             # EHIGN neural network architecture
└── utils.py             # Utility functions for geometric calculations
```

### Core Functions

#### 1. Graph Construction Functions

**From PDB Files:**
- `pdb_path_to_graph()` - Build single graph from PDB file paths
- `pdb_path_list_to_batch()` - Build batch graph from PDB file lists

**From PDB Strings:**
- `pdbstr_to_graph()` - Build single graph from PDB strings
- `pdbstr_list_to_batch()` - Build batch graph from PDB string lists

**From Complex Structures:**
- `complex_to_graph()` - Build graph from complex PDB with chain specifications

**Pocket Generation:**
- `_generate_pocket_from_receptor()` - Generate pocket from receptor using PyMOL

#### 2. Model Management Functions

- `build_model()` - Create EHIGN model instance
- `load_checkpoint()` - Load pre-trained model weights
- `evaluate()` - Run inference on batch graphs

#### 3. EHIGN Class API

The `EHIGN` class provides a unified interface for model inference:

```python
from lazydock.nn.EHIGN import EHIGN

# Initialize model
ehign = EHIGN(checkpoint_path='path/to/model.pt', device='cuda')

# Build batch from file paths
batch_graph, sample_ids = ehign.build_batch_from_paths(
    ligand_path_list=['lig1.pdb', 'lig2.pdb'],
    pocket_path_list=['pocket1.pdb', 'pocket2.pdb']
)

# Or use receptor to automatically generate pocket
batch_graph, sample_ids = ehign.build_batch_from_paths(
    ligand_path_list=['lig1.pdb', 'lig2.pdb'],
    receptor_path_list=['rec1.pdb', 'rec2.pdb']
)

# Run prediction
scores_lp, scores_pl = ehign.predict_from_paths(
    ligand_path_list=['lig1.pdb', 'lig2.pdb'],
    receptor_path_list=['rec1.pdb', 'rec2.pdb']
)

# Or use complex structures with chain specifications
scores_lp, scores_pl = ehign.predict_from_complexes(
    complex_path_list=['complex1.pdb', 'complex2.pdb'],
    ligand_chain_list=['B', 'B'],
    receptor_chain_list=['A', 'A']
)
```

### Usage Examples

#### Basic Usage with File Paths

```python
from lazydock.nn.EHIGN import EHIGN

# Initialize model
ehign = EHIGN(checkpoint_path='model.pt')

# Predict from file paths
scores_lp, scores_pl = ehign.predict_from_paths(
    ligand_path_list=['data/ligand1.pdb', 'data/ligand2.pdb'],
    receptor_path_list=['data/receptor1.pdb', 'data/receptor2.pdb'],
    dis_threshold=5.0
)

print(f"Ligand→Pocket scores: {scores_lp}")
print(f"Pocket→Ligand scores: {scores_pl}")
```

#### Using PDB Strings

```python
from lazydock.nn.EHIGN import EHIGN

# Initialize model
ehign = EHIGN()

# Load pre-trained weights
ehign.load_model('model.pt')

# Predict from PDB strings
scores_lp, scores_pl = ehign.predict_from_strs(
    ligand_pdbstr_list=[ligand_str1, ligand_str2],
    receptor_pdbstr_list=[receptor_str1, receptor_str2]
)
```

#### Using Complex Structures

```python
from lazydock.nn.EHIGN import EHIGN

# Initialize model
ehign = EHIGN(checkpoint_path='model.pt')

# Predict from complex PDB files with chain specifications
scores_lp, scores_pl = ehign.predict_from_complexes(
    complex_path_list=['complex1.pdb', 'complex2.pdb'],
    ligand_chain_list=['B', 'B'],      # Ligand chain IDs
    receptor_chain_list=['A', 'A'],    # Receptor chain IDs
    dis_threshold=5.0
)

print(f"Affinity scores: {scores_lp}")
```

#### Manual Graph Construction

```python
from lazydock.nn.EHIGN import pdb_path_to_graph, graphs_to_batch, evaluate
from lazydock.nn.EHIGN import build_model, load_checkpoint

# Build individual graphs
graph1 = pdb_path_to_graph('lig1.pdb', receptor_path='rec1.pdb')
graph2 = pdb_path_to_graph('lig2.pdb', receptor_path='rec2.pdb')

# Create batch
batch_graph = graphs_to_batch([graph1, graph2])

# Build and load model
model = build_model()
load_checkpoint(model, 'model.pt')

# Evaluate
scores_lp, scores_pl = evaluate(model, batch_graph)
```

### Model Architecture

EHIGN uses a heterogeneous graph neural network with:

- **Node Features**: 35-dimensional atomic features
- **Edge Features**: 17-dimensional intramolecular edge features
- **Inter-edge Features**: 11-dimensional intermolecular edge features
- **Graph Convolution Layers**: Custom CIGConv and NIGConv layers
- **Output**: Bidirectional affinity scores (ligand→pocket and pocket→ligand)

### Configuration Parameters

```python
# Default model configuration
model_config = {
    'node_feat_size': 35,      # Atomic feature dimension
    'edge_feat_size': 17,      # Intramolecular edge feature dimension
    'inter_edge_feat_size': 11, # Intermolecular edge feature dimension
    'hidden_feat_size': 256,   # Hidden layer dimension
    'layer_num': 3             # Number of graph convolution layers
}
```

### Pocket Generation

The module uses PyMOL for automatic pocket generation:

1. Loads receptor and ligand structures
2. Removes water molecules and inorganic compounds
3. Selects residues within specified distance (default: 5.0Å) from ligand
4. Generates pocket PDB structure

### Error Handling

The module includes robust error handling:
- Invalid PDB files are skipped with warnings
- Failed pocket generation returns None
- Batch processing continues even if individual samples fail
- Comprehensive error messages for debugging


### Citation

If you use EHIGN in your research, please cite their original paper:

```
Yang, Z., Zhong, W., Lv, Q., Dong, T., Chen, G., and Chen, C.Y.-C. (2024). Interaction-Based Inductive Bias in Graph Neural Networks: Enhancing Protein-Ligand Binding Affinity Predictions From 3D Structures. IEEE Transactions on Pattern Analysis and Machine Intelligence, 1-18. 10.1109/TPAMI.2024.3400515.
```

### License

This module is part of the LazyDock package. See the main LICENSE file for details.