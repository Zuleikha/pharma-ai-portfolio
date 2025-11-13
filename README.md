# Pharma AI Portfolio

A comprehensive collection of pharmaceutical and drug discovery tools using RDKit and machine learning.

## Overview

This repository contains tools and utilities for computational drug discovery, molecular analysis, and pharmaceutical research. It leverages RDKit for cheminformatics and includes various modules for drug design and molecular property prediction.

## Features

- **Molecular Property Analysis**: Calculate and analyze molecular descriptors and drug-likeness properties
- **Fragment-Based Drug Design**: Tools for fragment library generation and optimization
- **Molecular Docking Preparation**: Prepare molecules for docking studies with proper 3D conformations

## Repository Structure

```
pharma-ai-portfolio/
├── 01-rdkit-basics/       # Basic RDKit tutorials and examples
├── src/                   # Source code modules
│   ├── molecular_property_analyzer.py
│   ├── fragment_based_drug_design.py
│   └── molecular_docking_prep.py
├── images/                # Images and visualizations
└── README.md             # This file
```

## Installation

### Prerequisites

- Python 3.7+
- RDKit
- NumPy
- Pandas (optional, for data handling)

### Install RDKit

```bash
# Using conda (recommended)
conda install -c conda-forge rdkit

# Using pip
pip install rdkit
```

## Usage

### Molecular Property Analyzer

```python
from src.molecular_property_analyzer import MolecularPropertyAnalyzer

analyzer = MolecularPropertyAnalyzer()
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
properties = analyzer.analyze_molecule(smiles)
print(properties)
```

### Fragment-Based Drug Design

```python
from src.fragment_based_drug_design import FragmentBasedDesign

designer = FragmentBasedDesign()
fragments = designer.generate_fragments("CCO")
print(f"Generated {len(fragments)} fragments")
```

### Molecular Docking Preparation

```python
from src.molecular_docking_prep import MolecularDockingPrep

prep = MolecularDockingPrep()
mol = prep.prepare_ligand("CC(C)Cc1ccc(cc1)C(C)C(=O)O")  # Ibuprofen
prep.write_mol_file(mol, "ibuprofen_3d.mol")
```

## Modules

### molecular_property_analyzer.py
Calculates molecular properties including:
- Molecular weight
- LogP (lipophilicity)
- Number of hydrogen bond donors/acceptors
- Topological polar surface area (TPSA)
- Lipinski's Rule of Five compliance
- Quantitative Estimate of Drug-likeness (QED)

### fragment_based_drug_design.py
Tools for fragment-based drug design:
- RECAP fragmentation
- BRICS fragmentation
- Fragment library generation
- Scaffold extraction
- Fragment filtering and scoring

### molecular_docking_prep.py
Prepares molecules for molecular docking:
- 3D conformation generation
- Energy minimization
- Protonation state assignment
- File format conversion (PDB, MOL2, SDF)
- Multiple conformer generation

## Examples

See the `01-rdkit-basics/` directory for introductory examples and tutorials.

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This project is open source and available for educational and research purposes.

## References

- RDKit: Open-source cheminformatics - https://www.rdkit.org/
- Lipinski, C. A. et al. (2001). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings.
- Bickerton, G. R. et al. (2012). Quantifying the chemical beauty of drugs. Nature Chemistry, 4(2), 90-98.

## Contact

For questions or suggestions, please open an issue on GitHub.
