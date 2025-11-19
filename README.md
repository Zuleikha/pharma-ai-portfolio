# Pharma-AI Portfolio

A polished portfolio demonstrating practical skills in computational chemistry, structural bioinformatics, AlphaFold modeling, ligand preparation, docking, binding-site analysis, and data-driven drug-discovery workflows.

This repository documents a hands-on transition from software engineering into computational drug discovery and pharma AI/ML.

---

# Overview

This project combines:

- Protein structure prediction (AlphaFold / ColabFold)
- Structural analysis and visualization
- Ligand preparation and docking (RDKit + AutoDock Vina)
- Fragment-based drug design
- Molecular property analysis
- Practical decision-making workflows for hit identification
- Planned extensions: binding-site detection, ADMET, scoring, ML models

It functions as both a learning record and a working toolkit.

---

# Project Structure

```text
pharma-ai-portfolio/
│
├─ alphafold_target_pipeline/      # AlphaFold → analysis → docking workflow (DHFR)
│   ├─ data/
│   ├─ notebooks/
│   ├─ images/
│   ├─ output/
│   └─ src/
│
├─ notebooks/                      # Cheminformatics & structural analysis notebooks
│
├─ src/                            # Reusable Python modules
│
├─ images/                         # Molecule figures, conformers, visuals
│
├─ output/                         # Output tables (e.g., descriptor results)
│
├─ LEARNING_PATH.md
├─ FUTURE_WORK.md
└─ README.md
```

---

# Included Workflows

## 1. AlphaFold Target Pipeline (DHFR)
A full structural-biology and docking workflow:

- FASTA → AlphaFold structure prediction
- Model inspection & visualization
- Active-site analysis
- Ligand preparation and PDBQT conversion
- Docking using AutoDock Vina
- Binding pose and interaction examination

Detailed docs: `alphafold_target_pipeline/README.md`

---

## 2. Fragment-Based Drug Design (FFDD)

- Fragment generation
- Recombination
- Scoring and filtering
- Conformer generation
- Visualization

Notebook: `notebooks/fragment_based_drug_design.ipynb`

---

## 3. Molecular Property Analysis

Computes drug-likeness metrics:

- LogP, MW, TPSA
- HBD/HBA
- Rule-of-5 parameters
- Molecule visualization

Notebook: `notebooks/molecular_property_analysis.ipynb`

---

## 4. Docking Preparation

Automates:

- Molecule sanitization
- Protonation
- Geometry optimization
- Export to PDBQT

Notebook: `notebooks/Molecular Docking Preparation.ipynb`

---

# Example Images

### Molecular Figures
![Aspirin](images/aspirin.png)  
![Caffeine](images/best_2d.png)  
![Three Drugs](images/three_drugs.png)

---

# DHFR Structural Gallery

Images generated from AlphaFold + PyMOL:

![DHFR SS](alphafold_target_pipeline/images/structures/protein/dhfr_cartoon_secondary_structure.png)
![DHFR Gray](alphafold_target_pipeline/images/structures/protein/dhfr_cartoon_greyscale.png)
![DHFR Sticks](alphafold_target_pipeline/images/structures/protein/dhfr_all_atom_sticks.png)
![DHFR Active Site](alphafold_target_pipeline/images/structures/protein/dhfr_active_site_highlighted.png)
![DHFR Surface](alphafold_target_pipeline/images/structures/protein/dhfr_surface_exposed_residues.png)
![DHFR Surface 2](alphafold_target_pipeline/images/structures/protein/dhfr_active_site_surface_view.png)

---

# Why DHFR?

Dihydrofolate Reductase is a clinically validated target central to:

- DNA synthesis
- Cell proliferation
- Folate metabolism

It is targeted by:

- Antibacterial agents (trimethoprim, pyrimethamine)
- Anticancer therapeutics (methotrexate)

This makes it a suitable system for demonstrating structural modeling, docking, and ligand–protein analysis.

---

# Installation

```bash
git clone https://github.com/yourusername/pharma-ai-portfolio.git
cd pharma-ai-portfolio
pip install -r requirements.txt
```

(Optional) Conda environment:

```bash
conda env create -f environment.yml
conda activate pharma-ai-env
```

---

# Roadmap

See `FUTURE_WORK.md` for the full multi-phase roadmap, covering:

- Binding pocket detection
- Pose scoring and clustering
- Fragment growing and optimization
- ADMET filtering and prediction
- Machine-learning expansion

---

# Contributing

Contributions are welcome, especially improvements to documentation, examples, or workflows.

---

# License

MIT License.
