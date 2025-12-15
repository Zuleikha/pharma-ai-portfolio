# Pharma AI Portfolio

A focused portfolio demonstrating how machine learning, computational chemistry, and structural bioinformatics are applied in early stage drug discovery.

This repository documents a practical transition from software engineering into pharma focused AI and computational drug discovery, with emphasis on clear biological context, reproducible workflows, and decision making aligned with real industry pipelines.

---

## Portfolio Goal

This portfolio demonstrates how AI and machine learning integrate with structural biology and cheminformatics to support target analysis, ligand design, and hit identification in drug discovery. Each project is designed so a pharma recruiter can understand the problem, biology, method, and outcome within minutes.

---

## Core Projects

This portfolio is intentionally limited to three core areas:

1. Structure based target analysis and docking  
2. Fragment based drug design and molecular optimisation  
3. Data driven molecular property analysis and filtering  

Together, these represent key steps in early discovery workflows.

---

## Overview of Capabilities

The work in this repository covers:

- Protein structure prediction using AlphaFold and ColabFold  
- Structural inspection and binding site analysis  
- Ligand preparation and docking using RDKit and AutoDock Vina  
- Fragment based drug design workflows  
- Molecular property calculation and drug likeness filtering  
- End to end pipelines from biological input to actionable outputs  

This repository functions as both a learning record and a working technical toolkit.

---

## Project Structure

```text
pharma-ai-portfolio/
│
├─ alphafold_target_pipeline/
│   ├─ data/
│   ├─ notebooks/
│   ├─ images/
│   ├─ output/
│   └─ src/
│
├─ notebooks/
├─ src/
├─ images/
├─ output/
│
├─ LEARNING_PATH.md
├─ FUTURE_WORK.md
└─ README.md
```

---

## Included Workflows

### AlphaFold Target Pipeline (DHFR)

This workflow demonstrates a complete structure based drug discovery pipeline using Dihydrofolate Reductase as a validated biological target.

Steps covered:

- Selection of DHFR as a clinically relevant enzyme target  
- Retrieval and preparation of the protein FASTA sequence  
- Protein structure prediction using AlphaFold or ColabFold  
- Structural inspection and quality assessment  
- Identification of the active site and binding pocket  
- Ligand preparation using RDKit  
- Docking using AutoDock Vina  
- Binding pose and interaction analysis  

---

### Fragment Based Drug Design

- Fragment library generation  
- Fragment recombination and growing  
- Scoring and filtering  
- Conformer generation  
- Visualisation  

---

### Molecular Property Analysis

- Molecular weight  
- LogP  
- TPSA  
- Hydrogen bond donors and acceptors  
- Rule of five filtering  

---

### Docking Preparation Pipeline

- Molecule sanitisation  
- Protonation handling  
- Geometry optimisation  
- Export to PDBQT  

---

## Why DHFR

Dihydrofolate Reductase is a clinically validated drug target involved in DNA synthesis, cell proliferation, and folate metabolism. It is targeted by established therapeutics in infectious disease and oncology, making it a realistic system for demonstrating structure based drug discovery techniques.

---

## Installation

```bash
git clone https://github.com/Zuleikha/pharma-ai-portfolio.git
cd pharma-ai-portfolio
pip install -r requirements.txt
```

Optional conda environment:

```bash
conda env create -f environment.yml
conda activate pharma-ai-env
```

---

## Roadmap

Planned extensions include binding pocket detection, pose scoring and clustering, fragment growing and optimisation, ADMET filtering and prediction, and machine learning models for compound prioritisation.

---

## License

MIT License
