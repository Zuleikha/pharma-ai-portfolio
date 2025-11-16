# ✅ AlphaFold Target Pipeline — README

This module contains the complete pipeline used to:

✔ Run **AlphaFold/ColabFold**  
✔ Predict protein structure (DHFR in this case)  
✔ Visualize key structural features  
✔ Generate images for analysis & portfolio use  

It is part of the **pharma-ai-portfolio** project.

## 📁 Folder Structure
```
alphafold_target_pipeline/
│
├── notebooks/
│   ├── 01_run_alphafold_colab.ipynb
│   └── 02_structure_analysis.ipynb
│
├── data/
│   └── protein_sequence.fasta
│
├── output/
│   └── structures/
│       ├── DHFR_rank1.pdb
│       ├── DHFR_rank2.pdb
│       ├── DHFR_rank3.pdb
│       ├── DHFR_rank4.pdb
│       ├── DHFR_rank5.pdb
│
└── images/
    ├── dhfr_cartoon_secondary_structure.png
    ├── dhfr_cartoon_greyscale.png
    ├── dhfr_all_atom_sticks.png
    ├── dhfr_active_site_highlighted.png
    ├── dhfr_surface_exposed_residues.png
    └── dhfr_active_site_surface_view.png
```

---

# 📓 Notebooks Overview

## 1️⃣ 01_run_alphafold_colab.ipynb
Runs AlphaFold using **ColabFold**, generating:

- Multiple PDB models  
- MSA alignment via MMseqs2  
- Ranked models saved to `output/structures/`

This notebook is designed to be run in **Google Colab** for GPU acceleration.

---

## 2️⃣ 02_structure_analysis.ipynb
Loads AlphaFold-generated structures and creates publication-ready visualizations using **py3Dmol**.

Outputs include:

- 🌈 Rainbow-colored cartoon  
- ⚪ Greyscale cartoon  
- 🧱 Full-atom stick view  
- 🔴 Active-site residue highlights  
- ☁️ Surface representations  

Images saved inside **images/**.

---

# 🖼️ Image Gallery

## 🌈 Rainbow Colored Secondary Structure
![DHFR Cartoon Secondary Structure](dhfr_cartoon_secondary_structure.png)

## ⚪ Greyscale Cartoon
![DHFR Cartoon Greyscale](dhfr_cartoon_greyscale.png)

## 🧱 Full Atom Stick Representation
![DHFR All Atom Sticks](dhfr_all_atom_sticks.png)

## 🔴 Active Site Highlighted
![DHFR Active Site Highlighted](dhfr_active_site_highlighted.png)

## ☁️ Surface with Exposed Residues
![DHFR Surface Exposed Residues](dhfr_surface_exposed_residues.png)

## 🔵 Active Site Surface View
![DHFR Active Site Surface View](dhfr_active_site_surface_view.png)

---

# 🧬 DHFR Summary

**Dihydrofolate Reductase (DHFR)** is a key enzyme in:

- DNA synthesis  
- Folate metabolism  
- Cell replication  

It is a major drug target for:

- 🦠 antibiotics (e.g., trimethoprim)  
- 🎗️ cancer therapy (e.g., methotrexate)  

This makes it a perfect molecule for your drug-design portfolio.

---

# 🚀 Next Steps

This pipeline supports upcoming project phases:

1. **Binding pocket detection** (fpocket or PyMol)  
2. **Docking setup** (AutoDock Vina or DiffDock)  
3. **Ligand preparation & scoring**  
4. **End-to-end portfolio integration**

---
