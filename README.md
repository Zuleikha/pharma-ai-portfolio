# 🧬 AI-Powered Drug Discovery Portfolio  
A collection of machine-learning and structure-based drug discovery projects, including AlphaFold structure prediction, binding site analysis, ligand preparation, docking pipelines, and molecular visualization.

---

# 📁 Project Structure

```
pharma-ai-portfolio/
│
├── alphafold_target_pipeline/
│   ├── data/
│   ├── images/
│   │   └── structures/
│   ├── output/
│   ├── notebooks/
│   └── scripts
│
└── README.md  ← You are here
```

---

# 🧫 DHFR (Dihydrofolate Reductase) Structural Analysis

These images were generated using AlphaFold, PyMOL, and fpocket to visualize the predicted structure of DHFR and explore its binding pocket.

All images are stored here:  
`alphafold_target_pipeline/images/structures/`

---

## 🔍 DHFR Structural Visualizations

### 🔴 Active Site Highlighted  
<img src="alphafold_target_pipeline/images/structures/dhfr_active_site_highlighted.png" width="500"/>

---

### 🔴 Active Site Surface View  
<img src="alphafold_target_pipeline/images/structures/dhfr_active_site_surface_view.png" width="500"/>

---

### 🧱 All-Atom Stick Model  
<img src="alphafold_target_pipeline/images/structures/dhfr_all_atom_sticks.png" width="500"/>

---

### 🎨 Cartoon (Greyscale)  
<img src="alphafold_target_pipeline/images/structures/dhfr_cartoon_greyscale.png" width="500"/>

---

### 🌈 Cartoon by Secondary Structure  
<img src="alphafold_target_pipeline/images/structures/dhfr_cartoon_secondary_structure.png" width="500"/>

---

### 🟡 Surface-Exposed Residues  
<img src="alphafold_target_pipeline/images/structures/dhfr_surface_exposed_residues.png" width="500"/>

---

# 🔬 Docking Results (AutoDock Vina)

Ligands were generated with:

- RDKit (SMILES → 3D conformers)  
- OpenBabel (format conversion)  
- MGLTools `prepare_ligand4.py`  
- AutoDock Vina 1.2.7  

Docking outputs stored in:

```
alphafold_target_pipeline/output/docking/
```

---

## 🧡 Trimethoprim Docking Pose  
Affinity ≈ **−7.4 kcal/mol**

<img src="alphafold_target_pipeline/images/structures/dhfr_trimethoprim_docking.png" width="500"/>

---

## 🧡 Pyrimethamine Docking Pose  
Affinity ≈ **−8.0 kcal/mol**

<img src="alphafold_target_pipeline/images/structures/dhfr_pyrimethamine_docking.png" width="500"/>

---

# 🚀 Summary

This project demonstrates a complete structural biology & docking workflow:

- Protein structure prediction (AlphaFold)  
- Binding pocket detection (fpocket, PyMOL)  
- Ligand preparation (RDKit → OpenBabel → MGLTools)  
- Docking (AutoDock Vina)  
- Visualization and reporting (PyMOL)

Perfect for an ML + drug discovery portfolio.

---

# 📌 Next Steps
- Add DiffDock comparisons  
- Add scoring/affinity ranking  
- Add ligand fragmentation & linking  
- Add interactive 3D visualizations (NGL viewer)

---

# # 🔗 Contact / Research Profiles

- **GitHub:** https://github.com/Zuleikha/pharma-ai-portfolio
- **LinkedIn:** https://www.linkedin.com/in/zuleikha-k-45264a36a/
- **Email:** zuleikhak@gmail.com

