# 🧬 AlphaFold Target Pipeline — DHFR  
A complete, portfolio‑ready workflow showing how to:

✔ Run **AlphaFold/ColabFold**  
✔ Generate & analyze protein structures  
✔ Create publication‑quality images  
✔ Prepare ligands  
✔ Perform docking  
✔ Visualize binding interactions  

---

# 📁 Folder Structure

```text
alphafold_target_pipeline/
│
├── notebooks/
│   ├── 01_run_alphafold_colab.ipynb
│   └── 02_structure_analysis.ipynb
│
├── data/
│   ├── protein_sequence.fasta
│   ├── ligands/
│   └── protein/
│
├── images/
│   ├── protein_ligand_docking.png
│   └── structures/
│       ├── docking/
│       ├── interactions/
│       └── protein/
│
├── output/
│   ├── structures/
│   ├── docking/
│   └── vina/
│
└── src/
```

---

# 📓 Notebook Workflows

## **1️⃣ 01_run_alphafold_colab.ipynb**
Runs **ColabFold** to generate:

- 5 ranked PDB models  
- MSA via MMseqs2  
- PAE & pLDDT confidence plots  
- All saved to:  
  `output/structures/`

---

## **2️⃣ 02_structure_analysis.ipynb**
Visualizes AlphaFold output using **py3Dmol + PyMOL**:

Produces:

- Rainbow cartoon  
- Greyscale cartoon  
- Full‑atom sticks  
- Active‑site residues  
- Surface views  

Saved to:  
`images/structures/protein/`

---

# 🖼️ Protein Structure Gallery

![SS](images/structures/protein/dhfr_cartoon_secondary_structure.png)
![Gray](images/structures/protein/dhfr_cartoon_greyscale.png)
![Sticks](images/structures/protein/dhfr_all_atom_sticks.png)
![Active](images/structures/protein/dhfr_active_site_highlighted.png)
![Surface](images/structures/protein/dhfr_surface_exposed_residues.png)
![Surface2](images/structures/protein/dhfr_active_site_surface_view.png)

---

# 🔬 Docking Workflow

Ligands prepared in:

```
data/ligands/
```

Protein PDB & PDBQT:

```
data/protein/
```

Docking results:

```
output/docking/
output/vina/
images/structures/docking/
images/structures/interactions/
```

Interaction diagrams:

![Docking](images/structures/docking/dhfr_trimethoprim_docking.png)
![Docking2](images/structures/docking/dhfr_pyrimethamine_docking.png)
![Int](images/structures/interactions/dhfr_trimethoprim_interaction.png)
![Int2](images/structures/interactions/dhfr_pyrimethamine_interactions.png)

---

# 🧬 DHFR Summary

DHFR plays a central role in:

- DNA synthesis  
- Folate metabolism  
- Cell growth  

Targeted by:

- Methotrexate  
- Trimethoprim  
- Pyrimethamine  

This makes it an excellent molecule for demonstrating structure‑guided drug design.

---

# 🚀 Future Extensions

- Binding pocket mapping (fpocket / PyMOL)  
- Fragment elaboration  
- Interaction fingerprinting  
- Pose scoring  
- ADMET predictions  
- ML‑guided ligand ranking  

---

# 🧰 Notes

This pipeline intentionally mirrors workflows used in:

- Academic computational chemistry labs  
- Pharma / biotech structural biology teams  
- Structure‑based drug‑design pipelines  

Making it ideal for portfolio, interview, and learning use.
