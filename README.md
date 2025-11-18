# 🧬 Pharma‑AI Portfolio  
A polished, professional portfolio showcasing practical skills across **computational chemistry**, **structural bioinformatics**, **AlphaFold**, **docking**, **binding‑site analysis**, and **AI‑driven drug discovery**.  
Built end‑to‑end as a demonstration of real‑world drug discovery workflows.

---

# 🌟 Overview

This portfolio documents a hands‑on journey from software engineering into **computational drug discovery** and **pharma AI/ML**, combining:

- Protein structure prediction (AlphaFold / ColabFold)  
- Structural analysis & visualization  
- Ligand preparation & docking (RDKit + AutoDock Vina)  
- Fragment‑based design  
- Molecular property analysis  
- Data‑driven decision making  
- Future extensions: pocket detection, scoring, ADMET, and ML models  

This repository is both a **learning artifact** and a **working toolkit** for real drug‑discovery pipelines.

---

# 🗂️ Project Structure (Updated)

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
├─ notebooks/                      # Interactive cheminformatics & structural workflows
│
├─ src/                            # Core reusable Python modules (RDKit, docking, FBDD)
│
├─ images/                         # Portfolio figures (molecules, structures, conformers)
│
├─ output/                         # Output tables (e.g., molecular properties)
│
├─ LEARNING_PATH.md                # My structured journey into pharma/AI
├─ FUTURE_WORK.md                  # Detailed project roadmap
└─ README.md                       # <- You are here
```

---

# 🧩 Included Pipelines & Workflows

## **1️⃣ AlphaFold Target Pipeline (DHFR)**  
A complete structural‑biology workflow:

- Input FASTA → AlphaFold structure prediction  
- Structural visualization (py3Dmol + PyMOL)  
- Active site characterization  
- Ligand prep (RDKit → PDBQT)  
- Docking with AutoDock Vina  
- Interaction & binding pose analysis  

**Full documentation →** `alphafold_target_pipeline/README.md`

---

## **2️⃣ Fragment‑Based Drug Design (FBDD)**

Explores:

- Fragmentation  
- Recombination  
- Scoring & filtering  
- Conformer generation  
- Visualization  

Notebook → `notebooks/fragment_based_drug_design.ipynb`

---

## **3️⃣ Molecular Property Analysis**

Computes key drug‑likeness descriptors:

- LogP / MW / TPSA  
- HBD / HBA  
- Rule‑of‑5 flags  
- Molecular visualization  

Notebook → `notebooks/molecular_property_analysis.ipynb`

---

## **4️⃣ Docking Preparation**

Automates:

- Molecule sanitization  
- Geometry optimization  
- Protonation  
- Conversion to PDBQT  

Notebook → `notebooks/Molecular Docking Preparation.ipynb`

---

# 🖼️ Highlight Images

## Molecular Examples
![Aspirin](images/aspirin.png)
![Caffeine](images/best_2d.png)
![Three Drugs](images/three_drugs.png)

---

# 🖼️ DHFR Structural Gallery

These images were generated using your AlphaFold model + PyMOL:

![DHFR SS](alphafold_target_pipeline/images/structures/protein/dhfr_cartoon_secondary_structure.png)
![DHFR Gray](alphafold_target_pipeline/images/structures/protein/dhfr_cartoon_greyscale.png)
![DHFR Sticks](alphafold_target_pipeline/images/structures/protein/dhfr_all_atom_sticks.png)
![DHFR Active Site](alphafold_target_pipeline/images/structures/protein/dhfr_active_site_highlighted.png)
![DHFR Surface](alphafold_target_pipeline/images/structures/protein/dhfr_surface_exposed_residues.png)
![DHFR Surface 2](alphafold_target_pipeline/images/structures/protein/dhfr_active_site_surface_view.png)

---

# 🧬 Why DHFR?

**Dihydrofolate Reductase** is a clinically validated target involved in:

- DNA synthesis  
- Cell proliferation  
- Folate metabolism  

Drug classes targeting DHFR include:

- Antibiotics (trimethoprim, pyrimethamine)  
- Anticancer agents (methotrexate)

This makes DHFR ideal for showcasing structural modeling, docking, and ligand–protein interaction analysis.

---

# 🚀 Installation

```bash
git clone https://github.com/yourusername/pharma-ai-portfolio.git
cd pharma-ai-portfolio
pip install -r requirements.txt
```

---

# 📈 Roadmap

See **FUTURE_WORK.md** for the complete multi‑phase roadmap, including:

- Binding pocket detection (fpocket, ML models)  
- Pose scoring & clustering  
- Fragment growing & lead optimization  
- ADMET ML pipeline  
- End‑to‑end hit discovery workflow  

---

# 🤝 Contributing

Suggestions and contributions are welcome, especially improvements to documentation, workflows, and visualizations.

---

# 📄 License

MIT License.
