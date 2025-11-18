# 🧬 Pharma-AI Portfolio  
A modular, expanding portfolio demonstrating practical skills in **computational chemistry**, **structural bioinformatics**, and **AI-driven drug discovery**.  
Each section is a clean, self-contained workflow.

## 📚 Learning Journey: Tech to Pharma AI

This portfolio documents my transition from software engineering to pharmaceutical AI/ML, backed by structured learning in drug discovery fundamentals, computational chemistry, and AI applications in pharma.

**[View Complete Learning Path →](LEARNING_PATH.md)**

**Key Areas Covered:**
- Drug discovery process (target ID → clinical trials)
- Computational chemistry tools (RDKit, AlphaFold, AutoDock)
- ADMET prediction & molecular modeling
- AI/ML for drug design and property prediction

---

# 📁 Project Structure

```text
pharma-ai-portfolio/
│
├─ alphafold_target_pipeline/           # AlphaFold + DHFR structure & docking pipeline
│   ├─ data/
│   │   ├─ protein_sequence.fasta
│   │   ├─ ligands/                     # RDKit-generated ligands + PDBQT files
│   │   └─ protein/                     # DHFR PDB + PDBQT
│   ├─ notebooks/
│   │   ├─ 01_run_alphafold_colab.ipynb
│   │   └─ 02_structure_analysis.ipynb
│   ├─ images/
│   │   ├─ protein_ligand_docking.png   # Overview figure
│   │   └─ structures/                  # DHFR structural & docking renders (PNG)
│   ├─ output/
│   │   ├─ docking/                     # AutoDock Vina poses (PDBQT)
│   │   ├─ vina/                        # Legacy Vina tests
│   │   └─ structures/                  # AlphaFold PDBs + exported images
│   └─ src/                             # (future) pocket detection, scoring, utilities
│
├─ notebooks/                           # Small interactive workflows
│   ├─ fragment_based_drug_design.ipynb
│   ├─ molecular_property_analysis.ipynb
│   ├─ drug_analysis_demo.ipynb
│   ├─ Molecular Docking Preparation.ipynb
│
├─ src/                                 # Core reusable python modules
│   ├─ fragment_based_drug_design.py
│   ├─ molecular_docking_prep.py
│   ├─ molecular_property_analyzer.py
│
├─ archive/
│   └─ rdkit_basics/                    # Early RDKit experiments and scripts
│
├─ images/                              # General molecular structure figures
│   ├─ aspirin.png
│   ├─ best_2d.png
│   └─ three_drugs.png
│
├─ output/
│   └─ molecular_analysis_results.csv
│
├─ FUTURE_WORK.md                       # Detailed pipeline expansion roadmap
├─ LEARNING_PATH.md                     # Tech to Pharma AI learning journey
└─ README.md
```

---

# 🧩 Included Pipelines

## **1️⃣ AlphaFold Target Pipeline (DHFR)**

A full structural-biology mini-workflow:

✔ Input FASTA  
✔ Predict 3D structure via ColabFold/AlphaFold  
✔ Visualize via `py3Dmol` and PyMOL  
✔ Highlight the active site & exposed residues  
✔ Prepare protein + ligands for docking  
✔ Run AutoDock Vina and export portfolio-ready figures  

**Outputs:**  
- Structures → `alphafold_target_pipeline/output/structures/`  
- Docking poses → `alphafold_target_pipeline/output/docking/`  
- Images → `alphafold_target_pipeline/images/structures/`

---

## **2️⃣ Fragment-Based Drug Design**

Explores:

- Fragmentation  
- Recombination  
- Scoring  
- Filtering  
- Visualization  

Notebook → `notebooks/fragment_based_drug_design.ipynb`  
Core script → `src/fragment_based_drug_design.py`

---

## **3️⃣ Molecular Property Analysis**

Computes RDKit-based descriptors:

✔ LogP  
✔ Molecular weight  
✔ HBD / HBA  
✔ Polar surface area  
✔ Rule-of-5 metrics  

Notebook → `notebooks/molecular_property_analysis.ipynb`  
Output → `output/molecular_analysis_results.csv`

---

## **4️⃣ Docking Preparation**

Prepares molecules for docking:

✔ Protonation  
✔ Sanitization  
✔ Geometry optimization  
✔ SDF / MOL / PDBQT export  

Notebook → `notebooks/Molecular Docking Preparation.ipynb`  
Script → `src/molecular_docking_prep.py`

---

# 🧪 Molecular Structure Examples

## Drug Molecule Visualizations

![Aspirin Structure](images/aspirin.png)  
**Aspirin (Acetylsalicylic Acid)** — common NSAID showing ester and carboxylic acid functional groups.

![Caffeine Conformer](images/best_2d.png)  
**Caffeine Lowest Energy Conformer** — RDKit conformer analysis.

![Three Common Drugs](images/three_drugs.png)  
**Comparative Drug Structures** — Aspirin, Caffeine, and Ibuprofen side-by-side.

These visualizations demonstrate 2D molecular rendering and conformer analysis using RDKit, core tools in computational drug discovery.

---

# 🖼️ Structural Images (DHFR)

All DHFR structural and docking PNGs are stored here:

```text
alphafold_target_pipeline/images/structures/
```

Current set includes:

- `dhfr_cartoon_secondary_structure.png`  
- `dhfr_cartoon_greyscale.png`  
- `dhfr_all_atom_sticks.png`  
- `dhfr_active_site_highlighted.png`  
- `dhfr_surface_exposed_residues.png`  
- `dhfr_active_site_surface_view.png`  
- `dhfr_trimethoprim_docking.png`  
- `dhfr_pyrimethamine_docking.png`  

These images were generated from the AlphaFold DHFR model using PyMOL and are used throughout the portfolio to illustrate structure, active-site context, and ligand binding.

---

# 🧬 About DHFR (Dihydrofolate Reductase)

DHFR is essential for:

- DNA synthesis  
- Folate metabolism  
- Cell growth  

It is an important drug target for:

- Cancer therapy (e.g., methotrexate)  
- Antimicrobials (e.g., trimethoprim, pyrimethamine)  

This makes DHFR an ideal example for showcasing structural biology and drug-design skills in this portfolio.

---

# 🧱 Installation

```bash
git clone https://github.com/yourusername/pharma-ai-portfolio.git
cd pharma-ai-portfolio
```

Create a virtual environment:

```bash
python -m venv venv
# Windows (Git Bash / PowerShell)
source venv/Scripts/activate
```

Install required libraries:

```bash
pip install -r requirements.txt
```

---

# 🧪 Dependencies

- RDKit  
- NumPy / SciPy  
- Pandas  
- Matplotlib  
- py3Dmol  
- Biopython  
- scikit-learn  

*(AlphaFold itself is run externally in Colab/ColabFold; docking uses a local AutoDock Vina installation.)*

---

# 🚀 Roadmap (Next Steps)

**[View Detailed Pipeline Expansion Plan →](FUTURE_WORK.md)**

Planned next phases to complete the full discovery workflow:

1. Pose analysis & validation (RMSD clustering, interaction fingerprints)  
2. ADMET filtering and property-based triage  
3. Optional molecular dynamics for binding stability  
4. Hit optimization and fragment growing  
5. Final ranking and lead selection  
6. Portfolio packaging and write-up

---

# 🤝 Contributing

This is a personal learning portfolio, but suggestions and ideas are welcome.  
Pull requests focusing on documentation, testing, or new example workflows are appreciated.

---

# 📄 License

MIT License.
