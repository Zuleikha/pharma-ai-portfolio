# 🧬 Pharma-AI Portfolio  
**TL;DR:** A collection of small, focused projects in computational chemistry + structural bioinformatics. Each folder = one workflow. Minimal fluff, maximum clarity.

# 📁 Project Structure

```
pharma-ai-portfolio/
│
├─ alphafold_target_pipeline/      # Protein folding + structure analysis (AlphaFold)
│   ├─ data/
│   ├─ notebooks/
│   ├─ src/
│   └─ output/
│
├─ src/                            # Core Python scripts
│   ├─ fragment_based_drug_design.py
│   ├─ molecular_docking_prep.py
│   ├─ molecular_property_analyzer.py
│
├─ notebooks/                      # Interactive workflows
│   ├─ fragment_based_drug_design.ipynb
│   ├─ molecular_property_analysis.ipynb
│   ├─ drug_analysis_demo.ipynb
│   ├─ Molecular Docking Preparation.ipynb
│
├─ output/                         # Generated results
│   └─ molecular_analysis_results.csv
│
└─ README.md
```

# 📂 Folder Descriptions

### **src/**
Core logic scripts:
- `fragment_based_drug_design.py` — fragmenting, recombining, scoring  
- `molecular_docking_prep.py` — ligand cleaning, protonation, geometry prep  
- `molecular_property_analyzer.py` — RDKit property calculations  

### **notebooks/**
Quick interactive demos:
- Data exploration  
- Visualizations  
- End-to-end mini workflows  

### **alphafold_target_pipeline/**
Protein structure project:
- Sequence → AlphaFold prediction → pocket analysis → ligand exploration

### **output/**
Exported results:
- Property tables  
- Fragment outputs  
- Docking prep summaries  

# 🚀 Installation

```bash
git clone https://github.com/yourusername/pharma-ai-portfolio.git
cd pharma-ai-portfolio
```

Create & activate a virtual environment:

```bash
python -m venv venv
source venv/Scripts/activate  # Windows (Git Bash)
```

Install dependencies:

```bash
pip install -r requirements.txt
```

# 🧪 Requirements

Core libraries:
- RDKit  
- NumPy  
- Pandas  
- Matplotlib  
- Scikit-learn (optional)

# 📊 Example Workflows

**Fragment-Based Drug Design**  
Fragment → recombine → score → filter.

**Molecular Property Analysis**  
Compute LogP, MW, PSA, HBD/HBA → export.

**Docking Preparation**  
Clean → add H → optimize → SMILES → SDF.

**AlphaFold Pipeline**  
Protein sequence → predicted 3D structure → pocket → ligands.

# 📝 Notes  
Portfolio is actively growing — more workflows coming soon.  
Each project is intentionally small, clear, and self-contained.

# 🤝 Contributing  
PRs and suggestions welcome.

# 📄 License  
MIT License.
