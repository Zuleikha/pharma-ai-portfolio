# 📘 Pharma AI Portfolio

A collection of small, focused projects demonstrating practical applications of cheminformatics, AI, and computational drug design.  
This portfolio includes Python scripts, Jupyter notebooks, and analysis outputs related to:

- Fragment-Based Drug Design (FBDD)
- Molecular Property Analysis
- Molecular Docking Preparation
- RDKit cheminformatics workflows

## 📂 Project Structure

```
pharma-ai-portfolio/
│
├─ src/
│   ├─ fragment_based_drug_design.py
│   ├─ molecular_docking_prep.py
│   ├─ molecular_property_analyzer.py
│
├─ notebooks/
│   ├─ fragment_based_drug_design.ipynb
│   ├─ molecular_property_analysis.ipynb
│   ├─ drug_analysis_demo.ipynb
│   ├─ Molecular Docking Preparation.ipynb
│
├─ output/
│   └─ molecular_analysis_results.csv
│
└─ README.md
```

### Folder Descriptions

#### `src/`
Contains Python scripts for core logic:
- **fragment_based_drug_design.py** — performs fragment decomposition, recombination, and scoring.
- **molecular_docking_prep.py** — prepares ligands for docking (cleaning, protonation, etc.).
- **molecular_property_analyzer.py** — calculates RDKit-based physicochemical properties.

#### `notebooks/`
Interactive Jupyter notebooks used for:
- Demonstrations
- Data exploration
- Visualizations
- Step-by-step workflows

#### `output/`
Contains exported results from analyses:
- Molecular property tables
- Fragment analysis outputs
- Docking prep summaries

## 🚀 Installation

### 1. Clone the repository:
```bash
git clone https://github.com/yourusername/pharma-ai-portfolio.git
cd pharma-ai-portfolio
```

### 2. Create & activate a virtual environment:
```bash
python -m venv venv
source venv/Scripts/activate   # Git Bash on Windows
```

### 3. Install dependencies:
```bash
pip install -r requirements.txt
```

## 🧪 Requirements
Main libraries used:

- RDKit
- NumPy
- Pandas
- Matplotlib
- Scikit-learn (optional)

Add them to your `requirements.txt` if you plan to publish the repo.

## 📊 Example Workflows

### Fragment-Based Drug Design
- Decompose molecules into fragments
- Recombine fragments
- Score candidates using property filters

### Molecular Property Analysis
- Calculate LogP, MW, PSA, HBD/HBA, etc.
- Export tabular results

### Docking Preparation
- Clean molecular structures
- Add hydrogens and optimize geometry
- Prepare SMILES → SDF conversion pipeline

## 📝 Notes
This portfolio is a work-in-progress collection of computational chemistry workflows.  
More modules and notebook demos will be added over time.

## 🤝 Contributing
Pull requests and improvements are welcome.

## 📄 License
MIT License.
