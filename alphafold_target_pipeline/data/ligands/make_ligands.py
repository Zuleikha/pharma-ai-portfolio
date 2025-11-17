from rdkit import Chem
from rdkit.Chem import AllChem

# ---- LIGAND LIST ----
ligands = {
    "trimethoprim": "COc1cc(C)c(Nc2nc(N)nc(N)c2)cc1OC",
    "pyrimethamine": "CC1=NC(=NC(=N1)Cl)N(c2ccc(cc2)Cl)C",
    "fragment_core": "NC1=NC=NC(N)=N1"
}

# ---- GENERATE 3D STRUCTURES ----
for name, smiles in ligands.items():
    print(f"Processing {name}...")

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    Chem.MolToMolFile(mol, f"{name}.mol")

print("Done! 3D ligand files created.")
