# Project 3 – ML-assisted compound prioritisation
# Purpose: rank compounds using simple RDKit features

from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

def featurise_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
    }

def rank_compounds(smiles_list):
    rows = []
    for s in smiles_list:
        feats = featurise_molecule(s)
        if feats:
            feats["SMILES"] = s
            rows.append(feats)

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Simple heuristic score (placeholder for ML)
    df["Score"] = (
        -0.01 * df["MW"]
        -0.2 * df["LogP"].abs()
        -0.05 * df["TPSA"]
        +0.5 * (df["HBD"] + df["HBA"])
    )

    return df.sort_values("Score", ascending=False)
