from rdkit import Chem
from rdkit.Chem import rdMolAlign

ref = Chem.MolFromMolFile(
    r"alphafold_target_pipeline\data\ligands\trimethoprim.mol",
    removeHs=True
)

dock_raw = Chem.MolFromMolFile(
    r"alphafold_target_pipeline\output\docking\trimethoprim_docked_redock.sdf",
    removeHs=False,
    sanitize=False
)

if ref is None or dock_raw is None:
    raise ValueError("Failed to load one of the molecules")

dock = Chem.RemoveHs(dock_raw)

print("ref atoms:", ref.GetNumAtoms())
print("dock atoms:", dock.GetNumAtoms())

if ref.GetNumAtoms() != dock.GetNumAtoms():
    raise ValueError("Atom counts still differ – cannot compute RMSD safely")

atom_map = list(zip(range(ref.GetNumAtoms()), range(dock.GetNumAtoms())))

rmsd = rdMolAlign.AlignMol(dock, ref, atomMap=atom_map)
print(f"Trimethoprim redocking RMSD: {rmsd:.3f} Å")
