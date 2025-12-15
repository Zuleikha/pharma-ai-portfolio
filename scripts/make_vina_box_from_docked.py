from rdkit import Chem

mol = Chem.MolFromMolFile(
    r"alphafold_target_pipeline\output\docking\trimethoprim_docked.sdf",
    removeHs=False
)
if mol is None:
    raise ValueError("Could not load docked SDF")

conf = mol.GetConformer()

xs, ys, zs = [], [], []
for i in range(mol.GetNumAtoms()):
    pos = conf.GetAtomPosition(i)
    xs.append(pos.x)
    ys.append(pos.y)
    zs.append(pos.z)

def mid(a, b):
    return (a + b) / 2.0

center_x = mid(min(xs), max(xs))
center_y = mid(min(ys), max(ys))
center_z = mid(min(zs), max(zs))

# Box size = ligand span + 10 Å margin
size_x = (max(xs) - min(xs)) + 10.0
size_y = (max(ys) - min(ys)) + 10.0
size_z = (max(zs) - min(zs)) + 10.0

print("center_x =", center_x)
print("center_y =", center_y)
print("center_z =", center_z)
print("size_x =", size_x)
print("size_y =", size_y)
print("size_z =", size_z)
