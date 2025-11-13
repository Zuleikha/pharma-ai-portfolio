"""
Fragment-Based Drug Design
==========================

This module provides tools for fragment-based drug design (FBDD) using RDKit.
It includes methods for fragment generation, analysis, and optimization.

Fragment-based drug design is a powerful approach in drug discovery that involves
identifying small chemical fragments that bind to different parts of a biological
target, then linking or growing these fragments into larger, lead-like molecules.

Author: Pharma AI Portfolio
Date: 2025
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski
from rdkit.Chem import BRICS, Recap, Scaffolds
from rdkit.Chem import Fragments
from rdkit.Chem import rdMolDescriptors
from typing import List, Set, Dict, Tuple, Optional
import pandas as pd
from collections import Counter


class FragmentBasedDesign:
    """
    Tools for fragment-based drug design and analysis.

    This class provides methods for:
    - Fragment generation using RECAP and BRICS
    - Scaffold extraction
    - Fragment library creation
    - Fragment property analysis
    """

    def __init__(self):
        """Initialize FragmentBasedDesign."""
        self.fragment_mw_limit = 300  # Typical fragment MW limit
        self.fragment_hba_limit = 3
        self.fragment_hbd_limit = 3
        self.fragment_logp_limit = 3

    def generate_recap_fragments(self, smiles: str) -> List[str]:
        """
        Generate fragments using RECAP (Retrosynthetic Combinatorial Analysis Procedure).

        RECAP breaks molecules at chemical bonds that are commonly formed in
        synthetic chemistry.

        Args:
            smiles: SMILES string of molecule

        Returns:
            List of fragment SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error: Invalid SMILES string: {smiles}")
            return []

        recap_tree = Recap.RecapDecompose(mol)
        fragments = []

        for node in recap_tree.GetLeaves().values():
            fragment_smiles = node.mol
            if fragment_smiles:
                fragments.append(fragment_smiles)

        return list(set(fragments))  # Remove duplicates

    def generate_brics_fragments(self, smiles: str) -> List[str]:
        """
        Generate fragments using BRICS (Breaking of Retrosynthetically Interesting
        Chemical Substructures).

        BRICS is similar to RECAP but uses a different set of bond-breaking rules.

        Args:
            smiles: SMILES string of molecule

        Returns:
            List of fragment SMILES
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error: Invalid SMILES string: {smiles}")
            return []

        fragments = list(BRICS.BRICSDecompose(mol))
        # Remove dummy atoms from fragments
        cleaned_fragments = []
        for frag in fragments:
            # Remove BRICS dummy atoms (marked with [*])
            cleaned = frag.replace('[*]', '[H]')
            mol_frag = Chem.MolFromSmiles(cleaned)
            if mol_frag:
                cleaned_fragments.append(Chem.MolToSmiles(mol_frag))

        return list(set(cleaned_fragments))

    def extract_murcko_scaffold(self, smiles: str) -> Optional[str]:
        """
        Extract Murcko scaffold (framework) from molecule.

        The Murcko scaffold represents the core ring systems and linkers,
        removing all side chains.

        Args:
            smiles: SMILES string of molecule

        Returns:
            SMILES of scaffold or None if extraction fails
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        scaffold = Scaffolds.MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaffold) if scaffold else None

    def extract_generic_scaffold(self, smiles: str) -> Optional[str]:
        """
        Extract generic Murcko scaffold (all atoms as carbon, all bonds as single).

        Args:
            smiles: SMILES string of molecule

        Returns:
            SMILES of generic scaffold or None if extraction fails
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        scaffold = Scaffolds.MurckoScaffold.MakeScaffoldGeneric(mol)
        return Chem.MolToSmiles(scaffold) if scaffold else None

    def calculate_fragment_properties(self, smiles: str) -> Dict[str, float]:
        """
        Calculate properties relevant for fragment screening.

        Args:
            smiles: SMILES string of fragment

        Returns:
            Dictionary of fragment properties
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}

        properties = {
            'MolecularWeight': Descriptors.MolWt(mol),
            'LogP': Crippen.MolLogP(mol),
            'HBD': Lipinski.NumHDonors(mol),
            'HBA': Lipinski.NumHAcceptors(mol),
            'TPSA': Descriptors.TPSA(mol),
            'HeavyAtoms': mol.GetNumHeavyAtoms(),
            'RotatableBonds': Lipinski.NumRotatableBonds(mol),
            'AromaticRings': Lipinski.NumAromaticRings(mol),
            'SaturatedRings': Lipinski.NumSaturatedRings(mol),
            'Complexity': rdMolDescriptors.BertzCT(mol)
        }
        return properties

    def is_valid_fragment(self, smiles: str) -> bool:
        """
        Check if fragment meets typical fragment screening criteria.

        Rule of Three for fragments:
        - MW <= 300
        - LogP <= 3
        - HBD <= 3
        - HBA <= 3
        - Rotatable bonds <= 3
        - TPSA <= 60

        Args:
            smiles: SMILES string of fragment

        Returns:
            True if fragment passes filters
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rotatable = Lipinski.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)

        return (
            mw <= 300 and
            logp <= 3 and
            hbd <= 3 and
            hba <= 3 and
            rotatable <= 3 and
            tpsa <= 60
        )

    def filter_fragments(self, fragments: List[str]) -> List[str]:
        """
        Filter fragments based on Rule of Three criteria.

        Args:
            fragments: List of fragment SMILES

        Returns:
            List of valid fragment SMILES
        """
        valid_fragments = []
        for frag in fragments:
            if self.is_valid_fragment(frag):
                valid_fragments.append(frag)
        return valid_fragments

    def generate_fragment_library(self, smiles_list: List[str],
                                  method: str = 'brics') -> List[str]:
        """
        Generate a fragment library from a list of molecules.

        Args:
            smiles_list: List of molecule SMILES
            method: Fragmentation method ('brics' or 'recap')

        Returns:
            List of unique fragments
        """
        all_fragments = set()

        for smiles in smiles_list:
            if method.lower() == 'brics':
                fragments = self.generate_brics_fragments(smiles)
            elif method.lower() == 'recap':
                fragments = self.generate_recap_fragments(smiles)
            else:
                raise ValueError("Method must be 'brics' or 'recap'")

            all_fragments.update(fragments)

        return list(all_fragments)

    def score_fragment(self, smiles: str) -> float:
        """
        Calculate a simple fragment quality score.

        Score based on:
        - Ligand efficiency (arbitrary binding energy / heavy atoms)
        - Lipophilic efficiency
        - Molecular complexity

        Args:
            smiles: Fragment SMILES

        Returns:
            Fragment score (higher is better)
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0

        heavy_atoms = mol.GetNumHeavyAtoms()
        if heavy_atoms == 0:
            return 0.0

        logp = Crippen.MolLogP(mol)
        aromatic_rings = Lipinski.NumAromaticRings(mol)
        complexity = rdMolDescriptors.BertzCT(mol)

        # Simple scoring function (can be customized)
        score = (aromatic_rings * 2 + 5) / heavy_atoms
        score -= abs(logp) * 0.1  # Penalize high LogP
        score -= complexity / 1000  # Penalize high complexity

        return max(0, score)

    def analyze_fragment_library(self, fragments: List[str]) -> pd.DataFrame:
        """
        Analyze a fragment library and return statistics.

        Args:
            fragments: List of fragment SMILES

        Returns:
            DataFrame with fragment properties and scores
        """
        results = []

        for frag in fragments:
            props = self.calculate_fragment_properties(frag)
            if props:
                props['SMILES'] = frag
                props['IsValid'] = self.is_valid_fragment(frag)
                props['Score'] = self.score_fragment(frag)
                results.append(props)

        df = pd.DataFrame(results)
        return df.sort_values('Score', ascending=False)

    def find_common_scaffolds(self, smiles_list: List[str],
                             min_frequency: int = 2) -> List[Tuple[str, int]]:
        """
        Find common scaffolds in a set of molecules.

        Args:
            smiles_list: List of molecule SMILES
            min_frequency: Minimum number of occurrences

        Returns:
            List of (scaffold_smiles, count) tuples
        """
        scaffolds = []

        for smiles in smiles_list:
            scaffold = self.extract_murcko_scaffold(smiles)
            if scaffold:
                scaffolds.append(scaffold)

        counter = Counter(scaffolds)
        common_scaffolds = [
            (scaffold, count)
            for scaffold, count in counter.items()
            if count >= min_frequency
        ]

        return sorted(common_scaffolds, key=lambda x: x[1], reverse=True)

    def merge_fragments(self, frag1_smiles: str, frag2_smiles: str) -> List[str]:
        """
        Generate possible merged molecules from two fragments (simple concatenation).

        Note: This is a simplified example. Real fragment merging requires
        sophisticated chemistry-aware algorithms.

        Args:
            frag1_smiles: First fragment SMILES
            frag2_smiles: Second fragment SMILES

        Returns:
            List of merged molecule SMILES
        """
        frag1 = Chem.MolFromSmiles(frag1_smiles)
        frag2 = Chem.MolFromSmiles(frag2_smiles)

        if frag1 is None or frag2 is None:
            return []

        # Simple example: use BRICS to recombine
        merged = []
        try:
            # This is a placeholder - real implementation would use
            # sophisticated fragment linking algorithms
            combined_smiles = f"{frag1_smiles}.{frag2_smiles}"
            merged.append(combined_smiles)
        except Exception as e:
            print(f"Error merging fragments: {e}")

        return merged


def main():
    """Example usage of FragmentBasedDesign."""
    designer = FragmentBasedDesign()

    # Example drug molecules
    example_drugs = {
        'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'Aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
        'Paracetamol': 'CC(=O)Nc1ccc(O)cc1',
        'Warfarin': 'CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O'
    }

    print("="*70)
    print("FRAGMENT-BASED DRUG DESIGN ANALYSIS")
    print("="*70)

    # Example 1: BRICS Fragmentation
    print("\n1. BRICS Fragmentation of Ibuprofen:")
    print("-" * 70)
    brics_frags = designer.generate_brics_fragments(example_drugs['Ibuprofen'])
    for i, frag in enumerate(brics_frags, 1):
        print(f"   Fragment {i}: {frag}")

    # Example 2: RECAP Fragmentation
    print("\n2. RECAP Fragmentation of Warfarin:")
    print("-" * 70)
    recap_frags = designer.generate_recap_fragments(example_drugs['Warfarin'])
    for i, frag in enumerate(recap_frags, 1):
        print(f"   Fragment {i}: {frag}")

    # Example 3: Scaffold Extraction
    print("\n3. Scaffold Extraction:")
    print("-" * 70)
    for name, smiles in example_drugs.items():
        scaffold = designer.extract_murcko_scaffold(smiles)
        print(f"   {name:15s}: {scaffold}")

    # Example 4: Fragment Library Generation
    print("\n4. Generating Fragment Library:")
    print("-" * 70)
    fragment_library = designer.generate_fragment_library(
        list(example_drugs.values()),
        method='brics'
    )
    print(f"   Total unique fragments: {len(fragment_library)}")

    # Example 5: Filter valid fragments
    valid_fragments = designer.filter_fragments(fragment_library)
    print(f"   Valid fragments (Rule of 3): {len(valid_fragments)}")

    # Example 6: Analyze fragment library
    print("\n5. Fragment Analysis (Top 5 by score):")
    print("-" * 70)
    if valid_fragments:
        df = designer.analyze_fragment_library(valid_fragments[:10])
        print(df[['SMILES', 'MolecularWeight', 'LogP', 'Score', 'IsValid']].head().to_string(index=False))

    # Example 7: Find common scaffolds
    print("\n6. Common Scaffolds:")
    print("-" * 70)
    common_scaffolds = designer.find_common_scaffolds(list(example_drugs.values()))
    for scaffold, count in common_scaffolds:
        print(f"   Count: {count} - {scaffold}")

    print("\n" + "="*70)


if __name__ == "__main__":
    main()
