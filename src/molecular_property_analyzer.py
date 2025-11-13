"""
Molecular Property Analyzer
============================

This module provides tools for analyzing molecular properties and drug-likeness
using RDKit. It calculates various descriptors and evaluates molecules against
pharmaceutical property rules.

Author: Pharma AI Portfolio
Date: 2025
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, AllChem
from rdkit.Chem import rdMolDescriptors
import pandas as pd
from typing import Dict, List, Optional, Union


class MolecularPropertyAnalyzer:
    """
    A comprehensive molecular property analyzer for drug discovery.

    This class provides methods to calculate molecular descriptors,
    evaluate drug-likeness, and assess ADME properties.
    """

    def __init__(self):
        """Initialize the MolecularPropertyAnalyzer."""
        self.property_names = [
            'MolecularWeight', 'LogP', 'HBD', 'HBA', 'TPSA',
            'RotatableBonds', 'AromaticRings', 'QED'
        ]

    def smiles_to_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Convert SMILES string to RDKit molecule object.

        Args:
            smiles: SMILES string representation of molecule

        Returns:
            RDKit molecule object or None if invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error: Invalid SMILES string: {smiles}")
        return mol

    def calculate_basic_properties(self, mol: Chem.Mol) -> Dict[str, float]:
        """
        Calculate basic molecular properties.

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary of property names and values
        """
        properties = {
            'MolecularWeight': Descriptors.MolWt(mol),
            'LogP': Crippen.MolLogP(mol),
            'HBD': Lipinski.NumHDonors(mol),
            'HBA': Lipinski.NumHAcceptors(mol),
            'TPSA': Descriptors.TPSA(mol),
            'RotatableBonds': Lipinski.NumRotatableBonds(mol),
            'AromaticRings': Lipinski.NumAromaticRings(mol),
            'HeavyAtomCount': Lipinski.HeavyAtomCount(mol),
            'FractionCsp3': Lipinski.FractionCSP3(mol),
            'MolarRefractivity': Crippen.MolMR(mol)
        }
        return properties

    def calculate_lipinski_properties(self, mol: Chem.Mol) -> Dict[str, Union[float, bool]]:
        """
        Calculate Lipinski's Rule of Five properties.

        The Rule of Five states that, in general, an orally active drug has:
        - MW <= 500 Da
        - LogP <= 5
        - HBD <= 5
        - HBA <= 10

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary with Lipinski properties and pass/fail status
        """
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)

        lipinski_pass = (
            mw <= 500 and
            logp <= 5 and
            hbd <= 5 and
            hba <= 10
        )

        violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10
        ])

        return {
            'MW': mw,
            'LogP': logp,
            'HBD': hbd,
            'HBA': hba,
            'LipinskiPass': lipinski_pass,
            'Violations': violations
        }

    def calculate_qed(self, mol: Chem.Mol) -> float:
        """
        Calculate Quantitative Estimate of Drug-likeness (QED).

        QED is a measure of drug-likeness based on the concept of desirability.
        Values range from 0 (unfavorable) to 1 (favorable).

        Args:
            mol: RDKit molecule object

        Returns:
            QED score (0-1)
        """
        return QED.qed(mol)

    def calculate_veber_properties(self, mol: Chem.Mol) -> Dict[str, Union[float, bool]]:
        """
        Calculate Veber's rules for oral bioavailability.

        Veber's rules state that oral bioavailability is likely if:
        - Rotatable bonds <= 10
        - TPSA <= 140 Ų

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary with Veber properties and pass/fail status
        """
        rotatable_bonds = Lipinski.NumRotatableBonds(mol)
        tpsa = Descriptors.TPSA(mol)

        veber_pass = (rotatable_bonds <= 10 and tpsa <= 140)

        return {
            'RotatableBonds': rotatable_bonds,
            'TPSA': tpsa,
            'VeberPass': veber_pass
        }

    def calculate_ghose_filter(self, mol: Chem.Mol) -> Dict[str, Union[float, bool]]:
        """
        Calculate Ghose filter properties.

        Ghose filter ranges for drug-like molecules:
        - MW: 160-480
        - LogP: -0.4 to 5.6
        - Atom count: 20-70
        - Molar refractivity: 40-130

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary with Ghose properties and pass/fail status
        """
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        atoms = mol.GetNumHeavyAtoms()
        mr = Crippen.MolMR(mol)

        ghose_pass = (
            160 <= mw <= 480 and
            -0.4 <= logp <= 5.6 and
            20 <= atoms <= 70 and
            40 <= mr <= 130
        )

        return {
            'MW': mw,
            'LogP': logp,
            'Atoms': atoms,
            'MolarRefractivity': mr,
            'GhosePass': ghose_pass
        }

    def analyze_molecule(self, smiles: str, verbose: bool = True) -> Dict:
        """
        Perform comprehensive analysis of a molecule.

        Args:
            smiles: SMILES string representation of molecule
            verbose: If True, print detailed results

        Returns:
            Dictionary containing all calculated properties
        """
        mol = self.smiles_to_mol(smiles)
        if mol is None:
            return {}

        results = {
            'SMILES': smiles,
            'BasicProperties': self.calculate_basic_properties(mol),
            'Lipinski': self.calculate_lipinski_properties(mol),
            'QED': self.calculate_qed(mol),
            'Veber': self.calculate_veber_properties(mol),
            'Ghose': self.calculate_ghose_filter(mol)
        }

        if verbose:
            self._print_results(results)

        return results

    def _print_results(self, results: Dict):
        """Print formatted analysis results."""
        print("\n" + "="*60)
        print("MOLECULAR PROPERTY ANALYSIS")
        print("="*60)
        print(f"\nSMILES: {results['SMILES']}")

        print("\n--- Basic Properties ---")
        for key, value in results['BasicProperties'].items():
            print(f"{key:20s}: {value:.2f}")

        print("\n--- Lipinski's Rule of Five ---")
        lipinski = results['Lipinski']
        print(f"Molecular Weight    : {lipinski['MW']:.2f} (≤ 500)")
        print(f"LogP               : {lipinski['LogP']:.2f} (≤ 5)")
        print(f"H-Bond Donors      : {lipinski['HBD']} (≤ 5)")
        print(f"H-Bond Acceptors   : {lipinski['HBA']} (≤ 10)")
        print(f"Violations         : {lipinski['Violations']}")
        print(f"Pass               : {'✓ Yes' if lipinski['LipinskiPass'] else '✗ No'}")

        print(f"\n--- Drug-likeness ---")
        print(f"QED Score          : {results['QED']:.3f}")

        print("\n--- Veber's Rules ---")
        veber = results['Veber']
        print(f"Rotatable Bonds    : {veber['RotatableBonds']} (≤ 10)")
        print(f"TPSA              : {veber['TPSA']:.2f} (≤ 140)")
        print(f"Pass               : {'✓ Yes' if veber['VeberPass'] else '✗ No'}")

        print("\n--- Ghose Filter ---")
        ghose = results['Ghose']
        print(f"Pass               : {'✓ Yes' if ghose['GhosePass'] else '✗ No'}")
        print("="*60 + "\n")

    def analyze_batch(self, smiles_list: List[str]) -> pd.DataFrame:
        """
        Analyze multiple molecules and return results as DataFrame.

        Args:
            smiles_list: List of SMILES strings

        Returns:
            Pandas DataFrame with all properties
        """
        results = []
        for smiles in smiles_list:
            analysis = self.analyze_molecule(smiles, verbose=False)
            if analysis:
                flat_result = {'SMILES': smiles}
                flat_result.update(analysis['BasicProperties'])
                flat_result.update(analysis['Lipinski'])
                flat_result['QED'] = analysis['QED']
                results.append(flat_result)

        return pd.DataFrame(results)


def main():
    """Example usage of MolecularPropertyAnalyzer."""
    analyzer = MolecularPropertyAnalyzer()

    # Example molecules
    examples = {
        'Aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
        'Ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
        'Caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'Penicillin G': 'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O'
    }

    print("Analyzing example drug molecules...\n")

    for name, smiles in examples.items():
        print(f"\n{'='*60}")
        print(f"Analyzing: {name}")
        print(f"{'='*60}")
        analyzer.analyze_molecule(smiles)

    # Batch analysis example
    print("\n\nBatch Analysis Example:")
    print("-" * 60)
    smiles_list = list(examples.values())
    df = analyzer.analyze_batch(smiles_list)
    print(df[['SMILES', 'MolecularWeight', 'LogP', 'QED', 'LipinskiPass']].to_string())


if __name__ == "__main__":
    main()
