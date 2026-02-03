"""
ML-Based Compound Prioritization
=================================

This module demonstrates machine learning for drug discovery, specifically
predicting molecular toxicity using Random Forest models with molecular fingerprints.

Features:
- Molecular fingerprint generation (Morgan, MACCS)
- Binary classification (toxic vs non-toxic)
- Model training and evaluation
- Compound ranking by predicted safety

Author: Pharma AI Portfolio
Date: 2025
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, MACCSkeys
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from sklearn.metrics import classification_report, confusion_matrix
import pickle
import os


class MolecularFingerprinter:
    """Generate molecular fingerprints for ML models."""
    
    def __init__(self, fingerprint_type='morgan', radius=2, n_bits=2048):
        """
        Initialize fingerprinter.
        
        Args:
            fingerprint_type: 'morgan' or 'maccs'
            radius: Radius for Morgan fingerprints (default 2)
            n_bits: Number of bits for Morgan fingerprints (default 2048)
        """
        self.fingerprint_type = fingerprint_type
        self.radius = radius
        self.n_bits = n_bits
    
    def smiles_to_fingerprint(self, smiles):
        """Convert SMILES to fingerprint vector."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        if self.fingerprint_type == 'morgan':
            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol, self.radius, nBits=self.n_bits
            )
        elif self.fingerprint_type == 'maccs':
            fp = MACCSkeys.GenMACCSKeys(mol)
        else:
            raise ValueError(f"Unknown fingerprint type: {self.fingerprint_type}")
        
        return np.array(fp)
    
    def batch_fingerprints(self, smiles_list):
        """Generate fingerprints for a list of SMILES."""
        fingerprints = []
        valid_smiles = []
        
        for smiles in smiles_list:
            fp = self.smiles_to_fingerprint(smiles)
            if fp is not None:
                fingerprints.append(fp)
                valid_smiles.append(smiles)
        
        return np.array(fingerprints), valid_smiles


class ToxicityPredictor:
    """ML model for predicting molecular toxicity."""
    
    def __init__(self, fingerprint_type='morgan', radius=2, n_bits=2048):
        """
        Initialize the toxicity predictor.
        
        Args:
            fingerprint_type: Type of fingerprint to use
            radius: Morgan fingerprint radius
            n_bits: Number of bits in fingerprint
        """
        self.fingerprinter = MolecularFingerprinter(
            fingerprint_type=fingerprint_type,
            radius=radius,
            n_bits=n_bits
        )
        self.model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            n_jobs=-1
        )
        self.is_trained = False
    
    def train(self, smiles_list, labels):
        """
        Train the toxicity prediction model.
        
        Args:
            smiles_list: List of SMILES strings
            labels: Binary labels (1 = toxic, 0 = non-toxic)
        
        Returns:
            Dictionary with training metrics
        """
        # Generate fingerprints
        X, valid_smiles = self.fingerprinter.batch_fingerprints(smiles_list)
        y = np.array([labels[i] for i, s in enumerate(smiles_list) 
                      if s in valid_smiles])
        
        # Train-test split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        
        # Train model
        print("Training Random Forest model...")
        self.model.fit(X_train, y_train)
        self.is_trained = True
        
        # Evaluate
        y_pred = self.model.predict(X_test)
        y_proba = self.model.predict_proba(X_test)[:, 1]
        
        metrics = {
            'accuracy': accuracy_score(y_test, y_pred),
            'precision': precision_score(y_test, y_pred),
            'recall': recall_score(y_test, y_pred),
            'roc_auc': roc_auc_score(y_test, y_proba),
            'confusion_matrix': confusion_matrix(y_test, y_pred),
            'n_train': len(X_train),
            'n_test': len(X_test)
        }
        
        print("\nModel Training Complete!")
        print(f"Training set: {metrics['n_train']} compounds")
        print(f"Test set: {metrics['n_test']} compounds")
        print(f"\nPerformance Metrics:")
        print(f"  Accuracy:  {metrics['accuracy']:.3f}")
        print(f"  Precision: {metrics['precision']:.3f}")
        print(f"  Recall:    {metrics['recall']:.3f}")
        print(f"  ROC-AUC:   {metrics['roc_auc']:.3f}")
        
        return metrics
    
    def predict(self, smiles_list):
        """
        Predict toxicity for new compounds.
        
        Args:
            smiles_list: List of SMILES strings
        
        Returns:
            DataFrame with predictions and probabilities
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")
        
        # Generate fingerprints
        X, valid_smiles = self.fingerprinter.batch_fingerprints(smiles_list)
        
        # Predict
        predictions = self.model.predict(X)
        probabilities = self.model.predict_proba(X)[:, 1]
        
        # Create results DataFrame
        results = pd.DataFrame({
            'SMILES': valid_smiles,
            'Predicted_Toxic': predictions.astype(bool),
            'Toxicity_Probability': probabilities,
            'Safety_Score': 1 - probabilities
        })
        
        # Sort by safety (least toxic first)
        results = results.sort_values('Safety_Score', ascending=False)
        
        return results
    
    def rank_compounds(self, smiles_list):
        """
        Rank compounds by predicted safety.
        
        Args:
            smiles_list: List of SMILES strings
        
        Returns:
            DataFrame ranked by safety score
        """
        results = self.predict(smiles_list)
        results['Rank'] = range(1, len(results) + 1)
        
        print("\nCompound Ranking (by Safety):")
        print("=" * 70)
        for _, row in results.iterrows():
            status = "SAFE" if not row['Predicted_Toxic'] else "TOXIC"
            print(f"Rank {row['Rank']:2.0f} | {status:5s} | "
                  f"Safety: {row['Safety_Score']:.3f} | "
                  f"SMILES: {row['SMILES']}")
        
        return results
    
    def save_model(self, filepath='toxicity_model.pkl'):
        """Save the trained model to disk."""
        if not self.is_trained:
            raise ValueError("No trained model to save")
        
        model_data = {
            'model': self.model,
            'fingerprinter': self.fingerprinter
        }
        
        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)
        
        print(f"\nModel saved to: {filepath}")
    
    @classmethod
    def load_model(cls, filepath='toxicity_model.pkl'):
        """Load a trained model from disk."""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)
        
        predictor = cls()
        predictor.model = model_data['model']
        predictor.fingerprinter = model_data['fingerprinter']
        predictor.is_trained = True
        
        print(f"Model loaded from: {filepath}")
        return predictor


def generate_synthetic_training_data():
    """
    Generate synthetic training data for demonstration.
    
    In real applications, use actual toxicity databases like Tox21 or ToxCast.
    """
    # Example toxic compounds (known toxins/problematic structures)
    toxic_compounds = [
        'C1=CC=C(C=C1)N(=O)=O',  # Nitrobenzene
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O.O',  # Modified problematic structure
        'C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3C2=O',  # Anthraquinone
        'C1CCOC1',  # Tetrahydrofuran (solvent toxicity)
        'C1=CC=C(C=C1)Cl',  # Chlorobenzene
        'CC(=O)NC1=CC=C(C=C1)O',  # Acetaminophen derivative
        'C1=CC=C(C=C1)CCCl',  # Alkyl halide
        'O=C1C=CC(=O)C=C1',  # Benzoquinone
        'C1=CC=C(C=C1)N=NC2=CC=CC=C2',  # Azobenzene
        'CC(C)(C)OOC(C)(C)C',  # Peroxide
    ]
    
    # Example safe compounds (approved drugs/natural products)
    safe_compounds = [
        'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin
        'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
        'CC(C)NCC(COc1ccccc1)O',  # Propranolol
        'CN(C)CCOC(c1ccccc1)c2ccccc2',  # Diphenhydramine
        'CC(C)(C)NCC(c1ccc(c(c1)CO)O)O',  # Salbutamol
        'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O',  # Penicillin G
        'COc1ccc2c(c1)c(=O)c(cn2C)C(=O)O',  # Nalidixic acid
        'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen variant
        'Cc1ccc(cc1)S(=O)(=O)N',  # Toluenesulfonamide
    ]
    
    # Create labels (1 = toxic, 0 = safe)
    all_smiles = toxic_compounds + safe_compounds
    labels = [1] * len(toxic_compounds) + [0] * len(safe_compounds)
    
    return all_smiles, labels


def demo_toxicity_prediction():
    """Demonstrate the complete ML workflow for toxicity prediction."""
    
    print("=" * 70)
    print("ML-BASED COMPOUND PRIORITIZATION - TOXICITY PREDICTION")
    print("=" * 70)
    
    # Generate training data
    print("\n[1] Generating training data...")
    smiles_list, labels = generate_synthetic_training_data()
    print(f"    Total compounds: {len(smiles_list)}")
    print(f"    Toxic: {sum(labels)}, Safe: {len(labels) - sum(labels)}")
    
    # Initialize and train model
    print("\n[2] Training ML model...")
    predictor = ToxicityPredictor(fingerprint_type='morgan', radius=2, n_bits=2048)
    metrics = predictor.train(smiles_list, labels)
    
    # Test on new compounds
    print("\n[3] Predicting toxicity for new compounds...")
    test_compounds = [
        'CC(C)Cc1ccc(cc1)C(C)C(=O)O',  # Ibuprofen (should be safe)
        'C1=CC=C(C=C1)N(=O)=O',  # Nitrobenzene (should be toxic)
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine (should be safe)
        'O=C1C=CC(=O)C=C1',  # Quinone (should be toxic)
        'CC(=O)Oc1ccccc1C(=O)O',  # Aspirin (should be safe)
    ]
    
    results = predictor.rank_compounds(test_compounds)
    
    # Save model
    print("\n[4] Saving model...")
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
    os.makedirs(output_dir, exist_ok=True)
    model_path = os.path.join(output_dir, 'toxicity_model.pkl')
    predictor.save_model(model_path)
    
    print("\n" + "=" * 70)
    print("DEMO COMPLETE")
    print("=" * 70)
    
    return predictor, results


def main():
    """Run the demonstration."""
    predictor, results = demo_toxicity_prediction()
    
    # Show summary statistics
    print("\nSummary Statistics:")
    print(f"  Total compounds ranked: {len(results)}")
    print(f"  Predicted safe: {sum(~results['Predicted_Toxic'])}")
    print(f"  Predicted toxic: {sum(results['Predicted_Toxic'])}")
    print(f"  Average safety score: {results['Safety_Score'].mean():.3f}")


if __name__ == "__main__":
    main()