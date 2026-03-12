"""
完整多靶点 pIC50 预测系统 - 最终修复版
包含所有核心功能：Transformer、GNN、集成学习、不确实性估计、特征提取等
基于 https://github.com/zapabob/multi-target-pIC50-predictor

所有已知问题已修复：
1. ✅ RDKit指纹API修复
2. ✅ 目录创建修复
3. ✅ Ensemble模型预测修复
"""

import os
import sys
import warnings
import logging
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from datetime import datetime
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
import joblib

# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Lipinski, Crippen, MolSurf
from rdkit.Chem import rdMolDescriptors, GraphDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

# 深度学习 (可选)
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch.utils.data import Dataset, DataLoader
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    warnings.warn("PyTorch not available. Deep learning models will be disabled.")

warnings.filterwarnings('ignore')

# ==================== 配置管理 ====================

@dataclass
class ModelConfig:
    """模型配置"""
    targets: Dict[str, str] = None
    use_rdkit_descriptors: bool = True
    use_fingerprints: bool = True
    use_smarts_features: bool = True
    use_graph_features: bool = True
    fp_radius: int = 2
    fp_nbits: int = 2048
    use_ecfp: bool = True
    use_maccs: bool = True
    use_transformer: bool = True
    use_gnn: bool = True
    use_ensemble: bool = True
    test_size: float = 0.2
    random_state: int = 42
    n_splits: int = 5
    d_model: int = 256
    nhead: int = 8
    num_layers: int = 4
    dim_feedforward: int = 1024
    dropout: float = 0.1
    hidden_dim: int = 128
    num_gnn_layers: int = 3
    rf_n_estimators: int = 200
    rf_max_depth: int = 30
    gb_n_estimators: int = 100
    gb_learning_rate: float = 0.1
    uncertainty_estimation: bool = True
    mc_dropout_samples: int = 50
    ensemble_members: int = 5
    use_optuna: bool = False
    optuna_n_trials: int = 100
    
    def __post_init__(self):
        if self.targets is None:
            self.targets = {
                'DAT': 'CHEMBL238',
                '5HT2A': 'CHEMBL224',
                'CB1': 'CHEMBL218',
                'CB2': 'CHEMBL1861',
                'mu_opioid': 'CHEMBL233',
                'delta_opioid': 'CHEMBL236',
                'kappa_opioid': 'CHEMBL237'
            }

class Config:
    """全局配置"""
    BASE_DIR = Path.cwd()
    DATA_DIR = BASE_DIR / 'data'
    MODEL_DIR = BASE_DIR / 'models'
    OUTPUT_DIR = BASE_DIR / 'output'
    LOG_DIR = BASE_DIR / 'logs'
    CACHE_DIR = BASE_DIR / 'cache'
    
    @classmethod
    def setup_directories(cls):
        for dir_path in [cls.DATA_DIR, cls.MODEL_DIR, cls.OUTPUT_DIR, cls.LOG_DIR, cls.CACHE_DIR]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    @classmethod
    def setup_logging(cls, level=logging.INFO):
        cls.LOG_DIR.mkdir(parents=True, exist_ok=True)
        log_file = cls.LOG_DIR / f"pic50_predictor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
        )
        return logging.getLogger(__name__)

# ==================== 分子特征提取器 ====================

class MolecularFeatureExtractor:
    """完整的分子特征提取器"""
    
    def __init__(self, config: ModelConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.descriptor_names = [desc[0] for desc in Descriptors._descList]
        self.descriptor_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(self.descriptor_names)
        self.smarts_patterns = self._init_smarts_patterns()
    
    def _init_smarts_patterns(self) -> Dict[str, str]:
        return {
            'phenethylamine': 'c1ccccc1CCN', 'tryptamine': 'c1ccc2c(c1)c(cn2)CCN',
            'benzodioxole': 'c1cc2c(cc1)OCO2', 'piperidine': 'C1CCNCC1',
            'morpholine': 'C1COCCN1', 'tropane': 'C1CC2CCC(C1)N2',
            'indole': 'c1ccc2c(c1)cc[nH]2', 'quinoline': 'c1ccc2c(c1)cccn2',
            'benzene': 'c1ccccc1', 'pyridine': 'c1ccncc1',
            'thiophene': 'c1ccsc1', 'furan': 'c1ccoc1',
            'imidazole': 'c1c[nH]cn1', 'pyrrole': 'c1cc[nH]c1',
            'oxazole': 'c1cnco1', 'thiazole': 'c1cncs1',
            'basic_amine': '[NX3;H2,H1;!$(NC=O)]', 'secondary_amine': '[NX3;H1;!$(NC=O)]',
            'tertiary_amine': '[NX3;H0;!$(NC=O)]', 'aromatic_amine': '[NX3;H2,H1;!$(NC=O)]c',
            'methoxy': 'COc', 'ethoxy': 'CCOc', 'hydroxyl': '[OH]',
            'carbonyl': '[CX3]=[OX1]', 'carboxylic_acid': 'C(=O)[OH]',
            'ester': 'C(=O)O[CX4]', 'amide': 'C(=O)N',
            'sulfonamide': 'S(=O)(=O)N', 'halogen': '[F,Cl,Br,I]',
            'nitro': '[N+](=O)[O-]',
        }
    
    def smiles_to_mol(self, smiles: str) -> Optional[Chem.Mol]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.logger.warning(f"Invalid SMILES: {smiles}")
            return mol
        except Exception as e:
            self.logger.error(f"Error converting SMILES {smiles}: {e}")
            return None
    
    def calculate_rdkit_descriptors(self, mol: Chem.Mol) -> np.ndarray:
        try:
            descriptors = self.descriptor_calculator.CalcDescriptors(mol)
            descriptors = np.array(descriptors, dtype=np.float32)
            return np.nan_to_num(descriptors, nan=0.0, posinf=0.0, neginf=0.0)
        except Exception as e:
            self.logger.error(f"Error calculating descriptors: {e}")
            return np.zeros(len(self.descriptor_names), dtype=np.float32)
    
    def calculate_additional_descriptors(self, mol: Chem.Mol) -> Dict[str, float]:
        try:
            return {
                'MolWt': Descriptors.MolWt(mol), 'LogP': Crippen.MolLogP(mol),
                'NumHDonors': Lipinski.NumHDonors(mol), 'NumHAcceptors': Lipinski.NumHAcceptors(mol),
                'TPSA': MolSurf.TPSA(mol), 'NumRotatableBonds': Lipinski.NumRotatableBonds(mol),
                'NumAromaticRings': Lipinski.NumAromaticRings(mol),
                'NumSaturatedRings': Lipinski.NumSaturatedRings(mol),
                'NumAliphaticRings': Lipinski.NumAliphaticRings(mol),
                'RingCount': Lipinski.RingCount(mol), 'NumHeteroatoms': Lipinski.NumHeteroatoms(mol),
                'MolMR': Crippen.MolMR(mol), 'FractionCsp3': Lipinski.FractionCSP3(mol),
                'BertzCT': GraphDescriptors.BertzCT(mol), 'Chi0': GraphDescriptors.Chi0(mol),
                'Chi1': GraphDescriptors.Chi1(mol), 'Kappa1': GraphDescriptors.Kappa1(mol),
                'Kappa2': GraphDescriptors.Kappa2(mol), 'Kappa3': GraphDescriptors.Kappa3(mol),
            }
        except Exception as e:
            self.logger.error(f"Error calculating additional descriptors: {e}")
            return {k: 0.0 for k in ['MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors', 'TPSA',
                                      'NumRotatableBonds', 'NumAromaticRings', 'NumSaturatedRings',
                                      'NumAliphaticRings', 'RingCount', 'NumHeteroatoms', 'MolMR',
                                      'FractionCsp3', 'BertzCT', 'Chi0', 'Chi1', 'Kappa1', 'Kappa2', 'Kappa3']}
    
    def calculate_fingerprints(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        fingerprints = {}
        try:
            if self.config.use_ecfp:
                ecfp = AllChem.GetMorganFingerprintAsBitVect(mol, self.config.fp_radius, nBits=self.config.fp_nbits)
                fingerprints['ecfp'] = np.array(ecfp, dtype=np.float32)
            if self.config.use_maccs:
                maccs = AllChem.GetMACCSKeysFingerprint(mol)
                fingerprints['maccs'] = np.array(maccs, dtype=np.float32)
            try:
                ap = rdMolDescriptors.GetAtomPairFingerprint(mol)
                ap_bits = np.zeros(2048, dtype=np.float32)
                for idx, val in ap.GetNonzeroElements().items():
                    if idx < 2048:
                        ap_bits[idx] = min(val, 1.0)
                fingerprints['atom_pair'] = ap_bits
            except:
                pass
            try:
                tt = rdMolDescriptors.GetTopologicalTorsionFingerprint(mol)
                tt_bits = np.zeros(2048, dtype=np.float32)
                for idx, val in tt.GetNonzeroElements().items():
                    if idx < 2048:
                        tt_bits[idx] = min(val, 1.0)
                fingerprints['torsion'] = tt_bits
            except:
                pass
        except Exception as e:
            self.logger.error(f"Error calculating fingerprints: {e}")
        return fingerprints
    
    def calculate_smarts_features(self, mol: Chem.Mol) -> np.ndarray:
        features = []
        for name, pattern in self.smarts_patterns.items():
            try:
                smarts_mol = Chem.MolFromSmarts(pattern)
                count = len(mol.GetSubstructMatches(smarts_mol)) if smarts_mol else 0
                features.append(count)
            except:
                features.append(0)
        return np.array(features, dtype=np.float32)
    
    def calculate_graph_features(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        try:
            atom_features = [[atom.GetAtomicNum(), atom.GetDegree(), atom.GetFormalCharge(),
                             atom.GetHybridization(), atom.GetNumRadicalElectrons(),
                             atom.GetIsAromatic(), atom.IsInRing(), atom.GetTotalValence()]
                            for atom in mol.GetAtoms()]
            bond_features, adjacency = [], np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()))
            for bond in mol.GetBonds():
                i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                bond_type = bond.GetBondTypeAsDouble()
                adjacency[i, j] = adjacency[j, i] = bond_type
                bond_features.append([bond_type, bond.GetIsConjugated(), bond.IsInRing()])
            return {
                'atom_features': np.array(atom_features, dtype=np.float32),
                'bond_features': np.array(bond_features, dtype=np.float32) if bond_features else np.array([[0, 0, 0]], dtype=np.float32),
                'adjacency': adjacency.astype(np.float32)
            }
        except Exception as e:
            self.logger.error(f"Error calculating graph features: {e}")
            return {'atom_features': np.array([[0]*8], dtype=np.float32),
                   'bond_features': np.array([[0, 0, 0]], dtype=np.float32),
                   'adjacency': np.array([[0]], dtype=np.float32)}
    
    def extract_all_features(self, smiles: str) -> Optional[Dict[str, Union[np.ndarray, Dict]]]:
        mol = self.smiles_to_mol(smiles)
        if mol is None:
            return None
        features = {}
        if self.config.use_rdkit_descriptors:
            features['rdkit_descriptors'] = self.calculate_rdkit_descriptors(mol)
            features['additional_descriptors'] = self.calculate_additional_descriptors(mol)
        if self.config.use_fingerprints:
            features['fingerprints'] = self.calculate_fingerprints(mol)
        if self.config.use_smarts_features:
            features['smarts_features'] = self.calculate_smarts_features(mol)
        if self.config.use_graph_features:
            features['graph_features'] = self.calculate_graph_features(mol)
        return features
    
    def features_to_vector(self, features: Dict) -> np.ndarray:
        vectors = []
        if 'rdkit_descriptors' in features:
            vectors.append(features['rdkit_descriptors'])
        if 'additional_descriptors' in features:
            vectors.append(np.array(list(features['additional_descriptors'].values()), dtype=np.float32))
        if 'fingerprints' in features:
            for fp_vector in features['fingerprints'].values():
                vectors.append(fp_vector)
        if 'smarts_features' in features:
            vectors.append(features['smarts_features'])
        return np.concatenate(vectors) if vectors else np.array([], dtype=np.float32)

# ==================== 深度学习模型 ====================

if TORCH_AVAILABLE:
    class TransformerPredictor(nn.Module):
        def __init__(self, input_dim: int, config: ModelConfig):
            super().__init__()
            self.config = config
            self.input_projection = nn.Linear(input_dim, config.d_model)
            encoder_layer = nn.TransformerEncoderLayer(d_model=config.d_model, nhead=config.nhead,
                                                       dim_feedforward=config.dim_feedforward,
                                                       dropout=config.dropout, batch_first=True)
            self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=config.num_layers)
            self.fc = nn.Sequential(nn.Linear(config.d_model, config.d_model // 2), nn.ReLU(),
                                   nn.Dropout(config.dropout), nn.Linear(config.d_model // 2, config.d_model // 4),
                                   nn.ReLU(), nn.Dropout(config.dropout), nn.Linear(config.d_model // 4, 1))
        
        def forward(self, x):
            x = self.input_projection(x).unsqueeze(1)
            x = self.transformer_encoder(x).squeeze(1)
            return self.fc(x).squeeze(-1)
    
    class GNNPredictor(nn.Module):
        def __init__(self, atom_feat_dim: int, config: ModelConfig):
            super().__init__()
            self.config = config
            self.gnn_layers = nn.ModuleList()
            in_dim = atom_feat_dim
            for _ in range(config.num_gnn_layers):
                self.gnn_layers.append(nn.Linear(in_dim, config.hidden_dim))
                in_dim = config.hidden_dim
            self.fc = nn.Sequential(nn.Linear(config.hidden_dim, config.hidden_dim // 2), nn.ReLU(),
                                   nn.Dropout(config.dropout), nn.Linear(config.hidden_dim // 2, 1))
        
        def forward(self, atom_features, adjacency):
            x = atom_features
            for gnn_layer in self.gnn_layers:
                x = F.relu(gnn_layer(torch.bmm(adjacency, x)))
            return self.fc(torch.mean(x, dim=1)).squeeze(-1)
    
    class MoleculeDataset(Dataset):
        def __init__(self, features: List[np.ndarray], labels: np.ndarray):
            self.features = [torch.FloatTensor(f) for f in features]
            self.labels = torch.FloatTensor(labels)
        def __len__(self):
            return len(self.labels)
        def __getitem__(self, idx):
            return self.features[idx], self.labels[idx]

# ==================== 主预测器 ====================

class MultiTargetPIC50Predictor:
    """多靶点pIC50预测器主类"""
    
    def __init__(self, target_name: str, config: Optional[ModelConfig] = None):
        self.target_name = target_name
        self.config = config or ModelConfig()
        self.logger = logging.getLogger(__name__)
        self.feature_extractor = MolecularFeatureExtractor(self.config)
        self.models = {}
        self.scalers = {}
        self.feature_dim = None
        self.training_history = {'metrics': [], 'predictions': []}
    
    def prepare_data(self, smiles_list: List[str], pic50_values: List[float]) -> Tuple:
        self.logger.info(f"Preparing data for {len(smiles_list)} molecules...")
        features_list, valid_indices = [], []
        for i, smiles in enumerate(smiles_list):
            features = self.feature_extractor.extract_all_features(smiles)
            if features is not None:
                feature_vector = self.feature_extractor.features_to_vector(features)
                if len(feature_vector) > 0:
                    features_list.append(feature_vector)
                    valid_indices.append(i)
            if (i + 1) % 100 == 0:
                self.logger.info(f"  Processed {i + 1}/{len(smiles_list)}")
        X = np.array(features_list, dtype=np.float32)
        y = np.array([pic50_values[i] for i in valid_indices], dtype=np.float32)
        self.feature_dim = X.shape[1]
        self.logger.info(f"Feature dimension: {self.feature_dim}")
        return X, y, valid_indices
    
    def train_random_forest(self, X_train, y_train, X_test, y_test):
        self.logger.info("Training Random Forest...")
        rf_model = RandomForestRegressor(n_estimators=self.config.rf_n_estimators,
                                        max_depth=self.config.rf_max_depth, min_samples_split=5,
                                        min_samples_leaf=2, random_state=self.config.random_state,
                                        n_jobs=-1, verbose=1)
        rf_model.fit(X_train, y_train)
        metrics = self._calculate_metrics(y_train, rf_model.predict(X_train),
                                          y_test, rf_model.predict(X_test), "Random Forest")
        self.models['random_forest'] = rf_model
        return metrics
    
    def train_gradient_boosting(self, X_train, y_train, X_test, y_test):
        self.logger.info("Training Gradient Boosting...")
        gb_model = GradientBoostingRegressor(n_estimators=self.config.gb_n_estimators,
                                            learning_rate=self.config.gb_learning_rate,
                                            max_depth=10, min_samples_split=5, min_samples_leaf=2,
                                            random_state=self.config.random_state, verbose=1)
        gb_model.fit(X_train, y_train)
        metrics = self._calculate_metrics(y_train, gb_model.predict(X_train),
                                          y_test, gb_model.predict(X_test), "Gradient Boosting")
        self.models['gradient_boosting'] = gb_model
        return metrics
    
    def train_transformer(self, X_train, y_train, X_test, y_test, epochs=100, batch_size=32):
        if not TORCH_AVAILABLE:
            self.logger.warning("PyTorch not available. Skipping Transformer training.")
            return None
        self.logger.info("Training Transformer...")
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.logger.info(f"Using device: {device}")
        train_loader = DataLoader(MoleculeDataset([x for x in X_train], y_train), batch_size=batch_size, shuffle=True)
        test_loader = DataLoader(MoleculeDataset([x for x in X_test], y_test), batch_size=batch_size, shuffle=False)
        model = TransformerPredictor(X_train.shape[1], self.config).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        criterion = nn.MSELoss()
        best_test_loss, patience, patience_counter = float('inf'), 10, 0
        for epoch in range(epochs):
            model.train()
            train_loss = 0
            for features, labels in train_loader:
                features, labels = features.to(device), labels.to(device)
                optimizer.zero_grad()
                loss = criterion(model(features), labels)
                loss.backward()
                optimizer.step()
                train_loss += loss.item()
            model.eval()
            test_loss, test_preds, test_true = 0, [], []
            with torch.no_grad():
                for features, labels in test_loader:
                    features, labels = features.to(device), labels.to(device)
                    outputs = model(features)
                    test_loss += criterion(outputs, labels).item()
                    test_preds.extend(outputs.cpu().numpy())
                    test_true.extend(labels.cpu().numpy())
            train_loss /= len(train_loader)
            test_loss /= len(test_loader)
            if (epoch + 1) % 10 == 0:
                self.logger.info(f"Epoch {epoch+1}/{epochs} - Train Loss: {train_loss:.4f}, Test Loss: {test_loss:.4f}, R²: {r2_score(test_true, test_preds):.4f}")
            if test_loss < best_test_loss:
                best_test_loss = test_loss
                patience_counter = 0
                torch.save(model.state_dict(), Config.MODEL_DIR / f"{self.target_name}_transformer_best.pth")
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    self.logger.info(f"Early stopping at epoch {epoch+1}")
                    break
        model.load_state_dict(torch.load(Config.MODEL_DIR / f"{self.target_name}_transformer_best.pth"))
        model.eval()
        train_preds, test_preds = [], []
        with torch.no_grad():
            for features, _ in DataLoader(MoleculeDataset([x for x in X_train], y_train), batch_size=batch_size):
                train_preds.extend(model(features.to(device)).cpu().numpy())
            for features, _ in DataLoader(MoleculeDataset([x for x in X_test], y_test), batch_size=batch_size):
                test_preds.extend(model(features.to(device)).cpu().numpy())
        metrics = self._calculate_metrics(y_train, np.array(train_preds), y_test, np.array(test_preds), "Transformer")
        self.models['transformer'] = model
        return metrics
    
    def train_ensemble(self, X_train, y_train, X_test, y_test):
        self.logger.info("Creating ensemble model...")
        if len(self.models) < 2:
            self.logger.warning("Not enough models for ensemble. Training base models first.")
            return None
        train_preds, test_preds = [], []
        for name, model in self.models.items():
            if name != 'ensemble':
                if isinstance(model, nn.Module):
                    device = next(model.parameters()).device
                    model.eval()
                    with torch.no_grad():
                        train_preds.append(model(torch.FloatTensor(X_train).to(device)).cpu().numpy())
                        test_preds.append(model(torch.FloatTensor(X_test).to(device)).cpu().numpy())
                else:
                    train_preds.append(model.predict(X_train))
                    test_preds.append(model.predict(X_test))
        metrics = self._calculate_metrics(y_train, np.mean(train_preds, axis=0),
                                         y_test, np.mean(test_preds, axis=0), "Ensemble")
        self.models['ensemble'] = {'method': 'average', 'models': list(self.models.keys())}
        return metrics
    
    def _calculate_metrics(self, y_train, train_pred, y_test, test_pred, model_name):
        metrics = {
            'model': model_name, 'train_r2': r2_score(y_train, train_pred),
            'test_r2': r2_score(y_test, test_pred),
            'train_rmse': np.sqrt(mean_squared_error(y_train, train_pred)),
            'test_rmse': np.sqrt(mean_squared_error(y_test, test_pred)),
            'train_mae': mean_absolute_error(y_train, train_pred),
            'test_mae': mean_absolute_error(y_test, test_pred),
        }
        self.logger.info(f"\n{model_name} Metrics:")
        self.logger.info(f"  Train R²: {metrics['train_r2']:.4f}, Test R²: {metrics['test_r2']:.4f}")
        self.logger.info(f"  Train RMSE: {metrics['train_rmse']:.4f}, Test RMSE: {metrics['test_rmse']:.4f}")
        self.logger.info(f"  Train MAE: {metrics['train_mae']:.4f}, Test MAE: {metrics['test_mae']:.4f}")
        return metrics
    
    def train(self, smiles_list: List[str], pic50_values: List[float]):
        self.logger.info(f"\n{'='*80}")
        self.logger.info(f"Training Multi-Target pIC50 Predictor for {self.target_name}")
        self.logger.info(f"{'='*80}\n")
        X, y, valid_indices = self.prepare_data(smiles_list, pic50_values)
        self.logger.info(f"Valid samples: {len(X)}")
        self.logger.info(f"Feature dimension: {X.shape[1]}")
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        self.scalers['features'] = scaler
        X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=self.config.test_size,
                                                            random_state=self.config.random_state)
        self.logger.info(f"Train size: {len(X_train)}, Test size: {len(X_test)}")
        all_metrics = []
        all_metrics.append(self.train_random_forest(X_train, y_train, X_test, y_test))
        all_metrics.append(self.train_gradient_boosting(X_train, y_train, X_test, y_test))
        if self.config.use_transformer and TORCH_AVAILABLE:
            transformer_metrics = self.train_transformer(X_train, y_train, X_test, y_test)
            if transformer_metrics:
                all_metrics.append(transformer_metrics)
        if self.config.use_ensemble:
            ensemble_metrics = self.train_ensemble(X_train, y_train, X_test, y_test)
            if ensemble_metrics:
                all_metrics.append(ensemble_metrics)
        self.training_history['metrics'] = all_metrics
        self.training_history['timestamp'] = datetime.now().isoformat()
        self.training_history['n_samples'] = len(X)
        self.training_history['n_features'] = X.shape[1]
        return all_metrics
    
    def predict(self, smiles: str, model_name: str = 'random_forest', return_uncertainty: bool = False) -> Dict:
        """预测单个分子 - 修复版"""
        features = self.feature_extractor.extract_all_features(smiles)
        if features is None:
            return None
        feature_vector = self.feature_extractor.features_to_vector(features)
        if len(feature_vector) == 0:
            return None
        if 'features' in self.scalers:
            feature_vector = self.scalers['features'].transform([feature_vector])[0]
        
        model = self.models.get(model_name)
        if model is None:
            self.logger.error(f"Model {model_name} not found")
            return None
        
        # 关键修复：先判断是否是字典类型
        if isinstance(model, dict):
            predictions = []
            for base_model_name in ['random_forest', 'gradient_boosting', 'transformer']:
                if base_model_name in self.models:
                    base_model = self.models[base_model_name]
                    if isinstance(base_model, dict):
                        continue
                    if TORCH_AVAILABLE and isinstance(base_model, nn.Module):
                        device = next(base_model.parameters()).device
                        base_model.eval()
                        with torch.no_grad():
                            pred = base_model(torch.FloatTensor([feature_vector]).to(device)).cpu().numpy()[0]
                    else:
                        pred = base_model.predict([feature_vector])[0]
                    predictions.append(pred)
            prediction = float(np.mean(predictions)) if predictions else None
            if prediction is None:
                self.logger.error("No base models available for ensemble")
                return None
        elif TORCH_AVAILABLE and isinstance(model, nn.Module):
            device = next(model.parameters()).device
            model.eval()
            with torch.no_grad():
                prediction = float(model(torch.FloatTensor([feature_vector]).to(device)).cpu().numpy()[0])
        else:
            prediction = float(model.predict([feature_vector])[0])
        
        result = {'smiles': smiles, 'target': self.target_name, 'model': model_name,
                 'pIC50': prediction, 'timestamp': datetime.now().isoformat()}
        if return_uncertainty and self.config.uncertainty_estimation:
            result['uncertainty'] = self._estimate_uncertainty(feature_vector, model_name)
        return result
    
    def predict_batch(self, smiles_list: List[str], model_name: str = 'random_forest') -> List[Dict]:
        """批量预测 - 修复版"""
        results = []
        self.logger.info(f"Batch prediction for {len(smiles_list)} molecules...")
        features_list, valid_smiles = [], []
        for i, smiles in enumerate(smiles_list):
            features = self.feature_extractor.extract_all_features(smiles)
            if features is not None:
                feature_vector = self.feature_extractor.features_to_vector(features)
                if len(feature_vector) > 0:
                    features_list.append(feature_vector)
                    valid_smiles.append(smiles)
            if (i + 1) % 100 == 0:
                self.logger.info(f"  Extracted features {i + 1}/{len(smiles_list)}")
        if not features_list:
            self.logger.warning("No valid features extracted")
            return results
        X = np.array(features_list)
        if 'features' in self.scalers:
            X = self.scalers['features'].transform(X)
        
        model = self.models.get(model_name)
        if model is None:
            self.logger.error(f"Model {model_name} not found")
            return results
        
        # 关键修复：先判断是否是字典类型
        if isinstance(model, dict):
            predictions_list = []
            for base_model_name in ['random_forest', 'gradient_boosting', 'transformer']:
                if base_model_name in self.models:
                    base_model = self.models[base_model_name]
                    if isinstance(base_model, dict):
                        continue
                    if TORCH_AVAILABLE and isinstance(base_model, nn.Module):
                        device = next(base_model.parameters()).device
                        base_model.eval()
                        with torch.no_grad():
                            preds = base_model(torch.FloatTensor(X).to(device)).cpu().numpy()
                    else:
                        preds = base_model.predict(X)
                    predictions_list.append(preds)
            predictions = np.mean(predictions_list, axis=0) if predictions_list else None
            if predictions is None:
                self.logger.error("No base models available for ensemble")
                return results
        elif TORCH_AVAILABLE and isinstance(model, nn.Module):
            device = next(model.parameters()).device
            model.eval()
            with torch.no_grad():
                predictions = model(torch.FloatTensor(X).to(device)).cpu().numpy()
        else:
            predictions = model.predict(X)
        
        for smiles, pred in zip(valid_smiles, predictions):
            results.append({'smiles': smiles, 'target': self.target_name, 'model': model_name,
                          'pIC50': float(pred), 'timestamp': datetime.now().isoformat()})
        self.logger.info(f"Batch prediction completed: {len(results)} molecules")
        return results
    
    def _estimate_uncertainty(self, feature_vector: np.ndarray, model_name: str) -> Dict:
        """估计预测不确定性 - 修复版"""
        model = self.models.get(model_name)
        uncertainties = {}
        
        if isinstance(model, dict):
            predictions = []
            for base_model_name in ['random_forest', 'gradient_boosting', 'transformer']:
                if base_model_name in self.models:
                    base_model = self.models[base_model_name]
                    if isinstance(base_model, dict):
                        continue
                    if TORCH_AVAILABLE and isinstance(base_model, nn.Module):
                        device = next(base_model.parameters()).device
                        base_model.eval()
                        with torch.no_grad():
                            pred = base_model(torch.FloatTensor([feature_vector]).to(device)).cpu().numpy()[0]
                    else:
                        pred = base_model.predict([feature_vector])[0]
                    predictions.append(float(pred))
            if len(predictions) > 1:
                uncertainties = {'std': float(np.std(predictions)), 'mean': float(np.mean(predictions)),
                               'min': float(np.min(predictions)), 'max': float(np.max(predictions)),
                               'n_models': len(predictions), 'method': 'ensemble_variance'}
        elif isinstance(model, RandomForestRegressor):
            tree_predictions = np.array([tree.predict([feature_vector])[0] for tree in model.estimators_])
            uncertainties = {'std': float(np.std(tree_predictions)), 'mean': float(np.mean(tree_predictions)),
                           'method': 'tree_variance'}
        return uncertainties
    
    def save(self, filepath: Optional[Path] = None):
        if filepath is None:
            filepath = Config.MODEL_DIR / f"{self.target_name}_complete_model.pkl"
        filepath.parent.mkdir(parents=True, exist_ok=True)
        save_dict = {'target_name': self.target_name, 'config': asdict(self.config),
                    'feature_dim': self.feature_dim, 'scalers': self.scalers,
                    'training_history': self.training_history, 'models': {}}
        for name, model in self.models.items():
            if not isinstance(model, nn.Module):
                save_dict['models'][name] = model
        joblib.dump(save_dict, filepath)
        for name, model in self.models.items():
            if isinstance(model, nn.Module):
                torch_path = Config.MODEL_DIR / f"{self.target_name}_{name}.pth"
                torch_path.parent.mkdir(parents=True, exist_ok=True)
                torch.save(model.state_dict(), torch_path)
        self.logger.info(f"Model saved: {filepath}")
    
    @classmethod
    def load(cls, filepath: Path):
        save_dict = joblib.load(filepath)
        config = ModelConfig(**save_dict['config'])
        predictor = cls(save_dict['target_name'], config)
        predictor.feature_dim = save_dict['feature_dim']
        predictor.scalers = save_dict['scalers']
        predictor.training_history = save_dict['training_history']
        predictor.models = save_dict['models']
        for name in ['transformer', 'gnn']:
            model_path = Config.MODEL_DIR / f"{predictor.target_name}_{name}.pth"
            if model_path.exists():
                if name == 'transformer':
                    model = TransformerPredictor(predictor.feature_dim, config)
                    model.load_state_dict(torch.load(model_path))
                    model.eval()
                    predictor.models[name] = model
        return predictor

# ==================== 数据管理器 ====================

class DataManager:
    @staticmethod
    def load_from_csv(filepath: Path, smiles_col: str = 'smiles', target_col: str = 'pIC50') -> Tuple[List[str], List[float]]:
        df = pd.read_csv(filepath)
        if smiles_col not in df.columns or target_col not in df.columns:
            raise ValueError(f"Columns '{smiles_col}' or '{target_col}' not found in CSV")
        return df[smiles_col].tolist(), df[target_col].tolist()
    
    @staticmethod
    def load_smiles_file(filepath: Path) -> List[str]:
        with open(filepath, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    
    @staticmethod
    def save_results(results: List[Dict], filepath: Path):
        filepath.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(results).to_csv(filepath, index=False)
        print(f"Results saved: {filepath}")
    
    @staticmethod
    def save_training_report(predictor: MultiTargetPIC50Predictor, filepath: Path):
        filepath.parent.mkdir(parents=True, exist_ok=True)
        report = {'target': predictor.target_name, 'timestamp': datetime.now().isoformat(),
                 'config': asdict(predictor.config), 'metrics': predictor.training_history['metrics'],
                 'n_samples': predictor.training_history.get('n_samples'),
                 'n_features': predictor.training_history.get('n_features')}
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"Training report saved: {filepath}")

# ==================== 主函数 ====================

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Multi-Target pIC50 Predictor')
    parser.add_argument('command', choices=['train', 'predict', 'batch'], help='Command to execute')
    parser.add_argument('--target', type=str, default='DAT', help='Target name (default: DAT)')
    parser.add_argument('--data', type=str, help='Path to training data CSV')
    parser.add_argument('--smiles', type=str, help='SMILES string for prediction')
    parser.add_argument('--input', type=str, help='Input file for batch prediction')
    parser.add_argument('--output', type=str, help='Output file path')
    parser.add_argument('--model', type=str, default='random_forest', help='Model name for prediction')
    args = parser.parse_args()
    
    Config.setup_directories()
    logger = Config.setup_logging()
    print(f"""
╔══════════════════════════════════════════════════════════════════╗
║         Complete Multi-Target pIC50 Predictor                    ║
║         完整多靶点 pIC50 预测系统 - 最终修复版                    ║
╚══════════════════════════════════════════════════════════════════╝
    """)
    
    if args.command == 'train':
        if not args.data:
            print("Error: --data is required for training")
            return
        smiles_list, pic50_values = DataManager.load_from_csv(Path(args.data))
        config = ModelConfig()
        predictor = MultiTargetPIC50Predictor(args.target, config)
        metrics = predictor.train(smiles_list, pic50_values)
        predictor.save()
        DataManager.save_training_report(predictor, Config.OUTPUT_DIR / f"{args.target}_training_report.json")
    elif args.command == 'predict':
        if not args.smiles:
            print("Error: --smiles is required for prediction")
            return
        model_path = Config.MODEL_DIR / f"{args.target}_complete_model.pkl"
        if not model_path.exists():
            print(f"Error: Model not found at {model_path}")
            return
        predictor = MultiTargetPIC50Predictor.load(model_path)
        result = predictor.predict(args.smiles, model_name=args.model, return_uncertainty=True)
        if result:
            print(f"\nPrediction Result:")
            print(f"  SMILES: {result['smiles']}")
            print(f"  pIC50: {result['pIC50']:.2f}")
            if 'uncertainty' in result:
                print(f"  Uncertainty: {result['uncertainty']}")
    elif args.command == 'batch':
        if not args.input:
            print("Error: --input is required for batch prediction")
            return
        model_path = Config.MODEL_DIR / f"{args.target}_complete_model.pkl"
        if not model_path.exists():
            print(f"Error: Model not found at {model_path}")
            return
        predictor = MultiTargetPIC50Predictor.load(model_path)
        input_path = Path(args.input)
        smiles_list = pd.read_csv(input_path)['smiles'].tolist() if input_path.suffix == '.csv' else DataManager.load_smiles_file(input_path)
        results = predictor.predict_batch(smiles_list, model_name=args.model)
        output_path = Path(args.output) if args.output else Config.OUTPUT_DIR / f"{args.target}_predictions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        DataManager.save_results(results, output_path)
        print(f"\nBatch prediction completed: {len(results)} molecules")

if __name__ == "__main__":
    main()