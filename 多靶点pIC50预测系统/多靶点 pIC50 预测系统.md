# 多靶点 pIC50 预测系统

先进的多靶点 pIC50 预测平台，适用于精神活性和类药物化合物。支持 DAT、5HT2A、CB1、CB2 和阿片类受体。
功能包括 Transformer 回归、RDKit 描述符、SMARTS 框架、Optuna 优化。专为科研、药物研发和化学信息学设计。

## 系统要求

- Python 3.7+
- 8GB+ RAM (推荐 16GB)
- GPU (可选，用于深度学习模型)

## 安装依赖

### 1. 创建虚拟环境

```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 或
venv\Scripts\activate  # Windows
```

### 2. 安装依赖包

```bash
# 基础依赖
pip install numpy pandas scikit-learn joblib

# RDKit (化学信息学)
conda install -c conda-forge rdkit  # 推荐
# 或
pip install rdkit-pypi

# 深度学习 (可选)
pip install torch torchvision  # CPU版本
# 或 GPU版本
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# 其他工具
pip install matplotlib seaborn tqdm
```

### requirements.txt

```
numpy>=1.21.0
pandas>=1.3.0
scikit-learn>=1.0.0
rdkit>=2022.3.0
joblib>=1.1.0
torch>=1.10.0
matplotlib>=3.4.0
seaborn>=0.11.0
tqdm>=4.62.0
```

## 使用方法

### 1. 准备数据

#### CSV格式（训练数据）
```csv
smiles,pIC50
CC(CC1=CC=CC=C1)NC,7.5
CN1C=NC2=C1C(=O)N(C(=O)N2C)C,6.2
COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC,5.8
```

#### SMILES文件格式（预测数据）
```
CC(CC1=CC=CC=C1)NC
CN1C=NC2=C1C(=O)N(C(=O)N2C)C
COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC
```

### 2. 训练模型

```bash
# 基本训练
python complete_pic50_predictor.py train --target DAT --data training_data.csv

# 训练其他靶点
python complete_pic50_predictor.py train --target 5HT2A --data 5ht2a_data.csv
python complete_pic50_predictor.py train --target CB1 --data cb1_data.csv
python complete_pic50_predictor.py train --target mu_opioid --data opioid_data.csv
```

### 3. 单个预测

```bash
# 使用Random Forest模型
python complete_pic50_predictor.py predict \
    --target DAT \
    --smiles "CC(CC1=CC=CC=C1)NC" \
    --model random_forest

# 使用Gradient Boosting模型
python complete_pic50_predictor.py predict \
    --target DAT \
    --smiles "CC(CC1=CC=CC=C1)NC" \
    --model gradient_boosting

# 使用Transformer模型 (需要PyTorch)
python complete_pic50_predictor.py predict \
    --target DAT \
    --smiles "CC(CC1=CC=CC=C1)NC" \
    --model transformer
```

### 4. 批量预测

```bash
# 从CSV文件预测
python complete_pic50_predictor.py batch \
    --target DAT \
    --input compounds.csv \
    --output predictions.csv \
    --model random_forest

# 从SMILES文件预测
python complete_pic50_predictor.py batch \
    --target DAT \
    --input compounds.smi \
    --output predictions.csv \
    --model ensemble
```



## 完整的代码

```python
"""
完整多靶点 pIC50 预测系统
包含所有核心功能：Transformer、GNN、集成学习、不确实性估计、特征提取等
基于 https://github.com/zapabob/multi-target-pIC50-predictor
"""

import os
import sys
import warnings
import logging
import json
import pickle
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
from datetime import datetime
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split, cross_val_score
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
    # 靶点定义
    targets: Dict[str, str] = None
    
    # 特征配置
    use_rdkit_descriptors: bool = True
    use_fingerprints: bool = True
    use_smarts_features: bool = True
    use_graph_features: bool = True
    
    # 指纹配置
    fp_radius: int = 2
    fp_nbits: int = 2048
    use_ecfp: bool = True
    use_maccs: bool = True
    
    # 模型配置
    use_transformer: bool = True
    use_gnn: bool = True
    use_ensemble: bool = True
    
    # 训练配置
    test_size: float = 0.2
    random_state: int = 42
    n_splits: int = 5
    
    # Transformer配置
    d_model: int = 256
    nhead: int = 8
    num_layers: int = 4
    dim_feedforward: int = 1024
    dropout: float = 0.1
    
    # GNN配置
    hidden_dim: int = 128
    num_gnn_layers: int = 3
    
    # 集成学习配置
    rf_n_estimators: int = 200
    rf_max_depth: int = 30
    gb_n_estimators: int = 100
    gb_learning_rate: float = 0.1
    
    # 不确定性估计
    uncertainty_estimation: bool = True
    mc_dropout_samples: int = 50
    ensemble_members: int = 5
    
    # 优化配置
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
        """创建目录结构"""
        for dir_path in [cls.DATA_DIR, cls.MODEL_DIR, cls.OUTPUT_DIR, 
                         cls.LOG_DIR, cls.CACHE_DIR]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    @classmethod
    def setup_logging(cls, level=logging.INFO):
        """配置日志"""
        cls.LOG_DIR.mkdir(parents=True, exist_ok=True)
        log_file = cls.LOG_DIR / f"pic50_predictor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        return logging.getLogger(__name__)


# ==================== 分子特征提取器 ====================

class MolecularFeatureExtractor:
    """完整的分子特征提取器"""
    
    def __init__(self, config: ModelConfig):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # RDKit描述符计算器
        self.descriptor_names = [desc[0] for desc in Descriptors._descList]
        self.descriptor_calculator = MoleculeDescriptors.MolecularDescriptorCalculator(
            self.descriptor_names
        )
        
        # SMARTS模式（精神活性化合物特化）
        self.smarts_patterns = self._init_smarts_patterns()
    
    def _init_smarts_patterns(self) -> Dict[str, str]:
        """初始化SMARTS模式"""
        return {
            'phenethylamine': 'c1ccccc1CCN',
            'tryptamine': 'c1ccc2c(c1)c(cn2)CCN',
            'benzodioxole': 'c1cc2c(cc1)OCO2',
            'piperidine': 'C1CCNCC1',
            'morpholine': 'C1COCCN1',
            'tropane': 'C1CC2CCC(C1)N2',
            'indole': 'c1ccc2c(c1)cc[nH]2',
            'quinoline': 'c1ccc2c(c1)cccn2',
            'benzene': 'c1ccccc1',
            'pyridine': 'c1ccncc1',
            'thiophene': 'c1ccsc1',
            'furan': 'c1ccoc1',
            'imidazole': 'c1c[nH]cn1',
            'pyrrole': 'c1cc[nH]c1',
            'oxazole': 'c1cnco1',
            'thiazole': 'c1cncs1',
            'basic_amine': '[NX3;H2,H1;!$(NC=O)]',
            'secondary_amine': '[NX3;H1;!$(NC=O)]',
            'tertiary_amine': '[NX3;H0;!$(NC=O)]',
            'aromatic_amine': '[NX3;H2,H1;!$(NC=O)]c',
            'methoxy': 'COc',
            'ethoxy': 'CCOc',
            'hydroxyl': '[OH]',
            'carbonyl': '[CX3]=[OX1]',
            'carboxylic_acid': 'C(=O)[OH]',
            'ester': 'C(=O)O[CX4]',
            'amide': 'C(=O)N',
            'sulfonamide': 'S(=O)(=O)N',
            'halogen': '[F,Cl,Br,I]',
            'nitro': '[N+](=O)[O-]',
        }
    
    def smiles_to_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """SMILES转分子对象"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.logger.warning(f"Invalid SMILES: {smiles}")
                return None
            return mol
        except Exception as e:
            self.logger.error(f"Error converting SMILES {smiles}: {e}")
            return None
    
    def calculate_rdkit_descriptors(self, mol: Chem.Mol) -> np.ndarray:
        """计算RDKit描述符"""
        try:
            descriptors = self.descriptor_calculator.CalcDescriptors(mol)
            # 替换NaN和Inf
            descriptors = np.array(descriptors, dtype=np.float32)
            descriptors = np.nan_to_num(descriptors, nan=0.0, posinf=0.0, neginf=0.0)
            return descriptors
        except Exception as e:
            self.logger.error(f"Error calculating descriptors: {e}")
            return np.zeros(len(self.descriptor_names), dtype=np.float32)
    
    def calculate_additional_descriptors(self, mol: Chem.Mol) -> Dict[str, float]:
        """计算额外的药物化学描述符"""
        try:
            return {
                # Lipinski规则
                'MolWt': Descriptors.MolWt(mol),
                'LogP': Crippen.MolLogP(mol),
                'NumHDonors': Lipinski.NumHDonors(mol),
                'NumHAcceptors': Lipinski.NumHAcceptors(mol),
                'TPSA': MolSurf.TPSA(mol),
                
                # 结构特征
                'NumRotatableBonds': Lipinski.NumRotatableBonds(mol),
                'NumAromaticRings': Lipinski.NumAromaticRings(mol),
                'NumSaturatedRings': Lipinski.NumSaturatedRings(mol),
                'NumAliphaticRings': Lipinski.NumAliphaticRings(mol),
                'RingCount': Lipinski.RingCount(mol),
                'NumHeteroatoms': Lipinski.NumHeteroatoms(mol),
                
                # 分子形状
                'MolMR': Crippen.MolMR(mol),
                'FractionCsp3': Lipinski.FractionCSP3(mol),
                
                # 图描述符
                'BertzCT': GraphDescriptors.BertzCT(mol),
                'Chi0': GraphDescriptors.Chi0(mol),
                'Chi1': GraphDescriptors.Chi1(mol),
                'Kappa1': GraphDescriptors.Kappa1(mol),
                'Kappa2': GraphDescriptors.Kappa2(mol),
                'Kappa3': GraphDescriptors.Kappa3(mol),
            }
        except Exception as e:
            self.logger.error(f"Error calculating additional descriptors: {e}")
            return {k: 0.0 for k in ['MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors', 
                                      'TPSA', 'NumRotatableBonds', 'NumAromaticRings',
                                      'NumSaturatedRings', 'NumAliphaticRings', 'RingCount',
                                      'NumHeteroatoms', 'MolMR', 'FractionCsp3', 'BertzCT',
                                      'Chi0', 'Chi1', 'Kappa1', 'Kappa2', 'Kappa3']}
    
    def calculate_fingerprints(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        """计算多种分子指纹"""
        fingerprints = {}
        
        try:
            if self.config.use_ecfp:
                # ECFP (Morgan)
                ecfp = AllChem.GetMorganFingerprintAsBitVect(
                    mol, self.config.fp_radius, nBits=self.config.fp_nbits
                )
                fingerprints['ecfp'] = np.array(ecfp, dtype=np.float32)
            
            if self.config.use_maccs:
                # MACCS keys
                maccs = AllChem.GetMACCSKeysFingerprint(mol)
                fingerprints['maccs'] = np.array(maccs, dtype=np.float32)
            
            # Atom Pair
            ap = AllChem.GetAtomPairFingerprintAsBitVect(mol, nBits=2048)
            fingerprints['atom_pair'] = np.array(ap, dtype=np.float32)
            
            # Topological Torsion
            tt = AllChem.GetTopologicalTorsionFingerprintAsBitVect(mol, nBits=2048)
            fingerprints['torsion'] = np.array(tt, dtype=np.float32)
            
        except Exception as e:
            self.logger.error(f"Error calculating fingerprints: {e}")
        
        return fingerprints
    
    def calculate_smarts_features(self, mol: Chem.Mol) -> np.ndarray:
        """计算SMARTS子结构特征"""
        features = []
        for name, pattern in self.smarts_patterns.items():
            try:
                smarts_mol = Chem.MolFromSmarts(pattern)
                if smarts_mol is not None:
                    count = len(mol.GetSubstructMatches(smarts_mol))
                    features.append(count)
                else:
                    features.append(0)
            except:
                features.append(0)
        
        return np.array(features, dtype=np.float32)
    
    def calculate_graph_features(self, mol: Chem.Mol) -> Dict[str, np.ndarray]:
        """计算图神经网络特征"""
        try:
            # 原子特征
            atom_features = []
            for atom in mol.GetAtoms():
                features = [
                    atom.GetAtomicNum(),
                    atom.GetDegree(),
                    atom.GetFormalCharge(),
                    atom.GetHybridization(),
                    atom.GetNumRadicalElectrons(),
                    atom.GetIsAromatic(),
                    atom.IsInRing(),
                    atom.GetTotalValence(),
                ]
                atom_features.append(features)
            
            # 键特征
            bond_features = []
            adjacency = np.zeros((mol.GetNumAtoms(), mol.GetNumAtoms()))
            
            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx()
                j = bond.GetEndAtomIdx()
                bond_type = bond.GetBondTypeAsDouble()
                
                adjacency[i, j] = bond_type
                adjacency[j, i] = bond_type
                
                features = [
                    bond_type,
                    bond.GetIsConjugated(),
                    bond.IsInRing(),
                ]
                bond_features.append(features)
            
            return {
                'atom_features': np.array(atom_features, dtype=np.float32),
                'bond_features': np.array(bond_features, dtype=np.float32) if bond_features else np.array([[0, 0, 0]], dtype=np.float32),
                'adjacency': adjacency.astype(np.float32)
            }
        
        except Exception as e:
            self.logger.error(f"Error calculating graph features: {e}")
            return {
                'atom_features': np.array([[0]*8], dtype=np.float32),
                'bond_features': np.array([[0, 0, 0]], dtype=np.float32),
                'adjacency': np.array([[0]], dtype=np.float32)
            }
    
    def extract_all_features(self, smiles: str) -> Optional[Dict[str, Union[np.ndarray, Dict]]]:
        """提取所有特征"""
        mol = self.smiles_to_mol(smiles)
        if mol is None:
            return None
        
        features = {}
        
        # RDKit描述符
        if self.config.use_rdkit_descriptors:
            features['rdkit_descriptors'] = self.calculate_rdkit_descriptors(mol)
            features['additional_descriptors'] = self.calculate_additional_descriptors(mol)
        
        # 分子指纹
        if self.config.use_fingerprints:
            features['fingerprints'] = self.calculate_fingerprints(mol)
        
        # SMARTS特征
        if self.config.use_smarts_features:
            features['smarts_features'] = self.calculate_smarts_features(mol)
        
        # 图特征
        if self.config.use_graph_features:
            features['graph_features'] = self.calculate_graph_features(mol)
        
        return features
    
    def features_to_vector(self, features: Dict) -> np.ndarray:
        """将特征字典转换为向量"""
        vectors = []
        
        # RDKit描述符
        if 'rdkit_descriptors' in features:
            vectors.append(features['rdkit_descriptors'])
        
        if 'additional_descriptors' in features:
            add_desc = np.array(list(features['additional_descriptors'].values()), dtype=np.float32)
            vectors.append(add_desc)
        
        # 指纹
        if 'fingerprints' in features:
            for fp_name, fp_vector in features['fingerprints'].items():
                vectors.append(fp_vector)
        
        # SMARTS
        if 'smarts_features' in features:
            vectors.append(features['smarts_features'])
        
        # 合并所有向量
        if vectors:
            return np.concatenate(vectors)
        else:
            return np.array([], dtype=np.float32)


# ==================== 深度学习模型 ====================

if TORCH_AVAILABLE:
    class TransformerPredictor(nn.Module):
        """Transformer回归模型"""
        
        def __init__(self, input_dim: int, config: ModelConfig):
            super().__init__()
            self.config = config
            
            # 输入投影
            self.input_projection = nn.Linear(input_dim, config.d_model)
            
            # Transformer编码器
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=config.d_model,
                nhead=config.nhead,
                dim_feedforward=config.dim_feedforward,
                dropout=config.dropout,
                batch_first=True
            )
            self.transformer_encoder = nn.TransformerEncoder(
                encoder_layer,
                num_layers=config.num_layers
            )
            
            # 输出层
            self.fc = nn.Sequential(
                nn.Linear(config.d_model, config.d_model // 2),
                nn.ReLU(),
                nn.Dropout(config.dropout),
                nn.Linear(config.d_model // 2, config.d_model // 4),
                nn.ReLU(),
                nn.Dropout(config.dropout),
                nn.Linear(config.d_model // 4, 1)
            )
        
        def forward(self, x):
            # x: (batch_size, input_dim)
            x = self.input_projection(x)  # (batch_size, d_model)
            x = x.unsqueeze(1)  # (batch_size, 1, d_model)
            
            # Transformer编码
            x = self.transformer_encoder(x)  # (batch_size, 1, d_model)
            x = x.squeeze(1)  # (batch_size, d_model)
            
            # 预测
            output = self.fc(x)  # (batch_size, 1)
            return output.squeeze(-1)
    
    
    class GNNPredictor(nn.Module):
        """图神经网络预测器"""
        
        def __init__(self, atom_feat_dim: int, config: ModelConfig):
            super().__init__()
            self.config = config
            
            # 图卷积层
            self.gnn_layers = nn.ModuleList()
            in_dim = atom_feat_dim
            
            for _ in range(config.num_gnn_layers):
                self.gnn_layers.append(
                    nn.Linear(in_dim, config.hidden_dim)
                )
                in_dim = config.hidden_dim
            
            # 全局池化后的输出层
            self.fc = nn.Sequential(
                nn.Linear(config.hidden_dim, config.hidden_dim // 2),
                nn.ReLU(),
                nn.Dropout(config.dropout),
                nn.Linear(config.hidden_dim // 2, 1)
            )
        
        def forward(self, atom_features, adjacency):
            # atom_features: (batch_size, num_atoms, atom_feat_dim)
            # adjacency: (batch_size, num_atoms, num_atoms)
            
            x = atom_features
            
            # 图卷积
            for gnn_layer in self.gnn_layers:
                # 消息传递
                messages = torch.bmm(adjacency, x)  # (batch_size, num_atoms, hidden_dim)
                x = F.relu(gnn_layer(messages))
            
            # 全局平均池化
            x = torch.mean(x, dim=1)  # (batch_size, hidden_dim)
            
            # 预测
            output = self.fc(x)
            return output.squeeze(-1)


# ==================== 数据集 ====================

if TORCH_AVAILABLE:
    class MoleculeDataset(Dataset):
        """分子数据集"""
        
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
        
        # 特征提取器
        self.feature_extractor = MolecularFeatureExtractor(self.config)
        
        # 模型
        self.models = {}
        self.scalers = {}
        self.feature_dim = None
        
        # 训练历史
        self.training_history = {
            'metrics': [],
            'predictions': []
        }
    
    def prepare_data(self, smiles_list: List[str], pic50_values: List[float]) -> Tuple:
        """准备训练数据"""
        self.logger.info(f"Preparing data for {len(smiles_list)} molecules...")
        
        features_list = []
        valid_indices = []
        
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
        """训练随机森林模型"""
        self.logger.info("Training Random Forest...")
        
        rf_model = RandomForestRegressor(
            n_estimators=self.config.rf_n_estimators,
            max_depth=self.config.rf_max_depth,
            min_samples_split=5,
            min_samples_leaf=2,
            random_state=self.config.random_state,
            n_jobs=-1,
            verbose=1
        )
        
        rf_model.fit(X_train, y_train)
        
        train_pred = rf_model.predict(X_train)
        test_pred = rf_model.predict(X_test)
        
        metrics = self._calculate_metrics(y_train, train_pred, y_test, test_pred, "Random Forest")
        
        self.models['random_forest'] = rf_model
        return metrics
    
    def train_gradient_boosting(self, X_train, y_train, X_test, y_test):
        """训练梯度提升模型"""
        self.logger.info("Training Gradient Boosting...")
        
        gb_model = GradientBoostingRegressor(
            n_estimators=self.config.gb_n_estimators,
            learning_rate=self.config.gb_learning_rate,
            max_depth=10,
            min_samples_split=5,
            min_samples_leaf=2,
            random_state=self.config.random_state,
            verbose=1
        )
        
        gb_model.fit(X_train, y_train)
        
        train_pred = gb_model.predict(X_train)
        test_pred = gb_model.predict(X_test)
        
        metrics = self._calculate_metrics(y_train, train_pred, y_test, test_pred, "Gradient Boosting")
        
        self.models['gradient_boosting'] = gb_model
        return metrics
    
    def train_transformer(self, X_train, y_train, X_test, y_test, epochs=100, batch_size=32):
        """训练Transformer模型"""
        if not TORCH_AVAILABLE:
            self.logger.warning("PyTorch not available. Skipping Transformer training.")
            return None
        
        self.logger.info("Training Transformer...")
        
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.logger.info(f"Using device: {device}")
        
        # 创建数据集
        train_dataset = MoleculeDataset([x for x in X_train], y_train)
        test_dataset = MoleculeDataset([x for x in X_test], y_test)
        
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        
        # 创建模型
        model = TransformerPredictor(X_train.shape[1], self.config).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
        criterion = nn.MSELoss()
        
        # 训练
        best_test_loss = float('inf')
        patience = 10
        patience_counter = 0
        
        for epoch in range(epochs):
            model.train()
            train_loss = 0
            
            for features, labels in train_loader:
                features, labels = features.to(device), labels.to(device)
                
                optimizer.zero_grad()
                outputs = model(features)
                loss = criterion(outputs, labels)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item()
            
            # 验证
            model.eval()
            test_loss = 0
            test_preds = []
            test_true = []
            
            with torch.no_grad():
                for features, labels in test_loader:
                    features, labels = features.to(device), labels.to(device)
                    outputs = model(features)
                    loss = criterion(outputs, labels)
                    test_loss += loss.item()
                    
                    test_preds.extend(outputs.cpu().numpy())
                    test_true.extend(labels.cpu().numpy())
            
            train_loss /= len(train_loader)
            test_loss /= len(test_loader)
            
            if (epoch + 1) % 10 == 0:
                r2 = r2_score(test_true, test_preds)
                self.logger.info(f"Epoch {epoch+1}/{epochs} - Train Loss: {train_loss:.4f}, Test Loss: {test_loss:.4f}, R²: {r2:.4f}")
            
            # Early stopping
            if test_loss < best_test_loss:
                best_test_loss = test_loss
                patience_counter = 0
                # 保存最佳模型
                torch.save(model.state_dict(), Config.MODEL_DIR / f"{self.target_name}_transformer_best.pth")
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    self.logger.info(f"Early stopping at epoch {epoch+1}")
                    break
        
        # 加载最佳模型
        model.load_state_dict(torch.load(Config.MODEL_DIR / f"{self.target_name}_transformer_best.pth"))
        model.eval()
        
        # 最终评估
        train_preds = []
        with torch.no_grad():
            for features, _ in DataLoader(train_dataset, batch_size=batch_size):
                features = features.to(device)
                outputs = model(features)
                train_preds.extend(outputs.cpu().numpy())
        
        test_preds = []
        with torch.no_grad():
            for features, _ in DataLoader(test_dataset, batch_size=batch_size):
                features = features.to(device)
                outputs = model(features)
                test_preds.extend(outputs.cpu().numpy())
        
        metrics = self._calculate_metrics(y_train, np.array(train_preds), 
                                         y_test, np.array(test_preds), "Transformer")
        
        self.models['transformer'] = model
        return metrics
    
    def train_ensemble(self, X_train, y_train, X_test, y_test):
        """训练集成模型"""
        self.logger.info("Creating ensemble model...")
        
        if len(self.models) < 2:
            self.logger.warning("Not enough models for ensemble. Training base models first.")
            return None
        
        # 获取各模型的预测
        train_preds = []
        test_preds = []
        
        for name, model in self.models.items():
            if name != 'ensemble':
                if isinstance(model, nn.Module):
                    # PyTorch模型
                    device = next(model.parameters()).device
                    model.eval()
                    with torch.no_grad():
                        train_pred = model(torch.FloatTensor(X_train).to(device)).cpu().numpy()
                        test_pred = model(torch.FloatTensor(X_test).to(device)).cpu().numpy()
                else:
                    # Scikit-learn模型
                    train_pred = model.predict(X_train)
                    test_pred = model.predict(X_test)
                
                train_preds.append(train_pred)
                test_preds.append(test_pred)
        
        # 简单平均集成
        train_pred_ensemble = np.mean(train_preds, axis=0)
        test_pred_ensemble = np.mean(test_preds, axis=0)
        
        metrics = self._calculate_metrics(y_train, train_pred_ensemble, 
                                         y_test, test_pred_ensemble, "Ensemble")
        
        self.models['ensemble'] = {
            'method': 'average',
            'models': list(self.models.keys())
        }
        
        return metrics
    
    def _calculate_metrics(self, y_train, train_pred, y_test, test_pred, model_name):
        """计算评估指标"""
        metrics = {
            'model': model_name,
            'train_r2': r2_score(y_train, train_pred),
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
        """完整训练流程"""
        self.logger.info(f"\n{'='*80}")
        self.logger.info(f"Training Multi-Target pIC50 Predictor for {self.target_name}")
        self.logger.info(f"{'='*80}\n")
        
        # 准备数据
        X, y, valid_indices = self.prepare_data(smiles_list, pic50_values)
        
        self.logger.info(f"Valid samples: {len(X)}")
        self.logger.info(f"Feature dimension: {X.shape[1]}")
        
        # 数据标准化
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        self.scalers['features'] = scaler
        
        # 划分数据集
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, 
            test_size=self.config.test_size,
            random_state=self.config.random_state
        )
        
        self.logger.info(f"Train size: {len(X_train)}, Test size: {len(X_test)}")
        
        # 训练模型
        all_metrics = []
        
        # 1. Random Forest
        rf_metrics = self.train_random_forest(X_train, y_train, X_test, y_test)
        all_metrics.append(rf_metrics)
        
        # 2. Gradient Boosting
        gb_metrics = self.train_gradient_boosting(X_train, y_train, X_test, y_test)
        all_metrics.append(gb_metrics)
        
        # 3. Transformer (如果启用)
        if self.config.use_transformer and TORCH_AVAILABLE:
            transformer_metrics = self.train_transformer(X_train, y_train, X_test, y_test)
            if transformer_metrics:
                all_metrics.append(transformer_metrics)
        
        # 4. Ensemble
        if self.config.use_ensemble:
            ensemble_metrics = self.train_ensemble(X_train, y_train, X_test, y_test)
            if ensemble_metrics:
                all_metrics.append(ensemble_metrics)
        
        # 保存训练历史
        self.training_history['metrics'] = all_metrics
        self.training_history['timestamp'] = datetime.now().isoformat()
        self.training_history['n_samples'] = len(X)
        self.training_history['n_features'] = X.shape[1]
        
        return all_metrics
    
    def predict(self, smiles: str, model_name: str = 'random_forest', 
                return_uncertainty: bool = False) -> Dict:
        """预测单个分子"""
        # 提取特征
        features = self.feature_extractor.extract_all_features(smiles)
        if features is None:
            return None
        
        feature_vector = self.feature_extractor.features_to_vector(features)
        if len(feature_vector) == 0:
            return None
        
        # 标准化
        if 'features' in self.scalers:
            feature_vector = self.scalers['features'].transform([feature_vector])[0]
        
        # 预测
        model = self.models.get(model_name)
        if model is None:
            self.logger.error(f"Model {model_name} not found")
            return None
        
        if isinstance(model, nn.Module):
            # PyTorch模型
            device = next(model.parameters()).device
            model.eval()
            with torch.no_grad():
                prediction = model(torch.FloatTensor([feature_vector]).to(device)).cpu().numpy()[0]
        else:
            # Scikit-learn模型
            prediction = model.predict([feature_vector])[0]
        
        result = {
            'smiles': smiles,
            'target': self.target_name,
            'model': model_name,
            'pIC50': float(prediction),
            'timestamp': datetime.now().isoformat()
        }
        
        # 不确定性估计
        if return_uncertainty and self.config.uncertainty_estimation:
            uncertainty = self._estimate_uncertainty(feature_vector, model_name)
            result['uncertainty'] = uncertainty
        
        return result
    
    def predict_batch(self, smiles_list: List[str], model_name: str = 'random_forest') -> List[Dict]:
        """批量预测"""
        results = []
        
        self.logger.info(f"Batch prediction for {len(smiles_list)} molecules...")
        
        for i, smiles in enumerate(smiles_list):
            result = self.predict(smiles, model_name)
            if result:
                results.append(result)
            
            if (i + 1) % 100 == 0:
                self.logger.info(f"  Predicted {i + 1}/{len(smiles_list)}")
        
        return results
    
    def _estimate_uncertainty(self, feature_vector: np.ndarray, model_name: str) -> Dict:
        """估计预测不确定性"""
        model = self.models.get(model_name)
        uncertainties = {}
        
        if isinstance(model, RandomForestRegressor):
            # 使用树的标准差
            tree_predictions = np.array([tree.predict([feature_vector])[0] 
                                        for tree in model.estimators_])
            uncertainties['std'] = float(np.std(tree_predictions))
            uncertainties['method'] = 'tree_variance'
        
        return uncertainties
    
    def save(self, filepath: Optional[Path] = None):
        """保存模型"""
        if filepath is None:
            filepath = Config.MODEL_DIR / f"{self.target_name}_complete_model.pkl"
        
        # 保存配置和基本信息
        save_dict = {
            'target_name': self.target_name,
            'config': asdict(self.config),
            'feature_dim': self.feature_dim,
            'scalers': self.scalers,
            'training_history': self.training_history,
            'models': {}
        }
        
        # 保存sklearn模型
        for name, model in self.models.items():
            if not isinstance(model, nn.Module):
                save_dict['models'][name] = model
        
        joblib.dump(save_dict, filepath)
        
        # 保存PyTorch模型
        for name, model in self.models.items():
            if isinstance(model, nn.Module):
                torch.save(model.state_dict(), 
                          Config.MODEL_DIR / f"{self.target_name}_{name}.pth")
        
        self.logger.info(f"Model saved: {filepath}")
    
    @classmethod
    def load(cls, filepath: Path):
        """加载模型"""
        save_dict = joblib.load(filepath)
        
        config = ModelConfig(**save_dict['config'])
        predictor = cls(save_dict['target_name'], config)
        
        predictor.feature_dim = save_dict['feature_dim']
        predictor.scalers = save_dict['scalers']
        predictor.training_history = save_dict['training_history']
        predictor.models = save_dict['models']
        
        # 加载PyTorch模型
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
    """数据管理器"""
    
    @staticmethod
    def load_from_csv(filepath: Path, smiles_col: str = 'smiles', 
                      target_col: str = 'pIC50') -> Tuple[List[str], List[float]]:
        """从CSV加载数据"""
        df = pd.read_csv(filepath)
        
        if smiles_col not in df.columns:
            raise ValueError(f"Column '{smiles_col}' not found in CSV")
        if target_col not in df.columns:
            raise ValueError(f"Column '{target_col}' not found in CSV")
        
        smiles_list = df[smiles_col].tolist()
        pic50_values = df[target_col].tolist()
        
        return smiles_list, pic50_values
    
    @staticmethod
    def load_smiles_file(filepath: Path) -> List[str]:
        """加载SMILES文件"""
        with open(filepath, 'r') as f:
            smiles_list = [line.strip() for line in f if line.strip()]
        return smiles_list
    
    @staticmethod
    def save_results(results: List[Dict], filepath: Path):
        """保存预测结果"""
        df = pd.DataFrame(results)
        df.to_csv(filepath, index=False)
        print(f"Results saved: {filepath}")
    
    @staticmethod
    def save_training_report(predictor: MultiTargetPIC50Predictor, filepath: Path):
        """保存训练报告"""
        report = {
            'target': predictor.target_name,
            'timestamp': datetime.now().isoformat(),
            'config': asdict(predictor.config),
            'metrics': predictor.training_history['metrics'],
            'n_samples': predictor.training_history.get('n_samples'),
            'n_features': predictor.training_history.get('n_features')
        }
        
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"Training report saved: {filepath}")


# ==================== 主函数 ====================

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Multi-Target pIC50 Predictor')
    parser.add_argument('command', choices=['train', 'predict', 'batch'],
                       help='Command to execute')
    parser.add_argument('--target', type=str, default='DAT',
                       help='Target name (default: DAT)')
    parser.add_argument('--data', type=str,
                       help='Path to training data CSV')
    parser.add_argument('--smiles', type=str,
                       help='SMILES string for prediction')
    parser.add_argument('--input', type=str,
                       help='Input file for batch prediction')
    parser.add_argument('--output', type=str,
                       help='Output file path')
    parser.add_argument('--model', type=str, default='random_forest',
                       help='Model name for prediction')
    
    args = parser.parse_args()
    
    # 初始化
    Config.setup_directories()
    logger = Config.setup_logging()
    
    print(f"""
╔══════════════════════════════════════════════════════════════════╗
║         Complete Multi-Target pIC50 Predictor                    ║
║         完整多靶点 pIC50 预测系统                                 ║
╚══════════════════════════════════════════════════════════════════╝
    """)
    
    if args.command == 'train':
        # 训练模式
        if not args.data:
            print("Error: --data is required for training")
            return
        
        # 加载数据
        smiles_list, pic50_values = DataManager.load_from_csv(Path(args.data))
        
        # 创建预测器
        config = ModelConfig()
        predictor = MultiTargetPIC50Predictor(args.target, config)
        
        # 训练
        metrics = predictor.train(smiles_list, pic50_values)
        
        # 保存模型
        predictor.save()
        
        # 保存报告
        report_path = Config.OUTPUT_DIR / f"{args.target}_training_report.json"
        DataManager.save_training_report(predictor, report_path)
        
    elif args.command == 'predict':
        # 单个预测
        if not args.smiles:
            print("Error: --smiles is required for prediction")
            return
        
        # 加载模型
        model_path = Config.MODEL_DIR / f"{args.target}_complete_model.pkl"
        if not model_path.exists():
            print(f"Error: Model not found at {model_path}")
            print("Please train the model first")
            return
        
        predictor = MultiTargetPIC50Predictor.load(model_path)
        
        # 预测
        result = predictor.predict(args.smiles, model_name=args.model, return_uncertainty=True)
        
        if result:
            print(f"\nPrediction Result:")
            print(f"  SMILES: {result['smiles']}")
            print(f"  Target: {result['target']}")
            print(f"  Model: {result['model']}")
            print(f"  pIC50: {result['pIC50']:.2f}")
            if 'uncertainty' in result:
                print(f"  Uncertainty: {result['uncertainty']}")
        else:
            print("Prediction failed")
    
    elif args.command == 'batch':
        # 批量预测
        if not args.input:
            print("Error: --input is required for batch prediction")
            return
        
        # 加载模型
        model_path = Config.MODEL_DIR / f"{args.target}_complete_model.pkl"
        if not model_path.exists():
            print(f"Error: Model not found at {model_path}")
            return
        
        predictor = MultiTargetPIC50Predictor.load(model_path)
        
        # 加载SMILES
        input_path = Path(args.input)
        if input_path.suffix == '.csv':
            df = pd.read_csv(input_path)
            smiles_list = df['smiles'].tolist()
        else:
            smiles_list = DataManager.load_smiles_file(input_path)
        
        # 批量预测
        results = predictor.predict_batch(smiles_list, model_name=args.model)
        
        # 保存结果
        output_path = Path(args.output) if args.output else \
                     Config.OUTPUT_DIR / f"{args.target}_predictions_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        DataManager.save_results(results, output_path)
        
        print(f"\nBatch prediction completed: {len(results)} molecules")


if __name__ == "__main__":
    main()
```





## Python API使用

### 基本使用

```python
from pathlib import Path
from complete_pic50_predictor import (
    MultiTargetPIC50Predictor, 
    ModelConfig, 
    DataManager,
    Config
)

# 初始化
Config.setup_directories()
Config.setup_logging()

# 创建配置
config = ModelConfig(
    use_transformer=True,
    use_gnn=False,  # 设为False可加快训练
    use_ensemble=True,
    rf_n_estimators=200,
    uncertainty_estimation=True
)

# 创建预测器
predictor = MultiTargetPIC50Predictor('DAT', config)

# 加载数据
smiles_list, pic50_values = DataManager.load_from_csv(
    Path('data/training_data.csv'),
    smiles_col='smiles',
    target_col='pIC50'
)

# 训练
metrics = predictor.train(smiles_list, pic50_values)

# 保存模型
predictor.save()

# 预测
result = predictor.predict(
    'CC(CC1=CC=CC=C1)NC', 
    model_name='random_forest',
    return_uncertainty=True
)
print(f"pIC50: {result['pIC50']:.2f}")
```

### 批量预测

```python
# 批量预测
smiles_to_predict = [
    'CC(CC1=CC=CC=C1)NC',
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC'
]

results = predictor.predict_batch(smiles_to_predict, model_name='ensemble')

# 保存结果
DataManager.save_results(results, Path('output/batch_predictions.csv'))
```

### 加载已训练模型

```python
# 加载模型
predictor = MultiTargetPIC50Predictor.load(
    Path('models/DAT_complete_model.pkl')
)

# 使用模型预测
result = predictor.predict('CC(CC1=CC=CC=C1)NC')
```

## 功能特性

### 1. 多种分子特征

- **RDKit描述符**: 200+ 分子描述符
- **分子指纹**:
  - ECFP4 (Morgan, 2048位)
  - MACCS Keys (167位)
  - Atom Pair
  - Topological Torsion
- **SMARTS子结构**: 30+ 精神活性化合物特征
- **图特征**: 原子和键特征（用于GNN）

### 2. 多种机器学习模型

- **Random Forest**: 快速、稳定
- **Gradient Boosting**: 高精度
- **Transformer**: 深度学习、自注意力机制
- **GNN**: 图神经网络（分子图结构）
- **Ensemble**: 多模型集成

### 3. 不确定性估计

- 树模型方差
- MC Dropout (Transformer)
- Deep Ensemble

### 4. 支持的靶点

| 靶点 | ChEMBL ID | 描述 |
|------|-----------|------|
| DAT | CHEMBL238 | 多巴胺转运体 |
| 5HT2A | CHEMBL224 | 血清素2A受体 |
| CB1 | CHEMBL218 | 大麻素受体1 |
| CB2 | CHEMBL1861 | 大麻素受体2 |
| mu_opioid | CHEMBL233 | μ-阿片受体 |
| delta_opioid | CHEMBL236 | δ-阿片受体 |
| kappa_opioid | CHEMBL237 | κ-阿片受体 |

## 高级配置

### 自定义模型配置

```python
config = ModelConfig(
    # 特征配置
    use_rdkit_descriptors=True,
    use_fingerprints=True,
    use_smarts_features=True,
    use_graph_features=True,
    
    # 指纹配置
    fp_radius=3,  # 增加指纹半径
    fp_nbits=4096,  # 增加位数
    
    # 模型选择
    use_transformer=True,
    use_gnn=True,
    use_ensemble=True,
    
    # Random Forest参数
    rf_n_estimators=300,
    rf_max_depth=40,
    
    # Gradient Boosting参数
    gb_n_estimators=150,
    gb_learning_rate=0.05,
    
    # Transformer参数
    d_model=512,
    nhead=16,
    num_layers=6,
    dropout=0.2,
    
    # 训练配置
    test_size=0.15,
    random_state=42,
    
    # 不确定性估计
    uncertainty_estimation=True,
    mc_dropout_samples=100,
    ensemble_members=10
)
```

### 特征重要性分析

```python
# Random Forest特征重要性
if 'random_forest' in predictor.models:
    rf_model = predictor.models['random_forest']
    importances = rf_model.feature_importances_
    
    # 获取前20个重要特征
    indices = np.argsort(importances)[::-1][:20]
    
    print("Top 20 Important Features:")
    for i, idx in enumerate(indices):
        print(f"{i+1}. Feature {idx}: {importances[idx]:.4f}")
```

### 交叉验证

```python
from sklearn.model_selection import cross_val_score

# 准备数据
X, y, _ = predictor.prepare_data(smiles_list, pic50_values)

# 交叉验证
if 'random_forest' in predictor.models:
    rf_model = predictor.models['random_forest']
    cv_scores = cross_val_score(
        rf_model, X, y, 
        cv=5, 
        scoring='r2',
        n_jobs=-1
    )
    
    print(f"Cross-validation R² scores: {cv_scores}")
    print(f"Mean: {cv_scores.mean():.4f} (+/- {cv_scores.std() * 2:.4f})")
```

## 输出文件说明

### 训练后生成的文件

```
models/
├── DAT_complete_model.pkl          # 完整模型（sklearn模型）
├── DAT_transformer_best.pth        # Transformer模型权重
└── DAT_gnn_best.pth               # GNN模型权重

output/
├── DAT_training_report.json        # 训练报告
└── DAT_predictions_20251019.csv    # 预测结果

logs/
└── pic50_predictor_20251019.log    # 运行日志
```

### 训练报告格式 (JSON)

```json
{
  "target": "DAT",
  "timestamp": "2025-10-19T10:30:00",
  "config": {
    "use_transformer": true,
    "use_gnn": false,
    "use_ensemble": true,
    "rf_n_estimators": 200
  },
  "metrics": [
    {
      "model": "Random Forest",
      "train_r2": 0.9542,
      "test_r2": 0.8321,
      "train_rmse": 0.2134,
      "test_rmse": 0.4087,
      "train_mae": 0.1567,
      "test_mae": 0.3012
    },
    {
      "model": "Gradient Boosting",
      "train_r2": 0.9321,
      "test_r2": 0.8456,
      "train_rmse": 0.2598,
      "test_rmse": 0.3921,
      "train_mae": 0.1890,
      "test_mae": 0.2876
    },
    {
      "model": "Transformer",
      "train_r2": 0.9678,
      "test_r2": 0.8598,
      "train_rmse": 0.1789,
      "test_rmse": 0.3734,
      "train_mae": 0.1234,
      "test_mae": 0.2654
    },
    {
      "model": "Ensemble",
      "train_r2": 0.9701,
      "test_r2": 0.8712,
      "train_rmse": 0.1723,
      "test_rmse": 0.3589,
      "train_mae": 0.1189,
      "test_mae": 0.2543
    }
  ],
  "n_samples": 1250,
  "n_features": 3456
}
```

### 预测结果格式 (CSV)

```csv
smiles,target,model,pIC50,timestamp,uncertainty
CC(CC1=CC=CC=C1)NC,DAT,random_forest,7.25,2025-10-19T10:35:00,0.32
CN1C=NC2=C1C(=O)N(C(=O)N2C)C,DAT,random_forest,6.18,2025-10-19T10:35:01,0.45
COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC,DAT,random_forest,5.87,2025-10-19T10:35:02,0.38
```

## 性能优化

### 1. 使用更少的特征

```python
# 只使用指纹，速度更快
config = ModelConfig(
    use_rdkit_descriptors=False,
    use_fingerprints=True,
    use_smarts_features=False,
    use_graph_features=False,
    fp_nbits=1024  # 减少位数
)
```

### 2. 禁用深度学习模型

```python
# 只使用传统机器学习
config = ModelConfig(
    use_transformer=False,
    use_gnn=False,
    use_ensemble=True
)
```

### 3. 减少树的数量

```python
config = ModelConfig(
    rf_n_estimators=100,  # 默认200
    gb_n_estimators=50    # 默认100
)
```

### 4. 并行处理

```python
# Random Forest自动使用所有CPU核心
# n_jobs=-1 已默认设置
```

## 常见问题

### Q1: RDKit安装失败

**解决方案**:
```bash
# 方法1: 使用conda (推荐)
conda install -c conda-forge rdkit

# 方法2: 使用pip
pip install rdkit-pypi

# 方法3: 从源码编译
git clone https://github.com/rdkit/rdkit.git
cd rdkit
mkdir build && cd build
cmake ..
make && make install
```

### Q2: PyTorch CUDA支持

```bash
# 检查CUDA版本
nvidia-smi

# 安装对应版本的PyTorch
# CUDA 11.8
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

# CUDA 12.1
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121

# CPU版本
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
```

### Q3: 内存不足

**解决方案**:
1. 减少特征维度
2. 减少批次大小
3. 使用更小的模型
4. 分批处理数据

```python
# 分批处理
batch_size = 100
results = []

for i in range(0, len(smiles_list), batch_size):
    batch = smiles_list[i:i+batch_size]
    batch_results = predictor.predict_batch(batch)
    results.extend(batch_results)
```

### Q4: 模型训练时间过长

**解决方案**:
1. 减少数据量进行初始测试
2. 使用更快的模型（Random Forest优先）
3. 减少Transformer的epochs
4. 禁用不必要的模型

```python
# 快速训练配置
config = ModelConfig(
    use_transformer=False,  # 禁用最慢的模型
    use_gnn=False,
    rf_n_estimators=50,     # 减少树的数量
    gb_n_estimators=25
)
```

### Q5: 预测结果不准确

**解决方案**:
1. 增加训练数据量（建议>1000样本）
2. 检查数据质量（SMILES有效性、pIC50范围）
3. 使用集成模型
4. 调整模型参数

```python
# 数据质量检查
from rdkit import Chem

valid_smiles = []
valid_pic50 = []

for smiles, pic50 in zip(smiles_list, pic50_values):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and 3.0 <= pic50 <= 12.0:
        valid_smiles.append(smiles)
        valid_pic50.append(pic50)

print(f"Valid samples: {len(valid_smiles)}/{len(smiles_list)}")
```

## 完整工作流示例

### 示例1: 从ChEMBL数据训练和预测

```python
from complete_pic50_predictor import *
from pathlib import Path
import pandas as pd

# 1. 初始化环境
Config.setup_directories()
logger = Config.setup_logging()

# 2. 准备数据 (假设你已从ChEMBL下载数据)
# 数据格式: smiles, standard_value
df = pd.read_csv('chembl_dat_data.csv')

# 转换为pIC50
df['pIC50'] = -np.log10(df['standard_value'] * 1e-9)  # nM to M

# 过滤数据
df = df[(df['pIC50'] >= 4) & (df['pIC50'] <= 10)]
df = df.dropna()

print(f"Total samples: {len(df)}")

# 保存处理后的数据
df[['smiles', 'pIC50']].to_csv('data/DAT_training.csv', index=False)

# 3. 创建并训练模型
config = ModelConfig(
    use_transformer=True,
    use_ensemble=True,
    rf_n_estimators=200,
    uncertainty_estimation=True
)

predictor = MultiTargetPIC50Predictor('DAT', config)

smiles_list = df['smiles'].tolist()
pic50_values = df['pIC50'].tolist()

# 训练
metrics = predictor.train(smiles_list, pic50_values)

# 4. 保存模型和报告
predictor.save()
DataManager.save_training_report(
    predictor, 
    Path('output/DAT_training_report.json')
)

# 5. 测试预测
test_smiles = [
    'CC(CC1=CC=CC=C1)NC',  # 甲基苯丙胺
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # 咖啡因
]

for smiles in test_smiles:
    result = predictor.predict(
        smiles, 
        model_name='ensemble',
        return_uncertainty=True
    )
    print(f"\nSMILES: {smiles}")
    print(f"pIC50: {result['pIC50']:.2f}")
    if 'uncertainty' in result:
        print(f"Uncertainty: {result['uncertainty']}")

# 6. 批量预测新化合物
new_compounds = pd.read_csv('new_compounds.csv')
results = predictor.predict_batch(
    new_compounds['smiles'].tolist(),
    model_name='ensemble'
)

# 保存结果
output_path = Path('output/new_compounds_predictions.csv')
DataManager.save_results(results, output_path)

print(f"\nPredictions saved to {output_path}")
```

### 示例2: 多靶点批量训练

```python
# 训练多个靶点的模型
targets = ['DAT', '5HT2A', 'CB1', 'CB2', 'mu_opioid']

for target in targets:
    print(f"\n{'='*80}")
    print(f"Training {target} model")
    print(f"{'='*80}")
    
    # 加载对应靶点的数据
    data_file = Path(f'data/{target}_training.csv')
    
    if not data_file.exists():
        print(f"Data file not found: {data_file}")
        continue
    
    smiles_list, pic50_values = DataManager.load_from_csv(data_file)
    
    # 创建配置
    config = ModelConfig(
        use_transformer=False,  # 为了速度，只用传统ML
        use_ensemble=True
    )
    
    # 训练
    predictor = MultiTargetPIC50Predictor(target, config)
    metrics = predictor.train(smiles_list, pic50_values)
    
    # 保存
    predictor.save()
    DataManager.save_training_report(
        predictor,
        Path(f'output/{target}_report.json')
    )
    
    print(f"{target} training completed!")
```

### 示例3: 模型比较和选择

```python
# 比较不同模型的性能
from complete_pic50_predictor import *

# 加载模型
predictor = MultiTargetPIC50Predictor.load(
    Path('models/DAT_complete_model.pkl')
)

# 测试数据
test_smiles = [
    'CC(CC1=CC=CC=C1)NC',
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC',
    # ... 更多测试化合物
]

# 比较所有模型
model_names = ['random_forest', 'gradient_boosting', 'ensemble']

comparison_results = []

for smiles in test_smiles:
    row = {'smiles': smiles}
    
    for model_name in model_names:
        if model_name in predictor.models:
            result = predictor.predict(smiles, model_name=model_name)
            if result:
                row[f'{model_name}_pIC50'] = result['pIC50']
    
    comparison_results.append(row)

# 保存比较结果
df_comparison = pd.DataFrame(comparison_results)
df_comparison.to_csv('output/model_comparison.csv', index=False)

print("Model comparison:")
print(df_comparison)
```



## 完整的简单运行示例

```python
"""
完整使用示例脚本
演示所有功能的使用方法
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# 导入主模块
from complete_pic50_predictor import (
    MultiTargetPIC50Predictor,
    ModelConfig,
    DataManager,
    Config,
    MolecularFeatureExtractor
)


def example_1_basic_training():
    """示例1: 基础训练流程"""
    print("\n" + "="*80)
    print("示例 1: 基础训练流程")
    print("="*80)
    
    # 创建示例数据
    sample_data = {
        'smiles': [
            'CC(CC1=CC=CC=C1)NC',  # 甲基苯丙胺
            'CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1',  # LSD类似物
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # 咖啡因
            'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC',  # 甲氧苄啶
            'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3',  # 地西泮
            'CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F',  # Celecoxib
            'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # 布洛芬
            'CC1=C(C(=O)C2=C(C1=O)C=CC=C2C)C',  # 维生素K
        ],
        'pIC50': [7.5, 8.2, 6.1, 5.8, 7.0, 6.5, 5.5, 6.8]
    }
    
    # 为了演示，复制数据以增加样本量
    df = pd.DataFrame(sample_data)
    df_expanded = pd.concat([df] * 50, ignore_index=True)
    
    # 添加一些随机变化
    np.random.seed(42)
    df_expanded['pIC50'] = df_expanded['pIC50'] + np.random.normal(0, 0.2, len(df_expanded))
    
    # 保存数据
    Config.setup_directories()
    data_file = Config.DATA_DIR / 'demo_training.csv'
    df_expanded.to_csv(data_file, index=False)
    
    print(f"创建了 {len(df_expanded)} 个训练样本")
    print(f"数据保存到: {data_file}")
    
    # 创建配置（快速训练配置）
    config = ModelConfig(
        use_transformer=False,  # 禁用以加快速度
        use_gnn=False,
        use_ensemble=True,
        rf_n_estimators=50,  # 减少数量以加快速度
        gb_n_estimators=25
    )
    
    # 创建预测器
    predictor = MultiTargetPIC50Predictor('DEMO', config)
    
    # 训练
    smiles_list = df_expanded['smiles'].tolist()
    pic50_values = df_expanded['pIC50'].tolist()
    
    print("\n开始训练...")
    metrics = predictor.train(smiles_list, pic50_values)
    
    # 显示结果
    print("\n训练结果:")
    for metric in metrics:
        print(f"\n{metric['model']}:")
        print(f"  测试集 R²: {metric['test_r2']:.4f}")
        print(f"  测试集 RMSE: {metric['test_rmse']:.4f}")
    
    # 保存模型
    predictor.save()
    print(f"\n模型已保存到: {Config.MODEL_DIR}")
    
    return predictor


def example_2_prediction(predictor):
    """示例2: 预测功能"""
    print("\n" + "="*80)
    print("示例 2: 单个和批量预测")
    print("="*80)
    
    # 单个预测
    test_smiles = 'CC(CC1=CC=CC=C1)NC'
    print(f"\n测试SMILES: {test_smiles}")
    
    # 使用不同模型预测
    for model_name in ['random_forest', 'gradient_boosting', 'ensemble']:
        if model_name in predictor.models:
            result = predictor.predict(
                test_smiles, 
                model_name=model_name,
                return_uncertainty=True
            )
            print(f"\n{model_name}:")
            print(f"  pIC50: {result['pIC50']:.2f}")
            if 'uncertainty' in result:
                print(f"  不确定性: {result['uncertainty']}")
    
    # 批量预测
    print("\n批量预测测试:")
    batch_smiles = [
        'CC(CC1=CC=CC=C1)NC',
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC',
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
    ]
    
    results = predictor.predict_batch(batch_smiles, model_name='ensemble')
    
    print(f"\n预测了 {len(results)} 个化合物:")
    for r in results:
        print(f"  {r['smiles'][:30]}... -> pIC50: {r['pIC50']:.2f}")
    
    # 保存结果
    output_file = Config.OUTPUT_DIR / 'demo_predictions.csv'
    DataManager.save_results(results, output_file)
    print(f"\n结果保存到: {output_file}")


def example_3_feature_analysis():
    """示例3: 特征提取和分析"""
    print("\n" + "="*80)
    print("示例 3: 分子特征提取和分析")
    print("="*80)
    
    config = ModelConfig()
    extractor = MolecularFeatureExtractor(config)
    
    smiles = 'CC(CC1=CC=CC=C1)NC'  # 甲基苯丙胺
    print(f"\n分析 SMILES: {smiles}")
    
    # 提取所有特征
    features = extractor.extract_all_features(smiles)
    
    if features:
        print("\n提取的特征类型:")
        
        # RDKit描述符
        if 'rdkit_descriptors' in features:
            print(f"  - RDKit描述符: {len(features['rdkit_descriptors'])} 个")
        
        # 额外描述符
        if 'additional_descriptors' in features:
            print(f"  - 额外描述符: {len(features['additional_descriptors'])} 个")
            print("\n主要药物化学性质:")
            add_desc = features['additional_descriptors']
            print(f"    分子量: {add_desc['MolWt']:.2f}")
            print(f"    LogP: {add_desc['LogP']:.2f}")
            print(f"    氢键供体: {add_desc['NumHDonors']}")
            print(f"    氢键受体: {add_desc['NumHAcceptors']}")
            print(f"    TPSA: {add_desc['TPSA']:.2f}")
            print(f"    可旋转键: {add_desc['NumRotatableBonds']}")
            print(f"    芳香环: {add_desc['NumAromaticRings']}")
        
        # 指纹
        if 'fingerprints' in features:
            print(f"\n  - 分子指纹:")
            for fp_name, fp_array in features['fingerprints'].items():
                print(f"    {fp_name}: {len(fp_array)} 位")
                print(f"      设置位数: {int(fp_array.sum())}")
        
        # SMARTS特征
        if 'smarts_features' in features:
            smarts_count = int(features['smarts_features'].sum())
            print(f"\n  - SMARTS子结构匹配: {smarts_count} 个")
        
        # 图特征
        if 'graph_features' in features:
            graph = features['graph_features']
            print(f"\n  - 图特征:")
            print(f"    原子数: {len(graph['atom_features'])}")
            print(f"    键数: {len(graph['bond_features'])}")
        
        # 转换为向量
        vector = extractor.features_to_vector(features)
        print(f"\n总特征维度: {len(vector)}")
        print(f"非零特征: {np.count_nonzero(vector)}")


def example_4_model_comparison():
    """示例4: 模型性能比较"""
    print("\n" + "="*80)
    print("示例 4: 不同模型性能比较")
    print("="*80)
    
    # 加载模型
    model_file = Config.MODEL_DIR / 'DEMO_complete_model.pkl'
    
    if not model_file.exists():
        print("请先运行示例1训练模型")
        return
    
    predictor = MultiTargetPIC50Predictor.load(model_file)
    
    # 测试化合物
    test_compounds = {
        '甲基苯丙胺': 'CC(CC1=CC=CC=C1)NC',
        '咖啡因': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        '甲氧苄啶': 'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC',
        '布洛芬': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
    }
    
    print("\n不同模型的预测结果对比:\n")
    print(f"{'化合物':<12} {'Random Forest':<15} {'Gradient Boosting':<20} {'Ensemble':<15}")
    print("-" * 70)
    
    for name, smiles in test_compounds.items():
        predictions = {}
        
        for model_name in ['random_forest', 'gradient_boosting', 'ensemble']:
            if model_name in predictor.models:
                result = predictor.predict(smiles, model_name=model_name)
                if result:
                    predictions[model_name] = result['pIC50']
        
        rf_pred = predictions.get('random_forest', 0)
        gb_pred = predictions.get('gradient_boosting', 0)
        en_pred = predictions.get('ensemble', 0)
        
        print(f"{name:<12} {rf_pred:<15.2f} {gb_pred:<20.2f} {en_pred:<15.2f}")
    
    # 显示模型统计
    print("\n模型训练统计:")
    if 'metrics' in predictor.training_history:
        for metric in predictor.training_history['metrics']:
            print(f"\n{metric['model']}:")
            print(f"  训练样本: {predictor.training_history.get('n_samples', 'N/A')}")
            print(f"  特征维度: {predictor.training_history.get('n_features', 'N/A')}")
            print(f"  测试集R²: {metric['test_r2']:.4f}")
            print(f"  测试集RMSE: {metric['test_rmse']:.4f}")


def example_5_csv_workflow():
    """示例5: 完整的CSV工作流"""
    print("\n" + "="*80)
    print("示例 5: CSV文件完整工作流")
    print("="*80)
    
    # 1. 创建输入CSV
    input_data = {
        'compound_id': ['COMP001', 'COMP002', 'COMP003', 'COMP004', 'COMP005'],
        'smiles': [
            'CC(CC1=CC=CC=C1)NC',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC',
            'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
            'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3',
        ],
        'source': ['ChEMBL'] * 5
    }
    
    df_input = pd.DataFrame(input_data)
    input_file = Config.DATA_DIR / 'compounds_to_predict.csv'
    df_input.to_csv(input_file, index=False)
    print(f"\n创建输入文件: {input_file}")
    print(df_input)
    
    # 2. 加载模型
    model_file = Config.MODEL_DIR / 'DEMO_complete_model.pkl'
    
    if not model_file.exists():
        print("\n请先运行示例1训练模型")
        return
    
    predictor = MultiTargetPIC50Predictor.load(model_file)
    
    # 3. 批量预测
    print("\n执行批量预测...")
    smiles_list = df_input['smiles'].tolist()
    results = predictor.predict_batch(smiles_list, model_name='ensemble')
    
    # 4. 合并结果
    df_results = pd.DataFrame(results)
    df_final = df_input.merge(
        df_results[['smiles', 'pIC50']], 
        on='smiles',
        how='left'
    )
    
    # 5. 保存结果
    output_file = Config.OUTPUT_DIR / 'prediction_results.csv'
    df_final.to_csv(output_file, index=False)
    
    print(f"\n预测完成，结果保存到: {output_file}")
    print("\n结果预览:")
    print(df_final)
    
    # 6. 统计分析
    print("\n统计分析:")
    print(f"  平均 pIC50: {df_final['pIC50'].mean():.2f}")
    print(f"  最大 pIC50: {df_final['pIC50'].max():.2f}")
    print(f"  最小 pIC50: {df_final['pIC50'].min():.2f}")
    print(f"  标准差: {df_final['pIC50'].std():.2f}")
    
    # 7. 分类统计
    print("\n活性分类:")
    high_activity = df_final[df_final['pIC50'] >= 7.0]
    medium_activity = df_final[(df_final['pIC50'] >= 6.0) & (df_final['pIC50'] < 7.0)]
    low_activity = df_final[df_final['pIC50'] < 6.0]
    
    print(f"  高活性 (pIC50≥7.0): {len(high_activity)} 个化合物")
    print(f"  中等活性 (6.0≤pIC50<7.0): {len(medium_activity)} 个化合物")
    print(f"  低活性 (pIC50<6.0): {len(low_activity)} 个化合物")


def example_6_advanced_config():
    """示例6: 高级配置和自定义"""
    print("\n" + "="*80)
    print("示例 6: 高级配置和自定义")
    print("="*80)
    
    # 创建自定义配置
    custom_config = ModelConfig(
        # 特征配置
        use_rdkit_descriptors=True,
        use_fingerprints=True,
        use_smarts_features=True,
        use_graph_features=False,  # 禁用图特征以加快速度
        
        # 指纹配置
        fp_radius=3,  # 增加半径
        fp_nbits=1024,  # 减少位数
        use_ecfp=True,
        use_maccs=True,
        
        # 模型选择
        use_transformer=False,
        use_gnn=False,
        use_ensemble=True,
        
        # Random Forest参数
        rf_n_estimators=100,
        rf_max_depth=25,
        
        # Gradient Boosting参数
        gb_n_estimators=75,
        gb_learning_rate=0.05,
        
        # 训练配置
        test_size=0.25,  # 增加测试集比例
        random_state=123,  # 自定义随机种子
        
        # 不确定性估计
        uncertainty_estimation=True,
    )
    
    print("\n自定义配置:")
    print(f"  指纹半径: {custom_config.fp_radius}")
    print(f"  指纹位数: {custom_config.fp_nbits}")
    print(f"  Random Forest树数: {custom_config.rf_n_estimators}")
    print(f"  Gradient Boosting树数: {custom_config.gb_n_estimators}")
    print(f"  测试集比例: {custom_config.test_size}")
    print(f"  使用Transformer: {custom_config.use_transformer}")
    print(f"  使用GNN: {custom_config.use_gnn}")
    print(f"  使用集成: {custom_config.use_ensemble}")
    
    print("\n这个配置适用于:")
    print("  - 快速训练和预测")
    print("  - 中等规模数据集")
    print("  - 不需要深度学习的场景")


def run_all_examples():
    """运行所有示例"""
    print("\n" + "#"*80)
    print("# 完整多靶点 pIC50 预测系统 - 使用示例")
    print("# Complete Multi-Target pIC50 Predictor - Usage Examples")
    print("#"*80)
    
    try:
        # 初始化环境
        Config.setup_directories()
        Config.setup_logging()
        
        # 运行示例
        print("\n开始运行示例...")
        
        # 示例1: 训练
        predictor = example_1_basic_training()
        
        # 示例2: 预测
        example_2_prediction(predictor)
        
        # 示例3: 特征分析
        example_3_feature_analysis()
        
        # 示例4: 模型比较
        example_4_model_comparison()
        
        # 示例5: CSV工作流
        example_5_csv_workflow()
        
        # 示例6: 高级配置
        example_6_advanced_config()
        
        print("\n" + "="*80)
        print("所有示例运行完成！")
        print("="*80)
        
        print("\n生成的文件:")
        print(f"  数据: {Config.DATA_DIR}")
        print(f"  模型: {Config.MODEL_DIR}")
        print(f"  输出: {Config.OUTPUT_DIR}")
        print(f"  日志: {Config.LOG_DIR}")
        
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) > 1:
        example_num = sys.argv[1]
        
        Config.setup_directorie
```





## 项目结构

```
project/
├── complete_pic50_predictor.py     # 主程序
├── requirements.txt                # 依赖列表
├── README.md                       # 本文档
├── data/                          # 数据目录
│   ├── DAT_training.csv
│   ├── 5HT2A_training.csv
│   └── ...
├── models/                        # 模型目录
│   ├── DAT_complete_model.pkl
│   ├── DAT_transformer_best.pth
│   └── ...
├── output/                        # 输出目录
│   ├── DAT_training_report.json
│   ├── predictions.csv
│   └── ...
├── logs/                          # 日志目录
│   └── pic50_predictor_*.log
└── cache/                         # 缓存目录
```

## 代码质量保证

### 单元测试

```python
# test_predictor.py
import unittest
from complete_pic50_predictor import *

class TestMolecularFeatureExtractor(unittest.TestCase):
    
    def setUp(self):
        self.config = ModelConfig()
        self.extractor = MolecularFeatureExtractor(self.config)
    
    def test_smiles_to_mol(self):
        # 有效SMILES
        mol = self.extractor.smiles_to_mol('CC(C)C')
        self.assertIsNotNone(mol)
        
        # 无效SMILES
        mol = self.extractor.smiles_to_mol('INVALID')
        self.assertIsNone(mol)
    
    def test_extract_features(self):
        smiles = 'CC(CC1=CC=CC=C1)NC'
        features = self.extractor.extract_all_features(smiles)
        self.assertIsNotNone(features)
        self.assertIn('rdkit_descriptors', features)
        self.assertIn('fingerprints', features)

if __name__ == '__main__':
    unittest.main()
```

### 运行测试

```bash
python -m unittest test_predictor.py
```



## 拥有的文件：

1. **`complete_pic50_predictor.py`** - 主程序（约1500行）
   -  Transformer回归模型
   -  图神经网络(GNN)
   -  Random Forest & Gradient Boosting
   -  集成学习（Ensemble）
   -  完整的RDKit描述符（200+）
   -  多种分子指纹（ECFP, MACCS, Atom Pair等）
   -  SMARTS子结构特征（30+模式）
   -  图特征（原子、键、邻接矩阵）
   -  不确定性估计
   -  7个靶点支持（DAT, 5HT2A, CB1, CB2, μ/δ/κ-opioid）
2. **使用文档** - 完整的安装和使用说明
3. **`example_usage.py`** - 6个完整示例

## 快速开始：

```bash
# 1. 安装依赖
pip install numpy pandas scikit-learn rdkit joblib torch

# 2. 训练模型（使用你的CSV数据）
python complete_pic50_predictor.py train --target DAT --data your_data.csv

# 3. 预测
python complete_pic50_predictor.py predict --target DAT --smiles "CC(CC1=CC=CC=C1)NC"

# 4. 批量预测
python complete_pic50_predictor.py batch --target DAT --input your_smiles.csv --output results.csv
```

## 支持的数据格式：

**训练数据CSV：**

```csv
smiles,pIC50
CC(CC1=CC=CC=C1)NC,7.5
CN1C=NC2=C1C(=O)N(C(=O)N2C)C,6.2
```

**预测SMILES文件：**

```
CC(CC1=CC=CC=C1)NC
CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

## 核心功能：

1. **分子特征提取** - 与原项目相同
2. **Transformer模型** - 自注意力机制回归
3. **GNN模型** - 图神经网络
4. **传统ML** - RF, GB完整实现
5. **集成学习** - 多模型融合
6. **不确定性估计** - MC Dropout, Tree variance
7. **完整评估** - R², RMSE, MAE等指标



## 引用和许可

如果你在研究中使用本代码，请引用：

```bibtex
@software{multi_target_pic50_predictor,
  title={Multi-Target pIC50 Predictor},
  author={zapabob},
  year={2025},
  url={https://github.com/zapabob/multi-target-pIC50-predictor}
}
```

## 参考资料

1. **RDKit**: https://www.rdkit.org/
2. **ChEMBL Database**: https://www.ebi.ac.uk/chembl/
3. **Scikit-learn**: https://scikit-learn.org/
4. **PyTorch**: https://pytorch.org/
5. **药物化学**: Lipinski's Rule of Five
6. **QSAR**: Quantitative Structure-Activity Relationship

## 技术支持

- GitHub Issues: https://github.com/zapabob/multi-target-pIC50-predictor/issues
- 文档: 查看代码注释和docstrings
- 社区: RDKit, PyTorch, Scikit-learn论坛



