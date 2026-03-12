# Advanced Dual-Target Affinity Prediction System

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🚀 Overview

A **state-of-the-art** deep learning system for predicting drug-protein binding affinity across multiple targets. This system integrates cutting-edge techniques from graph neural networks, attention mechanisms, and contrastive learning to achieve superior prediction performance.

### Key Features

- ✅ **Graph Neural Networks (GNN)**: Process molecular graph structure for enhanced feature extraction
- ✅ **Multi-Head Attention**: Model complex drug-protein interactions
- ✅ **Contrastive Learning**: Improve feature representations through self-supervised learning
- ✅ **Transfer Learning**: Support for pre-trained encoders
- ✅ **Multi-Task Learning**: Joint prediction for multiple targets
- ✅ **Comprehensive Error Handling**: Robust fallback mechanisms at every step
- ✅ **SCI-Quality Visualization**: Publication-ready figures (300 DPI)
- ✅ **Cross-Validation**: Rigorous model validation
- ✅ **Detailed Logging**: Complete execution tracking

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Architecture](#architecture)
- [Configuration](#configuration)
- [Data Format](#data-format)
- [Training](#training)
- [Evaluation](#evaluation)
- [Visualization](#visualization)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## 🔧 Installation

### Option 1: Using pip

```bash
# Clone the repository
git clone https://github.com/your-repo/advanced-dual-target-mtl.git
cd advanced-dual-target-mtl

# Install dependencies
pip install -r requirements.txt
```

### Option 2: Using conda (Recommended)

```bash
# Create conda environment
conda create -n affinity python=3.10
conda activate affinity

# Install PyTorch (adjust for your CUDA version)
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# Install other dependencies
pip install -r requirements.txt
```

### System Requirements

- **Python**: 3.8 or higher
- **PyTorch**: 2.0 or higher
- **Memory**: At least 8GB RAM (16GB recommended)
- **GPU**: Optional but recommended (CUDA 11.0+)

## 🎯 Quick Start

### Basic Execution

```bash
python advanced_main.py
```

This will:
1. Generate synthetic demonstration data (2000 samples)
2. Train the model for 200 epochs with early stopping
3. Perform 5-fold cross-validation
4. Generate 7 publication-quality figures
5. Save results and reports

### Expected Output

```
================================================================================
Advanced Dual-Target Affinity Prediction System
================================================================================

Features:
  ✓ Graph Neural Networks (GNN)
  ✓ Multi-Head Attention Mechanism
  ✓ Contrastive Learning
  ✓ Transfer Learning Support
  ✓ Advanced Error Handling
  ✓ SCI-Quality Visualization
================================================================================

[Training progress...]

================================================================================
EXECUTION COMPLETED SUCCESSFULLY!
================================================================================

Generated Files:
  • Model: models/best_model.pth
  • Results: outputs/results.json
  • Report: outputs/report.md
  • Figures: figures/

Key Results:
  Target 1 R²: 0.8532
  Target 2 R²: 0.8421
  CV Target 1 RMSE: 0.3245 ± 0.0152
  CV Target 2 RMSE: 0.3312 ± 0.0167
================================================================================
```

## 🏗️ Architecture

### Model Components

```
┌─────────────────────────────────────────────────────────────┐
│                   Input Layer                                │
│  Drug Features (128) | Protein 1 (128) | Protein 2 (128)    │
└───────────────┬─────────────────────────────────────────────┘
                │
                ▼
┌───────────────────────────────────────────────────────────┐
│            Graph Neural Network (GNN)                      │
│  • 3 Graph Convolutional Layers                           │
│  • Message Passing Neural Network                         │
│  • Aggregation: Sum/Mean/Max                             │
└───────────────┬───────────────────────────────────────────┘
                │
                ▼
┌───────────────────────────────────────────────────────────┐
│              Feature Encoders                              │
│  Drug Encoder (→ 512)  |  Protein Encoders (→ 512)       │
│  • Batch Normalization                                    │
│  • Dropout (0.3)                                         │
└───────────────┬───────────────────────────────────────────┘
                │
                ▼
┌───────────────────────────────────────────────────────────┐
│         Multi-Head Attention (8 heads)                    │
│  • Drug-Protein Cross-Attention                          │
│  • Self-Attention Mechanism                              │
│  • Residual Connections                                  │
└───────────────┬───────────────────────────────────────────┘
                │
                ▼
┌───────────────────────────────────────────────────────────┐
│            Shared Network                                  │
│  [1024 → 512 → 256]                                       │
│  • Residual Connections                                   │
│  • Batch Normalization                                    │
└───────────────┬───────────────────────────────────────────┘
                │
                ├─────────────────┬─────────────────┐
                ▼                 ▼                 ▼
┌──────────────────┐  ┌──────────────────┐  ┌──────────────────┐
│  Task 1 Head     │  │  Task 2 Head     │  │ Contrastive Head │
│  [256→128→1]     │  │  [256→128→1]     │  │  [256→128]       │
└──────────────────┘  └──────────────────┘  └──────────────────┘
        │                     │                     │
        ▼                     ▼                     ▼
   Affinity 1           Affinity 2         Contrastive Loss
```

### Loss Function

The total loss is a weighted combination:

```
L_total = w₁·L_task1 + w₂·L_task2 + λ·L_contrastive

where:
- L_task1, L_task2: MSE losses for affinity prediction
- L_contrastive: NT-Xent contrastive loss
- w₁, w₂: Task weights (dynamic or fixed)
- λ: Contrastive weight (default: 0.1)
```

## ⚙️ Configuration

### Basic Configuration

Edit `AdvancedConfig` class in `advanced_dual_target_mtl.py`:

```python
class AdvancedConfig:
    # Model Architecture
    drug_embedding_dim = 512
    protein_embedding_dim = 512
    use_gnn = True
    use_attention = True
    use_contrastive = True
    
    # Training
    batch_size = 64
    epochs = 200
    learning_rate = 0.001
    
    # Optimizer
    optimizer_name = 'adamw'  # 'adam', 'adamw', 'sgd'
    scheduler_name = 'cosine_warm'  # 'cosine', 'plateau', 'cosine_warm'
```

### Advanced Options

```python
# GNN Configuration
config.gnn_hidden_dim = 256
config.gnn_num_layers = 3
config.gnn_dropout = 0.2

# Attention Configuration
config.num_attention_heads = 8
config.attention_dropout = 0.1

# Contrastive Learning
config.contrastive_temperature = 0.07
config.contrastive_weight = 0.1

# Early Stopping
config.patience = 20
config.min_delta = 1e-4
```

## 📊 Data Format

### Input Data Structure

```python
# Drug features: [N, drug_feature_dim]
drug_features = np.array([[...], [...], ...])  # Shape: (N, 128)

# Protein features: [N, protein_feature_dim]
protein1_features = np.array([[...], [...], ...])  # Shape: (N, 128)
protein2_features = np.array([[...], [...], ...])  # Shape: (N, 128)

# Labels (binding affinity in kcal/mol)
labels_task1 = np.array([...])  # Shape: (N,)
labels_task2 = np.array([...])  # Shape: (N,)

# Optional: Graph structure for GNN
edge_indices = [torch.tensor([[...], [...]]), ...]  # List of [2, E] tensors
```

### Using Real Data

Replace the synthetic data generation in `advanced_main.py`:

```python
import pandas as pd

# Load your data
df = pd.read_csv('your_data.csv')

# Extract features (example using RDKit for molecules)
from rdkit import Chem
from rdkit.Chem import AllChem

def extract_drug_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=128)
    return np.array(fp)

drug_features = np.array([
    extract_drug_features(smiles) for smiles in df['smiles']
])

# Similarly for protein features...
```

## 🎓 Training

### Standard Training

```python
from advanced_dual_target_mtl import *

# Create model
config = AdvancedConfig()
model = AdvancedDualTargetPredictor(config)

# Create trainer
trainer = AdvancedTrainer(model, config)

# Train
history = trainer.train(train_loader, val_loader)
```

### Training with Custom Parameters

```python
# Modify configuration
config.epochs = 300
config.learning_rate = 0.0005
config.batch_size = 128

# Enable more aggressive training
config.use_contrastive = True
config.contrastive_weight = 0.2

# Train
trainer = AdvancedTrainer(model, config)
history = trainer.train(train_loader, val_loader)
```

### Resume Training

```python
# Load checkpoint
trainer.load_checkpoint('best_model.pth')

# Continue training
history = trainer.train(train_loader, val_loader)
```

## 📈 Evaluation

### Model Evaluation

```python
# Load best model
checkpoint = torch.load('models/best_model.pth')
model.load_state_dict(checkpoint['model_state_dict'])

# Evaluate
metrics = evaluate_model(model, test_loader, config.device)

print(f"Task 1 R²: {metrics['task1_r2']:.4f}")
print(f"Task 2 R²: {metrics['task2_r2']:.4f}")
```

### Cross-Validation

```python
cv_results = cross_validation(
    drug_features, protein1_features, protein2_features,
    labels_task1, labels_task2, config
)

print(f"CV RMSE: {np.mean(cv_results['task1_rmse']):.4f} ± "
      f"{np.std(cv_results['task1_rmse']):.4f}")
```

## 🎨 Visualization

### Generated Figures

The system automatically generates 7 publication-quality figures:

1. **01_training_curves.png**
   - Training and validation loss over epochs
   - Separate curves for each task
   - Learning rate schedule

2. **02_prediction_scatter.png**
   - Predicted vs. true affinity scatter plots
   - Regression lines and perfect prediction lines
   - Performance metrics overlay

3. **03_residual_analysis.png**
   - Residual distribution histograms
   - Residuals vs. predicted values
   - Trend analysis

4. **04_metrics_comparison.png**
   - Bar chart comparing multiple metrics
   - Side-by-side comparison of both tasks

5. **05_error_distribution.png**
   - Box plots of absolute and relative errors
   - Statistical summaries

6. **06_qq_plots.png**
   - Q-Q plots for normality checking
   - Shapiro-Wilk test statistics

7. **07_comprehensive_summary.png** ⭐
   - **Most important figure for publications**
   - 8 subplots with complete overview
   - Ready for direct use in papers

### Custom Visualization

```python
from advanced_visualization import *

# Generate specific plot
plot_training_curves(history, figure_dir='my_figures')
plot_prediction_scatter(metrics, figure_dir='my_figures')

# Generate all plots
generate_all_plots(history, metrics, figure_dir='my_figures')
```

## 🔬 Advanced Usage

### 1. Transfer Learning

```python
# Load pre-trained encoder
config.use_pretrained = True
config.pretrained_path = 'pretrained_encoder.pth'
config.freeze_encoder = True  # Freeze encoder weights

model = AdvancedDualTargetPredictor(config)
```

### 2. Custom Loss Function

```python
class CustomTrainer(AdvancedTrainer):
    def train_epoch(self, train_loader):
        # Add custom loss terms
        for batch in train_loader:
            # ... standard forward pass ...
            
            # Add custom regularization
            custom_loss = torch.mean(model.some_weights ** 2)
            loss = loss + 0.01 * custom_loss
            
            # ... backward pass ...
```

### 3. Model Ensemble

```python
# Train multiple models
models = []
for i in range(5):
    model = AdvancedDualTargetPredictor(config)
    trainer = AdvancedTrainer(model, config)
    history = trainer.train(train_loader, val_loader)
    models.append(model)

# Ensemble prediction
def ensemble_predict(models, data_loader):
    predictions = []
    for model in models:
        model.eval()
        with torch.no_grad():
            for batch in data_loader:
                pred1, pred2, _ = model(batch['drug'], 
                                       batch['protein1'], 
                                       batch['protein2'])
                predictions.append((pred1, pred2))
    
    # Average predictions
    avg_pred1 = torch.mean(torch.stack([p[0] for p in predictions]), dim=0)
    avg_pred2 = torch.mean(torch.stack([p[1] for p in predictions]), dim=0)
    
    return avg_pred1, avg_pred2
```

### 4. Hyperparameter Tuning

```python
import optuna

def objective(trial):
    # Suggest hyperparameters
    config.learning_rate = trial.suggest_loguniform('lr', 1e-5, 1e-2)
    config.batch_size = trial.suggest_int('batch_size', 32, 128)
    config.dropout_rate = trial.suggest_uniform('dropout', 0.1, 0.5)
    
    # Train and evaluate
    model = AdvancedDualTargetPredictor(config)
    trainer = AdvancedTrainer(model, config)
    history = trainer.train(train_loader, val_loader)
    
    return trainer.best_val_loss

# Run optimization
study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=50)

print("Best parameters:", study.best_params)
```

## 🐛 Troubleshooting

### Common Issues

#### 1. Out of Memory (OOM)

**Problem**: CUDA out of memory error

**Solution**:
```python
# Reduce batch size
config.batch_size = 32  # or smaller

# Reduce model size
config.hidden_dims = [512, 256]
config.drug_embedding_dim = 256
config.protein_embedding_dim = 256
```

#### 2. Training Not Converging

**Problem**: Loss not decreasing or unstable

**Solution**:
```python
# Lower learning rate
config.learning_rate = 0.0001

# Increase dropout
config.dropout_rate = 0.5

# Add gradient clipping (already included)
torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
```

#### 3. Overfitting

**Problem**: Train loss decreasing but validation loss increasing

**Solution**:
```python
# Increase regularization
config.dropout_rate = 0.5
config.weight_decay = 1e-4

# Use data augmentation
config.augment_data = True
config.noise_level = 0.05

# Early stopping (already enabled)
config.patience = 15
```

#### 4. Underfitting

**Problem**: Both train and validation loss remain high

**Solution**:
```python
# Increase model capacity
config.hidden_dims = [2048, 1024, 512, 256]
config.drug_embedding_dim = 1024
config.protein_embedding_dim = 1024

# Train longer
config.epochs = 500

# Reduce regularization
config.dropout_rate = 0.2
config.weight_decay = 1e-6
```

### Error Handling

The system includes comprehensive error handling:

- **Data Validation**: Automatic NaN detection and filling
- **Model Failures**: Fallback to zero tensors with warnings
- **Training Interrupts**: Automatic checkpoint saving
- **Visualization Errors**: Graceful degradation with logging

All errors are logged to `training.log` with full stack traces.

## 📚 Citation

If you use this system in your research, please cite:

```bibtex
@software{advanced_dual_target_mtl,
  title = {Advanced Dual-Target Affinity Prediction with Multi-Task Learning},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/your-repo/advanced-dual-target-mtl},
  note = {Graph Neural Networks, Multi-Head Attention, Contrastive Learning}
}
```

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📧 Contact

For questions or issues:
- GitHub Issues: [Create an issue](https://github.com/your-repo/issues)
- Email: your.email@example.com

## 🙏 Acknowledgments

This project builds upon research in:
- Graph Neural Networks (Gilmer et al., 2017)
- Attention Mechanisms (Vaswani et al., 2017)
- Contrastive Learning (Chen et al., 2020)
- Multi-Task Learning (Caruana, 1997)

## 📝 Changelog

### Version 2.0 (2025-10-23)

- ✅ Added Graph Neural Network support
- ✅ Implemented Multi-Head Attention mechanism
- ✅ Integrated Contrastive Learning
- ✅ Added Transfer Learning support
- ✅ Enhanced error handling and robustness
- ✅ Improved SCI-quality visualization
- ✅ Added comprehensive documentation
- ✅ Optimized for publication-ready results

### Version 1.0 (2025-01-01)

- Initial release with basic MTL model

---

**Happy Researching! 🚀**

*For more information, visit our [documentation](https://your-docs-site.com)*
