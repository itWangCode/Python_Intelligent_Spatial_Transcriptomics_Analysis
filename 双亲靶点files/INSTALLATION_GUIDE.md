# Advanced Dual-Target MTL - Installation and Usage Guide

## 📦 Installation Instructions

### Step 1: Install Python and Dependencies

#### For Linux/Mac:

```bash
# Install PyTorch (CPU version)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# Or for CUDA 11.8 (GPU version)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Install other dependencies
pip install numpy pandas scikit-learn scipy matplotlib seaborn tqdm
```

#### For Windows:

```bash
# Install PyTorch (CPU version)
pip install torch torchvision torchaudio

# Install other dependencies
pip install numpy pandas scikit-learn scipy matplotlib seaborn tqdm
```

#### Using Conda (Recommended):

```bash
# Create new environment
conda create -n affinity python=3.10
conda activate affinity

# Install PyTorch with CUDA
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# Install other packages
conda install numpy pandas scikit-learn scipy matplotlib seaborn tqdm -c conda-forge
```

### Step 2: Verify Installation

```bash
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
```

## 🚀 Quick Start Guide

### 1. Basic Execution

```bash
# Run with default settings
python advanced_main.py
```

Expected output:
```
================================================================================
Advanced Dual-Target Affinity Prediction System
================================================================================

Step 1: Loading configuration...
Step 2: Generating/Loading data...
Step 3: Preparing data loaders...
Step 4: Creating and training model...
[Training progress...]
Step 5: Evaluating best model...
Step 6: Performing cross-validation...
Step 7: Generating visualizations...
Step 8: Saving results...

================================================================================
EXECUTION COMPLETED SUCCESSFULLY!
================================================================================
```

### 2. Using Your Own Data

#### Prepare Your Data File (CSV format):

```csv
drug_smiles,protein1_sequence,protein2_sequence,affinity_1,affinity_2
CCO,MKTAYIAKQ...,MKGDLRPEV...,-8.5,-9.2
CC(=O)OC1=CC...,MGAAASIQG...,MTEYKLVVV...,-7.8,-8.1
...
```

#### Modify `advanced_main.py`:

```python
# Replace the data generation section (around line 180)

# OLD CODE:
drug_features, protein1_features, protein2_features, \
    labels_task1, labels_task2 = generate_synthetic_data(...)

# NEW CODE:
import pandas as pd
from feature_extractor import extract_drug_features, extract_protein_features

# Load your data
df = pd.read_csv('your_data.csv')

# Extract features
drug_features = np.array([
    extract_drug_features(smiles) for smiles in df['drug_smiles']
])
protein1_features = np.array([
    extract_protein_features(seq) for seq in df['protein1_sequence']
])
protein2_features = np.array([
    extract_protein_features(seq) for seq in df['protein2_sequence']
])

labels_task1 = df['affinity_1'].values
labels_task2 = df['affinity_2'].values
```

### 3. Customizing Model Configuration

Edit the configuration section in `advanced_main.py`:

```python
# Load default config
config = AdvancedConfig()

# Customize architecture
config.drug_embedding_dim = 1024  # Increase for more capacity
config.protein_embedding_dim = 1024
config.hidden_dims = [2048, 1024, 512, 256]

# Customize GNN
config.use_gnn = True
config.gnn_hidden_dim = 512
config.gnn_num_layers = 4

# Customize attention
config.use_attention = True
config.num_attention_heads = 16

# Customize training
config.batch_size = 128
config.epochs = 300
config.learning_rate = 0.0005

# Save configuration
config.save('my_config.json')
```

## 📊 Understanding the Output

### Directory Structure After Execution

```
.
├── advanced_dual_target_mtl.py    # Core model code
├── advanced_visualization.py       # Visualization module
├── advanced_main.py               # Main execution script
├── requirements.txt               # Dependencies
├── README.md                      # Documentation
├── training.log                   # Detailed logs
│
├── models/                        # Saved models
│   └── best_model.pth            # Best model checkpoint
│
├── outputs/                       # Results
│   ├── config.json               # Configuration used
│   ├── results.json              # Detailed metrics
│   └── report.md                 # Markdown report
│
└── figures/                       # Visualizations
    ├── 01_training_curves.png
    ├── 02_prediction_scatter.png
    ├── 03_residual_analysis.png
    ├── 04_metrics_comparison.png
    ├── 05_error_distribution.png
    ├── 06_qq_plots.png
    └── 07_comprehensive_summary.png  ⭐ Main figure
```

### Understanding Results

#### results.json

```json
{
  "test_metrics": {
    "task1": {
      "rmse": 0.3245,      // Root Mean Square Error
      "mae": 0.2512,       // Mean Absolute Error
      "r2": 0.8532,        // R-squared (higher is better)
      "pcc": 0.9234,       // Pearson Correlation
      "spearman": 0.9187   // Spearman Correlation
    },
    "task2": {...}
  },
  "cross_validation": {
    "task1": {
      "rmse_mean": 0.3312,
      "rmse_std": 0.0167   // Standard deviation across folds
    }
  }
}
```

#### Interpreting Metrics

- **R² (R-squared)**: 
  - > 0.9: Excellent
  - 0.7-0.9: Good
  - 0.5-0.7: Moderate
  - < 0.5: Poor

- **RMSE**: Lower is better (depends on data scale)
  - For affinity (kcal/mol): < 0.5 is good

- **PCC (Pearson)**: Measures linear correlation
  - > 0.9: Strong correlation

## 🎨 Understanding the Figures

### Figure 1: Training Curves (01_training_curves.png)

Shows the learning progress:
- **Panel A**: Overall loss trajectory
- **Panel B**: Task 1 specific loss
- **Panel C**: Task 2 specific loss
- **Panel D**: Learning rate schedule

**What to look for:**
- Decreasing training and validation loss
- Small gap between train and val (no overfitting)
- Smooth convergence

### Figure 2: Prediction Scatter (02_prediction_scatter.png)

Shows prediction accuracy:
- Points should cluster around the diagonal (perfect prediction)
- R² value indicates fit quality
- Outliers indicate difficult cases

### Figure 7: Comprehensive Summary ⭐

**Most important for papers!**
- 8 panels with complete overview
- Ready for direct publication use
- 300 DPI resolution

## 🔧 Advanced Configuration

### 1. Small Dataset (< 500 samples)

```python
config = AdvancedConfig()
config.drug_embedding_dim = 256
config.protein_embedding_dim = 256
config.hidden_dims = [512, 256]
config.dropout_rate = 0.5
config.batch_size = 32
config.epochs = 300
```

### 2. Large Dataset (> 5000 samples)

```python
config = AdvancedConfig()
config.drug_embedding_dim = 1024
config.protein_embedding_dim = 1024
config.hidden_dims = [2048, 1024, 512, 256]
config.dropout_rate = 0.3
config.batch_size = 256
config.epochs = 200
```

### 3. CPU-Only Training

```python
config = AdvancedConfig()
config.device = torch.device('cpu')
config.batch_size = 32  # Smaller for CPU
config.num_workers = 0
```

### 4. GPU Optimization

```python
config = AdvancedConfig()
config.device = torch.device('cuda')
config.batch_size = 256  # Larger for GPU
config.pin_memory = True
config.num_workers = 4
```

## 🐛 Common Issues and Solutions

### Issue 1: "CUDA out of memory"

**Solution 1**: Reduce batch size
```python
config.batch_size = 32  # or 16
```

**Solution 2**: Reduce model size
```python
config.hidden_dims = [512, 256]
config.drug_embedding_dim = 256
```

**Solution 3**: Use CPU
```python
config.device = torch.device('cpu')
```

### Issue 2: Training too slow

**Solution 1**: Use GPU
```python
config.device = torch.device('cuda')
```

**Solution 2**: Reduce cross-validation folds
```python
config.cv_folds = 3  # instead of 5
```

**Solution 3**: Reduce epochs for CV
Modify in `advanced_main.py`:
```python
# In cross_validation function
config.epochs = min(30, original_epochs)  # Faster CV
```

### Issue 3: Poor performance

**Solution 1**: Increase model capacity
```python
config.hidden_dims = [2048, 1024, 512, 256]
config.epochs = 500
```

**Solution 2**: Adjust learning rate
```python
config.learning_rate = 0.0001  # Try different values
```

**Solution 3**: Enable all features
```python
config.use_gnn = True
config.use_attention = True
config.use_contrastive = True
```

### Issue 4: Overfitting

**Solution 1**: Increase regularization
```python
config.dropout_rate = 0.5
config.weight_decay = 1e-4
```

**Solution 2**: Use data augmentation
```python
config.augment_data = True
config.noise_level = 0.05
```

**Solution 3**: Early stopping (already enabled)
```python
config.patience = 15  # Stop earlier
```

## 📈 Performance Benchmarks

### Expected Performance (Synthetic Data)

| Metric | Target 1 | Target 2 |
|--------|----------|----------|
| R² | 0.85-0.90 | 0.84-0.89 |
| RMSE | 0.30-0.35 | 0.31-0.36 |
| Training Time (CPU) | ~20 min | ~20 min |
| Training Time (GPU) | ~5 min | ~5 min |

### System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| RAM | 8 GB | 16 GB |
| GPU | None (CPU) | NVIDIA GPU (4GB+) |
| Storage | 1 GB | 5 GB |
| Python | 3.8 | 3.10+ |

## 🎯 Best Practices

### 1. Data Preparation

✅ **DO:**
- Normalize features
- Check for NaN values
- Balance dataset if possible
- Use stratified splits for CV

❌ **DON'T:**
- Use raw unnormalized data
- Ignore missing values
- Mix different scales

### 2. Model Training

✅ **DO:**
- Start with default config
- Monitor validation loss
- Use early stopping
- Save best model

❌ **DON'T:**
- Overtrain without validation
- Skip cross-validation
- Ignore warning messages

### 3. Result Reporting

✅ **DO:**
- Report mean ± std from CV
- Include all metrics (R², RMSE, PCC)
- Show comprehensive figure
- Document configuration

❌ **DON'T:**
- Report only best fold
- Cherry-pick metrics
- Skip uncertainty estimates

## 📞 Getting Help

### Check the Log File

```bash
tail -n 100 training.log
```

### Enable Debug Mode

```python
import logging
logger.setLevel(logging.DEBUG)
```

### Report Issues

When reporting issues, include:
1. Python version: `python --version`
2. PyTorch version: `python -c "import torch; print(torch.__version__)"`
3. Command used
4. Error message
5. Last 50 lines of `training.log`

## 🎓 Learning Resources

### Understanding the Model

1. **Graph Neural Networks**
   - Paper: "Neural Message Passing for Quantum Chemistry" (Gilmer et al., 2017)
   - Tutorial: https://distill.pub/2021/gnn-intro/

2. **Attention Mechanisms**
   - Paper: "Attention Is All You Need" (Vaswani et al., 2017)
   - Tutorial: https://jalammar.github.io/illustrated-transformer/

3. **Contrastive Learning**
   - Paper: "A Simple Framework for Contrastive Learning" (Chen et al., 2020)
   - Tutorial: https://lilianweng.github.io/posts/2021-05-31-contrastive/

### PyTorch Resources

- Official Tutorial: https://pytorch.org/tutorials/
- Documentation: https://pytorch.org/docs/stable/index.html

## 📝 Next Steps

1. **Experiment with configurations**
   ```bash
   # Try different settings
   python advanced_main.py
   ```

2. **Use your own data**
   - Prepare CSV file
   - Implement feature extraction
   - Train model

3. **Optimize hyperparameters**
   - Use grid search
   - Try Optuna for automatic tuning

4. **Publish results**
   - Use Figure 07 (comprehensive summary)
   - Report CV results with uncertainties
   - Document configuration

---

**Questions?** Check the README.md or create an issue on GitHub.

**Good luck with your research! 🚀**
