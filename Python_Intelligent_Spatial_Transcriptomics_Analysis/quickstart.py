#!/usr/bin/env python3
"""
DeepST 快速开始演示
展示完整的分析流程

使用方法:
    python quickstart.py
"""

import os
import sys
import numpy as np
import warnings

warnings.filterwarnings('ignore')

print("""
╔════════════════════════════════════════════════════════════╗
║                                                            ║
║         DeepST 空间转录组学分析 - 快速开始演示              ║
║                                                            ║
║         macOS Intel 兼容版本                               ║
║         适用于博士研究和 AIDD 应用                          ║
║                                                            ║
╚════════════════════════════════════════════════════════════╝
""")

print("\n检查依赖...")

# 检查必要的库
required_modules = {
    'torch': 'PyTorch',
    'torch_geometric': 'PyTorch Geometric',
    'scanpy': 'Scanpy',
    'anndata': 'AnnData',
    'sklearn': 'Scikit-learn',
    'matplotlib': 'Matplotlib'
}

missing_modules = []

for module, name in required_modules.items():
    try:
        __import__(module)
        print(f"  ✓ {name}")
    except ImportError:
        print(f"  ✗ {name} - 未安装")
        missing_modules.append(name)

if missing_modules:
    print(f"\n❌ 缺少依赖: {', '.join(missing_modules)}")
    print("\n请先运行环境配置脚本:")
    print("  ./setup_deepst_macos.sh")
    print("\n或手动安装依赖:")
    print("  pip install -r requirements.txt")
    sys.exit(1)

print("\n✓ 所有依赖已安装\n")

# 导入模块
from deepst_core import DeviceManager
from deepst_analysis import DeepSTAnalyzer

print("="*60)
print("演示 1: 设备检测")
print("="*60)

dm = DeviceManager(prefer_gpu=True)
print(f"\n当前使用设备: {dm.device}")
print(f"CPU线程数: {dm.device.type == 'cpu' and 'N/A' or 'GPU模式'}")

print("\n" + "="*60)
print("演示 2: 生成模拟数据")
print("="*60)

from sklearn.datasets import make_blobs
import anndata as ad
import scanpy as sc

# 生成模拟空间转录组学数据
print("\n生成模拟数据...")
n_spots = 1000
n_genes = 2000
n_domains = 7

# 空间坐标 (模拟组织切片)
spatial_coords, true_labels = make_blobs(
    n_samples=n_spots,
    n_features=2,
    centers=n_domains,
    cluster_std=5.0,
    random_state=42
)

# 基因表达 (负二项分布,模拟RNA-seq数据)
gene_expression = np.random.negative_binomial(5, 0.3, (n_spots, n_genes))

# 创建 AnnData 对象
adata = ad.AnnData(X=gene_expression)
adata.obsm['spatial'] = spatial_coords
adata.obs['ground_truth'] = true_labels.astype(str)
adata.var_names = [f'Gene_{i}' for i in range(n_genes)]

print(f"✓ 数据生成完成:")
print(f"  - Spots: {n_spots}")
print(f"  - Genes: {n_genes}")
print(f"  - 真实域数: {n_domains}")

print("\n" + "="*60)
print("演示 3: DeepST 分析流程")
print("="*60)

# 创建分析器
analyzer = DeepSTAnalyzer(
    n_domains=n_domains,
    hidden_dim=256,      # 较小的模型用于演示
    latent_dim=64,
    use_gpu=True,
    output_dir='./quickstart_results'
)

# 设置数据
analyzer.adata = adata

print("\n[1/6] 数据预处理...")
analyzer.prepare_data(
    use_highly_variable=False,
    pca_components=100,
    spatial_smooth=True
)

print("\n[2/6] 构建模型...")
analyzer.build_model()

print("\n[3/6] 训练模型 (100 epochs - 快速演示)...")
history = analyzer.train(n_epochs=100, print_every=25)

print("\n[4/6] 提取嵌入表示...")
embeddings = analyzer.get_embeddings()

print("\n[5/6] 识别空间域...")
domain_labels = analyzer.identify_domains()

print("\n[6/6] 生成可视化...")
analyzer.visualize_results()

print("\n" + "="*60)
print("演示 4: 评估结果")
print("="*60)

metrics = analyzer.evaluate_performance('ground_truth')

print("\n性能指标:")
print(f"  - ARI (调整兰德指数): {metrics['ARI']:.4f}")
print(f"  - NMI (归一化互信息): {metrics['NMI']:.4f}")

print("\n" + "="*60)
print("演示 5: AIDD 应用")
print("="*60)

from deepst_aidd import (
    SpatialTargetIdentifier,
    TumorMicroenvironmentAnalyzer
)

print("\n识别空间域标志基因...")
target_id = SpatialTargetIdentifier(analyzer.adata)
markers = target_id.identify_domain_markers(top_n=10)

print("\n前10个域0的标志基因:")
domain_0_markers = markers[markers['domain'] == '0'].head(10)
print(domain_0_markers[['names', 'logfoldchanges', 'pvals_adj']].to_string(index=False))

print("\n注释肿瘤微环境...")
tme = TumorMicroenvironmentAnalyzer(analyzer.adata)

# 随机选择一些基因作为标志物用于演示
all_genes = analyzer.adata.var_names.tolist()
tumor_markers = all_genes[:3]
immune_markers = all_genes[3:6]
stromal_markers = all_genes[6:9]

annotation = tme.annotate_tumor_regions(
    tumor_markers,
    immune_markers,
    stromal_markers
)

print("\n" + "="*60)
print("演示完成!")
print("="*60)

print("\n结果保存在: ./quickstart_results/")
print("\n生成的文件:")
print("  ✓ deepst_results.h5ad - 完整分析结果")
print("  ✓ domain_labels.csv - 空间域标签")
print("  ✓ metrics.csv - 评估指标")
print("  ✓ figures/ - 可视化图表")

print("\n下一步:")
print("  1. 查看 README.md 了解详细使用方法")
print("  2. 查看 DeepST_Analysis.md 了解技术原理")
print("  3. 使用您自己的数据运行分析")
print("  4. 探索 AIDD 扩展功能")

print("\n" + "="*60)
print("Happy analyzing! 🔬")
print("="*60 + "\n")
