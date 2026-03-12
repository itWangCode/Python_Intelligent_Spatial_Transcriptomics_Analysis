# 📚 Advanced Dual-Target Affinity Prediction - File Index

## 🗂️ 文件导航

### 快速开始
- **START_HERE.md** (本文件) - 开始阅读这里！
- **quick_start.py** - 快速示例脚本
- **PROJECT_DELIVERY.md** - 项目交付说明（简版）

### 核心代码（必读）
1. **advanced_dual_target_mtl.py** (60KB, 2000+行)
   - 主要模型实现
   - 包含：GNN, 注意力, 对比学习, 训练器
   
2. **advanced_visualization.py** (30KB, 800+行)
   - SCI级别可视化
   - 7个高质量图表生成函数
   
3. **advanced_main.py** (21KB, 500+行)
   - 主执行脚本
   - 完整的训练、验证、测试流程

### 文档（按阅读顺序）
1. **README.md** (19KB) ⭐
   - **最重要的文档**
   - 完整的项目介绍
   - 安装、使用、示例

2. **INSTALLATION_GUIDE.md** (12KB)
   - 详细安装步骤
   - 常见问题解决
   - 配置指南

3. **TECHNICAL_DETAILS.md** (14KB)
   - 技术实现细节
   - 科学原理
   - 架构设计

### 配置文件
- **requirements.txt** - Python依赖包列表

---

## 🚀 3分钟快速开始

### 方法1: 检查依赖并运行演示

```bash
python quick_start.py
```

这将：
1. ✅ 检查所有依赖是否安装
2. ✅ 运行一个快速演示（10轮训练）
3. ✅ 显示基本结果

### 方法2: 完整训练

```bash
# 安装依赖（如果还没安装）
pip install torch numpy pandas scikit-learn scipy matplotlib seaborn tqdm

# 运行完整训练
python advanced_main.py
```

这将：
1. ✅ 生成2000个合成样本
2. ✅ 训练200轮（带早停）
3. ✅ 执行5折交叉验证
4. ✅ 生成7张SCI级别图表
5. ✅ 保存结果到outputs/

---

## 📖 详细阅读路径

### 如果你是新手
1. 先读 **README.md** 的"Quick Start"部分
2. 运行 `python quick_start.py` 看看效果
3. 阅读 **INSTALLATION_GUIDE.md** 了解配置
4. 阅读 **README.md** 的其他部分

### 如果你有经验
1. 直接运行 `python advanced_main.py`
2. 阅读 **TECHNICAL_DETAILS.md** 了解实现
3. 查看代码：**advanced_dual_target_mtl.py**
4. 根据需要自定义配置

### 如果你要发论文
1. 阅读 **README.md** 的"Citation"部分
2. 阅读 **TECHNICAL_DETAILS.md** 的"Publications"部分
3. 运行完整训练获取结果
4. 使用 `figures/07_comprehensive_summary.png` 作为主图

---

## 🎯 核心功能位置

### 图神经网络 (GNN)
- **文件**: advanced_dual_target_mtl.py
- **行数**: 104-281
- **类**: `GraphConvLayer`, `MolecularGNN`

### 多头注意力
- **文件**: advanced_dual_target_mtl.py
- **行数**: 286-441
- **类**: `MultiHeadAttention`, `DrugProteinAttention`

### 对比学习
- **文件**: advanced_dual_target_mtl.py
- **行数**: 446-494
- **类**: `ContrastiveLoss`

### 主模型
- **文件**: advanced_dual_target_mtl.py
- **行数**: 499-710
- **类**: `AdvancedDualTargetPredictor`

### 可视化
- **文件**: advanced_visualization.py
- **函数**: 
  - `plot_training_curves()` - 训练曲线
  - `plot_prediction_scatter()` - 预测散点图
  - `plot_residual_analysis()` - 残差分析
  - `plot_comprehensive_summary()` - 综合总结 ⭐

---

## 📁 运行后的目录结构

```
.
├── advanced_dual_target_mtl.py   # 核心模型
├── advanced_visualization.py     # 可视化
├── advanced_main.py              # 主程序
├── quick_start.py                # 快速开始
├── requirements.txt              # 依赖
├── README.md                     # 主文档
├── INSTALLATION_GUIDE.md         # 安装指南
├── TECHNICAL_DETAILS.md          # 技术细节
├── PROJECT_DELIVERY.md           # 交付说明
├── START_HERE.md                 # 本文件
│
├── training.log                  # 训练日志
│
├── models/                       # 模型保存
│   └── best_model.pth
│
├── outputs/                      # 结果输出
│   ├── config.json              # 配置
│   ├── results.json             # 详细结果
│   └── report.md                # Markdown报告
│
└── figures/                      # 图表
    ├── 01_training_curves.png
    ├── 02_prediction_scatter.png
    ├── 03_residual_analysis.png
    ├── 04_metrics_comparison.png
    ├── 05_error_distribution.png
    ├── 06_qq_plots.png
    └── 07_comprehensive_summary.png  ⭐
```

---

## 💡 常见问题快速解答

### Q1: 如何使用自己的数据？
**A**: 修改 `advanced_main.py` 第180行左右，替换 `generate_synthetic_data()` 为你的数据加载代码。参考 INSTALLATION_GUIDE.md 第2节。

### Q2: 如何调整模型参数？
**A**: 修改 `advanced_main.py` 中的 `config` 对象，或直接编辑 `advanced_dual_target_mtl.py` 中的 `AdvancedConfig` 类。

### Q3: 内存不足怎么办？
**A**: 减小 `config.batch_size` 和 `config.hidden_dims`。参考 INSTALLATION_GUIDE.md 的"常见问题"部分。

### Q4: 训练太慢怎么办？
**A**: 使用GPU（`config.device = torch.device('cuda')`），或减少 `config.epochs` 和 `config.cv_folds`。

### Q5: 如何获得更好的性能？
**A**: 
1. 增加模型容量（`config.hidden_dims = [2048, 1024, 512, 256]`）
2. 训练更长时间（`config.epochs = 500`）
3. 启用所有功能（`use_gnn=True`, `use_attention=True`, `use_contrastive=True`）

---

## 🔍 代码特性速查

### 错误处理
- ✅ 所有关键函数都有try-except
- ✅ 自动NaN处理
- ✅ 回退机制
- ✅ 详细日志

### 训练特性
- ✅ 早停 (Early Stopping)
- ✅ 学习率调度 (3种)
- ✅ 梯度裁剪
- ✅ 动态任务权重
- ✅ 检查点保存

### 评估特性
- ✅ K折交叉验证
- ✅ 多指标评估 (10+)
- ✅ 残差分析
- ✅ 统计检验

### 可视化特性
- ✅ 300 DPI分辨率
- ✅ 色盲友好配色
- ✅ 完整统计信息
- ✅ 7张SCI级别图表

---

## 📞 获取帮助

### 在线资源
- **GitHub**: 提交Issue
- **文档**: 阅读完整的README.md
- **日志**: 检查training.log

### 离线资源
- **示例**: 运行quick_start.py
- **指南**: 阅读INSTALLATION_GUIDE.md
- **技术**: 阅读TECHNICAL_DETAILS.md

---

## ✅ 快速检查清单

### 运行前
- [ ] Python 3.8+ 已安装
- [ ] PyTorch 已安装
- [ ] 其他依赖已安装（运行quick_start.py检查）
- [ ] 有足够的磁盘空间（至少1GB）

### 运行中
- [ ] 训练损失在下降
- [ ] 验证损失稳定
- [ ] 没有NaN或Inf错误
- [ ] 日志正常输出

### 运行后
- [ ] models/best_model.pth 已生成
- [ ] outputs/ 目录有结果
- [ ] figures/ 目录有图表
- [ ] training.log 无严重错误

---

## 🎉 准备好了吗？

### 立即开始！

```bash
# 方案A: 快速演示（5分钟）
python quick_start.py

# 方案B: 完整训练（20分钟 CPU / 5分钟 GPU）
python advanced_main.py
```

---

## 📚 推荐阅读顺序

1. **START_HERE.md** (本文件) ← 你在这里
2. **README.md** - 完整文档
3. **quick_start.py** - 运行示例
4. **INSTALLATION_GUIDE.md** - 详细指南
5. **TECHNICAL_DETAILS.md** - 深入了解

---

**开始你的科研之旅！🚀**

有问题？查看 README.md 或 INSTALLATION_GUIDE.md 的"常见问题"部分。
