# 📁 批量自动化药物发现系统 - 完整指南

## 🎯 功能概述

**自动读取PDB文件夹 → 批量处理 → 统一汇总**

```
PDB文件夹/
├── protein1.pdb  →  完整流程  →  结果1/
├── protein2.pdb  →  完整流程  →  结果2/
├── protein3.pdb  →  完整流程  →  结果3/
└── ...           →  ...       →  ...
                       ↓
                  统一汇总和排名
```

## 🚀 快速开始

### 准备PDB文件夹

```bash
# 创建PDB文件夹
mkdir my_pdbs

# 复制你的PDB文件到文件夹
cp protein1.pdb my_pdbs/
cp protein2.pdb my_pdbs/
cp protein3.pdb my_pdbs/
# ... 更多PDB文件

# 或者使用自动下载
python auto_pdb_downloader.py \
  --mode download \
  --target_type kinases \
  --num_targets 20 \
  --output_dir my_pdbs
```

### 运行批量处理

```bash
# 基础用法 - 自动处理文件夹中所有PDB
python batch_discovery.py --pdb_folder my_pdbs

# 完整参数
python batch_discovery.py \
  --pdb_folder pdb_structures \
  --output batch_results \
  --candidates 100 \
  --sampling_steps 500
```

**就这么简单！系统会自动:**
1. ✅ 扫描文件夹中所有PDB文件
2. ✅ 对每个靶点运行完整工业级流程
3. ✅ 自动保存进度(支持断点续传)
4. ✅ 汇总所有结果并排名
5. ✅ 生成HTML报告

## 📋 完整使用示例

### 示例1: 基础批量处理

```bash
# 处理my_pdbs文件夹中的所有PDB
# 每个靶点生成100个候选分子
python batch_discovery.py \
  --pdb_folder my_pdbs \
  --output batch_results
```

**预期输出:**
```
batch_results/
├── protein1/                    # 靶点1结果
│   ├── molecules/               # 生成的分子
│   ├── complete_results.csv     # 完整结果
│   └── top20_candidates.csv     # Top 20
├── protein2/                    # 靶点2结果
├── protein3/                    # 靶点3结果
├── all_targets_combined.csv     # 所有结果合并 ⭐
├── top50_all_targets.csv        # 全局Top 50 ⭐
├── best_per_target.csv          # 每靶点最佳 ⭐
├── batch_report.html            # HTML报告 ⭐
├── batch_statistics.json        # 统计信息
└── progress.json                # 进度追踪
```

### 示例2: 快速筛选(不含对接)

```bash
# 快速模式 - 跳过对接步骤
python batch_discovery.py \
  --pdb_folder my_pdbs \
  --candidates 50 \
  --sampling_steps 300 \
  --no_docking \
  --output quick_batch
  
  
  
  
python industrial_pipeline.py \
  --pdb pdb_structures/3WWS.pdb \
  --candidates 10 \
  --output test_results
```





**用时**: 约5-10分钟/靶点 (CPU)

### 示例3: 高质量GPU模式

```bash
# GPU加速 + 高质量采样
python batch_discovery.py \
  --pdb_folder my_pdbs \
  --candidates 500 \
  --sampling_steps 1000 \
  --device cuda \
  --output premium_batch
```

**用时**: 约30-60分钟/靶点 (GPU)

### 示例4: 限制处理数量(测试)

```bash
# 只处理前5个PDB文件
python batch_discovery.py \
  --pdb_folder my_pdbs \
  --max_targets 5 \
  --output test_batch
```

### 示例5: 断点续传

```bash
# 第一次运行(中断)
python batch_discovery.py --pdb_folder my_pdbs

# 继续运行(自动跳过已完成的)
python batch_discovery.py --pdb_folder my_pdbs

# 强制重新运行
python batch_discovery.py --pdb_folder my_pdbs --no_resume
```

## 📊 输出结果详解

### 1. all_targets_combined.csv

**所有靶点的所有分子汇总**

| 字段 | 说明 |
|------|------|
| target | 靶点名称 |
| mol_id | 分子ID |
| final_score | 综合评分 |
| admet_score | ADMET评分 |
| druglikeness_score | 类药性评分 |
| vina_score | 对接评分 |
| recommended | 是否推荐 |
| ... | 更多字段 |

**用法:**
```python
import pandas as pd

df = pd.read_csv('batch_results/all_targets_combined.csv')
print(f"总分子数: {len(df)}")
print(f"推荐数: {df['recommended'].sum()}")
```

### 2. top50_all_targets.csv

**全局Top 50最佳候选**

按 `final_score` 排序，推荐进行实验验证

**用法:**
```python
top50 = pd.read_csv('batch_results/top50_all_targets.csv')

# 查看最佳候选
print(top50[['target', 'mol_id', 'final_score', 
            'admet_score', 'vina_score']].head(10))
```

### 3. best_per_target.csv

**每个靶点的最佳分子**

快速查看每个靶点的最优候选

**用法:**
```python
best = pd.read_csv('batch_results/best_per_target.csv')

# 按靶点查看
for _, row in best.iterrows():
    print(f"{row['target']}: Score={row['final_score']:.1f}, "
          f"Vina={row['vina_score']:.2f}")
```

### 4. batch_report.html

**可视化HTML报告**

在浏览器中打开查看:
- 总体统计
- Top 20候选表格
- 各靶点汇总
- 图表和可视化

```bash
# 打开报告
open batch_results/batch_report.html
# 或
xdg-open batch_results/batch_report.html
```

### 5. batch_statistics.json

**详细统计信息**

```json
{
  "total_targets": 10,
  "completed": 9,
  "failed": 1,
  "success_rate": 90.0,
  "total_molecules": 850,
  "total_recommended": 45,
  "results": [...]
}
```

## 🎯 完整参数说明

| 参数 | 说明 | 默认值 | 示例 |
|------|------|--------|------|
| **必需参数** |
| `--pdb_folder` | PDB文件夹路径 | 必需 | `./my_pdbs` |
| **输出参数** |
| `--output` | 输出目录 | `batch_results` | `./results` |
| **生成参数** |
| `--candidates` | 每靶点候选数 | 100 | 500 |
| `--num_atoms` | 原子数范围 | [15, 35] | 20 40 |
| `--sampling_steps` | 采样步数 | 500 | 1000 |
| **处理参数** |
| `--no_docking` | 禁用对接 | False | - |
| `--max_targets` | 最大处理数 | None(全部) | 10 |
| `--no_resume` | 不续传 | False | - |
| **设备参数** |
| `--device` | 计算设备 | 自动 | cuda |

## 📈 性能和时间估算

### CPU模式

| 配置 | 靶点数 | 候选/靶点 | 总时间 | 推荐数 |
|------|--------|-----------|--------|--------|
| 快速 | 5 | 50 | ~30分钟 | 10-15 |
| 标准 | 10 | 100 | ~2小时 | 20-40 |
| 完整 | 20 | 100 | ~4小时 | 40-80 |

### GPU模式 (RTX 3090)

| 配置 | 靶点数 | 候选/靶点 | 总时间 | 推荐数 |
|------|--------|-----------|--------|--------|
| 快速 | 10 | 100 | ~20分钟 | 20-40 |
| 标准 | 50 | 200 | ~2小时 | 100-200 |
| 大规模 | 100 | 500 | ~8小时 | 300-500 |

## 🔥 实战案例

### 案例1: 癌症靶点批量筛选

```bash
# 1. 下载20个癌症相关激酶
python auto_pdb_downloader.py \
  --mode download \
  --target_type kinases \
  --disease cancer \
  --num_targets 20 \
  --output_dir cancer_kinases

# 2. 批量处理
python batch_discovery.py \
  --pdb_folder cancer_kinases \
  --candidates 200 \
  --sampling_steps 800 \
  --device cuda \
  --output cancer_screening

# 3. 查看结果
python -c "
import pandas as pd
df = pd.read_csv('cancer_screening/top50_all_targets.csv')
print(df[['target', 'mol_id', 'final_score', 'vina_score']].head(20))
"
```

**预期产出:**
- 20个靶点
- ~4000个初始分子
- ~200个推荐候选
- Top 50全局最佳

### 案例2: 多靶点组合筛选

```bash
# 准备混合靶点文件夹
mkdir multi_targets
cp kinase_targets/*.pdb multi_targets/
cp gpcr_targets/*.pdb multi_targets/
cp protease_targets/*.pdb multi_targets/

# 批量处理
python batch_discovery.py \
  --pdb_folder multi_targets \
  --candidates 150 \
  --device cuda \
  --output multi_target_screen
```

### 案例3: 快速预筛选

```bash
# 第一轮: 快速筛选(不对接)
python batch_discovery.py \
  --pdb_folder all_targets \
  --candidates 50 \
  --no_docking \
  --output round1

# 选择Top靶点
# (手动或脚本选择最佳靶点的PDB)

# 第二轮: 精细筛选(含对接)
python batch_discovery.py \
  --pdb_folder top_targets \
  --candidates 500 \
  --sampling_steps 1000 \
  --device cuda \
  --output round2_refined
```

## 🛠️ 高级用法

### Python API

```python
from batch_discovery import BatchDrugDiscovery

# 初始化
batch = BatchDrugDiscovery(
    pdb_folder='my_pdbs',
    output_base='results',
    device='cuda'
)

# 自定义流程参数
stats = batch.run_batch_processing(
    candidates_per_target=200,
    num_atoms_range=(20, 35),
    sampling_steps=800,
    enable_docking=True,
    max_targets=None  # 处理全部
)

# 查看统计
print(f"成功率: {stats['success_rate']:.1f}%")
print(f"推荐候选: {stats['total_recommended']}")
```

### 自定义筛选阈值

```python
from batch_discovery import BatchDrugDiscovery

batch = BatchDrugDiscovery(pdb_folder='my_pdbs')

# 修改筛选标准
batch.pipeline.thresholds = {
    'admet_score': 65.0,  # 更严格
    'druglikeness_score': 60.0,
    'lipinski_violations': 0,
    'pains_alerts': 0,
    'vina_score': -7.5,
}

stats = batch.run_batch_processing(candidates_per_target=100)
```

### 结果后处理

```python
import pandas as pd
import matplotlib.pyplot as plt

# 读取所有结果
df = pd.read_csv('batch_results/all_targets_combined.csv')

# 1. 按靶点统计
target_stats = df.groupby('target').agg({
    'final_score': 'mean',
    'recommended': 'sum',
    'mol_id': 'count'
}).rename(columns={'mol_id': 'total'})

print(target_stats.sort_values('recommended', ascending=False))

# 2. 筛选优质候选
premium = df[
    (df['final_score'] >= 70) &
    (df['admet_score'] >= 65) &
    (df['lipinski_pass'] == True) &
    (df['vina_score'] < -7.0) &
    (df['pains_alerts'] == 0)
]

print(f"优质候选数: {len(premium)}")
premium.to_csv('premium_candidates.csv', index=False)

# 3. 可视化
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# 靶点分布
df.groupby('target')['recommended'].sum().plot(kind='bar', ax=axes[0,0])
axes[0,0].set_title('Recommended per Target')

# 评分分布
axes[0,1].hist(df['final_score'], bins=30)
axes[0,1].set_title('Final Score Distribution')

# ADMET vs Vina
axes[1,0].scatter(df['admet_score'], df['vina_score'], alpha=0.5)
axes[1,0].set_xlabel('ADMET Score')
axes[1,0].set_ylabel('Vina Score')

# 类药性分布
axes[1,1].hist(df['druglikeness_score'], bins=30)
axes[1,1].set_title('DrugLikeness Distribution')

plt.tight_layout()
plt.savefig('batch_analysis.png', dpi=300)
```

## 💡 最佳实践

### 1. 渐进式筛选策略

```bash
# 阶段1: 快速探索 (10个靶点)
python batch_discovery.py \
  --pdb_folder targets \
  --max_targets 10 \
  --candidates 50 \
  --no_docking

# 阶段2: 选择性深入 (5个最佳靶点)
python batch_discovery.py \
  --pdb_folder top5_targets \
  --candidates 300 \
  --sampling_steps 1000 \
  --device cuda

# 阶段3: 精细优化 (2个靶点)
python batch_discovery.py \
  --pdb_folder final_targets \
  --candidates 1000 \
  --sampling_steps 1500 \
  --device cuda
```

### 2. 并行处理多GPU

```bash
# GPU 0
CUDA_VISIBLE_DEVICES=0 python batch_discovery.py \
  --pdb_folder batch1 --device cuda &

# GPU 1
CUDA_VISIBLE_DEVICES=1 python batch_discovery.py \
  --pdb_folder batch2 --device cuda &
```

### 3. 定时批量处理

```bash
# 创建批处理脚本
cat > batch_job.sh << 'EOF'
#!/bin/bash
python batch_discovery.py \
  --pdb_folder daily_targets \
  --output results_$(date +%Y%m%d) \
  --candidates 200 \
  --device cuda
EOF

# 设置定时任务
crontab -e
# 每天凌晨2点运行
0 2 * * * /path/to/batch_job.sh
```

## 🆘 故障排除

### Q1: 某个PDB处理失败?

**A**: 系统会自动跳过并继续处理其他PDB
- 检查 `progress.json` 中的 `failed` 列表
- 查看日志: `batch_processing.log`

### Q2: 断点续传不工作?

**A**: 
```bash
# 检查进度文件
cat batch_results/progress.json

# 清除进度重新运行
rm batch_results/progress.json
python batch_discovery.py --pdb_folder my_pdbs
```

### Q3: 内存不足?

**A**: 
1. 减少每批处理的候选数
2. 减少 `max_targets`
3. 使用 `--no_docking`

### Q4: 如何只处理特定PDB?

**A**: 
```bash
# 创建子文件夹
mkdir selected_pdbs
cp target1.pdb selected_pdbs/
cp target2.pdb selected_pdbs/

python batch_discovery.py --pdb_folder selected_pdbs
```

## 📞 技术支持

- 查看完整日志: `batch_processing.log`
- 检查进度: `progress.json`
- 统计信息: `batch_statistics.json`

---

**开始您的批量药物发现之旅! 🚀**
