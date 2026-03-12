# 🏭 工业级药物发现完整流程

## 系统架构

```
输入PDB → 分子生成 → ADMET筛选 → 类药性过滤 → 分子对接 → 综合评分 → 推荐候选
          ↓           ↓             ↓              ↓           ↓            ↓
        扩散模型    化学验证      药代动力学      Lipinski     Vina       Top-N
                                   毒性预测       PAINS       Scoring    Ranking
```

## 🚀 快速开始

### 完整流程 (一键运行)

```bash
# 基础用法 (100个候选 → ADMET → 类药性 → 对接)
python industrial_pipeline.py \
  --pdb pdb_structures/4U3Y.pdb \
  --candidates 100 \
  --output industrial_results

# 高质量模式 (500个候选, 1000步采样)
python industrial_pipeline.py \
  --pdb pdb_structures/4U3Y.pdb  \
  --candidates 500 \
  --sampling_steps 1000 \
  --device cpu \
  --output high_quality_results
```

### 不含对接 (加速)

```bash
python industrial_pipeline.py \
  --pdb protein.pdb \
  --candidates 200 \
  --no_docking \
  --output fast_results
```

## 📊 完整流程详解

### 阶段1: 分子生成

**目标**: 生成候选分子库

**方法**:
- E(3)等变扩散模型
- 蛋白质口袋引导生成
- 可配置原子数和采样步数

**输出**: 
- 有效3D分子结构
- SDF文件
- SMILES字符串

**质量控制**:
- 化学价验证
- 连通性检查
- 基本结构合理性

---

### 阶段2: ADMET预测

**目标**: 评估药代动力学和毒性

#### A - Absorption (吸收)
| 指标 | 标准 | 解释 |
|------|------|------|
| **Caco-2通透性** | > -5.5 | 肠道吸收能力 |
| **HIA** | > 70% | 人肠道吸收率 |
| **P-gp底物** | False | 不是外排泵底物 |

#### D - Distribution (分布)
| 指标 | 标准 | 解释 |
|------|------|------|
| **Vd** | 0.5-5.0 L/kg | 组织分布 |
| **BBB通透** | > -1 | 血脑屏障穿透 |
| **血浆蛋白结合** | 50-95% | 与蛋白结合率 |

#### M - Metabolism (代谢)
| 指标 | 标准 | 解释 |
|------|------|------|
| **CYP底物** | 检测5种CYP酶 | 代谢途径 |
| **CYP抑制** | 避免强抑制 | 药物相互作用风险 |

#### E - Excretion (排泄)
| 指标 | 标准 | 解释 |
|------|------|------|
| **清除率** | 5-30 mL/min/kg | 药物清除速度 |
| **半衰期** | 2-12 h | 药物作用时间 |

#### T - Toxicity (毒性)
| 指标 | 标准 | 解释 |
|------|------|------|
| **Ames致突变** | Negative | 无遗传毒性 |
| **hERG抑制** | Negative | 心脏安全性 |
| **肝毒性** | Negative | 肝脏安全性 |
| **LD50** | > 500 mg/kg | 急性毒性 |

**筛选标准**:
- ADMET综合评分 ≥ 50
- 无Ames致突变性
- 无hERG抑制

---

### 阶段3: 类药性评估

#### Lipinski五规则 (Ro5)
| 规则 | 标准 | 权重 |
|------|------|------|
| 分子量 | ≤ 500 Da | ⭐⭐⭐ |
| LogP | ≤ 5 | ⭐⭐⭐ |
| H键供体 | ≤ 5 | ⭐⭐ |
| H键受体 | ≤ 10 | ⭐⭐ |

**标准**: ≤ 1个违规

#### Veber规则
- 可旋转键 ≤ 10
- TPSA ≤ 140 Ų

#### Ghose规则
- MW: 160-480 Da
- LogP: -0.4 to 5.6
- 原子数: 20-70

#### PAINS过滤
- 泛活性化合物检测
- 消除假阳性
- 标准: 0个警报

**筛选标准**:
- 类药性评分 ≥ 50
- Lipinski违规 ≤ 1
- 无PAINS警报

---

### 阶段4: 分子对接

**工具**: AutoDock Vina

**参数**:
- Box大小: 20×20×20 Å
- Exhaustiveness: 8
- 对接模式: 9个

**评分**:
- Vina Score (kcal/mol)
- 越负越好
- 标准: < -6.0

**输出**:
- 结合亲和力
- 结合构象
- 结合位点

---

### 阶段5: 综合评分

**评分公式**:
```
Final Score = 
  ADMET Score × 0.35 +
  DrugLikeness Score × 0.25 +
  QED × 100 × 0.20 +
  |Vina Score| × 10 × 0.20
```

**权重说明**:
- **ADMET (35%)**: 药代动力学最重要
- **类药性 (25%)**: 成药可能性
- **QED (20%)**: 药物相似性
- **对接 (20%)**: 靶点亲和力

**推荐标准**:
- 综合评分 ≥ 60
- ADMET评分 ≥ 60
- Lipinski通过
- Vina评分 < -6.0

## 📈 完整示例

### 示例1: 癌症靶点筛选

```bash
# 针对EGFR激酶的药物筛选
python industrial_pipeline.py \
  --pdb EGFR_kinase.pdb \
  --candidates 500 \
  --sampling_steps 1000 \
  --device cuda \
  --output egfr_screening

# 预期结果:
# - 500个初始候选
# - ~200个通过ADMET
# - ~100个通过类药性
# - ~50个完成对接
# - Top 10推荐候选
```

### 示例2: 快速探索

```bash
# 快速评估(不含对接)
python industrial_pipeline.py \
  --pdb target.pdb \
  --candidates 200 \
  --sampling_steps 300 \
  --no_docking \
  --output quick_screening
```

### 示例3: 高质量生成

```bash
# 最高质量模式
python industrial_pipeline.py \
  --pdb target.pdb \
  --candidates 1000 \
  --num_atoms 20 40 \
  --sampling_steps 1500 \
  --device cuda \
  --output premium_results
```

## 📊 输出文件

```
industrial_results/
├── molecules/                    # 生成的分子SDF文件
│   ├── mol_0000.sdf
│   ├── mol_0001.sdf
│   └── ...
├── complete_results.json         # 完整结果(JSON)
├── complete_results.csv          # 完整结果(CSV)
├── top20_candidates.csv          # Top 20候选
└── analysis_report.json          # 分析报告
```

### complete_results.csv 字段

| 字段 | 说明 | 单位 |
|------|------|------|
| mol_id | 分子ID | - |
| smiles | SMILES字符串 | - |
| mw | 分子量 | Da |
| logp | 脂溶性 | - |
| qed | QED评分 | 0-1 |
| admet_score | ADMET评分 | 0-100 |
| caco2 | 通透性 | log cm/s |
| hia | 肠道吸收 | % |
| ames_toxic | 致突变性 | Bool |
| herg_risk | 心脏风险 | Bool |
| druglikeness_score | 类药性评分 | 0-100 |
| lipinski_pass | Lipinski通过 | Bool |
| pains_alerts | PAINS警报 | Int |
| vina_score | Vina评分 | kcal/mol |
| final_score | 综合评分 | 0-100 |
| rank | 排名 | Int |
| recommended | 推荐 | Bool |

## 🎯 结果分析

### 查看Top候选

```python
import pandas as pd

# 读取结果
df = pd.read_csv('industrial_results/top20_candidates.csv')

# 显示最佳候选
print(df[['mol_id', 'final_score', 'vina_score', 'admet_score', 
         'druglikeness_score', 'recommended']].head(10))
```

### 筛选推荐分子

```python
# 筛选推荐的分子
recommended = df[df['recommended'] == True]
print(f"推荐候选数: {len(recommended)}")

# 筛选特定条件
high_quality = df[
    (df['final_score'] >= 70) &
    (df['lipinski_pass'] == True) &
    (df['vina_score'] < -7.0) &
    (df['ames_toxic'] == False)
]
print(f"高质量候选数: {len(high_quality)}")
```

### 可视化分析

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. 综合评分分布
axes[0, 0].hist(df['final_score'], bins=20)
axes[0, 0].set_title('Final Score Distribution')
axes[0, 0].set_xlabel('Score')

# 2. ADMET vs 类药性
axes[0, 1].scatter(df['admet_score'], df['druglikeness_score'])
axes[0, 1].set_xlabel('ADMET Score')
axes[0, 1].set_ylabel('DrugLikeness Score')

# 3. Vina评分分布
axes[1, 0].hist(df['vina_score'], bins=20)
axes[1, 0].set_title('Vina Score Distribution')
axes[1, 0].set_xlabel('Vina Score (kcal/mol)')

# 4. MW vs LogP
axes[1, 1].scatter(df['mw'], df['logp'])
axes[1, 1].set_xlabel('Molecular Weight')
axes[1, 1].set_ylabel('LogP')

plt.tight_layout()
plt.savefig('analysis.png')
```

## ⚙️ 高级配置

### 自定义阈值

```python
from industrial_pipeline import IndustrialDrugDiscoveryPipeline

pipeline = IndustrialDrugDiscoveryPipeline()

# 修改阈值
pipeline.thresholds = {
    'admet_score': 60.0,  # 更严格
    'druglikeness_score': 55.0,
    'lipinski_violations': 0,  # 必须完全通过
    'pains_alerts': 0,
    'vina_score': -7.0,  # 更强结合
}

# 运行
profiles = pipeline.run_complete_pipeline(
    pdb_path='target.pdb',
    num_candidates=100
)
```

### Python API完整示例

```python
from industrial_pipeline import IndustrialDrugDiscoveryPipeline
from targetdiff_full import ModelConfig

# 自定义配置
config = ModelConfig()
config.num_diffusion_steps = 2000
config.sampling_type = 'ddim'

# 初始化
pipeline = IndustrialDrugDiscoveryPipeline(
    model_path='pretrained_model.pt',
    device='cuda',
    config=config
)

# 运行完整流程
profiles = pipeline.run_complete_pipeline(
    pdb_path='protein.pdb',
    num_candidates=500,
    num_atoms_range=(20, 40),
    sampling_steps=1000,
    output_dir='custom_results',
    enable_docking=True
)

# 分析结果
recommended = [p for p in profiles if p.recommended]
print(f"发现 {len(recommended)} 个推荐候选")

for i, p in enumerate(recommended[:5], 1):
    print(f"\n候选 {i}:")
    print(f"  SMILES: {p.smiles}")
    print(f"  综合评分: {p.final_score:.1f}")
    print(f"  ADMET评分: {p.admet_score:.1f}")
    print(f"  Vina评分: {p.vina_score:.2f} kcal/mol")
```

## 📊 性能基准

| 配置 | 候选数 | GPU | 时间 | 推荐数 | 成功率 |
|------|--------|-----|------|--------|--------|
| 快速 | 100 | RTX 3090 | ~30分钟 | 5-10 | 5-10% |
| 标准 | 500 | RTX 3090 | ~2小时 | 20-30 | 4-6% |
| 高质量 | 1000 | RTX 3090 | ~4小时 | 40-60 | 4-6% |

## 🎯 最佳实践

### 1. 初步筛选 (快速)

```bash
# 100个候选, 不含对接
python industrial_pipeline.py \
  --pdb target.pdb \
  --candidates 100 \
  --sampling_steps 300 \
  --no_docking
```

**用时**: ~10分钟  
**目的**: 快速评估可行性

### 2. 中等规模 (平衡)

```bash
# 500个候选, 含对接
python industrial_pipeline.py \
  --pdb target.pdb \
  --candidates 500 \
  --sampling_steps 500 \
  --device cuda
```

**用时**: ~2小时  
**目的**: 发现高质量候选

### 3. 深度挖掘 (彻底)

```bash
# 1000个候选, 高质量采样
python industrial_pipeline.py \
  --pdb target.pdb \
  --candidates 1000 \
  --sampling_steps 1000 \
  --device cuda
```

**用时**: ~4-6小时  
**目的**: 最大化候选数量和质量

## 🔬 科学验证

### 验证通过ADMET的分子

```bash
# 导出推荐分子的SDF
python -c "
import pandas as pd
from rdkit import Chem

df = pd.read_csv('industrial_results/complete_results.csv')
recommended = df[df['recommended'] == True]

for _, row in recommended.iterrows():
    print(f\"{row['mol_id']}: Score={row['final_score']:.1f}\")
"
```

### 实验验证建议

1. **先导化合物筛选**:
   - 选择Top 10-20候选
   - 进行体外活性测试
   - 验证结合亲和力

2. **ADMET实验验证**:
   - Caco-2细胞渗透实验
   - CYP酶抑制实验
   - Ames测试

3. **结构优化**:
   - 基于实验数据优化
   - 结构-活性关系分析
   - 再次迭代生成

## 🆘 故障排除

### Q: ADMET评分都很低?

**A**: 
1. 增加候选数 (`--candidates 500`)
2. 调整原子数范围
3. 降低ADMET阈值

### Q: 没有推荐的分子?

**A**: 
1. 放宽筛选标准
2. 检查PDB文件质量
3. 增加采样步数

### Q: 对接失败?

**A**: 
1. 检查Vina安装
2. 准备PDBQT文件
3. 或使用 `--no_docking`

## 📚 参考文献

1. **ADMET预测**: 
   - SwissADME
   - pkCSM
   - ADMETlab

2. **类药性规则**:
   - Lipinski's Rule of Five
   - Veber's Rule
   - PAINS filters

3. **分子对接**:
   - AutoDock Vina
   - Glide
   - Gold

---

**开始发现您的候选药物分子! 🚀**
