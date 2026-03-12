# BioSentVec 完整工具包

这是一个**易于使用**的BioSentVec（生物医学句子向量）工具包，支持自动GPU/CPU检测、多线程训练，适合小白用户快速上手。

## 主要特性

**自动设备检测**: 自动检测并使用GPU（如果可用），否则使用CPU多线程  
**多线程支持**: CPU模式自动启用多线程加速  
**简单易用**: 提供快速开始脚本，无需复杂配置  
**完整功能**: 支持词向量训练、句子相似度计算、词语相似度查询  
**示例数据**: 内置医学文本示例，开箱即用  

## 系统要求

- **Python**: 3.8 或更高版本
- **内存**: 最低 4GB RAM（推荐 8GB+）
- **磁盘**: 至少 2GB 可用空间
- **可选**: NVIDIA GPU（支持CUDA加速）

## 快速安装

### 步骤 1: 安装Python依赖

```bash
# 安装所有依赖
pip install -r requirements.txt

# 或者手动安装
pip install torch numpy gensim nltk tqdm
```

### 步骤 2: 验证安装

```bash
python -c "import torch; print('GPU可用:', torch.cuda.is_available())"
```

## 快速开始（三种方式）

### 方式 1: 菜单模式（最简单，推荐新手）

```bash
python quick_start.py
```

然后选择：
- 选项 1: 快速开始（自动运行完整流程）
- 选项 2: 交互模式（输入自定义句子）

### 方式 2: 快速开始模式

```bash
python quick_start.py quick
```

自动执行：
1. 准备示例数据
2. 训练词向量模型
3. 测试词语相似度
4. 测试句子相似度

### 方式 3: 交互模式

```bash
python quick_start.py interactive
```

实时输入句子，计算相似度。

## 详细使用说明

### 1. 基础使用示例

```python
from biosentvec_toolkit import BioSentVecToolkit

# 初始化工具包
toolkit = BioSentVecToolkit()

# 准备数据（使用示例数据）
train_file, test_file = toolkit.prepare_sample_data()

# 训练模型
toolkit.train_word_embeddings(
    corpus_file=train_file,
    vector_size=100,      # 词向量维度
    window=10,            # 上下文窗口大小
    epochs=10,            # 训练轮数
    save_name='my_model'  # 模型名称
)

# 计算句子相似度
similarity = toolkit.compute_similarity(
    "The patient has diabetes",
    "Diabetes was diagnosed in the patient"
)
print(f"相似度: {similarity:.4f}")
```

### 2. 使用自己的数据训练

```python
from biosentvec_toolkit import BioSentVecToolkit

toolkit = BioSentVecToolkit()

# 训练模型（使用您自己的文本文件）
toolkit.train_word_embeddings(
    corpus_file='path/to/your/medical_texts.txt',  # 每行一个句子
    vector_size=200,      # 更大的维度获得更好的效果
    window=20,            # PubMed推荐使用20
    min_count=5,          # 最小词频
    epochs=10,            # 训练轮数
    save_name='custom_biowordvec'
)
```

### 3. 加载已训练的模型

```python
toolkit = BioSentVecToolkit()

# 加载最新训练的模型
toolkit.load_word_model()

# 或加载指定模型
toolkit.load_word_model('biosentvec_workspace/models/my_model.model')

# 使用模型
similarity = toolkit.compute_similarity(
    "Cancer treatment",
    "Oncology therapy"
)
```

### 4. 查找相似词

```python
# 查找与某个词最相似的词
similar_words = toolkit.find_similar_words('diabetes', topn=10)

for word, score in similar_words:
    print(f"{word}: {score:.4f}")
```

### 5. 下载预训练模型（可选）

```python
# 下载NCBI官方预训练的BioWordVec模型（13GB）
toolkit.download_pretrained_model(model_type='word')

# 下载BioSentVec模型（21GB）
toolkit.download_pretrained_model(model_type='sentence')
```

**注意**: 预训练模型非常大，下载可能需要较长时间。

## 数据集说明

### 训练数据格式

训练数据应该是**纯文本文件**，每行一个句子。例如：

```text
The patient was diagnosed with type 2 diabetes mellitus.
Hypertension is a major risk factor for cardiovascular disease.
The study investigated the efficacy of novel cancer treatments.
...
```

### 推荐的医学文本数据源

1. **PubMed摘要**
   - 下载地址: https://pubmed.ncbi.nlm.nih.gov/
   - 包含数千万篇生物医学文献摘要

2. **MIMIC-III临床笔记**
   - 地址: https://physionet.org/content/mimiciii/
   - 需要申请访问权限

3. **医学论文和报告**
   - 您自己的医学文本数据
   - 临床记录、病历等

### 示例数据

工具包内置了15条训练句子和3条测试句子作为示例：

**训练数据示例**:
```
- The patient was diagnosed with type 2 diabetes mellitus.
- Hypertension is a major risk factor for cardiovascular disease.
- The study investigated the efficacy of novel cancer treatments.
- ...（共15条）
```

**测试数据示例**:
```
- Patients with diabetes require regular monitoring.
- High blood pressure increases cardiovascular risk.
- Novel therapies show promise in treating malignancies.
```

## 设备检测与性能优化

### 自动GPU检测

工具包会自动检测是否有可用的GPU：

```
检测到GPU: NVIDIA GeForce RTX 3080
设备: GPU (CUDA 11.8)
```

### CPU多线程

如果没有GPU，工具会自动启用CPU多线程：

```
未检测到GPU，将使用CPU
CPU核心数: 16
CPU模式：将使用 15 个线程进行并行处理
```

### 性能提示

- **GPU模式**: 适合大规模数据训练
- **CPU模式**: 自动使用 (CPU核心数 - 1) 个线程
- **推荐配置**: 
  - 小数据集（<10万句子）: CPU多线程即可
  - 大数据集（>100万句子）: 建议使用GPU

## 参数调优指南

### 词向量训练参数

| 参数 | 说明 | 推荐值 | 范围 |
|------|------|--------|------|
| `vector_size` | 词向量维度 | 200 | 50-700 |
| `window` | 上下文窗口大小 | 10-20 | 5-30 |
| `min_count` | 最小词频 | 5 | 1-10 |
| `epochs` | 训练轮数 | 5-10 | 3-20 |
| `sg` | 训练算法 | 0 (CBOW) | 0或1 |

### 训练算法选择

- **CBOW (sg=0)**: 
  - 更快，适合大数据集
  - 对高频词效果好
  
- **Skip-gram (sg=1)**: 
  - 更慢，但对低频词效果更好
  - 适合小数据集

## 文件结构

```
biosentvec_workspace/           # 工作目录
├── data/                       # 数据目录
│   ├── train_data.txt         # 训练数据
│   └── test_data.txt          # 测试数据
├── models/                     # 模型目录
│   ├── biowordvec_demo.model  # 训练的模型
│   └── biowordvec_demo.wv     # 词向量文件
└── outputs/                    # 输出目录
```

## 常见问题

### Q1: 如何判断训练效果好不好？

**答**: 观察以下指标：
1. 训练loss是否持续下降
2. 相似词是否合理（使用`find_similar_words`测试）
3. 句子相似度是否符合预期

### Q2: 训练需要多长时间？

**答**: 取决于数据量和硬件：
- 示例数据（15句）: 几秒钟
- 10万句子: 5-30分钟（CPU）/ 2-5分钟（GPU）
- 100万句子: 1-3小时（CPU）/ 10-30分钟（GPU）

### Q3: 内存不够怎么办？

**答**: 尝试：
1. 减小`vector_size`
2. 增大`min_count`（过滤低频词）
3. 分批训练数据

### Q4: 如何提高准确度？

**答**:
1. 使用更多训练数据
2. 增加`vector_size`（如200或300）
3. 调整`window`大小
4. 增加训练`epochs`
5. 使用领域相关的数据

### Q5: CPU训练太慢怎么办？

**答**: 工具已自动启用多线程。如果仍然慢：
1. 减小数据量
2. 减小`vector_size`
3. 减少`epochs`
4. 考虑使用GPU或云服务

## 完整示例代码

### 示例1: 医学文档相似度

```python
from biosentvec_toolkit import BioSentVecToolkit

# 初始化
toolkit = BioSentVecToolkit()

# 准备您的医学文本数据
# train_data.txt 内容示例:
# Type 2 diabetes is a chronic condition.
# Insulin resistance is common in diabetes.
# ...

# 训练模型
toolkit.train_word_embeddings(
    corpus_file='medical_corpus.txt',
    vector_size=200,
    window=15,
    epochs=10
)

# 计算文档相似度
doc1 = "Patient presents with elevated blood glucose levels."
doc2 = "The patient has high blood sugar."
doc3 = "Treatment involves dietary modifications."

sim_12 = toolkit.compute_similarity(doc1, doc2)
sim_13 = toolkit.compute_similarity(doc1, doc3)

print(f"文档1 vs 文档2: {sim_12:.4f}")  # 应该较高
print(f"文档1 vs 文档3: {sim_13:.4f}")  # 应该较低
```

### 示例2: 批量处理

```python
import pandas as pd
from biosentvec_toolkit import BioSentVecToolkit

toolkit = BioSentVecToolkit()
toolkit.load_word_model()

# 读取CSV数据
df = pd.read_csv('medical_questions.csv')

# 批量计算相似度
results = []
for i, row in df.iterrows():
    sim = toolkit.compute_similarity(
        row['question1'],
        row['question2']
    )
    results.append(sim)

df['similarity'] = results
df.to_csv('results.csv', index=False)
```

## 参考文献

1. Chen Q, Peng Y, Lu Z. "BioSentVec: creating sentence embeddings for biomedical texts." IEEE International Conference on Healthcare Informatics. 2019.

2. Zhang Y, Chen Q, Yang Z, Lin H, Lu Z. "BioWordVec, improving biomedical word embeddings with subword information and MeSH." Scientific Data. 2019.

## 开始使用

现在就开始吧！

```bash
# 方法1: 最简单
python quick_start.py

# 方法2: 直接快速开始
python quick_start.py quick

# 方法3: 交互模式
python quick_start.py interactive
```

