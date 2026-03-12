# BioSentVec 完整工具包使用指南

## 您获得的完整文件包

```
BioSentVec_Complete_Toolkit/
├──  核心文件
│   ├── biosentvec_toolkit.py      # 主工具包（16KB）
│   ├── requirements.txt            # 依赖列表
│   └── README.md                   # 完整文档（10KB）
│
├──  快速开始
│   ├── quick_start.py             # 快速开始脚本（7KB）
│   ├── install.sh                 # Linux/Mac安装脚本
│   └── install.bat                # Windows安装脚本
│
├──  高级功能
│   ├── advanced_examples.py       # 高级示例（11KB）
│   ├── data_preparation.py        # 数据准备工具（11KB）
│   └── test_toolkit.py           # 测试套件（10KB）
│
└──  文档
    ├── QUICK_REFERENCE.md         # 快速参考（5KB）
    └── START_HERE.md             # 本文件
```

## 10秒快速开始

### Linux / Mac
```bash
# 1. 安装
bash install.sh

# 2. 运行
python quick_start.py
```

### Windows
```cmd
# 1. 安装
install.bat

# 2. 运行
python quick_start.py
```

## 三种使用模式

### 模式1️⃣: 菜单模式（最简单）
```bash
python quick_start.py
```
然后按提示选择：
- 选项1: 自动演示所有功能
- 选项2: 交互式计算相似度

### 模式2️⃣: 自动演示
```bash
python quick_start.py quick
```
自动完成：数据准备 → 训练 → 测试 → 结果

### 模式3️⃣: 自定义编程
```python
from biosentvec_toolkit import BioSentVecToolkit

toolkit = BioSentVecToolkit()
train_file, _ = toolkit.prepare_sample_data()
toolkit.train_word_embeddings(corpus_file=train_file)
sim = toolkit.compute_similarity("句子1", "句子2")
```

## 特性清单

**自动GPU/CPU检测** - 有GPU用GPU，无GPU自动多线程  
**零配置启动** - 内置示例数据，开箱即用  
**完整功能** - 训练、相似度计算、词语分析  
**易于扩展** - 使用自己的医学文本数据  
**详细文档** - README + 快速参考 + 代码注释  
**测试套件** - 自动验证所有功能  

## 核心功能演示

### 1. 句子相似度计算
```python
toolkit = BioSentVecToolkit()
toolkit.load_word_model()  # 加载已训练模型

similarity = toolkit.compute_similarity(
    "The patient has diabetes",
    "Diabetes was diagnosed in the patient"
)
# 输出: 0.8234 (高度相似)
```

### 2. 医学术语分析
```python
similar_words = toolkit.find_similar_words('diabetes', topn=5)
# 输出:
# mellitus: 0.8921
# glucose: 0.8234
# insulin: 0.7845
# ...
```

### 3. 批量处理
```python
pairs = [
    ("text1", "text2"),
    ("text3", "text4")
]

for t1, t2 in pairs:
    sim = toolkit.compute_similarity(t1, t2)
    print(f"{sim:.4f}")
```

## 支持的数据源

### 推荐数据集

1. **PubMed文献摘要**
   - 官方来源，权威可靠
   - 包含数百万医学文献

2. **MIMIC-III临床笔记**
   - 真实临床数据
   - 需申请访问权限

3. **自定义医学文本**
   - 医院病历（需脱敏）
   - 医学论文
   - 临床指南

### 数据格式要求

**训练数据文件 (train.txt)**
```text
The patient was diagnosed with type 2 diabetes.
Hypertension is a major risk factor.
Cancer treatment includes chemotherapy.
...
```

- 每行一个句子
- UTF-8编码
- 建议：至少1000句以上

## 系统要求

### 最低配置
- Python 3.8+
- 4GB RAM
- 2GB 磁盘空间

### 推荐配置
- Python 3.10+
- 8GB+ RAM
- NVIDIA GPU（可选，但大幅提速）

### 依赖包
```
torch >= 2.0.0
numpy >= 1.24.0
gensim >= 4.3.0
nltk >= 3.8.0
tqdm >= 4.65.0
```

## 学习路径

### 新手（10分钟）
1. 运行 `python test_toolkit.py` 验证安装
2. 运行 `python quick_start.py` 体验功能
3. 查看 `QUICK_REFERENCE.md` 了解常用代码

### 进阶（30分钟）
1. 运行 `python advanced_examples.py` 学习高级用法
2. 阅读 `README.md` 了解参数调优
3. 尝试使用自己的数据训练

### 专家（1小时）
1. 阅读 `biosentvec_toolkit.py` 源代码
2. 自定义扩展功能
3. 集成到自己的项目中

## 典型应用场景

### 场景1: 医学问答系统
```python
# 构建知识库索引
knowledge_base = ["答案1", "答案2", ...]

# 用户提问
query = "什么是糖尿病症状？"

# 找最相似答案
scores = [(toolkit.compute_similarity(query, ans), ans) 
          for ans in knowledge_base]
best_answer = max(scores)[1]
```

### 场景2: 文档聚类
```python
# 计算所有文档对的相似度
docs = ["文档1", "文档2", ...]
matrix = [[toolkit.compute_similarity(d1, d2) 
           for d2 in docs] for d1 in docs]

# 聚类分析
# ...
```

### 场景3: 语义搜索
```python
# 搜索最相关的医学文献
query = "糖尿病治疗方案"
results = sorted(documents, 
    key=lambda d: toolkit.compute_similarity(query, d),
    reverse=True
)
```

## 专业提示

### 提高准确度
1. **使用更多数据** - 至少10万句以上
2. **领域专业数据** - 使用医学专业文本
3. **调整参数** - 增加 vector_size 到 200-300
4. **多轮训练** - epochs 设为 10-20

### 提升性能
1. **GPU加速** - 使用CUDA，速度提升5-10倍
2. **CPU多线程** - 已自动启用
3. **批处理** - 一次处理多个句子对
4. **模型缓存** - 避免重复加载模型

### 避免错误
1. **数据编码** - 确保UTF-8编码
2. **内存管理** - 大数据集分批处理
3. **模型保存** - 及时保存训练好的模型
4. **参数验证** - 使用测试脚本验证

## 故障排除

### 常见问题

**Q: 安装依赖失败？**
```bash
# 尝试使用国内镜像
pip install -r requirements.txt -i https://pypi.tuna.tsinghua.edu.cn/simple
```

**Q: 内存不足？**
```python
# 减小参数
toolkit.train_word_embeddings(
    corpus_file=file,
    vector_size=50,  # 从200降到50
    min_count=10     # 过滤低频词
)
```

**Q: 训练太慢？**
- 检查是否启用了多线程（自动）
- 考虑使用GPU
- 减少数据量或epochs

**Q: 相似度不准？**
- 增加训练数据
- 使用领域相关数据
- 调整 window 和 vector_size

## 技术支持

### 获取帮助
1. **查看文档**: `README.md`, `QUICK_REFERENCE.md`
2. **运行测试**: `python test_toolkit.py`
3. **查看示例**: `advanced_examples.py`

### 文件说明
- `biosentvec_toolkit.py` - 完整工具包实现
- `quick_start.py` - 快速开始和演示
- `test_toolkit.py` - 功能测试验证
- `README.md` - 完整使用文档
- `QUICK_REFERENCE.md` - 快速参考手册

## 开始您的旅程

### 第一步：验证安装
```bash
python test_toolkit.py
```

### 第二步：快速体验
```bash
python quick_start.py
```

### 第三步：开始创造
```python
from biosentvec_toolkit import BioSentVecToolkit
# 开始您的项目！
```
