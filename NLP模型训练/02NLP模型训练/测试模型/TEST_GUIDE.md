# 🧪 模型测试完全指南

训练完成后，使用以下4种方式测试模型！

---

## ✅ 你的训练已成功！

看到这些信息说明训练成功：
```
模型已保存: biosentvec_data/models/biowordvec_custom.model
词向量已保存: biosentvec_data/models/biowordvec_custom.wv
训练完成!
```

---

## 🎯 4种测试方式

### 方式1: 快速测试（推荐新手）

```bash
python tests/test_model.py
```

**包含：**
- ✅ 词语相似度测试
- ✅ 句子相似度测试
- ✅ 交互式输入测试

**输出示例：**
```
[2/4] 测试词语相似度
'diabetes' 的相似词:
  mellitus              相似度: 0.8921
  glucose               相似度: 0.8234
  insulin               相似度: 0.7845

[3/4] 测试句子相似度
配对 1:
  句子A: The patient was diagnosed with diabetes
  句子B: Diabetes was diagnosed in the patient
  相似度: 0.9151 - ⭐⭐⭐⭐⭐ 高度相似
```

---

### 方式2: 交互式测试（最灵活）

```bash
python tests/interactive_test.py
```

**功能菜单：**
1. 计算句子相似度（实时输入）
2. 查找相似词
3. 批量测试
4. 退出

**使用示例：**
```
请选择功能 (1/2/3/4): 1

--- 句子相似度计算 ---
句子1: patient has diabetes
句子2: diabetes diagnosed in patient

相似度: 0.9151
图示: [█████████████████████████████████████████████░] 91.5%
评价: ⭐⭐⭐⭐⭐ 高度相似
```

---

### 方式3: 批量测试（测试集评估）

```bash
python tests/batch_test.py
```

**自动测试多个句子对，输出准确率**

**输出示例：**
```
测试统计
总测试数: 8
正确数: 7
准确率: 87.5%

结果已保存至: test_results.csv
```

---

### 方式4: 性能评估（专业用户）

```bash
python tests/evaluate_model.py
```

**包含：**
- 词汇覆盖率
- 相关性分析
- 语义聚类测试
- 模型统计信息

---

## 💻 代码中使用

### 基础使用

```python
from biosentvec_toolkit import BioSentVecToolkit

# 加载模型
toolkit = BioSentVecToolkit()
toolkit.load_word_model()  # 自动加载最新模型

# 测试1: 句子相似度
sim = toolkit.compute_similarity(
    "patient has diabetes",
    "diabetes diagnosed"
)
print(f"相似度: {sim:.4f}")

# 测试2: 查找相似词
similar = toolkit.find_similar_words('diabetes', topn=5)
for word, score in similar:
    print(f"{word}: {score:.4f}")
```

---

## 📊 测试自己的数据

### 创建测试文件

```python
# test_my_data.py
from biosentvec_toolkit import BioSentVecToolkit

toolkit = BioSentVecToolkit()
toolkit.load_word_model()

# 你的句子对
my_pairs = [
    ("句子1", "句子2"),
    ("句子3", "句子4"),
]

for s1, s2 in my_pairs:
    sim = toolkit.compute_similarity(s1, s2)
    print(f"[{sim:.4f}] {s1} <-> {s2}")
```

---

## 🔍 查看模型信息

```python
from biosentvec_toolkit import BioSentVecToolkit

toolkit = BioSentVecToolkit()
model = toolkit.load_word_model()

print(f"词汇表大小: {len(model.wv)}")
print(f"向量维度: {model.vector_size}")
print(f"窗口大小: {model.window}")

# 检查某个词是否在词汇表中
if 'diabetes' in model.wv:
    print("✓ 'diabetes' 在词汇表中")
```

---

## 📈 性能基准

根据你的训练结果 `相似度: 0.9151`，这是**非常好的结果**！

**参考标准：**
- 0.9+ = 优秀 ⭐⭐⭐⭐⭐
- 0.8-0.9 = 良好 ⭐⭐⭐⭐
- 0.6-0.8 = 中等 ⭐⭐⭐
- <0.6 = 需要改进 ⭐⭐

---

## 💡 测试技巧

### 1. 测试不同类型的句子对

```python
# 高度相似
("patient has diabetes", "diabetes diagnosed")  # 应该 >0.8

# 中等相似（同领域不同概念）
("diabetes treatment", "cancer therapy")  # 应该 0.4-0.7

# 低相似（不相关）
("diabetes", "weather")  # 应该 <0.3
```

### 2. 测试医学术语

```python
medical_terms = [
    'diabetes', 'hypertension', 'cancer', 
    'treatment', 'diagnosis', 'therapy'
]

for term in medical_terms:
    similar = toolkit.find_similar_words(term, topn=3)
    print(f"{term}: {similar}")
```

### 3. 测试边界情况

```python
# 空句子
sim = toolkit.compute_similarity("", "test")

# 单词句子
sim = toolkit.compute_similarity("diabetes", "cancer")

# 长句子
sim = toolkit.compute_similarity(
    "The patient was diagnosed with type 2 diabetes...",
    "Another long sentence..."
)
```

---

## 🎯 快速开始

**最简单的方式：**

```bash
# 1. 运行快速测试
python tests/test_model.py

# 2. 如果想交互式测试
python tests/interactive_test.py
```

---

## 📝 保存测试结果

```python
# 保存到文件
with open('test_results.txt', 'w') as f:
    for s1, s2 in pairs:
        sim = toolkit.compute_similarity(s1, s2)
        f.write(f"{s1}\t{s2}\t{sim:.4f}\n")
```

---

## ⚡ 常见问题

**Q: 如何知道模型训练得好不好？**

A: 运行 `python tests/evaluate_model.py` 查看：
- 词汇覆盖率 >80% = 好
- 句子相似度相关性 >0.7 = 好

**Q: 相似度总是很低怎么办？**

A: 可能原因：
- 训练数据太少（增加数据）
- 词汇不在词汇表（检查是否在训练数据中）
- 需要更多训练轮数（增加epochs）

**Q: 如何测试特定领域？**

A: 创建该领域的测试句子对，然后运行批量测试。

---

**开始测试你的模型吧！** 🚀
