#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
模型性能评估
评估模型在不同任务上的表现
"""

from biosentvec_toolkit import BioSentVecToolkit
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

print("=" * 70)
print("模型性能评估")
print("=" * 70)

# 加载模型
print("\n加载模型...")
toolkit = BioSentVecToolkit()
toolkit.load_word_model()

# 1. 词汇覆盖率测试
print("\n" + "=" * 70)
print("1. 词汇覆盖率测试")
print("=" * 70)

test_words = [
    'diabetes', 'cancer', 'hypertension', 'asthma', 'arthritis',
    'patient', 'treatment', 'disease', 'symptom', 'diagnosis',
    'medication', 'therapy', 'chronic', 'acute', 'clinical',
    'medical', 'health', 'care', 'study', 'research'
]

covered = 0
for word in test_words:
    if word in toolkit.word_model.wv:
        covered += 1

coverage = covered / len(test_words) * 100
print(f"测试词汇: {len(test_words)}")
print(f"覆盖词汇: {covered}")
print(f"覆盖率: {coverage:.1f}%")

# 2. 句子相似度任务
print("\n" + "=" * 70)
print("2. 句子相似度任务")
print("=" * 70)

# 医学句子对（人工标注的相似度 0-1）
sentence_pairs = [
    ("patient has diabetes", "diabetes diagnosed", 0.9),
    ("cancer treatment", "chemotherapy", 0.8),
    ("high blood pressure", "hypertension", 0.9),
    ("heart attack", "myocardial infarction", 0.95),
    ("diabetes", "cancer", 0.3),
    ("treatment", "weather", 0.1),
]

human_scores = []
model_scores = []

print("\n句子对相似度比较:")
print("-" * 70)

for sent1, sent2, human_score in sentence_pairs:
    model_score = toolkit.compute_similarity(sent1, sent2)
    human_scores.append(human_score)
    model_scores.append(model_score)
    
    diff = abs(human_score - model_score)
    print(f"\n句子1: {sent1}")
    print(f"句子2: {sent2}")
    print(f"人工评分: {human_score:.2f}")
    print(f"模型评分: {model_score:.4f}")
    print(f"差异: {diff:.4f}")

# 计算相关性
if len(human_scores) > 1:
    correlation, p_value = spearmanr(human_scores, model_scores)
    print(f"\nSpearman相关系数: {correlation:.4f}")
    print(f"P值: {p_value:.4f}")

# 3. 语义聚类测试
print("\n" + "=" * 70)
print("3. 语义聚类测试")
print("=" * 70)

# 三组医学概念
groups = {
    '糖尿病': ['diabetes', 'glucose', 'insulin', 'hyperglycemia'],
    '心血管': ['hypertension', 'heart', 'cardiovascular', 'cardiac'],
    '癌症': ['cancer', 'tumor', 'oncology', 'malignancy']
}

print("\n组内平均相似度:")
for group_name, words in groups.items():
    similarities = []
    for i, w1 in enumerate(words):
        for w2 in words[i+1:]:
            if w1 in toolkit.word_model.wv and w2 in toolkit.word_model.wv:
                sim = toolkit.word_model.wv.similarity(w1, w2)
                similarities.append(sim)
    
    if similarities:
        avg_sim = np.mean(similarities)
        print(f"  {group_name}: {avg_sim:.4f}")

# 4. 模型统计信息
print("\n" + "=" * 70)
print("4. 模型统计信息")
print("=" * 70)

print(f"词汇表大小: {len(toolkit.word_model.wv)}")
print(f"向量维度: {toolkit.word_model.vector_size}")
print(f"窗口大小: {toolkit.word_model.window}")

# 词频统计
if hasattr(toolkit.word_model.wv, 'index_to_key'):
    top_words = toolkit.word_model.wv.index_to_key[:20]
    print(f"\n前20个高频词:")
    for i, word in enumerate(top_words, 1):
        print(f"  {i:2d}. {word}")

print("\n" + "=" * 70)
print("评估完成！")
print("=" * 70)
