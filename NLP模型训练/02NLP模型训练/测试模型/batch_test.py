#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量测试模型性能
从文件读取测试数据
"""

from biosentvec_toolkit import BioSentVecToolkit
import pandas as pd
from pathlib import Path

print("=" * 70)
print("批量测试模型")
print("=" * 70)

# 加载模型
print("\n加载模型...")
toolkit = BioSentVecToolkit()
toolkit.load_word_model()

# 测试数据
test_data = [
    # 句子1, 句子2, 预期相似度（高/中/低）
    ("patient has diabetes", "diabetes diagnosed in patient", "高"),
    ("cancer treatment started", "patient began chemotherapy", "高"),
    ("high blood pressure", "hypertension detected", "高"),
    ("kidney function declining", "renal impairment noted", "高"),
    
    ("diabetes management", "cancer treatment", "中"),
    ("heart disease", "lung disease", "中"),
    
    ("diabetes", "completely unrelated topic", "低"),
    ("cancer", "weather forecast", "低"),
]

print(f"\n测试 {len(test_data)} 个句子对...")
print("=" * 70)

results = []

for i, (sent1, sent2, expected) in enumerate(test_data, 1):
    sim = toolkit.compute_similarity(sent1, sent2)
    
    # 判断预测是否正确
    if expected == "高" and sim > 0.6:
        status = "✓"
    elif expected == "中" and 0.3 < sim <= 0.6:
        status = "✓"
    elif expected == "低" and sim <= 0.3:
        status = "✓"
    else:
        status = "✗"
    
    results.append({
        'ID': i,
        '句子1': sent1,
        '句子2': sent2,
        '预期': expected,
        '相似度': sim,
        '状态': status
    })
    
    print(f"\n测试 {i}:")
    print(f"  句子A: {sent1}")
    print(f"  句子B: {sent2}")
    print(f"  预期: {expected}相似, 实际: {sim:.4f} {status}")

# 统计
print("\n" + "=" * 70)
print("测试统计")
print("=" * 70)

correct = sum(1 for r in results if r['状态'] == '✓')
total = len(results)
accuracy = correct / total * 100

print(f"总测试数: {total}")
print(f"正确数: {correct}")
print(f"准确率: {accuracy:.1f}%")

# 保存结果
output_file = 'test_results.csv'
df = pd.DataFrame(results)
df.to_csv(output_file, index=False, encoding='utf-8')
print(f"\n结果已保存至: {output_file}")
