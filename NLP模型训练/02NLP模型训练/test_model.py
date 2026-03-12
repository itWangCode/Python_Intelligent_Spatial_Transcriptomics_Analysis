#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试已训练的BioSentVec模型
"""

from biosentvec_toolkit import BioSentVecToolkit
# 或使用简化版
# from simple_biosentvec import SimpleBioSentVec

print("=" * 70)
print("BioSentVec 模型测试")
print("=" * 70)

# 1. 加载模型
print("\n[1/4] 加载模型...")
toolkit = BioSentVecToolkit()
toolkit.load_word_model()  # 自动加载最新模型

# 2. 测试词语相似度
print("\n[2/4] 测试词语相似度")
print("-" * 70)

test_words = ['diabetes', 'patient', 'treatment', 'disease', 'cancer']

for word in test_words:
    print(f"\n'{word}' 的相似词:")
    similar = toolkit.find_similar_words(word, topn=5)
    if similar:
        for sim_word, score in similar:
            print(f"  {sim_word:20s} 相似度: {score:.4f}")
    else:
        print(f"  (词 '{word}' 不在词汇表中)")

# 3. 测试句子相似度
print("\n" + "=" * 70)
print("[3/4] 测试句子相似度")
print("-" * 70)

sentence_pairs = [
    # 高度相似的句子对
    (
        "The patient was diagnosed with diabetes",
        "Diabetes was diagnosed in the patient"
    ),
    (
        "Cancer treatment includes chemotherapy",
        "Chemotherapy is used for cancer treatment"
    ),
    # 中等相似
    (
        "Hypertension is a major risk factor",
        "High blood pressure can be dangerous"
    ),
    # 低相似
    (
        "The patient has diabetes",
        "Cancer treatment options"
    ),
]

for i, (sent1, sent2) in enumerate(sentence_pairs, 1):
    sim = toolkit.compute_similarity(sent1, sent2)
    
    # 评级
    if sim > 0.8:
        rating = "⭐⭐⭐⭐⭐ 高度相似"
    elif sim > 0.6:
        rating = "⭐⭐⭐⭐ 较为相似"
    elif sim > 0.4:
        rating = "⭐⭐⭐ 中等相似"
    elif sim > 0.2:
        rating = "⭐⭐ 稍微相关"
    else:
        rating = "⭐ 相关性低"
    
    print(f"\n配对 {i}:")
    print(f"  句子A: {sent1}")
    print(f"  句子B: {sent2}")
    print(f"  相似度: {sim:.4f} - {rating}")

# 4. 自定义测试
print("\n" + "=" * 70)
print("[4/4] 自定义测试")
print("-" * 70)
print("您可以输入自己的句子进行测试（输入 'q' 退出）")

while True:
    print("\n")
    sent1 = input("句子1: ").strip()
    if sent1.lower() == 'q':
        break
    
    sent2 = input("句子2: ").strip()
    if sent2.lower() == 'q':
        break
    
    if sent1 and sent2:
        sim = toolkit.compute_similarity(sent1, sent2)
        print(f"\n相似度: {sim:.4f}")
        
        if sim > 0.8:
            print("评价: 高度相似 ⭐⭐⭐⭐⭐")
        elif sim > 0.6:
            print("评价: 较为相似 ⭐⭐⭐⭐")
        elif sim > 0.4:
            print("评价: 中等相似 ⭐⭐⭐")
        elif sim > 0.2:
            print("评价: 稍微相关 ⭐⭐")
        else:
            print("评价: 相关性低 ⭐")

print("\n" + "=" * 70)
print("测试完成！")
print("=" * 70)
