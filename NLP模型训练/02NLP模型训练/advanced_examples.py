#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BioSentVec 高级使用示例
展示更多实用功能
"""

from biosentvec_toolkit import BioSentVecToolkit
import numpy as np
from pathlib import Path


def example_1_medical_document_clustering():
    """示例1: 医学文档聚类"""
    print("\n" + "=" * 70)
    print("示例 1: 医学文档聚类")
    print("=" * 70)
    
    toolkit = BioSentVecToolkit()
    
    # 准备示例医学文档
    documents = [
        "Type 2 diabetes mellitus is characterized by insulin resistance.",
        "Patients with diabetes often require blood glucose monitoring.",
        "Diabetic neuropathy is a common complication of diabetes.",
        "Hypertension increases the risk of stroke and heart disease.",
        "High blood pressure can be managed with lifestyle changes.",
        "Cardiovascular disease is a leading cause of mortality.",
        "Cancer cells exhibit uncontrolled growth and division.",
        "Chemotherapy is used to treat various types of cancer.",
        "Immunotherapy has shown promise in treating cancer.",
    ]
    
    # 创建临时训练数据
    train_file = Path('temp_train.txt')
    with open(train_file, 'w', encoding='utf-8') as f:
        for doc in documents:
            f.write(doc.lower() + '\n')
    
    # 训练模型
    print("\n训练模型...")
    toolkit.train_word_embeddings(
        corpus_file=train_file,
        vector_size=50,
        epochs=20,
        save_name='clustering_model'
    )
    
    # 计算文档向量
    print("\n计算文档向量...")
    vectors = []
    for doc in documents:
        vec = toolkit.get_sentence_vector(doc)
        vectors.append(vec)
    
    # 计算相似度矩阵
    print("\n文档相似度矩阵:")
    print("\n主题分类:")
    print("  文档 0-2: 糖尿病相关")
    print("  文档 3-5: 心血管疾病相关")
    print("  文档 6-8: 癌症相关")
    print()
    
    similarity_matrix = np.zeros((len(documents), len(documents)))
    for i in range(len(documents)):
        for j in range(len(documents)):
            if i == j:
                similarity_matrix[i][j] = 1.0
            else:
                vec_i = vectors[i]
                vec_j = vectors[j]
                sim = np.dot(vec_i, vec_j) / (np.linalg.norm(vec_i) * np.linalg.norm(vec_j))
                similarity_matrix[i][j] = sim
    
    # 打印矩阵
    print("     ", end="")
    for i in range(len(documents)):
        print(f"Doc{i:2d}", end="  ")
    print()
    
    for i in range(len(documents)):
        print(f"Doc{i:2d}", end="  ")
        for j in range(len(documents)):
            print(f"{similarity_matrix[i][j]:5.2f}", end="  ")
        print()
    
    # 清理
    train_file.unlink()
    
    print("\n观察: 同类别文档之间的相似度应该更高！")


def example_2_sentence_semantic_search():
    """示例2: 语义搜索"""
    print("\n" + "=" * 70)
    print("示例 2: 医学问答语义搜索")
    print("=" * 70)
    
    toolkit = BioSentVecToolkit()
    
    # 医学知识库
    knowledge_base = [
        "Diabetes is a metabolic disorder characterized by high blood sugar levels.",
        "Symptoms of diabetes include frequent urination, increased thirst, and fatigue.",
        "Type 1 diabetes is an autoimmune condition affecting insulin production.",
        "Type 2 diabetes is associated with insulin resistance and obesity.",
        "Hypertension is defined as persistently elevated blood pressure.",
        "Heart disease is the leading cause of death worldwide.",
        "Regular exercise can help prevent cardiovascular disease.",
        "Cancer involves the uncontrolled growth of abnormal cells.",
        "Chemotherapy uses drugs to kill cancer cells.",
        "Asthma is a chronic respiratory condition causing breathing difficulties.",
    ]
    
    # 创建训练文件
    train_file = Path('kb_train.txt')
    with open(train_file, 'w', encoding='utf-8') as f:
        for text in knowledge_base:
            f.write(text.lower() + '\n')
    
    # 训练
    print("\n建立知识库索引...")
    toolkit.train_word_embeddings(
        corpus_file=train_file,
        vector_size=100,
        epochs=15,
        save_name='search_model'
    )
    
    # 用户查询
    queries = [
        "What are the symptoms of high blood sugar?",
        "How to prevent heart problems?",
        "Tell me about cancer treatment options.",
    ]
    
    print("\n语义搜索结果:")
    print("=" * 70)
    
    for query in queries:
        print(f"\n查询: {query}")
        print("-" * 70)
        
        # 计算查询与知识库的相似度
        scores = []
        for kb_text in knowledge_base:
            sim = toolkit.compute_similarity(query, kb_text)
            scores.append((sim, kb_text))
        
        # 排序并显示前3个结果
        scores.sort(reverse=True)
        print("最相关的答案:")
        for i, (score, text) in enumerate(scores[:3], 1):
            print(f"  {i}. [相似度: {score:.4f}] {text}")
    
    train_file.unlink()


def example_3_medical_term_similarity():
    """示例3: 医学术语相似度分析"""
    print("\n" + "=" * 70)
    print("示例 3: 医学术语相似度分析")
    print("=" * 70)
    
    toolkit = BioSentVecToolkit()
    
    # 准备包含医学术语的句子
    medical_corpus = [
        "diabetes mellitus glucose insulin hyperglycemia",
        "hypertension blood pressure cardiovascular heart",
        "cancer tumor malignant metastasis oncology",
        "asthma respiratory bronchial inflammation",
        "arthritis joint inflammation rheumatoid",
        "diabetes insulin glucose pancreas",
        "hypertension arterial pressure systolic diastolic",
        "cancer chemotherapy radiation therapy",
        "asthma inhaler bronchodilator respiratory",
        "arthritis pain swelling joints",
    ] * 5  # 重复以增加训练数据
    
    train_file = Path('terms_train.txt')
    with open(train_file, 'w', encoding='utf-8') as f:
        for text in medical_corpus:
            f.write(text.lower() + '\n')
    
    print("\n训练医学术语模型...")
    toolkit.train_word_embeddings(
        corpus_file=train_file,
        vector_size=100,
        window=5,
        min_count=1,
        epochs=30,
        save_name='medical_terms_model'
    )
    
    # 测试术语相似度
    print("\n医学术语相似词分析:")
    print("=" * 70)
    
    test_terms = ['diabetes', 'cancer', 'hypertension', 'asthma']
    
    for term in test_terms:
        similar = toolkit.find_similar_words(term, topn=5)
        if similar:
            print(f"\n'{term}' 的相关术语:")
            for word, score in similar:
                print(f"  • {word:20s} 相似度: {score:.4f}")
    
    train_file.unlink()


def example_4_batch_similarity_computation():
    """示例4: 批量相似度计算"""
    print("\n" + "=" * 70)
    print("示例 4: 批量相似度计算与结果导出")
    print("=" * 70)
    
    toolkit = BioSentVecToolkit()
    train_file, _ = toolkit.prepare_sample_data()
    
    # 训练模型
    print("\n训练模型...")
    toolkit.train_word_embeddings(
        corpus_file=train_file,
        vector_size=100,
        epochs=10,
        save_name='batch_model'
    )
    
    # 批量句子对
    sentence_pairs = [
        ("The patient has diabetes", "Diabetes diagnosis in patient"),
        ("High blood pressure detected", "Hypertension observed"),
        ("Cancer treatment started", "Patient began chemotherapy"),
        ("Asthma symptoms worsened", "Respiratory problems increased"),
        ("Kidney function declining", "Renal impairment noted"),
    ]
    
    print("\n批量计算相似度...")
    results = []
    
    print("\n句子对相似度:")
    print("-" * 70)
    for i, (sent1, sent2) in enumerate(sentence_pairs, 1):
        sim = toolkit.compute_similarity(sent1, sent2)
        results.append({
            'pair_id': i,
            'sentence1': sent1,
            'sentence2': sent2,
            'similarity': sim
        })
        print(f"{i}. [{sim:.4f}] {sent1} <-> {sent2}")
    
    # 导出结果
    output_file = toolkit.output_dir / 'similarity_results.txt'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("句子对相似度分析结果\n")
        f.write("=" * 70 + "\n\n")
        for r in results:
            f.write(f"配对 {r['pair_id']}:\n")
            f.write(f"  句子1: {r['sentence1']}\n")
            f.write(f"  句子2: {r['sentence2']}\n")
            f.write(f"  相似度: {r['similarity']:.4f}\n\n")
    
    print(f"\n结果已保存至: {output_file}")


def main():
    """运行所有示例"""
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          BioSentVec 高级使用示例集                         ║
    ╚════════════════════════════════════════════════════════════╝
    
    包含以下示例:
    
    1. 医学文档聚类
       - 演示如何将相似的医学文档聚类
    
    2. 语义搜索
       - 构建医学知识库问答系统
    
    3. 医学术语相似度
       - 分析医学术语之间的关系
    
    4. 批量相似度计算
       - 批量处理句子对并导出结果
    """)
    
    choice = input("\n选择示例 (1/2/3/4) 或 'all' 运行全部: ").strip()
    
    if choice == '1':
        example_1_medical_document_clustering()
    elif choice == '2':
        example_2_sentence_semantic_search()
    elif choice == '3':
        example_3_medical_term_similarity()
    elif choice == '4':
        example_4_batch_similarity_computation()
    elif choice.lower() == 'all':
        example_1_medical_document_clustering()
        example_2_sentence_semantic_search()
        example_3_medical_term_similarity()
        example_4_batch_similarity_computation()
    else:
        print("无效选择！")
        return
    
    print("\n" + "=" * 70)
    print("示例运行完成!")
    print("=" * 70)


if __name__ == '__main__':
    main()
