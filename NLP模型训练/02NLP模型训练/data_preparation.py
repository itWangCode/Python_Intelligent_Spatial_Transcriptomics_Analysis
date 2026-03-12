#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
数据准备工具
帮助用户准备和验证训练数据
"""

import os
from pathlib import Path
import nltk
from collections import Counter
import matplotlib.pyplot as plt


class DataPreparationTool:
    """数据准备和分析工具"""
    
    def __init__(self):
        self.setup_nltk()
    
    def setup_nltk(self):
        """设置NLTK"""
        try:
            nltk.data.find('tokenizers/punkt')
        except LookupError:
            print("下载NLTK数据...")
            nltk.download('punkt', quiet=True)
    
    def validate_corpus_file(self, file_path):
        """
        验证语料文件
        
        Args:
            file_path: 文件路径
            
        Returns:
            统计信息字典
        """
        print(f"\n分析文件: {file_path}")
        print("=" * 70)
        
        if not Path(file_path).exists():
            print(f"错误: 文件不存在 - {file_path}")
            return None
        
        # 读取文件
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        # 统计信息
        stats = {
            'total_lines': len(lines),
            'empty_lines': 0,
            'total_tokens': 0,
            'unique_tokens': set(),
            'sentence_lengths': [],
            'token_frequencies': Counter()
        }
        
        # 分析每行
        for line in lines:
            line = line.strip()
            
            if not line:
                stats['empty_lines'] += 1
                continue
            
            # 分词
            tokens = nltk.word_tokenize(line.lower())
            stats['total_tokens'] += len(tokens)
            stats['unique_tokens'].update(tokens)
            stats['sentence_lengths'].append(len(tokens))
            stats['token_frequencies'].update(tokens)
        
        # 计算统计量
        stats['unique_tokens'] = len(stats['unique_tokens'])
        stats['avg_sentence_length'] = (
            sum(stats['sentence_lengths']) / len(stats['sentence_lengths'])
            if stats['sentence_lengths'] else 0
        )
        stats['min_sentence_length'] = min(stats['sentence_lengths']) if stats['sentence_lengths'] else 0
        stats['max_sentence_length'] = max(stats['sentence_lengths']) if stats['sentence_lengths'] else 0
        
        # 打印统计信息
        print(f"\n文件统计:")
        print(f"  总行数: {stats['total_lines']}")
        print(f"  空行数: {stats['empty_lines']}")
        print(f"  有效句子: {stats['total_lines'] - stats['empty_lines']}")
        print(f"\n词汇统计:")
        print(f"  总词数: {stats['total_tokens']}")
        print(f"  唯一词数: {stats['unique_tokens']}")
        print(f"  词汇丰富度: {stats['unique_tokens'] / stats['total_tokens']:.4f}")
        print(f"\n句子长度:")
        print(f"  平均: {stats['avg_sentence_length']:.2f} 词")
        print(f"  最短: {stats['min_sentence_length']} 词")
        print(f"  最长: {stats['max_sentence_length']} 词")
        
        # 显示高频词
        print(f"\n前20个高频词:")
        for word, freq in stats['token_frequencies'].most_common(20):
            print(f"  {word:20s} {freq:5d}")
        
        # 建议
        print(f"\n数据质量评估:")
        if stats['total_tokens'] < 10000:
            print("  ⚠ 警告: 数据量较小，建议至少10万词以上")
        elif stats['total_tokens'] < 100000:
            print("  ⚠ 提示: 数据量适中，可用于快速实验")
        else:
            print("  ✓ 数据量充足，适合训练")
        
        if stats['unique_tokens'] / stats['total_tokens'] < 0.1:
            print("  ⚠ 词汇丰富度较低，可能需要更多样化的文本")
        else:
            print("  ✓ 词汇丰富度良好")
        
        if stats['avg_sentence_length'] < 5:
            print("  ⚠ 句子平均长度较短")
        elif stats['avg_sentence_length'] > 50:
            print("  ⚠ 句子平均长度较长，建议分句")
        else:
            print("  ✓ 句子长度适中")
        
        return stats
    
    def create_sample_medical_corpus(self, output_file='medical_corpus.txt', num_sentences=100):
        """
        创建示例医学语料
        
        Args:
            output_file: 输出文件名
            num_sentences: 句子数量
        """
        print(f"\n创建示例医学语料...")
        
        # 医学模板
        templates = [
            "The patient presented with {symptom} and {symptom}.",
            "{disease} is characterized by {symptom}.",
            "Treatment for {disease} includes {treatment}.",
            "{symptom} is a common symptom of {disease}.",
            "The physician diagnosed the patient with {disease}.",
            "{treatment} is effective in managing {disease}.",
            "Risk factors for {disease} include {risk_factor}.",
            "The patient's {symptom} improved with {treatment}.",
            "{disease} can lead to complications such as {complication}.",
            "Regular {test} is recommended for monitoring {disease}.",
        ]
        
        # 医学术语
        symptoms = [
            "fever", "cough", "fatigue", "headache", "nausea", "pain",
            "shortness of breath", "dizziness", "chest pain", "abdominal pain"
        ]
        
        diseases = [
            "diabetes", "hypertension", "cancer", "asthma", "arthritis",
            "heart disease", "chronic kidney disease", "COPD", "pneumonia"
        ]
        
        treatments = [
            "medication", "physical therapy", "surgery", "chemotherapy",
            "lifestyle modification", "insulin therapy", "antibiotics"
        ]
        
        risk_factors = [
            "obesity", "smoking", "age", "family history", "sedentary lifestyle"
        ]
        
        complications = [
            "organ failure", "infection", "disability", "cardiovascular events"
        ]
        
        tests = [
            "blood glucose monitoring", "blood pressure measurement",
            "imaging studies", "laboratory tests"
        ]
        
        import random
        
        sentences = []
        for _ in range(num_sentences):
            template = random.choice(templates)
            sentence = template.format(
                symptom=random.choice(symptoms),
                disease=random.choice(diseases),
                treatment=random.choice(treatments),
                risk_factor=random.choice(risk_factors),
                complication=random.choice(complications),
                test=random.choice(tests)
            )
            sentences.append(sentence)
        
        # 保存文件
        with open(output_file, 'w', encoding='utf-8') as f:
            for sent in sentences:
                f.write(sent + '\n')
        
        print(f"✓ 已创建 {num_sentences} 条医学句子")
        print(f"✓ 保存至: {output_file}")
        
        # 验证
        self.validate_corpus_file(output_file)
    
    def split_train_test(self, input_file, train_ratio=0.8):
        """
        将数据分割为训练集和测试集
        
        Args:
            input_file: 输入文件
            train_ratio: 训练集比例
        """
        print(f"\n分割数据集...")
        
        # 读取数据
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        # 随机打乱
        import random
        random.shuffle(lines)
        
        # 分割
        split_idx = int(len(lines) * train_ratio)
        train_lines = lines[:split_idx]
        test_lines = lines[split_idx:]
        
        # 保存
        base_name = Path(input_file).stem
        train_file = f"{base_name}_train.txt"
        test_file = f"{base_name}_test.txt"
        
        with open(train_file, 'w', encoding='utf-8') as f:
            for line in train_lines:
                f.write(line + '\n')
        
        with open(test_file, 'w', encoding='utf-8') as f:
            for line in test_lines:
                f.write(line + '\n')
        
        print(f"✓ 训练集: {train_file} ({len(train_lines)} 条)")
        print(f"✓ 测试集: {test_file} ({len(test_lines)} 条)")


def main():
    """主函数"""
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          数据准备工具                                      ║
    ╚════════════════════════════════════════════════════════════╝
    
    功能:
    1. 验证语料文件
    2. 创建示例医学语料
    3. 分割训练/测试集
    """)
    
    tool = DataPreparationTool()
    
    while True:
        print("\n请选择操作:")
        print("  1. 验证现有语料文件")
        print("  2. 创建示例医学语料")
        print("  3. 分割训练/测试集")
        print("  4. 退出")
        
        choice = input("\n选择 (1/2/3/4): ").strip()
        
        if choice == '1':
            file_path = input("请输入文件路径: ").strip()
            tool.validate_corpus_file(file_path)
        
        elif choice == '2':
            output_file = input("输出文件名 (默认: medical_corpus.txt): ").strip()
            if not output_file:
                output_file = 'medical_corpus.txt'
            
            num_sentences = input("句子数量 (默认: 100): ").strip()
            num_sentences = int(num_sentences) if num_sentences else 100
            
            tool.create_sample_medical_corpus(output_file, num_sentences)
        
        elif choice == '3':
            input_file = input("请输入输入文件路径: ").strip()
            ratio = input("训练集比例 (默认: 0.8): ").strip()
            ratio = float(ratio) if ratio else 0.8
            
            tool.split_train_test(input_file, ratio)
        
        elif choice == '4':
            print("\n再见!")
            break
        
        else:
            print("无效选择，请重试")


if __name__ == '__main__':
    main()
