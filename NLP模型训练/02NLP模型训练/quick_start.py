#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BioSentVec 快速开始脚本
适合小白用户快速上手
"""

from biosentvec_toolkit import BioSentVecToolkit
import sys


def quick_start():
    """快速开始 - 最简单的使用方式"""
    
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          BioSentVec 快速开始工具                           ║
    ║      生物医学文本向量化和相似度计算                         ║
    ╚════════════════════════════════════════════════════════════╝
    """)
    
    # 初始化
    print("\n[1/5] 初始化工具包...")
    toolkit = BioSentVecToolkit()
    
    # 准备数据
    print("\n[2/5] 准备示例数据...")
    train_file, test_file = toolkit.prepare_sample_data()
    print(f"✓ 训练数据: {train_file}")
    print(f"✓ 测试数据: {test_file}")
    
    # 训练模型
    print("\n[3/5] 训练词向量模型...")
    print("(这可能需要几分钟，取决于数据量和硬件配置)")
    
    model = toolkit.train_word_embeddings(
        corpus_file=train_file,
        vector_size=100,
        window=10,
        epochs=10,
        save_name='quick_start_model'
    )
    
    # 测试相似词
    print("\n[4/5] 测试词语相似度...")
    print("\n查找与 'diabetes' 相似的词:")
    similar_words = toolkit.find_similar_words('diabetes', topn=5)
    if similar_words:
        for word, score in similar_words:
            print(f"  • {word:20s} 相似度: {score:.4f}")
    
    # 测试句子相似度
    print("\n[5/5] 测试句子相似度...")
    
    sentence_pairs = [
        (
            "The patient was diagnosed with type 2 diabetes mellitus.",
            "Patients with diabetes require regular monitoring."
        ),
        (
            "Hypertension is a major risk factor for cardiovascular disease.",
            "High blood pressure increases cardiovascular risk."
        ),
        (
            "The study investigated the efficacy of novel cancer treatments.",
            "Patients with diabetes require regular monitoring."
        )
    ]
    
    print("\n句子相似度计算结果:")
    print("-" * 70)
    
    for i, (sent1, sent2) in enumerate(sentence_pairs, 1):
        similarity = toolkit.compute_similarity(sent1, sent2)
        print(f"\n配对 {i}:")
        print(f"  句子A: {sent1}")
        print(f"  句子B: {sent2}")
        print(f"  相似度: {similarity:.4f} {'(高度相似)' if similarity > 0.7 else '(相关)' if similarity > 0.5 else '(相关性较低)'}")
    
    print("\n" + "=" * 70)
    print("✓ 快速开始完成!")
    print("=" * 70)
    
    print("\n您可以使用以下代码计算自己的句子:")
    print("""
    from biosentvec_toolkit import BioSentVecToolkit
    
    # 加载工具包
    toolkit = BioSentVecToolkit()
    toolkit.load_word_model()  # 加载已训练的模型
    
    # 计算相似度
    similarity = toolkit.compute_similarity(
        "Your first sentence here",
        "Your second sentence here"
    )
    print(f"相似度: {similarity:.4f}")
    """)


def interactive_mode():
    """交互模式"""
    
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          BioSentVec 交互模式                               ║
    ╚════════════════════════════════════════════════════════════╝
    """)
    
    toolkit = BioSentVecToolkit()
    
    # 检查是否有已训练的模型
    model_loaded = False
    try:
        toolkit.load_word_model()
        model_loaded = True
        print("✓ 已加载最新训练的模型")
    except:
        print("! 未找到已训练的模型，请先训练模型")
    
    if not model_loaded:
        print("\n是否使用示例数据训练一个新模型?")
        choice = input("输入 'y' 继续，'n' 退出: ").strip().lower()
        
        if choice == 'y':
            train_file, _ = toolkit.prepare_sample_data()
            toolkit.train_word_embeddings(
                corpus_file=train_file,
                vector_size=100,
                epochs=10,
                save_name='interactive_model'
            )
        else:
            print("退出程序")
            return
    
    print("\n" + "=" * 70)
    print("现在您可以输入两个句子来计算它们的相似度")
    print("(输入 'q' 退出)")
    print("=" * 70)
    
    while True:
        print("\n")
        sent1 = input("请输入第一个句子: ").strip()
        if sent1.lower() == 'q':
            break
        
        sent2 = input("请输入第二个句子: ").strip()
        if sent2.lower() == 'q':
            break
        
        if sent1 and sent2:
            similarity = toolkit.compute_similarity(sent1, sent2)
            print(f"\n相似度: {similarity:.4f}")
            
            if similarity > 0.8:
                print("评价: 高度相似 ⭐⭐⭐⭐⭐")
            elif similarity > 0.6:
                print("评价: 较为相似 ⭐⭐⭐⭐")
            elif similarity > 0.4:
                print("评价: 中等相似 ⭐⭐⭐")
            elif similarity > 0.2:
                print("评价: 稍微相关 ⭐⭐")
            else:
                print("评价: 相关性较低 ⭐")
        else:
            print("请输入有效的句子!")
    
    print("\n谢谢使用!")


def show_menu():
    """显示菜单"""
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          BioSentVec 工具包 - 主菜单                        ║
    ╚════════════════════════════════════════════════════════════╝
    
    请选择操作模式:
    
    1. 快速开始 (推荐新手)
       - 自动准备示例数据
       - 训练模型
       - 运行测试
    
    2. 交互模式
       - 输入自定义句子
       - 实时计算相似度
    
    3. 退出
    
    """)
    
    choice = input("请输入选项 (1/2/3): ").strip()
    return choice


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'quick':
            quick_start()
        elif sys.argv[1] == 'interactive':
            interactive_mode()
        else:
            print("未知参数，请使用: python quick_start.py [quick|interactive]")
    else:
        while True:
            choice = show_menu()
            
            if choice == '1':
                quick_start()
                break
            elif choice == '2':
                interactive_mode()
                break
            elif choice == '3':
                print("\n再见!")
                break
            else:
                print("\n无效选项，请重新选择\n")
