#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
交互式模型测试工具
实时输入句子，查看相似度
"""

from biosentvec_toolkit import BioSentVecToolkit

def main():
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          BioSentVec 交互式测试工具                         ║
    ╚════════════════════════════════════════════════════════════╝
    """)
    
    # 加载模型
    print("正在加载模型...")
    toolkit = BioSentVecToolkit()
    model = toolkit.load_word_model()
    
    if model is None:
        print("错误: 未找到模型！请先训练模型。")
        return
    
    print(f"✓ 模型加载成功")
    print(f"✓ 词汇表大小: {len(model.wv)}")
    print(f"✓ 向量维度: {model.vector_size}")
    
    print("\n" + "=" * 70)
    print("功能菜单:")
    print("  1. 计算句子相似度")
    print("  2. 查找相似词")
    print("  3. 批量测试")
    print("  4. 退出")
    print("=" * 70)
    
    while True:
        print("\n")
        choice = input("请选择功能 (1/2/3/4): ").strip()
        
        if choice == '1':
            # 句子相似度
            print("\n--- 句子相似度计算 ---")
            sent1 = input("句子1: ").strip()
            sent2 = input("句子2: ").strip()
            
            if sent1 and sent2:
                sim = toolkit.compute_similarity(sent1, sent2)
                print(f"\n相似度: {sim:.4f}")
                
                # 可视化
                bar_length = int(sim * 50)
                bar = "█" * bar_length + "░" * (50 - bar_length)
                print(f"图示: [{bar}] {sim*100:.1f}%")
                
                if sim > 0.8:
                    print("评价: ⭐⭐⭐⭐⭐ 高度相似")
                elif sim > 0.6:
                    print("评价: ⭐⭐⭐⭐ 较为相似")
                elif sim > 0.4:
                    print("评价: ⭐⭐⭐ 中等相似")
                elif sim > 0.2:
                    print("评价: ⭐⭐ 稍微相关")
                else:
                    print("评价: ⭐ 相关性低")
        
        elif choice == '2':
            # 相似词查询
            print("\n--- 相似词查询 ---")
            word = input("输入词语: ").strip().lower()
            topn = input("返回数量 (默认10): ").strip()
            topn = int(topn) if topn else 10
            
            similar = toolkit.find_similar_words(word, topn=topn)
            
            if similar:
                print(f"\n与 '{word}' 最相似的 {topn} 个词:")
                print("-" * 50)
                for i, (sim_word, score) in enumerate(similar, 1):
                    bar_length = int(score * 30)
                    bar = "█" * bar_length
                    print(f"{i:2d}. {sim_word:20s} {score:.4f} {bar}")
            else:
                print(f"词 '{word}' 不在词汇表中")
        
        elif choice == '3':
            # 批量测试
            print("\n--- 批量测试 ---")
            print("输入句子对，每行格式: 句子1 | 句子2")
            print("输入空行结束")
            
            pairs = []
            while True:
                line = input("> ").strip()
                if not line:
                    break
                
                if '|' in line:
                    sent1, sent2 = line.split('|', 1)
                    pairs.append((sent1.strip(), sent2.strip()))
            
            if pairs:
                print(f"\n测试 {len(pairs)} 个句子对:")
                print("-" * 70)
                
                for i, (s1, s2) in enumerate(pairs, 1):
                    sim = toolkit.compute_similarity(s1, s2)
                    print(f"{i}. [{sim:.4f}] {s1} <-> {s2}")
        
        elif choice == '4':
            print("\n再见!")
            break
        
        else:
            print("无效选项，请重新选择")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n程序被中断")
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
