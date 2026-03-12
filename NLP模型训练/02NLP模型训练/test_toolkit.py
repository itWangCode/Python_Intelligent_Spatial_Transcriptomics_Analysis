#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BioSentVec 工具包测试脚本
验证所有功能是否正常工作
"""

import sys
import traceback
from pathlib import Path


def test_imports():
    """测试1: 检查所有依赖包"""
    print("\n" + "=" * 70)
    print("测试 1: 检查依赖包")
    print("=" * 70)
    
    packages = {
        'torch': 'PyTorch',
        'numpy': 'NumPy',
        'gensim': 'Gensim',
        'nltk': 'NLTK',
        'tqdm': 'tqdm'
    }
    
    all_ok = True
    for pkg, name in packages.items():
        try:
            __import__(pkg)
            print(f"✓ {name:15s} - 已安装")
        except ImportError:
            print(f"✗ {name:15s} - 未安装")
            all_ok = False
    
    return all_ok


def test_device_detection():
    """测试2: 设备检测"""
    print("\n" + "=" * 70)
    print("测试 2: 设备检测")
    print("=" * 70)
    
    try:
        import torch
        import multiprocessing
        
        print(f"CPU核心数: {multiprocessing.cpu_count()}")
        
        if torch.cuda.is_available():
            print(f"✓ GPU可用: {torch.cuda.get_device_name(0)}")
            print(f"  CUDA版本: {torch.version.cuda}")
        else:
            print("⚠ GPU不可用，将使用CPU")
        
        return True
    except Exception as e:
        print(f"✗ 错误: {e}")
        return False


def test_toolkit_initialization():
    """测试3: 工具包初始化"""
    print("\n" + "=" * 70)
    print("测试 3: 工具包初始化")
    print("=" * 70)
    
    try:
        from biosentvec_toolkit import BioSentVecToolkit
        
        toolkit = BioSentVecToolkit(base_dir='./test_workspace')
        print("✓ 工具包初始化成功")
        print(f"  基础目录: {toolkit.base_dir}")
        print(f"  设备: {toolkit.device}")
        print(f"  工作线程: {toolkit.num_workers}")
        
        return True, toolkit
    except Exception as e:
        print(f"✗ 错误: {e}")
        traceback.print_exc()
        return False, None


def test_data_preparation(toolkit):
    """测试4: 数据准备"""
    print("\n" + "=" * 70)
    print("测试 4: 数据准备")
    print("=" * 70)
    
    try:
        train_file, test_file = toolkit.prepare_sample_data()
        
        # 检查文件
        if Path(train_file).exists() and Path(test_file).exists():
            print(f"✓ 训练文件: {train_file}")
            print(f"✓ 测试文件: {test_file}")
            
            # 读取并显示行数
            with open(train_file, 'r') as f:
                train_lines = len(f.readlines())
            with open(test_file, 'r') as f:
                test_lines = len(f.readlines())
            
            print(f"  训练句子数: {train_lines}")
            print(f"  测试句子数: {test_lines}")
            
            return True, train_file, test_file
        else:
            print("✗ 文件创建失败")
            return False, None, None
            
    except Exception as e:
        print(f"✗ 错误: {e}")
        traceback.print_exc()
        return False, None, None


def test_model_training(toolkit, train_file):
    """测试5: 模型训练"""
    print("\n" + "=" * 70)
    print("测试 5: 模型训练")
    print("=" * 70)
    print("(这可能需要几分钟...)")
    
    try:
        model = toolkit.train_word_embeddings(
            corpus_file=train_file,
            vector_size=50,  # 小维度加快测试
            window=5,
            min_count=1,
            epochs=5,
            save_name='test_model'
        )
        
        if model is not None:
            print("✓ 模型训练成功")
            print(f"  词汇表大小: {len(model.wv)}")
            print(f"  向量维度: {model.vector_size}")
            return True
        else:
            print("✗ 模型训练失败")
            return False
            
    except Exception as e:
        print(f"✗ 错误: {e}")
        traceback.print_exc()
        return False


def test_word_similarity(toolkit):
    """测试6: 词语相似度"""
    print("\n" + "=" * 70)
    print("测试 6: 词语相似度")
    print("=" * 70)
    
    try:
        test_words = ['diabetes', 'patient', 'disease']
        
        for word in test_words:
            similar = toolkit.find_similar_words(word, topn=3)
            if similar:
                print(f"\n'{word}' 的相似词:")
                for sim_word, score in similar:
                    print(f"  {sim_word}: {score:.4f}")
            else:
                print(f"⚠ 词 '{word}' 不在词汇表中")
        
        return True
        
    except Exception as e:
        print(f"✗ 错误: {e}")
        traceback.print_exc()
        return False


def test_sentence_similarity(toolkit):
    """测试7: 句子相似度"""
    print("\n" + "=" * 70)
    print("测试 7: 句子相似度")
    print("=" * 70)
    
    try:
        sent_pairs = [
            (
                "The patient has diabetes",
                "Diabetes was diagnosed in patient"
            ),
            (
                "Cancer treatment options",
                "The patient has diabetes"
            )
        ]
        
        for sent1, sent2 in sent_pairs:
            sim = toolkit.compute_similarity(sent1, sent2)
            print(f"\n句子1: {sent1}")
            print(f"句子2: {sent2}")
            print(f"相似度: {sim:.4f}")
        
        return True
        
    except Exception as e:
        print(f"✗ 错误: {e}")
        traceback.print_exc()
        return False


def test_model_save_load(toolkit):
    """测试8: 模型保存和加载"""
    print("\n" + "=" * 70)
    print("测试 8: 模型保存和加载")
    print("=" * 70)
    
    try:
        # 模型应该已经在训练时保存了
        model_files = list(toolkit.model_dir.glob('*.model'))
        
        if model_files:
            print(f"✓ 找到 {len(model_files)} 个模型文件")
            
            # 尝试加载
            loaded_model = toolkit.load_word_model(model_files[0])
            if loaded_model:
                print(f"✓ 模型加载成功: {model_files[0].name}")
                return True
            else:
                print("✗ 模型加载失败")
                return False
        else:
            print("⚠ 未找到模型文件")
            return False
            
    except Exception as e:
        print(f"✗ 错误: {e}")
        traceback.print_exc()
        return False


def cleanup_test_files():
    """清理测试文件"""
    print("\n" + "=" * 70)
    print("清理测试文件")
    print("=" * 70)
    
    import shutil
    
    test_dir = Path('./test_workspace')
    if test_dir.exists():
        choice = input(f"\n是否删除测试目录 {test_dir}? (y/n): ").strip().lower()
        if choice == 'y':
            shutil.rmtree(test_dir)
            print("✓ 测试目录已删除")
        else:
            print("⚠ 保留测试目录")


def run_all_tests():
    """运行所有测试"""
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          BioSentVec 工具包测试套件                         ║
    ╚════════════════════════════════════════════════════════════╝
    
    将运行以下测试:
    1. 检查依赖包
    2. 设备检测
    3. 工具包初始化
    4. 数据准备
    5. 模型训练
    6. 词语相似度
    7. 句子相似度
    8. 模型保存和加载
    """)
    
    input("按Enter开始测试...")
    
    results = []
    
    # 测试1: 依赖包
    results.append(('依赖包检查', test_imports()))
    
    if not results[-1][1]:
        print("\n✗ 依赖包检查失败，请先运行: pip install -r requirements.txt")
        return
    
    # 测试2: 设备检测
    results.append(('设备检测', test_device_detection()))
    
    # 测试3: 初始化
    success, toolkit = test_toolkit_initialization()
    results.append(('工具包初始化', success))
    
    if not success:
        print("\n✗ 工具包初始化失败")
        return
    
    # 测试4: 数据准备
    success, train_file, test_file = test_data_preparation(toolkit)
    results.append(('数据准备', success))
    
    if not success:
        print("\n✗ 数据准备失败")
        return
    
    # 测试5: 模型训练
    results.append(('模型训练', test_model_training(toolkit, train_file)))
    
    if not results[-1][1]:
        print("\n✗ 模型训练失败")
        return
    
    # 测试6: 词语相似度
    results.append(('词语相似度', test_word_similarity(toolkit)))
    
    # 测试7: 句子相似度
    results.append(('句子相似度', test_sentence_similarity(toolkit)))
    
    # 测试8: 保存加载
    results.append(('模型保存加载', test_model_save_load(toolkit)))
    
    # 总结
    print("\n" + "=" * 70)
    print("测试总结")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✓ 通过" if result else "✗ 失败"
        print(f"{test_name:20s} {status}")
    
    print(f"\n总计: {passed}/{total} 测试通过")
    
    if passed == total:
        print("\n🎉 所有测试通过！工具包已准备就绪！")
    else:
        print(f"\n⚠ {total - passed} 个测试失败，请检查错误信息")
    
    # 清理
    cleanup_test_files()


if __name__ == '__main__':
    try:
        run_all_tests()
    except KeyboardInterrupt:
        print("\n\n测试被用户中断")
    except Exception as e:
        print(f"\n\n测试过程中发生错误: {e}")
        traceback.print_exc()
