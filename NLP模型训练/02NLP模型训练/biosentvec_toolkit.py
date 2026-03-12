#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BioSentVec 完整工具包
支持自动GPU/CPU检测、多线程处理、模型训练和使用
"""

import os
import sys
import torch
import multiprocessing
from pathlib import Path
import numpy as np
import gensim
from gensim.models import FastText
from gensim.models.callbacks import CallbackAny2Vec
import logging
import urllib.request
from tqdm import tqdm
import zipfile
import nltk
from typing import List, Union
import warnings
warnings.filterwarnings('ignore')

# 设置日志
logging.basicConfig(
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)
logger = logging.getLogger(__name__)


class DownloadProgressBar(tqdm):
    """下载进度条"""
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)


class TrainingCallback(CallbackAny2Vec):
    """训练进度回调"""
    def __init__(self):
        self.epoch = 0
        self.losses = []

    def on_epoch_end(self, model):
        loss = model.get_latest_training_loss()
        self.losses.append(loss)
        logger.info(f'Epoch {self.epoch} 完成, Loss: {loss:.4f}')
        self.epoch += 1


class BioSentVecToolkit:
    """BioSentVec 工具包主类"""
    
    def __init__(self, base_dir='./biosentvec_data'):
        """
        初始化工具包
        
        Args:
            base_dir: 基础数据目录
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        
        # 子目录
        self.data_dir = self.base_dir / 'data'
        self.model_dir = self.base_dir / 'models'
        self.output_dir = self.base_dir / 'outputs'
        
        for dir_path in [self.data_dir, self.model_dir, self.output_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # 检测设备
        self.device = self._detect_device()
        self.num_workers = self._get_optimal_workers()
        
        # 模型
        self.word_model = None
        
        logger.info(f"设备: {self.device}")
        logger.info(f"CPU核心数: {multiprocessing.cpu_count()}")
        logger.info(f"使用线程数: {self.num_workers}")
        
        # 下载nltk数据
        self._setup_nltk()
    
    def _detect_device(self):
        """检测可用设备"""
        if torch.cuda.is_available():
            device = f"GPU (CUDA {torch.version.cuda})"
            logger.info(f"检测到GPU: {torch.cuda.get_device_name(0)}")
            return 'cuda'
        else:
            logger.info("未检测到GPU，将使用CPU")
            return 'cpu'
    
    def _get_optimal_workers(self):
        """获取最优工作线程数"""
        cpu_count = multiprocessing.cpu_count()
        if self.device == 'cpu':
            # CPU模式：使用所有核心减1
            workers = max(1, cpu_count - 1)
            logger.info(f"CPU模式：将使用 {workers} 个线程进行并行处理")
        else:
            # GPU模式：使用较少的线程
            workers = max(1, cpu_count // 2)
            logger.info(f"GPU模式：将使用 {workers} 个线程")
        return workers
    
    def _setup_nltk(self):
        """设置NLTK数据"""
        # 下载所需的NLTK数据包
        nltk_packages = ['punkt', 'punkt_tab']
        
        for package in nltk_packages:
            try:
                if package == 'punkt':
                    nltk.data.find('tokenizers/punkt')
                elif package == 'punkt_tab':
                    nltk.data.find('tokenizers/punkt_tab/english')
            except LookupError:
                logger.info(f"下载NLTK {package}数据...")
                nltk.download(package, quiet=True)
    
    def download_pretrained_model(self, model_type='word'):
        """
        下载预训练模型
        
        Args:
            model_type: 'word' 或 'sentence'
        """
        if model_type == 'word':
            url = 'https://ftp.ncbi.nlm.nih.gov/pub/lu/Suppl/BioSentVec/BioWordVec_PubMed_MIMICIII_d200.vec.bin'
            output_path = self.model_dir / 'BioWordVec_PubMed_MIMICIII_d200.vec.bin'
        else:
            url = 'https://ftp.ncbi.nlm.nih.gov/pub/lu/Suppl/BioSentVec/BioSentVec_PubMed_MIMICIII-bigram_d700.bin'
            output_path = self.model_dir / 'BioSentVec_PubMed_MIMICIII-bigram_d700.bin'
        
        if output_path.exists():
            logger.info(f"模型已存在: {output_path}")
            return output_path
        
        logger.info(f"开始下载预训练模型...")
        logger.info(f"URL: {url}")
        logger.info(f"保存至: {output_path}")
        
        with DownloadProgressBar(unit='B', unit_scale=True, miniters=1, desc='下载进度') as t:
            urllib.request.urlretrieve(url, output_path, reporthook=t.update_to)
        
        logger.info(f"下载完成!")
        return output_path
    
    def prepare_sample_data(self):
        """准备示例数据集"""
        # 训练数据
        train_texts = [
            "The patient was diagnosed with type 2 diabetes mellitus.",
            "Hypertension is a major risk factor for cardiovascular disease.",
            "The study investigated the efficacy of novel cancer treatments.",
            "Chronic kidney disease requires long-term management.",
            "Inflammatory markers were elevated in patients with rheumatoid arthritis.",
            "The mechanism of action involves inhibition of enzyme activity.",
            "Clinical trials demonstrated significant improvement in survival rates.",
            "Genetic mutations contribute to the development of various diseases.",
            "The treatment protocol includes combination therapy with multiple drugs.",
            "Biomarkers play a crucial role in early disease detection.",
            "Neurological symptoms were observed in COVID-19 patients.",
            "The pharmacokinetics of the drug were studied extensively.",
            "Metabolic syndrome increases the risk of heart disease.",
            "Immunotherapy has revolutionized cancer treatment approaches.",
            "The pathophysiology involves complex molecular mechanisms.",
        ]
        
        # 测试数据
        test_texts = [
            "Patients with diabetes require regular monitoring.",
            "High blood pressure increases cardiovascular risk.",
            "Novel therapies show promise in treating malignancies.",
        ]
        
        # 保存训练数据
        train_file = self.data_dir / 'train_data.txt'
        with open(train_file, 'w', encoding='utf-8') as f:
            for text in train_texts:
                f.write(text.lower() + '\n')
        
        # 保存测试数据
        test_file = self.data_dir / 'test_data.txt'
        with open(test_file, 'w', encoding='utf-8') as f:
            for text in test_texts:
                f.write(text.lower() + '\n')
        
        logger.info(f"示例训练数据: {train_file} ({len(train_texts)} 条)")
        logger.info(f"示例测试数据: {test_file} ({len(test_texts)} 条)")
        
        return train_file, test_file
    
    def preprocess_text(self, text):
        """
        文本预处理
        
        Args:
            text: 输入文本
            
        Returns:
            处理后的词列表
        """
        # 转小写
        text = text.lower()
        # 分词
        tokens = nltk.word_tokenize(text)
        return tokens
    
    def load_corpus(self, file_path):
        """
        加载语料库
        
        Args:
            file_path: 文件路径
            
        Yields:
            处理后的句子
        """
        logger.info(f"加载语料库: {file_path}")
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:
                    yield self.preprocess_text(line)
    
    def train_word_embeddings(
        self,
        corpus_file,
        vector_size=200,
        window=20,
        min_count=5,
        epochs=5,
        sg=0,  # 0=CBOW, 1=Skip-gram
        save_name='biowordvec_custom'
    ):
        """
        训练词向量模型
        
        Args:
            corpus_file: 语料文件路径
            vector_size: 向量维度
            window: 窗口大小
            min_count: 最小词频
            epochs: 训练轮数
            sg: 训练算法 (0=CBOW, 1=Skip-gram)
            save_name: 保存的模型名称
        """
        logger.info("=" * 50)
        logger.info("开始训练词向量模型")
        logger.info("=" * 50)
        logger.info(f"参数设置:")
        logger.info(f"  - 向量维度: {vector_size}")
        logger.info(f"  - 窗口大小: {window}")
        logger.info(f"  - 最小词频: {min_count}")
        logger.info(f"  - 训练轮数: {epochs}")
        logger.info(f"  - 算法: {'Skip-gram' if sg else 'CBOW'}")
        logger.info(f"  - 线程数: {self.num_workers}")
        
        # 加载语料
        sentences = list(self.load_corpus(corpus_file))
        logger.info(f"加载了 {len(sentences)} 条句子")
        
        # 训练回调
        callback = TrainingCallback()
        
        # 训练模型
        logger.info("开始训练...")
        self.word_model = FastText(
            sentences=sentences,
            vector_size=vector_size,
            window=window,
            min_count=min_count,
            workers=self.num_workers,
            sg=sg,
            epochs=epochs,
            callbacks=[callback]
        )
        
        # 保存模型
        model_path = self.model_dir / f'{save_name}.model'
        self.word_model.save(str(model_path))
        logger.info(f"模型已保存: {model_path}")
        
        # 保存词向量
        vector_path = self.model_dir / f'{save_name}.wv'
        self.word_model.wv.save(str(vector_path))
        logger.info(f"词向量已保存: {vector_path}")
        
        logger.info("训练完成!")
        return self.word_model
    
    def load_word_model(self, model_path=None):
        """
        加载词向量模型
        
        Args:
            model_path: 模型路径，如果为None则使用最新训练的模型
        """
        if model_path is None:
            # 查找最新的模型
            models = list(self.model_dir.glob('*.model'))
            if not models:
                logger.error("未找到模型文件!")
                return None
            model_path = max(models, key=os.path.getctime)
        
        logger.info(f"加载模型: {model_path}")
        self.word_model = FastText.load(str(model_path))
        logger.info("模型加载成功!")
        return self.word_model
    
    def get_sentence_vector(self, sentence):
        """
        获取句子向量（通过词向量平均）
        
        Args:
            sentence: 输入句子
            
        Returns:
            句子向量
        """
        if self.word_model is None:
            logger.error("请先加载或训练模型!")
            return None
        
        tokens = self.preprocess_text(sentence)
        vectors = []
        
        for token in tokens:
            if token in self.word_model.wv:
                vectors.append(self.word_model.wv[token])
        
        if not vectors:
            return np.zeros(self.word_model.vector_size)
        
        return np.mean(vectors, axis=0)
    
    def compute_similarity(self, text1, text2):
        """
        计算两个句子的相似度
        
        Args:
            text1: 句子1
            text2: 句子2
            
        Returns:
            相似度分数 (0-1)
        """
        vec1 = self.get_sentence_vector(text1)
        vec2 = self.get_sentence_vector(text2)
        
        if vec1 is None or vec2 is None:
            return None
        
        # 余弦相似度
        similarity = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
        return float(similarity)
    
    def find_similar_words(self, word, topn=10):
        """
        查找相似词
        
        Args:
            word: 目标词
            topn: 返回前N个相似词
            
        Returns:
            相似词列表
        """
        if self.word_model is None:
            logger.error("请先加载或训练模型!")
            return None
        
        try:
            similar = self.word_model.wv.most_similar(word, topn=topn)
            return similar
        except KeyError:
            logger.warning(f"词 '{word}' 不在词汇表中")
            return None
    
    def evaluate_on_test_data(self, test_file):
        """
        在测试数据上评估
        
        Args:
            test_file: 测试文件路径
        """
        logger.info("=" * 50)
        logger.info("测试数据评估")
        logger.info("=" * 50)
        
        with open(test_file, 'r', encoding='utf-8') as f:
            test_texts = [line.strip() for line in f if line.strip()]
        
        logger.info(f"测试句子数量: {len(test_texts)}")
        
        # 计算句子向量
        vectors = []
        for text in test_texts:
            vec = self.get_sentence_vector(text)
            vectors.append(vec)
        
        # 计算相似度矩阵
        logger.info("\n相似度矩阵:")
        print("\n" + " " * 10, end="")
        for i in range(len(test_texts)):
            print(f"句子{i+1:2d}", end="  ")
        print()
        
        for i, text1 in enumerate(test_texts):
            print(f"句子{i+1:2d}:", end="  ")
            for j, text2 in enumerate(test_texts):
                if i == j:
                    print("  1.00 ", end="  ")
                else:
                    sim = self.compute_similarity(text1, text2)
                    print(f"{sim:6.3f}", end="  ")
            print(f"  | {test_texts[i][:40]}...")
        
        return vectors


def main():
    """主函数 - 演示完整流程"""
    print("=" * 70)
    print("BioSentVec 工具包 - 完整演示")
    print("=" * 70)
    
    # 1. 初始化工具包
    toolkit = BioSentVecToolkit(base_dir='./biosentvec_workspace')
    
    # 2. 准备示例数据
    print("\n[步骤 1] 准备示例数据")
    train_file, test_file = toolkit.prepare_sample_data()
    
    # 3. 训练词向量模型
    print("\n[步骤 2] 训练词向量模型")
    print("提示: 在真实应用中，您需要提供更大的医学文本语料库")
    print("推荐数据源:")
    print("  - PubMed abstracts")
    print("  - MIMIC-III clinical notes")
    print("  - 医学论文和报告")
    print()
    
    user_input = input("是否训练新模型? (y/n, 默认y): ").strip().lower()
    if user_input != 'n':
        toolkit.train_word_embeddings(
            corpus_file=train_file,
            vector_size=100,  # 示例用较小维度
            window=10,
            min_count=1,
            epochs=10,
            save_name='biowordvec_demo'
        )
    else:
        toolkit.load_word_model()
    
    # 4. 测试词相似度
    print("\n[步骤 3] 测试词相似度")
    test_words = ['diabetes', 'disease', 'treatment', 'patient']
    for word in test_words:
        similar = toolkit.find_similar_words(word, topn=5)
        if similar:
            print(f"\n与 '{word}' 最相似的词:")
            for sim_word, score in similar:
                print(f"  {sim_word:20s} {score:.4f}")
    
    # 5. 测试句子相似度
    print("\n[步骤 4] 测试句子相似度")
    toolkit.evaluate_on_test_data(test_file)
    
    # 6. 自定义测试
    print("\n[步骤 5] 自定义句子相似度测试")
    test_pairs = [
        ("The patient has diabetes", "Diabetes was diagnosed in the patient"),
        ("Cancer treatment options", "Cardiovascular disease management"),
    ]
    
    for text1, text2 in test_pairs:
        sim = toolkit.compute_similarity(text1, text2)
        print(f"\n句子1: {text1}")
        print(f"句子2: {text2}")
        print(f"相似度: {sim:.4f}")
    
    print("\n" + "=" * 70)
    print("演示完成!")
    print("=" * 70)
    print("\n下一步:")
    print("1. 准备您自己的医学文本数据")
    print("2. 使用 toolkit.train_word_embeddings() 训练模型")
    print("3. 使用 toolkit.compute_similarity() 计算句子相似度")
    print("4. 或下载预训练模型: toolkit.download_pretrained_model()")


if __name__ == '__main__':
    main()