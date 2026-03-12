#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PubMed 医学文献爬虫
自动下载医学摘要用于训练BioSentVec模型
"""

import requests
from bs4 import BeautifulSoup
import time
import re
from pathlib import Path
import random
from tqdm import tqdm
import json


class PubMedCrawler:
    """PubMed爬虫"""
    
    def __init__(self, output_dir='./pubmed_data'):
        self.base_url = "https://pubmed.ncbi.nlm.nih.gov"
        self.api_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 请求头，模拟浏览器
        self.headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
            'Connection': 'keep-alive',
        }
        
        # 会话保持连接
        self.session = requests.Session()
        self.session.headers.update(self.headers)
        
        # 重试配置
        self.max_retries = 5
        self.retry_delay = 3
    
    def request_with_retry(self, url, params=None):
        """带重试的请求，自动重连"""
        for attempt in range(self.max_retries):
            try:
                response = self.session.get(
                    url, 
                    params=params, 
                    timeout=30,
                    verify=True
                )
                response.raise_for_status()
                return response
            
            except requests.exceptions.ConnectionError:
                if attempt < self.max_retries - 1:
                    wait_time = self.retry_delay * (attempt + 1)
                    print(f"⚠ 网络断开，{wait_time}秒后重试... ({attempt + 1}/{self.max_retries})")
                    time.sleep(wait_time)
                else:
                    print(f"✗ 网络连接失败，已重试{self.max_retries}次")
                    raise
            
            except requests.exceptions.Timeout:
                if attempt < self.max_retries - 1:
                    print(f"⚠ 请求超时，重试中... ({attempt + 1}/{self.max_retries})")
                    time.sleep(self.retry_delay)
                else:
                    raise
            
            except Exception as e:
                if attempt < self.max_retries - 1:
                    print(f"⚠ 请求错误: {e}，重试中...")
                    time.sleep(self.retry_delay)
                else:
                    raise
        
        return None
    
    def search_pubmed(self, keyword, max_results=1000):
        """
        搜索PubMed
        
        Args:
            keyword: 搜索关键词（如 "diabetes"）
            max_results: 最大结果数（1000, 2000, 或 0表示不限）
        
        Returns:
            文章ID列表
        """
        print(f"\n🔍 搜索关键词: {keyword}")
        print(f"📊 目标数量: {'不限' if max_results == 0 else max_results}")
        
        # 第一步：搜索获取ID列表
        search_url = f"{self.api_url}/esearch.fcgi"
        
        if max_results == 0:
            # 不限制，先获取总数
            params = {
                'db': 'pubmed',
                'term': keyword,
                'retmode': 'json',
                'retmax': 1
            }
            response = self.request_with_retry(search_url, params)
            data = response.json()
            total_count = int(data['esearchresult']['count'])
            print(f"✓ 找到 {total_count} 篇文章")
            max_results = total_count
        
        # 分批获取ID（每次最多10000）
        all_ids = []
        batch_size = 10000
        
        for start in range(0, max_results, batch_size):
            retmax = min(batch_size, max_results - start)
            
            params = {
                'db': 'pubmed',
                'term': keyword,
                'retstart': start,
                'retmax': retmax,
                'retmode': 'json'
            }
            
            response = self.request_with_retry(search_url, params)
            if response:
                data = response.json()
                ids = data['esearchresult']['idlist']
                all_ids.extend(ids)
                print(f"  已获取 {len(all_ids)}/{max_results} 个ID")
            
            time.sleep(0.5)  # 避免请求过快
        
        return all_ids
    
    def fetch_abstracts(self, id_list):
        """
        获取摘要内容
        
        Args:
            id_list: PubMed ID列表
        
        Returns:
            摘要文本列表
        """
        print(f"\n📥 开始下载 {len(id_list)} 篇摘要...")
        
        abstracts = []
        fetch_url = f"{self.api_url}/efetch.fcgi"
        
        # 分批获取（每次200篇）
        batch_size = 200
        
        with tqdm(total=len(id_list), desc="下载进度") as pbar:
            for i in range(0, len(id_list), batch_size):
                batch_ids = id_list[i:i + batch_size]
                
                params = {
                    'db': 'pubmed',
                    'id': ','.join(batch_ids),
                    'rettype': 'abstract',
                    'retmode': 'xml'
                }
                
                response = self.request_with_retry(fetch_url, params)
                if response:
                    # 解析XML
                    soup = BeautifulSoup(response.content, 'xml')
                    articles = soup.find_all('PubmedArticle')
                    
                    for article in articles:
                        abstract_elem = article.find('AbstractText')
                        if abstract_elem:
                            text = abstract_elem.get_text()
                            abstracts.append(text)
                
                pbar.update(len(batch_ids))
                time.sleep(0.5)  # 遵守NCBI请求限制
        
        return abstracts
    
    def clean_and_split_sentences(self, text):
        """
        清理文本并分句
        
        Args:
            text: 摘要文本
        
        Returns:
            句子列表
        """
        # 移除特殊字符
        text = re.sub(r'\s+', ' ', text)
        text = text.strip()
        
        # 简单分句（按句号、问号、感叹号分）
        sentences = re.split(r'[.!?]+', text)
        
        # 清理并过滤
        clean_sentences = []
        for sent in sentences:
            sent = sent.strip()
            # 只保留长度合适的句子（5-300词）
            word_count = len(sent.split())
            if 5 <= word_count <= 300:
                clean_sentences.append(sent)
        
        return clean_sentences
    
    def save_to_file(self, sentences, filename):
        """保存到文件"""
        output_path = self.output_dir / filename
        
        with open(output_path, 'w', encoding='utf-8') as f:
            for sent in sentences:
                f.write(sent + '\n')
        
        return output_path
    
    def crawl(self, keyword, max_results=1000, output_filename=None):
        """
        完整爬取流程
        
        Args:
            keyword: 搜索关键词
            max_results: 最大结果数（0表示不限）
            output_filename: 输出文件名
        
        Returns:
            输出文件路径
        """
        print("=" * 70)
        print("PubMed 医学文献爬虫")
        print("=" * 70)
        
        # 1. 搜索文章
        id_list = self.search_pubmed(keyword, max_results)
        
        if not id_list:
            print("✗ 未找到文章")
            return None
        
        # 2. 获取摘要
        abstracts = self.fetch_abstracts(id_list)
        print(f"✓ 成功下载 {len(abstracts)} 篇摘要")
        
        # 3. 分句和清理
        print(f"\n📝 处理文本...")
        all_sentences = []
        for abstract in tqdm(abstracts, desc="分句处理"):
            sentences = self.clean_and_split_sentences(abstract)
            all_sentences.extend(sentences)
        
        print(f"✓ 提取了 {len(all_sentences)} 个句子")
        
        # 4. 保存文件
        if output_filename is None:
            output_filename = f"pubmed_{keyword.replace(' ', '_')}.txt"
        
        output_path = self.save_to_file(all_sentences, output_filename)
        
        print(f"\n✓ 保存成功: {output_path}")
        print(f"✓ 总句子数: {len(all_sentences)}")
        
        # 保存统计信息
        stats = {
            'keyword': keyword,
            'articles': len(abstracts),
            'sentences': len(all_sentences),
            'output_file': str(output_path)
        }
        
        stats_file = self.output_dir / f"{output_filename}.stats.json"
        with open(stats_file, 'w', encoding='utf-8') as f:
            json.dump(stats, f, indent=2)
        
        return output_path


def main():
    """主函数"""
    print("""
    ╔════════════════════════════════════════════════════════════╗
    ║          PubMed 医学文献爬虫                               ║
    ║      自动下载医学摘要用于BioSentVec训练                    ║
    ╚════════════════════════════════════════════════════════════╝
    """)
    
    # 初始化爬虫
    crawler = PubMedCrawler()
    
    # 输入关键词
    keyword = input("请输入搜索关键词（如 diabetes）: ").strip()
    if not keyword:
        keyword = "diabetes"
        print(f"使用默认关键词: {keyword}")
    
    # 输入数量限制
    print("\n选择下载数量:")
    print("  1. 1000篇")
    print("  2. 2000篇")
    print("  3. 不限")
    
    choice = input("请选择 (1/2/3): ").strip()
    
    if choice == '1':
        max_results = 1000
    elif choice == '2':
        max_results = 2000
    elif choice == '3':
        max_results = 0  # 不限
    else:
        max_results = 1000
        print("默认选择: 1000篇")
    
    # 开始爬取
    try:
        output_path = crawler.crawl(keyword, max_results)
        
        if output_path:
            print("\n" + "=" * 70)
            print("🎉 爬取完成！")
            print("=" * 70)
            print(f"\n数据已保存至: {output_path}")
            print(f"\n下一步:")
            print(f"  1. 使用此文件训练模型:")
            print(f"     toolkit.train_word_embeddings(corpus_file='{output_path}')")
            
    except KeyboardInterrupt:
        print("\n\n⚠ 用户中断")
    except Exception as e:
        print(f"\n✗ 错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
