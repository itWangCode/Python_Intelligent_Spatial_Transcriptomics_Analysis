#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PubMed简单爬虫 - 单文件版本
直接运行即可
"""

import requests
from bs4 import BeautifulSoup
import time
import re
from pathlib import Path
from tqdm import tqdm


# ============ 配置区 ============
KEYWORD = "diabetes"              # 搜索关键词
MAX_RESULTS = 1000                # 1000, 2000, 或 0（不限）
OUTPUT_FILE = "pubmed_data.txt"   # 输出文件名
# ================================


def request_with_retry(url, params=None, max_retries=5):
    """带重试的请求"""
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, params=params, headers=headers, timeout=30)
            response.raise_for_status()
            return response
        except:
            if attempt < max_retries - 1:
                wait = 3 * (attempt + 1)
                print(f"⚠ 重试中... ({attempt + 1}/{max_retries})，等待{wait}秒")
                time.sleep(wait)
            else:
                raise
    return None


def search_pubmed(keyword, max_results):
    """搜索PubMed获取ID"""
    print(f"🔍 搜索: {keyword}")
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    if max_results == 0:
        params = {'db': 'pubmed', 'term': keyword, 'retmode': 'json', 'retmax': 1}
        r = request_with_retry(url, params)
        max_results = int(r.json()['esearchresult']['count'])
        print(f"✓ 找到 {max_results} 篇")
    
    all_ids = []
    batch = 10000
    
    for start in range(0, max_results, batch):
        params = {
            'db': 'pubmed',
            'term': keyword,
            'retstart': start,
            'retmax': min(batch, max_results - start),
            'retmode': 'json'
        }
        r = request_with_retry(url, params)
        ids = r.json()['esearchresult']['idlist']
        all_ids.extend(ids)
        print(f"  已获取 {len(all_ids)}/{max_results} ID")
        time.sleep(0.5)
    
    return all_ids


def fetch_abstracts(id_list):
    """获取摘要"""
    print(f"📥 下载 {len(id_list)} 篇摘要...")
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    abstracts = []
    batch = 200
    
    with tqdm(total=len(id_list)) as pbar:
        for i in range(0, len(id_list), batch):
            batch_ids = id_list[i:i + batch]
            params = {
                'db': 'pubmed',
                'id': ','.join(batch_ids),
                'rettype': 'abstract',
                'retmode': 'xml'
            }
            r = request_with_retry(url, params)
            soup = BeautifulSoup(r.content, 'xml')
            
            for article in soup.find_all('PubmedArticle'):
                abstract = article.find('AbstractText')
                if abstract:
                    abstracts.append(abstract.get_text())
            
            pbar.update(len(batch_ids))
            time.sleep(0.5)
    
    return abstracts


def split_sentences(text):
    """分句"""
    text = re.sub(r'\s+', ' ', text).strip()
    sentences = re.split(r'[.!?]+', text)
    
    result = []
    for s in sentences:
        s = s.strip()
        words = len(s.split())
        if 5 <= words <= 300:
            result.append(s)
    
    return result


def main():
    print("=" * 70)
    print("PubMed 爬虫")
    print("=" * 70)
    
    # 1. 搜索
    ids = search_pubmed(KEYWORD, MAX_RESULTS)
    
    # 2. 下载
    abstracts = fetch_abstracts(ids)
    print(f"✓ 下载了 {len(abstracts)} 篇")
    
    # 3. 分句
    print("📝 分句中...")
    all_sentences = []
    for abstract in tqdm(abstracts):
        all_sentences.extend(split_sentences(abstract))
    
    # 4. 保存
    Path(OUTPUT_FILE).parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        for sent in all_sentences:
            f.write(sent + '\n')
    
    print(f"\n✓ 完成!")
    print(f"✓ 文件: {OUTPUT_FILE}")
    print(f"✓ 句子数: {len(all_sentences)}")
    print("=" * 70)


if __name__ == '__main__':
    main()
