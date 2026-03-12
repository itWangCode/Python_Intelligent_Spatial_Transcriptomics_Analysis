#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
批量爬取多个疾病关键词
"""

from pubmed_crawler import PubMedCrawler
import time


def batch_crawl():
    """批量爬取"""
    
    # 初始化爬虫
    crawler = PubMedCrawler(output_dir='./pubmed_batch_data')
    
    # 定义要爬取的疾病关键词
    keywords = [
        "diabetes",
        "hypertension", 
        "cancer",
        "asthma",
        "arthritis",
        "alzheimer",
        "heart disease",
        "stroke",
        "COPD",
        "kidney disease"
    ]
    
    # 每个关键词爬取数量
    max_results_per_keyword = 1000  # 可改为 2000 或 0（不限）
    
    print("=" * 70)
    print(f"批量爬取 {len(keywords)} 个疾病关键词")
    print(f"每个关键词: {max_results_per_keyword} 篇")
    print("=" * 70)
    
    all_files = []
    
    for i, keyword in enumerate(keywords, 1):
        print(f"\n[{i}/{len(keywords)}] 处理关键词: {keyword}")
        
        try:
            output_path = crawler.crawl(
                keyword=keyword,
                max_results=max_results_per_keyword,
                output_filename=f"{keyword.replace(' ', '_')}.txt"
            )
            
            if output_path:
                all_files.append(output_path)
            
            # 间隔，避免请求过快
            if i < len(keywords):
                print(f"\n等待3秒后继续...")
                time.sleep(3)
                
        except Exception as e:
            print(f"✗ 错误: {e}")
            continue
    
    # 合并所有文件
    if all_files:
        print("\n" + "=" * 70)
        print("合并所有文件...")
        
        merged_file = crawler.output_dir / 'all_diseases_merged.txt'
        
        with open(merged_file, 'w', encoding='utf-8') as outfile:
            for file_path in all_files:
                with open(file_path, 'r', encoding='utf-8') as infile:
                    outfile.write(infile.read())
        
        print(f"✓ 合并完成: {merged_file}")
        
        # 统计
        with open(merged_file, 'r', encoding='utf-8') as f:
            total_lines = len(f.readlines())
        
        print(f"✓ 总句子数: {total_lines}")
        print("\n" + "=" * 70)
        print("🎉 批量爬取完成！")
        print("=" * 70)


if __name__ == '__main__':
    batch_crawl()
