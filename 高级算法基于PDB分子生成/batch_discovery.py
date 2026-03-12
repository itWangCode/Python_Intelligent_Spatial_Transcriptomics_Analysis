#!/usr/bin/env python3
"""
批量自动化药物发现系统
自动读取PDB文件夹，对每个靶点进行完整的工业级筛选

特性:
- 自动扫描PDB文件夹
- 并行处理多个靶点
- 完整的ADMET + 类药性 + 对接流程
- 统一的结果汇总和排名
- 进度追踪和断点续传
"""

import os
import sys
import glob
import json
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime
import pandas as pd
from tqdm import tqdm
import argparse

# 导入工业级流程
from industrial_pipeline import IndustrialDrugDiscoveryPipeline, CompleteMoleculeProfile

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('batch_processing.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class BatchDrugDiscovery:
    """批量药物发现系统"""
    
    def __init__(self, 
                 pdb_folder: str,
                 output_base: str = 'batch_results',
                 device: str = None,
                 resume: bool = True):
        """
        初始化批量处理系统
        
        Args:
            pdb_folder: PDB文件夹路径
            output_base: 输出基础目录
            device: 计算设备
            resume: 是否断点续传
        """
        self.pdb_folder = Path(pdb_folder)
        self.output_base = Path(output_base)
        self.device = device
        self.resume = resume
        
        # 创建输出目录
        self.output_base.mkdir(parents=True, exist_ok=True)
        
        # 进度文件
        self.progress_file = self.output_base / 'progress.json'
        self.progress = self._load_progress()
        
        # 初始化流程
        self.pipeline = IndustrialDrugDiscoveryPipeline(device=device)
        
        logger.info(f"批量处理系统初始化完成")
        logger.info(f"PDB文件夹: {self.pdb_folder}")
        logger.info(f"输出目录: {self.output_base}")
    
    def scan_pdb_files(self) -> List[Path]:
        """扫描PDB文件"""
        # 支持多种扩展名
        pdb_files = []
        for ext in ['*.pdb', '*.PDB', '*.ent']:
            pdb_files.extend(self.pdb_folder.glob(ext))
        
        pdb_files = sorted(set(pdb_files))
        logger.info(f"发现 {len(pdb_files)} 个PDB文件")
        
        return pdb_files
    
    def run_batch_processing(self,
                           candidates_per_target: int = 100,
                           num_atoms_range: tuple = (15, 35),
                           sampling_steps: int = 500,
                           enable_docking: bool = True,
                           max_targets: Optional[int] = None) -> Dict:
        """
        运行批量处理
        
        Args:
            candidates_per_target: 每个靶点生成的候选数
            num_atoms_range: 原子数范围
            sampling_steps: 采样步数
            enable_docking: 是否启用对接
            max_targets: 最大处理靶点数(None=全部)
        """
        # 扫描PDB文件
        pdb_files = self.scan_pdb_files()
        
        if max_targets:
            pdb_files = pdb_files[:max_targets]
        
        logger.info("=" * 80)
        logger.info(f"批量处理 {len(pdb_files)} 个靶点")
        logger.info("=" * 80)
        
        # 统计信息
        stats = {
            'total_targets': len(pdb_files),
            'completed': 0,
            'failed': 0,
            'total_molecules': 0,
            'total_recommended': 0,
            'start_time': datetime.now().isoformat(),
            'results': []
        }
        
        # 处理每个PDB
        for idx, pdb_file in enumerate(pdb_files, 1):
            target_name = pdb_file.stem
            
            # 检查是否已处理
            if self.resume and target_name in self.progress.get('completed', []):
                logger.info(f"[{idx}/{len(pdb_files)}] 跳过已完成: {target_name}")
                stats['completed'] += 1
                continue
            
            logger.info("\n" + "=" * 80)
            logger.info(f"[{idx}/{len(pdb_files)}] 处理靶点: {target_name}")
            logger.info("=" * 80)
            
            try:
                # 处理单个靶点
                result = self._process_single_target(
                    pdb_file=pdb_file,
                    target_name=target_name,
                    candidates=candidates_per_target,
                    num_atoms_range=num_atoms_range,
                    sampling_steps=sampling_steps,
                    enable_docking=enable_docking
                )
                
                # 更新统计
                stats['completed'] += 1
                stats['total_molecules'] += result['total_molecules']
                stats['total_recommended'] += result['recommended_count']
                stats['results'].append(result)
                
                # 保存进度
                self._save_progress(target_name, 'completed')
                
                logger.info(f"✓ {target_name} 完成: {result['recommended_count']} 个推荐候选")
                
            except Exception as e:
                logger.error(f"✗ {target_name} 失败: {e}")
                stats['failed'] += 1
                self._save_progress(target_name, 'failed')
                continue
        
        # 完成处理
        stats['end_time'] = datetime.now().isoformat()
        stats['success_rate'] = stats['completed'] / stats['total_targets'] * 100
        
        # 汇总结果
        self._generate_summary(stats)
        
        logger.info("\n" + "=" * 80)
        logger.info("批量处理完成!")
        logger.info("=" * 80)
        logger.info(f"总靶点数: {stats['total_targets']}")
        logger.info(f"成功: {stats['completed']}")
        logger.info(f"失败: {stats['failed']}")
        logger.info(f"成功率: {stats['success_rate']:.1f}%")
        logger.info(f"总分子数: {stats['total_molecules']}")
        logger.info(f"推荐候选: {stats['total_recommended']}")
        logger.info(f"结果保存在: {self.output_base}")
        
        return stats
    
    def _process_single_target(self,
                              pdb_file: Path,
                              target_name: str,
                              candidates: int,
                              num_atoms_range: tuple,
                              sampling_steps: int,
                              enable_docking: bool) -> Dict:
        """处理单个靶点"""
        # 输出目录
        target_output = self.output_base / target_name
        target_output.mkdir(parents=True, exist_ok=True)
        
        # 运行工业级流程
        start_time = time.time()
        
        profiles = self.pipeline.run_complete_pipeline(
            pdb_path=str(pdb_file),
            num_candidates=candidates,
            num_atoms_range=num_atoms_range,
            sampling_steps=sampling_steps,
            output_dir=str(target_output),
            enable_docking=enable_docking
        )
        
        elapsed_time = time.time() - start_time
        
        # 统计结果
        recommended = [p for p in profiles if p.recommended]
        
        result = {
            'target_name': target_name,
            'pdb_file': str(pdb_file),
            'total_molecules': len(profiles),
            'recommended_count': len(recommended),
            'elapsed_time': elapsed_time,
            'avg_final_score': sum(p.final_score for p in profiles) / len(profiles) if profiles else 0,
            'avg_admet_score': sum(p.admet_score for p in profiles) / len(profiles) if profiles else 0,
            'avg_vina_score': sum(p.vina_score for p in profiles) / len(profiles) if profiles else 0,
            'top_score': profiles[0].final_score if profiles else 0,
            'output_dir': str(target_output)
        }
        
        return result
    
    def _load_progress(self) -> Dict:
        """加载进度"""
        if self.progress_file.exists():
            with open(self.progress_file, 'r') as f:
                return json.load(f)
        return {'completed': [], 'failed': []}
    
    def _save_progress(self, target_name: str, status: str):
        """保存进度"""
        if status not in self.progress:
            self.progress[status] = []
        
        if target_name not in self.progress[status]:
            self.progress[status].append(target_name)
        
        with open(self.progress_file, 'w') as f:
            json.dump(self.progress, f, indent=2)
    
    def _generate_summary(self, stats: Dict):
        """生成汇总报告"""
        # 保存统计信息
        stats_file = self.output_base / 'batch_statistics.json'
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        # 汇总所有结果
        all_profiles = []
        
        for result in stats['results']:
            target_name = result['target_name']
            result_file = Path(result['output_dir']) / 'complete_results.csv'
            
            if result_file.exists():
                df = pd.read_csv(result_file)
                df['target'] = target_name
                all_profiles.append(df)
        
        if all_profiles:
            # 合并所有结果
            combined_df = pd.concat(all_profiles, ignore_index=True)
            
            # 保存合并结果
            combined_file = self.output_base / 'all_targets_combined.csv'
            combined_df.to_csv(combined_file, index=False)
            logger.info(f"合并结果已保存: {combined_file}")
            
            # 生成Top排名
            recommended = combined_df[combined_df['recommended'] == True]
            top_overall = recommended.nlargest(50, 'final_score')
            
            top_file = self.output_base / 'top50_all_targets.csv'
            top_overall.to_csv(top_file, index=False)
            logger.info(f"Top 50 候选已保存: {top_file}")
            
            # 每个靶点的最佳候选
            best_per_target = combined_df.loc[combined_df.groupby('target')['final_score'].idxmax()]
            best_file = self.output_base / 'best_per_target.csv'
            best_per_target.to_csv(best_file, index=False)
            logger.info(f"每个靶点最佳候选已保存: {best_file}")
            
            # 生成HTML报告
            self._generate_html_report(combined_df, stats)
    
    def _generate_html_report(self, df: pd.DataFrame, stats: Dict):
        """生成HTML报告"""
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>批量药物发现报告</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .stats {{ background: #ecf0f1; padding: 20px; border-radius: 5px; }}
        .recommended {{ color: green; font-weight: bold; }}
    </style>
</head>
<body>
    <h1>🏭 批量药物发现报告</h1>
    
    <div class="stats">
        <h2>总体统计</h2>
        <p><strong>处理靶点数:</strong> {stats['total_targets']}</p>
        <p><strong>成功完成:</strong> {stats['completed']}</p>
        <p><strong>失败:</strong> {stats['failed']}</p>
        <p><strong>成功率:</strong> {stats['success_rate']:.1f}%</p>
        <p><strong>总分子数:</strong> {stats['total_molecules']}</p>
        <p><strong>推荐候选:</strong> {stats['total_recommended']}</p>
        <p><strong>开始时间:</strong> {stats['start_time']}</p>
        <p><strong>结束时间:</strong> {stats['end_time']}</p>
    </div>
    
    <h2>Top 20 推荐候选</h2>
    <table>
        <tr>
            <th>排名</th>
            <th>靶点</th>
            <th>分子ID</th>
            <th>综合评分</th>
            <th>ADMET</th>
            <th>类药性</th>
            <th>Vina</th>
            <th>推荐</th>
        </tr>
"""
        
        # Top 20
        top20 = df.nlargest(20, 'final_score')
        for idx, row in top20.iterrows():
            recommended_tag = '<span class="recommended">✓</span>' if row['recommended'] else ''
            html_content += f"""
        <tr>
            <td>{row.get('rank', '-')}</td>
            <td>{row['target']}</td>
            <td>{row['mol_id']}</td>
            <td>{row['final_score']:.1f}</td>
            <td>{row['admet_score']:.1f}</td>
            <td>{row['druglikeness_score']:.1f}</td>
            <td>{row['vina_score']:.2f}</td>
            <td>{recommended_tag}</td>
        </tr>
"""
        
        html_content += """
    </table>
    
    <h2>各靶点汇总</h2>
    <table>
        <tr>
            <th>靶点</th>
            <th>总分子数</th>
            <th>推荐数</th>
            <th>平均评分</th>
            <th>最高评分</th>
        </tr>
"""
        
        # 各靶点统计
        for result in stats['results']:
            html_content += f"""
        <tr>
            <td>{result['target_name']}</td>
            <td>{result['total_molecules']}</td>
            <td>{result['recommended_count']}</td>
            <td>{result['avg_final_score']:.1f}</td>
            <td>{result['top_score']:.1f}</td>
        </tr>
"""
        
        html_content += """
    </table>
</body>
</html>
"""
        
        # 保存HTML
        html_file = self.output_base / 'batch_report.html'
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logger.info(f"HTML报告已保存: {html_file}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='批量自动化药物发现系统',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 处理文件夹中所有PDB
  python batch_discovery.py --pdb_folder ./pdbs --output batch_results
  
  # 快速模式(不含对接)
  python batch_discovery.py --pdb_folder ./pdbs --no_docking --candidates 50
  
  # 高质量模式
  python batch_discovery.py --pdb_folder ./pdbs --candidates 500 --sampling_steps 1000 --device cuda
  
  # 限制处理数量
  python batch_discovery.py --pdb_folder ./pdbs --max_targets 5
        """
    )
    
    # 必需参数
    parser.add_argument('--pdb_folder', type=str, required=True,
                       help='PDB文件夹路径')
    
    # 输出参数
    parser.add_argument('--output', type=str, default='batch_results',
                       help='输出目录 (默认: batch_results)')
    
    # 生成参数
    parser.add_argument('--candidates', type=int, default=100,
                       help='每个靶点生成的候选数 (默认: 100)')
    parser.add_argument('--num_atoms', type=int, nargs=2, default=[15, 35],
                       help='原子数范围 (默认: 15 35)')
    parser.add_argument('--sampling_steps', type=int, default=500,
                       help='采样步数 (默认: 500)')
    
    # 处理参数
    parser.add_argument('--no_docking', action='store_true',
                       help='禁用分子对接')
    parser.add_argument('--max_targets', type=int, default=None,
                       help='最大处理靶点数 (默认: 全部)')
    parser.add_argument('--no_resume', action='store_true',
                       help='不使用断点续传')
    
    # 设备参数
    parser.add_argument('--device', type=str, default=None,
                       help='计算设备 (cpu/cuda/cuda:0, 默认: 自动)')
    
    args = parser.parse_args()
    
    # 检查PDB文件夹
    if not os.path.exists(args.pdb_folder):
        print(f"错误: PDB文件夹不存在: {args.pdb_folder}")
        sys.exit(1)
    
    # 初始化系统
    print("=" * 80)
    print("批量自动化药物发现系统")
    print("=" * 80)
    print(f"PDB文件夹: {args.pdb_folder}")
    print(f"输出目录: {args.output}")
    print(f"每靶点候选数: {args.candidates}")
    print(f"采样步数: {args.sampling_steps}")
    print(f"对接: {'禁用' if args.no_docking else '启用'}")
    print(f"设备: {args.device or '自动'}")
    print("=" * 80)
    
    # 运行批量处理
    batch_system = BatchDrugDiscovery(
        pdb_folder=args.pdb_folder,
        output_base=args.output,
        device=args.device,
        resume=not args.no_resume
    )
    
    stats = batch_system.run_batch_processing(
        candidates_per_target=args.candidates,
        num_atoms_range=tuple(args.num_atoms),
        sampling_steps=args.sampling_steps,
        enable_docking=not args.no_docking,
        max_targets=args.max_targets
    )
    
    print("\n" + "=" * 80)
    print("✓ 批量处理完成!")
    print("=" * 80)
    print(f"查看结果: {args.output}/")
    print(f"  - all_targets_combined.csv  (所有结果)")
    print(f"  - top50_all_targets.csv     (Top 50候选)")
    print(f"  - best_per_target.csv       (每靶点最佳)")
    print(f"  - batch_report.html         (HTML报告)")
    print(f"  - batch_statistics.json     (统计信息)")


if __name__ == '__main__':
    main()
