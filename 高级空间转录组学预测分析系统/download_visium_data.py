"""
10x Genomics Visium数据下载和准备脚本
自动下载、解压和组织数据
"""

import os
import urllib.request
import tarfile
import gzip
import shutil
from pathlib import Path
import scanpy as sc
import pandas as pd


class VisiumDataDownloader:
    """10x Genomics Visium数据下载器"""
    
    # 10x Genomics公开数据集URL
    DATASETS = {
        'breast_cancer': {
            'name': 'Human Breast Cancer (Block A Section 1)',
            'base_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/',
            'files': {
                'filtered_matrix': 'V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5',
                'raw_matrix': 'V1_Breast_Cancer_Block_A_Section_1_raw_feature_bc_matrix.h5',
                'spatial': 'V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz',
                'molecule_info': 'V1_Breast_Cancer_Block_A_Section_1_molecule_info.h5',
                'analysis': 'V1_Breast_Cancer_Block_A_Section_1_analysis.tar.gz',
                'cloupe': 'V1_Breast_Cancer_Block_A_Section_1_cloupe.cloupe'
            }
        },
        'breast_cancer_ff': {
            'name': 'Human Breast Cancer Fresh Frozen',
            'base_url': 'https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FreshFrozen_Human_Breast_Cancer/',
            'files': {
                'filtered_matrix': 'CytAssist_FreshFrozen_Human_Breast_Cancer_filtered_feature_bc_matrix.h5',
                'raw_matrix': 'CytAssist_FreshFrozen_Human_Breast_Cancer_raw_feature_bc_matrix.h5',
                'spatial': 'CytAssist_FreshFrozen_Human_Breast_Cancer_spatial.tar.gz',
                'analysis': 'CytAssist_FreshFrozen_Human_Breast_Cancer_analysis.tar.gz',
            }
        },
        'mouse_brain': {
            'name': 'Mouse Brain Serial Section 1',
            'base_url': 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/',
            'files': {
                'filtered_matrix': 'V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5',
                'spatial': 'V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz',
                'analysis': 'V1_Mouse_Brain_Sagittal_Anterior_analysis.tar.gz',
            }
        }
    }
    
    def __init__(self, output_dir='./visium_data'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
    def download_file(self, url, output_path, desc="Downloading"):
        """下载单个文件并显示进度"""
        print(f"{desc}: {output_path.name}")
        
        try:
            # 检查文件是否已存在
            if output_path.exists():
                print(f"  文件已存在,跳过下载")
                return True
            
            # 下载文件
            def reporthook(block_num, block_size, total_size):
                downloaded = block_num * block_size
                if total_size > 0:
                    percent = min(downloaded * 100 / total_size, 100)
                    mb_downloaded = downloaded / (1024 * 1024)
                    mb_total = total_size / (1024 * 1024)
                    print(f"\r  进度: {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", 
                          end='', flush=True)
            
            urllib.request.urlretrieve(url, output_path, reporthook)
            print()  # 换行
            return True
            
        except Exception as e:
            print(f"\n  下载失败: {e}")
            return False
    
    def extract_tarball(self, tar_path, extract_to):
        """解压tar.gz文件"""
        print(f"解压: {tar_path.name}")
        
        try:
            with tarfile.open(tar_path, 'r:gz') as tar:
                tar.extractall(extract_to)
            print(f"  解压完成")
            return True
        except Exception as e:
            print(f"  解压失败: {e}")
            return False
    
    def download_dataset(self, dataset_key='breast_cancer', 
                        download_all=False, required_files=None):
        """
        下载指定的数据集
        
        参数:
            dataset_key: 数据集标识 ('breast_cancer', 'mouse_brain', etc.)
            download_all: 是否下载所有文件
            required_files: 必需文件列表 (例如 ['filtered_matrix', 'spatial'])
        """
        if dataset_key not in self.DATASETS:
            print(f"未知数据集: {dataset_key}")
            print(f"可用数据集: {list(self.DATASETS.keys())}")
            return False
        
        dataset = self.DATASETS[dataset_key]
        dataset_dir = self.output_dir / dataset_key
        dataset_dir.mkdir(exist_ok=True, parents=True)
        
        print("="*80)
        print(f"下载数据集: {dataset['name']}")
        print(f"保存位置: {dataset_dir}")
        print("="*80)
        
        # 确定要下载的文件
        if download_all:
            files_to_download = dataset['files'].keys()
        elif required_files:
            files_to_download = required_files
        else:
            # 默认下载最基本的文件
            files_to_download = ['filtered_matrix', 'spatial']
        
        success = True
        for file_key in files_to_download:
            if file_key not in dataset['files']:
                print(f"警告: 文件 '{file_key}' 不存在于此数据集")
                continue
            
            filename = dataset['files'][file_key]
            url = dataset['base_url'] + filename
            output_path = dataset_dir / filename
            
            # 下载文件
            if not self.download_file(url, output_path):
                success = False
                continue
            
            # 如果是压缩文件,解压
            if filename.endswith('.tar.gz'):
                if not self.extract_tarball(output_path, dataset_dir):
                    success = False
        
        if success:
            print("\n" + "="*80)
            print("下载完成!")
            print("="*80)
            print(f"数据位置: {dataset_dir}")
            
            # 列出下载的文件
            print("\n下载的文件:")
            for item in sorted(dataset_dir.rglob('*')):
                if item.is_file():
                    size_mb = item.stat().st_size / (1024 * 1024)
                    print(f"  {item.relative_to(dataset_dir)} ({size_mb:.2f} MB)")
        
        return success
    
    def load_dataset(self, dataset_key='breast_cancer', filtered=True):
        """
        加载已下载的数据集为AnnData对象
        
        参数:
            dataset_key: 数据集标识
            filtered: 使用过滤后的矩阵 (True) 或原始矩阵 (False)
        """
        dataset_dir = self.output_dir / dataset_key
        
        if not dataset_dir.exists():
            print(f"数据集目录不存在: {dataset_dir}")
            print("请先下载数据集")
            return None
        
        print(f"加载数据集: {dataset_key}")
        
        # 查找h5文件
        matrix_type = 'filtered' if filtered else 'raw'
        h5_files = list(dataset_dir.glob(f'*{matrix_type}*.h5'))
        
        if not h5_files:
            print(f"未找到 {matrix_type} h5文件")
            return None
        
        h5_file = h5_files[0]
        print(f"  矩阵文件: {h5_file.name}")
        
        # 加载数据
        adata = sc.read_10x_h5(h5_file)
        
        # 查找spatial文件夹
        spatial_dirs = list(dataset_dir.glob('**/spatial'))
        
        if spatial_dirs:
            spatial_dir = spatial_dirs[0]
            print(f"  空间数据: {spatial_dir}")
            
            # 加载空间坐标
            tissue_positions = None
            for pos_file in ['tissue_positions.csv', 'tissue_positions_list.csv']:
                pos_path = spatial_dir / pos_file
                if pos_path.exists():
                    tissue_positions = pd.read_csv(pos_path, header=None if 'list' in pos_file else 0)
                    break
            
            if tissue_positions is not None:
                # 提取坐标 (通常在列4和列5, 或named as 'pxl_col_in_fullres' 和 'pxl_row_in_fullres')
                if tissue_positions.shape[1] >= 6:
                    # 无header格式
                    coords = tissue_positions.iloc[:, [4, 5]].values
                    barcodes = tissue_positions.iloc[:, 0].values
                else:
                    # 有header格式
                    coords = tissue_positions[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values
                    barcodes = tissue_positions['barcode'].values
                
                # 匹配barcodes
                common_barcodes = adata.obs_names.intersection(barcodes)
                
                if len(common_barcodes) > 0:
                    # 创建坐标映射
                    coord_dict = dict(zip(barcodes, coords))
                    spatial_coords = np.array([coord_dict.get(bc, [0, 0]) for bc in adata.obs_names])
                    adata.obsm['spatial'] = spatial_coords
                    print(f"  成功加载 {len(common_barcodes)} 个spots的空间坐标")
                else:
                    print("  警告: barcodes不匹配")
            
            # 加载组织图像
            img_files = {
                'hires': spatial_dir / 'tissue_hires_image.png',
                'lowres': spatial_dir / 'tissue_lowres_image.png'
            }
            
            images = {}
            for img_type, img_path in img_files.items():
                if img_path.exists():
                    from PIL import Image
                    images[img_type] = np.array(Image.open(img_path))
                    print(f"  加载 {img_type} 图像: {images[img_type].shape}")
            
            if images:
                adata.uns['spatial'] = {
                    dataset_key: {
                        'images': images,
                        'scalefactors': {
                            'tissue_hires_scalef': 1.0,
                            'tissue_lowres_scalef': 1.0
                        }
                    }
                }
        
        print(f"\n数据集加载完成:")
        print(f"  Spots: {adata.n_obs}")
        print(f"  Genes: {adata.n_vars}")
        
        return adata


def prepare_analysis_ready_data(adata, output_path='./processed_data.h5ad'):
    """
    准备分析就绪的数据
    执行基本QC和预处理
    """
    print("\n准备分析数据...")
    
    # 1. 基本统计
    print("  计算QC指标...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    # 2. 过滤
    print("  过滤低质量spots...")
    print(f"    过滤前: {adata.n_obs} spots")
    
    # 过滤基因数太少的spots
    sc.pp.filter_cells(adata, min_genes=200)
    # 过滤在太少spots中表达的基因
    sc.pp.filter_genes(adata, min_cells=3)
    
    print(f"    过滤后: {adata.n_obs} spots, {adata.n_vars} genes")
    
    # 3. 归一化
    print("  归一化...")
    adata.layers['counts'] = adata.X.copy()  # 保存原始计数
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 4. 高变基因
    print("  识别高变基因...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    print(f"    发现 {adata.var['highly_variable'].sum()} 个高变基因")
    
    # 5. 保存
    print(f"  保存到: {output_path}")
    adata.write(output_path)
    
    print("  准备完成!")
    return adata


# 主函数
if __name__ == "__main__":
    import argparse
    import numpy as np
    
    parser = argparse.ArgumentParser(description='下载和准备10x Genomics Visium数据')
    parser.add_argument('--dataset', type=str, default='breast_cancer',
                       choices=['breast_cancer', 'breast_cancer_ff', 'mouse_brain'],
                       help='数据集名称')
    parser.add_argument('--output-dir', type=str, default='./visium_data',
                       help='输出目录')
    parser.add_argument('--download-all', action='store_true',
                       help='下载所有文件 (否则只下载基本文件)')
    parser.add_argument('--skip-download', action='store_true',
                       help='跳过下载,只加载数据')
    
    args = parser.parse_args()
    
    # 创建下载器
    downloader = VisiumDataDownloader(output_dir=args.output_dir)
    
    # 下载数据
    if not args.skip_download:
        success = downloader.download_dataset(
            dataset_key=args.dataset,
            download_all=args.download_all
        )
        
        if not success:
            print("下载失败!")
            exit(1)
    
    # 加载数据
    print("\n" + "="*80)
    adata = downloader.load_dataset(dataset_key=args.dataset)
    
    if adata is None:
        print("数据加载失败!")
        exit(1)
    
    # 准备数据
    output_path = Path(args.output_dir) / args.dataset / 'processed_data.h5ad'
    adata = prepare_analysis_ready_data(adata, output_path=output_path)
    
    print("\n" + "="*80)
    print("所有步骤完成!")
    print("="*80)
    print(f"处理后的数据已保存到: {output_path}")
    print("\n下一步:")
    print(f"  python spatial_transcriptomics_prediction.py --data {output_path}")
