"""
DeepST 扩展模块 - AI药物发现与设计 (AIDD)
Spatial Transcriptomics for Drug Discovery

功能:
1. 空间特异性靶点识别
2. 药物响应预测
3. 肿瘤微环境分析
4. 生物标志物发现
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from typing import List, Dict, Tuple
import warnings

warnings.filterwarnings('ignore')


class SpatialTargetIdentifier:
    """空间特异性靶点识别器"""
    
    def __init__(self, adata, domain_key: str = 'DeepST_domain'):
        """
        参数:
            adata: AnnData对象
            domain_key: 空间域标签的key
        """
        self.adata = adata
        self.domain_key = domain_key
        self.marker_genes = {}
    
    def identify_domain_markers(
        self,
        method: str = 'wilcoxon',
        top_n: int = 50,
        min_fold_change: float = 1.5,
        max_pval: float = 0.05
    ) -> pd.DataFrame:
        """
        识别每个空间域的标志基因
        
        参数:
            method: 统计检验方法 ('wilcoxon', 't-test')
            top_n: 每个域返回的顶部基因数量
            min_fold_change: 最小倍数变化
            max_pval: 最大p值
            
        返回:
            marker_df: 标志基因数据框
        """
        print(f"\n{'='*60}")
        print("识别空间域标志基因")
        print(f"{'='*60}")
        
        # 使用 Scanpy 进行差异表达分析
        sc.tl.rank_genes_groups(
            self.adata,
            groupby=self.domain_key,
            method=method,
            n_genes=top_n
        )
        
        # 提取结果
        marker_list = []
        
        domains = self.adata.obs[self.domain_key].unique()
        
        for domain in domains:
            genes = sc.get.rank_genes_groups_df(
                self.adata,
                group=domain
            )
            
            # 过滤
            genes = genes[
                (genes['logfoldchanges'] >= np.log2(min_fold_change)) &
                (genes['pvals_adj'] <= max_pval)
            ]
            
            genes['domain'] = domain
            marker_list.append(genes.head(top_n))
            
            self.marker_genes[domain] = genes['names'].head(top_n).tolist()
            
            print(f"域 {domain}: {len(genes)} 个标志基因")
        
        marker_df = pd.concat(marker_list, ignore_index=True)
        
        print(f"\n✓ 总计 {len(marker_df)} 个域特异性标志基因")
        
        return marker_df
    
    def identify_potential_targets(
        self,
        target_domain: str,
        druggable_genes: List[str] = None,
        expression_threshold: float = 1.0
    ) -> pd.DataFrame:
        """
        识别特定域的潜在药物靶点
        
        参数:
            target_domain: 目标域
            druggable_genes: 可成药基因列表
            expression_threshold: 表达阈值
            
        返回:
            targets_df: 潜在靶点数据框
        """
        print(f"\n寻找域 {target_domain} 的潜在药物靶点...")
        
        if target_domain not in self.marker_genes:
            raise ValueError(f"请先运行 identify_domain_markers()")
        
        markers = self.marker_genes[target_domain]
        
        # 计算域内平均表达
        domain_mask = self.adata.obs[self.domain_key] == target_domain
        domain_expression = self.adata[domain_mask].to_df().mean()
        
        # 过滤高表达基因
        high_expr_markers = [
            gene for gene in markers
            if gene in domain_expression.index and
            domain_expression[gene] >= expression_threshold
        ]
        
        # 如果提供可成药基因列表,取交集
        if druggable_genes is not None:
            druggable_set = set(druggable_genes)
            potential_targets = [
                gene for gene in high_expr_markers
                if gene in druggable_set
            ]
            print(f"  可成药靶点: {len(potential_targets)}")
        else:
            potential_targets = high_expr_markers
        
        # 创建结果数据框
        targets_df = pd.DataFrame({
            'gene': potential_targets,
            'mean_expression': [domain_expression[g] for g in potential_targets],
            'domain': target_domain
        })
        
        targets_df = targets_df.sort_values('mean_expression', ascending=False)
        
        print(f"✓ 识别到 {len(targets_df)} 个潜在靶点")
        
        return targets_df
    
    def visualize_target_expression(
        self,
        target_gene: str,
        save_path: str = None
    ):
        """可视化靶点基因的空间表达"""
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # 空间表达图
        sc.pl.spatial(
            self.adata,
            color=target_gene,
            spot_size=150,
            title=f'{target_gene} Expression',
            ax=axes[0],
            show=False,
            frameon=False,
            cmap='viridis'
        )
        
        # 小提琴图 - 按域
        sc.pl.violin(
            self.adata,
            keys=target_gene,
            groupby=self.domain_key,
            ax=axes[1],
            show=False
        )
        axes[1].set_title(f'{target_gene} Expression by Domain')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        plt.show()


class TumorMicroenvironmentAnalyzer:
    """肿瘤微环境分析器"""
    
    def __init__(self, adata, domain_key: str = 'DeepST_domain'):
        self.adata = adata
        self.domain_key = domain_key
    
    def annotate_tumor_regions(
        self,
        tumor_markers: List[str],
        immune_markers: List[str],
        stromal_markers: List[str]
    ) -> pd.DataFrame:
        """
        注释肿瘤、免疫和基质区域
        
        参数:
            tumor_markers: 肿瘤标志物列表
            immune_markers: 免疫标志物列表
            stromal_markers: 基质标志物列表
            
        返回:
            annotation_df: 注释结果
        """
        print(f"\n{'='*60}")
        print("注释肿瘤微环境")
        print(f"{'='*60}")
        
        # 计算每个类型的标志物得分
        def calc_signature_score(markers):
            available = [m for m in markers if m in self.adata.var_names]
            if len(available) == 0:
                return np.zeros(self.adata.n_obs)
            return self.adata[:, available].to_df().mean(axis=1).values
        
        tumor_score = calc_signature_score(tumor_markers)
        immune_score = calc_signature_score(immune_markers)
        stromal_score = calc_signature_score(stromal_markers)
        
        # 分配注释
        scores = np.column_stack([tumor_score, immune_score, stromal_score])
        annotation = np.argmax(scores, axis=1)
        
        annotation_labels = ['Tumor', 'Immune', 'Stromal']
        self.adata.obs['TME_annotation'] = [annotation_labels[i] for i in annotation]
        
        # 统计
        for label in annotation_labels:
            count = (self.adata.obs['TME_annotation'] == label).sum()
            print(f"{label}: {count} spots ({count/len(self.adata)*100:.1f}%)")
        
        # 创建结果数据框
        annotation_df = pd.DataFrame({
            'spot_id': self.adata.obs_names,
            'annotation': self.adata.obs['TME_annotation'],
            'tumor_score': tumor_score,
            'immune_score': immune_score,
            'stromal_score': stromal_score
        })
        
        return annotation_df
    
    def analyze_immune_infiltration(
        self,
        immune_markers: Dict[str, List[str]]
    ) -> pd.DataFrame:
        """
        分析免疫浸润
        
        参数:
            immune_markers: 免疫细胞类型标志物字典
                {'CD8_T': [...], 'CD4_T': [...], 'Macrophage': [...], ...}
                
        返回:
            infiltration_df: 浸润分析结果
        """
        print(f"\n分析免疫浸润...")
        
        infiltration_scores = {}
        
        for cell_type, markers in immune_markers.items():
            available = [m for m in markers if m in self.adata.var_names]
            if len(available) > 0:
                score = self.adata[:, available].to_df().mean(axis=1).values
                infiltration_scores[f'{cell_type}_score'] = score
                print(f"  {cell_type}: {len(available)} 个标志物")
        
        infiltration_df = pd.DataFrame(infiltration_scores)
        infiltration_df['spot_id'] = self.adata.obs_names
        infiltration_df['domain'] = self.adata.obs[self.domain_key].values
        
        return infiltration_df
    
    def identify_tumor_boundary(
        self,
        tumor_label: str = 'Tumor',
        k_neighbors: int = 6
    ) -> np.ndarray:
        """
        识别肿瘤边界
        
        参数:
            tumor_label: 肿瘤区域标签
            k_neighbors: 邻居数量
            
        返回:
            boundary_mask: 边界mask
        """
        print(f"\n识别肿瘤边界...")
        
        from sklearn.neighbors import NearestNeighbors
        
        spatial_coords = self.adata.obsm['spatial']
        is_tumor = (self.adata.obs['TME_annotation'] == tumor_label).values
        
        # 找邻居
        nbrs = NearestNeighbors(n_neighbors=k_neighbors).fit(spatial_coords)
        _, indices = nbrs.kneighbors(spatial_coords)
        
        # 边界定义:肿瘤spot但至少有一个非肿瘤邻居
        boundary_mask = np.zeros(len(is_tumor), dtype=bool)
        
        for i in range(len(is_tumor)):
            if is_tumor[i]:
                neighbors_tumor = is_tumor[indices[i]]
                if not neighbors_tumor.all():
                    boundary_mask[i] = True
        
        self.adata.obs['is_boundary'] = boundary_mask
        
        n_boundary = boundary_mask.sum()
        print(f"✓ 识别到 {n_boundary} 个边界spots")
        
        return boundary_mask


class DrugResponsePredictor:
    """药物响应预测器"""
    
    def __init__(self, adata, domain_key: str = 'DeepST_domain'):
        self.adata = adata
        self.domain_key = domain_key
        self.model = None
    
    def prepare_features(
        self,
        use_embeddings: bool = True,
        use_gene_expression: bool = True,
        gene_subset: List[str] = None
    ) -> np.ndarray:
        """
        准备预测特征
        
        参数:
            use_embeddings: 是否使用DeepST嵌入
            use_gene_expression: 是否使用基因表达
            gene_subset: 基因子集
            
        返回:
            features: 特征矩阵
        """
        feature_list = []
        
        if use_embeddings and 'DeepST_embed' in self.adata.obsm:
            feature_list.append(self.adata.obsm['DeepST_embed'])
        
        if use_gene_expression:
            if gene_subset is not None:
                expr = self.adata[:, gene_subset].to_df().values
            else:
                expr = self.adata.to_df().values
            feature_list.append(expr)
        
        features = np.hstack(feature_list) if len(feature_list) > 1 else feature_list[0]
        
        return features
    
    def train_response_model(
        self,
        response_labels: np.ndarray,
        features: np.ndarray = None,
        model_type: str = 'random_forest'
    ):
        """
        训练药物响应预测模型
        
        参数:
            response_labels: 响应标签 (0: 不响应, 1: 响应)
            features: 特征矩阵
            model_type: 模型类型
        """
        print(f"\n{'='*60}")
        print("训练药物响应预测模型")
        print(f"{'='*60}")
        
        if features is None:
            features = self.prepare_features()
        
        # 选择模型
        if model_type == 'random_forest':
            self.model = RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=42
            )
        elif model_type == 'logistic':
            self.model = LogisticRegression(
                max_iter=1000,
                random_state=42
            )
        else:
            raise ValueError(f"未知模型类型: {model_type}")
        
        # 训练
        self.model.fit(features, response_labels)
        
        # 评估
        from sklearn.metrics import accuracy_score, roc_auc_score
        
        pred_labels = self.model.predict(features)
        pred_proba = self.model.predict_proba(features)[:, 1]
        
        accuracy = accuracy_score(response_labels, pred_labels)
        auc = roc_auc_score(response_labels, pred_proba)
        
        print(f"✓ 训练完成")
        print(f"  准确率: {accuracy:.4f}")
        print(f"  AUC: {auc:.4f}")
        
        # 保存预测概率
        self.adata.obs['drug_response_prob'] = pred_proba
    
    def predict_spatial_response(
        self,
        features: np.ndarray = None
    ) -> np.ndarray:
        """
        预测空间药物响应
        
        返回:
            response_proba: 响应概率
        """
        if self.model is None:
            raise ValueError("请先训练模型")
        
        if features is None:
            features = self.prepare_features()
        
        response_proba = self.model.predict_proba(features)[:, 1]
        
        return response_proba
    
    def identify_resistance_regions(
        self,
        threshold: float = 0.5
    ) -> pd.DataFrame:
        """
        识别耐药区域
        
        参数:
            threshold: 响应概率阈值
            
        返回:
            resistance_df: 耐药区域信息
        """
        if 'drug_response_prob' not in self.adata.obs:
            raise ValueError("请先训练并预测")
        
        response_prob = self.adata.obs['drug_response_prob'].values
        is_resistant = response_prob < threshold
        
        self.adata.obs['is_resistant'] = is_resistant
        
        # 按域统计
        resistance_by_domain = []
        
        for domain in self.adata.obs[self.domain_key].unique():
            domain_mask = self.adata.obs[self.domain_key] == domain
            domain_resistant = is_resistant[domain_mask].sum()
            domain_total = domain_mask.sum()
            
            resistance_by_domain.append({
                'domain': domain,
                'resistant_spots': domain_resistant,
                'total_spots': domain_total,
                'resistance_rate': domain_resistant / domain_total
            })
        
        resistance_df = pd.DataFrame(resistance_by_domain)
        resistance_df = resistance_df.sort_values('resistance_rate', ascending=False)
        
        print(f"\n耐药区域统计:")
        print(resistance_df.to_string(index=False))
        
        return resistance_df


class BiomarkerDiscovery:
    """生物标志物发现"""
    
    def __init__(self, adata, domain_key: str = 'DeepST_domain'):
        self.adata = adata
        self.domain_key = domain_key
    
    def discover_prognostic_markers(
        self,
        survival_time: np.ndarray,
        event_observed: np.ndarray,
        top_n: int = 20
    ) -> pd.DataFrame:
        """
        发现预后标志物
        
        参数:
            survival_time: 生存时间
            event_observed: 事件发生标记
            top_n: 返回顶部基因数量
            
        返回:
            markers_df: 标志物数据框
        """
        print(f"\n{'='*60}")
        print("发现预后生物标志物")
        print(f"{'='*60}")
        
        from scipy.stats import spearmanr
        
        correlations = []
        
        gene_names = self.adata.var_names
        expression = self.adata.to_df()
        
        for gene in gene_names:
            corr, pval = spearmanr(expression[gene], survival_time)
            correlations.append({
                'gene': gene,
                'correlation': corr,
                'pvalue': pval
            })
        
        markers_df = pd.DataFrame(correlations)
        markers_df['abs_correlation'] = np.abs(markers_df['correlation'])
        markers_df = markers_df.sort_values('abs_correlation', ascending=False)
        
        # 多重检验校正
        from statsmodels.stats.multitest import multipletests
        
        _, markers_df['pvalue_adj'], _, _ = multipletests(
            markers_df['pvalue'],
            method='fdr_bh'
        )
        
        # 过滤显著基因
        significant = markers_df[markers_df['pvalue_adj'] < 0.05]
        
        print(f"✓ 发现 {len(significant)} 个显著预后标志物")
        
        return markers_df.head(top_n)
    
    def spatial_heterogeneity_score(self) -> pd.DataFrame:
        """
        计算空间异质性得分
        
        返回:
            heterogeneity_df: 异质性得分
        """
        print(f"\n计算空间异质性...")
        
        from sklearn.metrics import pairwise_distances
        from sklearn.neighbors import NearestNeighbors
        
        embeddings = self.adata.obsm['DeepST_embed']
        spatial_coords = self.adata.obsm['spatial']
        
        # 找邻居
        nbrs = NearestNeighbors(n_neighbors=7).fit(spatial_coords)
        _, indices = nbrs.kneighbors(spatial_coords)
        
        # 计算每个spot与邻居的差异
        heterogeneity_scores = []
        
        for i in range(len(embeddings)):
            neighbor_embeddings = embeddings[indices[i][1:]]  # 排除自己
            distances = pairwise_distances(
                [embeddings[i]],
                neighbor_embeddings
            )[0]
            heterogeneity_scores.append(distances.mean())
        
        self.adata.obs['heterogeneity_score'] = heterogeneity_scores
        
        # 按域统计
        heterogeneity_by_domain = []
        
        for domain in self.adata.obs[self.domain_key].unique():
            domain_mask = self.adata.obs[self.domain_key] == domain
            mean_score = np.array(heterogeneity_scores)[domain_mask].mean()
            
            heterogeneity_by_domain.append({
                'domain': domain,
                'mean_heterogeneity': mean_score
            })
        
        heterogeneity_df = pd.DataFrame(heterogeneity_by_domain)
        
        print(f"✓ 异质性计算完成")
        print(heterogeneity_df.to_string(index=False))
        
        return heterogeneity_df


def aidd_example_workflow():
    """AIDD 应用示例工作流"""
    
    print("\n" + "="*60)
    print("DeepST AIDD 应用示例")
    print("="*60 + "\n")
    
    # 假设已经运行了 DeepST 分析
    # adata = sc.read_h5ad('deepst_results.h5ad')
    
    # 生成示例数据
    from sklearn.datasets import make_blobs
    import anndata as ad
    
    n_spots = 1000
    spatial_coords, true_labels = make_blobs(
        n_samples=n_spots,
        n_features=2,
        centers=7,
        random_state=42
    )
    
    gene_expression = np.random.negative_binomial(5, 0.3, (n_spots, 100))
    
    adata = ad.AnnData(X=gene_expression)
    adata.obsm['spatial'] = spatial_coords
    adata.obs['DeepST_domain'] = true_labels.astype(str)
    adata.obsm['DeepST_embed'] = np.random.randn(n_spots, 128)
    
    # 添加模拟基因名
    adata.var_names = [f'Gene_{i}' for i in range(100)]
    
    print("1. 识别空间靶点")
    print("-" * 60)
    target_id = SpatialTargetIdentifier(adata)
    markers = target_id.identify_domain_markers(top_n=10)
    
    targets = target_id.identify_potential_targets(
        target_domain='0',
        expression_threshold=0.5
    )
    
    print("\n2. 肿瘤微环境分析")
    print("-" * 60)
    tme = TumorMicroenvironmentAnalyzer(adata)
    
    # 模拟标志物
    tumor_markers = ['Gene_0', 'Gene_1', 'Gene_2']
    immune_markers = ['Gene_3', 'Gene_4', 'Gene_5']
    stromal_markers = ['Gene_6', 'Gene_7', 'Gene_8']
    
    annotation = tme.annotate_tumor_regions(
        tumor_markers,
        immune_markers,
        stromal_markers
    )
    
    print("\n3. 药物响应预测")
    print("-" * 60)
    predictor = DrugResponsePredictor(adata)
    
    # 模拟响应标签
    response_labels = np.random.binomial(1, 0.7, n_spots)
    
    predictor.train_response_model(response_labels)
    
    resistance = predictor.identify_resistance_regions(threshold=0.5)
    
    print("\n4. 生物标志物发现")
    print("-" * 60)
    biomarker = BiomarkerDiscovery(adata)
    
    # 模拟生存数据
    survival_time = np.random.exponential(10, n_spots)
    event_observed = np.random.binomial(1, 0.6, n_spots)
    
    prognostic = biomarker.discover_prognostic_markers(
        survival_time,
        event_observed,
        top_n=10
    )
    
    heterogeneity = biomarker.spatial_heterogeneity_score()
    
    print(f"\n{'='*60}")
    print("AIDD 分析完成!")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    aidd_example_workflow()
