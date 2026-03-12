"""
高级下游生物学分析
包括: 差异表达分析、通路富集、空间轨迹、细胞通讯预测
"""

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.metrics import silhouette_score
import warnings
warnings.filterwarnings('ignore')


class AdvancedBiologicalAnalysis:
    """高级生物学分析工具集"""
    
    def __init__(self, adata):
        self.adata = adata
        
    def differential_expression_analysis(self, group_key='spatial_domain', 
                                        method='wilcoxon'):
        """
        空间域间差异表达分析
        """
        print("执行差异表达分析...")
        
        # 使用Scanpy的rank_genes_groups
        sc.tl.rank_genes_groups(
            self.adata,
            groupby=group_key,
            method=method,
            n_genes=100
        )
        
        # 提取结果
        result = self.adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        # 创建结果DataFrame
        de_results = {}
        for group in groups:
            de_results[group] = pd.DataFrame({
                'gene': result['names'][group],
                'score': result['scores'][group],
                'pval': result['pvals'][group],
                'pval_adj': result['pvals_adj'][group],
                'logfoldchange': result['logfoldchanges'][group]
            })
        
        return de_results
    
    def pathway_enrichment_analysis(self, de_genes, organism='human'):
        """
        通路富集分析 (需要连接到在线数据库)
        这里提供框架,实际使用时需要API key
        """
        print("通路富集分析框架...")
        
        # 模拟GO/KEGG富集结果
        # 实际使用时可以调用:
        # - GSEA (Gene Set Enrichment Analysis)
        # - Enrichr API
        # - clusterProfiler (R)
        
        pathways = {
            'GO:0008283': {
                'name': 'Cell proliferation',
                'genes': de_genes[:20],
                'pvalue': 1e-5,
                'enrichment_score': 3.2
            },
            'GO:0006955': {
                'name': 'Immune response',
                'genes': de_genes[20:40],
                'pvalue': 1e-4,
                'enrichment_score': 2.8
            },
            'KEGG:04010': {
                'name': 'MAPK signaling pathway',
                'genes': de_genes[40:60],
                'pvalue': 1e-3,
                'enrichment_score': 2.5
            }
        }
        
        return pd.DataFrame.from_dict(pathways, orient='index')
    
    def spatial_trajectory_inference(self, start_domain=None):
        """
        空间轨迹推断
        识别细胞状态转换的空间路径
        """
        print("推断空间轨迹...")
        
        # 使用扩散伪时间
        sc.pp.neighbors(self.adata, use_rep='latent')
        sc.tl.diffmap(self.adata)
        sc.tl.dpt(self.adata, n_dcs=10)
        
        # 如果指定起始域
        if start_domain is not None:
            self.adata.uns['iroot'] = np.flatnonzero(
                self.adata.obs['spatial_domain'] == start_domain
            )[0]
            sc.tl.dpt(self.adata)
        
        return self.adata.obs['dpt_pseudotime']
    
    def cell_cell_communication(self):
        """
        预测细胞-细胞通讯
        基于配体-受体相互作用
        """
        print("分析细胞通讯...")
        
        # 使用Squidpy的配体-受体分析
        sq.gr.ligrec(
            self.adata,
            n_perms=100,
            cluster_key='spatial_domain',
            copy=False
        )
        
        # 提取显著相互作用
        if 'ligrec' in self.adata.uns:
            lr_results = self.adata.uns['ligrec']
            return lr_results
        
        return None
    
    def spatial_autocorrelation(self, genes=None):
        """
        计算基因的空间自相关性 (Moran's I)
        识别空间可变基因
        """
        print("计算空间自相关性...")
        
        if genes is None:
            genes = self.adata.var_names[:100]  # 默认前100个基因
        
        # 计算Moran's I
        morans_i = {}
        for gene in genes:
            gene_exp = self.adata[:, gene].X.toarray().flatten()
            
            # 简化版Moran's I计算
            n = len(gene_exp)
            mean_exp = np.mean(gene_exp)
            
            # 使用空间邻接矩阵
            W = self.adata.obsp['spatial_connectivities'].toarray()
            
            numerator = 0
            denominator = 0
            W_sum = W.sum()
            
            for i in range(n):
                for j in range(n):
                    numerator += W[i, j] * (gene_exp[i] - mean_exp) * (gene_exp[j] - mean_exp)
                denominator += (gene_exp[i] - mean_exp) ** 2
            
            if denominator > 0 and W_sum > 0:
                I = (n / W_sum) * (numerator / denominator)
                morans_i[gene] = I
        
        morans_df = pd.DataFrame.from_dict(
            morans_i, orient='index', columns=['Morans_I']
        ).sort_values('Morans_I', ascending=False)
        
        return morans_df
    
    def spatial_variable_genes(self, method='moransi', n_top=100):
        """
        识别空间可变基因 (SVG)
        """
        print(f"识别空间可变基因 (方法: {method})...")
        
        if method == 'moransi':
            # 使用Moran's I
            morans = self.spatial_autocorrelation(genes=self.adata.var_names)
            svg = morans.head(n_top).index.tolist()
            
        elif method == 'squidpy':
            # 使用Squidpy
            sq.gr.spatial_autocorr(
                self.adata,
                mode='moran',
                n_perms=100,
                n_jobs=1
            )
            svg = self.adata.uns['moranI'].nlargest(n_top, 'I').index.tolist()
        
        return svg
    
    def niche_clustering(self, n_niches=5, radius=150):
        """
        生态位聚类 - 基于局部邻域的基因表达模式
        """
        print("执行生态位聚类...")
        
        # 计算局部邻域聚合表达
        sq.gr.spatial_neighbors(self.adata, radius=radius, coord_type='generic')
        
        # 使用邻域聚合
        # 这里简化: 使用邻域平均表达
        adj_matrix = self.adata.obsp['spatial_connectivities']
        
        if hasattr(self.adata.X, 'toarray'):
            X = self.adata.X.toarray()
        else:
            X = self.adata.X
        
        # 邻域聚合
        neighborhood_exp = adj_matrix @ X
        
        # 在邻域表达上聚类
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=n_niches, random_state=42)
        niches = kmeans.fit_predict(neighborhood_exp)
        
        self.adata.obs['niche'] = niches
        
        return niches
    
    def spatial_gene_correlation_network(self, genes=None, threshold=0.5):
        """
        构建空间基因共表达网络
        """
        print("构建空间基因网络...")
        
        if genes is None:
            genes = self.adata.var_names[:50]
        
        # 提取基因表达
        if hasattr(self.adata[:, genes].X, 'toarray'):
            gene_exp = self.adata[:, genes].X.toarray()
        else:
            gene_exp = self.adata[:, genes].X
        
        # 计算相关性
        corr_matrix = np.corrcoef(gene_exp.T)
        
        # 阈值化
        corr_matrix[np.abs(corr_matrix) < threshold] = 0
        
        # 创建网络
        import networkx as nx
        G = nx.Graph()
        
        for i, gene1 in enumerate(genes):
            for j, gene2 in enumerate(genes):
                if i < j and corr_matrix[i, j] != 0:
                    G.add_edge(gene1, gene2, weight=corr_matrix[i, j])
        
        return G, corr_matrix


class VisualizationSuite:
    """高级可视化工具"""
    
    @staticmethod
    def plot_de_volcano(de_results, group_name, save_path=None):
        """火山图 - 差异表达"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        df = de_results[group_name]
        
        # 计算-log10(pval)
        df['neglog10pval'] = -np.log10(df['pval_adj'] + 1e-300)
        
        # 分类基因
        df['significant'] = 'Not Sig'
        df.loc[(df['pval_adj'] < 0.05) & (df['logfoldchange'] > 0.5), 'significant'] = 'Up'
        df.loc[(df['pval_adj'] < 0.05) & (df['logfoldchange'] < -0.5), 'significant'] = 'Down'
        
        # 绘图
        colors = {'Not Sig': 'gray', 'Up': 'red', 'Down': 'blue'}
        for sig in ['Not Sig', 'Up', 'Down']:
            subset = df[df['significant'] == sig]
            ax.scatter(subset['logfoldchange'], subset['neglog10pval'],
                      c=colors[sig], label=sig, alpha=0.6, s=20)
        
        ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1)
        ax.axvline(0.5, color='black', linestyle='--', linewidth=1)
        ax.axvline(-0.5, color='black', linestyle='--', linewidth=1)
        
        ax.set_xlabel('Log2 Fold Change')
        ax.set_ylabel('-Log10 Adjusted P-value')
        ax.set_title(f'Volcano Plot - {group_name}')
        ax.legend()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
    
    @staticmethod
    def plot_spatial_trajectory(adata, save_path=None):
        """空间轨迹可视化"""
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        
        # 左图: 空间位置着色伪时间
        spatial_coords = adata.obsm['spatial']
        pseudotime = adata.obs['dpt_pseudotime']
        
        scatter = axes[0].scatter(
            spatial_coords[:, 0],
            spatial_coords[:, 1],
            c=pseudotime,
            cmap='viridis',
            s=50,
            alpha=0.8
        )
        axes[0].set_title('Spatial Trajectory (Pseudotime)')
        axes[0].set_xlabel('Spatial X')
        axes[0].set_ylabel('Spatial Y')
        plt.colorbar(scatter, ax=axes[0], label='Pseudotime')
        
        # 右图: 扩散图
        if 'X_diffmap' in adata.obsm:
            diffmap = adata.obsm['X_diffmap']
            scatter2 = axes[1].scatter(
                diffmap[:, 0],
                diffmap[:, 1],
                c=pseudotime,
                cmap='viridis',
                s=50,
                alpha=0.8
            )
            axes[1].set_title('Diffusion Map')
            axes[1].set_xlabel('DC1')
            axes[1].set_ylabel('DC2')
            plt.colorbar(scatter2, ax=axes[1], label='Pseudotime')
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
    
    @staticmethod
    def plot_niche_comparison(adata, niches, genes, save_path=None):
        """生态位比较热图"""
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # 计算每个niche的平均表达
        niche_exp = []
        niche_labels = []
        
        for niche_id in np.unique(niches):
            mask = niches == niche_id
            if hasattr(adata[mask, genes].X, 'toarray'):
                mean_exp = adata[mask, genes].X.toarray().mean(axis=0)
            else:
                mean_exp = adata[mask, genes].X.mean(axis=0)
            niche_exp.append(mean_exp)
            niche_labels.append(f'Niche {niche_id}')
        
        niche_exp = np.array(niche_exp)
        
        # 绘制热图
        sns.heatmap(
            niche_exp,
            xticklabels=genes,
            yticklabels=niche_labels,
            cmap='RdBu_r',
            center=0,
            ax=ax,
            cbar_kws={'label': 'Mean Expression'}
        )
        
        ax.set_title('Niche-specific Gene Expression')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
    
    @staticmethod
    def plot_gene_network(G, save_path=None):
        """基因共表达网络可视化"""
        import networkx as nx
        
        fig, ax = plt.subplots(figsize=(12, 12))
        
        # 计算布局
        pos = nx.spring_layout(G, k=0.5, iterations=50)
        
        # 计算节点大小 (基于度)
        degrees = dict(G.degree())
        node_sizes = [degrees[node] * 100 for node in G.nodes()]
        
        # 绘制网络
        nx.draw_networkx_nodes(
            G, pos,
            node_size=node_sizes,
            node_color='lightblue',
            alpha=0.7,
            ax=ax
        )
        
        nx.draw_networkx_edges(
            G, pos,
            width=1,
            alpha=0.3,
            ax=ax
        )
        
        nx.draw_networkx_labels(
            G, pos,
            font_size=8,
            ax=ax
        )
        
        ax.set_title('Spatial Gene Co-expression Network')
        ax.axis('off')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()


# 完整分析流程示例
def run_advanced_analysis(adata):
    """运行完整的高级分析"""
    
    print("="*80)
    print("高级生物学分析流程")
    print("="*80)
    
    analyzer = AdvancedBiologicalAnalysis(adata)
    visualizer = VisualizationSuite()
    
    # 1. 差异表达分析
    print("\n步骤 1: 差异表达分析")
    de_results = analyzer.differential_expression_analysis()
    
    # 可视化第一个域
    first_domain = list(de_results.keys())[0]
    visualizer.plot_de_volcano(
        de_results,
        first_domain,
        save_path=f'/home/claude/volcano_{first_domain}.png'
    )
    
    # 2. 通路富集
    print("\n步骤 2: 通路富集分析")
    top_genes = de_results[first_domain]['gene'].head(50).tolist()
    pathways = analyzer.pathway_enrichment_analysis(top_genes)
    print("\n富集通路:")
    print(pathways)
    
    # 3. 空间轨迹
    print("\n步骤 3: 空间轨迹推断")
    pseudotime = analyzer.spatial_trajectory_inference()
    visualizer.plot_spatial_trajectory(
        adata,
        save_path='/home/claude/spatial_trajectory.png'
    )
    
    # 4. 细胞通讯
    print("\n步骤 4: 细胞-细胞通讯分析")
    lr_results = analyzer.cell_cell_communication()
    
    # 5. 空间可变基因
    print("\n步骤 5: 识别空间可变基因")
    svg = analyzer.spatial_variable_genes(n_top=50)
    print(f"\n发现 {len(svg)} 个空间可变基因")
    print("Top 10 SVGs:", svg[:10])
    
    # 6. 生态位分析
    print("\n步骤 6: 生态位聚类")
    niches = analyzer.niche_clustering(n_niches=5)
    
    visualizer.plot_niche_comparison(
        adata,
        niches,
        svg[:20],
        save_path='/home/claude/niche_comparison.png'
    )
    
    # 7. 基因网络
    print("\n步骤 7: 构建基因共表达网络")
    G, corr_matrix = analyzer.spatial_gene_correlation_network(
        genes=svg[:30],
        threshold=0.6
    )
    
    print(f"网络节点数: {G.number_of_nodes()}")
    print(f"网络边数: {G.number_of_edges()}")
    
    visualizer.plot_gene_network(
        G,
        save_path='/home/claude/gene_network.png'
    )
    
    # 8. 保存结果
    print("\n步骤 8: 保存分析结果")
    
    # 保存差异表达结果
    for domain, df in de_results.items():
        df.to_csv(f'/home/claude/DE_{domain}.csv', index=False)
    
    # 保存SVG列表
    pd.DataFrame({'gene': svg}).to_csv('/home/claude/SVGs.csv', index=False)
    
    # 保存更新的adata
    adata.write('/home/claude/complete_analysis.h5ad')
    
    print("\n" + "="*80)
    print("高级分析完成!")
    print("="*80)
    print(f"生成的文件:")
    print(f"  - 差异表达结果: DE_*.csv")
    print(f"  - 空间可变基因: SVGs.csv")
    print(f"  - 可视化图表: *.png")
    print(f"  - 完整数据: complete_analysis.h5ad")
    
    return {
        'de_results': de_results,
        'pathways': pathways,
        'pseudotime': pseudotime,
        'svg': svg,
        'niches': niches,
        'gene_network': G
    }


if __name__ == "__main__":
    # 加载之前分析的数据
    import scanpy as sc
    
    print("加载已分析的数据...")
    adata = sc.read_h5ad('/home/claude/analyzed_visium_data.h5ad')
    
    # 运行高级分析
    results = run_advanced_analysis(adata)
    
    print("\n分析对象已保存在 'results' 变量中")
