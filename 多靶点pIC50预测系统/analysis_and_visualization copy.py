"""
完整数据分析和SCI级可视化
生成论文级别的图表和统计分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import json
from scipy import stats
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

from complete_pic50_predictor import (
    MultiTargetPIC50Predictor,
    Config
)

# 设置SCI论文级别的绘图风格
plt.style.use('seaborn-v0_8-paper')
sns.set_palette("husl")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 13
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'


class ComprehensiveAnalysis:
    """完整的数据分析和可视化"""
    
    def __init__(self, model_path: Path, output_dir: Path = None):
        self.model_path = model_path
        self.output_dir = output_dir or Config.OUTPUT_DIR
        self.figures_dir = self.output_dir / 'figures'
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
        # 加载模型
        print("加载模型...")
        self.predictor = MultiTargetPIC50Predictor.load(model_path)
        
        # 加载训练报告
        report_path = self.output_dir / f"{self.predictor.target_name}_training_report.json"
        with open(report_path, 'r') as f:
            self.training_report = json.load(f)
        
        print(f"✅ 模型加载成功: {self.predictor.target_name}")
        print(f"   训练样本: {self.training_report.get('n_samples', 'N/A')}")
        print(f"   特征维度: {self.training_report.get('n_features', 'N/A')}")
    
    def load_training_data(self):
        """加载训练数据用于分析"""
        data_file = Config.DATA_DIR / 'training_data.csv'
        if data_file.exists():
            df = pd.read_csv(data_file)
            print(f"✅ 加载训练数据: {len(df)} 样本")
            return df
        else:
            print("⚠️  训练数据文件不存在")
            return None
    
    def generate_all_plots(self):
        """生成所有论文级图表"""
        print("\n" + "="*80)
        print("生成SCI论文级图表")
        print("="*80)
        
        # 1. 模型性能对比图
        print("\n1. 生成模型性能对比图...")
        self.plot_model_comparison()
        
        # 2. 数据分布图
        print("2. 生成数据分布分析图...")
        train_df = self.load_training_data()
        if train_df is not None:
            self.plot_data_distribution(train_df)
        
        # 3. 预测散点图
        print("3. 生成预测结果散点图...")
        self.plot_prediction_scatter(train_df)
        
        # 4. 残差分析图
        print("4. 生成残差分析图...")
        self.plot_residual_analysis(train_df)
        
        # 5. 不确定性分析
        print("5. 生成不确定性分析图...")
        self.plot_uncertainty_analysis(train_df)
        
        # 6. 特征重要性
        print("6. 生成特征重要性图...")
        self.plot_feature_importance()
        
        # 7. 综合仪表板
        print("7. 生成综合分析仪表板...")
        self.plot_comprehensive_dashboard(train_df)
        
        print(f"\n✅ 所有图表已保存到: {self.figures_dir}")
    
    def plot_model_comparison(self):
        """图1: 模型性能对比（多指标柱状图）"""
        metrics = self.training_report['metrics']
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Model Performance Comparison', fontsize=14, fontweight='bold')
        
        models = [m['model'] for m in metrics]
        
        # R² 比较
        train_r2 = [m['train_r2'] for m in metrics]
        test_r2 = [m['test_r2'] for m in metrics]
        
        x = np.arange(len(models))
        width = 0.35
        
        axes[0, 0].bar(x - width/2, train_r2, width, label='Training', alpha=0.8, color='#3498db')
        axes[0, 0].bar(x + width/2, test_r2, width, label='Testing', alpha=0.8, color='#e74c3c')
        axes[0, 0].set_ylabel('R² Score', fontweight='bold')
        axes[0, 0].set_title('(A) Coefficient of Determination (R²)')
        axes[0, 0].set_xticks(x)
        axes[0, 0].set_xticklabels(models, rotation=45, ha='right')
        axes[0, 0].legend()
        axes[0, 0].grid(axis='y', alpha=0.3)
        axes[0, 0].set_ylim([0, 1])
        
        # RMSE 比较
        train_rmse = [m['train_rmse'] for m in metrics]
        test_rmse = [m['test_rmse'] for m in metrics]
        
        axes[0, 1].bar(x - width/2, train_rmse, width, label='Training', alpha=0.8, color='#3498db')
        axes[0, 1].bar(x + width/2, test_rmse, width, label='Testing', alpha=0.8, color='#e74c3c')
        axes[0, 1].set_ylabel('RMSE', fontweight='bold')
        axes[0, 1].set_title('(B) Root Mean Square Error')
        axes[0, 1].set_xticks(x)
        axes[0, 1].set_xticklabels(models, rotation=45, ha='right')
        axes[0, 1].legend()
        axes[0, 1].grid(axis='y', alpha=0.3)
        
        # MAE 比较
        train_mae = [m['train_mae'] for m in metrics]
        test_mae = [m['test_mae'] for m in metrics]
        
        axes[1, 0].bar(x - width/2, train_mae, width, label='Training', alpha=0.8, color='#3498db')
        axes[1, 0].bar(x + width/2, test_mae, width, label='Testing', alpha=0.8, color='#e74c3c')
        axes[1, 0].set_ylabel('MAE', fontweight='bold')
        axes[1, 0].set_title('(C) Mean Absolute Error')
        axes[1, 0].set_xticks(x)
        axes[1, 0].set_xticklabels(models, rotation=45, ha='right')
        axes[1, 0].legend()
        axes[1, 0].grid(axis='y', alpha=0.3)
        
        # 过拟合分析（Train R² - Test R²）
        overfitting = [m['train_r2'] - m['test_r2'] for m in metrics]
        colors = ['#2ecc71' if x < 0.15 else '#f39c12' if x < 0.25 else '#e74c3c' for x in overfitting]
        
        axes[1, 1].bar(x, overfitting, color=colors, alpha=0.8)
        axes[1, 1].set_ylabel('R² Difference', fontweight='bold')
        axes[1, 1].set_title('(D) Overfitting Analysis (Train R² - Test R²)')
        axes[1, 1].set_xticks(x)
        axes[1, 1].set_xticklabels(models, rotation=45, ha='right')
        axes[1, 1].axhline(y=0.15, color='orange', linestyle='--', alpha=0.5, label='Acceptable (0.15)')
        axes[1, 1].axhline(y=0.25, color='red', linestyle='--', alpha=0.5, label='High (0.25)')
        axes[1, 1].legend()
        axes[1, 1].grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'Fig1_Model_Comparison.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig1_Model_Comparison.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_data_distribution(self, df):
        """图2: 数据分布分析"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Data Distribution Analysis', fontsize=14, fontweight='bold')
        
        # pIC50分布直方图
        axes[0, 0].hist(df['pIC50'], bins=50, color='#3498db', alpha=0.7, edgecolor='black')
        axes[0, 0].set_xlabel('pIC50', fontweight='bold')
        axes[0, 0].set_ylabel('Frequency', fontweight='bold')
        axes[0, 0].set_title('(A) pIC50 Distribution')
        axes[0, 0].axvline(df['pIC50'].mean(), color='red', linestyle='--', 
                          label=f'Mean: {df["pIC50"].mean():.2f}')
        axes[0, 0].axvline(df['pIC50'].median(), color='green', linestyle='--', 
                          label=f'Median: {df["pIC50"].median():.2f}')
        axes[0, 0].legend()
        axes[0, 0].grid(axis='y', alpha=0.3)
        
        # 累积分布函数
        sorted_pic50 = np.sort(df['pIC50'])
        cumulative = np.arange(1, len(sorted_pic50) + 1) / len(sorted_pic50)
        axes[0, 1].plot(sorted_pic50, cumulative, color='#3498db', linewidth=2)
        axes[0, 1].set_xlabel('pIC50', fontweight='bold')
        axes[0, 1].set_ylabel('Cumulative Probability', fontweight='bold')
        axes[0, 1].set_title('(B) Cumulative Distribution Function')
        axes[0, 1].grid(alpha=0.3)
        
        # 正态性检验 Q-Q图
        stats.probplot(df['pIC50'], dist="norm", plot=axes[1, 0])
        axes[1, 0].set_title('(C) Q-Q Plot (Normality Test)')
        axes[1, 0].grid(alpha=0.3)
        
        # 箱线图和小提琴图
        parts = axes[1, 1].violinplot([df['pIC50']], positions=[0], widths=0.7,
                                      showmeans=True, showmedians=True)
        axes[1, 1].boxplot([df['pIC50']], positions=[0], widths=0.3)
        axes[1, 1].set_ylabel('pIC50', fontweight='bold')
        axes[1, 1].set_title('(D) Violin Plot with Box Plot')
        axes[1, 1].set_xticks([0])
        axes[1, 1].set_xticklabels(['Dataset'])
        axes[1, 1].grid(axis='y', alpha=0.3)
        
        # 添加统计信息
        stats_text = f"n = {len(df)}\n"
        stats_text += f"Mean = {df['pIC50'].mean():.3f}\n"
        stats_text += f"Std = {df['pIC50'].std():.3f}\n"
        stats_text += f"Min = {df['pIC50'].min():.3f}\n"
        stats_text += f"Max = {df['pIC50'].max():.3f}\n"
        stats_text += f"Skewness = {stats.skew(df['pIC50']):.3f}\n"
        stats_text += f"Kurtosis = {stats.kurtosis(df['pIC50']):.3f}"
        
        axes[1, 1].text(0.5, df['pIC50'].max(), stats_text,
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                       verticalalignment='top', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'Fig2_Data_Distribution.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig2_Data_Distribution.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_prediction_scatter(self, df):
        """图3: 预测vs实际值散点图（所有模型）"""
        if df is None:
            return
        
        # 使用ensemble模型预测所有数据
        print("   正在生成预测...")
        smiles_list = df['smiles'].tolist()
        
        # 分批预测以避免内存问题
        batch_size = 500
        all_predictions = {'random_forest': [], 'gradient_boosting': [], 'ensemble': []}
        
        for i in range(0, len(smiles_list), batch_size):
            batch = smiles_list[i:i+batch_size]
            for model_name in ['random_forest', 'gradient_boosting', 'ensemble']:
                results = self.predictor.predict_batch(batch, model_name=model_name)
                all_predictions[model_name].extend([r['pIC50'] for r in results])
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('Predicted vs Actual pIC50 Values', fontsize=14, fontweight='bold')
        
        for idx, (model_name, title) in enumerate([
            ('random_forest', '(A) Random Forest'),
            ('gradient_boosting', '(B) Gradient Boosting'),
            ('ensemble', '(C) Ensemble')
        ]):
            y_pred = all_predictions[model_name]
            y_true = df['pIC50'].values[:len(y_pred)]
            
            # 计算指标
            r2 = r2_score(y_true, y_pred)
            rmse = np.sqrt(mean_squared_error(y_true, y_pred))
            mae = mean_absolute_error(y_true, y_pred)
            
            # 散点图
            axes[idx].scatter(y_true, y_pred, alpha=0.5, s=20, color='#3498db', edgecolors='none')
            
            # 理想线
            min_val = min(y_true.min(), min(y_pred))
            max_val = max(y_true.max(), max(y_pred))
            axes[idx].plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label='Ideal')
            
            # 回归线
            z = np.polyfit(y_true, y_pred, 1)
            p = np.poly1d(z)
            axes[idx].plot(y_true, p(y_true), "g-", linewidth=2, alpha=0.8, label='Fit')
            
            axes[idx].set_xlabel('Actual pIC50', fontweight='bold')
            axes[idx].set_ylabel('Predicted pIC50', fontweight='bold')
            axes[idx].set_title(title)
            axes[idx].legend()
            axes[idx].grid(alpha=0.3)
            
            # 添加统计信息
            text = f'R² = {r2:.3f}\nRMSE = {rmse:.3f}\nMAE = {mae:.3f}\nn = {len(y_true)}'
            axes[idx].text(0.05, 0.95, text, transform=axes[idx].transAxes,
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                          verticalalignment='top', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'Fig3_Prediction_Scatter.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig3_Prediction_Scatter.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 保存预测结果
        prediction_df = df.copy()
        prediction_df['predicted_RF'] = all_predictions['random_forest']
        prediction_df['predicted_GB'] = all_predictions['gradient_boosting']
        prediction_df['predicted_Ensemble'] = all_predictions['ensemble']
        prediction_df.to_csv(self.output_dir / 'detailed_predictions.csv', index=False)
        print(f"   ✅ 详细预测结果已保存")
    
    def plot_residual_analysis(self, df):
        """图4: 残差分析"""
        if df is None:
            return
        
        # 使用ensemble模型
        results = self.predictor.predict_batch(df['smiles'].tolist(), model_name='ensemble')
        y_pred = np.array([r['pIC50'] for r in results])
        y_true = df['pIC50'].values[:len(y_pred)]
        residuals = y_true - y_pred
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Residual Analysis (Ensemble Model)', fontsize=14, fontweight='bold')
        
        # 残差散点图
        axes[0, 0].scatter(y_pred, residuals, alpha=0.5, s=20, color='#3498db', edgecolors='none')
        axes[0, 0].axhline(y=0, color='r', linestyle='--', linewidth=2)
        axes[0, 0].set_xlabel('Predicted pIC50', fontweight='bold')
        axes[0, 0].set_ylabel('Residuals', fontweight='bold')
        axes[0, 0].set_title('(A) Residuals vs Predicted Values')
        axes[0, 0].grid(alpha=0.3)
        
        # 残差直方图
        axes[0, 1].hist(residuals, bins=50, color='#3498db', alpha=0.7, edgecolor='black')
        axes[0, 1].set_xlabel('Residuals', fontweight='bold')
        axes[0, 1].set_ylabel('Frequency', fontweight='bold')
        axes[0, 1].set_title('(B) Residuals Distribution')
        axes[0, 1].axvline(0, color='r', linestyle='--', linewidth=2)
        axes[0, 1].grid(axis='y', alpha=0.3)
        
        # Q-Q图
        stats.probplot(residuals, dist="norm", plot=axes[1, 0])
        axes[1, 0].set_title('(C) Q-Q Plot of Residuals')
        axes[1, 0].grid(alpha=0.3)
        
        # 标准化残差
        std_residuals = residuals / np.std(residuals)
        axes[1, 1].scatter(y_pred, std_residuals, alpha=0.5, s=20, color='#3498db', edgecolors='none')
        axes[1, 1].axhline(y=0, color='r', linestyle='--', linewidth=2)
        axes[1, 1].axhline(y=2, color='orange', linestyle='--', alpha=0.5)
        axes[1, 1].axhline(y=-2, color='orange', linestyle='--', alpha=0.5)
        axes[1, 1].set_xlabel('Predicted pIC50', fontweight='bold')
        axes[1, 1].set_ylabel('Standardized Residuals', fontweight='bold')
        axes[1, 1].set_title('(D) Standardized Residuals')
        axes[1, 1].grid(alpha=0.3)
        
        # 添加统计信息
        outliers = np.sum(np.abs(std_residuals) > 2)
        stats_text = f"Mean = {np.mean(residuals):.3f}\n"
        stats_text += f"Std = {np.std(residuals):.3f}\n"
        stats_text += f"Outliers (|z|>2) = {outliers} ({100*outliers/len(residuals):.1f}%)"
        
        axes[1, 1].text(0.05, 0.95, stats_text, transform=axes[1, 1].transAxes,
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                       verticalalignment='top', fontsize=9)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'Fig4_Residual_Analysis.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig4_Residual_Analysis.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_uncertainty_analysis(self, df):
        """图5: 不确定性分析"""
        if df is None or len(df) > 1000:
            # 如果数据太多，随机采样
            df_sample = df.sample(min(1000, len(df)), random_state=42)
        else:
            df_sample = df
        
        print("   正在计算不确定性...")
        uncertainties = []
        predictions = []
        actuals = []
        
        for idx, row in df_sample.iterrows():
            result = self.predictor.predict(row['smiles'], model_name='ensemble', 
                                           return_uncertainty=True)
            if result and 'uncertainty' in result:
                predictions.append(result['pIC50'])
                actuals.append(row['pIC50'])
                unc = result['uncertainty']
                uncertainties.append(unc.get('std', 0))
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle('Prediction Uncertainty Analysis', fontsize=14, fontweight='bold')
        
        # 不确定性 vs 预测值
        scatter = axes[0].scatter(predictions, uncertainties, c=actuals, 
                                 cmap='viridis', alpha=0.6, s=50, edgecolors='black')
        axes[0].set_xlabel('Predicted pIC50', fontweight='bold')
        axes[0].set_ylabel('Prediction Uncertainty (Std)', fontweight='bold')
        axes[0].set_title('(A) Uncertainty vs Predicted Value')
        axes[0].grid(alpha=0.3)
        cbar = plt.colorbar(scatter, ax=axes[0])
        cbar.set_label('Actual pIC50', fontweight='bold')
        
        # 误差 vs 不确定性
        errors = np.abs(np.array(actuals) - np.array(predictions))
        axes[1].scatter(uncertainties, errors, alpha=0.6, s=50, 
                       color='#3498db', edgecolors='black')
        
        # 相关性分析
        if len(uncertainties) > 0:
            corr, p_value = stats.pearsonr(uncertainties, errors)
            z = np.polyfit(uncertainties, errors, 1)
            p = np.poly1d(z)
            axes[1].plot(uncertainties, p(uncertainties), "r-", linewidth=2, alpha=0.8)
            
            axes[1].text(0.05, 0.95, f'Correlation: {corr:.3f}\np-value: {p_value:.3e}',
                        transform=axes[1].transAxes,
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                        verticalalignment='top', fontsize=9)
        
        axes[1].set_xlabel('Prediction Uncertainty (Std)', fontweight='bold')
        axes[1].set_ylabel('Prediction Error (|Actual - Predicted|)', fontweight='bold')
        axes[1].set_title('(B) Error vs Uncertainty')
        axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'Fig5_Uncertainty_Analysis.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig5_Uncertainty_Analysis.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_feature_importance(self):
        """图6: 特征重要性"""
        if 'random_forest' not in self.predictor.models:
            return
        
        rf_model = self.predictor.models['random_forest']
        importances = rf_model.feature_importances_
        
        # 获取前20个重要特征
        indices = np.argsort(importances)[::-1][:20]
        top_importances = importances[indices]
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('Feature Importance Analysis (Random Forest)', fontsize=14, fontweight='bold')
        
        # 柱状图
        colors = plt.cm.viridis(np.linspace(0, 1, len(top_importances)))
        axes[0].barh(range(len(top_importances)), top_importances, color=colors, alpha=0.8)
        axes[0].set_yticks(range(len(top_importances)))
        axes[0].set_yticklabels([f'Feature {i}' for i in indices])
        axes[0].set_xlabel('Importance Score', fontweight='bold')
        axes[0].set_title('(A) Top 20 Important Features')
        axes[0].invert_yaxis()
        axes[0].grid(axis='x', alpha=0.3)
        
        # 累积重要性
        sorted_importances = np.sort(importances)[::-1]
        cumulative_importances = np.cumsum(sorted_importances)
        
        axes[1].plot(range(1, len(cumulative_importances)+1), cumulative_importances, 
                    linewidth=2, color='#3498db')
        axes[1].axhline(y=0.95, color='r', linestyle='--', label='95% Threshold')
        axes[1].set_xlabel('Number of Features', fontweight='bold')
        axes[1].set_ylabel('Cumulative Importance', fontweight='bold')
        axes[1].set_title('(B) Cumulative Feature Importance')
        axes[1].legend()
        axes[1].grid(alpha=0.3)
        
        # 找到达到95%所需的特征数
        n_features_95 = np.argmax(cumulative_importances >= 0.95) + 1
        axes[1].axvline(x=n_features_95, color='orange', linestyle='--', alpha=0.5)
        axes[1].text(n_features_95, 0.5, f'{n_features_95} features\nfor 95%',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                    fontsize=9)
        
        plt.tight_layout()
        plt.savefig(self.figures_dir / 'Fig6_Feature_Importance.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig6_Feature_Importance.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_comprehensive_dashboard(self, df):
        """图7: 综合分析仪表板"""
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        fig.suptitle('Comprehensive Analysis Dashboard', fontsize=16, fontweight='bold')
        
        # 模型性能雷达图
        ax1 = fig.add_subplot(gs[0, 0], projection='polar')
        metrics = self.training_report['metrics']
        ensemble_metrics = [m for m in metrics if m['model'] == 'Ensemble'][0]
        
        categories = ['R²', 'RMSE\n(inv)', 'MAE\n(inv)', 'Train-Test\nGap (inv)']
        values = [
            ensemble_metrics['test_r2'],
            1 - min(ensemble_metrics['test_rmse'] / 2, 1),  # 归一化
            1 - min(ensemble_metrics['test_mae'] / 2, 1),
            1 - min((ensemble_metrics['train_r2'] - ensemble_metrics['test_r2']) * 4, 1)
        ]
        
        angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
        values += values[:1]
        angles += angles[:1]
        
        ax1.plot(angles, values, 'o-', linewidth=2, color='#3498db')
        ax1.fill(angles, values, alpha=0.25, color='#3498db')
        ax1.set_xticks(angles[:-1])
        ax1.set_xticklabels(categories, fontsize=8)
        ax1.set_ylim(0, 1)
        ax1.set_title('Model Performance\nRadar Chart', fontsize=10, pad=20)
        ax1.grid(True)
        
        # 其他子图...
        # 这里可以继续添加更多分析图表
        
        plt.savefig(self.figures_dir / 'Fig7_Dashboard.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.figures_dir / 'Fig7_Dashboard.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    def generate_summary_report(self):
        """生成文字总结报告"""
        print("\n生成总结报告...")
        
        report_lines = []
        report_lines.append("="*80)
        report_lines.append("COMPREHENSIVE ANALYSIS REPORT")
        report_lines.append("="*80)
        report_lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"Model: {self.predictor.target_name}")
        report_lines.append(f"\n1. DATASET INFORMATION")
        report_lines.append(f"   - Total samples: {self.training_report.get('n_samples', 'N/A')}")
        report_lines.append(f"   - Feature dimension: {self.training_report.get('n_features', 'N/A')}")
        
        report_lines.append(f"\n2. MODEL PERFORMANCE")
        for metric in self.training_report['metrics']:
            report_lines.append(f"\n   {metric['model']}:")
            report_lines.append(f"     - Training R²: {metric['train_r2']:.4f}")
            report_lines.append(f"     - Testing R²: {metric['test_r2']:.4f}")
            report_lines.append(f"     - Testing RMSE: {metric['test_rmse']:.4f}")
            report_lines.append(f"     - Testing MAE: {metric['test_mae']:.4f}")
            overfitting = metric['train_r2'] - metric['test_r2']
            status = "Good" if overfitting < 0.15 else "Moderate" if overfitting < 0.25 else "High"
            report_lines.append(f"     - Overfitting: {overfitting:.4f} ({status})")
        
        report_lines.append(f"\n3. GENERATED FILES")
        report_lines.append(f"   - Figures directory: {self.figures_dir}")
        report_lines.append(f"   - Detailed predictions: {self.output_dir / 'detailed_predictions.csv'}")
        
        report_lines.append(f"\n{'='*80}\n")
        
        report_text = "\n".join(report_lines)
        
        # 保存到文件
        report_file = self.output_dir / 'Analysis_Report.txt'
        with open(report_file, 'w') as f:
            f.write(report_text)
        
        print(report_text)
        print(f"✅ 报告已保存: {report_file}")


def main():
    """主函数"""
    print("""
╔══════════════════════════════════════════════════════════════════╗
║           完整数据分析和SCI级可视化                               ║
║           Comprehensive Analysis & SCI-Level Visualization       ║
╚══════════════════════════════════════════════════════════════════╝
    """)
    
    # 加载模型
    model_path = Config.MODEL_DIR / 'USER_DATA_complete_model.pkl'
    
    if not model_path.exists():
        print(f"❌ 模型文件不存在: {model_path}")
        print("请先运行训练脚本生成模型")
        return
    
    # 创建分析对象
    analysis = ComprehensiveAnalysis(model_path)
    
    # 生成所有图表
    analysis.generate_all_plots()
    
    # 生成总结报告
    analysis.generate_summary_report()
    
    print("\n" + "🎉"*30)
    print("分析完成！")
    print(f"所有图表和结果已保存到: {analysis.output_dir}")
    print("🎉"*30)


if __name__ == "__main__":
    Config.setup_directories()
    main()
