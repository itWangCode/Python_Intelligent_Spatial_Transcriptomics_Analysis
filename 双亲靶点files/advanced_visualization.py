#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
SCI-Quality Visualization Module
=============================================================================
Generate publication-ready figures for scientific papers

Features:
- High-resolution (300 DPI) output
- Professional color schemes
- Comprehensive metrics visualization
- All text in English
- Journal-compliant formatting
=============================================================================
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional
from scipy import stats
import logging

logger = logging.getLogger('AdvancedMTL')

# =============================================================================
# Configuration
# =============================================================================

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titlesize'] = 13

# Color palette
COLORS = {
    'primary': '#2E86AB',      # Blue
    'secondary': '#A23B72',    # Purple
    'accent1': '#F18F01',      # Orange
    'accent2': '#C73E1D',      # Red
    'success': '#06A77D',      # Green
    'neutral': '#6C757D'       # Gray
}

# =============================================================================
# Visualization Functions
# =============================================================================

def save_figure(fig, filename: str, figure_dir: str = 'figures', 
                dpi: int = 300):
    """
    Save figure in publication quality
    
    Args:
        fig: Matplotlib figure object
        filename: Output filename
        figure_dir: Directory to save figures
        dpi: Resolution (300 DPI for journals)
    """
    try:
        Path(figure_dir).mkdir(parents=True, exist_ok=True)
        filepath = Path(figure_dir) / filename
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        logger.info(f"Figure saved: {filepath}")
        plt.close(fig)
    except Exception as e:
        logger.error(f"Error saving figure {filename}: {e}")


def plot_training_curves(history: Dict, figure_dir: str = 'figures'):
    """
    Plot training and validation curves
    
    Figure 1: Training Dynamics (4 subplots)
    """
    try:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Training Dynamics and Learning Progress', 
                    fontsize=14, fontweight='bold', y=0.995)
        
        epochs = range(1, len(history['train_loss']) + 1)
        
        # Subplot 1: Overall Loss
        ax = axes[0, 0]
        ax.plot(epochs, history['train_loss'], 
               label='Training', color=COLORS['primary'], 
               linewidth=2, alpha=0.8)
        ax.plot(epochs, history['val_loss'], 
               label='Validation', color=COLORS['secondary'], 
               linewidth=2, alpha=0.8)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Total Loss')
        ax.set_title('(A) Overall Training Loss')
        ax.legend(frameon=True, fancybox=True)
        ax.grid(True, alpha=0.3)
        
        # Subplot 2: Task 1 Loss
        ax = axes[0, 1]
        ax.plot(epochs, history['train_task1_loss'], 
               label='Training', color=COLORS['primary'], 
               linewidth=2, alpha=0.8)
        ax.plot(epochs, history['val_task1_loss'], 
               label='Validation', color=COLORS['secondary'], 
               linewidth=2, alpha=0.8)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Task 1 Loss (MSE)')
        ax.set_title('(B) Target 1 Affinity Prediction')
        ax.legend(frameon=True, fancybox=True)
        ax.grid(True, alpha=0.3)
        
        # Subplot 3: Task 2 Loss
        ax = axes[1, 0]
        ax.plot(epochs, history['train_task2_loss'], 
               label='Training', color=COLORS['primary'], 
               linewidth=2, alpha=0.8)
        ax.plot(epochs, history['val_task2_loss'], 
               label='Validation', color=COLORS['secondary'], 
               linewidth=2, alpha=0.8)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Task 2 Loss (MSE)')
        ax.set_title('(C) Target 2 Affinity Prediction')
        ax.legend(frameon=True, fancybox=True)
        ax.grid(True, alpha=0.3)
        
        # Subplot 4: Learning Rate
        ax = axes[1, 1]
        ax.plot(epochs, history['learning_rate'], 
               color=COLORS['accent1'], linewidth=2, alpha=0.8)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Learning Rate')
        ax.set_title('(D) Learning Rate Schedule')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        save_figure(fig, '01_training_curves.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting training curves: {e}")


def plot_prediction_scatter(metrics: Dict, figure_dir: str = 'figures'):
    """
    Plot predicted vs actual values for both tasks
    
    Figure 2: Prediction Performance
    """
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('Prediction Performance on Test Set', 
                    fontsize=14, fontweight='bold')
        
        predictions = metrics['predictions']
        
        # Task 1
        ax = axes[0]
        true_task1 = predictions['task1_true']
        pred_task1 = predictions['task1_pred']
        
        ax.scatter(true_task1, pred_task1, alpha=0.5, s=30, 
                  color=COLORS['primary'], edgecolors='white', linewidth=0.5)
        
        # Perfect prediction line
        min_val = min(true_task1.min(), pred_task1.min())
        max_val = max(true_task1.max(), pred_task1.max())
        ax.plot([min_val, max_val], [min_val, max_val], 
               'r--', linewidth=2, alpha=0.7, label='Perfect Prediction')
        
        # Linear regression line
        slope, intercept, r_value, _, _ = stats.linregress(true_task1, pred_task1)
        line_x = np.array([min_val, max_val])
        line_y = slope * line_x + intercept
        ax.plot(line_x, line_y, 'b-', linewidth=2, alpha=0.7, 
               label=f'Linear Fit (R²={r_value**2:.3f})')
        
        ax.set_xlabel('True Affinity (kcal/mol)')
        ax.set_ylabel('Predicted Affinity (kcal/mol)')
        ax.set_title('(A) Target 1 Prediction')
        
        # Add metrics text box
        textstr = '\n'.join([
            f"RMSE: {metrics['task1_rmse']:.3f}",
            f"MAE: {metrics['task1_mae']:.3f}",
            f"R²: {metrics['task1_r2']:.3f}",
            f"PCC: {metrics['task1_pcc']:.3f}"
        ])
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, 
               fontsize=9, verticalalignment='top', bbox=props)
        
        ax.legend(loc='lower right')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        # Task 2
        ax = axes[1]
        true_task2 = predictions['task2_true']
        pred_task2 = predictions['task2_pred']
        
        ax.scatter(true_task2, pred_task2, alpha=0.5, s=30, 
                  color=COLORS['secondary'], edgecolors='white', linewidth=0.5)
        
        # Perfect prediction line
        min_val = min(true_task2.min(), pred_task2.min())
        max_val = max(true_task2.max(), pred_task2.max())
        ax.plot([min_val, max_val], [min_val, max_val], 
               'r--', linewidth=2, alpha=0.7, label='Perfect Prediction')
        
        # Linear regression line
        slope, intercept, r_value, _, _ = stats.linregress(true_task2, pred_task2)
        line_x = np.array([min_val, max_val])
        line_y = slope * line_x + intercept
        ax.plot(line_x, line_y, 'b-', linewidth=2, alpha=0.7, 
               label=f'Linear Fit (R²={r_value**2:.3f})')
        
        ax.set_xlabel('True Affinity (kcal/mol)')
        ax.set_ylabel('Predicted Affinity (kcal/mol)')
        ax.set_title('(B) Target 2 Prediction')
        
        # Add metrics text box
        textstr = '\n'.join([
            f"RMSE: {metrics['task2_rmse']:.3f}",
            f"MAE: {metrics['task2_mae']:.3f}",
            f"R²: {metrics['task2_r2']:.3f}",
            f"PCC: {metrics['task2_pcc']:.3f}"
        ])
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, 
               fontsize=9, verticalalignment='top', bbox=props)
        
        ax.legend(loc='lower right')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        plt.tight_layout()
        save_figure(fig, '02_prediction_scatter.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting prediction scatter: {e}")


def plot_residual_analysis(metrics: Dict, figure_dir: str = 'figures'):
    """
    Plot residual analysis for both tasks
    
    Figure 3: Residual Analysis
    """
    try:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        fig.suptitle('Residual Analysis and Error Distribution', 
                    fontsize=14, fontweight='bold', y=0.995)
        
        predictions = metrics['predictions']
        
        # Calculate residuals
        residuals_task1 = predictions['task1_true'] - predictions['task1_pred']
        residuals_task2 = predictions['task2_true'] - predictions['task2_pred']
        
        # Task 1 - Residual Histogram
        ax = axes[0, 0]
        ax.hist(residuals_task1, bins=50, alpha=0.7, color=COLORS['primary'], 
               edgecolor='black', linewidth=0.5)
        ax.axvline(0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.set_xlabel('Residuals (kcal/mol)')
        ax.set_ylabel('Frequency')
        ax.set_title('(A) Target 1 Residual Distribution')
        
        # Add normal distribution overlay
        mu, sigma = residuals_task1.mean(), residuals_task1.std()
        x = np.linspace(residuals_task1.min(), residuals_task1.max(), 100)
        ax2 = ax.twinx()
        ax2.plot(x, stats.norm.pdf(x, mu, sigma) * len(residuals_task1) * 
                (residuals_task1.max() - residuals_task1.min()) / 50,
                'r-', linewidth=2, label='Normal Fit')
        ax2.set_ylabel('Probability Density', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
        textstr = f"μ = {mu:.3f}\nσ = {sigma:.3f}"
        ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax.grid(True, alpha=0.3)
        
        # Task 2 - Residual Histogram
        ax = axes[0, 1]
        ax.hist(residuals_task2, bins=50, alpha=0.7, color=COLORS['secondary'], 
               edgecolor='black', linewidth=0.5)
        ax.axvline(0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.set_xlabel('Residuals (kcal/mol)')
        ax.set_ylabel('Frequency')
        ax.set_title('(B) Target 2 Residual Distribution')
        
        # Add normal distribution overlay
        mu, sigma = residuals_task2.mean(), residuals_task2.std()
        x = np.linspace(residuals_task2.min(), residuals_task2.max(), 100)
        ax2 = ax.twinx()
        ax2.plot(x, stats.norm.pdf(x, mu, sigma) * len(residuals_task2) * 
                (residuals_task2.max() - residuals_task2.min()) / 50,
                'r-', linewidth=2, label='Normal Fit')
        ax2.set_ylabel('Probability Density', color='r')
        ax2.tick_params(axis='y', labelcolor='r')
        
        textstr = f"μ = {mu:.3f}\nσ = {sigma:.3f}"
        ax.text(0.95, 0.95, textstr, transform=ax.transAxes,
               fontsize=9, verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax.grid(True, alpha=0.3)
        
        # Task 1 - Residuals vs Predicted
        ax = axes[1, 0]
        ax.scatter(predictions['task1_pred'], residuals_task1, 
                  alpha=0.5, s=30, color=COLORS['primary'], 
                  edgecolors='white', linewidth=0.5)
        ax.axhline(0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.set_xlabel('Predicted Affinity (kcal/mol)')
        ax.set_ylabel('Residuals (kcal/mol)')
        ax.set_title('(C) Target 1 Residual vs Predicted')
        ax.grid(True, alpha=0.3)
        
        # Add trend line
        z = np.polyfit(predictions['task1_pred'], residuals_task1, 1)
        p = np.poly1d(z)
        x_line = np.linspace(predictions['task1_pred'].min(), 
                           predictions['task1_pred'].max(), 100)
        ax.plot(x_line, p(x_line), "b--", linewidth=2, alpha=0.7, 
               label=f'Trend (slope={z[0]:.4f})')
        ax.legend()
        
        # Task 2 - Residuals vs Predicted
        ax = axes[1, 1]
        ax.scatter(predictions['task2_pred'], residuals_task2, 
                  alpha=0.5, s=30, color=COLORS['secondary'], 
                  edgecolors='white', linewidth=0.5)
        ax.axhline(0, color='red', linestyle='--', linewidth=2, alpha=0.7)
        ax.set_xlabel('Predicted Affinity (kcal/mol)')
        ax.set_ylabel('Residuals (kcal/mol)')
        ax.set_title('(D) Target 2 Residual vs Predicted')
        ax.grid(True, alpha=0.3)
        
        # Add trend line
        z = np.polyfit(predictions['task2_pred'], residuals_task2, 1)
        p = np.poly1d(z)
        x_line = np.linspace(predictions['task2_pred'].min(), 
                           predictions['task2_pred'].max(), 100)
        ax.plot(x_line, p(x_line), "b--", linewidth=2, alpha=0.7, 
               label=f'Trend (slope={z[0]:.4f})')
        ax.legend()
        
        plt.tight_layout()
        save_figure(fig, '03_residual_analysis.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting residual analysis: {e}")


def plot_metrics_comparison(metrics: Dict, figure_dir: str = 'figures'):
    """
    Plot comprehensive metrics comparison
    
    Figure 4: Metrics Comparison
    """
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        fig.suptitle('Performance Metrics Comparison', 
                    fontsize=14, fontweight='bold')
        
        # Prepare data
        metrics_names = ['RMSE', 'MAE', 'R²', 'PCC', 'Spearman']
        task1_values = [
            metrics['task1_rmse'],
            metrics['task1_mae'],
            metrics['task1_r2'],
            metrics['task1_pcc'],
            metrics['task1_spearman']
        ]
        task2_values = [
            metrics['task2_rmse'],
            metrics['task2_mae'],
            metrics['task2_r2'],
            metrics['task2_pcc'],
            metrics['task2_spearman']
        ]
        
        x = np.arange(len(metrics_names))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, task1_values, width, 
                      label='Target 1', color=COLORS['primary'], 
                      alpha=0.8, edgecolor='black', linewidth=0.5)
        bars2 = ax.bar(x + width/2, task2_values, width, 
                      label='Target 2', color=COLORS['secondary'], 
                      alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Add value labels on bars
        def autolabel(bars):
            for bar in bars:
                height = bar.get_height()
                ax.annotate(f'{height:.3f}',
                          xy=(bar.get_x() + bar.get_width() / 2, height),
                          xytext=(0, 3),
                          textcoords="offset points",
                          ha='center', va='bottom', fontsize=8)
        
        autolabel(bars1)
        autolabel(bars2)
        
        ax.set_xlabel('Metrics')
        ax.set_ylabel('Value')
        ax.set_title('Quantitative Performance Comparison')
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_names)
        ax.legend(frameon=True, fancybox=True)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, '04_metrics_comparison.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting metrics comparison: {e}")


def plot_error_distribution(metrics: Dict, figure_dir: str = 'figures'):
    """
    Plot error distribution analysis
    
    Figure 5: Error Distribution
    """
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('Error Distribution Analysis', 
                    fontsize=14, fontweight='bold')
        
        predictions = metrics['predictions']
        
        # Calculate errors
        abs_errors_task1 = np.abs(predictions['task1_true'] - predictions['task1_pred'])
        abs_errors_task2 = np.abs(predictions['task2_true'] - predictions['task2_pred'])
        
        rel_errors_task1 = np.abs((predictions['task1_true'] - predictions['task1_pred']) / 
                                  (predictions['task1_true'] + 1e-8)) * 100
        rel_errors_task2 = np.abs((predictions['task2_true'] - predictions['task2_pred']) / 
                                  (predictions['task2_true'] + 1e-8)) * 100
        
        # Absolute Error Box Plot
        ax = axes[0]
        data = [abs_errors_task1, abs_errors_task2]
        bp = ax.boxplot(data, labels=['Target 1', 'Target 2'],
                       patch_artist=True, widths=0.6,
                       boxprops=dict(facecolor='lightblue', alpha=0.7),
                       medianprops=dict(color='red', linewidth=2),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5))
        
        # Color boxes
        colors = [COLORS['primary'], COLORS['secondary']]
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Absolute Error (kcal/mol)')
        ax.set_title('(A) Absolute Prediction Error')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add mean markers
        means = [abs_errors_task1.mean(), abs_errors_task2.mean()]
        ax.plot([1, 2], means, 'D', color='green', markersize=8, 
               label='Mean', zorder=3)
        ax.legend()
        
        # Relative Error Box Plot
        ax = axes[1]
        # Cap relative errors at 100% for visualization
        rel_errors_task1_capped = np.minimum(rel_errors_task1, 100)
        rel_errors_task2_capped = np.minimum(rel_errors_task2, 100)
        
        data = [rel_errors_task1_capped, rel_errors_task2_capped]
        bp = ax.boxplot(data, labels=['Target 1', 'Target 2'],
                       patch_artist=True, widths=0.6,
                       boxprops=dict(facecolor='lightblue', alpha=0.7),
                       medianprops=dict(color='red', linewidth=2),
                       whiskerprops=dict(linewidth=1.5),
                       capprops=dict(linewidth=1.5))
        
        # Color boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Relative Error (%)')
        ax.set_title('(B) Relative Prediction Error')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Add mean markers
        means = [rel_errors_task1_capped.mean(), rel_errors_task2_capped.mean()]
        ax.plot([1, 2], means, 'D', color='green', markersize=8, 
               label='Mean', zorder=3)
        ax.legend()
        
        plt.tight_layout()
        save_figure(fig, '05_error_distribution.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting error distribution: {e}")


def plot_qq_plots(metrics: Dict, figure_dir: str = 'figures'):
    """
    Plot Q-Q plots to check normality of residuals
    
    Figure 6: Q-Q Plots
    """
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('Quantile-Quantile Plots (Normality Check)', 
                    fontsize=14, fontweight='bold')
        
        predictions = metrics['predictions']
        
        # Calculate residuals
        residuals_task1 = predictions['task1_true'] - predictions['task1_pred']
        residuals_task2 = predictions['task2_true'] - predictions['task2_pred']
        
        # Task 1 Q-Q Plot
        ax = axes[0]
        stats.probplot(residuals_task1, dist="norm", plot=ax)
        ax.get_lines()[0].set_markerfacecolor(COLORS['primary'])
        ax.get_lines()[0].set_markeredgecolor('white')
        ax.get_lines()[0].set_markersize(6)
        ax.get_lines()[0].set_alpha(0.7)
        ax.get_lines()[1].set_color('red')
        ax.get_lines()[1].set_linewidth(2)
        ax.set_title('(A) Target 1 Residuals Q-Q Plot')
        ax.set_xlabel('Theoretical Quantiles')
        ax.set_ylabel('Sample Quantiles')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        shapiro_stat, shapiro_p = stats.shapiro(residuals_task1)
        textstr = f'Shapiro-Wilk Test:\nW = {shapiro_stat:.4f}\np-value = {shapiro_p:.4f}'
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Task 2 Q-Q Plot
        ax = axes[1]
        stats.probplot(residuals_task2, dist="norm", plot=ax)
        ax.get_lines()[0].set_markerfacecolor(COLORS['secondary'])
        ax.get_lines()[0].set_markeredgecolor('white')
        ax.get_lines()[0].set_markersize(6)
        ax.get_lines()[0].set_alpha(0.7)
        ax.get_lines()[1].set_color('red')
        ax.get_lines()[1].set_linewidth(2)
        ax.set_title('(B) Target 2 Residuals Q-Q Plot')
        ax.set_xlabel('Theoretical Quantiles')
        ax.set_ylabel('Sample Quantiles')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        shapiro_stat, shapiro_p = stats.shapiro(residuals_task2)
        textstr = f'Shapiro-Wilk Test:\nW = {shapiro_stat:.4f}\np-value = {shapiro_p:.4f}'
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
               fontsize=9, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        save_figure(fig, '06_qq_plots.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting Q-Q plots: {e}")


def plot_comprehensive_summary(history: Dict, metrics: Dict, 
                              figure_dir: str = 'figures'):
    """
    Plot comprehensive summary figure
    
    Figure 7: Comprehensive Summary (Most Important for Papers)
    """
    try:
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        fig.suptitle('Comprehensive Model Performance Summary', 
                    fontsize=16, fontweight='bold', y=0.995)
        
        predictions = metrics['predictions']
        
        # 1. Training Loss (Top Left)
        ax = fig.add_subplot(gs[0, 0])
        epochs = range(1, len(history['train_loss']) + 1)
        ax.plot(epochs, history['train_loss'], label='Training', 
               color=COLORS['primary'], linewidth=2)
        ax.plot(epochs, history['val_loss'], label='Validation', 
               color=COLORS['secondary'], linewidth=2)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Loss')
        ax.set_title('(A) Learning Curves')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 2. Task 1 Prediction (Top Middle)
        ax = fig.add_subplot(gs[0, 1])
        ax.scatter(predictions['task1_true'], predictions['task1_pred'], 
                  alpha=0.5, s=20, color=COLORS['primary'])
        min_val = min(predictions['task1_true'].min(), predictions['task1_pred'].min())
        max_val = max(predictions['task1_true'].max(), predictions['task1_pred'].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2)
        ax.set_xlabel('True Affinity')
        ax.set_ylabel('Predicted Affinity')
        ax.set_title(f'(B) Target 1 (R²={metrics["task1_r2"]:.3f})')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        # 3. Task 2 Prediction (Top Right)
        ax = fig.add_subplot(gs[0, 2])
        ax.scatter(predictions['task2_true'], predictions['task2_pred'], 
                  alpha=0.5, s=20, color=COLORS['secondary'])
        min_val = min(predictions['task2_true'].min(), predictions['task2_pred'].min())
        max_val = max(predictions['task2_true'].max(), predictions['task2_pred'].max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2)
        ax.set_xlabel('True Affinity')
        ax.set_ylabel('Predicted Affinity')
        ax.set_title(f'(C) Target 2 (R²={metrics["task2_r2"]:.3f})')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal', adjustable='box')
        
        # 4. Task 1 Residuals (Middle Left)
        ax = fig.add_subplot(gs[1, 0])
        residuals_task1 = predictions['task1_true'] - predictions['task1_pred']
        ax.hist(residuals_task1, bins=40, alpha=0.7, color=COLORS['primary'], 
               edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel('Residuals (kcal/mol)')
        ax.set_ylabel('Frequency')
        ax.set_title(f'(D) Target 1 Residuals (σ={residuals_task1.std():.3f})')
        ax.grid(True, alpha=0.3)
        
        # 5. Task 2 Residuals (Middle Middle)
        ax = fig.add_subplot(gs[1, 1])
        residuals_task2 = predictions['task2_true'] - predictions['task2_pred']
        ax.hist(residuals_task2, bins=40, alpha=0.7, color=COLORS['secondary'], 
               edgecolor='black')
        ax.axvline(0, color='red', linestyle='--', linewidth=2)
        ax.set_xlabel('Residuals (kcal/mol)')
        ax.set_ylabel('Frequency')
        ax.set_title(f'(E) Target 2 Residuals (σ={residuals_task2.std():.3f})')
        ax.grid(True, alpha=0.3)
        
        # 6. Metrics Comparison (Middle Right)
        ax = fig.add_subplot(gs[1, 2])
        metrics_names = ['RMSE', 'MAE', 'R²', 'PCC']
        task1_values = [metrics['task1_rmse'], metrics['task1_mae'], 
                       metrics['task1_r2'], metrics['task1_pcc']]
        task2_values = [metrics['task2_rmse'], metrics['task2_mae'], 
                       metrics['task2_r2'], metrics['task2_pcc']]
        x = np.arange(len(metrics_names))
        width = 0.35
        ax.bar(x - width/2, task1_values, width, label='Target 1', 
              color=COLORS['primary'], alpha=0.8)
        ax.bar(x + width/2, task2_values, width, label='Target 2', 
              color=COLORS['secondary'], alpha=0.8)
        ax.set_ylabel('Value')
        ax.set_title('(F) Performance Metrics')
        ax.set_xticks(x)
        ax.set_xticklabels(metrics_names, rotation=45)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        # 7. Learning Rate Schedule (Bottom Left)
        ax = fig.add_subplot(gs[2, 0])
        ax.plot(epochs, history['learning_rate'], color=COLORS['accent1'], linewidth=2)
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Learning Rate')
        ax.set_title('(G) Learning Rate Schedule')
        ax.set_yscale('log')
        ax.grid(True, alpha=0.3)
        
        # 8. Error Distribution (Bottom Middle and Right)
        ax = fig.add_subplot(gs[2, 1:])
        abs_errors_task1 = np.abs(residuals_task1)
        abs_errors_task2 = np.abs(residuals_task2)
        data = [abs_errors_task1, abs_errors_task2]
        bp = ax.boxplot(data, labels=['Target 1', 'Target 2'],
                       patch_artist=True, widths=0.5)
        colors = [COLORS['primary'], COLORS['secondary']]
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_ylabel('Absolute Error (kcal/mol)')
        ax.set_title('(H) Error Distribution Comparison')
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, '07_comprehensive_summary.png', figure_dir)
        
    except Exception as e:
        logger.error(f"Error plotting comprehensive summary: {e}")


def generate_all_plots(history: Dict, metrics: Dict, 
                      figure_dir: str = 'figures'):
    """
    Generate all publication-quality plots
    
    Args:
        history: Training history dictionary
        metrics: Evaluation metrics dictionary
        figure_dir: Directory to save figures
    """
    logger.info("\n" + "="*70)
    logger.info("Generating SCI-Quality Figures")
    logger.info("="*70)
    
    try:
        # Create figure directory
        Path(figure_dir).mkdir(parents=True, exist_ok=True)
        
        # Generate all plots
        logger.info("1/7 Plotting training curves...")
        plot_training_curves(history, figure_dir)
        
        logger.info("2/7 Plotting prediction scatter...")
        plot_prediction_scatter(metrics, figure_dir)
        
        logger.info("3/7 Plotting residual analysis...")
        plot_residual_analysis(metrics, figure_dir)
        
        logger.info("4/7 Plotting metrics comparison...")
        plot_metrics_comparison(metrics, figure_dir)
        
        logger.info("5/7 Plotting error distribution...")
        plot_error_distribution(metrics, figure_dir)
        
        logger.info("6/7 Plotting Q-Q plots...")
        plot_qq_plots(metrics, figure_dir)
        
        logger.info("7/7 Plotting comprehensive summary...")
        plot_comprehensive_summary(history, metrics, figure_dir)
        
        logger.info("\n" + "="*70)
        logger.info("All figures generated successfully!")
        logger.info(f"Figures saved in: {Path(figure_dir).absolute()}")
        logger.info("="*70)
        
    except Exception as e:
        logger.error(f"Error generating plots: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")


# =============================================================================
# End of advanced_visualization.py
# =============================================================================
