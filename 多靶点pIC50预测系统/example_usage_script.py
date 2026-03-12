"""
真实数据测试 - 使用用户提供的数据文件
训练数据：data_ubstructures_matches.csv
测试数据：predicted_pIC50_df.csv
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

from complete_pic50_predictor import (
    MultiTargetPIC50Predictor,
    ModelConfig,
    DataManager,
    Config
)

def load_and_validate_data():
    """加载并验证数据文件"""
    print("="*80)
    print("数据加载和验证")
    print("="*80)
    
    # 训练数据
    train_file = Path('data_ubstructures_matches.csv')
    if not train_file.exists():
        print(f"❌ 训练数据文件不存在: {train_file}")
        print("请确保 data_ubstructures_matches.csv 在当前目录下")
        return None, None
    
    print(f"\n✅ 找到训练数据: {train_file}")
    train_df = pd.read_csv(train_file)
    print(f"   训练数据形状: {train_df.shape}")
    print(f"   列名: {list(train_df.columns)}")
    
    # 检查必需的列
    if 'smiles' not in train_df.columns:
        print("   ❌ 缺少 'smiles' 列")
        if 'SMILES' in train_df.columns:
            print("   发现 'SMILES' 列，将其重命名为 'smiles'")
            train_df = train_df.rename(columns={'SMILES': 'smiles'})
        else:
            print("   可用列:", list(train_df.columns))
            return None, None
    
    if 'pIC50' not in train_df.columns:
        print("   ❌ 缺少 'pIC50' 列")
        # 尝试查找可能的列名
        possible_cols = [col for col in train_df.columns if 'pic50' in col.lower() or 'activity' in col.lower()]
        if possible_cols:
            print(f"   可能的活性列: {possible_cols}")
            print(f"   请手动指定活性列名，或将列重命名为 'pIC50'")
        return None, None
    
    # 清理数据
    train_df = train_df.dropna(subset=['smiles', 'pIC50'])
    train_df = train_df[train_df['pIC50'].between(3, 12)]  # 过滤异常值
    
    print(f"   清理后样本数: {len(train_df)}")
    print(f"   pIC50 范围: [{train_df['pIC50'].min():.2f}, {train_df['pIC50'].max():.2f}]")
    print(f"   pIC50 均值: {train_df['pIC50'].mean():.2f} ± {train_df['pIC50'].std():.2f}")
    
    # 测试数据
    test_file = Path('predicted_pIC50_df.csv')
    if not test_file.exists():
        print(f"\n⚠️  测试数据文件不存在: {test_file}")
        print("将使用训练数据的一部分作为测试")
        test_df = None
    else:
        print(f"\n✅ 找到测试数据: {test_file}")
        test_df = pd.read_csv(test_file)
        print(f"   测试数据形状: {test_df.shape}")
        print(f"   列名: {list(test_df.columns)}")
        
        if 'smiles' not in test_df.columns:
            if 'SMILES' in test_df.columns:
                test_df = test_df.rename(columns={'SMILES': 'smiles'})
            else:
                print("   ❌ 缺少 'smiles' 列")
                test_df = None
        
        if test_df is not None:
            test_df = test_df.dropna(subset=['smiles'])
            print(f"   清理后样本数: {len(test_df)}")
    
    return train_df, test_df


def train_model_with_real_data(train_df):
    """使用真实数据训练模型"""
    print("\n" + "="*80)
    print("模型训练")
    print("="*80)
    
    # 保存到标准位置
    Config.setup_directories()
    data_file = Config.DATA_DIR / 'training_data.csv'
    train_df[['smiles', 'pIC50']].to_csv(data_file, index=False)
    print(f"\n数据已保存到: {data_file}")
    
    # 配置模型（根据数据量调整）
    n_samples = len(train_df)
    
    if n_samples < 100:
        print(f"⚠️  样本数较少 ({n_samples})，使用快速配置")
        config = ModelConfig(
            use_transformer=False,
            use_gnn=False,
            use_ensemble=True,
            rf_n_estimators=50,
            gb_n_estimators=25
        )
    elif n_samples < 500:
        print(f"✅ 样本数适中 ({n_samples})，使用标准配置")
        config = ModelConfig(
            use_transformer=False,
            use_gnn=False,
            use_ensemble=True,
            rf_n_estimators=100,
            gb_n_estimators=50
        )
    else:
        print(f"✅ 样本数充足 ({n_samples})，使用完整配置")
        config = ModelConfig(
            use_transformer=False,  # 可设为True如果有PyTorch和GPU
            use_gnn=False,
            use_ensemble=True,
            rf_n_estimators=200,
            gb_n_estimators=100
        )
    
    print("\n开始训练...")
    predictor = MultiTargetPIC50Predictor('USER_DATA', config)
    
    smiles_list = train_df['smiles'].tolist()
    pic50_values = train_df['pIC50'].tolist()
    
    try:
        metrics = predictor.train(smiles_list, pic50_values)
        
        print("\n✅ 训练完成！")
        print(f"   可用模型: {list(predictor.models.keys())}")
        
        # 显示性能
        print("\n模型性能:")
        for metric in metrics:
            print(f"\n  {metric['model']}:")
            print(f"    训练集 R²: {metric['train_r2']:.4f}")
            print(f"    测试集 R²: {metric['test_r2']:.4f}")
            print(f"    测试集 RMSE: {metric['test_rmse']:.4f}")
            print(f"    测试集 MAE: {metric['test_mae']:.4f}")
        
        # 保存模型
        predictor.save()
        print(f"\n✅ 模型已保存")
        
        # 保存训练报告
        report_path = Config.OUTPUT_DIR / 'USER_DATA_training_report.json'
        DataManager.save_training_report(predictor, report_path)
        
        return predictor
        
    except Exception as e:
        print(f"\n❌ 训练失败: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_predictions(predictor, test_df):
    """测试预测功能"""
    print("\n" + "="*80)
    print("预测测试")
    print("="*80)
    
    # 测试单个预测
    print("\n1. 单个预测测试")
    test_smiles = test_df['smiles'].iloc[0] if test_df is not None else 'CC(CC1=CC=CC=C1)NC'
    print(f"   测试 SMILES: {test_smiles}")
    
    for model_name in ['random_forest', 'gradient_boosting', 'ensemble']:
        if model_name in predictor.models:
            try:
                result = predictor.predict(
                    test_smiles,
                    model_name=model_name,
                    return_uncertainty=True
                )
                
                if result:
                    print(f"\n   {model_name}:")
                    print(f"     ✅ pIC50: {result['pIC50']:.2f}")
                    if 'uncertainty' in result and result['uncertainty']:
                        unc = result['uncertainty']
                        if 'std' in unc:
                            print(f"     不确定性 (std): {unc['std']:.3f}")
                else:
                    print(f"   {model_name}: ❌ 预测失败")
                    return False
                    
            except Exception as e:
                print(f"   {model_name}: ❌ 错误 - {e}")
                import traceback
                traceback.print_exc()
                return False
    
    # 测试批量预测
    print("\n2. 批量预测测试")
    
    if test_df is not None and len(test_df) > 0:
        # 使用测试数据
        batch_size = min(10, len(test_df))
        batch_smiles = test_df['smiles'].head(batch_size).tolist()
        print(f"   使用测试数据的前 {batch_size} 个样本")
    else:
        # 使用示例数据
        batch_smiles = [
            'CC(CC1=CC=CC=C1)NC',
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
            'COC1=C(C=C2C(=C1)C(=NC(=N2)N)N)OC',
        ]
        print(f"   使用 {len(batch_smiles)} 个示例 SMILES")
    
    try:
        results = predictor.predict_batch(batch_smiles, model_name='ensemble')
        
        if results:
            print(f"\n   ✅ 批量预测成功: {len(results)} 个化合物")
            print("\n   预测结果:")
            for i, r in enumerate(results[:5], 1):  # 显示前5个
                print(f"     {i}. {r['smiles'][:40]}... -> pIC50: {r['pIC50']:.2f}")
            
            if len(results) > 5:
                print(f"     ... (还有 {len(results)-5} 个结果)")
            
            # 保存结果
            output_file = Config.OUTPUT_DIR / 'batch_predictions.csv'
            DataManager.save_results(results, output_file)
            print(f"\n   ✅ 结果已保存: {output_file}")
            
            # 统计信息
            pIC50_values = [r['pIC50'] for r in results]
            print(f"\n   预测统计:")
            print(f"     均值: {np.mean(pIC50_values):.2f}")
            print(f"     标准差: {np.std(pIC50_values):.2f}")
            print(f"     范围: [{np.min(pIC50_values):.2f}, {np.max(pIC50_values):.2f}]")
            
        else:
            print("   ❌ 批量预测失败：返回空结果")
            return False
            
    except Exception as e:
        print(f"   ❌ 批量预测错误: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True


def main():
    """主函数"""
    print("""
╔══════════════════════════════════════════════════════════════════╗
║              真实数据测试 - 完整工作流                            ║
║              Real Data Test - Complete Workflow                  ║
╚══════════════════════════════════════════════════════════════════╝
    """)
    
    # 初始化环境
    Config.setup_directories()
    Config.setup_logging()
    
    # 1. 加载数据
    train_df, test_df = load_and_validate_data()
    
    if train_df is None:
        print("\n❌ 数据加载失败，测试终止")
        print("\n请确保以下文件存在:")
        print("  - data_ubstructures_matches.csv (训练数据)")
        print("  - predicted_pIC50_df.csv (测试数据，可选)")
        print("\n文件必须包含以下列:")
        print("  - smiles 或 SMILES (化合物结构)")
        print("  - pIC50 (活性数据)")
        return False
    
    print(f"\n{'='*80}")
    print(f"数据加载成功！开始训练...")
    print(f"{'='*80}")
    
    # 2. 训练模型
    predictor = train_model_with_real_data(train_df)
    
    if predictor is None:
        print("\n❌ 模型训练失败")
        return False
    
    # 3. 测试预测
    success = test_predictions(predictor, test_df)
    
    if not success:
        print("\n❌ 预测测试失败")
        return False
    
    # 4. 总结
    print("\n" + "="*80)
    print("测试总结")
    print("="*80)
    
    print("\n✅ 所有测试通过！")
    print("\n生成的文件:")
    print(f"  模型: {Config.MODEL_DIR / 'USER_DATA_complete_model.pkl'}")
    print(f"  训练报告: {Config.OUTPUT_DIR / 'USER_DATA_training_report.json'}")
    print(f"  预测结果: {Config.OUTPUT_DIR / 'batch_predictions.csv'}")
    print(f"  日志: {Config.LOG_DIR}")
    
    print("\n" + "🎉"*30)
    print("系统已准备就绪！")
    print("你现在可以使用训练好的模型进行预测了。")
    print("🎉"*30)
    
    print("\n使用方法:")
    print("  1. 单个预测:")
    print("     python complete_pic50_predictor.py predict --target USER_DATA --smiles 'YOUR_SMILES'")
    print("\n  2. 批量预测:")
    print("     python complete_pic50_predictor.py batch --target USER_DATA --input your_file.csv")
    print("\n  3. 在Python中使用:")
    print("     from complete_pic50_predictor import MultiTargetPIC50Predictor")
    print("     predictor = MultiTargetPIC50Predictor.load('models/USER_DATA_complete_model.pkl')")
    print("     result = predictor.predict('CC(CC1=CC=CC=C1)NC', model_name='ensemble')")
    
    return True


if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\n⚠️  用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ 未预期的错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)