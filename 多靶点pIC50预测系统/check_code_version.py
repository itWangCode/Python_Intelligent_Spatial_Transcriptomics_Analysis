"""
检查代码版本 - 验证文件是否已更新
"""

import re
from pathlib import Path

def check_file_version():
    """检查文件是否包含修复"""
    
    file_path = Path('complete_pic50_predictor.py')
    
    if not file_path.exists():
        print("❌ 找不到 complete_pic50_predictor.py")
        return False
    
    print("检查文件版本...")
    print("="*80)
    
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 检查1: predict方法中是否有正确的判断顺序
    print("\n1. 检查 predict() 方法...")
    
    # 找到predict方法
    predict_match = re.search(
        r'def predict\(self.*?\n(.*?)(?=\n    def |\nclass |\Z)',
        content,
        re.DOTALL
    )
    
    if predict_match:
        predict_code = predict_match.group(1)
        
        # 检查是否先判断isinstance(model, dict)
        if 'if isinstance(model, dict):' in predict_code:
            # 检查这个判断是否在其他判断之前
            dict_pos = predict_code.find('if isinstance(model, dict):')
            module_pos = predict_code.find('elif TORCH_AVAILABLE and isinstance(model, nn.Module):')
            else_pos = predict_code.find('else:\n            # Scikit-learn')
            
            if dict_pos < module_pos and dict_pos < else_pos:
                print("   ✅ predict() 方法判断顺序正确")
                print(f"      - isinstance(model, dict) 在位置: {dict_pos}")
            else:
                print("   ❌ predict() 方法判断顺序错误")
                print(f"      - isinstance(model, dict) 在位置: {dict_pos}")
                print(f"      - 应该在所有其他判断之前")
                return False
        else:
            print("   ❌ predict() 方法缺少 isinstance(model, dict) 判断")
            return False
    else:
        print("   ❌ 找不到 predict() 方法")
        return False
    
    # 检查2: predict_batch方法
    print("\n2. 检查 predict_batch() 方法...")
    
    batch_match = re.search(
        r'def predict_batch\(self.*?\n(.*?)(?=\n    def |\nclass |\Z)',
        content,
        re.DOTALL
    )
    
    if batch_match:
        batch_code = batch_match.group(1)
        
        if 'if isinstance(model, dict):' in batch_code:
            print("   ✅ predict_batch() 方法包含正确判断")
        else:
            print("   ❌ predict_batch() 方法缺少正确判断")
            return False
    else:
        print("   ❌ 找不到 predict_batch() 方法")
        return False
    
    # 检查3: _estimate_uncertainty方法
    print("\n3. 检查 _estimate_uncertainty() 方法...")
    
    uncertainty_match = re.search(
        r'def _estimate_uncertainty\(self.*?\n(.*?)(?=\n    def |\nclass |\Z)',
        content,
        re.DOTALL
    )
    
    if uncertainty_match:
        uncertainty_code = uncertainty_match.group(1)
        
        if 'if isinstance(model, dict):' in uncertainty_code:
            print("   ✅ _estimate_uncertainty() 方法包含正确判断")
        else:
            print("   ❌ _estimate_uncertainty() 方法缺少正确判断")
            return False
    else:
        print("   ❌ 找不到 _estimate_uncertainty() 方法")
        return False
    
    # 检查4: 查找旧的错误代码
    print("\n4. 检查是否有旧的错误代码...")
    
    # 搜索第895行附近的代码
    lines = content.split('\n')
    
    if len(lines) > 895:
        line_895 = lines[894]  # 0-indexed
        print(f"   第895行: {line_895.strip()}")
        
        if 'prediction = model.predict([feature_vector])[0]' in line_895:
            print("   ❌ 第895行仍然是旧代码！")
            print("   这行应该不会被执行到（应该被 isinstance(model, dict) 拦截）")
            
            # 检查上下文
            print(f"\n   上下文 (第890-900行):")
            for i in range(max(0, 889), min(len(lines), 900)):
                marker = " >>> " if i == 894 else "     "
                print(f"{marker}{i+1}: {lines[i]}")
            
            return False
    
    print("\n" + "="*80)
    print("✅ 文件版本检查通过！")
    print("="*80)
    
    return True


def show_fix_instructions():
    """显示修复说明"""
    print("\n" + "="*80)
    print("修复说明")
    print("="*80)
    
    print("""
你的文件没有正确更新。请按以下步骤操作：

方法1：手动修改文件
--------------------
在 complete_pic50_predictor.py 文件中找到 predict() 方法（大约第850-900行）

找到这段代码：
    model = self.models.get(model_name)
    if model is None:
        self.logger.error(f"Model {model_name} not found")
        return None
    
    # 旧代码 - 删除这些 ❌
    if model_name == 'ensemble' and isinstance(model, dict):
        ...
    elif isinstance(model, nn.Module):
        ...
    else:
        prediction = model.predict([feature_vector])[0]  # ← 这行导致错误！

替换为：
    model = self.models.get(model_name)
    if model is None:
        self.logger.error(f"Model {model_name} not found")
        return None
    
    # 新代码 - 先判断类型 ✅
    if isinstance(model, dict):  # ← 注意：先判断dict！
        # Ensemble模型：使用所有基模型的平均预测
        predictions = []
        for base_model_name in ['random_forest', 'gradient_boosting', 'transformer']:
            if base_model_name in self.models:
                base_model = self.models[base_model_name]
                
                # 跳过ensemble自身（字典类型）
                if isinstance(base_model, dict):
                    continue
                
                # PyTorch模型
                if TORCH_AVAILABLE and isinstance(base_model, nn.Module):
                    device = next(base_model.parameters()).device
                    base_model.eval()
                    with torch.no_grad():
                        pred = base_model(torch.FloatTensor([feature_vector]).to(device)).cpu().numpy()[0]
                # Scikit-learn模型
                else:
                    pred = base_model.predict([feature_vector])[0]
                
                predictions.append(pred)
        
        if predictions:
            prediction = float(np.mean(predictions))
        else:
            self.logger.error("No base models available for ensemble")
            return None
    
    elif TORCH_AVAILABLE and isinstance(model, nn.Module):
        # PyTorch模型
        device = next(model.parameters()).device
        model.eval()
        with torch.no_grad():
            prediction = float(model(torch.FloatTensor([feature_vector]).to(device)).cpu().numpy()[0])
    
    else:
        # Scikit-learn模型
        prediction = float(model.predict([feature_vector])[0])

关键点：
1. ⚠️  第一个判断必须是 if isinstance(model, dict):
2. ⚠️  不要用 if model_name == 'ensemble' and isinstance(model, dict):
3. ⚠️  dict判断必须在最前面，不能在 elif 中

同样修改 predict_batch() 和 _estimate_uncertainty() 方法！

方法2：重新下载文件
--------------------
从artifact重新复制整个 complete_pic50_predictor.py 文件的内容。
    """)


if __name__ == "__main__":
    print("""
╔══════════════════════════════════════════════════════════════════╗
║                    代码版本检查工具                               ║
║                Code Version Checker                              ║
╚══════════════════════════════════════════════════════════════════╝
    """)
    
    is_correct = check_file_version()
    
    if not is_correct:
        show_fix_instructions()
        print("\n⚠️  请修复后重新运行测试")
    else:
        print("\n✅ 代码版本正确，可以正常使用！")
