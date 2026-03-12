import torch
import numpy
import pandas
import scipy
import sklearn
import matplotlib

import requests
import bs4
import tqdm


print("=== 安装验证 ===")

print(f"PyTorch  : {torch.__version__}")
print(f"NumPy    : {numpy.__version__}")

print(f"GPU 可用 : {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"GPU 型号 : {torch.cuda.get_device_name(0)}")
    mem = torch.cuda.get_device_properties(0).total_memory
    print(f"显存大小 : {mem / 1024**3:.1f} GB")

print("所有依赖安装成功！")
