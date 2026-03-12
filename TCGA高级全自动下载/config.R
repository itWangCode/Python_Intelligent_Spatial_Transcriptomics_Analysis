# TCGA分析配置文件
# 复制这个文件并重命名为 config.R，然后修改参数

# ========================== 工作目录设置 ==========================
WORK_DIR <- "/Users/wangyang/Desktop/BCBM_TCGA_download"

# ========================== 项目设置 ==========================
# TCGA项目ID
# 常见项目:
#   TCGA-BRCA: 乳腺癌
#   TCGA-LUAD: 肺腺癌
#   TCGA-LUSC: 肺鳞癌
#   TCGA-PRAD: 前列腺癌
#   TCGA-COAD: 结肠腺癌
#   TCGA-READ: 直肠腺癌
PROJECT_ID <- "TCGA-BRCA"

# ========================== 差异分析参数 ==========================
# log2 Fold Change阈值（绝对值）
THRESHOLD_LOGFC <- 1.0

# 校正后p值阈值
THRESHOLD_ADJP <- 0.05

# ========================== 可视化参数 ==========================
# 热图展示的最大基因数（每个方向）
MAX_DISPLAY_GENES <- 50

# 火山图标注的top基因数
TOP_GENES_VOLCANO <- 20

# ========================== 过滤参数 ==========================
# 基因过滤：至少在N个样本中表达
MIN_SAMPLE_COUNT <- 10

# 基因过滤：最小count值
MIN_COUNT <- 10

# ========================== 生存分析参数 ==========================
# 用于生存分析的基因（留空则使用top差异基因）
SURVIVAL_GENE <- ""

# 表达分组方法: "median" 或 "quartile"
SURVIVAL_GROUP_METHOD <- "median"

# ========================== 下载设置 ==========================
# 每批下载的文件数
FILES_PER_CHUNK <- 10

# 下载方法: "api" 或 "client"
DOWNLOAD_METHOD <- "api"

# ========================== 其他设置 ==========================
# 是否强制重新运行所有步骤（忽略断点续传）
FORCE_RERUN <- FALSE

# 随机种子（用于可重复性）
RANDOM_SEED <- 42

# CPU核心数（用于并行计算）
N_CORES <- 4

# ========================== 颜色方案 ==========================
# 肿瘤样本颜色
COLOR_TUMOR <- "#E74C3C"

# 正常样本颜色  
COLOR_NORMAL <- "#3498DB"

# 上调基因颜色
COLOR_UP <- "#E74C3C"

# 下调基因颜色
COLOR_DOWN <- "#3498DB"

# 不显著颜色
COLOR_NS <- "#95A5A6"

# ========================== 高级设置 ==========================
# 是否使用DESeq2（FALSE则使用limma）
USE_DESEQ2 <- TRUE

# 是否进行批次效应校正
CORRECT_BATCH <- FALSE

# 批次信息列名（如果CORRECT_BATCH为TRUE）
BATCH_COLUMN <- "plate"

# 是否保存中间结果
SAVE_INTERMEDIATE <- TRUE

# 日志级别: "INFO", "DEBUG", "WARNING", "ERROR"
LOG_LEVEL <- "INFO"
