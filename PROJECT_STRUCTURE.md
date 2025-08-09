# AWS Omics 奶牛基因组学项目结构

## 📁 目录结构

```
aws-omics/
├── 📄 项目根文件
│   ├── README.md                          # 项目主要说明
│   ├── PROJECT_STRUCTURE.md               # 项目结构说明 (本文件)
│   ├── .gitignore                         # Git忽略文件
│   ├── .gitattributes                     # Git属性配置
│   └── download.log                       # 下载日志
│
├── 🛠️ scripts/                           # 所有脚本文件
│   ├── simulation/                        # 基因组模拟脚本
│   │   ├── simulate_full_genome_gpu.py    # GPU全基因组模拟
│   │   ├── simulate_single_chromosome_gpu.py # GPU单染色体模拟
│   │   ├── parallel_genome_simulation.py  # 并行基因组模拟
│   │   ├── generate_mutated_fna_gpu.py    # GPU生成变异基因组 ⭐
│   │   ├── generate_mutated_genome.py     # 生成变异基因组
│   │   └── generate_sequencing_data_gpu.py # GPU生成测序数据
│   ├── analysis/                          # 数据分析脚本
│   │   └── analyze_genomic_data.py        # 基因组数据分析
│   ├── utils/                             # 工具脚本
│   │   ├── consolidate_simulation_data.py # 数据整合工具
│   │   └── aws_omics_workflow_example.py  # AWS Omics示例
│   ├── run_mutated_fna_generation.sh      # 运行变异基因组生成 ⭐
│   ├── run_sequencing_simulation.sh       # 运行测序模拟
│   ├── run_parallel_simulation.sh         # 运行并行模拟
│   ├── download_genomic_data.sh           # 下载基因组数据
│   ├── deploy_omics_simulation.sh         # 部署Omics模拟
│   └── 部署牛基因组学.sh                   # 中文部署脚本
│
├── 📚 docs/                              # 文档目录
│   ├── user_guides/                       # 用户指南
│   │   ├── 快速入门指南.md                 # 快速开始
│   │   ├── 牛基因组下载指南.md             # 数据下载指南
│   │   ├── GPU测序数据生成说明.md          # GPU测序说明 ⭐
│   │   └── GPU性能优化说明.md              # GPU优化说明 ⭐
│   ├── technical/                         # 技术文档
│   │   ├── 基因工程术语.md                 # 术语解释
│   │   ├── 牛基因组学AWS_Omics演示.md      # AWS Omics演示
│   │   ├── 基因组学分析工作流程.md         # 分析流程
│   │   └── 基因组数据下载总结.md           # 下载总结
│   └── archive/                           # 过时文档存档
│       ├── 并行基因组模拟说明.md           # 已过时
│       ├── GPU模拟脚本使用说明.md          # 已过时
│       ├── 模拟测序数据使用指南.md         # 已过时
│       ├── 项目完整总结.md                 # 已过时
│       ├── 完整基因组Mock数据生成报告.md   # 已过时
│       └── 数据合理性分析.md               # 已过时
│
├── 🔄 workflows/                         # 工作流文件
│   ├── wdl/                              # WDL工作流
│   │   ├── variant_calling_workflow.wdl   # 变异检测工作流
│   │   └── multi_file_alignment_workflow.wdl # 多文件比对工作流
│   └── examples/                          # 工作流示例
│
├── 💾 data/                              # 数据目录
│   ├── reference/                         # 参考数据 (软链接到genomic_data)
│   └── variants/                          # 变异数据
│       └── consolidated_simulation_data/  # 整合的模拟数据
│
├── 📊 outputs/                           # 输出结果
│   ├── simulation/                        # 模拟结果
│   ├── sequencing/                        # 测序数据
│   │   └── sequencing_test/               # 测序测试结果
│   └── mutated_genome/                    # 变异基因组
│       ├── mutated_genome_data/           # 变异基因组数据
│       └── mutated_genome_test/           # 变异基因组测试
│
├── 🧬 genomic_data/                      # 原始基因组数据 (Git LFS)
│   ├── reference/                         # 参考基因组
│   ├── annotation/                        # 基因注释
│   ├── protein/                           # 蛋白质数据
│   └── ...                               # 其他基因组数据
│
└── 🐍 gpu_env/                           # Python GPU环境
    └── ...                               # 虚拟环境文件
```

## 🎯 主要功能模块

### 1. 核心脚本 (⭐ 推荐使用)
- **`scripts/simulation/generate_mutated_fna_gpu.py`**: GPU加速生成变异基因组FNA文件
- **`scripts/run_mutated_fna_generation.sh`**: 一键运行变异基因组生成

### 2. 数据流程
```
参考基因组 (genomic_data/reference/)
    ↓
变异数据 (data/variants/consolidated_simulation_data/)
    ↓
GPU并行处理 (scripts/simulation/)
    ↓
变异基因组 (outputs/mutated_genome/)
```

### 3. 文档分类
- **用户指南**: 面向用户的操作说明
- **技术文档**: 深入的技术细节
- **存档文档**: 过时但保留的历史文档

## 🚀 快速开始

### 1. 生成变异基因组 (推荐)
```bash
# 使用自动化脚本
./scripts/run_mutated_fna_generation.sh

# 或直接运行
cd scripts/simulation
python generate_mutated_fna_gpu.py \
    --reference ../../genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    --variants ../../data/variants/consolidated_simulation_data/consolidated_variants.vcf \
    --output ../../outputs/mutated_genome/test
```

### 2. 查看文档
```bash
# 用户指南
cat docs/user_guides/GPU测序数据生成说明.md

# 技术文档
cat docs/technical/牛基因组学AWS_Omics演示.md
```

## 🧹 清理说明

### 已整理的改进
1. **脚本分类**: 按功能分为simulation、analysis、utils
2. **文档分层**: 用户指南、技术文档、存档分离
3. **数据组织**: 输入数据和输出结果分离
4. **过时内容**: 移至archive目录，不删除以保留历史

### 推荐使用的文件
- ✅ `scripts/simulation/generate_mutated_fna_gpu.py` - 最新GPU优化脚本
- ✅ `scripts/run_mutated_fna_generation.sh` - 自动化运行脚本
- ✅ `docs/user_guides/GPU测序数据生成说明.md` - 最新使用说明
- ✅ `docs/user_guides/GPU性能优化说明.md` - 性能优化指南

### 可以忽略的文件
- 📦 `docs/archive/` - 过时文档，仅作历史参考
- 📦 旧版本的模拟脚本 - 已被GPU优化版本替代

## 📋 维护建议

1. **新增脚本**: 放入对应的scripts子目录
2. **新增文档**: 根据用途放入docs的相应子目录
3. **输出数据**: 统一放入outputs目录
4. **过时内容**: 移至archive而非删除

---

**最后更新**: 2025年8月8日  
**项目状态**: ✅ 结构化整理完成
