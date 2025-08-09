# AWS Omics 奶牛基因组学演示项目

本项目展示如何使用AWS Omics服务构建完整的奶牛基因组分析流水线，从原始测序数据到生物学解释的端到端解决方案。

## 🚀 快速开始

## 0-下载数据

https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR16760538&display=metadata

### 1. 生成变异基因组 (推荐)
```bash
# 激活GPU环境并运行
source gpu_env/bin/activate
./scripts/run_mutated_fna_generation.sh
```

### 2. 查看项目结构
```bash
# 查看详细的项目结构说明
cat PROJECT_STRUCTURE.md
```

## 📁 项目结构概览

```
aws-omics/
├── scripts/           # 所有脚本文件
│   ├── simulation/    # 基因组模拟脚本 ⭐
│   ├── analysis/      # 数据分析脚本
│   └── utils/         # 工具脚本
├── docs/              # 文档目录
│   ├── user_guides/   # 用户指南 ⭐
│   ├── technical/     # 技术文档
│   └── archive/       # 过时文档存档
├── workflows/         # WDL工作流文件
├── data/              # 输入数据
├── outputs/           # 输出结果
└── genomic_data/      # 原始基因组数据 (Git LFS)
```

## 🎯 核心功能

### 1. GPU加速变异基因组生成
- **脚本**: `scripts/simulation/generate_mutated_fna_gpu.py`
- **功能**: 基于3,154,016个变异生成完整的变异基因组
- **性能**: 充分利用Apple M4 Max GPU，~100,000变异/秒

### 2. 测序数据模拟
- **脚本**: `scripts/simulation/generate_sequencing_data_gpu.py`
- **功能**: 生成高质量的配对端测序数据
- **格式**: 标准FASTQ格式，支持不同覆盖深度

### 3. AWS Omics集成
- **工作流**: `workflows/wdl/`
- **功能**: 完整的变异检测和分析流水线

## 📊 数据说明

### 输入数据
- **参考基因组**: 奶牛 ARS-UCD1.2 版本 (~2.75GB)
- **变异数据**: 3,154,016个高质量变异 (230MB VCF)

### 输出数据
- **变异基因组**: 完整的FNA格式基因组文件
- **测序数据**: 标准FASTQ格式测序reads
- **分析报告**: 详细的处理统计信息

## 🔧 系统要求

- **操作系统**: macOS (支持MPS)
- **Python**: 3.8+ with PyTorch
- **GPU**: Apple Silicon (M1/M2/M4) 推荐
- **内存**: 16GB+ 推荐
- **存储**: 10GB+ 可用空间

## 📖 使用指南

### 快速测试
```bash
# 测试3个染色体的变异基因组生成
source gpu_env/bin/activate
python scripts/simulation/generate_mutated_fna_gpu.py \
    --reference genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    --variants data/variants/consolidated_simulation_data/consolidated_variants.vcf \
    --output outputs/mutated_genome/test \
    --chromosomes NC_037328.1 NC_037329.1 NC_037330.1 \
    --max-workers 16
```

### 完整运行
```bash
# 处理所有30个主要染色体
./scripts/run_mutated_fna_generation.sh
```

## 📚 文档

- **用户指南**: `docs/user_guides/` - 操作说明和快速开始
- **技术文档**: `docs/technical/` - 深入的技术细节
- **项目结构**: `PROJECT_STRUCTURE.md` - 完整的项目结构说明

## ⚡ 性能特性

- **GPU加速**: 充分利用Apple M4 Max GPU
- **大规模并行**: 50,000变异/批次处理
- **内存优化**: 智能GPU内存管理
- **高吞吐量**: ~100,000变异/秒处理速度

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个项目！

## 📄 许可证

本项目遵循开源许可证。数据使用请遵循NCBI的使用条款。

---

**最后更新**: 2025年8月8日  
**项目状态**: ✅ 结构化整理完成，GPU优化就绪
