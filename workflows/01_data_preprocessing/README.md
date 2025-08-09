# 01_数据预处理 (Data Preprocessing)

## 📋 目录说明

本目录包含基因组测序数据预处理相关的工作流和工具，是生物信息学分析流程的第一步。

## 🎯 主要功能

### 质量控制 (Quality Control)
- **FastQC**: 测序数据质量评估
- **MultiQC**: 多样本质量报告汇总
- 质量分数分布分析
- 序列长度分布检查
- GC含量分析

### 数据清理 (Data Cleaning)
- **Trimmomatic**: 接头序列去除
- **FastP**: 高效的质量过滤和接头去除
- 低质量碱基修剪
- 最小长度过滤
- N碱基处理

### 去重复 (Deduplication)
- PCR重复序列识别和去除
- 光学重复序列处理
- 序列复杂度评估

## 📁 文件结构

```
01_data_preprocessing/
├── README.md                           # 本说明文档
├── preprocessing_workflow.wdl          # 主预处理工作流
├── 预处理工作流完整指南.md              # 详细使用指南
├── tasks/                              # 任务定义目录
│   ├── fastqc_task.wdl                # FastQC质量控制任务
│   ├── trimmomatic_task.wdl           # Trimmomatic接头去除任务
│   ├── fastp_task.wdl                 # FastP质量过滤任务
│   ├── deduplication_task.wdl         # 去重复任务
│   └── multiqc_task.wdl               # MultiQC报告生成任务
├── inputs/                             # 输入参数配置
│   └── preprocessing_inputs.json       # 预处理参数配置文件
└── examples/                           # 示例文件
    └── example_workflow.wdl           # 示例工作流
```

## 🚀 使用方法

### 1. 基本预处理流程
```bash
# 使用AWS Omics运行预处理工作流
aws omics start-run \
    --workflow-id <workflow-id> \
    --parameters file://inputs/preprocessing_inputs.json \
    --name "cattle-preprocessing-$(date +%Y%m%d)"
```

### 2. 输入数据要求
- **FASTQ文件**: 原始测序数据（支持单端和双端测序）
- **质量格式**: Phred+33质量编码
- **文件格式**: 支持压缩格式（.gz, .bz2）

### 3. 输出结果
- **清理后的FASTQ文件**: 去除接头和低质量序列
- **质量控制报告**: HTML格式的详细质量分析
- **统计信息**: 处理前后的数据统计对比

## ⚙️ 参数配置

### 质量过滤参数
- **最小质量分数**: 20 (Phred score)
- **最小序列长度**: 50 bp
- **滑动窗口**: 4:20 (窗口大小:平均质量)

### 接头序列
- **Illumina通用接头**: TruSeq3-PE-2.fa
- **自定义接头**: 支持用户定义的接头序列

## 📊 质量指标

### 处理前检查
- 总序列数量
- 平均序列长度
- 质量分数分布
- GC含量分布

### 处理后验证
- 保留序列比例
- 质量改善程度
- 接头去除效率
- 重复序列去除率

## 🔗 下一步流程

预处理完成后，清理的数据将用于：
- **02_sequence_alignment**: 序列比对分析
- **质量验证**: 确保数据质量满足下游分析要求

## 📚 相关文档

- [预处理工作流完整指南](./预处理工作流完整指南.md)
- [AWS Omics 用户指南](https://docs.aws.amazon.com/omics/)
- [生物信息学最佳实践](../docs/bioinformatics_best_practices.md)

---

**注意**: 数据预处理是整个基因组分析流程的基础，高质量的预处理直接影响后续分析结果的准确性。
