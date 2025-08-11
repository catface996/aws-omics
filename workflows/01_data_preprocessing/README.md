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

---

# 🚀 生产就绪工作流：cow-preprocessing-fixed-final

## 工作流概述

`cow-preprocessing-fixed-final` (工作流ID: 3947992) 是一个专为奶牛基因组测序数据设计的预处理工作流，集成了质量控制、数据过滤和去重功能。该工作流已在AWS HealthOmics平台上成功验证，能够高效处理大规模基因组数据。

**状态**: ✅ ACTIVE | **创建时间**: 2025-08-10T02:50:40Z

## 工作流架构

```
输入FASTQ文件
    ↓
[InitialQC] ← FastQC质量控制
    ↓
[RunFastp] ← 数据过滤和修剪
    ↓
[RemoveDuplicates] ← 去除重复序列
    ↓
输出处理后的数据
```

## 核心功能

### 1. 初始质量控制 (InitialQC)
- **工具**: FastQC v0.12.1
- **功能**: 生成原始数据质量报告
- **资源配置**: 8 CPU, 8GB 内存
- **实例类型**: omics.c.2xlarge

### 2. 数据过滤和修剪 (RunFastp)
- **工具**: Fastp v0.23.4
- **功能**: 
  - 低质量序列过滤
  - 接头序列去除
  - PolyG/PolyX修剪
  - 长度过滤
- **资源配置**: 16 CPU, 32GB 内存
- **实例类型**: omics.c.4xlarge

### 3. 去重处理 (RemoveDuplicates)
- **工具**: SeqKit v2.5.1
- **功能**: 精确去重 (exact method)
- **资源配置**: 16 CPU, 24GB 内存
- **实例类型**: omics.c.4xlarge

## 参数配置

### 必需参数
- `input_fastq`: 输入FASTQ文件路径 (支持AWS Omics Sequence Store)
- `sample_name`: 样本名称

### 可选参数
| 参数名 | 默认值 | 描述 |
|--------|--------|------|
| `min_quality` | 20 | 最低质量阈值 |
| `min_length` | 50 | 最短序列长度 |
| `max_length` | 500 | 最长序列长度 |
| `complexity_threshold` | 30 | 复杂度阈值 |
| `enable_polyg_trimming` | true | 启用PolyG修剪 |
| `enable_polyx_trimming` | true | 启用PolyX修剪 |
| `dedup_method` | "exact" | 去重方法 |
| `fastp_cpu` | 16 | Fastp CPU核数 |
| `fastp_memory_gb` | 32 | Fastp内存大小 |
| `fastqc_cpu` | 8 | FastQC CPU核数 |
| `fastqc_memory_gb` | 8 | FastQC内存大小 |
| `dedup_cpu` | 16 | 去重CPU核数 |
| `dedup_memory_gb` | 24 | 去重内存大小 |

## 🎯 成功运行案例

### 运行信息
- **运行ID**: 5969296
- **运行名称**: cow-preprocessing-fixed-final-20250810-025200
- **状态**: ✅ COMPLETED
- **总运行时间**: 2小时32分钟 (7,232秒)
- **存储类型**: STATIC (2,400 GiB)

### 输入数据质量
基于SRR16760538样本的处理结果：
- **样本**: SRR16760538 (奶牛测序数据)
- **总读段数**: 197.04M
- **总碱基数**: 29.38G
- **Q20/Q30比例**: 99.998%
- **GC含量**: 43.81%
- **重复率**: 6.72% (单端数据可能高估)

### 任务执行详情

#### 1. InitialQC 任务
- **执行时间**: 21分钟3秒
- **CPU利用率**: 平均0.91核 (11.4%效率)
- **内存使用**: 平均1.64GB

#### 2. RunFastp 任务
- **执行时间**: 9分钟57秒
- **CPU利用率**: 平均4.52核 (28.3%效率)
- **内存使用**: 平均2.98GB

#### 3. RemoveDuplicates 任务
- **执行时间**: 1小时36分钟33秒 (最长任务)
- **CPU利用率**: 平均1.16核 (7.2%效率)
- **内存使用**: 平均3.63GB (峰值17.10GB)

## ⚡ 性能优化建议

基于实际运行数据的优化建议：

### CPU资源优化
- **RemoveDuplicates任务**: CPU利用率仅7.2%，建议减少到4-8核
- **InitialQC任务**: CPU利用率11.4%，建议减少到4核
- **RunFastp任务**: CPU利用率28.3%，配置合理

### 存储优化
- 当前STATIC存储利用率仅8.4%，建议改用DYNAMIC存储模式以节省成本

## 🚀 使用方法

### 1. 通过AWS CLI启动
```bash
aws omics start-run \
    --workflow-id 3947992 \
    --role-arn arn:aws:iam::YOUR-ACCOUNT:role/OmicsServiceRole \
    --name "cow-preprocessing-$(date +%Y%m%d-%H%M%S)" \
    --output-uri "s3://your-bucket/omics-outputs/preprocessing/" \
    --parameters '{
        "input_fastq": "s3://your-sequence-store-path/sample.fastq.gz",
        "sample_name": "your-sample-name",
        "min_quality": 20,
        "min_length": 50,
        "dedup_method": "exact"
    }' \
    --storage-type DYNAMIC
```

### 2. 推荐参数配置
```json
{
    "input_fastq": "s3://your-sequence-store/sample.fastq.gz",
    "sample_name": "sample_name",
    "min_quality": 20,
    "min_length": 50,
    "max_length": 500,
    "complexity_threshold": 30,
    "enable_polyg_trimming": true,
    "enable_polyx_trimming": true,
    "dedup_method": "exact",
    "fastp_cpu": 16,
    "fastp_memory_gb": 32,
    "fastqc_cpu": 4,
    "fastqc_memory_gb": 8,
    "dedup_cpu": 8,
    "dedup_memory_gb": 24
}
```

## 📊 输出文件

工作流输出存储在S3路径：
```
s3://catface996-genomic/omics-outputs/preprocessing-fixed-final/{run-id}/
```

### 主要输出文件
- `{sample_name}_filtered.fastq.gz` - 过滤后的FASTQ文件
- `{sample_name}_dedup.fastq.gz` - 去重后的FASTQ文件
- `{sample_name}_fastqc.html` - FastQC质量报告
- `{sample_name}_fastp.html` - Fastp处理报告
- `{sample_name}_fastp.json` - Fastp统计数据

## 🔗 下一步流程

预处理完成后，清理的数据将用于：
- **02_sequence_alignment**: 序列比对分析
- **质量验证**: 确保数据质量满足下游分析要求

## 📚 相关文档

- [预处理工作流完整指南](./预处理工作流完整指南.md)
- [AWS HealthOmics 用户指南](https://docs.aws.amazon.com/omics/)
- [FastQC 文档](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Fastp 文档](https://github.com/OpenGene/fastp)
- [SeqKit 文档](https://bioinf.shenwei.me/seqkit/)

---

**最后更新**: 2025年8月10日  
**状态**: ✅ 生产就绪  
**注意**: 数据预处理是整个基因组分析流程的基础，高质量的预处理直接影响后续分析结果的准确性。
