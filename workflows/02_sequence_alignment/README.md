# 02_序列比对 (Sequence Alignment)

## 📋 目录说明

本目录包含基因组序列比对相关的工作流和工具，负责将预处理后的测序reads比对到参考基因组上。

## 🔄 序列比对工作流任务详解

### 任务1: BuildReferenceIndex - 参考基因组索引构建 🔧

**用途**: 为参考基因组构建必要的索引文件，以支持高效的序列比对

**详细功能**:
- **BWA索引构建**: 创建BWA-MEM算法所需的索引文件(.amb, .ann, .bwt, .pac, .sa)
- **SAMtools索引**: 生成FASTA索引文件(.fai)，用于快速随机访问参考序列
- **序列字典**: 创建序列字典文件(.dict)，包含染色体信息和长度

**技术细节**:
- **工具**: BWA index + SAMtools faidx
- **输入**: 参考基因组FASTA文件 (ARS-UCD1.2, ~2.75GB)
- **输出**: 多个索引文件 (总计约5-6GB)
- **资源配置**: 8 CPU, 16GB内存, 200GB存储
- **运行时间**: 约50分钟 (已完成 ✅)

**重要性**: 
- 索引构建是一次性操作，后续比对可重复使用
- 索引质量直接影响比对速度和准确性
- 大型基因组(如奶牛2.75GB)需要充足的内存和存储空间

### 任务2: BWAAlignment - BWA序列比对 🎯

**用途**: 将预处理后的测序reads比对到参考基因组上，生成比对结果

**详细功能**:
- **序列比对**: 使用BWA-MEM算法进行高精度比对
- **读组信息**: 添加样本、文库、测序平台等元数据
- **比对参数优化**: 针对短读长测序数据优化参数

**技术细节**:
- **工具**: BWA-MEM v0.7.17
- **输入**: 清洁FASTQ文件 (~1.3GB压缩) + 参考基因组索引
- **输出**: SAM格式比对文件 (~15-20GB未压缩)
- **资源配置**: 16 CPU, 32GB内存, 400GB存储
- **当前状态**: 正在运行 🔄 (已运行约1小时)

**关键参数**:
- **最小种子长度**: 19bp (平衡敏感性和特异性)
- **带宽**: 100 (允许的插入缺失范围)
- **比对评分**: A=1, B=4, O=6, E=1 (匹配/错配/gap惩罚)
- **读组标签**: @RG\tID:样本ID\tSM:样本名\tPL:ILLUMINA

**性能特点**:
- **高通量**: 处理数百万条reads
- **高精度**: 支持复杂基因组区域比对
- **内存效率**: 流式处理减少内存占用

### 任务3: ProcessAlignment - 比对后处理和质量评估 📊

**用途**: 将SAM文件转换为标准BAM格式，并进行质量过滤和统计分析

**详细功能**:
- **格式转换**: SAM转BAM，减少存储空间
- **质量过滤**: 移除低质量比对 (MAPQ < 20)
- **排序索引**: 按基因组坐标排序并创建索引
- **统计分析**: 生成详细的比对质量报告

**技术细节**:
- **工具**: SAMtools v1.17
- **输入**: SAM比对文件 (~15-20GB)
- **输出**: 排序BAM文件 (~3-5GB) + 索引 + 统计报告
- **资源配置**: 8 CPU, 16GB内存, 500GB存储
- **当前状态**: 等待中 ⏳

**处理步骤**:
1. **SAM转BAM**: `samtools view -b` (压缩比约4:1)
2. **质量过滤**: 移除未比对reads (-F 4) 和低质量比对 (-q 20)
3. **坐标排序**: `samtools sort` 按染色体位置排序
4. **索引创建**: `samtools index` 生成.bai索引文件
5. **统计生成**: 
   - `samtools stats`: 详细比对统计
   - `samtools flagstat`: 比对标志统计
   - `samtools idxstats`: 每条染色体比对统计

**质量指标**:
- **总比对率**: 期望 >95%
- **唯一比对率**: 期望 >90%
- **重复率**: 通常 <20%
- **插入片段大小**: 期望均值约300-500bp

## 🎯 工作流整体架构

```
预处理FASTQ → [索引构建] → [BWA比对] → [BAM处理] → 最终BAM文件
     ↓              ↓           ↓           ↓           ↓
  1.3GB压缩      索引文件     SAM文件     BAM文件     统计报告
                 (~6GB)      (~20GB)     (~5GB)      (数MB)
```

## 📈 性能基准和预期

### 当前运行状态 (运行ID: 1809888)
- **总预计运行时间**: 2.5-3小时
- **已完成**: BuildReferenceIndex (50分钟)
- **进行中**: BWAAlignment (1小时+，预计还需30-60分钟)
- **待执行**: ProcessAlignment (预计20-30分钟)

### 资源利用效率
- **CPU利用率**: BWA比对阶段约80-90%
- **内存使用**: 峰值约30GB (BWA比对阶段)
- **I/O吞吐**: 读取约100MB/s，写入约50MB/s
- **网络传输**: S3读写约20-30MB/s

## 🔧 工作流特性

### 容错和监控
- **自动重试**: 任务失败自动重试3次
- **CloudWatch日志**: 详细的执行日志记录
- **状态监控**: 实时任务状态和资源使用监控
- **错误诊断**: 自动错误分类和建议

### 数据管理
- **输入验证**: 自动检查FASTQ和参考基因组完整性
- **中间文件**: 自动清理临时文件节省存储
- **输出组织**: 结构化输出目录便于后续分析

## 📁 文件结构

```
02_sequence_alignment/
├── README.md                           # 本说明文档 (已更新)
├── sequence_alignment_v2_ecr.wdl      # 主要工作流 (使用私有ECR)
├── sequence_alignment_v1_complete.wdl # 完整版工作流
├── multi_file_alignment_workflow.wdl  # 多文件比对工作流
└── sequence_alignment_v2_ecr.zip      # 部署包
```

## 🚀 使用方法

### 启动序列比对工作流
```bash
aws omics start-run \
    --workflow-id 2495995 \
    --role-arn arn:aws:iam::ACCOUNT:role/OmicsServiceRole \
    --name "cow-sequence-alignment-$(date +%Y%m%d-%H%M%S)" \
    --output-uri "s3://your-bucket/omics-outputs/sequence-alignment/" \
    --parameters '{
        "sample_name": "SRR16760538",
        "cleaned_fastq": "s3://path/to/cleaned.fastq.gz",
        "reference_genome": "s3://path/to/reference.fna",
        "bwa_cpu": 16,
        "bwa_memory_gb": 32,
        "samtools_cpu": 8,
        "samtools_memory_gb": 16
    }'
```

### 监控运行状态
```bash
# 检查运行状态
aws omics get-run --id <run-id>

# 查看任务详情
aws omics list-run-tasks --id <run-id>

# 获取运行日志
aws omics get-run-logs --id <run-id>
```

## 📊 输出结果

### 主要输出文件
- **final_bam**: 排序后的比对结果 (约3-5GB)
- **final_bam_index**: BAM文件索引 (.bai)
- **alignment_stats**: 详细比对统计 (samtools stats)
- **flagstat_report**: 比对标志统计 (samtools flagstat)
- **idxstats_report**: 染色体比对统计 (samtools idxstats)

### 输出位置
```
s3://catface996-genomic/omics-outputs/sequence-alignment-v2-final/1809888/
├── out/
│   ├── final_bam/
│   │   └── SRR16760538.sorted.bam
│   ├── final_bam_index/
│   │   └── SRR16760538.sorted.bam.bai
│   ├── alignment_stats/
│   │   └── SRR16760538_alignment_stats.txt
│   ├── flagstat_report/
│   │   └── SRR16760538_flagstat.txt
│   └── idxstats_report/
│       └── SRR16760538_idxstats.txt
└── reference_index_files/
    └── [BWA和SAMtools索引文件]
```

## 🔗 数据流向

### 输入来源
- **预处理工作流**: 清洁FASTQ文件 (运行ID: 5969296)
- **参考基因组**: ARS-UCD1.2 奶牛基因组

### 输出去向
- **变异检测工作流**: BAM文件用于SNP/INDEL检测
- **结构变异分析**: BAM文件用于大片段变异检测
- **覆盖度分析**: 基因组覆盖度统计

## 📚 相关文档

- [BWA用户手册](http://bio-bwa.sourceforge.net/bwa.shtml)
- [SAMtools文档](http://www.htslib.org/doc/samtools.html)
- [AWS Omics最佳实践](https://docs.aws.amazon.com/omics/latest/dev/workflows-best-practices.html)
- [SAM/BAM格式规范](https://samtools.github.io/hts-specs/SAMv1.pdf)

---

**最后更新**: 2025年8月10日  
**当前状态**: BWA比对任务正在运行中 🔄  
**预计完成**: 约1-1.5小时后
