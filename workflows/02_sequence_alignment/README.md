# 02_序列比对 (Sequence Alignment)

## 📋 目录说明

本目录包含基因组序列比对相关的工作流和工具，负责将预处理后的测序reads比对到参考基因组上。

## 🎯 主要功能

### 序列比对 (Read Alignment)
- **BWA-MEM**: 高效的短序列比对算法
- **Bowtie2**: 替代比对工具
- **STAR**: RNA-seq数据比对（如需要）
- 多文件并行比对处理

### 比对后处理 (Post-alignment Processing)
- **SAM/BAM格式转换**: 标准化比对结果格式
- **排序和索引**: 提高后续分析效率
- **比对质量评估**: 比对率和覆盖度统计
- **重复标记**: 标记PCR和光学重复

### 质量控制 (Quality Assessment)
- 比对率统计
- 插入片段大小分布
- 覆盖深度分析
- 比对质量分数分布

## 📁 文件结构

```
02_sequence_alignment/
├── README.md                           # 本说明文档
├── multi_file_alignment_workflow.wdl  # 多文件比对主工作流
├── tasks/                              # 任务定义目录（待创建）
│   ├── bwa_alignment_task.wdl         # BWA比对任务
│   ├── sam_processing_task.wdl        # SAM/BAM处理任务
│   ├── alignment_qc_task.wdl          # 比对质量控制任务
│   └── duplicate_marking_task.wdl     # 重复标记任务
├── inputs/                             # 输入参数配置（待创建）
│   └── alignment_inputs.json          # 比对参数配置文件
└── scripts/                            # 辅助脚本（待创建）
    ├── prepare_reference.sh           # 参考基因组索引准备
    └── alignment_stats.py             # 比对统计分析
```

## 🚀 使用方法

### 1. 参考基因组准备
```bash
# 创建BWA索引
bwa index reference_genome.fasta

# 创建samtools索引
samtools faidx reference_genome.fasta

# 创建序列字典
picard CreateSequenceDictionary \
    R=reference_genome.fasta \
    O=reference_genome.dict
```

### 2. 多文件比对流程
```bash
# 使用AWS Omics运行多文件比对工作流
aws omics start-run \
    --workflow-id <workflow-id> \
    --parameters file://inputs/alignment_inputs.json \
    --name "cattle-alignment-$(date +%Y%m%d)"
```

### 3. 输入数据要求
- **预处理后的FASTQ文件**: 来自01_data_preprocessing步骤
- **参考基因组**: 奶牛ARS-UCD1.2基因组
- **读组信息**: 样本ID、文库信息、测序平台

## ⚙️ 比对参数配置

### BWA-MEM参数
- **最小种子长度**: 19
- **带宽**: 100
- **Z-dropoff**: 100
- **3'端剪切惩罚**: 5

### SAM/BAM处理
- **排序方式**: 按坐标排序
- **压缩级别**: 6 (平衡速度和存储)
- **索引创建**: 自动生成BAI索引

### 质量过滤
- **最小比对质量**: MAPQ ≥ 20
- **去除未比对reads**: 是
- **去除次要比对**: 是

## 📊 输出结果

### 主要输出文件
- **sorted.bam**: 排序后的比对结果
- **sorted.bam.bai**: BAM文件索引
- **alignment_stats.txt**: 比对统计信息
- **coverage_report.html**: 覆盖度分析报告

### 质量指标
- **总比对率**: >95% (高质量数据)
- **唯一比对率**: >90%
- **平均覆盖深度**: 30X (推荐)
- **覆盖均匀性**: CV < 0.3

## 🔧 工作流特性

### 多文件支持
- 同时处理多个FASTQ文件对
- 自动合并同一样本的多个文库
- 保持读组信息完整性

### 并行处理
- 染色体级别并行比对
- 多线程BWA比对
- 分布式计算优化

### 错误处理
- 自动重试机制
- 中间文件清理
- 详细日志记录

## 📈 性能优化

### 计算资源
- **CPU**: 16-32核心推荐
- **内存**: 64GB推荐
- **存储**: 高IOPS SSD存储

### 优化策略
- 参考基因组预加载
- 流式处理减少I/O
- 智能任务调度

## 🔗 数据流向

### 输入来源
- **01_data_preprocessing**: 清理后的FASTQ文件
- **genomic_data/02_reference_genome**: 参考基因组文件

### 输出去向
- **03_variant_detection**: BAM文件用于变异检测
- **genomic_data/03_alignment_results**: 存储比对结果

## 📚 相关文档

- [BWA用户手册](http://bio-bwa.sourceforge.net/bwa.shtml)
- [SAMtools文档](http://www.htslib.org/doc/samtools.html)
- [AWS Omics最佳实践](https://docs.aws.amazon.com/omics/latest/dev/workflows-best-practices.html)

---

**注意**: 序列比对是计算密集型步骤，合理的参数配置和资源分配对性能至关重要。
