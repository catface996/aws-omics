# 奶牛基因组数据下载摘要

## 下载信息
- 下载时间: Fri  8 Aug 2025 10:51:19 CST
- 数据源: NCBI Genome Database
- 基因组版本: ARS-UCD1.2
- 物种: Bos taurus (奶牛)

## 目录结构

目录按照生物信息学分析流程顺序组织：**数据预处理** → **序列比对** → **变异检测** → **变异注释** → **统计分析** → **生物学解释**

```
genomic_data/
├── 01_raw_data/              # 原始测序数据 + 质量控制 (数据预处理阶段)
├── 02_reference_genome/      # 参考基因组序列 (序列比对阶段)
├── 03_alignment_results/     # 序列比对结果 (序列比对阶段)
├── 04_genome_assembly/       # 基因组组装信息 (序列比对阶段)
├── 05_genome_annotation/     # 基因组注释数据 (序列比对 + 变异注释阶段)
├── 06_genomic_features/      # 基因组特征数据 (序列比对 + 变异注释阶段)
├── 07_rna_sequences/         # RNA序列数据 (生物学解释阶段)
└── 08_protein_sequences/     # 蛋白质序列数据 (生物学解释阶段)
```

### 目录详细说明

#### 🔬 数据预处理阶段
- **01_raw_data/**: 从NCBI SRA下载的原始FASTQ测序数据
  - 包含配对端测序reads (SRR16760538_1.fastq.gz, SRR16760538_2.fastq.gz)
  - 包含MD5校验和文件 (md5checksums.txt)
  - 用于质量控制、去接头、过滤等预处理步骤

#### 🧬 序列比对阶段  
- **02_reference_genome/**: 奶牛参考基因组 ARS-UCD1.2
  - 用于reads比对的标准基因组序列
  - 变异检测的参考标准

- **03_alignment_results/**: 序列比对输出文件
  - SAM/BAM格式的比对结果文件
  - BAM索引文件 (.bai)
  - 比对统计和质量报告

- **04_genome_assembly/**: 基因组组装统计和报告
  - 染色体/scaffold结构信息
  - 组装质量统计数据

#### 📝 序列比对 + 变异注释阶段
- **05_genome_annotation/**: 基因组功能注释
  - GFF/GTF格式的基因注释文件
  - 用于指导比对策略和变异功能注释

- **06_genomic_features/**: 基因组特征数据
  - 基因、外显子、内含子边界信息
  - 用于变异影响预测和功能分析

#### 🔍 生物学解释阶段 (遵循中心法则: DNA → RNA → Protein)
- **07_rna_sequences/**: RNA序列数据 (转录产物)
  - mRNA、tRNA、rRNA等序列信息
  - 转录组水平的功能分析
  - 蛋白质翻译的直接模板

- **08_protein_sequences/**: 蛋白质序列数据 (翻译产物)
  - 用于预测变异对蛋白质功能的影响
  - 功能域和结构分析
  - 最终的基因表达产物

## 主要文件说明

### 原始测序数据 (01_raw_data/)
- `SRR16760538_1.fastq.gz`: 配对端测序数据 R1
- `SRR16760538_2.fastq.gz`: 配对端测序数据 R2
- `SRR16760538.fastq.gz`: 单端测序数据
- `md5checksums.txt`: MD5校验和文件

### 参考基因组 (02_reference_genome/)
- `GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz`: 参考基因组FASTA格式
- `GCF_002263795.1_ARS-UCD1.2_genomic.fna`: 解压后的参考基因组

### 序列比对结果 (03_alignment_results/)
- `aligned.sam`: 原始比对结果 (SAM格式，可选)
- `aligned_sorted.bam`: 排序后的BAM文件
- `aligned_sorted.bam.bai`: BAM索引文件
- `alignment_stats.txt`: 比对统计信息
- `alignment_summary.html`: 比对质量报告

### 基因组组装信息 (04_genome_assembly/)
- `GCF_002263795.1_ARS-UCD1.2_assembly_report.txt`: 组装报告
- `GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt`: 组装统计信息

### 基因组注释 (05_genome_annotation/)
- `GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz`: GTF格式注释
- `GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz`: GFF3格式注释
- `GCF_002263795.1_ARS-UCD1.2_cds_from_genomic.fna.gz`: CDS序列

## 分析流程建议

### 1. 数据预处理阶段
```bash
# 验证数据完整性
cd 01_raw_data/
md5sum -c md5checksums.txt

# 质量控制检查
fastqc *.fastq.gz
```

### 2. 序列比对阶段
```bash
# 解压参考基因组
cd 02_reference_genome/
gunzip GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz

# 建立BWA索引
bwa index GCF_002263795.1_ARS-UCD1.2_genomic.fna

# 序列比对并保存到alignment_results目录
bwa mem 02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    01_raw_data/SRR16760538_1.fastq.gz \
    01_raw_data/SRR16760538_2.fastq.gz > 03_alignment_results/aligned.sam

# SAM转BAM并排序
cd 03_alignment_results/
samtools view -bS aligned.sam | samtools sort -o aligned_sorted.bam
samtools index aligned_sorted.bam

# 生成比对统计
samtools flagstat aligned_sorted.bam > alignment_stats.txt
```

### 3. 变异检测阶段
```bash
# 变异检测
bcftools mpileup -f 02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    03_alignment_results/aligned_sorted.bam | bcftools call -mv -Oz -o variants.vcf.gz
```

### 4. 变异注释阶段
```bash
# 使用基因组注释进行变异注释
snpEff ann -v ARS-UCD1.2 variants.vcf.gz > annotated_variants.vcf
```

## 引用信息
Rosen, B.D., Bickhart, D.M., Schnabel, R.D. et al. 
De novo assembly of the cattle reference genome with single-molecule sequencing. 
GigaScience 9, giaa021 (2020).
