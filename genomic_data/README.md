# 奶牛基因组数据下载摘要

## 下载信息
- 下载时间: Fri  8 Aug 2025 10:51:19 CST
- 数据源: NCBI Genome Database
- 基因组版本: ARS-UCD1.2
- 物种: Bos taurus (奶牛)

## 目录结构
```
genomic_data/
├── reference/          # 参考基因组序列
├── annotation/         # 基因注释文件
├── protein/           # 蛋白质序列
├── rna/              # RNA序列数据
├── features/         # 功能特征表
├── assembly/         # 组装信息
└── checksums/        # 校验文件
```

## 主要文件说明

### 参考基因组
- `GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz`: 参考基因组FASTA格式

### 基因注释
- `GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz`: GTF格式注释
- `GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz`: GFF3格式注释

### 蛋白质数据
- `GCF_002263795.1_ARS-UCD1.2_protein.faa.gz`: 蛋白质序列

## 使用建议
1. 解压文件: `gunzip *.gz`
2. 建立索引: 使用samtools、bwa等工具为参考基因组建立索引
3. 质量控制: 验证文件完整性和格式正确性

## 引用信息
Rosen, B.D., Bickhart, D.M., Schnabel, R.D. et al. 
De novo assembly of the cattle reference genome with single-molecule sequencing. 
GigaScience 9, giaa021 (2020).
