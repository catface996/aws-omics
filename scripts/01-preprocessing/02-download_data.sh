#!/bin/bash

# 数据下载脚本
# 用于下载项目所需的大型数据文件

set -e

echo "开始下载基因组数据..."

# 创建必要的目录
mkdir -p genomic_data/reference
mkdir -p genomic_data/annotation
mkdir -p genomic_data/protein
mkdir -p data/public_fastq

# 下载参考基因组
echo "下载参考基因组..."
if [ ! -f "genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz" ]; then
    wget -O genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz"
fi

# 下载注释文件
echo "下载注释文件..."
if [ ! -f "genomic_data/annotation/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz" ]; then
    wget -O genomic_data/annotation/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2/GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz"
fi

# 解压参考基因组（如果需要）
if [ ! -f "genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna" ]; then
    echo "解压参考基因组..."
    gunzip -k genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz
fi

echo "数据下载完成！"
echo "注意：大型FASTQ文件需要从NCBI SRA单独下载"
echo "使用命令：fastq-dump --split-files SRR16760538"
