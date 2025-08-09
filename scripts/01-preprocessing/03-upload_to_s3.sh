#!/bin/bash

# S3上传脚本 - 上传模拟测序数据到 catface996-genomic bucket
echo "🚀 开始上传模拟测序数据到 S3..."

# 设置变量
BUCKET_NAME="catface996-genomic"
REGION="us-east-1"
LOCAL_DATA_DIR="outputs/mutated_genome"

# 检查AWS CLI配置
echo "🔧 检查AWS配置..."
aws sts get-caller-identity --region $REGION

if [ $? -ne 0 ]; then
    echo "❌ AWS CLI未正确配置，请先配置AWS凭证"
    exit 1
fi

echo ""
echo "📁 准备上传的数据:"
echo "  本地目录: $LOCAL_DATA_DIR"
echo "  S3 Bucket: s3://$BUCKET_NAME"
echo "  区域: $REGION"

# 检查本地数据
if [ ! -d "$LOCAL_DATA_DIR" ]; then
    echo "❌ 本地数据目录不存在: $LOCAL_DATA_DIR"
    exit 1
fi

echo ""
echo "📊 数据大小统计:"
du -sh $LOCAL_DATA_DIR

echo ""
echo "🔄 开始上传..."

# 1. 上传完整基因组文件
echo ""
echo "1️⃣ 上传完整基因组文件..."
aws s3 cp "$LOCAL_DATA_DIR/complete_mutated_genome.fna" \
    "s3://$BUCKET_NAME/reference-genome/complete_mutated_genome.fna" \
    --region $REGION \
    --storage-class STANDARD \
    --metadata "data-type=reference-genome,source=gpu-intensive-simulation,chromosomes=30"

if [ $? -eq 0 ]; then
    echo "✅ 完整基因组文件上传成功"
else
    echo "❌ 完整基因组文件上传失败"
fi

# 2. 上传分散染色体文件
echo ""
echo "2️⃣ 上传分散染色体文件..."
aws s3 sync "$LOCAL_DATA_DIR/complete_gpu_intensive/" \
    "s3://$BUCKET_NAME/chromosomes/" \
    --region $REGION \
    --storage-class STANDARD \
    --exclude "*.json" \
    --metadata "data-type=chromosome-files,source=gpu-intensive-simulation"

if [ $? -eq 0 ]; then
    echo "✅ 分散染色体文件上传成功"
else
    echo "❌ 分散染色体文件上传失败"
fi

# 3. 上传报告文件
echo ""
echo "3️⃣ 上传报告文件..."
aws s3 cp "$LOCAL_DATA_DIR/complete_mutated_genome_merge_report.json" \
    "s3://$BUCKET_NAME/reports/complete_mutated_genome_merge_report.json" \
    --region $REGION \
    --storage-class STANDARD \
    --metadata "data-type=report,content=merge-validation"

aws s3 cp "$LOCAL_DATA_DIR/complete_gpu_intensive/gpu_intensive_generation_report.json" \
    "s3://$BUCKET_NAME/reports/gpu_intensive_generation_report.json" \
    --region $REGION \
    --storage-class STANDARD \
    --metadata "data-type=report,content=generation-stats"

if [ $? -eq 0 ]; then
    echo "✅ 报告文件上传成功"
else
    echo "❌ 报告文件上传失败"
fi

echo ""
echo "📋 验证上传结果..."
aws s3 ls "s3://$BUCKET_NAME/" --region $REGION --recursive --human-readable --summarize

echo ""
echo "🎉 上传完成!"
echo ""
echo "📁 S3 目录结构:"
echo "s3://$BUCKET_NAME/"
echo "├── reference-genome/"
echo "│   └── complete_mutated_genome.fna      # 完整基因组 (2.5GB)"
echo "├── chromosomes/"
echo "│   ├── NC_037328.1_mutated.fna         # 染色体文件 x30"
echo "│   └── ..."
echo "└── reports/"
echo "    ├── complete_mutated_genome_merge_report.json"
echo "    └── gpu_intensive_generation_report.json"
echo ""
echo "🔗 访问链接:"
echo "  AWS Console: https://s3.console.aws.amazon.com/s3/buckets/$BUCKET_NAME"
echo "  CLI访问: aws s3 ls s3://$BUCKET_NAME/ --region $REGION"
