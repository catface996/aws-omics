#!/bin/bash

# 修复版工作流v9部署脚本
# 解决FastQC Perl模块依赖问题

set -e

echo "=== 部署修复版工作流v9 ==="

# 配置变量
WORKFLOW_NAME="cow-preprocessing-readset-v9-fixed"
WORKFLOW_FILE="workflows/01_data_preprocessing/preprocessing_readset_v9_fixed.wdl"
REGION="us-east-1"
OUTPUT_BUCKET="s3://catface996-genomic/omics-outputs/preprocessing-readset-v9-fixed/"

# 检查工作流文件是否存在
if [ ! -f "$WORKFLOW_FILE" ]; then
    echo "错误: 工作流文件 $WORKFLOW_FILE 不存在"
    exit 1
fi

echo "1. 创建工作流ZIP包..."
cd workflows/01_data_preprocessing
zip -q preprocessing_readset_v9_fixed.zip preprocessing_readset_v9_fixed.wdl
cd ../..

echo "2. 创建AWS Omics工作流..."
WORKFLOW_ID=$(aws omics create-workflow \
    --region $REGION \
    --name "$WORKFLOW_NAME" \
    --description "修复版奶牛基因组预处理工作流v9 - 使用公共生物信息学镜像解决Perl依赖问题" \
    --engine WDL \
    --definition-zip fileb://workflows/01_data_preprocessing/preprocessing_readset_v9_fixed.zip \
    --main main.wdl \
    --query 'id' \
    --output text)

echo "工作流创建成功，ID: $WORKFLOW_ID"

echo "3. 准备运行参数..."
cat > /tmp/run_params_v9.json << EOF
{
    "sample_name": "SRR16760538",
    "input_fastq": "s3://864899854573-4118666-ij76xyhmqz6rwbszcqsw5wb7m655huse1b-s3alias/864899854573/sequenceStore/4118666176/readSet/5825575116/SRR16760538.fastq.gz",
    "fastqc_memory_gb": 8
}
EOF

echo "4. 启动工作流运行..."
RUN_NAME="cow-preprocessing-readset-v9-fixed-$(date +%Y%m%d-%H%M%S)"

RUN_ID=$(aws omics start-run \
    --region $REGION \
    --workflow-id $WORKFLOW_ID \
    --name "$RUN_NAME" \
    --role-arn "arn:aws:iam::864899854573:role/OmicsServiceRole" \
    --parameters file:///tmp/run_params_v9.json \
    --output-uri "$OUTPUT_BUCKET" \
    --storage-capacity 1200 \
    --query 'id' \
    --output text)

echo "工作流运行启动成功！"
echo "运行ID: $RUN_ID"
echo "运行名称: $RUN_NAME"
echo "输出位置: $OUTPUT_BUCKET$RUN_ID"

echo "5. 监控运行状态..."
echo "使用以下命令监控运行状态:"
echo "aws omics get-run --region $REGION --id $RUN_ID"

# 清理临时文件
rm -f /tmp/run_params_v9.json

echo "=== 部署完成 ==="
