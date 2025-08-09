#!/bin/bash

# AWS Omics 奶牛基因组分析部署脚本
# 使用模拟测序数据

set -e

# 配置变量
PROJECT_NAME="cattle-genomics-simulation"
REGION="us-east-1"
BUCKET_NAME="cattle-genomics-sim-$(date +%s)"
WORKFLOW_NAME="cattle-variant-calling"
REFERENCE_STORE_NAME="cattle-reference-store"
SEQUENCE_STORE_NAME="cattle-sequence-store"

# 颜色输出
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 检查AWS CLI配置
check_aws_config() {
    log_info "检查AWS CLI配置..."
    
    if ! command -v aws &> /dev/null; then
        log_error "AWS CLI未安装"
        exit 1
    fi
    
    if ! aws sts get-caller-identity &> /dev/null; then
        log_error "AWS凭证未配置"
        exit 1
    fi
    
    ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
    log_info "当前AWS账户: $ACCOUNT_ID"
}

# 创建S3存储桶
create_s3_bucket() {
    log_info "创建S3存储桶: $BUCKET_NAME"
    
    aws s3 mb s3://$BUCKET_NAME --region $REGION
    
    # 设置生命周期策略
    cat > lifecycle.json << EOF
{
    "Rules": [
        {
            "ID": "GenomicsDataLifecycle",
            "Status": "Enabled",
            "Transitions": [
                {
                    "Days": 30,
                    "StorageClass": "STANDARD_IA"
                },
                {
                    "Days": 90,
                    "StorageClass": "GLACIER"
                }
            ]
        }
    ]
}
EOF
    
    aws s3api put-bucket-lifecycle-configuration \
        --bucket $BUCKET_NAME \
        --lifecycle-configuration file://lifecycle.json
    
    rm lifecycle.json
}

# 上传数据到S3
upload_data() {
    log_info "上传基因组数据到S3..."
    
    # 上传参考基因组
    aws s3 cp genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz \
        s3://$BUCKET_NAME/reference/genome.fna.gz
    
    # 上传模拟测序数据
    if [ -f "simulated_data/simulated_reads.fastq" ]; then
        aws s3 cp simulated_data/simulated_reads.fastq \
            s3://$BUCKET_NAME/reads/sample1.fastq
    else
        log_warn "模拟测序数据不存在，请先运行 simulate_sequencing_data.py"
    fi
    
    # 上传工作流定义
    aws s3 cp variant_calling_workflow.wdl \
        s3://$BUCKET_NAME/workflows/variant_calling.wdl
}

# 创建IAM角色
create_iam_role() {
    log_info "创建IAM角色..."
    
    # Omics服务角色
    cat > omics-trust-policy.json << EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Principal": {
                "Service": "omics.amazonaws.com"
            },
            "Action": "sts:AssumeRole"
        }
    ]
}
EOF

    aws iam create-role \
        --role-name OmicsExecutionRole \
        --assume-role-policy-document file://omics-trust-policy.json || true
    
    # 附加策略
    aws iam attach-role-policy \
        --role-name OmicsExecutionRole \
        --policy-arn arn:aws:iam::aws:policy/service-role/AmazonOmicsServiceRolePolicy || true
    
    rm omics-trust-policy.json
}

# 创建Omics资源
create_omics_resources() {
    log_info "创建Omics资源..."
    
    # 创建参考存储
    REFERENCE_STORE_ID=$(aws omics create-reference-store \
        --name $REFERENCE_STORE_NAME \
        --description "Cattle genome reference store" \
        --query 'id' --output text)
    
    log_info "参考存储ID: $REFERENCE_STORE_ID"
    
    # 创建序列存储
    SEQUENCE_STORE_ID=$(aws omics create-sequence-store \
        --name $SEQUENCE_STORE_NAME \
        --description "Cattle sequence data store" \
        --query 'id' --output text)
    
    log_info "序列存储ID: $SEQUENCE_STORE_ID"
    
    # 导入参考基因组
    aws omics create-reference \
        --reference-store-id $REFERENCE_STORE_ID \
        --name "ARS-UCD1.2" \
        --source-uri s3://$BUCKET_NAME/reference/genome.fna.gz
}

# 主函数
main() {
    echo "=========================================="
    echo "  AWS Omics 奶牛基因组分析部署"
    echo "=========================================="
    
    check_aws_config
    create_s3_bucket
    upload_data
    create_iam_role
    create_omics_resources
    
    log_info "部署完成！"
    echo "S3存储桶: $BUCKET_NAME"
    echo "参考存储: $REFERENCE_STORE_ID"
    echo "序列存储: $SEQUENCE_STORE_ID"
}

main "$@"
