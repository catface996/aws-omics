#!/bin/bash

# AWS Omics 存储创建脚本
# 创建Reference Store和Sequence Store用于奶牛基因组项目

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 项目配置
PROJECT_NAME="cow-genomics"
REGION="us-east-1"
REFERENCE_STORE_NAME="cow-reference-genomes"
SEQUENCE_STORE_NAME="cow-genomics-sequences"

echo -e "${BLUE}🧬 AWS Omics 存储创建工具${NC}"
echo "=================================================="
echo "项目: $PROJECT_NAME"
echo "区域: $REGION"
echo ""

# 检查AWS CLI配置
echo "🔍 检查AWS CLI配置..."
if ! aws sts get-caller-identity >/dev/null 2>&1; then
    echo -e "${RED}❌ AWS CLI未配置或无权限${NC}"
    echo "请运行: aws configure"
    exit 1
fi

ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
echo -e "${GREEN}✅ AWS账号: $ACCOUNT_ID${NC}"

# 检查IAM角色
echo ""
echo "🔍 检查IAM服务角色..."
ROLE_NAME="OmicsServiceRole"
if aws iam get-role --role-name $ROLE_NAME >/dev/null 2>&1; then
    echo -e "${GREEN}✅ IAM角色已存在: $ROLE_NAME${NC}"
else
    echo -e "${YELLOW}⚠️ 创建IAM服务角色...${NC}"
    
    # 创建信任策略
    cat > /tmp/omics-trust-policy.json << EOF
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

    # 创建角色
    aws iam create-role \
        --role-name $ROLE_NAME \
        --assume-role-policy-document file:///tmp/omics-trust-policy.json \
        --description "Service role for AWS Omics workflows"
    
    # 附加权限策略
    aws iam attach-role-policy \
        --role-name $ROLE_NAME \
        --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
    
    aws iam attach-role-policy \
        --role-name $ROLE_NAME \
        --policy-arn arn:aws:iam::aws:policy/CloudWatchLogsFullAccess
    
    # 清理临时文件
    rm -f /tmp/omics-trust-policy.json
    
    echo -e "${GREEN}✅ IAM角色创建成功${NC}"
fi

ROLE_ARN="arn:aws:iam::$ACCOUNT_ID:role/$ROLE_NAME"

# 创建Reference Store
echo ""
echo "🏗️ 创建Reference Store..."
REFERENCE_STORE_RESPONSE=$(aws omics create-reference-store \
    --name $REFERENCE_STORE_NAME \
    --description "奶牛参考基因组存储，包含ARS-UCD1.2等标准参考基因组" \
    --region $REGION \
    --output json 2>/dev/null || echo "")

if [[ -n "$REFERENCE_STORE_RESPONSE" ]]; then
    REFERENCE_STORE_ID=$(echo "$REFERENCE_STORE_RESPONSE" | jq -r '.id')
    echo -e "${GREEN}✅ Reference Store创建成功${NC}"
    echo "   ID: $REFERENCE_STORE_ID"
    echo "   名称: $REFERENCE_STORE_NAME"
else
    # 检查是否已存在
    REFERENCE_STORE_ID=$(aws omics list-reference-stores \
        --region $REGION \
        --query "referenceStores[?name=='$REFERENCE_STORE_NAME'].id" \
        --output text 2>/dev/null || echo "")
    
    if [[ -n "$REFERENCE_STORE_ID" && "$REFERENCE_STORE_ID" != "None" ]]; then
        echo -e "${YELLOW}⚠️ Reference Store已存在${NC}"
        echo "   ID: $REFERENCE_STORE_ID"
    else
        echo -e "${RED}❌ Reference Store创建失败${NC}"
        exit 1
    fi
fi

# 创建Sequence Store
echo ""
echo "🏗️ 创建Sequence Store..."
SEQUENCE_STORE_RESPONSE=$(aws omics create-sequence-store \
    --name $SEQUENCE_STORE_NAME \
    --description "奶牛基因组测序数据存储，用于存储FASTQ文件和相关元数据" \
    --region $REGION \
    --output json 2>/dev/null || echo "")

if [[ -n "$SEQUENCE_STORE_RESPONSE" ]]; then
    SEQUENCE_STORE_ID=$(echo "$SEQUENCE_STORE_RESPONSE" | jq -r '.id')
    echo -e "${GREEN}✅ Sequence Store创建成功${NC}"
    echo "   ID: $SEQUENCE_STORE_ID"
    echo "   名称: $SEQUENCE_STORE_NAME"
    
    # 显示S3访问点信息
    S3_ACCESS_POINT=$(echo "$SEQUENCE_STORE_RESPONSE" | jq -r '.s3Access.s3AccessPointArn')
    echo "   S3访问点: $S3_ACCESS_POINT"
else
    # 检查是否已存在
    SEQUENCE_STORE_ID=$(aws omics list-sequence-stores \
        --region $REGION \
        --query "sequenceStores[?name=='$SEQUENCE_STORE_NAME'].id" \
        --output text 2>/dev/null || echo "")
    
    if [[ -n "$SEQUENCE_STORE_ID" && "$SEQUENCE_STORE_ID" != "None" ]]; then
        echo -e "${YELLOW}⚠️ Sequence Store已存在${NC}"
        echo "   ID: $SEQUENCE_STORE_ID"
    else
        echo -e "${RED}❌ Sequence Store创建失败${NC}"
        exit 1
    fi
fi

# 保存配置信息
echo ""
echo "💾 保存配置信息..."
CONFIG_FILE="$(dirname "$0")/omics_stores_config.json"
cat > "$CONFIG_FILE" << EOF
{
  "project_name": "$PROJECT_NAME",
  "region": "$REGION",
  "account_id": "$ACCOUNT_ID",
  "role_arn": "$ROLE_ARN",
  "reference_store": {
    "id": "$REFERENCE_STORE_ID",
    "name": "$REFERENCE_STORE_NAME"
  },
  "sequence_store": {
    "id": "$SEQUENCE_STORE_ID",
    "name": "$SEQUENCE_STORE_NAME"
  },
  "created_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF

echo -e "${GREEN}✅ 配置已保存到: $CONFIG_FILE${NC}"

# 显示总结
echo ""
echo -e "${BLUE}📊 创建总结${NC}"
echo "=================================================="
echo "Reference Store ID: $REFERENCE_STORE_ID"
echo "Sequence Store ID:  $SEQUENCE_STORE_ID"
echo "IAM角色ARN:        $ROLE_ARN"
echo ""
echo -e "${GREEN}✅ AWS Omics存储创建完成！${NC}"
echo ""
echo "🔗 AWS控制台链接:"
echo "Reference Store: https://console.aws.amazon.com/omics/home?region=$REGION#/reference-stores/$REFERENCE_STORE_ID"
echo "Sequence Store:  https://console.aws.amazon.com/omics/home?region=$REGION#/sequence-stores/$SEQUENCE_STORE_ID"
echo ""
echo "📋 下一步:"
echo "1. 运行 02-import_reference_genome.sh 导入参考基因组"
echo "2. 运行 03-import_sequencing_data.sh 导入测序数据"
