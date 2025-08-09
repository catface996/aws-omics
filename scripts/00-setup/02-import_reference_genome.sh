#!/bin/bash

# AWS Omics 参考基因组导入脚本
# 将奶牛参考基因组ARS-UCD1.2导入到Reference Store

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 脚本目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/omics_stores_config.json"

echo -e "${BLUE}🧬 AWS Omics 参考基因组导入工具${NC}"
echo "=================================================="

# 检查配置文件
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo -e "${RED}❌ 配置文件不存在: $CONFIG_FILE${NC}"
    echo "请先运行 01-create_omics_stores.sh"
    exit 1
fi

# 读取配置
REGION=$(jq -r '.region' "$CONFIG_FILE")
REFERENCE_STORE_ID=$(jq -r '.reference_store.id' "$CONFIG_FILE")
ROLE_ARN=$(jq -r '.role_arn' "$CONFIG_FILE")

echo "区域: $REGION"
echo "Reference Store ID: $REFERENCE_STORE_ID"
echo ""

# 参考基因组配置
REFERENCE_NAME="ARS-UCD1.2"
REFERENCE_DESCRIPTION="奶牛参考基因组 ARS-UCD1.2 版本，包含所有染色体和支架序列"
REFERENCE_S3_PATH="s3://catface996-genomic/genomic_data/02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna"

# 检查参考基因组文件是否存在
echo "🔍 检查参考基因组文件..."
if aws s3api head-object --bucket catface996-genomic --key genomic_data/02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna >/dev/null 2>&1; then
    echo -e "${GREEN}✅ 参考基因组文件存在${NC}"
    FILE_SIZE=$(aws s3api head-object --bucket catface996-genomic --key genomic_data/02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna --query ContentLength --output text)
    echo "   文件大小: $(numfmt --to=iec $FILE_SIZE)"
else
    echo -e "${RED}❌ 参考基因组文件不存在: $REFERENCE_S3_PATH${NC}"
    echo "请确保文件已上传到S3"
    exit 1
fi

# 检查是否已经导入
echo ""
echo "🔍 检查现有参考基因组..."
EXISTING_REFERENCES=$(aws omics list-references \
    --reference-store-id $REFERENCE_STORE_ID \
    --region $REGION \
    --query "references[?name=='$REFERENCE_NAME'].id" \
    --output text 2>/dev/null || echo "")

if [[ -n "$EXISTING_REFERENCES" && "$EXISTING_REFERENCES" != "None" ]]; then
    echo -e "${YELLOW}⚠️ 参考基因组已存在: $REFERENCE_NAME${NC}"
    echo "   ID: $EXISTING_REFERENCES"
    read -p "是否重新导入? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "跳过导入"
        exit 0
    fi
fi

# 开始导入
echo ""
echo "🚀 开始导入参考基因组..."
echo "名称: $REFERENCE_NAME"
echo "描述: $REFERENCE_DESCRIPTION"
echo "源文件: $REFERENCE_S3_PATH"

IMPORT_JOB_RESPONSE=$(aws omics start-reference-import-job \
    --reference-store-id $REFERENCE_STORE_ID \
    --role-arn $ROLE_ARN \
    --sources name="$REFERENCE_NAME",description="$REFERENCE_DESCRIPTION",sourceFile="$REFERENCE_S3_PATH" \
    --region $REGION \
    --output json)

IMPORT_JOB_ID=$(echo "$IMPORT_JOB_RESPONSE" | jq -r '.id')
echo -e "${GREEN}✅ 导入任务已启动${NC}"
echo "   任务ID: $IMPORT_JOB_ID"

# 监控导入进度
echo ""
echo "⏳ 监控导入进度..."
while true; do
    STATUS_RESPONSE=$(aws omics get-reference-import-job \
        --reference-store-id $REFERENCE_STORE_ID \
        --id $IMPORT_JOB_ID \
        --region $REGION \
        --output json)
    
    STATUS=$(echo "$STATUS_RESPONSE" | jq -r '.status')
    STATUS_MESSAGE=$(echo "$STATUS_RESPONSE" | jq -r '.statusMessage')
    
    case $STATUS in
        "SUBMITTED"|"IN_PROGRESS")
            echo "   状态: $STATUS - $STATUS_MESSAGE"
            sleep 30
            ;;
        "COMPLETED")
            echo -e "${GREEN}🎉 导入成功完成！${NC}"
            
            # 获取参考基因组ID
            REFERENCE_ID=$(echo "$STATUS_RESPONSE" | jq -r '.sources[0].referenceId')
            COMPLETION_TIME=$(echo "$STATUS_RESPONSE" | jq -r '.completionTime')
            
            echo "   参考基因组ID: $REFERENCE_ID"
            echo "   完成时间: $COMPLETION_TIME"
            
            # 更新配置文件
            jq --arg ref_id "$REFERENCE_ID" --arg import_job_id "$IMPORT_JOB_ID" \
               '.reference_store.reference_id = $ref_id | .reference_store.import_job_id = $import_job_id' \
               "$CONFIG_FILE" > "$CONFIG_FILE.tmp" && mv "$CONFIG_FILE.tmp" "$CONFIG_FILE"
            
            break
            ;;
        "FAILED"|"CANCELLED")
            echo -e "${RED}❌ 导入失败${NC}"
            echo "   状态: $STATUS"
            echo "   错误信息: $STATUS_MESSAGE"
            
            # 显示详细错误
            echo "$STATUS_RESPONSE" | jq '.sources[0]'
            exit 1
            ;;
        *)
            echo -e "${YELLOW}⚠️ 未知状态: $STATUS${NC}"
            sleep 30
            ;;
    esac
done

# 显示总结
echo ""
echo -e "${BLUE}📊 导入总结${NC}"
echo "=================================================="
echo "参考基因组名称: $REFERENCE_NAME"
echo "参考基因组ID:   $REFERENCE_ID"
echo "导入任务ID:     $IMPORT_JOB_ID"
echo "Reference Store: $REFERENCE_STORE_ID"
echo ""
echo -e "${GREEN}✅ 参考基因组导入完成！${NC}"
echo ""
echo "🔗 AWS控制台链接:"
echo "https://console.aws.amazon.com/omics/home?region=$REGION#/reference-stores/$REFERENCE_STORE_ID/references/$REFERENCE_ID"
echo ""
echo "📋 下一步:"
echo "运行 03-import_sequencing_data.sh 导入测序数据"
