#!/bin/bash

# AWS Omics 测序数据导入脚本
# 将奶牛基因组测序数据导入到Sequence Store

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

echo -e "${BLUE}🧬 AWS Omics 测序数据导入工具${NC}"
echo "=================================================="

# 检查配置文件
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo -e "${RED}❌ 配置文件不存在: $CONFIG_FILE${NC}"
    echo "请先运行 01-create_omics_stores.sh"
    exit 1
fi

# 读取配置
REGION=$(jq -r '.region' "$CONFIG_FILE")
SEQUENCE_STORE_ID=$(jq -r '.sequence_store.id' "$CONFIG_FILE")
ROLE_ARN=$(jq -r '.role_arn' "$CONFIG_FILE")

echo "区域: $REGION"
echo "Sequence Store ID: $SEQUENCE_STORE_ID"
echo ""

# 显示数据源选项
echo "选择测序数据源:"
echo "1. 本地S3存储的双端数据 (推荐用于测试)"
echo "2. 本地S3存储的单端数据"
echo "3. NCBI SRA原始双端数据 (推荐用于生产)"
echo ""

read -p "请选择 (1-3): " -n 1 -r
echo

case $REPLY in
    1)
        echo "🔬 使用本地S3双端测序数据"
        DATA_SOURCE="local-paired"
        SAMPLE_ID="SRR16760538-local-paired"
        SOURCE_FILE_1="s3://catface996-genomic/genomic_data/01_raw_data/SRR16760538_1.fastq.gz"
        SOURCE_FILE_2="s3://catface996-genomic/genomic_data/01_raw_data/SRR16760538_2.fastq.gz"
        DESCRIPTION="奶牛基因组双端测序数据，来自本地S3存储"
        ;;
    2)
        echo "🔬 使用本地S3单端测序数据"
        DATA_SOURCE="local-single"
        SAMPLE_ID="SRR16760538-local-single"
        SOURCE_FILE_1="s3://catface996-genomic/genomic_data/01_raw_data/SRR16760538.fastq.gz"
        SOURCE_FILE_2=""
        DESCRIPTION="奶牛基因组单端测序数据，来自本地S3存储"
        ;;
    3)
        echo "🔬 使用NCBI SRA原始双端数据"
        DATA_SOURCE="ncbi-paired"
        SAMPLE_ID="SRR16760538-ncbi-original"
        SOURCE_FILE_1="s3://sra-pub-src-4/SRR16760538/X4M_clean_1.fq.gz.1"
        SOURCE_FILE_2="s3://sra-pub-src-4/SRR16760538/X4M_clean_2.fq.gz.1"
        DESCRIPTION="奶牛基因组双端测序数据，来自NCBI SRA原始数据"
        ;;
    *)
        echo "❌ 无效选择"
        exit 1
        ;;
esac

SUBJECT_ID="cow-sample-001"

# 检查现有Read Sets
echo ""
echo "🔍 检查现有Read Sets..."
EXISTING_READSETS=$(aws omics list-read-sets \
    --sequence-store-id $SEQUENCE_STORE_ID \
    --region $REGION \
    --query "readSets[?sampleId=='$SAMPLE_ID'].id" \
    --output text 2>/dev/null || echo "")

if [[ -n "$EXISTING_READSETS" && "$EXISTING_READSETS" != "None" ]]; then
    echo -e "${YELLOW}⚠️ Read Set已存在: $SAMPLE_ID${NC}"
    echo "   ID: $EXISTING_READSETS"
    read -p "是否重新导入? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "跳过导入"
        exit 0
    fi
fi

# 准备导入参数
echo ""
echo "🚀 开始导入测序数据..."
echo "样本ID: $SAMPLE_ID"
echo "主体ID: $SUBJECT_ID"
echo "描述: $DESCRIPTION"

if [[ -n "$SOURCE_FILE_2" ]]; then
    echo "源文件1: $SOURCE_FILE_1"
    echo "源文件2: $SOURCE_FILE_2"
    echo "类型: 双端测序"
    
    IMPORT_JOB_RESPONSE=$(aws omics start-read-set-import-job \
        --sequence-store-id $SEQUENCE_STORE_ID \
        --role-arn $ROLE_ARN \
        --sources sourceFileType=FASTQ,sourceFiles="{source1=$SOURCE_FILE_1,source2=$SOURCE_FILE_2}",sampleId="$SAMPLE_ID",subjectId="$SUBJECT_ID",description="$DESCRIPTION" \
        --region $REGION \
        --output json)
else
    echo "源文件: $SOURCE_FILE_1"
    echo "类型: 单端测序"
    
    IMPORT_JOB_RESPONSE=$(aws omics start-read-set-import-job \
        --sequence-store-id $SEQUENCE_STORE_ID \
        --role-arn $ROLE_ARN \
        --sources sourceFileType=FASTQ,sourceFiles="{source1=$SOURCE_FILE_1}",sampleId="$SAMPLE_ID",subjectId="$SUBJECT_ID",description="$DESCRIPTION" \
        --region $REGION \
        --output json)
fi

IMPORT_JOB_ID=$(echo "$IMPORT_JOB_RESPONSE" | jq -r '.id')
echo -e "${GREEN}✅ 导入任务已启动${NC}"
echo "   任务ID: $IMPORT_JOB_ID"

# 监控导入进度
echo ""
echo "⏳ 监控导入进度..."
START_TIME=$(date +%s)

while true; do
    STATUS_RESPONSE=$(aws omics get-read-set-import-job \
        --sequence-store-id $SEQUENCE_STORE_ID \
        --id $IMPORT_JOB_ID \
        --region $REGION \
        --output json)
    
    STATUS=$(echo "$STATUS_RESPONSE" | jq -r '.status')
    STATUS_MESSAGE=$(echo "$STATUS_RESPONSE" | jq -r '.statusMessage')
    
    # 计算运行时间
    CURRENT_TIME=$(date +%s)
    ELAPSED=$((CURRENT_TIME - START_TIME))
    ELAPSED_MIN=$((ELAPSED / 60))
    ELAPSED_SEC=$((ELAPSED % 60))
    
    case $STATUS in
        "SUBMITTED"|"IN_PROGRESS")
            echo "   状态: $STATUS - 运行时间: ${ELAPSED_MIN}m${ELAPSED_SEC}s"
            sleep 30
            ;;
        "COMPLETED")
            echo -e "${GREEN}🎉 导入成功完成！${NC}"
            
            # 获取Read Set ID
            READ_SET_ID=$(echo "$STATUS_RESPONSE" | jq -r '.sources[0].readSetId // empty')
            COMPLETION_TIME=$(echo "$STATUS_RESPONSE" | jq -r '.completionTime')
            
            if [[ -n "$READ_SET_ID" ]]; then
                echo "   Read Set ID: $READ_SET_ID"
            fi
            echo "   完成时间: $COMPLETION_TIME"
            echo "   总运行时间: ${ELAPSED_MIN}m${ELAPSED_SEC}s"
            
            # 更新配置文件
            jq --arg data_source "$DATA_SOURCE" --arg sample_id "$SAMPLE_ID" --arg readset_id "$READ_SET_ID" --arg import_job_id "$IMPORT_JOB_ID" \
               '.sequence_store.data_source = $data_source | .sequence_store.sample_id = $sample_id | .sequence_store.readset_id = $readset_id | .sequence_store.import_job_id = $import_job_id' \
               "$CONFIG_FILE" > "$CONFIG_FILE.tmp" && mv "$CONFIG_FILE.tmp" "$CONFIG_FILE"
            
            break
            ;;
        "COMPLETED_WITH_FAILURES")
            echo -e "${RED}❌ 导入完成但有失败${NC}"
            echo "   状态: $STATUS"
            echo "   消息: $STATUS_MESSAGE"
            
            # 显示详细错误
            SOURCE_STATUS=$(echo "$STATUS_RESPONSE" | jq -r '.sources[0].status')
            SOURCE_MESSAGE=$(echo "$STATUS_RESPONSE" | jq -r '.sources[0].statusMessage')
            echo "   源文件状态: $SOURCE_STATUS"
            echo "   源文件错误: $SOURCE_MESSAGE"
            
            echo ""
            echo "💡 常见解决方案:"
            echo "1. 检查FASTQ文件格式是否标准"
            echo "2. 尝试使用NCBI原始数据 (选项3)"
            echo "3. 检查文件是否损坏"
            exit 1
            ;;
        "FAILED"|"CANCELLED")
            echo -e "${RED}❌ 导入失败${NC}"
            echo "   状态: $STATUS"
            echo "   错误信息: $STATUS_MESSAGE"
            exit 1
            ;;
        *)
            echo -e "${YELLOW}⚠️ 未知状态: $STATUS${NC}"
            sleep 30
            ;;
    esac
done

# 显示Read Set信息
if [[ -n "$READ_SET_ID" ]]; then
    echo ""
    echo "📊 Read Set详细信息:"
    aws omics get-read-set-metadata \
        --sequence-store-id $SEQUENCE_STORE_ID \
        --id $READ_SET_ID \
        --region $REGION \
        --output table
fi

# 显示总结
echo ""
echo -e "${BLUE}📊 导入总结${NC}"
echo "=================================================="
echo "样本ID:         $SAMPLE_ID"
echo "Read Set ID:    $READ_SET_ID"
echo "导入任务ID:     $IMPORT_JOB_ID"
echo "Sequence Store: $SEQUENCE_STORE_ID"
echo "数据源:         $DATA_SOURCE"
echo ""
echo -e "${GREEN}✅ 测序数据导入完成！${NC}"
echo ""
echo "🔗 AWS控制台链接:"
echo "https://console.aws.amazon.com/omics/home?region=$REGION#/sequence-stores/$SEQUENCE_STORE_ID"
echo ""
echo "📋 下一步:"
echo "现在可以在工作流中使用Read Set ID: $READ_SET_ID"
