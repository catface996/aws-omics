#!/bin/bash

# AWS Omics Read Set 预处理工作流运行脚本
# 直接处理Sequence Store中的Read Set数据

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 脚本目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CONFIG_FILE="$PROJECT_ROOT/scripts/00-setup/omics_stores_config.json"

echo -e "${BLUE}🧬 AWS Omics Read Set 预处理工作流${NC}"
echo "=================================================="

# 检查配置文件
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo -e "${RED}❌ 配置文件不存在: $CONFIG_FILE${NC}"
    echo "请先运行 scripts/00-setup/00-setup_omics_environment.sh"
    exit 1
fi

# 读取配置
REGION=$(jq -r '.region' "$CONFIG_FILE")
REFERENCE_STORE_ID=$(jq -r '.reference_store.id' "$CONFIG_FILE")
REFERENCE_ID=$(jq -r '.reference_store.reference_id' "$CONFIG_FILE")
SEQUENCE_STORE_ID=$(jq -r '.sequence_store.id' "$CONFIG_FILE")
READSET_ID=$(jq -r '.sequence_store.readset_id' "$CONFIG_FILE")
SAMPLE_ID=$(jq -r '.sequence_store.sample_id' "$CONFIG_FILE")
ROLE_ARN=$(jq -r '.role_arn' "$CONFIG_FILE")

echo "区域: $REGION"
echo "Reference Store: $REFERENCE_STORE_ID"
echo "Sequence Store: $SEQUENCE_STORE_ID"
echo "Read Set ID: $READSET_ID"
echo "样本ID: $SAMPLE_ID"
echo ""

# 检查Read Set是否存在
echo "🔍 验证Read Set..."
if aws omics get-read-set-metadata \
    --sequence-store-id "$SEQUENCE_STORE_ID" \
    --id "$READSET_ID" \
    --region "$REGION" >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Read Set验证成功${NC}"
else
    echo -e "${RED}❌ Read Set不存在或无法访问${NC}"
    echo "Read Set ID: $READSET_ID"
    echo "Sequence Store ID: $SEQUENCE_STORE_ID"
    exit 1
fi

# 检查参考基因组
echo "🔍 验证参考基因组..."
if aws omics get-reference-metadata \
    --reference-store-id "$REFERENCE_STORE_ID" \
    --id "$REFERENCE_ID" \
    --region "$REGION" >/dev/null 2>&1; then
    echo -e "${GREEN}✅ 参考基因组验证成功${NC}"
else
    echo -e "${RED}❌ 参考基因组不存在或无法访问${NC}"
    echo "Reference ID: $REFERENCE_ID"
    echo "Reference Store ID: $REFERENCE_STORE_ID"
    exit 1
fi

# 检查工作流是否存在
WORKFLOW_NAME="cow-preprocessing-readset"
echo ""
echo "🔍 检查预处理工作流..."

WORKFLOW_ID=$(aws omics list-workflows \
    --region "$REGION" \
    --query "items[?name=='$WORKFLOW_NAME'].id" \
    --output text 2>/dev/null || echo "")

if [[ -z "$WORKFLOW_ID" || "$WORKFLOW_ID" == "None" ]]; then
    echo -e "${YELLOW}⚠️ 工作流不存在，正在创建...${NC}"
    
    # 创建工作流
    WORKFLOW_RESPONSE=$(aws omics create-workflow \
        --name "$WORKFLOW_NAME" \
        --description "奶牛基因组Read Set预处理工作流，包含质量控制、接头去除、质量过滤和去重复" \
        --definition-zip "fileb://$PROJECT_ROOT/workflows/01_data_preprocessing/preprocessing_workflow_readset.wdl" \
        --engine WDL \
        --region "$REGION" \
        --tags Project=CowGenomics,Type=Preprocessing,Version=ReadSet \
        --output json)
    
    WORKFLOW_ID=$(echo "$WORKFLOW_RESPONSE" | jq -r '.id')
    echo -e "${GREEN}✅ 工作流创建成功${NC}"
    echo "   工作流ID: $WORKFLOW_ID"
else
    echo -e "${GREEN}✅ 工作流已存在${NC}"
    echo "   工作流ID: $WORKFLOW_ID"
fi

# 构建Read Set ARN
READSET_ARN="arn:aws:omics:$REGION:$(aws sts get-caller-identity --query Account --output text):sequenceStore/$SEQUENCE_STORE_ID/readSet/$READSET_ID"
REFERENCE_ARN="arn:aws:omics:$REGION:$(aws sts get-caller-identity --query Account --output text):referenceStore/$REFERENCE_STORE_ID/reference/$REFERENCE_ID"

# 准备工作流参数
echo ""
echo "🚀 准备启动预处理工作流..."
echo "Read Set ARN: $READSET_ARN"
echo "Reference ARN: $REFERENCE_ARN"

# 创建参数文件
PARAMS_FILE="/tmp/readset_preprocessing_params.json"
cat > "$PARAMS_FILE" << EOF
{
  "PreprocessingWorkflowReadSet.sample_name": "$SAMPLE_ID",
  "PreprocessingWorkflowReadSet.input_readset": "$READSET_ARN",
  "PreprocessingWorkflowReadSet.reference_genome": "$REFERENCE_ARN",
  "PreprocessingWorkflowReadSet.min_length": 50,
  "PreprocessingWorkflowReadSet.min_quality": 20,
  "PreprocessingWorkflowReadSet.threads": 8,
  "PreprocessingWorkflowReadSet.paired_end": true,
  "PreprocessingWorkflowReadSet.max_length": 500,
  "PreprocessingWorkflowReadSet.complexity_threshold": 30,
  "PreprocessingWorkflowReadSet.enable_polyg_trimming": true,
  "PreprocessingWorkflowReadSet.enable_polyx_trimming": true,
  "PreprocessingWorkflowReadSet.dedup_method": "fastuniq",
  "PreprocessingWorkflowReadSet.fastqc_memory_gb": 8,
  "PreprocessingWorkflowReadSet.trimmomatic_memory_gb": 16,
  "PreprocessingWorkflowReadSet.fastp_memory_gb": 16,
  "PreprocessingWorkflowReadSet.dedup_memory_gb": 16,
  "PreprocessingWorkflowReadSet.multiqc_memory_gb": 8
}
EOF

echo "参数文件已创建: $PARAMS_FILE"

# 启动工作流
echo ""
echo "🚀 启动预处理工作流..."

RUN_NAME="cow-preprocessing-readset-$(date +%Y%m%d-%H%M%S)"
OUTPUT_URI="s3://catface996-genomic/omics-outputs/preprocessing-readset/"

RUN_RESPONSE=$(aws omics start-run \
    --workflow-id "$WORKFLOW_ID" \
    --role-arn "$ROLE_ARN" \
    --name "$RUN_NAME" \
    --parameters "file://$PARAMS_FILE" \
    --output-uri "$OUTPUT_URI" \
    --storage-capacity 1200 \
    --log-level ALL \
    --tags Project=CowGenomics,Type=Preprocessing,Sample="$SAMPLE_ID",DataSource=ReadSet \
    --region "$REGION" \
    --output json)

RUN_ID=$(echo "$RUN_RESPONSE" | jq -r '.id')

echo -e "${GREEN}🎉 工作流启动成功！${NC}"
echo "=================================================="
echo "运行ID: $RUN_ID"
echo "运行名称: $RUN_NAME"
echo "工作流ID: $WORKFLOW_ID"
echo "输出位置: $OUTPUT_URI$RUN_ID"
echo ""
echo "🔗 AWS控制台链接:"
echo "https://console.aws.amazon.com/omics/home?region=$REGION#/runs/$RUN_ID"
echo ""
echo "📋 监控命令:"
echo "aws omics get-run --id $RUN_ID --region $REGION"
echo ""
echo "⏳ 预计运行时间: 30-60分钟"

# 清理临时文件
rm -f "$PARAMS_FILE"

# 保存运行信息
RUN_INFO_FILE="$PROJECT_ROOT/outputs/preprocessing/readset_run_info.json"
mkdir -p "$(dirname "$RUN_INFO_FILE")"

cat > "$RUN_INFO_FILE" << EOF
{
  "run_id": "$RUN_ID",
  "run_name": "$RUN_NAME",
  "workflow_id": "$WORKFLOW_ID",
  "readset_arn": "$READSET_ARN",
  "reference_arn": "$REFERENCE_ARN",
  "sample_id": "$SAMPLE_ID",
  "output_uri": "$OUTPUT_URI$RUN_ID",
  "region": "$REGION",
  "started_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "console_url": "https://console.aws.amazon.com/omics/home?region=$REGION#/runs/$RUN_ID"
}
EOF

echo "运行信息已保存到: $RUN_INFO_FILE"
