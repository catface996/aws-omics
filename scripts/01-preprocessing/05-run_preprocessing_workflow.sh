#!/bin/bash

# AWS Omics 数据预处理工作流运行脚本
# 用于在AWS Omics上运行完整的测序数据预处理流水线

set -euo pipefail

echo "🧬 AWS Omics 数据预处理工作流启动器"
echo "=================================================="

# 配置参数
WORKFLOW_NAME="cow-genome-preprocessing"
WORKFLOW_VERSION="1.0"
REGION="us-east-1"
ROLE_ARN="arn:aws:iam::$(aws sts get-caller-identity --query Account --output text):role/OmicsWorkflowRole"

# 文件路径
WORKFLOW_DIR="workflows/01_data_preprocessing"
WORKFLOW_FILE="$WORKFLOW_DIR/preprocessing_workflow.wdl"
INPUTS_DIR="$WORKFLOW_DIR/inputs"

# 检查必要文件
echo "🔍 检查工作流文件..."
if [[ ! -f "$WORKFLOW_FILE" ]]; then
    echo "❌ 工作流文件不存在: $WORKFLOW_FILE"
    exit 1
fi

echo "✅ 工作流文件: $WORKFLOW_FILE"

# 显示菜单
echo ""
echo "选择运行模式:"
echo "1. 双端测序数据预处理 (默认)"
echo "2. 单端测序数据预处理"
echo "3. 批量处理多个样本"
echo "4. 创建/更新工作流定义"
echo "5. 查看现有工作流"
echo "6. 查看运行历史"

read -p "请选择 (1-6) [默认: 1]: " -n 1 -r
echo

# 如果没有输入，默认选择1
if [[ -z "$REPLY" ]]; then
    REPLY="1"
fi

case $REPLY in
    1)
        echo "🔬 双端测序数据预处理模式"
        INPUTS_FILE="$INPUTS_DIR/preprocessing_inputs.json"
        RUN_NAME="cow-preprocessing-pe-$(date +%Y%m%d-%H%M%S)"
        ;;
    2)
        echo "🔬 单端测序数据预处理模式"
        INPUTS_FILE="$INPUTS_DIR/preprocessing_inputs_single_end.json"
        RUN_NAME="cow-preprocessing-se-$(date +%Y%m%d-%H%M%S)"
        ;;
    3)
        echo "🔬 批量处理模式"
        echo "请确保已准备好批量输入文件"
        read -p "输入批量配置文件路径: " BATCH_FILE
        if [[ ! -f "$BATCH_FILE" ]]; then
            echo "❌ 批量配置文件不存在: $BATCH_FILE"
            exit 1
        fi
        INPUTS_FILE="$BATCH_FILE"
        RUN_NAME="cow-preprocessing-batch-$(date +%Y%m%d-%H%M%S)"
        ;;
    4)
        echo "🛠️ 创建/更新工作流定义"
        create_or_update_workflow
        exit 0
        ;;
    5)
        echo "📋 查看现有工作流"
        aws omics list-workflows --region $REGION --output table
        exit 0
        ;;
    6)
        echo "📊 查看运行历史"
        aws omics list-runs --region $REGION --output table
        exit 0
        ;;
    *)
        echo "❌ 无效选择"
        exit 1
        ;;
esac

# 创建或更新工作流定义的函数
create_or_update_workflow() {
    echo "🛠️ 创建/更新工作流定义..."
    
    # 检查工作流是否已存在
    EXISTING_WORKFLOW=$(aws omics list-workflows \
        --region $REGION \
        --query "items[?name=='$WORKFLOW_NAME'].id" \
        --output text 2>/dev/null || echo "")
    
    if [[ -n "$EXISTING_WORKFLOW" && "$EXISTING_WORKFLOW" != "None" ]]; then
        echo "📝 更新现有工作流: $EXISTING_WORKFLOW"
        
        # 创建新版本
        aws omics create-workflow \
            --region $REGION \
            --name "${WORKFLOW_NAME}-v$(date +%Y%m%d-%H%M%S)" \
            --description "奶牛基因组数据预处理工作流 - 更新版本" \
            --engine WDL \
            --definition-zip fileb://<(cd $WORKFLOW_DIR && zip -r - .) \
            --parameter-template file://$INPUTS_FILE \
            --tags Project=CowGenomics,Version=$WORKFLOW_VERSION,Type=Preprocessing
    else
        echo "🆕 创建新工作流定义..."
        
        # 打包工作流文件
        echo "📦 打包工作流文件..."
        cd $WORKFLOW_DIR
        zip -r ../preprocessing_workflow.zip .
        cd - > /dev/null
        
        # 创建工作流
        WORKFLOW_ID=$(aws omics create-workflow \
            --region $REGION \
            --name $WORKFLOW_NAME \
            --description "奶牛基因组测序数据预处理工作流，包含质量评估、接头去除、质量过滤、长度过滤和去重复处理" \
            --engine WDL \
            --definition-zip fileb://$WORKFLOW_DIR/../preprocessing_workflow.zip \
            --parameter-template file://$INPUTS_FILE \
            --tags Project=CowGenomics,Version=$WORKFLOW_VERSION,Type=Preprocessing \
            --query 'id' \
            --output text)
        
        echo "✅ 工作流创建成功: $WORKFLOW_ID"
        
        # 清理临时文件
        rm -f $WORKFLOW_DIR/../preprocessing_workflow.zip
    fi
}

# 检查输入文件
if [[ ! -f "$INPUTS_FILE" ]]; then
    echo "❌ 输入参数文件不存在: $INPUTS_FILE"
    echo "请先编辑输入参数文件，设置正确的S3路径和参数"
    exit 1
fi

echo "✅ 输入参数文件: $INPUTS_FILE"

# 显示输入参数预览
echo ""
echo "📋 输入参数预览:"
echo "----------------------------------------"
head -20 "$INPUTS_FILE"
echo "----------------------------------------"

# 支持自动确认模式
if [[ "${AUTO_CONFIRM:-false}" == "true" ]]; then
    echo "🤖 自动确认模式，跳过确认提示"
    REPLY="y"
else
    read -p "确认运行工作流? (Y/n) [默认: Y]: " -n 1 -r
    echo
    # 如果没有输入或输入为空，默认为y
    if [[ -z "$REPLY" ]]; then
        REPLY="y"
    fi
fi

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ 取消运行"
    exit 1
fi

# 确保工作流存在
echo "🔍 检查工作流定义..."
WORKFLOW_ID=$(aws omics list-workflows \
    --region $REGION \
    --query "items[?name=='$WORKFLOW_NAME'].id" \
    --output text 2>/dev/null || echo "")

if [[ -z "$WORKFLOW_ID" || "$WORKFLOW_ID" == "None" ]]; then
    echo "⚠️ 工作流不存在，正在创建..."
    create_or_update_workflow
    
    # 重新获取工作流ID
    WORKFLOW_ID=$(aws omics list-workflows \
        --region $REGION \
        --query "items[?name=='$WORKFLOW_NAME'].id" \
        --output text)
fi

echo "✅ 使用工作流: $WORKFLOW_ID"

# 创建运行组（如果不存在）
echo "🏗️ 检查运行组..."
RUN_GROUP_ID=$(aws omics list-run-groups \
    --region $REGION \
    --query "items[?name=='cow-genomics-preprocessing'].id" \
    --output text 2>/dev/null || echo "")

if [[ -z "$RUN_GROUP_ID" || "$RUN_GROUP_ID" == "None" ]]; then
    echo "🆕 创建运行组..."
    RUN_GROUP_ID=$(aws omics create-run-group \
        --region $REGION \
        --name "cow-genomics-preprocessing" \
        --max-cpus 256 \
        --max-runs 10 \
        --max-duration 7200 \
        --tags Project=CowGenomics,Type=Preprocessing \
        --query 'id' \
        --output text)
    echo "✅ 运行组创建成功: $RUN_GROUP_ID"
else
    echo "✅ 使用现有运行组: $RUN_GROUP_ID"
fi

# 启动工作流运行
echo ""
echo "🚀 启动工作流运行..."
echo "运行名称: $RUN_NAME"
echo "工作流ID: $WORKFLOW_ID"
echo "运行组ID: $RUN_GROUP_ID"

# 创建运行时参数文件（使用预定义的参数值）
RUNTIME_PARAMS=$(mktemp)
echo "🔄 准备运行参数..."
cat > "$RUNTIME_PARAMS" << EOF
{
  "PreprocessingWorkflow.sample_name": "SRR16760538",
  "PreprocessingWorkflow.fastq_r1": "s3://catface996-genomic/genomic_data/01_raw_data/SRR16760538_1.fastq.gz",
  "PreprocessingWorkflow.fastq_r2": "s3://catface996-genomic/genomic_data/01_raw_data/SRR16760538_2.fastq.gz",
  "PreprocessingWorkflow.reference_genome": "s3://catface996-genomic/genomic_data/02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna",
  "PreprocessingWorkflow.min_length": 50,
  "PreprocessingWorkflow.min_quality": 20,
  "PreprocessingWorkflow.threads": 8,
  "PreprocessingWorkflow.paired_end": true,
  "PreprocessingWorkflow.RunFastPPE.max_length": 500,
  "PreprocessingWorkflow.RunFastPPE.complexity_threshold": 30,
  "PreprocessingWorkflow.RunFastPPE.enable_polyg_trimming": true,
  "PreprocessingWorkflow.RunFastPPE.enable_polyx_trimming": true,
  "PreprocessingWorkflow.RemoveDuplicatesPE.method": "fastuniq",
  "PreprocessingWorkflow.RunFastQC.memory_gb": 8,
  "PreprocessingWorkflow.RunTrimmomaticPE.memory_gb": 16,
  "PreprocessingWorkflow.RunFastPPE.memory_gb": 16,
  "PreprocessingWorkflow.RemoveDuplicatesPE.memory_gb": 16,
  "PreprocessingWorkflow.RunMultiQC.memory_gb": 8
}
EOF

RUN_ID=$(aws omics start-run \
    --region $REGION \
    --workflow-id $WORKFLOW_ID \
    --workflow-type PRIVATE \
    --run-group-id $RUN_GROUP_ID \
    --name $RUN_NAME \
    --role-arn "arn:aws:iam::$(aws sts get-caller-identity --query Account --output text):role/OmicsServiceRole" \
    --parameters file://$RUNTIME_PARAMS \
    --output-uri "s3://catface996-genomic/omics-outputs/preprocessing/" \
    --log-level ALL \
    --tags Project=CowGenomics,Sample=SRR16760538,Type=Preprocessing \
    --query 'id' \
    --output text)

# 清理临时文件
rm -f "$RUNTIME_PARAMS"

echo ""
echo "🎉 工作流运行已启动!"
echo "运行ID: $RUN_ID"
echo ""
echo "📊 监控命令:"
echo "aws omics get-run --region $REGION --id $RUN_ID"
echo ""
echo "📋 查看日志:"
echo "aws omics list-run-tasks --region $REGION --id $RUN_ID"
echo ""
echo "🌐 AWS控制台:"
echo "https://console.aws.amazon.com/omics/home?region=$REGION#/runs/$RUN_ID"

# 可选：等待运行完成
read -p "是否等待运行完成? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "⏳ 等待运行完成..."
    
    while true; do
        STATUS=$(aws omics get-run \
            --region $REGION \
            --id $RUN_ID \
            --query 'status' \
            --output text)
        
        echo "当前状态: $STATUS ($(date))"
        
        case $STATUS in
            "COMPLETED")
                echo "🎉 运行成功完成!"
                break
                ;;
            "FAILED"|"CANCELLED")
                echo "❌ 运行失败或被取消"
                echo "查看详细信息:"
                aws omics get-run --region $REGION --id $RUN_ID
                exit 1
                ;;
            "RUNNING"|"PENDING"|"STARTING")
                sleep 30
                ;;
            *)
                echo "⚠️ 未知状态: $STATUS"
                sleep 30
                ;;
        esac
    done
    
    # 显示运行结果
    echo ""
    echo "📊 运行结果:"
    aws omics get-run --region $REGION --id $RUN_ID --output table
fi

echo ""
echo "✅ 脚本执行完成!"
