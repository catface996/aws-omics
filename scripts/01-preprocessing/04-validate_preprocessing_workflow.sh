#!/bin/bash

# 验证预处理工作流WDL语法的脚本

set -euo pipefail

echo "🔍 验证预处理工作流WDL语法"
echo "=================================="

WORKFLOW_DIR="workflows/wdl/preprocessing"
MAIN_WORKFLOW="$WORKFLOW_DIR/preprocessing_workflow.wdl"

# 检查文件是否存在
echo "📁 检查工作流文件..."
if [[ ! -f "$MAIN_WORKFLOW" ]]; then
    echo "❌ 主工作流文件不存在: $MAIN_WORKFLOW"
    exit 1
fi

echo "✅ 主工作流文件: $MAIN_WORKFLOW"

# 检查任务文件
TASK_FILES=(
    "$WORKFLOW_DIR/tasks/fastqc_task.wdl"
    "$WORKFLOW_DIR/tasks/trimmomatic_task.wdl"
    "$WORKFLOW_DIR/tasks/fastp_task.wdl"
    "$WORKFLOW_DIR/tasks/deduplication_task.wdl"
    "$WORKFLOW_DIR/tasks/multiqc_task.wdl"
)

echo ""
echo "📋 检查任务文件..."
for task_file in "${TASK_FILES[@]}"; do
    if [[ -f "$task_file" ]]; then
        echo "✅ $(basename $task_file)"
    else
        echo "❌ 缺失: $(basename $task_file)"
        exit 1
    fi
done

# 检查输入文件
INPUT_FILES=(
    "$WORKFLOW_DIR/inputs/preprocessing_inputs.json"
    "$WORKFLOW_DIR/inputs/preprocessing_inputs_single_end.json"
)

echo ""
echo "📝 检查输入参数文件..."
for input_file in "${INPUT_FILES[@]}"; do
    if [[ -f "$input_file" ]]; then
        echo "✅ $(basename $input_file)"
        # 验证JSON语法
        if jq empty "$input_file" 2>/dev/null; then
            echo "  ✅ JSON语法正确"
        else
            echo "  ❌ JSON语法错误"
            exit 1
        fi
    else
        echo "❌ 缺失: $(basename $input_file)"
        exit 1
    fi
done

# 基本WDL语法检查
echo ""
echo "🔍 基本WDL语法检查..."

# 检查版本声明
if grep -q "version 1.0" "$MAIN_WORKFLOW"; then
    echo "✅ WDL版本声明正确"
else
    echo "❌ 缺少或错误的WDL版本声明"
    exit 1
fi

# 检查工作流定义
if grep -q "workflow PreprocessingWorkflow" "$MAIN_WORKFLOW"; then
    echo "✅ 工作流定义正确"
else
    echo "❌ 工作流定义错误"
    exit 1
fi

# 检查导入语句
echo ""
echo "📦 检查导入语句..."
imports=(
    "import \"tasks/fastqc_task.wdl\""
    "import \"tasks/trimmomatic_task.wdl\""
    "import \"tasks/fastp_task.wdl\""
    "import \"tasks/deduplication_task.wdl\""
    "import \"tasks/multiqc_task.wdl\""
)

for import_stmt in "${imports[@]}"; do
    if grep -q "$import_stmt" "$MAIN_WORKFLOW"; then
        echo "✅ $(echo $import_stmt | cut -d'"' -f2)"
    else
        echo "❌ 缺失导入: $import_stmt"
        exit 1
    fi
done

# 检查任务文件中的任务定义
echo ""
echo "🔧 检查任务定义..."
task_checks=(
    "tasks/fastqc_task.wdl:task RunFastQC"
    "tasks/trimmomatic_task.wdl:task RunTrimmomaticPE"
    "tasks/trimmomatic_task.wdl:task RunTrimmomaticSE"
    "tasks/fastp_task.wdl:task RunFastPPE"
    "tasks/fastp_task.wdl:task RunFastPSE"
    "tasks/deduplication_task.wdl:task RemoveDuplicatesPE"
    "tasks/deduplication_task.wdl:task RemoveDuplicatesSE"
    "tasks/multiqc_task.wdl:task RunMultiQC"
)

for check in "${task_checks[@]}"; do
    file=$(echo $check | cut -d':' -f1)
    task=$(echo $check | cut -d':' -f2)
    
    if grep -q "$task" "$WORKFLOW_DIR/$file"; then
        echo "✅ $task in $(basename $file)"
    else
        echo "❌ 缺失任务: $task in $file"
        exit 1
    fi
done

# 检查Docker镜像引用
echo ""
echo "🐳 检查Docker镜像引用..."
docker_images=(
    "quay.io/biocontainers/fastqc"
    "quay.io/biocontainers/trimmomatic"
    "quay.io/biocontainers/fastp"
    "quay.io/biocontainers/multiqc"
)

for image in "${docker_images[@]}"; do
    if grep -r "$image" "$WORKFLOW_DIR/tasks/" >/dev/null; then
        echo "✅ $image"
    else
        echo "⚠️  未找到: $image (可能使用了不同的镜像)"
    fi
done

# 生成工作流摘要
echo ""
echo "📊 工作流摘要"
echo "=================================="
echo "主工作流: $(basename $MAIN_WORKFLOW)"
echo "任务数量: ${#TASK_FILES[@]}"
echo "输入配置: ${#INPUT_FILES[@]}"

# 统计代码行数
total_lines=0
for file in "$MAIN_WORKFLOW" "${TASK_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        lines=$(wc -l < "$file")
        total_lines=$((total_lines + lines))
        echo "  $(basename $file): $lines 行"
    fi
done
echo "总代码行数: $total_lines"

# 检查是否可以创建工作流包
echo ""
echo "📦 创建工作流包测试..."
cd "$WORKFLOW_DIR"
if zip -r ../preprocessing_workflow_test.zip . >/dev/null 2>&1; then
    package_size=$(ls -lh ../preprocessing_workflow_test.zip | awk '{print $5}')
    echo "✅ 工作流包创建成功: $package_size"
    rm -f ../preprocessing_workflow_test.zip
else
    echo "❌ 工作流包创建失败"
    exit 1
fi
cd - >/dev/null

echo ""
echo "🎉 工作流验证完成!"
echo ""
echo "📋 下一步操作:"
echo "1. 编辑输入参数文件，设置正确的S3路径"
echo "2. 运行: ./scripts/run_preprocessing_workflow.sh"
echo "3. 或直接使用AWS CLI创建工作流"
echo ""
echo "💡 提示:"
echo "- 确保AWS凭证已配置"
echo "- 确保有足够的IAM权限"
echo "- 检查S3存储桶访问权限"
