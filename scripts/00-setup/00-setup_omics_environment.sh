#!/bin/bash

# AWS Omics 环境完整设置脚本
# 一键创建Reference Store、Sequence Store并导入数据

set -e

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 脚本目录
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "${BLUE}🧬 AWS Omics 环境完整设置${NC}"
echo "=================================================="
echo "这个脚本将完成以下操作:"
echo "1. 创建IAM服务角色"
echo "2. 创建Reference Store和Sequence Store"
echo "3. 导入奶牛参考基因组 ARS-UCD1.2"
echo "4. 导入测序数据到Sequence Store"
echo ""

read -p "确认开始设置? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ 取消设置"
    exit 1
fi

# 步骤1: 创建存储
echo ""
echo -e "${BLUE}步骤 1/3: 创建AWS Omics存储${NC}"
echo "=================================================="
if [[ -x "$SCRIPT_DIR/01-create_omics_stores.sh" ]]; then
    "$SCRIPT_DIR/01-create_omics_stores.sh"
else
    echo -e "${RED}❌ 脚本不存在或不可执行: 01-create_omics_stores.sh${NC}"
    exit 1
fi

# 步骤2: 导入参考基因组
echo ""
echo -e "${BLUE}步骤 2/3: 导入参考基因组${NC}"
echo "=================================================="
if [[ -x "$SCRIPT_DIR/02-import_reference_genome.sh" ]]; then
    "$SCRIPT_DIR/02-import_reference_genome.sh"
else
    echo -e "${RED}❌ 脚本不存在或不可执行: 02-import_reference_genome.sh${NC}"
    exit 1
fi

# 步骤3: 导入测序数据
echo ""
echo -e "${BLUE}步骤 3/3: 导入测序数据${NC}"
echo "=================================================="
if [[ -x "$SCRIPT_DIR/03-import_sequencing_data.sh" ]]; then
    "$SCRIPT_DIR/03-import_sequencing_data.sh"
else
    echo -e "${RED}❌ 脚本不存在或不可执行: 03-import_sequencing_data.sh${NC}"
    exit 1
fi

# 完成总结
echo ""
echo -e "${GREEN}🎉 AWS Omics环境设置完成！${NC}"
echo "=================================================="

# 读取最终配置
CONFIG_FILE="$SCRIPT_DIR/omics_stores_config.json"
if [[ -f "$CONFIG_FILE" ]]; then
    REGION=$(jq -r '.region' "$CONFIG_FILE")
    REFERENCE_STORE_ID=$(jq -r '.reference_store.id' "$CONFIG_FILE")
    SEQUENCE_STORE_ID=$(jq -r '.sequence_store.id' "$CONFIG_FILE")
    REFERENCE_ID=$(jq -r '.reference_store.reference_id // "N/A"' "$CONFIG_FILE")
    READSET_ID=$(jq -r '.sequence_store.readset_id // "N/A"' "$CONFIG_FILE")
    
    echo "✅ Reference Store: $REFERENCE_STORE_ID"
    echo "✅ Sequence Store:  $SEQUENCE_STORE_ID"
    echo "✅ 参考基因组ID:    $REFERENCE_ID"
    echo "✅ Read Set ID:     $READSET_ID"
    echo ""
    echo "🔗 AWS控制台链接:"
    echo "Reference Store: https://console.aws.amazon.com/omics/home?region=$REGION#/reference-stores/$REFERENCE_STORE_ID"
    echo "Sequence Store:  https://console.aws.amazon.com/omics/home?region=$REGION#/sequence-stores/$SEQUENCE_STORE_ID"
fi

echo ""
echo "📋 下一步建议:"
echo "1. 在工作流中使用新的Reference Store和Sequence Store"
echo "2. 更新工作流参数以引用Read Set ID"
echo "3. 享受AWS Omics原生存储的性能和成本优势"
