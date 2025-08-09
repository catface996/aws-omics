#!/bin/bash

# 奶牛基因组数据下载脚本
# 基于项目文档中的数据源信息
# 作者: AWS Omics Demo
# 日期: $(date +%Y-%m-%d)

set -e

# 配置变量
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.1_ARS-UCD1.2"
DOWNLOAD_DIR="./genomic_data"
LOG_FILE="download.log"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"
}

log_progress() {
    echo -e "${BLUE}[PROGRESS]${NC} $1" | tee -a "$LOG_FILE"
}

# 检查依赖
check_dependencies() {
    log_info "检查系统依赖..."
    
    if ! command -v curl &> /dev/null; then
        log_error "curl未安装，请先安装curl"
        exit 1
    fi
    
    if ! command -v md5sum &> /dev/null && ! command -v md5 &> /dev/null; then
        log_warn "MD5校验工具未找到，将跳过文件完整性验证"
    fi
    
    log_info "依赖检查完成"
}

# 创建下载目录
create_directories() {
    log_info "创建下载目录结构..."
    
    mkdir -p "$DOWNLOAD_DIR"/{reference,annotation,protein,rna,features,assembly,checksums}
    
    log_info "目录结构创建完成"
}

# 下载文件函数
download_file() {
    local url="$1"
    local filename="$2"
    local description="$3"
    local target_dir="$4"
    
    log_progress "下载 $description..."
    
    if [ -f "$DOWNLOAD_DIR/$target_dir/$filename" ]; then
        log_warn "$filename 已存在，跳过下载"
        return 0
    fi
    
    if curl -L --fail --show-error --progress-bar \
        -o "$DOWNLOAD_DIR/$target_dir/$filename" \
        "$url/$filename"; then
        log_info "$description 下载完成"
        return 0
    else
        log_error "$description 下载失败"
        return 1
    fi
}

# 下载核心基因组数据
download_core_data() {
    log_info "开始下载核心基因组数据..."
    
    # 1. 参考基因组序列 (最重要)
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz" \
        "参考基因组序列 (803MB)" \
        "reference"
    
    # 2. GTF注释文件
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz" \
        "GTF基因注释文件 (21MB)" \
        "annotation"
    
    # 3. GFF注释文件
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz" \
        "GFF基因注释文件 (23MB)" \
        "annotation"
    
    # 4. 蛋白质序列
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_protein.faa.gz" \
        "蛋白质序列文件 (15MB)" \
        "protein"
    
    # 5. 编码序列
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_cds_from_genomic.fna.gz" \
        "编码序列文件 (18MB)" \
        "annotation"
    
    # 6. 功能特征表
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_feature_table.txt.gz" \
        "功能特征表 (4MB)" \
        "features"
    
    log_info "核心基因组数据下载完成"
}

# 下载扩展数据
download_extended_data() {
    log_info "开始下载扩展数据..."
    
    # RNA序列数据
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_rna.fna.gz" \
        "RNA序列文件 (51MB)" \
        "rna"
    
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_rna_from_genomic.fna.gz" \
        "基因组提取RNA序列 (41MB)" \
        "rna"
    
    # 重复序列信息
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_rm.out.gz" \
        "重复序列掩码结果 (171MB)" \
        "annotation"
    
    log_info "扩展数据下载完成"
}

# 下载辅助文件
download_auxiliary_files() {
    log_info "开始下载辅助文件..."
    
    # 组装信息
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_assembly_report.txt" \
        "基因组组装报告" \
        "assembly"
    
    download_file "$BASE_URL" \
        "GCF_002263795.1_ARS-UCD1.2_assembly_stats.txt" \
        "组装统计信息" \
        "assembly"
    
    # MD5校验文件
    download_file "$BASE_URL" \
        "md5checksums.txt" \
        "MD5校验文件" \
        "checksums"
    
    log_info "辅助文件下载完成"
}

# 验证文件完整性
verify_integrity() {
    log_info "开始验证文件完整性..."
    
    if [ -f "$DOWNLOAD_DIR/checksums/md5checksums.txt" ]; then
        cd "$DOWNLOAD_DIR"
        
        # 提取相关文件的MD5值进行验证
        while IFS= read -r line; do
            if [[ $line == *"genomic.fna.gz"* ]] || \
               [[ $line == *"genomic.gtf.gz"* ]] || \
               [[ $line == *"protein.faa.gz"* ]]; then
                
                expected_md5=$(echo "$line" | awk '{print $1}')
                filename=$(echo "$line" | awk '{print $2}' | sed 's|^./||')
                
                # 查找文件在哪个子目录
                found_file=$(find . -name "$(basename "$filename")" -type f | head -1)
                
                if [ -n "$found_file" ]; then
                    if command -v md5sum &> /dev/null; then
                        actual_md5=$(md5sum "$found_file" | awk '{print $1}')
                    elif command -v md5 &> /dev/null; then
                        actual_md5=$(md5 -q "$found_file")
                    else
                        log_warn "无法验证 $filename 的完整性"
                        continue
                    fi
                    
                    if [ "$expected_md5" = "$actual_md5" ]; then
                        log_info "✓ $filename 完整性验证通过"
                    else
                        log_error "✗ $filename 完整性验证失败"
                    fi
                fi
            fi
        done < "checksums/md5checksums.txt"
        
        cd - > /dev/null
    else
        log_warn "MD5校验文件不存在，跳过完整性验证"
    fi
}

# 生成数据摘要
generate_summary() {
    log_info "生成数据下载摘要..."
    
    cat > "$DOWNLOAD_DIR/README.md" << EOF
# 奶牛基因组数据下载摘要

## 下载信息
- 下载时间: $(date)
- 数据源: NCBI Genome Database
- 基因组版本: ARS-UCD1.2
- 物种: Bos taurus (奶牛)

## 目录结构
\`\`\`
genomic_data/
├── reference/          # 参考基因组序列
├── annotation/         # 基因注释文件
├── protein/           # 蛋白质序列
├── rna/              # RNA序列数据
├── features/         # 功能特征表
├── assembly/         # 组装信息
└── checksums/        # 校验文件
\`\`\`

## 主要文件说明

### 参考基因组
- \`GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz\`: 参考基因组FASTA格式

### 基因注释
- \`GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz\`: GTF格式注释
- \`GCF_002263795.1_ARS-UCD1.2_genomic.gff.gz\`: GFF3格式注释

### 蛋白质数据
- \`GCF_002263795.1_ARS-UCD1.2_protein.faa.gz\`: 蛋白质序列

## 使用建议
1. 解压文件: \`gunzip *.gz\`
2. 建立索引: 使用samtools、bwa等工具为参考基因组建立索引
3. 质量控制: 验证文件完整性和格式正确性

## 引用信息
Rosen, B.D., Bickhart, D.M., Schnabel, R.D. et al. 
De novo assembly of the cattle reference genome with single-molecule sequencing. 
GigaScience 9, giaa021 (2020).
EOF

    # 计算总下载大小
    total_size=$(du -sh "$DOWNLOAD_DIR" | cut -f1)
    log_info "数据下载完成！总大小: $total_size"
    log_info "数据摘要已保存到: $DOWNLOAD_DIR/README.md"
}

# 主函数
main() {
    echo "=========================================="
    echo "    奶牛基因组数据下载工具"
    echo "=========================================="
    echo
    
    # 初始化日志
    echo "下载开始时间: $(date)" > "$LOG_FILE"
    
    # 执行下载流程
    check_dependencies
    create_directories
    
    log_info "开始下载奶牛基因组数据..."
    log_info "数据源: $BASE_URL"
    log_info "下载目录: $DOWNLOAD_DIR"
    
    download_core_data
    
    # 询问是否下载扩展数据
    echo
    read -p "是否下载扩展数据 (RNA序列、重复序列等)? [y/N]: " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        download_extended_data
    fi
    
    download_auxiliary_files
    verify_integrity
    generate_summary
    
    echo
    echo "=========================================="
    log_info "所有数据下载完成！"
    echo "=========================================="
    echo
    echo "下载目录: $DOWNLOAD_DIR"
    echo "日志文件: $LOG_FILE"
    echo "数据摘要: $DOWNLOAD_DIR/README.md"
    echo
}

# 脚本入口
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
