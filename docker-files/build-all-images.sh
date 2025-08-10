#!/bin/bash

# AWS Omics Docker镜像构建脚本
# 用于构建和推送所有生物信息学工具镜像到ECR

set -e

# 配置变量
ECR_REGISTRY="864899854573.dkr.ecr.us-east-1.amazonaws.com"
REGION="us-east-1"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 检查Docker是否运行
check_docker() {
    if ! docker info >/dev/null 2>&1; then
        log_error "Docker is not running. Please start Docker first."
        exit 1
    fi
    log_success "Docker is running"
}

# ECR登录
ecr_login() {
    log_info "Logging into ECR..."
    aws ecr get-login-password --region $REGION | docker login --username AWS --password-stdin $ECR_REGISTRY
    if [ $? -eq 0 ]; then
        log_success "Successfully logged into ECR"
    else
        log_error "Failed to login to ECR"
        exit 1
    fi
}

# 构建镜像函数
build_image() {
    local dockerfile=$1
    local image_name=$2
    local tag=$3
    
    log_info "Building $image_name:$tag..."
    docker build -f $dockerfile -t $image_name:$tag .
    
    if [ $? -eq 0 ]; then
        log_success "Successfully built $image_name:$tag"
        return 0
    else
        log_error "Failed to build $image_name:$tag"
        return 1
    fi
}

# 推送镜像函数
push_image() {
    local local_image=$1
    local ecr_image=$2
    
    log_info "Tagging and pushing $local_image to $ecr_image..."
    
    # 标记镜像
    docker tag $local_image $ecr_image
    
    # 推送镜像
    docker push $ecr_image
    
    if [ $? -eq 0 ]; then
        log_success "Successfully pushed $ecr_image"
        return 0
    else
        log_error "Failed to push $ecr_image"
        return 1
    fi
}

# 主函数
main() {
    log_info "Starting AWS Omics Docker image build process..."
    
    # 检查Docker
    check_docker
    
    # ECR登录
    ecr_login
    
    # 定义镜像列表
    declare -A images=(
        ["Dockerfile.fastqc"]="omics/fastqc:0.12.1"
        ["Dockerfile.fastp"]="omics/fastp:0.23.4"
        ["Dockerfile.seqkit"]="omics/seqkit:2.5.1"
        ["Dockerfile.bwa"]="omics/bwa:0.7.17"
        ["Dockerfile.samtools"]="omics/samtools:1.17"
        ["Dockerfile.bwa-samtools"]="omics/bioinformatics:bwa-samtools-1.17"
    )
    
    # 构建计数器
    local built=0
    local failed=0
    local total=${#images[@]}
    
    # 构建所有镜像
    for dockerfile in "${!images[@]}"; do
        local image_tag="${images[$dockerfile]}"
        local image_name=$(echo $image_tag | cut -d':' -f1)
        local tag=$(echo $image_tag | cut -d':' -f2)
        
        if build_image $dockerfile $image_name $tag; then
            ((built++))
            
            # 推送到ECR
            local ecr_image="$ECR_REGISTRY/$image_tag"
            if push_image $image_tag $ecr_image; then
                log_success "✅ $image_tag -> ECR"
            else
                log_error "❌ Failed to push $image_tag"
                ((failed++))
            fi
        else
            log_error "❌ Failed to build $image_tag"
            ((failed++))
        fi
        
        echo "----------------------------------------"
    done
    
    # 总结
    log_info "Build Summary:"
    log_info "Total images: $total"
    log_success "Successfully built and pushed: $built"
    if [ $failed -gt 0 ]; then
        log_error "Failed: $failed"
        exit 1
    else
        log_success "🎉 All images built and pushed successfully!"
    fi
}

# 显示帮助
show_help() {
    echo "AWS Omics Docker Image Build Script"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help     Show this help message"
    echo "  --build-only   Only build images, don't push to ECR"
    echo "  --push-only    Only push existing images to ECR"
    echo ""
    echo "Examples:"
    echo "  $0                 # Build and push all images"
    echo "  $0 --build-only    # Only build images locally"
    echo "  $0 --push-only     # Only push existing images"
}

# 解析命令行参数
case "${1:-}" in
    -h|--help)
        show_help
        exit 0
        ;;
    --build-only)
        log_info "Build-only mode enabled"
        # TODO: Implement build-only mode
        ;;
    --push-only)
        log_info "Push-only mode enabled"
        # TODO: Implement push-only mode
        ;;
    "")
        main
        ;;
    *)
        log_error "Unknown option: $1"
        show_help
        exit 1
        ;;
esac
