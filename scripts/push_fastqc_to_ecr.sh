#!/bin/bash

# 将公共FastQC镜像推送到私有ECR
# 处理Mac到Linux的架构兼容性问题

set -e

echo "=== 推送FastQC镜像到私有ECR ==="

# 配置变量
REGION="us-east-1"
ACCOUNT_ID="864899854573"
ECR_REPO="omics/fastqc"
PUBLIC_IMAGE="quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
PRIVATE_IMAGE="$ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com/$ECR_REPO:0.12.1"

echo "1. 登录到ECR..."
aws ecr get-login-password --region $REGION | docker login --username AWS --password-stdin $ACCOUNT_ID.dkr.ecr.$REGION.amazonaws.com

echo "2. 检查ECR仓库是否存在..."
if ! aws ecr describe-repositories --region $REGION --repository-names $ECR_REPO >/dev/null 2>&1; then
    echo "创建ECR仓库: $ECR_REPO"
    aws ecr create-repository --region $REGION --repository-name $ECR_REPO
else
    echo "ECR仓库已存在: $ECR_REPO"
fi

echo "3. 拉取公共镜像（指定Linux/AMD64架构）..."
# 使用--platform参数确保拉取Linux/AMD64版本
docker pull --platform linux/amd64 $PUBLIC_IMAGE

echo "4. 重新标记镜像..."
docker tag $PUBLIC_IMAGE $PRIVATE_IMAGE

echo "5. 推送到私有ECR..."
docker push $PRIVATE_IMAGE

echo "6. 验证推送结果..."
aws ecr describe-images --region $REGION --repository-name $ECR_REPO --query 'imageDetails[0].{Digest:imageDigest,Size:imageSizeInBytes,PushedAt:imagePushedAt}'

echo "7. 清理本地镜像（可选）..."
read -p "是否删除本地镜像？(y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    docker rmi $PUBLIC_IMAGE $PRIVATE_IMAGE
    echo "本地镜像已删除"
fi

echo "=== 推送完成 ==="
echo "私有ECR镜像地址: $PRIVATE_IMAGE"
echo ""
echo "现在可以在WDL工作流中使用:"
echo "docker: \"$PRIVATE_IMAGE\""
