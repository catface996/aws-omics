# Docker镜像构建文件

本目录包含用于AWS Omics工作流的所有Docker镜像构建文件。这些Dockerfile基于项目中已成功执行的任务使用的镜像。

## 📋 镜像列表

### 数据预处理镜像

#### 1. FastQC v0.12.1
- **文件**: `Dockerfile.fastqc`
- **用途**: 基因组测序数据质量评估
- **ECR镜像**: `864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastqc:0.12.1`
- **摘要**: `sha256:86bc4b39fde28868e624bfcd2a153ebfa2fe0f0c206e78c686d820d9c4984dd6`

#### 2. fastp v0.23.4
- **文件**: `Dockerfile.fastp`
- **用途**: 高速测序数据预处理和质量控制
- **ECR镜像**: `864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastp:0.23.4`
- **摘要**: `sha256:8a33155af2ba1b32bb3b67707f88a061aaa35565f6bca5b6216453b85a0693c1`

#### 3. seqkit v2.5.1
- **文件**: `Dockerfile.seqkit`
- **用途**: FASTA/Q文件处理和序列去重
- **ECR镜像**: `864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/seqkit:2.5.1`
- **摘要**: `sha256:4a0a239e967c754149a6bf679ae7d514a418091f916b6e3d76959e1d731fa9fb`

### 序列比对镜像

#### 4. BWA v0.7.17
- **文件**: `Dockerfile.bwa`
- **用途**: 基因组序列比对
- **ECR镜像**: `864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bwa:0.7.17`
- **摘要**: `sha256:673b7c758b98a850e7a2face544089636fc948433647f65f9f9b31757712ff7f`

#### 5. SAMtools v1.17
- **文件**: `Dockerfile.samtools`
- **用途**: SAM/BAM文件处理和操作
- **ECR镜像**: `864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/samtools:1.17`
- **摘要**: `sha256:fb9237b5cf803fed459e366c36ed487d09cb2aa9a8a76853fe9554b2fe47d969`

#### 6. BWA + SAMtools 组合镜像
- **文件**: `Dockerfile.bwa-samtools`
- **用途**: 序列比对工作流，包含BWA和SAMtools
- **ECR镜像**: `864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bioinformatics:bwa-samtools-1.17`
- **摘要**: `sha256:c510aafa23bac3849bb20607efdd98aa225460180d346d8da2377763fb2fe478`

## 🚀 构建和推送命令

### 构建镜像
```bash
# 构建FastQC镜像
docker build -f Dockerfile.fastqc -t omics/fastqc:0.12.1 .

# 构建fastp镜像
docker build -f Dockerfile.fastp -t omics/fastp:0.23.4 .

# 构建seqkit镜像
docker build -f Dockerfile.seqkit -t omics/seqkit:2.5.1 .

# 构建BWA镜像
docker build -f Dockerfile.bwa -t omics/bwa:0.7.17 .

# 构建SAMtools镜像
docker build -f Dockerfile.samtools -t omics/samtools:1.17 .

# 构建BWA+SAMtools组合镜像
docker build -f Dockerfile.bwa-samtools -t omics/bioinformatics:bwa-samtools-1.17 .
```

### 标记和推送到ECR
```bash
# 获取ECR登录令牌
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin 864899854573.dkr.ecr.us-east-1.amazonaws.com

# 标记并推送FastQC
docker tag omics/fastqc:0.12.1 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastqc:0.12.1
docker push 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastqc:0.12.1

# 标记并推送fastp
docker tag omics/fastp:0.23.4 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastp:0.23.4
docker push 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastp:0.23.4

# 标记并推送seqkit
docker tag omics/seqkit:2.5.1 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/seqkit:2.5.1
docker push 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/seqkit:2.5.1

# 标记并推送BWA
docker tag omics/bwa:0.7.17 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bwa:0.7.17
docker push 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bwa:0.7.17

# 标记并推送SAMtools
docker tag omics/samtools:1.17 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/samtools:1.17
docker push 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/samtools:1.17

# 标记并推送BWA+SAMtools组合镜像
docker tag omics/bioinformatics:bwa-samtools-1.17 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bioinformatics:bwa-samtools-1.17
docker push 864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bioinformatics:bwa-samtools-1.17
```

## 📊 使用状态

### ✅ 已成功使用的镜像
这些镜像都已在AWS Omics工作流中成功使用：

1. **数据预处理工作流** (运行ID: 5969296) - ✅ 已完成
   - FastQC: 初始和最终质量评估
   - fastp: 数据清洗和质量过滤
   - seqkit: 序列去重处理

2. **序列比对工作流** (运行ID: 1809888) - 🔄 正在运行
   - BWA+SAMtools组合镜像: 索引构建和序列比对

## 🔧 技术特性

### 镜像优化
- **基础镜像**: Ubuntu 20.04 LTS (稳定性和兼容性)
- **多阶段构建**: 减少最终镜像大小
- **依赖管理**: 精确的版本控制和依赖安装
- **验证机制**: 每个镜像都包含安装验证步骤

### 安全性
- **最小权限**: 只安装必要的依赖
- **清理机制**: 构建后清理临时文件和缓存
- **版本固定**: 使用特定版本避免不确定性

## 📚 相关文档

- [Docker最佳实践](https://docs.docker.com/develop/dev-best-practices/)
- [AWS ECR用户指南](https://docs.aws.amazon.com/ecr/latest/userguide/)
- [AWS Omics容器镜像](https://docs.aws.amazon.com/omics/latest/dev/workflows-containers.html)

---

**最后更新**: 2025年8月10日  
**状态**: 所有镜像已成功构建并推送到ECR  
**用途**: AWS Omics基因组分析工作流
