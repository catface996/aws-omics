# AWS Omics 存储设置指南

本指南介绍如何为奶牛基因组项目设置AWS Omics的Reference Store和Sequence Store。

## 📋 概述

AWS Omics提供两种专用存储服务：
- **Reference Store**: 存储参考基因组，提供优化的索引和访问
- **Sequence Store**: 存储测序数据，提供高压缩率和快速访问

## 🎯 优势

### Reference Store
- ✅ 自动索引创建和维护
- ✅ 多工作流共享访问
- ✅ 版本管理支持
- ✅ 快速随机访问优化

### Sequence Store
- ✅ 比gzip高3-5倍的压缩率
- ✅ 显著降低存储成本
- ✅ 针对基因组分析优化的访问模式
- ✅ 丰富的元数据管理
- ✅ 与Omics工作流原生集成

## 🚀 快速开始

### 方法1: 一键设置（推荐）
```bash
cd scripts/00-setup
chmod +x *.sh
./00-setup_omics_environment.sh
```

### 方法2: 分步设置
```bash
cd scripts/00-setup
chmod +x *.sh

# 步骤1: 创建存储
./01-create_omics_stores.sh

# 步骤2: 导入参考基因组
./02-import_reference_genome.sh

# 步骤3: 导入测序数据
./03-import_sequencing_data.sh
```

## 📁 脚本说明

### 00-setup_omics_environment.sh
**功能**: 完整的环境设置脚本
- 调用所有子脚本
- 提供统一的设置体验
- 生成完整的配置总结

### 01-create_omics_stores.sh
**功能**: 创建AWS Omics存储
- 检查和创建IAM服务角色
- 创建Reference Store
- 创建Sequence Store
- 生成配置文件

**输出**:
- `omics_stores_config.json`: 存储配置信息

### 02-import_reference_genome.sh
**功能**: 导入参考基因组
- 导入奶牛参考基因组ARS-UCD1.2
- 监控导入进度
- 更新配置文件

**要求**:
- 参考基因组文件已上传到S3
- 路径: `s3://catface996-genomic/genomic_data/02_reference_genome/GCF_002263795.1_ARS-UCD1.2_genomic.fna`

### 03-import_sequencing_data.sh
**功能**: 导入测序数据
- 支持多种数据源选择
- 监控导入进度
- 生成Read Set

**数据源选项**:
1. **本地S3双端数据**: 适合测试
2. **本地S3单端数据**: 备选方案
3. **NCBI SRA原始数据**: 推荐用于生产（格式更标准）

## 📊 配置文件格式

`omics_stores_config.json` 包含以下信息：
```json
{
  "project_name": "cow-genomics",
  "region": "us-east-1",
  "account_id": "864899854573",
  "role_arn": "arn:aws:iam::864899854573:role/OmicsServiceRole",
  "reference_store": {
    "id": "4443887886",
    "name": "cow-reference-genomes",
    "reference_id": "6274639166",
    "import_job_id": "6758280994"
  },
  "sequence_store": {
    "id": "4118666176",
    "name": "cow-genomics-sequences",
    "data_source": "ncbi-paired",
    "sample_id": "SRR16760538-ncbi-original",
    "readset_id": "1234567890",
    "import_job_id": "7758569225"
  },
  "created_at": "2025-08-09T06:40:00Z"
}
```

## 🔧 故障排除

### 常见问题

#### 1. IAM权限不足
**症状**: 创建存储时权限错误
**解决**: 确保AWS CLI配置了足够权限的用户

#### 2. 参考基因组文件不存在
**症状**: 导入参考基因组时文件未找到
**解决**: 检查S3路径和文件是否存在

#### 3. 测序数据导入失败
**症状**: `File validation failed`
**解决**: 
- 尝试使用NCBI原始数据（选项3）
- 检查FASTQ文件格式
- 验证文件完整性

#### 4. 跨账号访问问题
**症状**: 访问其他账号S3 bucket失败
**解决**: 配置跨账号IAM权限和bucket policy

### 日志查看
```bash
# 查看导入任务详情
aws omics get-reference-import-job --reference-store-id <store-id> --id <job-id>
aws omics get-read-set-import-job --sequence-store-id <store-id> --id <job-id>

# 查看存储信息
aws omics get-reference-store --id <store-id>
aws omics get-sequence-store --id <store-id>
```

## 📈 性能优化

### 存储成本对比
| 存储方式 | 原始大小 | 压缩后大小 | 成本节省 |
|----------|----------|------------|----------|
| S3 (gzip) | 23.4 GB | 23.4 GB | 基准 |
| Sequence Store | 23.4 GB | ~6-8 GB | 65-70% |

### 访问性能
- **Reference Store**: 针对随机访问优化，支持快速区域查询
- **Sequence Store**: 针对批量读取优化，支持流式访问

## 🔗 相关链接

- [AWS Omics 官方文档](https://docs.aws.amazon.com/omics/)
- [Reference Store API参考](https://docs.aws.amazon.com/omics/latest/api/API_Operations_Amazon_Omics_Storage.html)
- [Sequence Store API参考](https://docs.aws.amazon.com/omics/latest/api/API_Operations_Amazon_Omics_Storage.html)

## 📋 下一步

设置完成后，您可以：
1. 在工作流中使用Reference Store和Sequence Store
2. 更新现有工作流参数
3. 享受更好的性能和更低的成本
4. 利用AWS Omics的原生集成功能

## 🤝 支持

如遇问题，请：
1. 查看本文档的故障排除部分
2. 检查AWS CloudWatch日志
3. 联系AWS支持团队
