# AWS Omics 存储设置脚本

本目录包含用于设置AWS Omics Reference Store和Sequence Store的脚本。

## 📁 文件说明

| 文件 | 功能 | 说明 |
|------|------|------|
| `00-setup_omics_environment.sh` | 一键完整设置 | 调用所有子脚本，完成完整环境设置 |
| `01-create_omics_stores.sh` | 创建存储 | 创建Reference Store和Sequence Store |
| `02-import_reference_genome.sh` | 导入参考基因组 | 导入ARS-UCD1.2奶牛参考基因组 |
| `03-import_sequencing_data.sh` | 导入测序数据 | 导入FASTQ测序数据到Sequence Store |
| `omics_stores_config.json` | 配置文件 | 存储创建的资源ID和配置信息 |

## 🚀 快速开始

### 一键设置（推荐）
```bash
./00-setup_omics_environment.sh
```

### 分步设置
```bash
# 1. 创建存储
./01-create_omics_stores.sh

# 2. 导入参考基因组
./02-import_reference_genome.sh

# 3. 导入测序数据
./03-import_sequencing_data.sh
```

## 📋 前提条件

1. **AWS CLI配置**: 确保已配置AWS CLI并有足够权限
2. **参考基因组文件**: 确保参考基因组已上传到S3
3. **测序数据**: 确保FASTQ文件可访问（本地S3或NCBI）

## 🎯 预期结果

设置完成后，您将获得：
- ✅ Reference Store (用于存储参考基因组)
- ✅ Sequence Store (用于存储测序数据)
- ✅ IAM服务角色 (OmicsServiceRole)
- ✅ 导入的参考基因组 (ARS-UCD1.2)
- ✅ 导入的测序数据 (Read Set)
- ✅ 配置文件 (omics_stores_config.json)

## 📊 优势

### 成本节省
- Sequence Store比S3存储节省65-70%成本
- 高压缩率减少存储需求

### 性能提升
- 针对基因组分析优化的访问模式
- 更快的数据读取和处理速度

### 集成优势
- 与AWS Omics工作流原生集成
- 自动元数据管理
- 版本控制支持

## 🔧 故障排除

常见问题请参考：`../../docs/user_guides/AWS_Omics_存储设置指南.md`

## 📞 支持

如需帮助，请查看项目文档或联系维护团队。
