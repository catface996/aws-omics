# 01-预处理 (Preprocessing)

## 目录说明
此目录包含基因组数据预处理相关的脚本和工具，按执行顺序编号。

## 主要功能
- 原始测序数据质量控制
- 接头去除和质量过滤
- 数据格式转换和标准化
- 数据下载和上传管理

## 脚本列表（按执行顺序）

### 01. 参考数据准备
- `01-download_genomic_data.sh` - 下载基因组参考数据（参考基因组、注释文件等）

### 02. 测序数据获取
- `02-download_data.sh` - 从 NCBI SRA 下载测序数据（FASTQ 文件）

### 03. 数据上传
- `03-upload_to_s3.sh` - 将 FASTQ 文件上传到 S3 存储

### 04. 工作流验证
- `04-validate_preprocessing_workflow.sh` - 验证预处理工作流配置和依赖

### 05. 预处理执行
- `05-run_preprocessing_workflow.sh` - 启动 AWS Omics 预处理工作流

## 执行流程

### 第一阶段：数据准备
1. **参考数据下载** (`01-download_genomic_data.sh`)
   - 下载奶牛参考基因组
   - 下载基因注释文件
   - 准备比对索引

2. **测序数据下载** (`02-download_data.sh`)
   - 从 NCBI SRA 下载 FASTQ 文件
   - 支持 SRA toolkit (fasterq-dump)
   - 自动分离双端测序数据

3. **数据上传** (`03-upload_to_s3.sh`)
   - 上传 FASTQ 文件到 S3
   - 上传参考基因组到 S3
   - 设置正确的 S3 权限

### 第二阶段：工作流执行
4. **配置验证** (`04-validate_preprocessing_workflow.sh`)
   - 检查 WDL 工作流语法
   - 验证输入参数配置
   - 确认 S3 文件路径有效性

5. **预处理运行** (`05-run_preprocessing_workflow.sh`)
   - 启动 AWS Omics 工作流
   - 执行质量控制 (FastQC)
   - 去接头和质量过滤 (Trimmomatic/fastp)
   - 去重复处理
   - 生成质量报告 (MultiQC)

## 输入数据源
- NCBI SRA 数据库的真实测序数据
- 支持双端测序 (Paired-end) 数据
- 奶牛参考基因组 (ARS-UCD1.2)

## 输出
预处理后的高质量测序数据，存储在 S3 中，准备用于 02-alignment 步骤。

## 使用说明
按顺序执行以下脚本：
```bash
# 1. 下载参考数据
./01-download_genomic_data.sh

# 2. 下载测序数据
./02-download_data.sh

# 3. 上传数据到 S3
./03-upload_to_s3.sh

# 4. 验证工作流配置
./04-validate_preprocessing_workflow.sh

# 5. 运行预处理工作流
./05-run_preprocessing_workflow.sh
```

## 配置文件
- `../workflows/wdl/preprocessing/inputs/preprocessing_inputs.json` - 预处理工作流输入参数
