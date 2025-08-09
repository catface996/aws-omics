# 奶牛基因组数据预处理工作流

这是一个完整的测序数据预处理工作流，专为AWS Omics设计，包含了从原始测序数据到分析就绪数据的所有必要步骤。

## 🎯 功能特性

### 核心处理步骤
1. **质量评估 (FastQC)** - 原始数据和最终数据的质量控制
2. **接头去除 (Trimmomatic)** - 去除测序接头和低质量区域
3. **高级质量过滤 (FastP)** - 精细的质量和长度过滤
4. **去重复序列** - 使用多种算法去除重复reads
5. **综合报告 (MultiQC)** - 生成统一的质量控制报告

### 支持的数据类型
- ✅ 双端测序数据 (Paired-end)
- ✅ 单端测序数据 (Single-end)
- ✅ 批量样本处理
- ✅ 自定义参数配置

## 📁 文件结构

```
preprocessing/
├── preprocessing_workflow.wdl     # 主工作流文件
├── tasks/                         # 任务定义目录
│   ├── fastqc_task.wdl           # FastQC质量评估
│   ├── trimmomatic_task.wdl      # 接头去除和质量过滤
│   ├── fastp_task.wdl            # 高级质量过滤
│   ├── deduplication_task.wdl    # 去重复处理
│   └── multiqc_task.wdl          # 综合报告生成
├── inputs/                        # 输入参数示例
│   ├── preprocessing_inputs.json          # 双端数据配置
│   └── preprocessing_inputs_single_end.json # 单端数据配置
└── README.md                      # 本文档
```

## 🚀 快速开始

### 1. 准备输入数据

确保你的测序数据已上传到S3：
```bash
# 示例S3路径
s3://catface996-genomic/sequencing-data/
├── sample_001_R1.fastq.gz
├── sample_001_R2.fastq.gz
└── sample_002_SE.fastq.gz
```

### 2. 配置输入参数

编辑 `inputs/preprocessing_inputs.json`：
```json
{
  "PreprocessingWorkflow.sample_name": "your_sample_name",
  "PreprocessingWorkflow.fastq_r1": "s3://your-bucket/path/to/R1.fastq.gz",
  "PreprocessingWorkflow.fastq_r2": "s3://your-bucket/path/to/R2.fastq.gz",
  "PreprocessingWorkflow.min_length": 50,
  "PreprocessingWorkflow.min_quality": 20,
  "PreprocessingWorkflow.threads": 8
}
```

### 3. 运行工作流

使用提供的脚本：
```bash
./scripts/run_preprocessing_workflow.sh
```

或者直接使用AWS CLI：
```bash
aws omics start-run \
    --workflow-id <workflow-id> \
    --run-group-id <run-group-id> \
    --parameters file://inputs/preprocessing_inputs.json \
    --output-uri s3://catface996-genomic/omics-outputs/preprocessing/
```

## ⚙️ 参数配置

### 基本参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `sample_name` | String | 必需 | 样本名称 |
| `fastq_r1` | File | 必需 | 正向测序文件 |
| `fastq_r2` | File | 可选 | 反向测序文件 |
| `min_length` | Int | 50 | 最小读长阈值 |
| `min_quality` | Int | 20 | 最小质量分数 |
| `threads` | Int | 8 | 并行线程数 |

### 高级参数

#### FastP 过滤参数
```json
{
  "PreprocessingWorkflow.RunFastPPE.max_length": 500,
  "PreprocessingWorkflow.RunFastPPE.complexity_threshold": 30,
  "PreprocessingWorkflow.RunFastPPE.enable_polyg_trimming": true,
  "PreprocessingWorkflow.RunFastPPE.enable_polyx_trimming": true
}
```

#### 去重复方法选择
```json
{
  "PreprocessingWorkflow.RemoveDuplicatesPE.method": "fastuniq"
}
```
支持的方法：`fastuniq`, `clumpify`, `seqkit`

#### 计算资源配置
```json
{
  "PreprocessingWorkflow.RunFastQC.memory_gb": 8,
  "PreprocessingWorkflow.RunTrimmomaticPE.memory_gb": 16,
  "PreprocessingWorkflow.RunFastPPE.memory_gb": 16,
  "PreprocessingWorkflow.RemoveDuplicatesPE.memory_gb": 16
}
```

## 📊 输出结果

### 主要输出文件

1. **处理后的FASTQ文件**
   - `processed_r1` - 最终处理的R1文件
   - `processed_r2` - 最终处理的R2文件（双端数据）

2. **质量控制报告**
   - `initial_qc_r1_html` - 原始数据FastQC报告
   - `final_qc_r1_html` - 最终数据FastQC报告
   - `multiqc_html` - 综合质量报告

3. **处理统计**
   - `trimmomatic_log` - 接头去除日志
   - `fastp_json` - FastP处理统计
   - `dedup_log` - 去重复处理日志

### 输出目录结构

```
s3://catface996-genomic/omics-outputs/preprocessing/
└── <run-id>/
    ├── processed_data/
    │   ├── sample_R1_processed.fastq.gz
    │   └── sample_R2_processed.fastq.gz
    ├── qc_reports/
    │   ├── initial_fastqc_reports/
    │   ├── final_fastqc_reports/
    │   └── multiqc_report.html
    └── processing_logs/
        ├── trimmomatic.log
        ├── fastp.json
        └── deduplication.log
```

## 🔧 工作流管理

### 创建工作流
```bash
aws omics create-workflow \
    --name cow-genome-preprocessing \
    --description "奶牛基因组数据预处理工作流" \
    --engine WDL \
    --definition-zip fileb://preprocessing_workflow.zip
```

### 查看工作流状态
```bash
aws omics list-workflows --query "items[?name=='cow-genome-preprocessing']"
```

### 监控运行状态
```bash
aws omics get-run --id <run-id>
aws omics list-run-tasks --id <run-id>
```

## 📈 性能优化

### 计算资源建议

| 数据大小 | CPU | 内存 | 磁盘 | 预计时间 |
|----------|-----|------|------|----------|
| < 5GB | 8核 | 16GB | 100GB | 1-2小时 |
| 5-20GB | 16核 | 32GB | 200GB | 2-4小时 |
| > 20GB | 32核 | 64GB | 500GB | 4-8小时 |

### 成本优化建议

1. **使用Spot实例** - 在任务定义中设置 `preemptible: 2`
2. **合理设置资源** - 根据数据大小调整CPU和内存
3. **批量处理** - 多个样本一起处理可以降低固定成本
4. **存储优化** - 及时清理中间文件

## 🐛 故障排除

### 常见问题

1. **内存不足错误**
   - 增加任务的 `memory_gb` 参数
   - 减少并行线程数

2. **磁盘空间不足**
   - 增加 `disk_gb` 参数
   - 检查输入文件大小

3. **权限错误**
   - 确保IAM角色有S3访问权限
   - 检查Omics服务权限

4. **文件路径错误**
   - 验证S3路径是否正确
   - 确保文件存在且可访问

### 调试命令

```bash
# 查看运行详情
aws omics get-run --id <run-id> --output json

# 查看任务日志
aws omics get-run-task --id <run-id> --task-id <task-id>

# 下载日志文件
aws s3 cp s3://omics-logs/<run-id>/ ./logs/ --recursive
```

## 🤝 贡献指南

1. Fork 项目
2. 创建功能分支
3. 提交更改
4. 创建Pull Request

## 📄 许可证

本项目遵循开源许可证。

## 📞 支持

如有问题，请：
1. 查看故障排除部分
2. 检查AWS Omics文档
3. 提交Issue到项目仓库

---

**最后更新**: 2025年8月8日  
**版本**: 1.0  
**兼容性**: AWS Omics, WDL 1.0
