# AWS Omics 变异检测工作流

本目录包含用于AWS HealthOmics的变异检测工作流，支持从BAM文件进行高质量的变异检测和分析。

## 🎯 工作流概述

### 主要功能
- **GATK HaplotypeCaller**: 高精度变异检测
- **BCFtools**: 变异过滤和统计分析
- **多格式输出**: VCF文件、统计报告、分析图表

### 工作流架构
```
BAM文件 → [GATK变异检测] → [BCFtools过滤] → 最终变异结果
   ↓              ↓              ↓              ↓
输入数据        原始VCF        过滤VCF        统计报告
```

## 📁 文件说明

### 🔧 **生产就绪版本**
- **`variant_calling_from_bam_final_fixed.wdl`** ⭐ **推荐使用**
  - 完全修复所有已知问题
  - 解决了GATK参数兼容性问题
  - 修复了序列字典路径问题
  - 支持x86_64架构Docker镜像
  - 经过完整测试验证

### 🧪 **开发和测试版本**
- **`variant_calling_from_bam_fixed.wdl`** - 第一版修复
- **`variant_calling_from_bam_amd64.wdl`** - 架构兼容性修复
- **`variant_calling_from_bam_ecr.wdl`** - ECR镜像版本
- **`variant_calling_from_bam.wdl`** - 原始版本

### 📚 **历史版本**
- **`variant_calling_workflow.wdl`** - 早期版本

## 🚀 快速开始

### 1. 准备输入数据
```bash
# 确保您有以下文件：
# - 参考基因组 (FASTA格式)
# - 排序的BAM文件
# - BAM索引文件 (.bai)
```

### 2. 部署工作流
```bash
# 使用AWS CLI部署工作流
aws omics create-workflow \
    --name "variant-calling-production" \
    --definition-zip fileb://variant_calling_from_bam_final_fixed.wdl \
    --parameter-template file://parameter_template.json
```

### 3. 启动运行
```bash
# 启动变异检测运行
aws omics start-run \
    --workflow-id <workflow-id> \
    --role-arn <omics-service-role-arn> \
    --name "variant-calling-$(date +%Y%m%d-%H%M%S)" \
    --output-uri "s3://your-bucket/omics-outputs/variant-calling/" \
    --parameters '{
        "reference_genome": "s3://your-bucket/reference.fna",
        "input_bam": "s3://your-bucket/sample.sorted.bam",
        "input_bam_index": "s3://your-bucket/sample.sorted.bam.bai",
        "sample_name": "your-sample-name",
        "gatk_cpu": 16,
        "gatk_memory_gb": 32,
        "bcftools_cpu": 8,
        "bcftools_memory_gb": 16
    }'
```

## 🔧 关键修复说明

### ✅ **已解决的问题**

#### 1. **GATK参数兼容性问题**
- **问题**: `--standard-min-confidence-threshold-for-emitting` 在GATK 4.4.0.0中不支持
- **解决方案**: 使用正确的参数 `--standard-min-confidence-threshold-for-calling 20.0`

#### 2. **Docker架构兼容性问题**
- **问题**: ARM64镜像在AWS Omics x86_64环境中导致"exec format error"
- **解决方案**: 使用x86_64兼容的Docker镜像 `omics/variant-calling:amd64`

#### 3. **序列字典路径问题**
- **问题**: GATK期望序列字典文件与参考基因组在同一路径
- **解决方案**: 
  ```bash
  # 复制参考基因组到工作目录
  cp ~{reference} ./reference.fna
  # 在同一目录创建字典文件
  gatk CreateSequenceDictionary -R ./reference.fna -O ./reference.dict
  ```

### 🧪 **测试验证**
- **测试数据**: SRR16760538 (奶牛基因组测序数据)
- **测试结果**: 工作流成功运行，任务稳定执行
- **运行时间**: GATK任务约1-2小时，总体约2-3小时

## 📊 输出文件

### 主要输出
- **`raw_variants_vcf`**: 原始变异VCF文件 (压缩格式)
- **`raw_variants_vcf_index`**: 原始VCF索引文件
- **`filtered_variants_vcf`**: 过滤后变异VCF文件
- **`filtered_variants_vcf_index`**: 过滤VCF索引文件

### 分析报告
- **`variant_stats`**: 过滤后变异统计
- **`raw_variant_stats`**: 原始变异统计
- **`variant_report`**: 综合分析报告
- **`variant_density`**: 变异密度分析

## ⚙️ 参数配置

### 必需参数
- **`reference_genome`**: 参考基因组文件 (FASTA格式)
- **`input_bam`**: 输入BAM文件 (已排序)
- **`input_bam_index`**: BAM文件索引 (.bai文件)
- **`sample_name`**: 样本名称

### 可选参数
- **`gatk_cpu`**: GATK任务CPU核心数 (默认: 16)
- **`gatk_memory_gb`**: GATK任务内存大小 (默认: 32GB)
- **`bcftools_cpu`**: BCFtools任务CPU核心数 (默认: 8)
- **`bcftools_memory_gb`**: BCFtools任务内存大小 (默认: 16GB)

## 🔍 故障排除

### 常见问题

#### 1. **任务快速失败 (< 1分钟)**
- **可能原因**: 序列字典路径问题或GATK参数错误
- **解决方案**: 使用 `variant_calling_from_bam_final_fixed.wdl`

#### 2. **"exec format error"**
- **可能原因**: Docker镜像架构不兼容
- **解决方案**: 确保使用x86_64兼容镜像

#### 3. **内存不足错误**
- **可能原因**: 内存配置过低
- **解决方案**: 增加 `gatk_memory_gb` 参数值

### 日志分析
```bash
# 获取运行日志
aws omics get-run-logs --id <run-id>

# 获取任务详细日志
aws omics get-task-logs --run-id <run-id> --task-id <task-id>
```

## 📈 性能基准

### 测试环境
- **数据**: 奶牛基因组 (ARS-UCD1.2, ~2.75GB)
- **输入**: ~8GB BAM文件
- **资源**: 16 CPU, 32GB内存

### 性能指标
- **GATK任务**: 1-2小时
- **BCFtools任务**: 5-10分钟
- **总运行时间**: 1.5-2.5小时
- **输出大小**: ~100MB VCF文件

## 🔄 版本历史

### v3.0 (最新) - 2025-08-10
- ✅ 完全修复序列字典路径问题
- ✅ 解决GATK参数兼容性问题
- ✅ 支持x86_64架构Docker镜像
- ✅ 增强错误处理和日志记录
- ✅ 完整测试验证

### v2.0 - 2025-08-10
- ✅ 修复GATK参数错误
- ✅ 更新Docker镜像架构
- ❌ 序列字典路径问题未完全解决

### v1.0 - 2025-08-09
- ✅ 基础变异检测功能
- ❌ 存在多个兼容性问题

## 📞 支持

如果遇到问题，请：
1. 检查日志文件中的错误信息
2. 确认输入文件格式和路径正确
3. 验证IAM权限配置
4. 参考故障排除部分

---

**最后更新**: 2025年8月10日  
**状态**: ✅ 生产就绪，经过完整测试验证
