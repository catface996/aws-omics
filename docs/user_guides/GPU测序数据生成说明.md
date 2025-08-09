# GPU加速测序数据生成器使用说明

## 🎯 功能概述

基于你已有的 `consolidated_variants.vcf`（3,154,016个变异）和奶牛参考基因组，生成高质量的模拟测序数据，充分利用Apple M4 Max GPU加速。

## 📊 输入数据状态

- ✅ **参考基因组**: `genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna` (2.75GB)
- ✅ **变异数据**: `consolidated_simulation_data/consolidated_variants.vcf` (230MB, 3,154,016个变异)
- ✅ **GPU设备**: Apple M4 Max (40核心) with MPS支持

## 🚀 快速开始

### 方法1: 使用自动化脚本（推荐）
```bash
# 运行自动化脚本，包含多种模式选择
./run_sequencing_simulation.sh
```

### 方法2: 手动运行
```bash
# 快速测试 (10x覆盖度, 前3个染色体)
python3 generate_sequencing_data_gpu.py \
    --reference genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    --variants consolidated_simulation_data/consolidated_variants.vcf \
    --output sequencing_data_test \
    --coverage 10 \
    --chromosomes NC_037328.1 NC_037329.1 NC_037330.1

# 标准模式 (30x覆盖度, 所有主要染色体)
python3 generate_sequencing_data_gpu.py \
    --reference genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    --variants consolidated_simulation_data/consolidated_variants.vcf \
    --output sequencing_data_30x \
    --coverage 30
```

## 🔧 参数说明

| 参数 | 说明 | 默认值 | 示例 |
|------|------|--------|------|
| `--reference` | 参考基因组FASTA文件 | 必需 | `genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna` |
| `--variants` | 变异VCF文件 | 必需 | `consolidated_simulation_data/consolidated_variants.vcf` |
| `--output` | 输出目录 | `./sequencing_data` | `my_sequencing_data` |
| `--coverage` | 测序覆盖深度 | `30.0` | `10`, `30`, `50` |
| `--chromosomes` | 指定染色体列表 | 所有主要染色体 | `NC_037328.1 NC_037329.1` |

## 📈 性能预期

基于你的Apple M4 Max GPU和之前的测试结果：

### 快速测试模式 (10x覆盖度, 3个染色体)
- **预计时间**: 2-3分钟
- **生成数据**: ~50万对reads
- **输出大小**: ~150MB FASTQ文件

### 标准模式 (30x覆盖度, 30个染色体)
- **预计时间**: 15-20分钟
- **生成数据**: ~1500万对reads
- **输出大小**: ~4.5GB FASTQ文件

### 高覆盖度模式 (50x覆盖度, 30个染色体)
- **预计时间**: 25-30分钟
- **生成数据**: ~2500万对reads
- **输出大小**: ~7.5GB FASTQ文件

## 📁 输出文件结构

```
sequencing_data/
├── NC_037328.1_R1.fastq.gz          # 染色体1正向reads
├── NC_037328.1_R2.fastq.gz          # 染色体1反向reads
├── NC_037329.1_R1.fastq.gz          # 染色体2正向reads
├── NC_037329.1_R2.fastq.gz          # 染色体2反向reads
├── ...
└── sequencing_generation_report.json # 详细生成报告
```

## 📊 生成的测序数据特性

### 技术参数
- **读长**: 150bp (配对端测序)
- **插入片段**: 300±50bp
- **测序错误率**: 0.1%
- **质量分数**: Phred+33编码
- **格式**: 标准FASTQ.gz压缩格式

### 生物学特性
- **变异密度**: 基于真实的3,154,016个变异
- **变异类型**: SNP和InDel
- **染色体覆盖**: 30个主要染色体 (NC_037328.1 - NC_037357.1)
- **基因组版本**: ARS-UCD1.2

## 🔍 质量验证

生成的数据包含：
1. **已知变异位点**: 所有变异都有准确的真值标注
2. **测序错误**: 模拟真实的测序错误模式
3. **覆盖度分布**: 符合泊松分布的覆盖度
4. **配对信息**: 正确的配对端关系

## 🧬 后续分析建议

生成的测序数据可用于：

### 1. 变异检测验证
```bash
# 使用GATK进行变异检测
gatk HaplotypeCaller \
    -R reference.fna \
    -I aligned_reads.bam \
    -O detected_variants.vcf
```

### 2. AWS Omics工作流测试
```bash
# 上传到S3并创建Read Set
aws s3 sync sequencing_data/ s3://your-bucket/sequencing-data/
aws omics create-read-set \
    --sequence-store-id <store-id> \
    --sources sourceFiles=s3://your-bucket/sequencing-data/
```

### 3. 准确性评估
```bash
# 比较检测到的变异与真值
bcftools isec \
    -p comparison_results \
    consolidated_variants.vcf \
    detected_variants.vcf
```

## ⚡ GPU优化特性

### 并行处理策略
- **染色体级并行**: 多个染色体同时处理
- **批量变异应用**: GPU批量处理变异
- **向量化操作**: 利用MPS进行向量化计算
- **内存优化**: 智能内存管理避免OOM

### 性能监控
脚本会实时显示：
- GPU利用率
- 处理速度 (reads/秒)
- 内存使用情况
- 预计完成时间

## 🐛 故障排除

### 常见问题

1. **GPU不可用**
   ```
   解决: 确保macOS支持MPS，更新到最新版本
   ```

2. **内存不足**
   ```
   解决: 减少覆盖度或分批处理染色体
   ```

3. **文件路径错误**
   ```
   解决: 检查参考基因组和VCF文件路径
   ```

### 性能调优

1. **提高并行度**
   ```python
   # 修改ThreadPoolExecutor的max_workers
   with ThreadPoolExecutor(max_workers=8) as executor:
   ```

2. **调整批处理大小**
   ```python
   # 修改batch_size参数
   batch_size = 2000  # 增加批处理大小
   ```

## 📋 检查清单

运行前确认：
- [ ] 参考基因组文件存在且完整
- [ ] VCF文件格式正确
- [ ] 有足够的磁盘空间 (至少10GB)
- [ ] GPU驱动正常工作
- [ ] Python环境包含所需依赖

## 🎯 下一步计划

1. **数据质量分析**: 使用FastQC分析生成的数据
2. **比对测试**: 使用BWA-MEM进行基因组比对
3. **变异检测**: 运行GATK HaplotypeCaller
4. **AWS Omics集成**: 上传数据并运行工作流
5. **准确性评估**: 计算变异检测的敏感性和特异性

---

**提示**: 首次运行建议使用快速测试模式，确认一切正常后再运行完整的标准模式。
