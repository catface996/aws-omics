# AWS Omics 奶牛基因组学演示项目

本项目展示如何使用AWS Omics服务构建完整的奶牛基因组分析流水线，从原始测序数据到生物学解释的端到端解决方案。

## 🔬 完整基因组分析流水线

### 高层次工作流架构
```
原始测序数据 → [数据预处理] → [序列比对] → [变异检测] → [变异注释] → [统计分析] → [生物学解释]
      ↓              ↓           ↓           ↓           ↓           ↓           ↓
   FASTQ文件      清洁数据      BAM文件      VCF文件    注释VCF    统计报告    功能解释
```

### 流水线各阶段详解

#### 1. **数据预处理 (Data Preprocessing)** 🧹
- **输入**: 原始FASTQ测序文件
- **工具**: FastQC + fastp + seqkit
- **功能**: 
  - 质量评估和控制
  - 接头去除和质量过滤
  - 序列去重复
- **输出**: 高质量清洁测序数据
- **AWS服务**: HealthOmics Workflows

#### 2. **序列比对 (Sequence Alignment)** 🎯
- **输入**: 清洁FASTQ文件 + 参考基因组
- **工具**: BWA-MEM2 / Minimap2
- **功能**:
  - 将测序reads比对到参考基因组
  - 生成比对结果和统计信息
- **输出**: 排序的BAM文件 + 比对统计
- **AWS服务**: HealthOmics Workflows + EC2

#### 3. **变异检测 (Variant Calling)** 🔍
- **输入**: BAM文件 + 参考基因组
- **工具**: GATK HaplotypeCaller / FreeBayes
- **功能**:
  - 识别SNPs (单核苷酸多态性)
  - 检测INDELs (插入缺失)
  - 结构变异检测
- **输出**: VCF变异文件
- **AWS服务**: HealthOmics Workflows + Batch

#### 4. **变异注释 (Variant Annotation)** 📝
- **输入**: VCF文件 + 注释数据库
- **工具**: VEP / SnpEff / ANNOVAR
- **功能**:
  - 基因功能注释
  - 蛋白质影响预测
  - 临床意义评估
- **输出**: 注释VCF文件
- **AWS服务**: HealthOmics Analytics + S3

#### 5. **统计分析 (Statistical Analysis)** 📊
- **输入**: 注释变异数据
- **工具**: R/Python + Bioconductor
- **功能**:
  - 变异频率统计
  - 群体遗传学分析
  - 关联性分析 (GWAS)
- **输出**: 统计报告和可视化
- **AWS服务**: SageMaker + QuickSight

#### 6. **生物学解释 (Biological Interpretation)** 🧬
- **输入**: 统计分析结果
- **工具**: GO分析 + KEGG通路分析
- **功能**:
  - 功能富集分析
  - 通路影响评估
  - 表型关联预测
- **输出**: 生物学意义报告
- **AWS服务**: Comprehend Medical + Lambda

### 🎯 奶牛基因组学应用场景

#### **育种改良**
- 产奶量相关基因变异识别
- 抗病性基因标记发现
- 肉质性状遗传标记

#### **疾病研究**
- 遗传性疾病变异检测
- 药物代谢基因分析
- 免疫相关基因研究

#### **群体遗传学**
- 品种间遗传差异分析
- 进化历史重构
- 遗传多样性评估

### 🏗️ AWS架构优势

#### **可扩展性**
- 自动资源调度和扩展
- 支持大规模并行处理
- 弹性计算资源管理

#### **成本优化**
- 按需付费模式
- 存储成本优化 (S3 + Omics Storage)
- 计算资源智能调度

#### **数据安全**
- 端到端加密
- 访问控制和审计
- 合规性保证

#### **集成生态**
- 与生物信息学工具深度集成
- 机器学习和AI能力
- 可视化和报告生成

## 🚀 快速开始

### 0. AWS Omics 存储设置 (新增 - 推荐)
```bash
# 一键设置Reference Store和Sequence Store
cd scripts/00-setup
./00-setup_omics_environment.sh
```

**优势**:
- ✅ 存储成本节省65-70%
- ✅ 针对基因组分析优化的性能
- ✅ 与AWS Omics工作流原生集成
- ✅ 自动元数据管理

### 1. 下载数据

https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR16760538&display=metadata

### 2. 生成变异基因组 (推荐)
```bash
# 激活GPU环境并运行
source gpu_env/bin/activate
./scripts/run_mutated_fna_generation.sh
```

### 3. 查看项目结构
```bash
# 查看详细的项目结构说明
cat PROJECT_STRUCTURE.md
```

## 📁 项目结构概览

```
aws-omics/
├── scripts/           # 所有脚本文件
│   ├── simulation/    # 基因组模拟脚本 ⭐
│   ├── analysis/      # 数据分析脚本
│   └── utils/         # 工具脚本
├── docs/              # 文档目录
│   ├── user_guides/   # 用户指南 ⭐
│   ├── technical/     # 技术文档
│   └── archive/       # 过时文档存档
├── workflows/         # WDL工作流文件
├── data/              # 输入数据
├── outputs/           # 输出结果
└── genomic_data/      # 原始基因组数据 (Git LFS)
```

## 🎯 核心功能

### 1. GPU加速变异基因组生成
- **脚本**: `scripts/simulation/generate_mutated_fna_gpu.py`
- **功能**: 基于3,154,016个变异生成完整的变异基因组
- **性能**: 充分利用Apple M4 Max GPU，~100,000变异/秒

### 2. 测序数据模拟
- **脚本**: `scripts/simulation/generate_sequencing_data_gpu.py`
- **功能**: 生成高质量的配对端测序数据
- **格式**: 标准FASTQ格式，支持不同覆盖深度

### 3. AWS Omics集成
- **工作流**: `workflows/wdl/`
- **功能**: 完整的变异检测和分析流水线

## 🧬 AWS Omics 数据预处理工作流

### 工作流架构
```
原始数据 → [InitialQC] → [RunFastp] → [RemoveDuplicates] → [FinalQC] → 最终清洁数据
                ↑              ↑              ↑              ↑
            质量评估      数据清洗(关键步骤)    去重复处理    最终质量评估
```

### 工作流步骤详解

#### 1. **InitialQC** - 初始质量评估
- **工具**: FastQC v0.12.1
- **功能**: 对原始测序数据进行质量评估
- **输出**: HTML质量报告 + ZIP详细数据
- **资源**: 8 CPU, 8GB 内存
- **运行时间**: ~15-25分钟

#### 2. **RunFastp** - 数据清洗 (关键步骤)
- **工具**: fastp v0.23.4
- **功能**: 
  - 🔧 接头序列自动检测和去除
  - 🎯 质量过滤 (Q20阈值)
  - 📏 长度过滤 (50-500bp)
  - 🧹 多聚体修剪 (poly-G/poly-X)
  - 🔄 滑动窗口质量修剪
  - 🧬 复杂度过滤 (30%阈值)
  - ✅ 序列错误校正
- **输出**: 清洁FASTQ + HTML/JSON报告
- **资源**: 16 CPU, 32GB 内存
- **运行时间**: ~10-15分钟

#### 3. **RemoveDuplicates** - 去重复处理
- **工具**: seqkit v2.5.1
- **功能**: 
  - 🔍 精确序列匹配去重 (exact method)
  - 📊 重复率统计分析
  - 💾 内存优化批处理
- **输出**: 去重FASTQ + 统计报告
- **资源**: 16 CPU, 24GB 内存
- **运行时间**: ~20-40分钟 (取决于数据量)

#### 4. **FinalQC** - 最终质量评估
- **工具**: FastQC v0.12.1
- **功能**: 对最终清洁数据进行质量验证
- **输出**: 最终质量报告
- **资源**: 8 CPU, 8GB 内存
- **运行时间**: ~10-20分钟

### 工作流特性

#### **并行执行优化**
- InitialQC 和 RunFastp 可并行启动
- 智能依赖管理确保数据流正确性
- 资源动态分配优化

#### **容错和监控**
- 详细的CloudWatch日志记录
- 任务级别的状态监控
- 自动错误诊断和报告

#### **数据质量保证**
- 多层质量控制检查
- 前后质量对比分析
- 详细的处理统计报告

### 使用示例

#### 启动预处理工作流
```bash
# 使用AWS Omics CLI启动工作流
aws omics start-run \
    --workflow-id 4976488 \
    --role-arn arn:aws:iam::ACCOUNT:role/OmicsServiceRole \
    --name "cow-preprocessing-$(date +%Y%m%d-%H%M%S)" \
    --output-uri "s3://your-bucket/omics-outputs/" \
    --parameters '{
        "sample_name": "SRR16760538",
        "input_fastq": "s3://your-sequence-store/path/to/data.fastq.gz",
        "fastp_cpu": 16,
        "fastp_memory_gb": 32,
        "min_quality": 20,
        "complexity_threshold": 30
    }'
```

#### 监控运行状态
```bash
# 检查运行状态
aws omics get-run --id <run-id>

# 查看任务详情
aws omics list-run-tasks --id <run-id>

# 获取运行日志
aws omics get-run-logs --id <run-id>
```

### 性能基准

**测试数据**: SRR16760538 (奶牛基因组测序数据, ~1.6GB压缩)

| 任务 | 运行时间 | CPU使用率 | 内存使用 | 输出大小 |
|------|----------|-----------|----------|----------|
| InitialQC | 22分钟 | ~60% | 1.6GB | 2MB报告 |
| RunFastp | 11分钟 | ~25% | 1.6GB | ~1.4GB清洁数据 |
| RemoveDuplicates | 32分钟 | ~0.1% | 1.6GB | ~1.3GB去重数据 |
| FinalQC | 15分钟 | ~60% | 1.6GB | 2MB报告 |

**总处理时间**: ~1.5小时  
**数据压缩率**: ~20% (去除低质量和重复序列)  
**质量提升**: Q20+ reads比例从85%提升到95%+

## 📊 数据说明

### 输入数据
- **参考基因组**: 奶牛 ARS-UCD1.2 版本 (~2.75GB)
- **变异数据**: 3,154,016个高质量变异 (230MB VCF)

### 输出数据
- **变异基因组**: 完整的FNA格式基因组文件
- **测序数据**: 标准FASTQ格式测序reads
- **分析报告**: 详细的处理统计信息

## 🔧 系统要求

- **操作系统**: macOS (支持MPS)
- **Python**: 3.8+ with PyTorch
- **GPU**: Apple Silicon (M1/M2/M4) 推荐
- **内存**: 16GB+ 推荐
- **存储**: 10GB+ 可用空间

## 📖 使用指南

### 快速测试
```bash
# 测试3个染色体的变异基因组生成
source gpu_env/bin/activate
python scripts/simulation/generate_mutated_fna_gpu.py \
    --reference genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna \
    --variants data/variants/consolidated_simulation_data/consolidated_variants.vcf \
    --output outputs/mutated_genome/test \
    --chromosomes NC_037328.1 NC_037329.1 NC_037330.1 \
    --max-workers 16
```

### 完整运行
```bash
# 处理所有30个主要染色体
./scripts/run_mutated_fna_generation.sh
```

## 📚 文档

- **用户指南**: `docs/user_guides/` - 操作说明和快速开始
- **技术文档**: `docs/technical/` - 深入的技术细节
- **项目结构**: `PROJECT_STRUCTURE.md` - 完整的项目结构说明

## ⚡ 性能特性

- **GPU加速**: 充分利用Apple M4 Max GPU
- **大规模并行**: 50,000变异/批次处理
- **内存优化**: 智能GPU内存管理
- **高吞吐量**: ~100,000变异/秒处理速度

## 🤝 贡献

欢迎提交Issue和Pull Request来改进这个项目！

## 📄 许可证

本项目遵循开源许可证。数据使用请遵循NCBI的使用条款。

---

**最后更新**: 2025年8月8日  
**项目状态**: ✅ 结构化整理完成，GPU优化就绪
