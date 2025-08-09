# 03_变异检测 (Variant Detection)

## 📋 目录说明

本目录包含基因组变异检测相关的工作流和工具，负责从比对结果中识别SNP、InDel等基因组变异。

## 🎯 主要功能

### 变异调用 (Variant Calling)
- **GATK HaplotypeCaller**: 高精度变异检测
- **BCFtools**: 快速变异调用
- **FreeBayes**: 贝叶斯变异检测
- **DeepVariant**: 深度学习变异检测

### 变异过滤 (Variant Filtering)
- **硬过滤**: 基于质量分数的严格过滤
- **VQSR**: 变异质量分数重校准
- **深度过滤**: 基于覆盖深度的过滤
- **等位基因频率过滤**: 基于群体频率的过滤

### 质量控制 (Quality Control)
- 变异调用统计
- Ti/Tv比率分析
- 变异密度分布
- 质量指标评估

## 📁 文件结构

```
03_variant_detection/
├── README.md                           # 本说明文档
├── variant_calling_workflow.wdl       # 变异检测主工作流
├── tasks/                              # 任务定义目录（待创建）
│   ├── gatk_haplotypecaller_task.wdl  # GATK变异调用任务
│   ├── bcftools_calling_task.wdl      # BCFtools变异调用任务
│   ├── variant_filtering_task.wdl     # 变异过滤任务
│   ├── variant_qc_task.wdl            # 变异质量控制任务
│   └── vcf_processing_task.wdl        # VCF文件处理任务
├── inputs/                             # 输入参数配置（待创建）
│   ├── variant_calling_inputs.json    # 变异检测参数配置
│   └── filtering_parameters.json      # 过滤参数配置
├── resources/                          # 资源文件（待创建）
│   ├── known_sites.vcf                # 已知变异位点
│   └── filtering_thresholds.txt       # 过滤阈值配置
└── scripts/                            # 辅助脚本（待创建）
    ├── variant_stats.py               # 变异统计分析
    └── vcf_validation.sh              # VCF文件验证
```

## 🚀 使用方法

### 1. 变异检测流程
```bash
# 使用AWS Omics运行变异检测工作流
aws omics start-run \
    --workflow-id <workflow-id> \
    --parameters file://inputs/variant_calling_inputs.json \
    --name "cattle-variant-calling-$(date +%Y%m%d)"
```

### 2. 输入数据要求
- **BAM文件**: 来自02_sequence_alignment步骤的比对结果
- **参考基因组**: 与比对使用的相同参考基因组
- **已知变异位点**: 用于基线重校准（可选）

### 3. 多样本变异检测
```bash
# 联合变异检测（多个样本）
aws omics start-run \
    --workflow-id <joint-calling-workflow-id> \
    --parameters file://inputs/joint_calling_inputs.json \
    --name "cattle-joint-calling-$(date +%Y%m%d)"
```

## ⚙️ 变异检测参数

### GATK HaplotypeCaller参数
- **最小基础质量**: 20
- **最小比对质量**: 20
- **标准置信度阈值**: 10.0
- **发射置信度阈值**: 30.0

### 变异过滤标准
- **质量分数**: QUAL > 30
- **覆盖深度**: DP > 10 && DP < 200
- **等位基因平衡**: AB > 0.2 && AB < 0.8
- **链偏向**: SOR < 3.0

### BCFtools参数
- **最小覆盖深度**: 10
- **最小变异质量**: 20
- **最小等位基因频率**: 0.1

## 📊 输出结果

### 主要输出文件
- **raw_variants.vcf**: 原始变异调用结果
- **filtered_variants.vcf**: 过滤后的高质量变异
- **variant_stats.txt**: 变异统计信息
- **variant_report.html**: 详细的变异分析报告

### 变异类型统计
- **SNP数量**: 单核苷酸多态性
- **InDel数量**: 插入缺失变异
- **Ti/Tv比率**: 转换/颠换比率
- **杂合/纯合比率**: 基因型分布

## 🔍 质量指标

### 变异调用质量
- **Ti/Tv比率**: 2.0-2.1 (全基因组)
- **InDel/SNP比率**: ~0.1-0.15
- **杂合率**: 0.1-0.3% (取决于物种)
- **新变异比例**: <5% (与已知数据库比较)

### 过滤效果
- **过滤前变异数**: 总变异数量
- **过滤后变异数**: 高质量变异数量
- **过滤率**: 通常20-40%
- **假阳性估计**: <5%

## 🧬 变异类型分析

### SNP分析
- 转换 (A↔G, C↔T)
- 颠换 (A/T↔C/G)
- 同义/非同义变异
- 启动子区域变异

### InDel分析
- 插入变异 (Insertions)
- 缺失变异 (Deletions)
- 框内/框移变异
- 重复序列区域InDel

### 结构变异 (可选)
- 大片段缺失 (>50bp)
- 重复序列扩增
- 倒位和易位
- 拷贝数变异 (CNV)

## 🔧 工作流特性

### 并行处理
- 染色体级别并行调用
- 区间分割并行处理
- 分布式计算优化

### 质量保证
- 多重验证机制
- 交叉验证不同算法
- 统计学质量评估

### 可扩展性
- 支持单样本和多样本
- 适配不同测序深度
- 灵活的过滤策略

## 📈 性能优化

### 计算资源
- **CPU**: 32-64核心推荐
- **内存**: 128GB推荐 (大基因组)
- **存储**: 高速SSD存储

### 优化策略
- 智能区间分割
- 内存映射文件访问
- 并行I/O操作

## 🔗 数据流向

### 输入来源
- **02_sequence_alignment**: BAM比对文件
- **genomic_data/02_reference_genome**: 参考基因组

### 输出去向
- **04_variant_annotation**: VCF文件用于功能注释
- **genomic_data/variants**: 存储变异检测结果

## 📚 相关文档

- [GATK最佳实践](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651)
- [BCFtools手册](http://samtools.github.io/bcftools/bcftools.html)
- [VCF格式规范](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
- [变异检测最佳实践](../docs/variant_calling_best_practices.md)

---

**注意**: 变异检测的准确性直接影响后续的功能分析，建议使用多种算法进行交叉验证。
