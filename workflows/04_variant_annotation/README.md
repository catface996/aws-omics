# 04_变异注释 (Variant Annotation)

## 📋 目录说明

本目录包含基因组变异功能注释相关的工作流和工具，负责为检测到的变异添加生物学功能信息。

## 🎯 主要功能

### 功能注释 (Functional Annotation)
- **SnpEff**: 基于基因组注释的变异效应预测
- **VEP**: Ensembl变异效应预测器
- **ANNOVAR**: 综合变异注释工具
- **ClinVar**: 临床相关变异注释

### 群体遗传学注释 (Population Genetics)
- **等位基因频率**: 群体中的变异频率
- **Hardy-Weinberg平衡**: 群体遗传学检验
- **连锁不平衡**: LD分析和单倍型构建
- **选择压力**: 正选择和负选择信号

### 功能影响预测 (Impact Prediction)
- **蛋白质功能预测**: SIFT, PolyPhen-2
- **保守性评分**: PhyloP, PhastCons
- **调控元件注释**: 启动子、增强子、转录因子结合位点
- **表观遗传修饰**: DNA甲基化、组蛋白修饰

## 📁 文件结构

```
04_variant_annotation/
├── README.md                           # 本说明文档
├── annotation_workflow.wdl            # 变异注释主工作流（待创建）
├── tasks/                              # 任务定义目录（待创建）
│   ├── snpeff_annotation_task.wdl     # SnpEff注释任务
│   ├── vep_annotation_task.wdl        # VEP注释任务
│   ├── annovar_task.wdl               # ANNOVAR注释任务
│   ├── population_stats_task.wdl      # 群体统计任务
│   └── functional_prediction_task.wdl # 功能预测任务
├── inputs/                             # 输入参数配置（待创建）
│   ├── annotation_inputs.json         # 注释参数配置
│   └── database_config.json           # 数据库配置文件
├── databases/                          # 注释数据库（待创建）
│   ├── snpeff_db/                     # SnpEff数据库
│   ├── vep_cache/                     # VEP缓存数据
│   ├── annovar_db/                    # ANNOVAR数据库
│   └── custom_annotations/            # 自定义注释数据
├── resources/                          # 资源文件（待创建）
│   ├── gene_ontology.obo              # Gene Ontology数据
│   ├── pathway_databases/             # 通路数据库
│   └── regulatory_elements.bed        # 调控元件注释
└── scripts/                            # 辅助脚本（待创建）
    ├── annotation_summary.py          # 注释结果汇总
    ├── functional_enrichment.R        # 功能富集分析
    └── variant_prioritization.py      # 变异优先级排序
```

## 🚀 使用方法

### 1. 变异注释流程
```bash
# 使用AWS Omics运行变异注释工作流
aws omics start-run \
    --workflow-id <workflow-id> \
    --parameters file://inputs/annotation_inputs.json \
    --name "cattle-annotation-$(date +%Y%m%d)"
```

### 2. 输入数据要求
- **VCF文件**: 来自03_variant_detection步骤的过滤变异
- **参考基因组**: 与变异检测使用的相同版本
- **基因组注释**: GTF/GFF格式的基因注释文件

### 3. 数据库准备
```bash
# 下载和准备SnpEff数据库
snpEff download -v ARS-UCD1.2.105

# 准备VEP缓存数据
vep_install -a cf -s bos_taurus -y ARS-UCD1.2

# 下载ANNOVAR数据库
annotate_variation.pl -buildver bosTau9 -downdb -webfrom annovar refGene bovinedb/
```

## ⚙️ 注释参数配置

### SnpEff参数
- **基因组版本**: ARS-UCD1.2.105
- **上游/下游区域**: 5000bp
- **剪接位点区域**: ±3bp
- **输出格式**: VCF + HTML报告

### VEP参数
- **缓存版本**: 最新版本
- **插件**: SIFT, PolyPhen, Condel
- **自定义注释**: 启用
- **输出格式**: VCF + 制表符分隔

### ANNOVAR参数
- **协议**: refGene,cytoBand,genomicSuperDups
- **操作**: g,r,r
- **构建版本**: bosTau9

## 📊 输出结果

### 主要输出文件
- **annotated_variants.vcf**: 完整注释的VCF文件
- **annotation_summary.html**: 注释结果可视化报告
- **functional_classes.txt**: 按功能分类的变异统计
- **high_impact_variants.vcf**: 高影响变异子集

### 注释信息类别
- **基因信息**: 基因名称、转录本ID、外显子/内含子
- **功能后果**: 同义/非同义、无义突变、剪接变异
- **蛋白质影响**: 氨基酸改变、结构域影响
- **调控影响**: 启动子、UTR、miRNA结合位点

## 🧬 功能分类

### 变异影响等级
- **高影响 (HIGH)**: 无义突变、框移变异、剪接位点
- **中等影响 (MODERATE)**: 非同义变异、密码子插入/缺失
- **低影响 (LOW)**: 同义变异、内含子变异
- **修饰符 (MODIFIER)**: 基因间区域、上游/下游变异

### 功能注释类别
- **编码序列变异**: 影响蛋白质编码的变异
- **调控区域变异**: 影响基因表达调控的变异
- **剪接变异**: 影响mRNA剪接的变异
- **非编码RNA变异**: 影响miRNA、lncRNA等的变异

## 📈 统计分析

### 变异分布统计
- **按染色体分布**: 各染色体变异密度
- **按基因分布**: 基因内变异数量统计
- **按功能分类**: 不同影响等级变异比例
- **按变异类型**: SNP vs InDel功能分布

### 功能富集分析
- **Gene Ontology富集**: 生物过程、分子功能、细胞组分
- **KEGG通路富集**: 代谢通路、信号通路
- **疾病关联**: 已知疾病相关基因富集
- **组织特异性**: 组织特异性表达基因富集

## 🔍 质量控制

### 注释完整性
- **注释覆盖率**: >95%变异获得功能注释
- **数据库版本**: 使用最新版本注释数据库
- **一致性检查**: 多个工具结果一致性验证
- **缺失注释**: 识别和报告未注释变异

### 功能预测准确性
- **已知变异验证**: 与已知功能变异比较
- **保守性一致性**: 保守区域变异影响预测
- **实验验证**: 与实验数据的一致性
- **文献支持**: 文献报道的功能变异验证

## 🔧 工作流特性

### 多工具集成
- 并行运行多个注释工具
- 结果整合和冲突解决
- 优先级排序算法

### 可定制化
- 支持自定义注释数据库
- 灵活的过滤和排序策略
- 用户定义的功能分类

### 可扩展性
- 支持大规模变异数据集
- 分布式并行处理
- 增量注释更新

## 🔗 数据流向

### 输入来源
- **03_variant_detection**: 过滤后的VCF文件
- **外部数据库**: 功能注释数据库

### 输出去向
- **05_statistical_analysis**: 注释变异用于统计分析
- **06_biological_interpretation**: 功能注释用于生物学解释

## 📚 相关文档

- [SnpEff用户手册](http://pcingola.github.io/SnpEff/)
- [VEP用户指南](https://www.ensembl.org/info/docs/tools/vep/index.html)
- [ANNOVAR文档](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
- [变异注释最佳实践](../docs/variant_annotation_best_practices.md)

---

**注意**: 变异注释的质量直接影响后续的生物学解释，建议使用多个工具进行交叉验证并定期更新注释数据库。
