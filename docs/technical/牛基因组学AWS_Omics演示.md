# 奶牛基因组分析 AWS Omics Demo方案

## 项目概述

本demo展示如何使用AWS Omics服务构建完整的奶牛基因组分析流水线，从原始测序数据到生物学解释的端到端解决方案。

### 目标
- 演示AWS Omics在大规模基因组分析中的应用
- 建立标准化的奶牛基因组分析流程
- 实现从基因型到表型的关联分析
- 提供可扩展的云端基因组学解决方案

### 技术架构
```
数据上传 → 预处理 → 参考基因组 → 序列比对 → 变异检测 → 变异注释 → 统计分析 → 生物学解释
    ↓         ↓         ↓          ↓         ↓         ↓         ↓         ↓
  S3存储   Omics     Reference   Omics    Omics    Omics    Omics    可视化
          Workflows   Store     Workflows Workflows Analytics Analytics  报告
```

## 数据准备

### 样本信息
- **物种**: 奶牛 (Bos taurus)
- **样本数量**: 100个个体
- **测序类型**: 全基因组重测序 (WGS)
- **测序平台**: Illumina NovaSeq 6000
- **测序深度**: 30X
- **读长**: 150bp paired-end

### 数据规模估算
- 每个样本原始数据: ~90GB (FASTQ)
- 总数据量: ~9TB
- 比对后数据: ~3TB (BAM文件)
- 变异数据: ~500MB (VCF文件)

## 阶段1: 数据上传与存储

### 1.1 S3存储架构设计

```
s3://cattle-genomics-demo/
├── raw-data/
│   ├── sample001/
│   │   ├── sample001_R1.fastq.gz
│   │   └── sample001_R2.fastq.gz
│   └── sample002/
├── reference/
│   ├── Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
│   ├── Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.fai
│   └── annotations/
├── processed/
│   ├── aligned/
│   ├── variants/
│   └── annotations/
└── results/
    ├── statistics/
    └── reports/
```

### 1.2 数据上传策略

**使用AWS CLI批量上传:**
```bash
# 配置多部分上传和并行传输
aws configure set default.s3.max_concurrent_requests 20
aws configure set default.s3.max_bandwidth 1GB/s
aws configure set default.s3.multipart_threshold 64MB

# 批量上传FASTQ文件
aws s3 sync ./local-fastq-files/ s3://cattle-genomics-demo/raw-data/ \
  --storage-class STANDARD_IA \
  --metadata project=cattle-genomics,stage=raw-data
```

**使用AWS DataSync (大规模数据):**
```bash
# 创建DataSync任务用于大规模数据传输
aws datasync create-task \
  --source-location-arn arn:aws:datasync:region:account:location/loc-xxx \
  --destination-location-arn arn:aws:datasync:region:account:location/loc-yyy \
  --cloud-watch-log-group-arn arn:aws:logs:region:account:log-group:datasync
```

### 1.3 存储优化配置

**生命周期管理:**
```json
{
  "Rules": [
    {
      "ID": "CattleGenomicsLifecycle",
      "Status": "Enabled",
      "Transitions": [
        {
          "Days": 30,
          "StorageClass": "STANDARD_IA"
        },
        {
          "Days": 90,
          "StorageClass": "GLACIER"
        },
        {
          "Days": 365,
          "StorageClass": "DEEP_ARCHIVE"
        }
      ]
    }
  ]
}
```

## 阶段2: 数据预处理

### 2.1 质量控制工作流

**创建Omics工作流:**
```bash
# 创建预处理工作流
aws omics create-workflow \
  --name "cattle-qc-preprocessing" \
  --description "Quality control and preprocessing for cattle genomics data" \
  --definition-uri s3://cattle-genomics-demo/workflows/qc-preprocessing.wdl \
  --parameter-template file://qc-parameters.json
```

**WDL工作流定义 (qc-preprocessing.wdl):**
```wdl
version 1.0

workflow CattleQCPreprocessing {
  input {
    File fastq_r1
    File fastq_r2
    String sample_id
    Int threads = 8
  }

  call FastQC {
    input:
      fastq_r1 = fastq_r1,
      fastq_r2 = fastq_r2,
      sample_id = sample_id
  }

  call Trimmomatic {
    input:
      fastq_r1 = fastq_r1,
      fastq_r2 = fastq_r2,
      sample_id = sample_id,
      threads = threads
  }

  call FastQCPost {
    input:
      fastq_r1 = Trimmomatic.trimmed_r1,
      fastq_r2 = Trimmomatic.trimmed_r2,
      sample_id = sample_id
  }

  output {
    File trimmed_r1 = Trimmomatic.trimmed_r1
    File trimmed_r2 = Trimmomatic.trimmed_r2
    File qc_report_pre = FastQC.qc_report
    File qc_report_post = FastQCPost.qc_report
  }
}

task FastQC {
  input {
    File fastq_r1
    File fastq_r2
    String sample_id
  }

  command <<<
    fastqc ~{fastq_r1} ~{fastq_r2} -o ./
    multiqc . -n ~{sample_id}_pre_qc_report
  >>>

  output {
    File qc_report = "${sample_id}_pre_qc_report.html"
  }

  runtime {
    docker: "quay.io/biocontainers/fastqc:0.11.9--0"
    cpu: 2
    memory: "4 GB"
  }
}

task Trimmomatic {
  input {
    File fastq_r1
    File fastq_r2
    String sample_id
    Int threads
  }

  command <<<
    trimmomatic PE -threads ~{threads} \
      ~{fastq_r1} ~{fastq_r2} \
      ~{sample_id}_trimmed_R1.fastq.gz ~{sample_id}_unpaired_R1.fastq.gz \
      ~{sample_id}_trimmed_R2.fastq.gz ~{sample_id}_unpaired_R2.fastq.gz \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  >>>

  output {
    File trimmed_r1 = "${sample_id}_trimmed_R1.fastq.gz"
    File trimmed_r2 = "${sample_id}_trimmed_R2.fastq.gz"
  }

  runtime {
    docker: "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
    cpu: threads
    memory: "8 GB"
  }
}
```

### 2.2 批量处理配置

**参数模板 (qc-parameters.json):**
```json
{
  "CattleQCPreprocessing.fastq_r1": {
    "description": "Forward reads FASTQ file"
  },
  "CattleQCPreprocessing.fastq_r2": {
    "description": "Reverse reads FASTQ file"
  },
  "CattleQCPreprocessing.sample_id": {
    "description": "Sample identifier"
  },
  "CattleQCPreprocessing.threads": {
    "description": "Number of CPU threads",
    "optional": true
  }
}
```

**批量运行脚本:**
```bash
#!/bin/bash
# 批量提交预处理任务

WORKFLOW_ID="12345678"  # 从create-workflow命令获取
ROLE_ARN="arn:aws:iam::account:role/OmicsWorkflowRole"

# 读取样本列表
while IFS=',' read -r sample_id fastq_r1 fastq_r2; do
  echo "Processing sample: $sample_id"
  
  # 创建参数文件
  cat > ${sample_id}_params.json << EOF
{
  "CattleQCPreprocessing.fastq_r1": "$fastq_r1",
  "CattleQCPreprocessing.fastq_r2": "$fastq_r2",
  "CattleQCPreprocessing.sample_id": "$sample_id",
  "CattleQCPreprocessing.threads": 8
}
EOF

  # 提交工作流运行
  aws omics start-run \
    --workflow-id $WORKFLOW_ID \
    --workflow-type PRIVATE \
    --role-arn $ROLE_ARN \
    --name "qc-${sample_id}" \
    --parameters file://${sample_id}_params.json \
    --output-uri s3://cattle-genomics-demo/processed/qc/${sample_id}/
    
done < sample_list.csv
```

## 阶段3: 参考基因组准备

### 3.1 Reference Store配置

**创建参考基因组存储:**
```bash
# 创建参考基因组存储
aws omics create-reference-store \
  --name "cattle-reference-store" \
  --description "Reference genomes for cattle genomics analysis"
```

**导入奶牛参考基因组:**
```bash
# 导入ARS-UCD1.2参考基因组
aws omics start-reference-import-job \
  --reference-store-id "1234567890123456" \
  --role-arn "arn:aws:iam::account:role/OmicsReferenceRole" \
  --sources '[
    {
      "sourceFile": "s3://cattle-genomics-demo/reference/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
      "name": "Bos_taurus_ARS-UCD1.2",
      "description": "Cattle reference genome ARS-UCD1.2"
    }
  ]'
```

### 3.2 参考基因组索引

**BWA索引构建工作流:**
```wdl
version 1.0

workflow BuildReferenceIndex {
  input {
    File reference_fasta
    String reference_name
  }

  call BWAIndex {
    input:
      reference_fasta = reference_fasta,
      reference_name = reference_name
  }

  call SamtoolsIndex {
    input:
      reference_fasta = reference_fasta
  }

  call CreateDict {
    input:
      reference_fasta = reference_fasta,
      reference_name = reference_name
  }

  output {
    Array[File] bwa_index = BWAIndex.index_files
    File fasta_index = SamtoolsIndex.fai_file
    File sequence_dict = CreateDict.dict_file
  }
}

task BWAIndex {
  input {
    File reference_fasta
    String reference_name
  }

  command <<<
    bwa index ~{reference_fasta}
    ls ~{reference_fasta}*
  >>>

  output {
    Array[File] index_files = glob("${reference_fasta}*")
  }

  runtime {
    docker: "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 SSD"
  }
}
```

## 阶段4: 序列比对

### 4.1 比对工作流设计

**BWA-MEM比对工作流:**
```wdl
version 1.0

workflow CattleAlignment {
  input {
    File fastq_r1
    File fastq_r2
    File reference_fasta
    Array[File] reference_index
    String sample_id
    String read_group_info
  }

  call BWAMem {
    input:
      fastq_r1 = fastq_r1,
      fastq_r2 = fastq_r2,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      sample_id = sample_id,
      read_group_info = read_group_info
  }

  call SortSam {
    input:
      input_sam = BWAMem.aligned_sam,
      sample_id = sample_id
  }

  call MarkDuplicates {
    input:
      input_bam = SortSam.sorted_bam,
      sample_id = sample_id
  }

  call BaseRecalibrator {
    input:
      input_bam = MarkDuplicates.dedup_bam,
      reference_fasta = reference_fasta,
      sample_id = sample_id
  }

  output {
    File final_bam = BaseRecalibrator.recalibrated_bam
    File final_bai = BaseRecalibrator.recalibrated_bai
    File alignment_metrics = MarkDuplicates.metrics_file
  }
}

task BWAMem {
  input {
    File fastq_r1
    File fastq_r2
    File reference_fasta
    Array[File] reference_index
    String sample_id
    String read_group_info
  }

  command <<<
    set -e
    bwa mem -t 16 -R "~{read_group_info}" \
      ~{reference_fasta} ~{fastq_r1} ~{fastq_r2} > ~{sample_id}.sam
  >>>

  output {
    File aligned_sam = "${sample_id}.sam"
  }

  runtime {
    docker: "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    cpu: 16
    memory: "32 GB"
    disks: "local-disk 200 SSD"
  }
}

task SortSam {
  input {
    File input_sam
    String sample_id
  }

  command <<<
    samtools sort -@ 8 -o ~{sample_id}_sorted.bam ~{input_sam}
    samtools index ~{sample_id}_sorted.bam
  >>>

  output {
    File sorted_bam = "${sample_id}_sorted.bam"
    File sorted_bai = "${sample_id}_sorted.bam.bai"
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.15.1--h1170115_0"
    cpu: 8
    memory: "16 GB"
  }
}
```

### 4.2 质量控制指标

**比对质量评估:**
```wdl
task AlignmentQC {
  input {
    File input_bam
    File reference_fasta
    String sample_id
  }

  command <<<
    # 计算比对统计
    samtools flagstat ~{input_bam} > ~{sample_id}_flagstat.txt
    samtools stats ~{input_bam} > ~{sample_id}_stats.txt
    
    # 计算覆盖度
    samtools depth ~{input_bam} | awk '{sum+=$3; count++} END {print "Average depth:", sum/count}' > ~{sample_id}_coverage.txt
    
    # Picard CollectAlignmentSummaryMetrics
    java -jar /picard.jar CollectAlignmentSummaryMetrics \
      I=~{input_bam} \
      O=~{sample_id}_alignment_summary.txt \
      R=~{reference_fasta}
  >>>

  output {
    File flagstat = "${sample_id}_flagstat.txt"
    File stats = "${sample_id}_stats.txt"
    File coverage = "${sample_id}_coverage.txt"
    File alignment_summary = "${sample_id}_alignment_summary.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/picard:2.27.4--hdfd78af_0"
    cpu: 4
    memory: "8 GB"
  }
}
```

## 阶段5: 变异检测

### 5.1 GATK变异检测工作流

**HaplotypeCaller变异检测:**
```wdl
version 1.0

workflow CattleVariantCalling {
  input {
    Array[File] input_bams
    Array[File] input_bais
    File reference_fasta
    File reference_fai
    File reference_dict
    Array[String] sample_ids
    File known_sites_vcf  # 已知变异位点
  }

  # 单样本变异检测
  scatter (i in range(length(input_bams))) {
    call HaplotypeCaller {
      input:
        input_bam = input_bams[i],
        input_bai = input_bais[i],
        reference_fasta = reference_fasta,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        sample_id = sample_ids[i]
    }
  }

  # 合并GVCF文件
  call CombineGVCFs {
    input:
      input_gvcfs = HaplotypeCaller.output_gvcf,
      reference_fasta = reference_fasta,
      reference_fai = reference_fai,
      reference_dict = reference_dict
  }

  # 联合基因型检测
  call GenotypeGVCFs {
    input:
      combined_gvcf = CombineGVCFs.combined_gvcf,
      reference_fasta = reference_fasta,
      reference_fai = reference_fai,
      reference_dict = reference_dict
  }

  # 变异质量过滤
  call VariantFiltration {
    input:
      input_vcf = GenotypeGVCFs.output_vcf,
      reference_fasta = reference_fasta
  }

  output {
    File final_vcf = VariantFiltration.filtered_vcf
    File final_vcf_index = VariantFiltration.filtered_vcf_index
    Array[File] individual_gvcfs = HaplotypeCaller.output_gvcf
  }
}

task HaplotypeCaller {
  input {
    File input_bam
    File input_bai
    File reference_fasta
    File reference_fai
    File reference_dict
    String sample_id
  }

  command <<<
    gatk HaplotypeCaller \
      -R ~{reference_fasta} \
      -I ~{input_bam} \
      -O ~{sample_id}.g.vcf.gz \
      -ERC GVCF \
      --native-pair-hmm-threads 4 \
      --max-alternate-alleles 3
  >>>

  output {
    File output_gvcf = "${sample_id}.g.vcf.gz"
    File output_gvcf_index = "${sample_id}.g.vcf.gz.tbi"
  }

  runtime {
    docker: "broadinstitute/gatk:4.3.0.0"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 SSD"
  }
}

task GenotypeGVCFs {
  input {
    File combined_gvcf
    File reference_fasta
    File reference_fai
    File reference_dict
  }

  command <<<
    gatk GenotypeGVCFs \
      -R ~{reference_fasta} \
      -V ~{combined_gvcf} \
      -O cattle_cohort.vcf.gz \
      --max-alternate-alleles 3
  >>>

  output {
    File output_vcf = "cattle_cohort.vcf.gz"
    File output_vcf_index = "cattle_cohort.vcf.gz.tbi"
  }

  runtime {
    docker: "broadinstitute/gatk:4.3.0.0"
    cpu: 8
    memory: "32 GB"
    disks: "local-disk 200 SSD"
  }
}

task VariantFiltration {
  input {
    File input_vcf
    File reference_fasta
  }

  command <<<
    # SNP过滤
    gatk SelectVariants \
      -R ~{reference_fasta} \
      -V ~{input_vcf} \
      -select-type SNP \
      -O snps.vcf.gz

    gatk VariantFiltration \
      -R ~{reference_fasta} \
      -V snps.vcf.gz \
      -O snps_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "SOR > 3.0" --filter-name "SOR3" \
      --filter-expression "FS > 60.0" --filter-name "FS60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ40"

    # INDEL过滤
    gatk SelectVariants \
      -R ~{reference_fasta} \
      -V ~{input_vcf} \
      -select-type INDEL \
      -O indels.vcf.gz

    gatk VariantFiltration \
      -R ~{reference_fasta} \
      -V indels.vcf.gz \
      -O indels_filtered.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "FS > 200.0" --filter-name "FS200" \
      --filter-expression "SOR > 10.0" --filter-name "SOR10"

    # 合并过滤后的变异
    gatk MergeVcfs \
      -I snps_filtered.vcf.gz \
      -I indels_filtered.vcf.gz \
      -O cattle_filtered.vcf.gz
  >>>

  output {
    File filtered_vcf = "cattle_filtered.vcf.gz"
    File filtered_vcf_index = "cattle_filtered.vcf.gz.tbi"
  }

  runtime {
    docker: "broadinstitute/gatk:4.3.0.0"
    cpu: 4
    memory: "16 GB"
  }
}
```

### 5.2 结构变异检测

**使用Manta检测结构变异:**
```wdl
task MantaSV {
  input {
    Array[File] input_bams
    Array[File] input_bais
    File reference_fasta
    File reference_fai
  }

  command <<<
    # 配置Manta
    configManta.py \
      --bam ~{sep=' --bam ' input_bams} \
      --referenceFasta ~{reference_fasta} \
      --runDir ./manta_analysis

    # 运行Manta
    cd manta_analysis
    python runWorkflow.py -m local -j 8
  >>>

  output {
    File sv_vcf = "manta_analysis/results/variants/diploidSV.vcf.gz"
    File sv_vcf_index = "manta_analysis/results/variants/diploidSV.vcf.gz.tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/manta:1.6.0--py27_0"
    cpu: 8
    memory: "16 GB"
  }
}
```

## 阶段6: 变异注释

### 6.1 功能注释工作流

**使用VEP进行变异注释:**
```wdl
version 1.0

workflow VariantAnnotation {
  input {
    File input_vcf
    File vep_cache_tar
    File reference_fasta
    String species = "bos_taurus"
    String assembly = "ARS-UCD1.2"
  }

  call VEPAnnotation {
    input:
      input_vcf = input_vcf,
      vep_cache_tar = vep_cache_tar,
      reference_fasta = reference_fasta,
      species = species,
      assembly = assembly
  }

  call SnpEffAnnotation {
    input:
      input_vcf = input_vcf,
      reference_fasta = reference_fasta
  }

  call CombineAnnotations {
    input:
      vep_vcf = VEPAnnotation.annotated_vcf,
      snpeff_vcf = SnpEffAnnotation.annotated_vcf
  }

  output {
    File final_annotated_vcf = CombineAnnotations.combined_vcf
    File annotation_summary = VEPAnnotation.summary_html
  }
}

task VEPAnnotation {
  input {
    File input_vcf
    File vep_cache_tar
    File reference_fasta
    String species
    String assembly
  }

  command <<<
    # 解压VEP缓存
    tar -xzf ~{vep_cache_tar}
    
    # 运行VEP注释
    vep \
      --input_file ~{input_vcf} \
      --output_file cattle_vep_annotated.vcf \
      --format vcf \
      --vcf \
      --species ~{species} \
      --assembly ~{assembly} \
      --cache \
      --dir_cache ./vep_cache \
      --fasta ~{reference_fasta} \
      --gene_phenotype \
      --regulatory \
      --protein \
      --symbol \
      --numbers \
      --biotype \
      --total_length \
      --canonical \
      --ccds \
      --domains \
      --pubmed \
      --variant_class \
      --stats_file cattle_vep_summary.html \
      --warning_file cattle_vep_warnings.txt \
      --fork 4
  >>>

  output {
    File annotated_vcf = "cattle_vep_annotated.vcf"
    File summary_html = "cattle_vep_summary.html"
    File warnings = "cattle_vep_warnings.txt"
  }

  runtime {
    docker: "ensemblorg/ensembl-vep:release_107.0"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 100 SSD"
  }
}

task SnpEffAnnotation {
  input {
    File input_vcf
    File reference_fasta
  }

  command <<<
    # 使用SnpEff进行注释
    snpEff -Xmx8g \
      -v Bos_taurus \
      ~{input_vcf} \
      > cattle_snpeff_annotated.vcf \
      2> cattle_snpeff.log
      
    # 生成统计报告
    snpEff -Xmx8g \
      -stats cattle_snpeff_stats.html \
      Bos_taurus \
      ~{input_vcf} \
      > /dev/null
  >>>

  output {
    File annotated_vcf = "cattle_snpeff_annotated.vcf"
    File stats_html = "cattle_snpeff_stats.html"
    File log_file = "cattle_snpeff.log"
  }

  runtime {
    docker: "quay.io/biocontainers/snpeff:5.1--hdfd78af_2"
    cpu: 4
    memory: "8 GB"
  }
}
```

### 6.2 数据库注释

**添加外部数据库注释:**
```wdl
task DatabaseAnnotation {
  input {
    File input_vcf
    File dbsnp_vcf
    File clinvar_vcf
    File cosmic_vcf
  }

  command <<<
    # 添加dbSNP注释
    bcftools annotate \
      -a ~{dbsnp_vcf} \
      -c ID \
      ~{input_vcf} \
      -O z \
      -o cattle_dbsnp.vcf.gz

    # 添加ClinVar注释
    bcftools annotate \
      -a ~{clinvar_vcf} \
      -c INFO/CLNSIG,INFO/CLNDN \
      cattle_dbsnp.vcf.gz \
      -O z \
      -o cattle_clinvar.vcf.gz

    # 添加COSMIC注释
    bcftools annotate \
      -a ~{cosmic_vcf} \
      -c INFO/CNT \
      cattle_clinvar.vcf.gz \
      -O z \
      -o cattle_annotated_final.vcf.gz

    # 创建索引
    tabix -p vcf cattle_annotated_final.vcf.gz
  >>>

  output {
    File final_vcf = "cattle_annotated_final.vcf.gz"
    File final_vcf_index = "cattle_annotated_final.vcf.gz.tbi"
  }

  runtime {
    docker: "quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0"
    cpu: 2
    memory: "4 GB"
  }
}
```

## 阶段7: 统计分析

### 7.1 群体遗传学分析

**使用PLINK进行群体分析:**
```wdl
task PopulationAnalysis {
  input {
    File input_vcf
    File sample_info
    String output_prefix = "cattle_population"
  }

  command <<<
    # 转换VCF到PLINK格式
    plink2 \
      --vcf ~{input_vcf} \
      --make-bed \
      --out ~{output_prefix}

    # 计算等位基因频率
    plink2 \
      --bfile ~{output_prefix} \
      --freq \
      --out ~{output_prefix}_freq

    # Hardy-Weinberg平衡检验
    plink2 \
      --bfile ~{output_prefix} \
      --hardy \
      --out ~{output_prefix}_hwe

    # 连锁不平衡分析
    plink2 \
      --bfile ~{output_prefix} \
      --r2 \
      --ld-window-kb 1000 \
      --ld-window 99999 \
      --ld-window-r2 0.2 \
      --out ~{output_prefix}_ld

    # 主成分分析
    plink2 \
      --bfile ~{output_prefix} \
      --pca 10 \
      --out ~{output_prefix}_pca

    # 近交系数计算
    plink2 \
      --bfile ~{output_prefix} \
      --het \
      --out ~{output_prefix}_het
  >>>

  output {
    File allele_freq = "${output_prefix}_freq.afreq"
    File hwe_results = "${output_prefix}_hwe.hardy"
    File ld_results = "${output_prefix}_ld.ld"
    File pca_eigenval = "${output_prefix}_pca.eigenval"
    File pca_eigenvec = "${output_prefix}_pca.eigenvec"
    File het_results = "${output_prefix}_het.het"
  }

  runtime {
    docker: "quay.io/biocontainers/plink2:2.00a3.3--h712d239_0"
    cpu: 4
    memory: "8 GB"
  }
}
```

### 7.2 关联分析

**GWAS分析工作流:**
```wdl
task GWASAnalysis {
  input {
    File genotype_bed
    File genotype_bim
    File genotype_fam
    File phenotype_file
    String trait_name
    Array[String] covariates
  }

  command <<<
    # 质量控制过滤
    plink2 \
      --bfile ~{basename(genotype_bed, ".bed")} \
      --maf 0.05 \
      --geno 0.1 \
      --mind 0.1 \
      --hwe 1e-6 \
      --make-bed \
      --out qc_filtered

    # GWAS分析
    plink2 \
      --bfile qc_filtered \
      --pheno ~{phenotype_file} \
      --pheno-name ~{trait_name} \
      --covar ~{phenotype_file} \
      --covar-name ~{sep=',' covariates} \
      --glm \
      --adjust \
      --out gwas_~{trait_name}

    # 曼哈顿图数据准备
    awk 'NR>1 {print $1, $3, $12}' gwas_~{trait_name}.~{trait_name}.glm.linear > manhattan_data.txt
  >>>

  output {
    File gwas_results = "gwas_${trait_name}.${trait_name}.glm.linear"
    File gwas_adjusted = "gwas_${trait_name}.${trait_name}.glm.linear.adjusted"
    File manhattan_data = "manhattan_data.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/plink2:2.00a3.3--h712d239_0"
    cpu: 8
    memory: "16 GB"
  }
}
```

## 阶段8: 生物学解释

### 8.1 功能富集分析

**基因本体论(GO)富集分析:**
```python
# Python脚本用于功能富集分析
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import requests
import json

def go_enrichment_analysis(gene_list, background_genes, go_annotations):
    """
    执行GO富集分析
    """
    enrichment_results = []
    
    for go_term, go_genes in go_annotations.items():
        # 计算重叠基因
        overlap = set(gene_list) & set(go_genes) & set(background_genes)
        
        # 超几何检验
        M = len(background_genes)  # 总基因数
        n = len(set(go_genes) & set(background_genes))  # GO term中的基因数
        N = len(gene_list)  # 候选基因数
        k = len(overlap)  # 重叠基因数
        
        if k > 0:
            p_value = hypergeom.sf(k-1, M, n, N)
            
            enrichment_results.append({
                'GO_term': go_term,
                'p_value': p_value,
                'overlap_genes': len(overlap),
                'go_genes': n,
                'candidate_genes': N,
                'genes': list(overlap)
            })
    
    # 多重检验校正
    df = pd.DataFrame(enrichment_results)
    df['fdr'] = df['p_value'] * len(df) / (df['p_value'].rank())
    
    return df.sort_values('p_value')

# 使用示例
significant_genes = ['GENE1', 'GENE2', 'GENE3']  # GWAS显著基因
go_results = go_enrichment_analysis(significant_genes, all_genes, go_annotations)
```

**KEGG通路分析:**
```wdl
task KEGGPathwayAnalysis {
  input {
    File gene_list
    String species_code = "bta"  # 奶牛KEGG代码
  }

  command <<<
    # 使用R进行KEGG分析
    Rscript << 'EOF'
    library(clusterProfiler)
    library(org.Bt.eg.db)
    library(KEGGREST)
    
    # 读取基因列表
    genes <- read.table("~{gene_list}", header=FALSE)$V1
    
    # 基因ID转换
    gene_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Bt.eg.db)
    
    # KEGG富集分析
    kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID,
                              organism = "~{species_code}",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)
    
    # 保存结果
    write.csv(as.data.frame(kegg_enrich), "kegg_enrichment_results.csv")
    
    # 生成可视化图
    pdf("kegg_dotplot.pdf", width=10, height=8)
    dotplot(kegg_enrich, showCategory=20)
    dev.off()
    
    pdf("kegg_barplot.pdf", width=10, height=8)
    barplot(kegg_enrich, showCategory=20)
    dev.off()
    EOF
  >>>

  output {
    File kegg_results = "kegg_enrichment_results.csv"
    File dotplot = "kegg_dotplot.pdf"
    File barplot = "kegg_barplot.pdf"
  }

  runtime {
    docker: "bioconductor/bioconductor_docker:RELEASE_3_16"
    cpu: 2
    memory: "8 GB"
  }
}
```

### 8.2 表型关联解释

**候选基因功能注释:**
```wdl
task CandidateGeneAnnotation {
  input {
    File gwas_results
    File gene_annotations
    Float pvalue_threshold = 5e-8
  }

  command <<<
    # Python脚本进行候选基因注释
    python3 << 'EOF'
    import pandas as pd
    import numpy as np
    
    # 读取GWAS结果
    gwas = pd.read_csv("~{gwas_results}", sep='\t')
    
    # 筛选显著位点
    significant = gwas[gwas['P'] < ~{pvalue_threshold}]
    
    # 读取基因注释
    annotations = pd.read_csv("~{gene_annotations}", sep='\t')
    
    # 基于位置匹配基因
    candidate_genes = []
    for _, snp in significant.iterrows():
        chr_genes = annotations[annotations['CHR'] == snp['CHR']]
        nearby_genes = chr_genes[
            (chr_genes['START'] <= snp['BP'] + 50000) & 
            (chr_genes['END'] >= snp['BP'] - 50000)
        ]
        
        for _, gene in nearby_genes.iterrows():
            candidate_genes.append({
                'SNP': snp['SNP'],
                'CHR': snp['CHR'],
                'BP': snp['BP'],
                'P': snp['P'],
                'GENE': gene['GENE_SYMBOL'],
                'GENE_NAME': gene['GENE_NAME'],
                'FUNCTION': gene['FUNCTION'],
                'DISTANCE': abs(gene['START'] - snp['BP'])
            })
    
    # 保存候选基因列表
    candidate_df = pd.DataFrame(candidate_genes)
    candidate_df.to_csv("candidate_genes.csv", index=False)
    
    # 生成基因功能总结
    with open("gene_function_summary.txt", "w") as f:
        f.write("Candidate Genes Functional Summary\n")
        f.write("="*50 + "\n\n")
        
        for gene in candidate_df['GENE'].unique():
            gene_info = candidate_df[candidate_df['GENE'] == gene].iloc[0]
            f.write(f"Gene: {gene}\n")
            f.write(f"Name: {gene_info['GENE_NAME']}\n")
            f.write(f"Function: {gene_info['FUNCTION']}\n")
            f.write(f"Associated SNPs: {len(candidate_df[candidate_df['GENE'] == gene])}\n")
            f.write("-" * 30 + "\n")
    EOF
  >>>

  output {
    File candidate_genes = "candidate_genes.csv"
    File function_summary = "gene_function_summary.txt"
  }

  runtime {
    docker: "python:3.9-slim"
    cpu: 2
    memory: "4 GB"
  }
}
```

### 8.3 可视化报告生成

**综合分析报告:**
```wdl
task GenerateReport {
  input {
    File gwas_results
    File candidate_genes
    File kegg_results
    File population_stats
    String trait_name
  }

  command <<<
    # R Markdown报告生成
    Rscript << 'EOF'
    library(rmarkdown)
    library(ggplot2)
    library(plotly)
    library(DT)
    library(knitr)
    
    # 创建R Markdown文档
    rmd_content <- '
---
title: "奶牛基因组关联分析报告"
subtitle: "性状: ~{trait_name}"
date: "`r Sys.Date()`"
output: 
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(ggplot2)
library(plotly)
library(DT)
library(dplyr)
```

## 分析概述

本报告展示了奶牛~{trait_name}性状的全基因组关联分析(GWAS)结果，包括：

- 显著关联位点识别
- 候选基因功能注释
- 生物学通路富集分析
- 群体遗传学参数

## GWAS结果

```{r gwas-manhattan}
# 读取GWAS结果
gwas_data <- read.table("~{gwas_results}", header=TRUE, stringsAsFactors=FALSE)

# 曼哈顿图
manhattan_plot <- ggplot(gwas_data, aes(x=BP, y=-log10(P), color=as.factor(CHR))) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=rep(c("darkblue", "orange"), 15)) +
  geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed") +
  labs(title="Manhattan Plot", x="Chromosome Position", y="-log10(P-value)") +
  theme_minimal() +
  theme(legend.position="none")

ggplotly(manhattan_plot)
```

```{r gwas-qq}
# QQ图
expected <- -log10(ppoints(nrow(gwas_data)))
observed <- -log10(sort(gwas_data$P))

qq_plot <- ggplot(data.frame(expected, observed), aes(x=expected, y=observed)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red") +
  labs(title="QQ Plot", x="Expected -log10(P)", y="Observed -log10(P)") +
  theme_minimal()

ggplotly(qq_plot)
```

## 候选基因

```{r candidate-genes}
# 候选基因表格
candidates <- read.csv("~{candidate_genes}")
datatable(candidates, options = list(pageLength = 10, scrollX = TRUE))
```

## 功能富集分析

```{r kegg-analysis}
# KEGG结果
kegg_data <- read.csv("~{kegg_results}")
if(nrow(kegg_data) > 0) {
  datatable(kegg_data, options = list(pageLength = 10, scrollX = TRUE))
} else {
  cat("No significant KEGG pathways found.")
}
```

## 群体遗传学统计

```{r population-stats}
# 群体统计信息
pop_stats <- read.table("~{population_stats}", header=TRUE)
summary(pop_stats)
```

## 结论

基于GWAS分析结果，我们识别出了与~{trait_name}性状显著关联的遗传变异位点。这些位点附近的候选基因主要参与以下生物学过程：

1. [根据实际结果填写主要生物学过程]
2. [相关代谢通路]
3. [可能的分子机制]

## 建议

1. **验证研究**: 建议在独立群体中验证显著关联位点
2. **功能研究**: 对候选基因进行功能验证实验
3. **育种应用**: 评估显著位点在分子标记辅助选择中的应用价值

'

    # 写入Rmd文件
    writeLines(rmd_content, "cattle_gwas_report.Rmd")
    
    # 渲染报告
    render("cattle_gwas_report.Rmd", output_file="cattle_gwas_report.html")
    EOF
  >>>

  output {
    File html_report = "cattle_gwas_report.html"
    File rmd_source = "cattle_gwas_report.Rmd"
  }

  runtime {
    docker: "rocker/tidyverse:4.2.0"
    cpu: 4
    memory: "8 GB"
  }
}
```

## 项目管理与监控

### 9.1 工作流编排

**主工作流定义:**
```wdl
version 1.0

workflow CattleGenomicsAnalysis {
  input {
    Array[File] fastq_r1_files
    Array[File] fastq_r2_files
    Array[String] sample_ids
    File reference_fasta
    File phenotype_file
    String trait_name
  }

  # 数据预处理
  scatter (i in range(length(sample_ids))) {
    call CattleQCPreprocessing {
      input:
        fastq_r1 = fastq_r1_files[i],
        fastq_r2 = fastq_r2_files[i],
        sample_id = sample_ids[i]
    }
  }

  # 序列比对
  scatter (i in range(length(sample_ids))) {
    call CattleAlignment {
      input:
        fastq_r1 = CattleQCPreprocessing.trimmed_r1[i],
        fastq_r2 = CattleQCPreprocessing.trimmed_r2[i],
        reference_fasta = reference_fasta,
        sample_id = sample_ids[i]
    }
  }

  # 变异检测
  call CattleVariantCalling {
    input:
      input_bams = CattleAlignment.final_bam,
      sample_ids = sample_ids,
      reference_fasta = reference_fasta
  }

  # 变异注释
  call VariantAnnotation {
    input:
      input_vcf = CattleVariantCalling.final_vcf
  }

  # 统计分析
  call GWASAnalysis {
    input:
      input_vcf = VariantAnnotation.final_annotated_vcf,
      phenotype_file = phenotype_file,
      trait_name = trait_name
  }

  # 生物学解释
  call KEGGPathwayAnalysis {
    input:
      gene_list = GWASAnalysis.candidate_genes
    }

  # 报告生成
  call GenerateReport {
    input:
      gwas_results = GWASAnalysis.gwas_results,
      kegg_results = KEGGPathwayAnalysis.kegg_results,
      trait_name = trait_name
  }

  output {
    File final_report = GenerateReport.html_report
    File annotated_vcf = VariantAnnotation.final_annotated_vcf
    File gwas_results = GWASAnalysis.gwas_results
  }
}
```

### 9.2 成本优化策略

**Spot实例使用:**
```json
{
  "computeEnvironment": {
    "type": "MANAGED",
    "state": "ENABLED",
    "computeResources": {
      "type": "EC2",
      "allocationStrategy": "SPOT_CAPACITY_OPTIMIZED",
      "minvCpus": 0,
      "maxvCpus": 1000,
      "desiredvCpus": 100,
      "instanceTypes": ["m5.large", "m5.xlarge", "c5.xlarge"],
      "spotIamFleetRequestRole": "arn:aws:iam::account:role/aws-ec2-spot-fleet-role",
      "bidPercentage": 50
    }
  }
}
```

**存储层级管理:**
```bash
# 自动化存储层级转换
aws s3api put-bucket-lifecycle-configuration \
  --bucket cattle-genomics-demo \
  --lifecycle-configuration file://lifecycle-policy.json
```

### 9.3 监控和告警

**CloudWatch监控设置:**
```bash
# 创建自定义指标
aws cloudwatch put-metric-data \
  --namespace "CattleGenomics" \
  --metric-data MetricName=WorkflowRuntime,Value=3600,Unit=Seconds

# 设置告警
aws cloudwatch put-metric-alarm \
  --alarm-name "LongRunningWorkflow" \
  --alarm-description "Workflow running too long" \
  --metric-name WorkflowRuntime \
  --namespace CattleGenomics \
  --statistic Average \
  --period 300 \
  --threshold 7200 \
  --comparison-operator GreaterThanThreshold
```

## 部署和运行指南

### 10.1 环境准备

**IAM角色配置:**
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "omics:*",
        "s3:GetObject",
        "s3:PutObject",
        "s3:ListBucket",
        "batch:SubmitJob",
        "batch:DescribeJobs",
        "logs:CreateLogGroup",
        "logs:CreateLogStream",
        "logs:PutLogEvents"
      ],
      "Resource": "*"
    }
  ]
}
```

**部署脚本:**
```bash
#!/bin/bash
# 部署奶牛基因组分析流水线

# 1. 创建S3存储桶
aws s3 mb s3://cattle-genomics-demo --region us-east-1

# 2. 上传工作流定义
aws s3 cp workflows/ s3://cattle-genomics-demo/workflows/ --recursive

# 3. 创建Omics工作流
WORKFLOW_ID=$(aws omics create-workflow \
  --name "cattle-genomics-pipeline" \
  --definition-uri s3://cattle-genomics-demo/workflows/main.wdl \
  --query 'id' --output text)

echo "Workflow created with ID: $WORKFLOW_ID"

# 4. 创建参考基因组存储
REF_STORE_ID=$(aws omics create-reference-store \
  --name "cattle-reference" \
  --query 'id' --output text)

echo "Reference store created with ID: $REF_STORE_ID"

# 5. 导入参考基因组
aws omics start-reference-import-job \
  --reference-store-id $REF_STORE_ID \
  --role-arn arn:aws:iam::account:role/OmicsRole \
  --sources file://reference-sources.json
```

### 10.2 运行示例

**批量提交任务:**
```bash
#!/bin/bash
# 批量运行奶牛基因组分析

WORKFLOW_ID="your-workflow-id"
ROLE_ARN="arn:aws:iam::account:role/OmicsRole"

# 读取样本清单
while IFS=',' read -r sample_id fastq_r1 fastq_r2 phenotype; do
  echo "Submitting analysis for sample: $sample_id"
  
  # 创建运行参数
  cat > ${sample_id}_run_params.json << EOF
{
  "CattleGenomicsAnalysis.sample_ids": ["$sample_id"],
  "CattleGenomicsAnalysis.fastq_r1_files": ["$fastq_r1"],
  "CattleGenomicsAnalysis.fastq_r2_files": ["$fastq_r2"],
  "CattleGenomicsAnalysis.reference_fasta": "s3://cattle-genomics-demo/reference/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
  "CattleGenomicsAnalysis.phenotype_file": "s3://cattle-genomics-demo/phenotypes/cattle_phenotypes.txt",
  "CattleGenomicsAnalysis.trait_name": "milk_yield"
}
EOF

  # 提交运行
  aws omics start-run \
    --workflow-id $WORKFLOW_ID \
    --workflow-type PRIVATE \
    --role-arn $ROLE_ARN \
    --name "cattle-analysis-${sample_id}" \
    --parameters file://${sample_id}_run_params.json \
    --output-uri s3://cattle-genomics-demo/results/${sample_id}/
    
done < sample_manifest.csv
```

## 预期结果和性能指标

### 11.1 分析结果

**预期输出文件:**
- 质控报告: `sample_qc_report.html`
- 比对统计: `alignment_metrics.txt`
- 变异文件: `cattle_filtered_annotated.vcf.gz`
- GWAS结果: `gwas_milk_yield.assoc.linear`
- 候选基因: `candidate_genes.csv`
- 功能富集: `kegg_enrichment_results.csv`
- 综合报告: `cattle_gwas_report.html`

### 11.2 性能基准

**处理时间估算:**
- 数据预处理: 2-4小时/样本
- 序列比对: 4-8小时/样本
- 变异检测: 2-6小时(100样本)
- 变异注释: 1-2小时
- 统计分析: 30分钟-2小时
- 总计: 24-48小时(100样本)

**成本估算:**
- 计算成本: $500-1000 (100样本)
- 存储成本: $50-100/月
- 数据传输: $20-50
- 总成本: $570-1150

### 11.3 质量控制指标

**数据质量标准:**
- 测序质量: Q30 > 85%
- 比对率: > 95%
- 重复率: < 20%
- 覆盖深度: 25-35X
- 变异检出率: > 99%

## 总结

本demo方案展示了如何使用AWS Omics构建完整的奶牛基因组分析流水线，实现了：

1. **可扩展性**: 支持从小规模到大规模样本分析
2. **标准化**: 使用WDL工作流确保分析的可重复性
3. **成本效益**: 通过Spot实例和存储优化降低成本
4. **自动化**: 端到端自动化分析流程
5. **可视化**: 生成交互式分析报告

该方案可作为农业基因组学研究的参考模板，并可根据具体需求进行定制和扩展。
