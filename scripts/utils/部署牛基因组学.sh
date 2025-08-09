#!/bin/bash

# AWS Omics 奶牛基因组分析流水线部署脚本
# 作者: AWS Omics Demo Team
# 日期: $(date +%Y-%m-%d)

set -e

# 配置变量
PROJECT_NAME="cattle-genomics-demo"
REGION="us-east-1"
BUCKET_NAME="cattle-genomics-demo-$(date +%s)"
WORKFLOW_NAME="cattle-genomics-pipeline"
REFERENCE_STORE_NAME="cattle-reference-store"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 检查AWS CLI配置
check_aws_config() {
    log_info "检查AWS CLI配置..."
    
    if ! command -v aws &> /dev/null; then
        log_error "AWS CLI未安装，请先安装AWS CLI"
        exit 1
    fi
    
    if ! aws sts get-caller-identity &> /dev/null; then
        log_error "AWS凭证未配置，请运行 'aws configure'"
        exit 1
    fi
    
    ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
    log_info "当前AWS账户: $ACCOUNT_ID"
}

# 创建IAM角色
create_iam_roles() {
    log_info "创建IAM角色..."
    
    # Omics工作流执行角色
    cat > omics-workflow-role-trust-policy.json << EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "omics.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

    cat > omics-workflow-role-policy.json << EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:ListBucket",
        "s3:DeleteObject"
      ],
      "Resource": [
        "arn:aws:s3:::${BUCKET_NAME}",
        "arn:aws:s3:::${BUCKET_NAME}/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "batch:SubmitJob",
        "batch:DescribeJobs",
        "batch:TerminateJob"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "logs:CreateLogGroup",
        "logs:CreateLogStream",
        "logs:PutLogEvents",
        "logs:DescribeLogGroups",
        "logs:DescribeLogStreams"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "ecr:GetAuthorizationToken",
        "ecr:BatchCheckLayerAvailability",
        "ecr:GetDownloadUrlForLayer",
        "ecr:BatchGetImage"
      ],
      "Resource": "*"
    }
  ]
}
EOF

    # 创建角色
    aws iam create-role \
        --role-name OmicsWorkflowExecutionRole \
        --assume-role-policy-document file://omics-workflow-role-trust-policy.json \
        --description "Role for AWS Omics workflow execution" \
        2>/dev/null || log_warn "角色可能已存在"

    # 附加策略
    aws iam put-role-policy \
        --role-name OmicsWorkflowExecutionRole \
        --policy-name OmicsWorkflowExecutionPolicy \
        --policy-document file://omics-workflow-role-policy.json

    WORKFLOW_ROLE_ARN="arn:aws:iam::${ACCOUNT_ID}:role/OmicsWorkflowExecutionRole"
    log_info "工作流执行角色ARN: $WORKFLOW_ROLE_ARN"
    
    # 清理临时文件
    rm -f omics-workflow-role-trust-policy.json omics-workflow-role-policy.json
}

# 创建S3存储桶
create_s3_bucket() {
    log_info "创建S3存储桶: $BUCKET_NAME"
    
    if aws s3 mb s3://$BUCKET_NAME --region $REGION; then
        log_info "S3存储桶创建成功"
    else
        log_error "S3存储桶创建失败"
        exit 1
    fi
    
    # 配置生命周期策略
    cat > lifecycle-policy.json << EOF
{
  "Rules": [
    {
      "ID": "CattleGenomicsLifecycle",
      "Status": "Enabled",
      "Filter": {
        "Prefix": "processed/"
      },
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
EOF

    aws s3api put-bucket-lifecycle-configuration \
        --bucket $BUCKET_NAME \
        --lifecycle-configuration file://lifecycle-policy.json
    
    log_info "生命周期策略配置完成"
    rm -f lifecycle-policy.json
}

# 创建工作流定义文件
create_workflow_files() {
    log_info "创建工作流定义文件..."
    
    mkdir -p workflows
    
    # 主工作流文件
    cat > workflows/main.wdl << 'EOF'
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
    call QCPreprocessing {
      input:
        fastq_r1 = fastq_r1_files[i],
        fastq_r2 = fastq_r2_files[i],
        sample_id = sample_ids[i]
    }
  }

  # 序列比对
  scatter (i in range(length(sample_ids))) {
    call BWAAlignment {
      input:
        fastq_r1 = QCPreprocessing.trimmed_r1[i],
        fastq_r2 = QCPreprocessing.trimmed_r2[i],
        reference_fasta = reference_fasta,
        sample_id = sample_ids[i]
    }
  }

  # 变异检测
  call VariantCalling {
    input:
      input_bams = BWAAlignment.final_bam,
      sample_ids = sample_ids,
      reference_fasta = reference_fasta
  }

  output {
    File final_vcf = VariantCalling.filtered_vcf
    Array[File] aligned_bams = BWAAlignment.final_bam
    Array[File] qc_reports = QCPreprocessing.qc_report
  }
}

task QCPreprocessing {
  input {
    File fastq_r1
    File fastq_r2
    String sample_id
  }

  command <<<
    # FastQC质量控制
    fastqc ~{fastq_r1} ~{fastq_r2} -o ./
    
    # Trimmomatic质量修剪
    trimmomatic PE -threads 4 \
      ~{fastq_r1} ~{fastq_r2} \
      ~{sample_id}_trimmed_R1.fastq.gz ~{sample_id}_unpaired_R1.fastq.gz \
      ~{sample_id}_trimmed_R2.fastq.gz ~{sample_id}_unpaired_R2.fastq.gz \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # 生成质控报告
    multiqc . -n ~{sample_id}_qc_report
  >>>

  output {
    File trimmed_r1 = "${sample_id}_trimmed_R1.fastq.gz"
    File trimmed_r2 = "${sample_id}_trimmed_R2.fastq.gz"
    File qc_report = "${sample_id}_qc_report.html"
  }

  runtime {
    docker: "quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0"
    cpu: 4
    memory: "8 GB"
  }
}

task BWAAlignment {
  input {
    File fastq_r1
    File fastq_r2
    File reference_fasta
    String sample_id
  }

  command <<<
    # BWA比对
    bwa mem -t 8 -R "@RG\tID:~{sample_id}\tSM:~{sample_id}\tPL:ILLUMINA" \
      ~{reference_fasta} ~{fastq_r1} ~{fastq_r2} | \
      samtools sort -@ 4 -o ~{sample_id}_sorted.bam -
    
    # 创建索引
    samtools index ~{sample_id}_sorted.bam
    
    # 标记重复
    picard MarkDuplicates \
      I=~{sample_id}_sorted.bam \
      O=~{sample_id}_dedup.bam \
      M=~{sample_id}_metrics.txt
    
    samtools index ~{sample_id}_dedup.bam
  >>>

  output {
    File final_bam = "${sample_id}_dedup.bam"
    File final_bai = "${sample_id}_dedup.bam.bai"
    File metrics = "${sample_id}_metrics.txt"
  }

  runtime {
    docker: "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    cpu: 8
    memory: "16 GB"
  }
}

task VariantCalling {
  input {
    Array[File] input_bams
    Array[String] sample_ids
    File reference_fasta
  }

  command <<<
    # GATK HaplotypeCaller
    gatk HaplotypeCaller \
      -R ~{reference_fasta} \
      ~{sep=' -I ' input_bams} \
      -O raw_variants.vcf.gz \
      --native-pair-hmm-threads 4
    
    # 变异过滤
    gatk VariantFiltration \
      -R ~{reference_fasta} \
      -V raw_variants.vcf.gz \
      -O filtered_variants.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "FS > 60.0" --filter-name "FS60"
  >>>

  output {
    File raw_vcf = "raw_variants.vcf.gz"
    File filtered_vcf = "filtered_variants.vcf.gz"
  }

  runtime {
    docker: "broadinstitute/gatk:4.3.0.0"
    cpu: 8
    memory: "32 GB"
  }
}
EOF

    # 参数模板文件
    cat > workflows/parameters.json << EOF
{
  "CattleGenomicsAnalysis.fastq_r1_files": {
    "description": "Array of forward reads FASTQ files"
  },
  "CattleGenomicsAnalysis.fastq_r2_files": {
    "description": "Array of reverse reads FASTQ files"
  },
  "CattleGenomicsAnalysis.sample_ids": {
    "description": "Array of sample identifiers"
  },
  "CattleGenomicsAnalysis.reference_fasta": {
    "description": "Reference genome FASTA file"
  },
  "CattleGenomicsAnalysis.phenotype_file": {
    "description": "Phenotype data file"
  },
  "CattleGenomicsAnalysis.trait_name": {
    "description": "Name of the trait to analyze"
  }
}
EOF

    log_info "工作流定义文件创建完成"
}

# 上传工作流文件到S3
upload_workflows() {
    log_info "上传工作流文件到S3..."
    
    aws s3 cp workflows/ s3://$BUCKET_NAME/workflows/ --recursive
    log_info "工作流文件上传完成"
}

# 创建Omics工作流
create_omics_workflow() {
    log_info "创建AWS Omics工作流..."
    
    WORKFLOW_ID=$(aws omics create-workflow \
        --name "$WORKFLOW_NAME" \
        --description "Cattle genomics analysis pipeline using AWS Omics" \
        --definition-uri "s3://$BUCKET_NAME/workflows/main.wdl" \
        --parameter-template file://workflows/parameters.json \
        --query 'id' --output text)
    
    if [ $? -eq 0 ]; then
        log_info "工作流创建成功，ID: $WORKFLOW_ID"
        echo "WORKFLOW_ID=$WORKFLOW_ID" > .env
    else
        log_error "工作流创建失败"
        exit 1
    fi
}

# 创建参考基因组存储
create_reference_store() {
    log_info "创建参考基因组存储..."
    
    REF_STORE_ID=$(aws omics create-reference-store \
        --name "$REFERENCE_STORE_NAME" \
        --description "Reference genomes for cattle genomics analysis" \
        --query 'id' --output text)
    
    if [ $? -eq 0 ]; then
        log_info "参考基因组存储创建成功，ID: $REF_STORE_ID"
        echo "REF_STORE_ID=$REF_STORE_ID" >> .env
    else
        log_error "参考基因组存储创建失败"
        exit 1
    fi
}

# 创建示例数据文件
create_sample_data() {
    log_info "创建示例数据文件..."
    
    mkdir -p sample_data
    
    # 样本清单文件
    cat > sample_data/sample_manifest.csv << EOF
sample_id,fastq_r1,fastq_r2,phenotype
sample001,s3://$BUCKET_NAME/raw-data/sample001/sample001_R1.fastq.gz,s3://$BUCKET_NAME/raw-data/sample001/sample001_R2.fastq.gz,milk_yield
sample002,s3://$BUCKET_NAME/raw-data/sample002/sample002_R1.fastq.gz,s3://$BUCKET_NAME/raw-data/sample002/sample002_R2.fastq.gz,milk_yield
sample003,s3://$BUCKET_NAME/raw-data/sample003/sample003_R1.fastq.gz,s3://$BUCKET_NAME/raw-data/sample003/sample003_R2.fastq.gz,milk_yield
EOF

    # 表型数据文件
    cat > sample_data/phenotypes.txt << EOF
FID	IID	milk_yield	fat_content	protein_content
sample001	sample001	8500	3.8	3.2
sample002	sample002	9200	4.1	3.4
sample003	sample003	7800	3.6	3.0
EOF

    # 运行参数示例
    cat > sample_data/run_parameters.json << EOF
{
  "CattleGenomicsAnalysis.fastq_r1_files": [
    "s3://$BUCKET_NAME/raw-data/sample001/sample001_R1.fastq.gz",
    "s3://$BUCKET_NAME/raw-data/sample002/sample002_R1.fastq.gz",
    "s3://$BUCKET_NAME/raw-data/sample003/sample003_R1.fastq.gz"
  ],
  "CattleGenomicsAnalysis.fastq_r2_files": [
    "s3://$BUCKET_NAME/raw-data/sample001/sample001_R2.fastq.gz",
    "s3://$BUCKET_NAME/raw-data/sample002/sample002_R2.fastq.gz",
    "s3://$BUCKET_NAME/raw-data/sample003/sample003_R2.fastq.gz"
  ],
  "CattleGenomicsAnalysis.sample_ids": [
    "sample001",
    "sample002", 
    "sample003"
  ],
  "CattleGenomicsAnalysis.reference_fasta": "s3://$BUCKET_NAME/reference/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
  "CattleGenomicsAnalysis.phenotype_file": "s3://$BUCKET_NAME/phenotypes/phenotypes.txt",
  "CattleGenomicsAnalysis.trait_name": "milk_yield"
}
EOF

    # 上传示例数据到S3
    aws s3 cp sample_data/phenotypes.txt s3://$BUCKET_NAME/phenotypes/
    
    log_info "示例数据文件创建完成"
}

# 创建运行脚本
create_run_script() {
    log_info "创建运行脚本..."
    
    cat > run_analysis.sh << EOF
#!/bin/bash

# 奶牛基因组分析运行脚本
# 使用方法: ./run_analysis.sh [参数文件]

set -e

# 加载环境变量
if [ -f .env ]; then
    source .env
else
    echo "错误: .env文件不存在，请先运行部署脚本"
    exit 1
fi

PARAMS_FILE=\${1:-sample_data/run_parameters.json}
RUN_NAME="cattle-analysis-\$(date +%Y%m%d-%H%M%S)"

echo "提交AWS Omics运行..."
echo "工作流ID: \$WORKFLOW_ID"
echo "参数文件: \$PARAMS_FILE"
echo "运行名称: \$RUN_NAME"

RUN_ID=\$(aws omics start-run \\
    --workflow-id \$WORKFLOW_ID \\
    --workflow-type PRIVATE \\
    --role-arn arn:aws:iam::\$(aws sts get-caller-identity --query Account --output text):role/OmicsWorkflowExecutionRole \\
    --name "\$RUN_NAME" \\
    --parameters file://\$PARAMS_FILE \\
    --output-uri s3://$BUCKET_NAME/results/\$RUN_NAME/ \\
    --query 'id' --output text)

if [ \$? -eq 0 ]; then
    echo "运行提交成功！"
    echo "运行ID: \$RUN_ID"
    echo "监控运行状态: aws omics get-run --id \$RUN_ID"
    echo "查看日志: aws logs describe-log-groups --log-group-name-prefix /aws/omics/WorkflowLog"
else
    echo "运行提交失败"
    exit 1
fi
EOF

    chmod +x run_analysis.sh
    log_info "运行脚本创建完成: run_analysis.sh"
}

# 创建监控脚本
create_monitoring_script() {
    log_info "创建监控脚本..."
    
    cat > monitor_runs.sh << EOF
#!/bin/bash

# 监控AWS Omics运行状态

echo "AWS Omics 运行状态监控"
echo "========================"

# 列出所有运行
echo "最近的运行:"
aws omics list-runs --query 'items[0:5].[id,name,status,creationTime]' --output table

echo ""
echo "输入运行ID查看详细状态 (按Enter跳过):"
read -r RUN_ID

if [ ! -z "\$RUN_ID" ]; then
    echo "运行详细信息:"
    aws omics get-run --id \$RUN_ID --query '{ID:id,Name:name,Status:status,CreationTime:creationTime,StartTime:startTime,StopTime:stopTime}' --output table
    
    echo ""
    echo "运行任务状态:"
    aws omics list-run-tasks --id \$RUN_ID --query 'items[].[taskId,name,status,creationTime]' --output table
fi
EOF

    chmod +x monitor_runs.sh
    log_info "监控脚本创建完成: monitor_runs.sh"
}

# 生成部署总结
generate_summary() {
    log_info "生成部署总结..."
    
    cat > DEPLOYMENT_SUMMARY.md << EOF
# AWS Omics 奶牛基因组分析流水线部署总结

## 部署信息
- **项目名称**: $PROJECT_NAME
- **AWS区域**: $REGION
- **S3存储桶**: $BUCKET_NAME
- **工作流名称**: $WORKFLOW_NAME
- **部署时间**: $(date)

## 创建的资源
1. **IAM角色**: OmicsWorkflowExecutionRole
2. **S3存储桶**: $BUCKET_NAME
3. **Omics工作流**: $WORKFLOW_NAME
4. **参考基因组存储**: $REFERENCE_STORE_NAME

## 文件结构
\`\`\`
.
├── workflows/
│   ├── main.wdl              # 主工作流定义
│   └── parameters.json       # 参数模板
├── sample_data/
│   ├── sample_manifest.csv   # 样本清单
│   ├── phenotypes.txt        # 表型数据
│   └── run_parameters.json   # 运行参数示例
├── run_analysis.sh           # 运行脚本
├── monitor_runs.sh           # 监控脚本
├── .env                      # 环境变量
└── DEPLOYMENT_SUMMARY.md     # 本文件
\`\`\`

## 使用方法

### 1. 上传测序数据
\`\`\`bash
# 上传FASTQ文件到指定位置
aws s3 cp sample001_R1.fastq.gz s3://$BUCKET_NAME/raw-data/sample001/
aws s3 cp sample001_R2.fastq.gz s3://$BUCKET_NAME/raw-data/sample001/
\`\`\`

### 2. 上传参考基因组
\`\`\`bash
# 下载奶牛参考基因组
wget http://ftp.ensembl.org/pub/release-107/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz
gunzip Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz

# 上传到S3
aws s3 cp Bos_taurus.ARS-UCD1.2.dna.toplevel.fa s3://$BUCKET_NAME/reference/
\`\`\`

### 3. 运行分析
\`\`\`bash
# 使用默认参数运行
./run_analysis.sh

# 使用自定义参数文件运行
./run_analysis.sh my_parameters.json
\`\`\`

### 4. 监控运行状态
\`\`\`bash
./monitor_runs.sh
\`\`\`

## 成本估算
- **计算成本**: \$10-50 per sample (取决于数据大小和分析复杂度)
- **存储成本**: \$0.023 per GB per month (S3 Standard)
- **数据传输**: \$0.09 per GB (出站传输)

## 注意事项
1. 确保有足够的AWS服务限额
2. 监控成本使用情况
3. 定期清理不需要的中间文件
4. 备份重要的分析结果

## 支持和故障排除
- 查看CloudWatch日志获取详细错误信息
- 检查IAM权限配置
- 验证S3存储桶访问权限
- 确认工作流定义语法正确

## 清理资源
如需删除所有创建的资源，请运行：
\`\`\`bash
# 删除S3存储桶内容
aws s3 rm s3://$BUCKET_NAME --recursive

# 删除S3存储桶
aws s3 rb s3://$BUCKET_NAME

# 删除IAM角色
aws iam delete-role-policy --role-name OmicsWorkflowExecutionRole --policy-name OmicsWorkflowExecutionPolicy
aws iam delete-role --role-name OmicsWorkflowExecutionRole

# 删除Omics工作流和参考存储需要通过控制台操作
\`\`\`
EOF

    log_info "部署总结已保存到 DEPLOYMENT_SUMMARY.md"
}

# 主函数
main() {
    echo "=========================================="
    echo "AWS Omics 奶牛基因组分析流水线部署脚本"
    echo "=========================================="
    echo ""
    
    check_aws_config
    create_iam_roles
    create_s3_bucket
    create_workflow_files
    upload_workflows
    create_omics_workflow
    create_reference_store
    create_sample_data
    create_run_script
    create_monitoring_script
    generate_summary
    
    echo ""
    log_info "=========================================="
    log_info "部署完成！"
    log_info "=========================================="
    log_info "S3存储桶: $BUCKET_NAME"
    log_info "工作流ID: $(cat .env | grep WORKFLOW_ID | cut -d'=' -f2)"
    log_info "参考存储ID: $(cat .env | grep REF_STORE_ID | cut -d'=' -f2)"
    log_info ""
    log_info "下一步操作:"
    log_info "1. 上传测序数据到 s3://$BUCKET_NAME/raw-data/"
    log_info "2. 上传参考基因组到 s3://$BUCKET_NAME/reference/"
    log_info "3. 运行 ./run_analysis.sh 开始分析"
    log_info "4. 使用 ./monitor_runs.sh 监控运行状态"
    log_info ""
    log_info "详细信息请查看 DEPLOYMENT_SUMMARY.md"
    echo ""
}

# 执行主函数
main "$@"
