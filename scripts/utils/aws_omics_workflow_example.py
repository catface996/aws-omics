#!/usr/bin/env python3
"""
AWS Omics 工作流示例
展示如何使用模拟的测序数据创建和运行AWS Omics工作流
"""

import json
import argparse
from pathlib import Path
import boto3
from botocore.exceptions import ClientError

class OmicsWorkflowManager:
    def __init__(self, region='us-east-1'):
        self.region = region
        self.omics_client = boto3.client('omics', region_name=region)
        self.s3_client = boto3.client('s3', region_name=region)
        
    def create_reference_store(self, name, description):
        """创建参考基因组存储"""
        try:
            response = self.omics_client.create_reference_store(
                name=name,
                description=description
            )
            store_id = response['id']
            print(f"创建参考存储成功: {store_id}")
            return store_id
        except ClientError as e:
            print(f"创建参考存储失败: {e}")
            return None
    
    def import_reference_genome(self, store_id, name, source_uri):
        """导入参考基因组"""
        try:
            response = self.omics_client.create_reference(
                referenceStoreId=store_id,
                name=name,
                sourceUri=source_uri
            )
            reference_id = response['id']
            print(f"导入参考基因组成功: {reference_id}")
            return reference_id
        except ClientError as e:
            print(f"导入参考基因组失败: {e}")
            return None
    
    def create_sequence_store(self, name, description):
        """创建序列数据存储"""
        try:
            response = self.omics_client.create_sequence_store(
                name=name,
                description=description
            )
            store_id = response['id']
            print(f"创建序列存储成功: {store_id}")
            return store_id
        except ClientError as e:
            print(f"创建序列存储失败: {e}")
            return None
    
    def upload_data_to_s3(self, local_file, bucket, key):
        """上传数据到S3"""
        try:
            self.s3_client.upload_file(local_file, bucket, key)
            s3_uri = f"s3://{bucket}/{key}"
            print(f"上传文件成功: {s3_uri}")
            return s3_uri
        except ClientError as e:
            print(f"上传文件失败: {e}")
            return None
    
    def create_variant_calling_workflow(self):
        """创建变异检测工作流的WDL定义"""
        wdl_content = '''
version 1.0

workflow VariantCallingWorkflow {
    input {
        File reference_genome
        File input_fastq
        String sample_name
        String output_prefix
    }
    
    call AlignReads {
        input:
            reference = reference_genome,
            fastq = input_fastq,
            sample_name = sample_name
    }
    
    call CallVariants {
        input:
            reference = reference_genome,
            aligned_bam = AlignReads.output_bam,
            sample_name = sample_name
    }
    
    call AnnotateVariants {
        input:
            variants_vcf = CallVariants.output_vcf,
            sample_name = sample_name,
            output_prefix = output_prefix
    }
    
    output {
        File aligned_bam = AlignReads.output_bam
        File variants_vcf = CallVariants.output_vcf
        File annotated_vcf = AnnotateVariants.output_vcf
        File analysis_report = AnnotateVariants.report
    }
}

task AlignReads {
    input {
        File reference
        File fastq
        String sample_name
    }
    
    command <<<
        # 建立BWA索引
        bwa index ~{reference}
        
        # 序列比对
        bwa mem -t 4 -R "@RG\\tID:~{sample_name}\\tSM:~{sample_name}\\tPL:ILLUMINA" \\
            ~{reference} ~{fastq} | \\
            samtools view -Sb - | \\
            samtools sort -o ~{sample_name}.sorted.bam -
        
        # 建立索引
        samtools index ~{sample_name}.sorted.bam
    >>>
    
    output {
        File output_bam = "~{sample_name}.sorted.bam"
        File output_bai = "~{sample_name}.sorted.bam.bai"
    }
    
    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
        cpu: 4
        memory: "8 GB"
        disks: "local-disk 100 SSD"
    }
}

task CallVariants {
    input {
        File reference
        File aligned_bam
        String sample_name
    }
    
    command <<<
        # 使用GATK进行变异检测
        gatk HaplotypeCaller \\
            -R ~{reference} \\
            -I ~{aligned_bam} \\
            -O ~{sample_name}.variants.vcf \\
            --emit-ref-confidence GVCF
    >>>
    
    output {
        File output_vcf = "~{sample_name}.variants.vcf"
    }
    
    runtime {
        docker: "broadinstitute/gatk:latest"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 50 SSD"
    }
}

task AnnotateVariants {
    input {
        File variants_vcf
        String sample_name
        String output_prefix
    }
    
    command <<<
        # 变异注释和统计
        bcftools stats ~{variants_vcf} > ~{output_prefix}.stats.txt
        
        # 生成注释VCF
        cp ~{variants_vcf} ~{output_prefix}.annotated.vcf
        
        # 生成分析报告
        echo "Variant Analysis Report for ~{sample_name}" > ~{output_prefix}.report.txt
        echo "Generated on: $(date)" >> ~{output_prefix}.report.txt
        echo "" >> ~{output_prefix}.report.txt
        echo "Variant Statistics:" >> ~{output_prefix}.report.txt
        bcftools stats ~{variants_vcf} | grep "^SN" >> ~{output_prefix}.report.txt
    >>>
    
    output {
        File output_vcf = "~{output_prefix}.annotated.vcf"
        File report = "~{output_prefix}.report.txt"
        File stats = "~{output_prefix}.stats.txt"
    }
    
    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 SSD"
    }
}
'''
        return wdl_content
    
    def create_workflow(self, name, definition_uri, parameter_template=None):
        """创建Omics工作流"""
        try:
            params = {
                'name': name,
                'definitionUri': definition_uri
            }
            
            if parameter_template:
                params['parameterTemplate'] = parameter_template
            
            response = self.omics_client.create_workflow(**params)
            workflow_id = response['id']
            print(f"创建工作流成功: {workflow_id}")
            return workflow_id
        except ClientError as e:
            print(f"创建工作流失败: {e}")
            return None
    
    def start_workflow_run(self, workflow_id, name, role_arn, parameters, output_uri):
        """启动工作流运行"""
        try:
            response = self.omics_client.start_run(
                workflowId=workflow_id,
                name=name,
                roleArn=role_arn,
                parameters=parameters,
                outputUri=output_uri
            )
            run_id = response['id']
            print(f"启动工作流运行成功: {run_id}")
            return run_id
        except ClientError as e:
            print(f"启动工作流运行失败: {e}")
            return None
    
    def get_run_status(self, run_id):
        """获取工作流运行状态"""
        try:
            response = self.omics_client.get_run(id=run_id)
            return response['status']
        except ClientError as e:
            print(f"获取运行状态失败: {e}")
            return None

def create_deployment_script(output_file):
    """创建部署脚本"""
    script_content = '''#!/bin/bash

# AWS Omics 奶牛基因组分析部署脚本
# 使用模拟测序数据

set -e

# 配置变量
PROJECT_NAME="cattle-genomics-simulation"
REGION="us-east-1"
BUCKET_NAME="cattle-genomics-sim-$(date +%s)"
WORKFLOW_NAME="cattle-variant-calling"
REFERENCE_STORE_NAME="cattle-reference-store"
SEQUENCE_STORE_NAME="cattle-sequence-store"

# 颜色输出
GREEN='\\033[0;32m'
YELLOW='\\033[1;33m'
RED='\\033[0;31m'
NC='\\033[0m'

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
        log_error "AWS CLI未安装"
        exit 1
    fi
    
    if ! aws sts get-caller-identity &> /dev/null; then
        log_error "AWS凭证未配置"
        exit 1
    fi
    
    ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
    log_info "当前AWS账户: $ACCOUNT_ID"
}

# 创建S3存储桶
create_s3_bucket() {
    log_info "创建S3存储桶: $BUCKET_NAME"
    
    aws s3 mb s3://$BUCKET_NAME --region $REGION
    
    # 设置生命周期策略
    cat > lifecycle.json << EOF
{
    "Rules": [
        {
            "ID": "GenomicsDataLifecycle",
            "Status": "Enabled",
            "Transitions": [
                {
                    "Days": 30,
                    "StorageClass": "STANDARD_IA"
                },
                {
                    "Days": 90,
                    "StorageClass": "GLACIER"
                }
            ]
        }
    ]
}
EOF
    
    aws s3api put-bucket-lifecycle-configuration \\
        --bucket $BUCKET_NAME \\
        --lifecycle-configuration file://lifecycle.json
    
    rm lifecycle.json
}

# 上传数据到S3
upload_data() {
    log_info "上传基因组数据到S3..."
    
    # 上传参考基因组
    aws s3 cp genomic_data/reference/GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz \\
        s3://$BUCKET_NAME/reference/genome.fna.gz
    
    # 上传模拟测序数据
    if [ -f "simulated_data/simulated_reads.fastq" ]; then
        aws s3 cp simulated_data/simulated_reads.fastq \\
            s3://$BUCKET_NAME/reads/sample1.fastq
    else
        log_warn "模拟测序数据不存在，请先运行 simulate_sequencing_data.py"
    fi
    
    # 上传工作流定义
    aws s3 cp variant_calling_workflow.wdl \\
        s3://$BUCKET_NAME/workflows/variant_calling.wdl
}

# 创建IAM角色
create_iam_role() {
    log_info "创建IAM角色..."
    
    # Omics服务角色
    cat > omics-trust-policy.json << EOF
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

    aws iam create-role \\
        --role-name OmicsExecutionRole \\
        --assume-role-policy-document file://omics-trust-policy.json || true
    
    # 附加策略
    aws iam attach-role-policy \\
        --role-name OmicsExecutionRole \\
        --policy-arn arn:aws:iam::aws:policy/service-role/AmazonOmicsServiceRolePolicy || true
    
    rm omics-trust-policy.json
}

# 创建Omics资源
create_omics_resources() {
    log_info "创建Omics资源..."
    
    # 创建参考存储
    REFERENCE_STORE_ID=$(aws omics create-reference-store \\
        --name $REFERENCE_STORE_NAME \\
        --description "Cattle genome reference store" \\
        --query 'id' --output text)
    
    log_info "参考存储ID: $REFERENCE_STORE_ID"
    
    # 创建序列存储
    SEQUENCE_STORE_ID=$(aws omics create-sequence-store \\
        --name $SEQUENCE_STORE_NAME \\
        --description "Cattle sequence data store" \\
        --query 'id' --output text)
    
    log_info "序列存储ID: $SEQUENCE_STORE_ID"
    
    # 导入参考基因组
    aws omics create-reference \\
        --reference-store-id $REFERENCE_STORE_ID \\
        --name "ARS-UCD1.2" \\
        --source-uri s3://$BUCKET_NAME/reference/genome.fna.gz
}

# 主函数
main() {
    echo "=========================================="
    echo "  AWS Omics 奶牛基因组分析部署"
    echo "=========================================="
    
    check_aws_config
    create_s3_bucket
    upload_data
    create_iam_role
    create_omics_resources
    
    log_info "部署完成！"
    echo "S3存储桶: $BUCKET_NAME"
    echo "参考存储: $REFERENCE_STORE_ID"
    echo "序列存储: $SEQUENCE_STORE_ID"
}

main "$@"
'''
    
    with open(output_file, 'w') as f:
        f.write(script_content)
    
    print(f"部署脚本已创建: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='AWS Omics工作流管理')
    parser.add_argument('--action', choices=['create-deployment-script', 'create-workflow-def'],
                       default='create-deployment-script',
                       help='执行的操作')
    parser.add_argument('--region', default='us-east-1',
                       help='AWS区域')
    
    args = parser.parse_args()
    
    if args.action == 'create-deployment-script':
        create_deployment_script('deploy_omics_simulation.sh')
        
    elif args.action == 'create-workflow-def':
        manager = OmicsWorkflowManager(args.region)
        wdl_content = manager.create_variant_calling_workflow()
        
        with open('variant_calling_workflow.wdl', 'w') as f:
            f.write(wdl_content)
        
        print("WDL工作流定义已创建: variant_calling_workflow.wdl")
    
    print("\n使用说明:")
    print("1. 运行模拟数据生成: python3 simulate_sequencing_data.py")
    print("2. 创建工作流定义: python3 aws_omics_workflow_example.py --action create-workflow-def")
    print("3. 运行部署脚本: chmod +x deploy_omics_simulation.sh && ./deploy_omics_simulation.sh")
    print("4. 在AWS控制台中监控工作流执行")

if __name__ == "__main__":
    main()
