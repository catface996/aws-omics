version 1.0

## 奶牛基因组数据预处理工作流 - Read Set版本
## 直接处理AWS Omics Sequence Store中的Read Set
## 包含质量评估、接头去除、质量过滤、长度过滤和去重复序列
##
## 作者: AWS Omics 奶牛基因组学演示项目
## 版本: 2.0 (Read Set支持)
## 日期: 2025-08-09

import "tasks/fastqc_task.wdl" as FastQC
import "tasks/trimmomatic_task.wdl" as Trimmomatic
import "tasks/fastp_task.wdl" as FastP
import "tasks/deduplication_task.wdl" as Dedup
import "tasks/multiqc_task.wdl" as MultiQC

workflow PreprocessingWorkflowReadSet {
    meta {
        description: "完整的测序数据预处理工作流，直接处理Read Set，包含质量控制、接头去除、质量过滤、长度过滤和去重复"
        version: "2.0"
        author: "AWS Omics Demo Project"
        keywords: ["preprocessing", "quality-control", "readset", "omics"]
    }

    parameter_meta {
        sample_name: "样本名称"
        input_readset: "输入的Read Set ARN或ID"
        reference_genome: "参考基因组ARN或S3路径"
        adapter_fasta: "接头序列文件（可选）"
        min_length: "最小读长阈值"
        min_quality: "最小质量分数"
        threads: "使用的CPU线程数"
        paired_end: "是否为双端测序数据"
    }

    input {
        String sample_name
        String input_readset
        String reference_genome
        File? adapter_fasta
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Boolean paired_end = true
        
        # 高级参数
        Int max_length = 500
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
        String dedup_method = "fastuniq"
        
        # 资源配置
        Int fastqc_memory_gb = 8
        Int trimmomatic_memory_gb = 16
        Int fastp_memory_gb = 16
        Int dedup_memory_gb = 16
        Int multiqc_memory_gb = 8
    }

    # 步骤1: 从Read Set提取FASTQ文件
    call ExtractReadSet {
        input:
            readset_arn = input_readset,
            sample_name = sample_name,
            paired_end = paired_end
    }

    # 步骤2: 原始数据质量评估
    call FastQC.RunFastQC as InitialQC {
        input:
            fastq_files = ExtractReadSet.fastq_files,
            sample_name = sample_name + "_initial",
            memory_gb = fastqc_memory_gb
    }

    # 步骤3: 使用FastP进行质量过滤和接头去除
    if (paired_end) {
        call FastP.RunFastPPE {
            input:
                fastq_r1 = ExtractReadSet.fastq_files[0],
                fastq_r2 = ExtractReadSet.fastq_files[1],
                sample_name = sample_name,
                min_length = min_length,
                min_quality = min_quality,
                max_length = max_length,
                complexity_threshold = complexity_threshold,
                enable_polyg_trimming = enable_polyg_trimming,
                enable_polyx_trimming = enable_polyx_trimming,
                threads = threads,
                memory_gb = fastp_memory_gb
        }
    }

    if (!paired_end) {
        call FastP.RunFastPSE {
            input:
                fastq = ExtractReadSet.fastq_files[0],
                sample_name = sample_name,
                min_length = min_length,
                min_quality = min_quality,
                max_length = max_length,
                complexity_threshold = complexity_threshold,
                enable_polyg_trimming = enable_polyg_trimming,
                enable_polyx_trimming = enable_polyx_trimming,
                threads = threads,
                memory_gb = fastp_memory_gb
        }
    }

    # 选择FastP输出
    Array[File] fastp_output = if paired_end then 
        select_all([RunFastPPE.cleaned_fastq_r1, RunFastPPE.cleaned_fastq_r2]) 
        else [select_first([RunFastPSE.cleaned_fastq])]

    # 步骤4: 去重复序列
    if (paired_end) {
        call Dedup.RemoveDuplicatesPE {
            input:
                fastq_r1 = fastp_output[0],
                fastq_r2 = fastp_output[1],
                sample_name = sample_name,
                method = dedup_method,
                memory_gb = dedup_memory_gb
        }
    }

    if (!paired_end) {
        call Dedup.RemoveDuplicatesSE {
            input:
                fastq = fastp_output[0],
                sample_name = sample_name,
                method = dedup_method,
                memory_gb = dedup_memory_gb
        }
    }

    # 选择去重输出
    Array[File] final_fastq = if paired_end then 
        select_all([RemoveDuplicatesPE.dedup_fastq_r1, RemoveDuplicatesPE.dedup_fastq_r2])
        else [select_first([RemoveDuplicatesSE.dedup_fastq])]

    # 步骤5: 最终质量评估
    call FastQC.RunFastQC as FinalQC {
        input:
            fastq_files = final_fastq,
            sample_name = sample_name + "_final",
            memory_gb = fastqc_memory_gb
    }

    # 步骤6: 生成综合质量报告
    call MultiQC.RunMultiQC {
        input:
            input_files = flatten([
                InitialQC.fastqc_reports,
                select_all([RunFastPPE.fastp_report, RunFastPSE.fastp_report]),
                select_all([RemoveDuplicatesPE.dedup_stats, RemoveDuplicatesSE.dedup_stats]),
                FinalQC.fastqc_reports
            ]),
            sample_name = sample_name,
            memory_gb = multiqc_memory_gb
    }

    output {
        # 最终处理后的FASTQ文件
        Array[File] processed_fastq = final_fastq
        
        # 质量控制报告
        Array[File] initial_qc_reports = InitialQC.fastqc_reports
        Array[File] final_qc_reports = FinalQC.fastqc_reports
        File multiqc_report = RunMultiQC.multiqc_report
        
        # 处理统计
        File? fastp_report_pe = RunFastPPE.fastp_report
        File? fastp_report_se = RunFastPSE.fastp_report
        File? dedup_stats_pe = RemoveDuplicatesPE.dedup_stats
        File? dedup_stats_se = RemoveDuplicatesSE.dedup_stats
        
        # 原始提取的FASTQ文件（用于比较）
        Array[File] original_fastq = ExtractReadSet.fastq_files
    }
}

# 从Read Set提取FASTQ文件的任务
task ExtractReadSet {
    input {
        String readset_arn
        String sample_name
        Boolean paired_end
    }

    command <<<
        set -e
        
        # 创建输出目录
        mkdir -p extracted_fastq
        
        # 从Read Set ARN中提取必要信息
        if [[ "~{readset_arn}" =~ ^arn:aws:omics:.*:readSet/(.*)$ ]]; then
            READSET_ID="${BASH_REMATCH[1]}"
            SEQUENCE_STORE_ID=$(echo "~{readset_arn}" | sed 's/.*sequenceStore\/\([^\/]*\).*/\1/')
        else
            # 如果直接提供ID
            READSET_ID="~{readset_arn}"
            # 需要从环境变量或配置中获取Sequence Store ID
            SEQUENCE_STORE_ID="${OMICS_SEQUENCE_STORE_ID}"
        fi
        
        echo "Extracting Read Set: $READSET_ID from Sequence Store: $SEQUENCE_STORE_ID"
        
        # 使用AWS CLI提取FASTQ文件
        if [ "~{paired_end}" = "true" ]; then
            # 双端测序
            aws omics get-read-set \
                --sequence-store-id "$SEQUENCE_STORE_ID" \
                --id "$READSET_ID" \
                --part-number 1 \
                --file SOURCE1 \
                extracted_fastq/~{sample_name}_R1.fastq.gz
            
            aws omics get-read-set \
                --sequence-store-id "$SEQUENCE_STORE_ID" \
                --id "$READSET_ID" \
                --part-number 1 \
                --file SOURCE2 \
                extracted_fastq/~{sample_name}_R2.fastq.gz
            
            echo "extracted_fastq/~{sample_name}_R1.fastq.gz" > fastq_files.txt
            echo "extracted_fastq/~{sample_name}_R2.fastq.gz" >> fastq_files.txt
        else
            # 单端测序
            aws omics get-read-set \
                --sequence-store-id "$SEQUENCE_STORE_ID" \
                --id "$READSET_ID" \
                --part-number 1 \
                --file SOURCE1 \
                extracted_fastq/~{sample_name}.fastq.gz
            
            echo "extracted_fastq/~{sample_name}.fastq.gz" > fastq_files.txt
        fi
        
        echo "Read Set extraction completed successfully"
    >>>

    output {
        Array[File] fastq_files = read_lines("fastq_files.txt")
    }

    runtime {
        docker: "public.ecr.aws/aws-cli/aws-cli:latest"
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 100 SSD"
    }

    meta {
        description: "从AWS Omics Sequence Store中提取Read Set为FASTQ文件"
    }
}
