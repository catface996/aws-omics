version 1.0

## 奶牛基因组数据预处理工作流
## 包含质量评估、接头去除、质量过滤、长度过滤和去重复序列
##
## 作者: AWS Omics 奶牛基因组学演示项目
## 版本: 1.0
## 日期: 2025-08-08

import "tasks/fastqc_task.wdl" as FastQC
import "tasks/trimmomatic_task.wdl" as Trimmomatic
import "tasks/fastp_task.wdl" as FastP
import "tasks/deduplication_task.wdl" as Dedup
import "tasks/multiqc_task.wdl" as MultiQC

workflow PreprocessingWorkflow {
    meta {
        description: "完整的测序数据预处理工作流，包含质量控制、接头去除、质量过滤、长度过滤和去重复"
        version: "1.0"
        author: "AWS Omics Demo Project"
    }

    parameter_meta {
        sample_name: "样本名称"
        fastq_r1: "正向测序文件 (R1)"
        fastq_r2: "反向测序文件 (R2，可选)"
        reference_genome: "参考基因组文件"
        adapter_fasta: "接头序列文件"
        min_length: "最小读长阈值"
        min_quality: "最小质量分数"
        threads: "并行线程数"
    }

    input {
        String sample_name
        File fastq_r1
        File? fastq_r2
        File? reference_genome
        File? adapter_fasta
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Boolean paired_end = defined(fastq_r2)
    }

    # 步骤1: 原始数据质量评估
    call FastQC.RunFastQC as InitialQC_R1 {
        input:
            fastq_file = fastq_r1,
            sample_name = sample_name + "_R1_initial",
            threads = threads
    }

    if (paired_end) {
        call FastQC.RunFastQC as InitialQC_R2 {
            input:
                fastq_file = select_first([fastq_r2]),
                sample_name = sample_name + "_R2_initial",
                threads = threads
        }
    }

    # 步骤2: 接头去除和质量过滤 (使用Trimmomatic)
    if (paired_end) {
        call Trimmomatic.RunTrimmomaticPE as TrimPE {
            input:
                fastq_r1 = fastq_r1,
                fastq_r2 = select_first([fastq_r2]),
                sample_name = sample_name,
                adapter_fasta = adapter_fasta,
                min_length = min_length,
                min_quality = min_quality,
                threads = threads
        }
    }

    if (!paired_end) {
        call Trimmomatic.RunTrimmomaticSE as TrimSE {
            input:
                fastq = fastq_r1,
                sample_name = sample_name,
                adapter_fasta = adapter_fasta,
                min_length = min_length,
                min_quality = min_quality,
                threads = threads
        }
    }

    # 步骤3: 高级质量过滤和长度过滤 (使用fastp)
    if (paired_end) {
        call FastP.RunFastPPE as FastPPE {
            input:
                fastq_r1 = select_first([TrimPE.trimmed_r1]),
                fastq_r2 = select_first([TrimPE.trimmed_r2]),
                sample_name = sample_name,
                min_length = min_length,
                min_quality = min_quality,
                threads = threads
        }
    }

    if (!paired_end) {
        call FastP.RunFastPSE as FastPSE {
            input:
                fastq = select_first([TrimSE.trimmed]),
                sample_name = sample_name,
                min_length = min_length,
                min_quality = min_quality,
                threads = threads
        }
    }

    # 步骤4: 去重复序列
    if (paired_end) {
        call Dedup.RemoveDuplicatesPE as DedupPE {
            input:
                fastq_r1 = select_first([FastPPE.filtered_r1]),
                fastq_r2 = select_first([FastPPE.filtered_r2]),
                sample_name = sample_name,
                threads = threads
        }
    }

    if (!paired_end) {
        call Dedup.RemoveDuplicatesSE as DedupSE {
            input:
                fastq = select_first([FastPSE.filtered]),
                sample_name = sample_name,
                threads = threads
        }
    }

    # 步骤5: 最终质量评估
    File final_r1 = select_first([DedupPE.dedup_r1, DedupSE.dedup])
    call FastQC.RunFastQC as FinalQC_R1 {
        input:
            fastq_file = final_r1,
            sample_name = sample_name + "_R1_final",
            threads = threads
    }

    if (paired_end) {
        call FastQC.RunFastQC as FinalQC_R2 {
            input:
                fastq_file = select_first([DedupPE.dedup_r2]),
                sample_name = sample_name + "_R2_final",
                threads = threads
        }
    }

    # 步骤6: 生成综合质量报告
    Array[File] all_qc_reports = flatten([
        [InitialQC_R1.fastqc_html, InitialQC_R1.fastqc_zip],
        select_all([InitialQC_R2.fastqc_html, InitialQC_R2.fastqc_zip]),
        [FinalQC_R1.fastqc_html, FinalQC_R1.fastqc_zip],
        select_all([FinalQC_R2.fastqc_html, FinalQC_R2.fastqc_zip]),
        select_all([FastPPE.fastp_html, FastPPE.fastp_json, FastPSE.fastp_html, FastPSE.fastp_json])
    ])

    call MultiQC.RunMultiQC as GenerateReport {
        input:
            input_files = all_qc_reports,
            sample_name = sample_name,
            output_name = sample_name + "_preprocessing_report"
    }

    # 输出结果
    output {
        # 最终处理后的FASTQ文件
        File processed_r1 = final_r1
        File? processed_r2 = DedupPE.dedup_r2
        
        # 质量控制报告
        File initial_qc_r1_html = InitialQC_R1.fastqc_html
        File initial_qc_r1_zip = InitialQC_R1.fastqc_zip
        File? initial_qc_r2_html = InitialQC_R2.fastqc_html
        File? initial_qc_r2_zip = InitialQC_R2.fastqc_zip
        
        File final_qc_r1_html = FinalQC_R1.fastqc_html
        File final_qc_r1_zip = FinalQC_R1.fastqc_zip
        File? final_qc_r2_html = FinalQC_R2.fastqc_html
        File? final_qc_r2_zip = FinalQC_R2.fastqc_zip
        
        # 处理统计报告
        File? trimmomatic_log = if paired_end then TrimPE.log_file else TrimSE.log_file
        File? fastp_html = if paired_end then FastPPE.fastp_html else FastPSE.fastp_html
        File? fastp_json = if paired_end then FastPPE.fastp_json else FastPSE.fastp_json
        File? dedup_log = if paired_end then DedupPE.log_file else DedupSE.log_file
        
        # 综合报告
        File multiqc_html = GenerateReport.multiqc_html
        File multiqc_data = GenerateReport.multiqc_data
        
        # 处理统计信息
        String preprocessing_summary = "样本 ${sample_name} 预处理完成"
    }
}
