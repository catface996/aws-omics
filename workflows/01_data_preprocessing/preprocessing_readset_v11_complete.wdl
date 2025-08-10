version 1.0

workflow CompletePreprocessingWorkflowV11 {
    meta {
        description: "完整的测序数据预处理工作流v11 - 包含质量评估、接头去除、质量过滤、长度过滤和去重复"
        version: "11.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File input_fastq
        
        # FastQC参数
        Int fastqc_memory_gb = 8
        
        # fastp参数
        Int fastp_memory_gb = 16
        Int min_length = 50
        Int max_length = 500
        Int min_quality = 20
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
        
        # 去重复参数
        Int dedup_memory_gb = 16
        String dedup_method = "exact"  # exact, prefix, suffix
    }

    # 步骤1: 原始数据质量评估
    call RunFastQC as InitialQC {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name + "_initial",
            memory_gb = fastqc_memory_gb
    }

    # 步骤2: 接头去除 + 质量过滤 + 长度过滤
    call RunFastp {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name,
            memory_gb = fastp_memory_gb,
            min_length = min_length,
            max_length = max_length,
            min_quality = min_quality,
            complexity_threshold = complexity_threshold,
            enable_polyg_trimming = enable_polyg_trimming,
            enable_polyx_trimming = enable_polyx_trimming
    }

    # 步骤3: 去重复序列
    call RemoveDuplicates {
        input:
            fastq_file = RunFastp.cleaned_fastq,
            sample_name = sample_name,
            memory_gb = dedup_memory_gb,
            method = dedup_method
    }

    # 步骤4: 最终质量评估
    call RunFastQC as FinalQC {
        input:
            fastq_file = RemoveDuplicates.dedup_fastq,
            sample_name = sample_name + "_final",
            memory_gb = fastqc_memory_gb
    }

    output {
        # 质量报告
        File initial_qc_report = InitialQC.fastqc_report
        File initial_qc_zip = InitialQC.fastqc_zip
        File final_qc_report = FinalQC.fastqc_report
        File final_qc_zip = FinalQC.fastqc_zip
        
        # 处理统计
        File fastp_report = RunFastp.fastp_report
        File dedup_stats = RemoveDuplicates.dedup_stats
        
        # 最终清洁数据
        File cleaned_fastq = RemoveDuplicates.dedup_fastq
        
        # 原始数据
        File original_fastq = input_fastq
    }
}

# FastQC质量评估任务
task RunFastQC {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
    }

    command <<<
        set -e
        
        echo "开始FastQC质量评估: ~{sample_name}"
        
        # 创建输出目录
        mkdir -p fastqc_output
        
        # 运行FastQC
        fastqc \
            --outdir fastqc_output \
            --threads 8 \
            --format fastq \
            --extract \
            ~{fastq_file}
        
        # 重命名输出文件
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
        
        echo "FastQC分析完成: ~{sample_name}"
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastqc:0.12.1"
        memory: "${memory_gb} GB"
        cpu: 8
        disks: "local-disk 50 SSD"
    }
}

# fastp预处理任务 (接头去除 + 质量过滤 + 长度过滤)
task RunFastp {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
        Int min_length
        Int max_length
        Int min_quality
        Int complexity_threshold
        Boolean enable_polyg_trimming
        Boolean enable_polyx_trimming
    }

    command <<<
        set -e
        
        echo "开始fastp预处理: ~{sample_name}"
        
        # 运行fastp
        fastp \
            --in1 ~{fastq_file} \
            --out1 ~{sample_name}_cleaned.fastq.gz \
            --json ~{sample_name}_fastp.json \
            --html ~{sample_name}_fastp.html \
            --thread 16 \
            --length_required ~{min_length} \
            --length_limit ~{max_length} \
            --qualified_quality_phred ~{min_quality} \
            --complexity_threshold ~{complexity_threshold} \
            ~{if enable_polyg_trimming then "--trim_poly_g" else ""} \
            ~{if enable_polyx_trimming then "--trim_poly_x" else ""} \
            --detect_adapter_for_se \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --verbose
        
        echo "fastp预处理完成: ~{sample_name}"
        
        # 显示统计信息
        echo "输出文件大小:"
        ls -lh ~{sample_name}_cleaned.fastq.gz
    >>>

    output {
        File cleaned_fastq = "${sample_name}_cleaned.fastq.gz"
        File fastp_report = "${sample_name}_fastp.html"
        File fastp_json = "${sample_name}_fastp.json"
    }

    runtime {
        docker: "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
        memory: "${memory_gb} GB"
        cpu: 16
        disks: "local-disk 100 SSD"
    }
}

# 去重复序列任务
task RemoveDuplicates {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
        String method
    }

    command <<<
        set -e
        
        echo "开始去重复处理: ~{sample_name}"
        echo "去重复方法: ~{method}"
        
        # 使用seqkit进行去重复
        seqkit rmdup \
            --by-seq \
            --ignore-case \
            ~{fastq_file} \
            --out-file ~{sample_name}_dedup.fastq.gz \
            --threads 16
        
        # 生成统计报告
        echo "=== 去重复统计报告 ===" > ~{sample_name}_dedup_stats.txt
        echo "样本名称: ~{sample_name}" >> ~{sample_name}_dedup_stats.txt
        echo "去重复方法: ~{method}" >> ~{sample_name}_dedup_stats.txt
        echo "处理时间: $(date)" >> ~{sample_name}_dedup_stats.txt
        echo "" >> ~{sample_name}_dedup_stats.txt
        
        # 统计原始和去重复后的序列数量
        echo "原始序列数量:" >> ~{sample_name}_dedup_stats.txt
        seqkit stats ~{fastq_file} >> ~{sample_name}_dedup_stats.txt
        echo "" >> ~{sample_name}_dedup_stats.txt
        
        echo "去重复后序列数量:" >> ~{sample_name}_dedup_stats.txt
        seqkit stats ~{sample_name}_dedup.fastq.gz >> ~{sample_name}_dedup_stats.txt
        
        echo "去重复处理完成: ~{sample_name}"
    >>>

    output {
        File dedup_fastq = "${sample_name}_dedup.fastq.gz"
        File dedup_stats = "${sample_name}_dedup_stats.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/seqkit:2.5.1--h9ee0642_0"
        memory: "${memory_gb} GB"
        cpu: 16
        disks: "local-disk 100 SSD"
    }
}
