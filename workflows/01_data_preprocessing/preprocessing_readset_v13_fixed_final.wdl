version 1.0

workflow CowPreprocessingFixed {
    meta {
        description: "奶牛基因组数据预处理工作流 - 包含所有修复"
        version: "13.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File input_fastq
        
        # FastQC参数
        Int fastqc_memory_gb = 8
        Int fastqc_cpu = 8
        
        # fastp参数
        Int fastp_memory_gb = 32
        Int fastp_cpu = 16
        Int min_length = 50
        Int max_length = 500
        Int min_quality = 20
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
        
        # 去重复参数
        Int dedup_memory_gb = 24
        Int dedup_cpu = 16
        String dedup_method = "exact"
    }

    # 步骤1: 原始数据质量评估
    call RunFastQC as InitialQC {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name + "_initial",
            memory_gb = fastqc_memory_gb,
            cpu = fastqc_cpu
    }

    # 步骤2: 数据清洗 (修复了fastp参数错误)
    call RunFastp {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name,
            memory_gb = fastp_memory_gb,
            cpu = fastp_cpu,
            min_length = min_length,
            max_length = max_length,
            min_quality = min_quality,
            complexity_threshold = complexity_threshold,
            enable_polyg_trimming = enable_polyg_trimming,
            enable_polyx_trimming = enable_polyx_trimming
    }

    # 步骤3: 去重复处理 (修复了bc命令依赖问题)
    call RemoveDuplicates {
        input:
            fastq_file = RunFastp.cleaned_fastq,
            sample_name = sample_name,
            memory_gb = dedup_memory_gb,
            cpu = dedup_cpu,
            method = dedup_method
    }

    output {
        # 质量报告
        File initial_qc_report = InitialQC.fastqc_report
        File initial_qc_zip = InitialQC.fastqc_zip
        
        # 处理统计
        File fastp_report = RunFastp.fastp_report
        File fastp_json = RunFastp.fastp_json
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
        Int cpu
    }

    command <<<
        set -e
        
        echo "=== FastQC质量评估开始 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{fastq_file}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "开始时间: $(date)"
        
        # 检查FastQC版本
        fastqc --version
        
        # 创建输出目录
        mkdir -p fastqc_output
        
        # 运行FastQC
        fastqc \
            --outdir fastqc_output \
            --threads ~{cpu} \
            --format fastq \
            --extract \
            --nogroup \
            ~{fastq_file}
        
        # 重命名输出文件
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
        
        echo "FastQC分析完成: ~{sample_name}"
        echo "完成时间: $(date)"
        echo "输出文件:"
        ls -lh ~{sample_name}_fastqc.*
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastqc:0.12.1"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 100 SSD"
    }
}

# fastp预处理任务 - 修复了参数错误
task RunFastp {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
        Int cpu
        Int min_length
        Int max_length
        Int min_quality
        Int complexity_threshold
        Boolean enable_polyg_trimming
        Boolean enable_polyx_trimming
    }

    command <<<
        set -e
        
        echo "=== fastp预处理开始 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{fastq_file}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "最小长度: ~{min_length}"
        echo "最大长度: ~{max_length}"
        echo "最小质量: ~{min_quality}"
        echo "复杂度阈值: ~{complexity_threshold}"
        echo "开始时间: $(date)"
        
        # 检查fastp版本
        fastp --version
        
        # 显示输入文件信息
        echo "输入文件大小:"
        ls -lh ~{fastq_file}
        
        # 运行fastp - 已修复参数错误，移除了不存在的 --detect_adapter_for_se
        fastp \
            --in1 ~{fastq_file} \
            --out1 ~{sample_name}_cleaned.fastq.gz \
            --json ~{sample_name}_fastp.json \
            --html ~{sample_name}_fastp.html \
            --thread ~{cpu} \
            --length_required ~{min_length} \
            --length_limit ~{max_length} \
            --qualified_quality_phred ~{min_quality} \
            --complexity_threshold ~{complexity_threshold} \
            ~{if enable_polyg_trimming then "--trim_poly_g" else ""} \
            ~{if enable_polyx_trimming then "--trim_poly_x" else ""} \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20 \
            --overrepresentation_analysis \
            --verbose
        
        echo "fastp预处理完成: ~{sample_name}"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "输出文件大小:"
        ls -lh ~{sample_name}_cleaned.fastq.gz
        
        # 显示处理统计
        echo "=== 处理统计 ==="
        if [ -f ~{sample_name}_fastp.json ]; then
            echo "JSON报告已生成"
        fi
        if [ -f ~{sample_name}_fastp.html ]; then
            echo "HTML报告已生成"
        fi
    >>>

    output {
        File cleaned_fastq = "${sample_name}_cleaned.fastq.gz"
        File fastp_report = "${sample_name}_fastp.html"
        File fastp_json = "${sample_name}_fastp.json"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastp:0.23.4"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 200 SSD"
    }
}

# 去重复序列任务 - 修复了bc命令依赖问题
task RemoveDuplicates {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
        Int cpu
        String method
    }

    command <<<
        set -e
        
        echo "=== 去重复处理开始 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{fastq_file}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "去重复方法: ~{method}"
        echo "开始时间: $(date)"
        
        # 检查seqkit版本
        seqkit version
        
        # 显示输入文件信息
        echo "输入文件大小:"
        ls -lh ~{fastq_file}
        
        # 统计原始序列数量
        echo "=== 原始序列统计 ==="
        seqkit stats ~{fastq_file}
        
        # 使用seqkit进行去重复
        seqkit rmdup \
            --by-seq \
            --ignore-case \
            ~{fastq_file} \
            --out-file ~{sample_name}_dedup.fastq.gz \
            --threads ~{cpu}
        
        # 统计去重复后的序列数量
        echo "=== 去重复后序列统计 ==="
        seqkit stats ~{sample_name}_dedup.fastq.gz
        
        # 生成详细统计报告 - 使用awk替代bc进行数学计算
        echo "=== 去重复统计报告 ===" > ~{sample_name}_dedup_stats.txt
        echo "样本名称: ~{sample_name}" >> ~{sample_name}_dedup_stats.txt
        echo "去重复方法: ~{method}" >> ~{sample_name}_dedup_stats.txt
        echo "处理时间: $(date)" >> ~{sample_name}_dedup_stats.txt
        echo "CPU核心数: ~{cpu}" >> ~{sample_name}_dedup_stats.txt
        echo "内存大小: ~{memory_gb}GB" >> ~{sample_name}_dedup_stats.txt
        echo "" >> ~{sample_name}_dedup_stats.txt
        
        echo "=== 原始序列统计 ===" >> ~{sample_name}_dedup_stats.txt
        seqkit stats ~{fastq_file} >> ~{sample_name}_dedup_stats.txt
        echo "" >> ~{sample_name}_dedup_stats.txt
        
        echo "=== 去重复后序列统计 ===" >> ~{sample_name}_dedup_stats.txt
        seqkit stats ~{sample_name}_dedup.fastq.gz >> ~{sample_name}_dedup_stats.txt
        echo "" >> ~{sample_name}_dedup_stats.txt
        
        # 使用awk计算去重复率，避免bc依赖
        ORIGINAL_COUNT=$(seqkit stats ~{fastq_file} | tail -n 1 | awk '{gsub(/,/, "", $4); print $4}')
        DEDUP_COUNT=$(seqkit stats ~{sample_name}_dedup.fastq.gz | tail -n 1 | awk '{gsub(/,/, "", $4); print $4}')
        
        if [ "$ORIGINAL_COUNT" -gt 0 ]; then
            DUPLICATE_COUNT=$((ORIGINAL_COUNT - DEDUP_COUNT))
            DUPLICATE_RATE=$(awk "BEGIN {printf \"%.2f\", ($DUPLICATE_COUNT * 100 / $ORIGINAL_COUNT)}")
            echo "重复序列数量: $DUPLICATE_COUNT" >> ~{sample_name}_dedup_stats.txt
            echo "重复率: ${DUPLICATE_RATE}%" >> ~{sample_name}_dedup_stats.txt
        fi
        
        echo "去重复处理完成: ~{sample_name}"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "输出文件大小:"
        ls -lh ~{sample_name}_dedup.fastq.gz
        
        # 显示最终统计
        echo "=== 最终统计摘要 ==="
        echo "原始序列数量: $ORIGINAL_COUNT"
        echo "去重复后数量: $DEDUP_COUNT"
        if [ "$ORIGINAL_COUNT" -gt 0 ]; then
            echo "重复序列数量: $DUPLICATE_COUNT"
            echo "重复率: ${DUPLICATE_RATE}%"
        fi
    >>>

    output {
        File dedup_fastq = "${sample_name}_dedup.fastq.gz"
        File dedup_stats = "${sample_name}_dedup_stats.txt"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/seqkit:2.5.1"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 200 SSD"
    }
}