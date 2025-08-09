version 1.0

## 修复版 Read Set 预处理工作流 v3
## 使用AWS Omics兼容的ECR镜像

workflow PreprocessingWorkflowReadSetV3 {
    meta {
        description: "修复版测序数据预处理工作流v3 - ECR镜像修复"
        version: "3.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File input_fastq
        File reference_genome
        Int min_length = 50
        Int min_quality = 20
        Int threads = 4
        Int max_length = 500
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
        String dedup_method = "simple"
        
        # 资源配置
        Int fastqc_memory_gb = 4
        Int fastp_memory_gb = 8
        Int dedup_memory_gb = 8
    }

    # 步骤1: 原始数据质量评估
    call RunFastQC as InitialQC {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name + "_initial",
            memory_gb = fastqc_memory_gb
    }

    # 步骤2: 使用FastP进行质量过滤和接头去除
    call RunFastP {
        input:
            fastq = input_fastq,
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

    # 步骤3: 去重复序列
    call RemoveDuplicates {
        input:
            fastq = RunFastP.cleaned_fastq,
            sample_name = sample_name,
            method = dedup_method,
            memory_gb = dedup_memory_gb
    }

    # 步骤4: 最终质量评估
    call RunFastQC as FinalQC {
        input:
            fastq_file = RemoveDuplicates.dedup_fastq,
            sample_name = sample_name + "_final",
            memory_gb = fastqc_memory_gb
    }

    output {
        # 最终处理后的FASTQ文件
        File processed_fastq = RemoveDuplicates.dedup_fastq
        
        # 质量控制报告
        File initial_qc_report = InitialQC.fastqc_report
        File final_qc_report = FinalQC.fastqc_report
        
        # 处理统计
        File fastp_report = RunFastP.fastp_report
        File dedup_stats = RemoveDuplicates.dedup_stats
        
        # 原始输入文件（用于比较）
        File original_fastq = input_fastq
    }
}

# FastQC质量评估任务 - 使用AWS ECR公共镜像
task RunFastQC {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
    }

    command <<<
        set -e
        
        # 更新包管理器并安装依赖
        apt-get update
        apt-get install -y openjdk-8-jre-headless wget unzip
        
        # 下载并安装FastQC
        cd /tmp
        wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
        unzip -q fastqc_v0.11.9.zip
        chmod +x FastQC/fastqc
        
        # 创建输出目录
        mkdir -p fastqc_output
        
        # 运行FastQC
        /tmp/FastQC/fastqc \
            --outdir fastqc_output \
            --threads 2 \
            --format fastq \
            ~{fastq_file}
        
        # 重命名输出文件
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "public.ecr.aws/ubuntu/ubuntu:20.04"
        memory: "${memory_gb} GB"
        cpu: 2
        disks: "local-disk 20 SSD"
    }
}

# FastP质量过滤任务 - 使用AWS ECR公共镜像
task RunFastP {
    input {
        File fastq
        String sample_name
        Int min_length
        Int min_quality
        Int max_length
        Int complexity_threshold
        Boolean enable_polyg_trimming
        Boolean enable_polyx_trimming
        Int threads
        Int memory_gb
    }

    command <<<
        set -e
        
        # 设置非交互模式
        export DEBIAN_FRONTEND=noninteractive
        
        # 安装依赖
        apt-get update
        apt-get install -y build-essential zlib1g-dev wget
        
        # 下载并编译FastP
        cd /tmp
        wget -q https://github.com/OpenGene/fastp/archive/v0.23.2.tar.gz
        tar -xzf v0.23.2.tar.gz
        cd fastp-0.23.2
        make -j~{threads}
        
        # 运行FastP
        ./fastp \
            --in1 ~{fastq} \
            --out1 ~{sample_name}_cleaned.fastq.gz \
            --json ~{sample_name}_fastp.json \
            --html ~{sample_name}_fastp.html \
            --thread ~{threads} \
            --length_required ~{min_length} \
            --qualified_quality_phred ~{min_quality} \
            --length_limit ~{max_length} \
            --complexity_threshold ~{complexity_threshold} \
            ~{if enable_polyg_trimming then "--trim_poly_g" else ""} \
            ~{if enable_polyx_trimming then "--trim_poly_x" else ""} \
            --detect_adapter_for_pe \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20
        
        # 移动输出文件到工作目录
        mv ~{sample_name}_cleaned.fastq.gz ./
        mv ~{sample_name}_fastp.html ./
        mv ~{sample_name}_fastp.json ./
    >>>

    output {
        File cleaned_fastq = "${sample_name}_cleaned.fastq.gz"
        File fastp_report = "${sample_name}_fastp.html"
        File fastp_json = "${sample_name}_fastp.json"
    }

    runtime {
        docker: "public.ecr.aws/ubuntu/ubuntu:20.04"
        memory: "${memory_gb} GB"
        cpu: "${threads}"
        disks: "local-disk 50 SSD"
    }
}

# 去重复序列任务 - 使用AWS ECR公共镜像
task RemoveDuplicates {
    input {
        File fastq
        String sample_name
        String method
        Int memory_gb
    }

    command <<<
        set -e
        
        # 安装gzip工具
        apt-get update
        apt-get install -y gzip bc
        
        # 使用简单的基于序列的去重复
        zcat ~{fastq} | \
        awk 'BEGIN{RS="@"; ORS=""} NR>1 {
            getline seq; getline plus; getline qual;
            if (!seen[seq]) {
                print "@" $0 "\n" seq "\n" plus "\n" qual "\n";
                seen[seq] = 1;
            }
        }' | gzip > ~{sample_name}_dedup.fastq.gz
        
        # 生成统计信息
        echo "Method: Simple deduplication" > ~{sample_name}_dedup_stats.txt
        echo "Input reads: $(zcat ~{fastq} | wc -l | awk '{print $1/4}')" >> ~{sample_name}_dedup_stats.txt
        echo "Output reads: $(zcat ~{sample_name}_dedup.fastq.gz | wc -l | awk '{print $1/4}')" >> ~{sample_name}_dedup_stats.txt
        
        # 计算去重复率
        INPUT_READS=$(zcat ~{fastq} | wc -l | awk '{print $1/4}')
        OUTPUT_READS=$(zcat ~{sample_name}_dedup.fastq.gz | wc -l | awk '{print $1/4}')
        DEDUP_RATE=$(echo "scale=2; (1 - $OUTPUT_READS / $INPUT_READS) * 100" | bc -l)
        echo "Deduplication rate: ${DEDUP_RATE}%" >> ~{sample_name}_dedup_stats.txt
    >>>

    output {
        File dedup_fastq = "${sample_name}_dedup.fastq.gz"
        File dedup_stats = "${sample_name}_dedup_stats.txt"
    }

    runtime {
        docker: "public.ecr.aws/ubuntu/ubuntu:20.04"
        memory: "${memory_gb} GB"
        cpu: 2
        disks: "local-disk 50 SSD"
    }
}
