version 1.0

## 简化版 FQ2BAM 工作流
## 模拟NVIDIA Parabricks FQ2BAM功能

workflow FQ2BAMWorkflow {
    meta {
        description: "简化版FQ2BAM工作流，包含预处理和比对"
        version: "1.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File input_fastq
        File reference_genome
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        
        # 资源配置
        Int fastqc_memory_gb = 8
        Int fastp_memory_gb = 16
        Int bwa_memory_gb = 32
        Int samtools_memory_gb = 16
    }

    # 步骤1: 原始数据质量评估
    call RunFastQC as InitialQC {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name + "_initial",
            memory_gb = fastqc_memory_gb
    }

    # 步骤2: 质量过滤和预处理
    call RunFastP {
        input:
            fastq = input_fastq,
            sample_name = sample_name,
            min_length = min_length,
            min_quality = min_quality,
            threads = threads,
            memory_gb = fastp_memory_gb
    }

    # 步骤3: BWA比对
    call BWAAlignment {
        input:
            fastq = RunFastP.cleaned_fastq,
            reference_genome = reference_genome,
            sample_name = sample_name,
            threads = threads,
            memory_gb = bwa_memory_gb
    }

    # 步骤4: SAM到BAM转换和排序
    call SamtoolsProcessing {
        input:
            sam_file = BWAAlignment.sam_file,
            sample_name = sample_name,
            threads = threads,
            memory_gb = samtools_memory_gb
    }

    # 步骤5: 最终质量评估
    call RunFastQC as FinalQC {
        input:
            fastq_file = RunFastP.cleaned_fastq,
            sample_name = sample_name + "_final",
            memory_gb = fastqc_memory_gb
    }

    output {
        # 主要输出：BAM文件
        File output_bam = SamtoolsProcessing.sorted_bam
        File output_bam_index = SamtoolsProcessing.bam_index
        
        # 预处理后的FASTQ
        File processed_fastq = RunFastP.cleaned_fastq
        
        # 质量控制报告
        File initial_qc_report = InitialQC.fastqc_report
        File final_qc_report = FinalQC.fastqc_report
        File fastp_report = RunFastP.fastp_report
        
        # 比对统计
        File alignment_stats = SamtoolsProcessing.alignment_stats
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
        
        # 安装FastQC
        yum update -y
        yum install -y java-1.8.0-openjdk wget unzip
        
        cd /tmp
        wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
        unzip fastqc_v0.11.9.zip
        chmod +x FastQC/fastqc
        
        mkdir -p fastqc_output
        
        /tmp/FastQC/fastqc \
            --outdir fastqc_output \
            --threads 2 \
            --format fastq \
            ~{fastq_file}
        
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "amazonlinux:2"
        memory: "${memory_gb} GB"
        cpu: 2
        disks: "local-disk 50 SSD"
    }
}

# FastP质量过滤任务
task RunFastP {
    input {
        File fastq
        String sample_name
        Int min_length
        Int min_quality
        Int threads
        Int memory_gb
    }

    command <<<
        set -e
        
        # 安装依赖
        yum update -y
        yum install -y gcc-c++ make zlib-devel wget tar
        
        # 下载并编译FastP
        cd /tmp
        wget https://github.com/OpenGene/fastp/archive/v0.23.2.tar.gz
        tar -xzf v0.23.2.tar.gz
        cd fastp-0.23.2
        make
        
        # 运行FastP
        ./fastp \
            --in1 ~{fastq} \
            --out1 ~{sample_name}_cleaned.fastq.gz \
            --json ~{sample_name}_fastp.json \
            --html ~{sample_name}_fastp.html \
            --thread ~{threads} \
            --length_required ~{min_length} \
            --qualified_quality_phred ~{min_quality} \
            --detect_adapter_for_pe \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20
        
        mv ~{sample_name}_cleaned.fastq.gz /tmp/
        mv ~{sample_name}_fastp.html /tmp/
        mv ~{sample_name}_fastp.json /tmp/
    >>>

    output {
        File cleaned_fastq = "${sample_name}_cleaned.fastq.gz"
        File fastp_report = "${sample_name}_fastp.html"
        File fastp_json = "${sample_name}_fastp.json"
    }

    runtime {
        docker: "amazonlinux:2"
        memory: "${memory_gb} GB"
        cpu: "${threads}"
        disks: "local-disk 100 SSD"
    }
}

# BWA比对任务
task BWAAlignment {
    input {
        File fastq
        File reference_genome
        String sample_name
        Int threads
        Int memory_gb
    }

    command <<<
        set -e
        
        # 安装BWA和samtools
        yum update -y
        yum install -y gcc make zlib-devel bzip2-devel xz-devel ncurses-devel wget tar
        
        # 安装BWA
        cd /tmp
        wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
        tar -xjf bwa-0.7.17.tar.bz2
        cd bwa-0.7.17
        make
        
        # 创建BWA索引
        ./bwa index ~{reference_genome}
        
        # 执行比对
        ./bwa mem -t ~{threads} ~{reference_genome} ~{fastq} > ~{sample_name}.sam
        
        mv ~{sample_name}.sam /tmp/
    >>>

    output {
        File sam_file = "${sample_name}.sam"
    }

    runtime {
        docker: "amazonlinux:2"
        memory: "${memory_gb} GB"
        cpu: "${threads}"
        disks: "local-disk 200 SSD"
    }
}

# Samtools处理任务
task SamtoolsProcessing {
    input {
        File sam_file
        String sample_name
        Int threads
        Int memory_gb
    }

    command <<<
        set -e
        
        # 安装samtools
        yum update -y
        yum install -y gcc make zlib-devel bzip2-devel xz-devel ncurses-devel wget tar
        
        cd /tmp
        wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2
        tar -xjf samtools-1.17.tar.bz2
        cd samtools-1.17
        make
        
        # SAM到BAM转换
        ./samtools view -@ ~{threads} -bS ~{sam_file} > ~{sample_name}.bam
        
        # BAM排序
        ./samtools sort -@ ~{threads} ~{sample_name}.bam -o ~{sample_name}_sorted.bam
        
        # 创建索引
        ./samtools index ~{sample_name}_sorted.bam
        
        # 生成统计信息
        ./samtools flagstat ~{sample_name}_sorted.bam > ~{sample_name}_alignment_stats.txt
        
        mv ~{sample_name}_sorted.bam /tmp/
        mv ~{sample_name}_sorted.bam.bai /tmp/
        mv ~{sample_name}_alignment_stats.txt /tmp/
    >>>

    output {
        File sorted_bam = "${sample_name}_sorted.bam"
        File bam_index = "${sample_name}_sorted.bam.bai"
        File alignment_stats = "${sample_name}_alignment_stats.txt"
    }

    runtime {
        docker: "amazonlinux:2"
        memory: "${memory_gb} GB"
        cpu: "${threads}"
        disks: "local-disk 200 SSD"
    }
}
