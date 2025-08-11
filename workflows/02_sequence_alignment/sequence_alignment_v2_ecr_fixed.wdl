version 1.0

workflow SequenceAlignmentWorkflowV2 {
    meta {
        description: "序列比对工作流v2 - 使用私有ECR镜像，包含BWA比对和SAMtools处理"
        version: "2.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File cleaned_fastq  # 来自预处理步骤的清洁FASTQ文件
        File reference_genome  # 参考基因组FASTA文件
        
        # BWA比对参数
        Int bwa_cpu = 16
        Int bwa_memory_gb = 32
        Int min_seed_length = 19
        Int bandwidth = 100
        
        # SAMtools参数
        Int samtools_cpu = 8
        Int samtools_memory_gb = 16
        Int min_mapq = 20
        
        # 索引构建参数
        Int index_cpu = 8
        Int index_memory_gb = 16
    }

    # 步骤1: 构建参考基因组索引 (使用SAMtools镜像，包含BWA)
    call BuildReferenceIndex {
        input:
            reference_genome = reference_genome,
            cpu = index_cpu,
            memory_gb = index_memory_gb
    }

    # 步骤2: BWA比对
    call BWAAlignment {
        input:
            fastq_file = cleaned_fastq,
            reference_genome = reference_genome,
            reference_index = BuildReferenceIndex.bwa_index,
            sample_name = sample_name,
            cpu = bwa_cpu,
            memory_gb = bwa_memory_gb,
            min_seed_length = min_seed_length,
            bandwidth = bandwidth
    }

    # 步骤3: SAM/BAM处理和基本质量评估
    call ProcessAlignment {
        input:
            sam_file = BWAAlignment.sam_file,
            sample_name = sample_name,
            cpu = samtools_cpu,
            memory_gb = samtools_memory_gb,
            min_mapq = min_mapq
    }

    output {
        # 最终比对结果
        File final_bam = ProcessAlignment.sorted_bam
        File final_bam_index = ProcessAlignment.sorted_bam_index
        
        # 质量评估报告
        File alignment_stats = ProcessAlignment.alignment_stats
        File flagstat_report = ProcessAlignment.flagstat_report
        
        # 索引文件（可重用）
        Array[File] reference_index_files = BuildReferenceIndex.bwa_index
    }
}

# 构建参考基因组索引
task BuildReferenceIndex {
    input {
        File reference_genome
        Int cpu
        Int memory_gb
    }

    command <<<
        set -e
        
        echo "=== 构建参考基因组索引 ==="
        echo "参考基因组: ~{reference_genome}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "开始时间: $(date)"
        
        # 复制参考基因组到工作目录
        cp ~{reference_genome} reference.fasta
        
        # 检查工具版本
        echo "=== 工具版本信息 ==="
        bwa 2>&1 | head -3 || echo "BWA not found in this image"
        samtools --version | head -2
        
        # 构建BWA索引
        echo "构建BWA索引..."
        bwa index reference.fasta
        
        # 构建SAMtools索引
        echo "构建SAMtools索引..."
        samtools faidx reference.fasta
        
        echo "索引构建完成"
        echo "完成时间: $(date)"
        
        # 列出所有索引文件
        echo "=== 生成的索引文件 ==="
        ls -la reference.*
    >>>

    output {
        Array[File] bwa_index = glob("reference.*")
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bioinformatics:bwa-samtools-1.17"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 200 SSD"
    }
}

# BWA序列比对
task BWAAlignment {
    input {
        File fastq_file
        File reference_genome
        Array[File] reference_index
        String sample_name
        Int cpu
        Int memory_gb
        Int min_seed_length
        Int bandwidth
    }

    command <<<
        set -e
        
        echo "=== BWA序列比对 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{fastq_file}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "最小种子长度: ~{min_seed_length}"
        echo "带宽: ~{bandwidth}"
        echo "开始时间: $(date)"
        
        # 复制参考基因组和索引文件
        cp ~{reference_genome} reference.fasta
        
        # 复制所有索引文件
        for index_file in ~{sep=' ' reference_index}; do
            cp "$index_file" ./
        done
        
        # 检查BWA版本
        echo "=== BWA版本信息 ==="
        bwa 2>&1 | head -3
        
        # 显示输入文件信息
        echo "=== 输入文件信息 ==="
        ls -lh ~{fastq_file}
        
        # 运行BWA比对
        echo "开始BWA-MEM比对..."
        bwa mem \
            -t ~{cpu} \
            -k ~{min_seed_length} \
            -w ~{bandwidth} \
            -r 1.5 \
            -A 1 \
            -B 4 \
            -O 6 \
            -E 1 \
            -L 5 \
            -R "@RG\tID:~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA\tLB:~{sample_name}\tPU:~{sample_name}" \
            reference.fasta \
            ~{fastq_file} \
            > ~{sample_name}.sam
        
        echo "BWA比对完成"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "=== 输出文件信息 ==="
        ls -lh ~{sample_name}.sam
        
        # 显示比对统计预览
        echo "=== 比对结果预览 ==="
        echo "总行数: $(wc -l < ~{sample_name}.sam)"
        echo "头部信息行数: $(grep -c '^@' ~{sample_name}.sam)"
        echo "比对记录数: $(grep -c -v '^@' ~{sample_name}.sam)"
    >>>

    output {
        File sam_file = "${sample_name}.sam"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bioinformatics:bwa-samtools-1.17"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 400 SSD"
    }
}

# SAM/BAM处理和质量评估
task ProcessAlignment {
    input {
        File sam_file
        String sample_name
        Int cpu
        Int memory_gb
        Int min_mapq
    }

    command <<<
        set -e
        
        echo "=== SAM/BAM处理和质量评估 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{sam_file}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "最小MAPQ: ~{min_mapq}"
        echo "开始时间: $(date)"
        
        # 检查samtools版本
        echo "=== SAMtools版本信息 ==="
        samtools --version
        
        # SAM转BAM并过滤低质量比对
        echo "SAM转BAM并过滤低质量比对..."
        samtools view \
            -@ ~{cpu} \
            -b \
            -q ~{min_mapq} \
            -F 4 \
            ~{sam_file} \
            | samtools sort \
                -@ ~{cpu} \
                -o ~{sample_name}.sorted.bam
        
        # 创建BAM索引
        echo "创建BAM索引..."
        samtools index ~{sample_name}.sorted.bam
        
        # 生成详细统计信息
        echo "生成比对统计信息..."
        samtools stats ~{sample_name}.sorted.bam > ~{sample_name}_alignment_stats.txt
        
        # 生成flagstat报告
        echo "生成flagstat报告..."
        samtools flagstat ~{sample_name}.sorted.bam > ~{sample_name}_flagstat.txt
        
        # 生成idxstats报告
        echo "生成idxstats报告..."
        samtools idxstats ~{sample_name}.sorted.bam > ~{sample_name}_idxstats.txt
        
        echo "SAM/BAM处理完成"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "=== 输出文件信息 ==="
        ls -lh ~{sample_name}.sorted.bam*
        ls -lh ~{sample_name}_*.txt
        
        # 显示关键统计信息
        echo "=== 关键统计信息 ==="
        echo "--- Flagstat报告 ---"
        cat ~{sample_name}_flagstat.txt
        
        echo ""
        echo "--- 比对质量摘要 ---"
        grep -E "^SN.*reads mapped|^SN.*reads unmapped|^SN.*reads properly paired|^SN.*average length|^SN.*average quality" ~{sample_name}_alignment_stats.txt || true
        
        echo ""
        echo "--- 染色体比对统计 (前10个) ---"
        head -10 ~{sample_name}_idxstats.txt
    >>>

    output {
        File sorted_bam = "${sample_name}.sorted.bam"
        File sorted_bam_index = "${sample_name}.sorted.bam.bai"
        File alignment_stats = "${sample_name}_alignment_stats.txt"
        File flagstat_report = "${sample_name}_flagstat.txt"
        File idxstats_report = "${sample_name}_idxstats.txt"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/samtools:1.17"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 500 SSD"
    }
}
