version 1.0

workflow SequenceAlignmentWorkflowV2Fixed {
    meta {
        description: "奶牛基因组序列比对工作流v2 - 修复版本，使用包含BWA和samtools的镜像"
        version: "2.1"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File cleaned_fastq
        File reference_genome
        
        # BWA比对参数
        Int bwa_cpu = 16
        Int bwa_memory_gb = 32
        Int min_seed_length = 19
        Int bandwidth = 100
        
        # SAMtools处理参数
        Int samtools_cpu = 8
        Int samtools_memory_gb = 16
        Int min_mapq = 20
        
        # 索引构建参数
        Int index_cpu = 8
        Int index_memory_gb = 16
    }

    # 步骤1: 构建参考基因组索引
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
        bwa 2>&1 | head -3
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
        echo "输入FASTQ: ~{fastq_file}"
        echo "参考基因组: ~{reference_genome}"
        echo "样本名称: ~{sample_name}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "开始时间: $(date)"
        
        # 复制索引文件到工作目录
        echo "复制索引文件..."
        for file in ~{sep=' ' reference_index}; do
            cp "$file" ./
        done
        
        # 检查索引文件
        echo "=== 索引文件检查 ==="
        ls -la reference.*
        
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
        echo "输入SAM文件: ~{sam_file}"
        echo "样本名称: ~{sample_name}"
        echo "最小MAPQ: ~{min_mapq}"
        echo "开始时间: $(date)"
        
        # 转换SAM到BAM并排序
        echo "转换SAM到BAM并排序..."
        samtools view -@ ~{cpu} -bS ~{sam_file} | \
        samtools sort -@ ~{cpu} -o ~{sample_name}_sorted.bam -
        
        # 创建BAM索引
        echo "创建BAM索引..."
        samtools index ~{sample_name}_sorted.bam
        
        # 过滤低质量比对
        echo "过滤低质量比对 (MAPQ >= ~{min_mapq})..."
        samtools view -@ ~{cpu} -b -q ~{min_mapq} ~{sample_name}_sorted.bam > ~{sample_name}_filtered.bam
        samtools index ~{sample_name}_filtered.bam
        
        # 生成比对统计
        echo "生成比对统计..."
        samtools stats ~{sample_name}_sorted.bam > ~{sample_name}_alignment_stats.txt
        samtools flagstat ~{sample_name}_sorted.bam > ~{sample_name}_flagstat.txt
        
        # 显示基本统计信息
        echo "=== 比对统计摘要 ==="
        echo "排序后BAM文件大小:"
        ls -lh ~{sample_name}_sorted.bam
        echo "过滤后BAM文件大小:"
        ls -lh ~{sample_name}_filtered.bam
        
        echo "Flagstat报告:"
        cat ~{sample_name}_flagstat.txt
        
        echo "处理完成"
        echo "完成时间: $(date)"
    >>>

    output {
        File sorted_bam = "${sample_name}_sorted.bam"
        File sorted_bam_index = "${sample_name}_sorted.bam.bai"
        File filtered_bam = "${sample_name}_filtered.bam"
        File filtered_bam_index = "${sample_name}_filtered.bam.bai"
        File alignment_stats = "${sample_name}_alignment_stats.txt"
        File flagstat_report = "${sample_name}_flagstat.txt"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/bioinformatics:bwa-samtools-1.17"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 300 SSD"
    }
}
