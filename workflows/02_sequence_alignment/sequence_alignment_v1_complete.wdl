version 1.0

workflow SequenceAlignmentWorkflowV1 {
    meta {
        description: "完整的序列比对工作流v1 - 包含索引构建、比对、后处理和质量评估"
        version: "1.0"
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
        
        # Picard参数
        Int picard_memory_gb = 24
        
        # Qualimap参数
        Int qualimap_memory_gb = 16
    }

    # 步骤1: 构建参考基因组索引
    call BuildReferenceIndex {
        input:
            reference_genome = reference_genome,
            cpu = 8,
            memory_gb = 16
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

    # 步骤3: SAM/BAM处理
    call ProcessAlignment {
        input:
            sam_file = BWAAlignment.sam_file,
            sample_name = sample_name,
            cpu = samtools_cpu,
            memory_gb = samtools_memory_gb,
            min_mapq = min_mapq
    }

    # 步骤4: 重复标记
    call MarkDuplicates {
        input:
            bam_file = ProcessAlignment.sorted_bam,
            sample_name = sample_name,
            memory_gb = picard_memory_gb
    }

    # 步骤5: 比对质量评估
    call AlignmentQC {
        input:
            bam_file = MarkDuplicates.marked_bam,
            bam_index = MarkDuplicates.marked_bam_index,
            sample_name = sample_name,
            memory_gb = qualimap_memory_gb
    }

    output {
        # 最终比对结果
        File final_bam = MarkDuplicates.marked_bam
        File final_bam_index = MarkDuplicates.marked_bam_index
        
        # 质量评估报告
        File alignment_stats = AlignmentQC.alignment_stats
        File qualimap_report = AlignmentQC.qualimap_report
        
        # 处理统计
        File duplicate_metrics = MarkDuplicates.duplicate_metrics
        File processing_stats = ProcessAlignment.processing_stats
        
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
        
        # 构建BWA索引
        echo "构建BWA索引..."
        bwa index reference.fasta
        
        # 构建SAMtools索引
        echo "构建SAMtools索引..."
        samtools faidx reference.fasta
        
        # 构建序列字典
        echo "构建序列字典..."
        picard CreateSequenceDictionary \
            R=reference.fasta \
            O=reference.dict
        
        echo "索引构建完成"
        echo "完成时间: $(date)"
        
        # 列出所有索引文件
        ls -la reference.*
    >>>

    output {
        Array[File] bwa_index = glob("reference.*")
    }

    runtime {
        docker: "quay.io/biocontainers/mulled-v2-002f51ea92721407ef440b921fb5940f424be842:43ec6124f9f4f875515f9548733b8b4e5fe9f506-0"
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
        echo "开始时间: $(date)"
        
        # 复制参考基因组和索引文件
        cp ~{reference_genome} reference.fasta
        
        # 复制所有索引文件
        for index_file in ~{sep=' ' reference_index}; do
            cp "$index_file" ./
        done
        
        # 检查BWA版本
        bwa 2>&1 | head -3
        
        # 运行BWA比对
        echo "开始BWA比对..."
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
            -R "@RG\tID:~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA\tLB:~{sample_name}" \
            reference.fasta \
            ~{fastq_file} \
            > ~{sample_name}.sam
        
        echo "BWA比对完成"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "输出文件大小:"
        ls -lh ~{sample_name}.sam
        
        # 显示比对统计
        echo "比对reads数量:"
        wc -l ~{sample_name}.sam
    >>>

    output {
        File sam_file = "${sample_name}.sam"
    }

    runtime {
        docker: "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 300 SSD"
    }
}

# SAM/BAM处理
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
        
        echo "=== SAM/BAM处理 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{sam_file}"
        echo "CPU核心数: ~{cpu}"
        echo "内存大小: ~{memory_gb}GB"
        echo "最小MAPQ: ~{min_mapq}"
        echo "开始时间: $(date)"
        
        # 检查samtools版本
        samtools --version
        
        # SAM转BAM并过滤低质量比对
        echo "SAM转BAM并过滤..."
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
        
        # 生成基本统计信息
        echo "生成统计信息..."
        samtools stats ~{sample_name}.sorted.bam > ~{sample_name}_processing_stats.txt
        
        echo "SAM/BAM处理完成"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "输出文件:"
        ls -lh ~{sample_name}.sorted.bam*
        
        # 显示基本统计
        echo "=== 基本统计 ==="
        samtools flagstat ~{sample_name}.sorted.bam
    >>>

    output {
        File sorted_bam = "${sample_name}.sorted.bam"
        File sorted_bam_index = "${sample_name}.sorted.bam.bai"
        File processing_stats = "${sample_name}_processing_stats.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
        memory: "${memory_gb} GB"
        cpu: cpu
        disks: "local-disk 400 SSD"
    }
}

# 重复标记
task MarkDuplicates {
    input {
        File bam_file
        String sample_name
        Int memory_gb
    }

    command <<<
        set -e
        
        echo "=== 重复标记 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{bam_file}"
        echo "内存大小: ~{memory_gb}GB"
        echo "开始时间: $(date)"
        
        # 检查Java版本
        java -version
        
        # 运行Picard MarkDuplicates
        echo "标记重复序列..."
        picard MarkDuplicates \
            I=~{bam_file} \
            O=~{sample_name}.marked.bam \
            M=~{sample_name}_duplicate_metrics.txt \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT \
            REMOVE_DUPLICATES=false \
            ASSUME_SORTED=true
        
        echo "重复标记完成"
        echo "完成时间: $(date)"
        
        # 显示输出文件信息
        echo "输出文件:"
        ls -lh ~{sample_name}.marked.bam*
        
        # 显示重复率统计
        echo "=== 重复率统计 ==="
        head -8 ~{sample_name}_duplicate_metrics.txt
    >>>

    output {
        File marked_bam = "${sample_name}.marked.bam"
        File marked_bam_index = "${sample_name}.marked.bai"
        File duplicate_metrics = "${sample_name}_duplicate_metrics.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/picard:3.0.0--hdfd78af_1"
        memory: "${memory_gb} GB"
        cpu: 4
        disks: "local-disk 500 SSD"
    }
}

# 比对质量评估
task AlignmentQC {
    input {
        File bam_file
        File bam_index
        String sample_name
        Int memory_gb
    }

    command <<<
        set -e
        
        echo "=== 比对质量评估 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入文件: ~{bam_file}"
        echo "内存大小: ~{memory_gb}GB"
        echo "开始时间: $(date)"
        
        # 生成详细的比对统计
        echo "生成比对统计..."
        samtools stats ~{bam_file} > ~{sample_name}_alignment_stats.txt
        samtools flagstat ~{bam_file} >> ~{sample_name}_alignment_stats.txt
        samtools idxstats ~{bam_file} >> ~{sample_name}_alignment_stats.txt
        
        # 运行Qualimap
        echo "运行Qualimap分析..."
        qualimap bamqc \
            -bam ~{bam_file} \
            -outdir ~{sample_name}_qualimap \
            -outformat HTML \
            -nt 4 \
            --java-mem-size=${memory_gb}G
        
        # 打包Qualimap结果
        tar -czf ~{sample_name}_qualimap_report.tar.gz ~{sample_name}_qualimap/
        
        echo "质量评估完成"
        echo "完成时间: $(date)"
        
        # 显示关键统计信息
        echo "=== 关键统计信息 ==="
        grep -E "^SN|reads mapped|reads unmapped|reads properly paired" ~{sample_name}_alignment_stats.txt
    >>>

    output {
        File alignment_stats = "${sample_name}_alignment_stats.txt"
        File qualimap_report = "${sample_name}_qualimap_report.tar.gz"
    }

    runtime {
        docker: "quay.io/biocontainers/qualimap:2.2.2d--1"
        memory: "${memory_gb} GB"
        cpu: 4
        disks: "local-disk 200 SSD"
    }
}
