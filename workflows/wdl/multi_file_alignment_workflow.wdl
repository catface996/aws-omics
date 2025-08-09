version 1.0

# AWS Omics多文件基因组比对工作流
# 支持将一头牛的多个基因文件进行统一比对分析

workflow CattleMultiFileAlignment {
    input {
        # 参考基因组
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        
        # 多个测序文件（同一头牛的不同染色体或不同测序批次）
        Array[File] fastq_r1_files
        Array[File] fastq_r2_files
        Array[String] read_group_names
        
        # 样本信息
        String sample_name
        String subject_id
        
        # 比对参数
        Int cpu_count = 16
        Int memory_gb = 32
        String docker_image = "public.ecr.aws/aws-genomics/bwa:0.7.17"
    }
    
    # 并行处理每对FASTQ文件
    scatter (i in range(length(fastq_r1_files))) {
        call AlignReads {
            input:
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                fastq_r1 = fastq_r1_files[i],
                fastq_r2 = fastq_r2_files[i],
                read_group_name = read_group_names[i],
                sample_name = sample_name,
                cpu_count = cpu_count,
                memory_gb = memory_gb,
                docker_image = docker_image
        }
    }
    
    # 合并所有比对结果
    call MergeBamFiles {
        input:
            bam_files = AlignReads.aligned_bam,
            bam_indexes = AlignReads.aligned_bam_index,
            sample_name = sample_name,
            cpu_count = cpu_count,
            memory_gb = memory_gb
    }
    
    # 标记重复序列
    call MarkDuplicates {
        input:
            input_bam = MergeBamFiles.merged_bam,
            sample_name = sample_name,
            cpu_count = cpu_count,
            memory_gb = memory_gb
    }
    
    # 变异检测
    call CallVariants {
        input:
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_dict = reference_dict,
            input_bam = MarkDuplicates.output_bam,
            input_bam_index = MarkDuplicates.output_bam_index,
            sample_name = sample_name,
            cpu_count = cpu_count,
            memory_gb = memory_gb
    }
    
    output {
        File final_bam = MarkDuplicates.output_bam
        File final_bam_index = MarkDuplicates.output_bam_index
        File variants_vcf = CallVariants.output_vcf
        File variants_vcf_index = CallVariants.output_vcf_index
        File alignment_metrics = MarkDuplicates.metrics_file
    }
}

# 比对单个FASTQ文件对
task AlignReads {
    input {
        File reference_fasta
        File reference_fasta_fai
        File fastq_r1
        File fastq_r2
        String read_group_name
        String sample_name
        Int cpu_count
        Int memory_gb
        String docker_image
    }
    
    String base_name = basename(fastq_r1, "_R1.fastq.gz")
    
    command <<<
        set -euo pipefail
        
        # 创建BWA索引（如果不存在）
        if [ ! -f "~{reference_fasta}.bwt" ]; then
            bwa index ~{reference_fasta}
        fi
        
        # BWA比对
        bwa mem -t ~{cpu_count} \
            -R "@RG\tID:~{read_group_name}\tSM:~{sample_name}\tPL:ILLUMINA\tLB:~{sample_name}_lib" \
            ~{reference_fasta} \
            ~{fastq_r1} \
            ~{fastq_r2} | \
        samtools sort -@ ~{cpu_count} -o ~{base_name}.sorted.bam -
        
        # 创建索引
        samtools index ~{base_name}.sorted.bam
    >>>
    
    output {
        File aligned_bam = "~{base_name}.sorted.bam"
        File aligned_bam_index = "~{base_name}.sorted.bam.bai"
    }
    
    runtime {
        docker: docker_image
        cpu: cpu_count
        memory: "~{memory_gb} GB"
        disks: "local-disk 100 SSD"
    }
}

# 合并多个BAM文件
task MergeBamFiles {
    input {
        Array[File] bam_files
        Array[File] bam_indexes
        String sample_name
        Int cpu_count
        Int memory_gb
    }
    
    command <<<
        set -euo pipefail
        
        # 使用samtools merge合并BAM文件
        samtools merge -@ ~{cpu_count} \
            ~{sample_name}.merged.bam \
            ~{sep=' ' bam_files}
        
        # 创建索引
        samtools index ~{sample_name}.merged.bam
    >>>
    
    output {
        File merged_bam = "~{sample_name}.merged.bam"
        File merged_bam_index = "~{sample_name}.merged.bam.bai"
    }
    
    runtime {
        docker: "public.ecr.aws/aws-genomics/samtools:1.15.1"
        cpu: cpu_count
        memory: "~{memory_gb} GB"
        disks: "local-disk 200 SSD"
    }
}

# 标记重复序列
task MarkDuplicates {
    input {
        File input_bam
        String sample_name
        Int cpu_count
        Int memory_gb
    }
    
    command <<<
        set -euo pipefail
        
        # 使用Picard MarkDuplicates
        java -Xmx~{memory_gb-4}g -jar /usr/picard/picard.jar MarkDuplicates \
            INPUT=~{input_bam} \
            OUTPUT=~{sample_name}.dedup.bam \
            METRICS_FILE=~{sample_name}.dedup_metrics.txt \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT
    >>>
    
    output {
        File output_bam = "~{sample_name}.dedup.bam"
        File output_bam_index = "~{sample_name}.dedup.bai"
        File metrics_file = "~{sample_name}.dedup_metrics.txt"
    }
    
    runtime {
        docker: "public.ecr.aws/aws-genomics/picard:2.27.4"
        cpu: cpu_count
        memory: "~{memory_gb} GB"
        disks: "local-disk 300 SSD"
    }
}

# 变异检测
task CallVariants {
    input {
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        File input_bam
        File input_bam_index
        String sample_name
        Int cpu_count
        Int memory_gb
    }
    
    command <<<
        set -euo pipefail
        
        # 使用GATK HaplotypeCaller进行变异检测
        gatk --java-options "-Xmx~{memory_gb-4}g" HaplotypeCaller \
            -R ~{reference_fasta} \
            -I ~{input_bam} \
            -O ~{sample_name}.variants.vcf.gz \
            --native-pair-hmm-threads ~{cpu_count}
        
        # 创建索引
        tabix -p vcf ~{sample_name}.variants.vcf.gz
    >>>
    
    output {
        File output_vcf = "~{sample_name}.variants.vcf.gz"
        File output_vcf_index = "~{sample_name}.variants.vcf.gz.tbi"
    }
    
    runtime {
        docker: "public.ecr.aws/aws-genomics/gatk:4.3.0.0"
        cpu: cpu_count
        memory: "~{memory_gb} GB"
        disks: "local-disk 100 SSD"
    }
}
