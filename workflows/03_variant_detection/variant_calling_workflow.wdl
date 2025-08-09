
version 1.0

workflow VariantCallingWorkflow {
    input {
        File reference_genome
        File input_fastq
        String sample_name
        String output_prefix
    }
    
    call AlignReads {
        input:
            reference = reference_genome,
            fastq = input_fastq,
            sample_name = sample_name
    }
    
    call CallVariants {
        input:
            reference = reference_genome,
            aligned_bam = AlignReads.output_bam,
            sample_name = sample_name
    }
    
    call AnnotateVariants {
        input:
            variants_vcf = CallVariants.output_vcf,
            sample_name = sample_name,
            output_prefix = output_prefix
    }
    
    output {
        File aligned_bam = AlignReads.output_bam
        File variants_vcf = CallVariants.output_vcf
        File annotated_vcf = AnnotateVariants.output_vcf
        File analysis_report = AnnotateVariants.report
    }
}

task AlignReads {
    input {
        File reference
        File fastq
        String sample_name
    }
    
    command <<<
        # 建立BWA索引
        bwa index ~{reference}
        
        # 序列比对
        bwa mem -t 4 -R "@RG\tID:~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA" \
            ~{reference} ~{fastq} | \
            samtools view -Sb - | \
            samtools sort -o ~{sample_name}.sorted.bam -
        
        # 建立索引
        samtools index ~{sample_name}.sorted.bam
    >>>
    
    output {
        File output_bam = "~{sample_name}.sorted.bam"
        File output_bai = "~{sample_name}.sorted.bam.bai"
    }
    
    runtime {
        docker: "biocontainers/bwa:v0.7.17_cv1"
        cpu: 4
        memory: "8 GB"
        disks: "local-disk 100 SSD"
    }
}

task CallVariants {
    input {
        File reference
        File aligned_bam
        String sample_name
    }
    
    command <<<
        # 使用GATK进行变异检测
        gatk HaplotypeCaller \
            -R ~{reference} \
            -I ~{aligned_bam} \
            -O ~{sample_name}.variants.vcf \
            --emit-ref-confidence GVCF
    >>>
    
    output {
        File output_vcf = "~{sample_name}.variants.vcf"
    }
    
    runtime {
        docker: "broadinstitute/gatk:latest"
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 50 SSD"
    }
}

task AnnotateVariants {
    input {
        File variants_vcf
        String sample_name
        String output_prefix
    }
    
    command <<<
        # 变异注释和统计
        bcftools stats ~{variants_vcf} > ~{output_prefix}.stats.txt
        
        # 生成注释VCF
        cp ~{variants_vcf} ~{output_prefix}.annotated.vcf
        
        # 生成分析报告
        echo "Variant Analysis Report for ~{sample_name}" > ~{output_prefix}.report.txt
        echo "Generated on: $(date)" >> ~{output_prefix}.report.txt
        echo "" >> ~{output_prefix}.report.txt
        echo "Variant Statistics:" >> ~{output_prefix}.report.txt
        bcftools stats ~{variants_vcf} | grep "^SN" >> ~{output_prefix}.report.txt
    >>>
    
    output {
        File output_vcf = "~{output_prefix}.annotated.vcf"
        File report = "~{output_prefix}.report.txt"
        File stats = "~{output_prefix}.stats.txt"
    }
    
    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 20 SSD"
    }
}
