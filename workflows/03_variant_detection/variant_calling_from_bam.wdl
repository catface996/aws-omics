version 1.0

workflow VariantCallingFromBAM {
    input {
        File reference_genome
        File input_bam
        File input_bam_index
        String sample_name
        Int gatk_cpu = 8
        Int gatk_memory_gb = 16
        Int bcftools_cpu = 4
        Int bcftools_memory_gb = 8
    }
    
    call CallVariantsGATK {
        input:
            reference = reference_genome,
            input_bam = input_bam,
            input_bam_index = input_bam_index,
            sample_name = sample_name,
            cpu = gatk_cpu,
            memory_gb = gatk_memory_gb
    }
    
    call ProcessVariants {
        input:
            variants_vcf = CallVariantsGATK.output_vcf,
            sample_name = sample_name,
            cpu = bcftools_cpu,
            memory_gb = bcftools_memory_gb
    }
    
    output {
        File raw_variants_vcf = CallVariantsGATK.output_vcf
        File filtered_variants_vcf = ProcessVariants.filtered_vcf
        File variant_stats = ProcessVariants.stats_file
        File variant_report = ProcessVariants.report_file
    }
}

task CallVariantsGATK {
    input {
        File reference
        File input_bam
        File input_bam_index
        String sample_name
        Int cpu
        Int memory_gb
    }
    
    command <<<
        set -euo pipefail
        
        # 创建参考基因组索引（如果不存在）
        if [ ! -f "~{reference}.fai" ]; then
            samtools faidx ~{reference}
        fi
        
        # 创建序列字典（如果不存在）
        if [ ! -f "$(basename ~{reference} .fna).dict" ]; then
            gatk CreateSequenceDictionary \
                -R ~{reference} \
                -O $(basename ~{reference} .fna).dict
        fi
        
        # 使用GATK HaplotypeCaller进行变异检测
        gatk --java-options "-Xmx~{memory_gb-2}g" HaplotypeCaller \
            -R ~{reference} \
            -I ~{input_bam} \
            -O ~{sample_name}.raw.variants.vcf \
            --native-pair-hmm-threads ~{cpu} \
            --standard-min-confidence-threshold-for-calling 20.0 \
            --standard-min-confidence-threshold-for-emitting 10.0
        
        # 压缩VCF文件
        bgzip ~{sample_name}.raw.variants.vcf
        tabix -p vcf ~{sample_name}.raw.variants.vcf.gz
    >>>
    
    output {
        File output_vcf = "~{sample_name}.raw.variants.vcf.gz"
        File output_vcf_index = "~{sample_name}.raw.variants.vcf.gz.tbi"
    }
    
    runtime {
        docker: "broadinstitute/gatk:4.4.0.0"
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk 200 SSD"
    }
}

task ProcessVariants {
    input {
        File variants_vcf
        String sample_name
        Int cpu
        Int memory_gb
    }
    
    command <<<
        set -euo pipefail
        
        # 基本变异过滤
        bcftools filter \
            -i 'QUAL>=20 && DP>=10' \
            -O z \
            -o ~{sample_name}.filtered.variants.vcf.gz \
            ~{variants_vcf}
        
        # 创建索引
        tabix -p vcf ~{sample_name}.filtered.variants.vcf.gz
        
        # 生成详细统计
        bcftools stats ~{variants_vcf} > ~{sample_name}.raw.stats.txt
        bcftools stats ~{sample_name}.filtered.variants.vcf.gz > ~{sample_name}.filtered.stats.txt
        
        # 生成分析报告
        cat > ~{sample_name}.variant.report.txt << EOF
# 变异检测分析报告
样本名称: ~{sample_name}
分析时间: $(date)
分析工具: GATK HaplotypeCaller + BCFtools

## 原始变异统计
$(bcftools stats ~{variants_vcf} | grep "^SN" | head -20)

## 过滤后变异统计  
$(bcftools stats ~{sample_name}.filtered.variants.vcf.gz | grep "^SN" | head -20)

## 变异类型分布
$(bcftools stats ~{sample_name}.filtered.variants.vcf.gz | grep "^ST")

## 质量分布
$(bcftools stats ~{sample_name}.filtered.variants.vcf.gz | grep "^QUAL")
EOF
        
        # 生成变异密度统计
        bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\n' ~{sample_name}.filtered.variants.vcf.gz | \
            head -1000 > ~{sample_name}.variant.density.txt
    >>>
    
    output {
        File filtered_vcf = "~{sample_name}.filtered.variants.vcf.gz"
        File filtered_vcf_index = "~{sample_name}.filtered.variants.vcf.gz.tbi"
        File stats_file = "~{sample_name}.filtered.stats.txt"
        File raw_stats_file = "~{sample_name}.raw.stats.txt"
        File report_file = "~{sample_name}.variant.report.txt"
        File density_file = "~{sample_name}.variant.density.txt"
    }
    
    runtime {
        docker: "biocontainers/bcftools:v1.17_cv1"
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk 100 SSD"
    }
}
