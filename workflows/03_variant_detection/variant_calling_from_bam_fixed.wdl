version 1.0

workflow VariantCallingFromBAM {
    input {
        File reference_genome
        File input_bam
        File input_bam_index
        String sample_name
        Int gatk_cpu = 16
        Int gatk_memory_gb = 32
        Int bcftools_cpu = 8
        Int bcftools_memory_gb = 16
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
        File raw_variants_vcf_index = CallVariantsGATK.output_vcf_index
        File filtered_variants_vcf = ProcessVariants.filtered_vcf
        File filtered_variants_vcf_index = ProcessVariants.filtered_vcf_index
        File variant_stats = ProcessVariants.stats_file
        File raw_variant_stats = ProcessVariants.raw_stats_file
        File variant_report = ProcessVariants.report_file
        File variant_density = ProcessVariants.density_file
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
        
        echo "=== GATK变异检测开始 ==="
        echo "样本名称: ~{sample_name}"
        echo "参考基因组: ~{reference}"
        echo "输入BAM: ~{input_bam}"
        echo "CPU核心数: ~{cpu}"
        echo "内存: ~{memory_gb}GB"
        echo "GATK版本: $(gatk --version | head -1)"
        echo "Java版本: $(java -version 2>&1 | head -1)"
        
        # 创建参考基因组索引（如果不存在）
        if [ ! -f "~{reference}.fai" ]; then
            echo "创建参考基因组索引..."
            samtools faidx ~{reference}
        else
            echo "参考基因组索引已存在"
        fi
        
        # 创建序列字典（如果不存在）
        DICT_FILE="$(basename ~{reference} .fna).dict"
        if [ ! -f "$DICT_FILE" ]; then
            echo "创建序列字典..."
            gatk CreateSequenceDictionary \
                -R ~{reference} \
                -O $DICT_FILE
        else
            echo "序列字典已存在"
        fi
        
        echo "=== 开始GATK HaplotypeCaller变异检测 ==="
        # 使用GATK HaplotypeCaller进行变异检测 - 修复参数错误
        gatk --java-options "-Xmx~{memory_gb-4}g" HaplotypeCaller \
            -R ~{reference} \
            -I ~{input_bam} \
            -O ~{sample_name}.raw.variants.vcf \
            --native-pair-hmm-threads ~{cpu} \
            --standard-min-confidence-threshold-for-calling 20.0 \
            --verbosity INFO \
            --create-output-variant-index true
        
        echo "=== 压缩和索引VCF文件 ==="
        # 压缩VCF文件
        bgzip ~{sample_name}.raw.variants.vcf
        tabix -p vcf ~{sample_name}.raw.variants.vcf.gz
        
        echo "=== GATK变异检测完成 ==="
        echo "输出文件: ~{sample_name}.raw.variants.vcf.gz"
        
        # 显示基本统计信息
        echo "=== 变异数量统计 ==="
        VARIANT_COUNT=$(zcat ~{sample_name}.raw.variants.vcf.gz | grep -v "^#" | wc -l)
        echo "检测到的变异总数: $VARIANT_COUNT"
        
        # 显示变异类型分布
        echo "=== 变异类型分布 ==="
        zcat ~{sample_name}.raw.variants.vcf.gz | grep -v "^#" | cut -f4,5 | \
        awk '{
            ref_len = length($1)
            alt_len = length($2)
            if (ref_len == 1 && alt_len == 1) print "SNP"
            else if (ref_len > alt_len) print "DEL"
            else if (ref_len < alt_len) print "INS"
            else print "COMPLEX"
        }' | sort | uniq -c | sort -nr
    >>>
    
    output {
        File output_vcf = "~{sample_name}.raw.variants.vcf.gz"
        File output_vcf_index = "~{sample_name}.raw.variants.vcf.gz.tbi"
    }
    
    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/variant-calling:amd64"
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk 400 SSD"
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
        
        echo "=== 变异后处理开始 ==="
        echo "样本名称: ~{sample_name}"
        echo "输入VCF: ~{variants_vcf}"
        echo "BCFtools版本: $(bcftools --version | head -1)"
        
        # 基本变异过滤
        echo "=== 应用质量过滤 ==="
        echo "过滤条件: QUAL>=20 && DP>=10"
        bcftools filter \
            -i 'QUAL>=20 && DP>=10' \
            -O z \
            -o ~{sample_name}.filtered.variants.vcf.gz \
            ~{variants_vcf}
        
        # 创建索引
        echo "=== 创建索引 ==="
        tabix -p vcf ~{sample_name}.filtered.variants.vcf.gz
        
        # 生成详细统计
        echo "=== 生成统计信息 ==="
        bcftools stats ~{variants_vcf} > ~{sample_name}.raw.stats.txt
        bcftools stats ~{sample_name}.filtered.variants.vcf.gz > ~{sample_name}.filtered.stats.txt
        
        # 生成分析报告
        echo "=== 生成分析报告 ==="
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

## 染色体分布
$(bcftools stats ~{sample_name}.filtered.variants.vcf.gz | grep "^CHR")

## 过滤效果对比
原始变异数量: $(zcat ~{variants_vcf} | grep -v "^#" | wc -l)
过滤后变异数量: $(zcat ~{sample_name}.filtered.variants.vcf.gz | grep -v "^#" | wc -l)
过滤率: $(echo "scale=2; ($(zcat ~{variants_vcf} | grep -v "^#" | wc -l) - $(zcat ~{sample_name}.filtered.variants.vcf.gz | grep -v "^#" | wc -l)) * 100 / $(zcat ~{variants_vcf} | grep -v "^#" | wc -l)" | bc)%
EOF
        
        # 生成变异密度统计
        echo "=== 生成变异密度统计 ==="
        bcftools query -f '%CHROM\t%POS\t%QUAL\t%DP\t%REF\t%ALT\t%FILTER\n' ~{sample_name}.filtered.variants.vcf.gz | \
            head -10000 > ~{sample_name}.variant.density.txt
        
        echo "=== 变异后处理完成 ==="
        FILTERED_COUNT=$(zcat ~{sample_name}.filtered.variants.vcf.gz | grep -v "^#" | wc -l)
        echo "过滤后变异数量: $FILTERED_COUNT"
        
        # 显示过滤后变异类型分布
        echo "=== 过滤后变异类型分布 ==="
        zcat ~{sample_name}.filtered.variants.vcf.gz | grep -v "^#" | cut -f4,5 | \
        awk '{
            ref_len = length($1)
            alt_len = length($2)
            if (ref_len == 1 && alt_len == 1) print "SNP"
            else if (ref_len > alt_len) print "DEL"
            else if (ref_len < alt_len) print "INS"
            else print "COMPLEX"
        }' | sort | uniq -c | sort -nr
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
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/variant-calling:amd64"
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk 200 SSD"
    }
}
