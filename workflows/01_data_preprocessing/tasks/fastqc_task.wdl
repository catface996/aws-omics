version 1.0

## FastQC 质量评估任务
## 用于测序数据的质量控制和评估

task RunFastQC {
    meta {
        description: "使用FastQC进行测序数据质量评估"
        version: "0.12.1"
    }

    parameter_meta {
        fastq_file: "输入的FASTQ文件"
        sample_name: "样本名称"
        threads: "并行线程数"
        memory_gb: "内存需求(GB)"
        disk_gb: "磁盘空间需求(GB)"
    }

    input {
        File fastq_file
        String sample_name
        Int threads = 4
        Int memory_gb = 8
        Int disk_gb = 50
    }

    String fastq_basename = basename(fastq_file, ".fastq.gz")
    String output_prefix = "${sample_name}_${fastq_basename}"

    command <<<
        set -euo pipefail
        
        echo "🔍 开始FastQC质量评估..."
        echo "输入文件: ~{fastq_file}"
        echo "样本名称: ~{sample_name}"
        echo "线程数: ~{threads}"
        
        # 创建输出目录
        mkdir -p fastqc_output
        
        # 运行FastQC
        fastqc \
            --threads ~{threads} \
            --outdir fastqc_output \
            --format fastq \
            --noextract \
            ~{fastq_file}
        
        # 重命名输出文件以包含样本名称
        original_html=$(find fastqc_output -name "*.html" | head -1)
        original_zip=$(find fastqc_output -name "*.zip" | head -1)
        
        if [[ -f "$original_html" ]]; then
            mv "$original_html" "fastqc_output/~{output_prefix}_fastqc.html"
        fi
        
        if [[ -f "$original_zip" ]]; then
            mv "$original_zip" "fastqc_output/~{output_prefix}_fastqc.zip"
        fi
        
        # 生成简要统计信息
        echo "FastQC分析完成" > fastqc_output/~{output_prefix}_summary.txt
        echo "输入文件: $(basename ~{fastq_file})" >> fastqc_output/~{output_prefix}_summary.txt
        echo "分析时间: $(date)" >> fastqc_output/~{output_prefix}_summary.txt
        
        # 提取关键质量指标
        if [[ -f "fastqc_output/~{output_prefix}_fastqc.zip" ]]; then
            unzip -q "fastqc_output/~{output_prefix}_fastqc.zip" -d temp_extract
            
            # 提取总序列数
            total_sequences=$(grep "Total Sequences" temp_extract/*/fastqc_data.txt | cut -f2)
            echo "总序列数: $total_sequences" >> fastqc_output/~{output_prefix}_summary.txt
            
            # 提取序列长度
            sequence_length=$(grep "Sequence length" temp_extract/*/fastqc_data.txt | cut -f2)
            echo "序列长度: $sequence_length" >> fastqc_output/~{output_prefix}_summary.txt
            
            # 提取GC含量
            gc_content=$(grep "%GC" temp_extract/*/fastqc_data.txt | cut -f2)
            echo "GC含量: $gc_content%" >> fastqc_output/~{output_prefix}_summary.txt
            
            rm -rf temp_extract
        fi
        
        echo "✅ FastQC分析完成"
        ls -la fastqc_output/
    >>>

    runtime {
        docker: "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File fastqc_html = "fastqc_output/${output_prefix}_fastqc.html"
        File fastqc_zip = "fastqc_output/${output_prefix}_fastqc.zip"
        File summary = "fastqc_output/${output_prefix}_summary.txt"
    }
}
