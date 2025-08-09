version 1.0

## Trimmomatic 接头去除和质量过滤任务
## 支持单端和双端测序数据

task RunTrimmomaticPE {
    meta {
        description: "使用Trimmomatic进行双端测序数据的接头去除和质量过滤"
        version: "0.39"
    }

    parameter_meta {
        fastq_r1: "正向测序文件 (R1)"
        fastq_r2: "反向测序文件 (R2)"
        sample_name: "样本名称"
        adapter_fasta: "接头序列文件"
        min_length: "最小读长阈值"
        min_quality: "最小质量分数"
        threads: "并行线程数"
    }

    input {
        File fastq_r1
        File fastq_r2
        String sample_name
        File? adapter_fasta
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Int memory_gb = 16
        Int disk_gb = 100
    }

    String adapter_file = if defined(adapter_fasta) then 
        "~{adapter_fasta}" else 
        "/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa"

    command <<<
        set -euo pipefail
        
        echo "✂️ 开始Trimmomatic双端数据处理..."
        echo "输入文件R1: ~{fastq_r1}"
        echo "输入文件R2: ~{fastq_r2}"
        echo "样本名称: ~{sample_name}"
        echo "接头文件: ~{adapter_file}"
        echo "最小长度: ~{min_length}"
        echo "最小质量: ~{min_quality}"
        echo "线程数: ~{threads}"
        
        # 创建输出目录
        mkdir -p trimmed_output
        
        # 如果没有提供接头文件，使用默认的TruSeq3-PE接头
        if [[ ! -f "~{adapter_file}" ]]; then
            echo "使用默认TruSeq3-PE接头序列"
            adapter_path="/usr/local/share/trimmomatic/adapters/TruSeq3-PE.fa"
        else
            adapter_path="~{adapter_file}"
        fi
        
        # 运行Trimmomatic PE模式
        trimmomatic PE \
            -threads ~{threads} \
            -phred33 \
            ~{fastq_r1} ~{fastq_r2} \
            trimmed_output/~{sample_name}_R1_paired.fastq.gz \
            trimmed_output/~{sample_name}_R1_unpaired.fastq.gz \
            trimmed_output/~{sample_name}_R2_paired.fastq.gz \
            trimmed_output/~{sample_name}_R2_unpaired.fastq.gz \
            ILLUMINACLIP:${adapter_path}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:~{min_quality} \
            MINLEN:~{min_length} \
            2> trimmed_output/~{sample_name}_trimmomatic.log
        
        # 生成统计报告
        echo "Trimmomatic处理完成" > trimmed_output/~{sample_name}_trimming_summary.txt
        echo "处理时间: $(date)" >> trimmed_output/~{sample_name}_trimming_summary.txt
        echo "参数设置:" >> trimmed_output/~{sample_name}_trimming_summary.txt
        echo "  最小长度: ~{min_length}" >> trimmed_output/~{sample_name}_trimming_summary.txt
        echo "  最小质量: ~{min_quality}" >> trimmed_output/~{sample_name}_trimming_summary.txt
        echo "  线程数: ~{threads}" >> trimmed_output/~{sample_name}_trimming_summary.txt
        
        # 统计输出文件
        echo "" >> trimmed_output/~{sample_name}_trimming_summary.txt
        echo "输出文件统计:" >> trimmed_output/~{sample_name}_trimming_summary.txt
        
        for file in trimmed_output/*.fastq.gz; do
            if [[ -f "$file" ]]; then
                reads=$(zcat "$file" | wc -l | awk '{print $1/4}')
                size=$(ls -lh "$file" | awk '{print $5}')
                echo "  $(basename $file): $reads reads, $size" >> trimmed_output/~{sample_name}_trimming_summary.txt
            fi
        done
        
        echo "✅ Trimmomatic处理完成"
        ls -la trimmed_output/
    >>>

    runtime {
        docker: "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File trimmed_r1 = "trimmed_output/${sample_name}_R1_paired.fastq.gz"
        File trimmed_r2 = "trimmed_output/${sample_name}_R2_paired.fastq.gz"
        File unpaired_r1 = "trimmed_output/${sample_name}_R1_unpaired.fastq.gz"
        File unpaired_r2 = "trimmed_output/${sample_name}_R2_unpaired.fastq.gz"
        File log_file = "trimmed_output/${sample_name}_trimmomatic.log"
        File summary = "trimmed_output/${sample_name}_trimming_summary.txt"
    }
}

task RunTrimmomaticSE {
    meta {
        description: "使用Trimmomatic进行单端测序数据的接头去除和质量过滤"
        version: "0.39"
    }

    input {
        File fastq
        String sample_name
        File? adapter_fasta
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Int memory_gb = 8
        Int disk_gb = 50
    }

    String adapter_file = if defined(adapter_fasta) then 
        "~{adapter_fasta}" else 
        "/usr/local/share/trimmomatic/adapters/TruSeq3-SE.fa"

    command <<<
        set -euo pipefail
        
        echo "✂️ 开始Trimmomatic单端数据处理..."
        echo "输入文件: ~{fastq}"
        echo "样本名称: ~{sample_name}"
        
        mkdir -p trimmed_output
        
        # 运行Trimmomatic SE模式
        trimmomatic SE \
            -threads ~{threads} \
            -phred33 \
            ~{fastq} \
            trimmed_output/~{sample_name}_trimmed.fastq.gz \
            ILLUMINACLIP:~{adapter_file}:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:~{min_quality} \
            MINLEN:~{min_length} \
            2> trimmed_output/~{sample_name}_trimmomatic.log
        
        # 生成统计报告
        reads=$(zcat trimmed_output/~{sample_name}_trimmed.fastq.gz | wc -l | awk '{print $1/4}')
        echo "Trimmomatic SE处理完成，输出 $reads reads" > trimmed_output/~{sample_name}_trimming_summary.txt
        
        echo "✅ Trimmomatic SE处理完成"
    >>>

    runtime {
        docker: "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File trimmed = "trimmed_output/${sample_name}_trimmed.fastq.gz"
        File log_file = "trimmed_output/${sample_name}_trimmomatic.log"
        File summary = "trimmed_output/${sample_name}_trimming_summary.txt"
    }
}
