version 1.0

## 去重复序列任务
## 使用多种策略去除重复的测序reads

task RemoveDuplicatesPE {
    meta {
        description: "使用FastUniq去除双端测序数据中的重复序列"
        version: "1.1"
    }

    parameter_meta {
        fastq_r1: "正向测序文件 (R1)"
        fastq_r2: "反向测序文件 (R2)"
        sample_name: "样本名称"
        threads: "并行线程数"
        method: "去重复方法: fastuniq, clumpify, or seqkit"
    }

    input {
        File fastq_r1
        File fastq_r2
        String sample_name
        Int threads = 8
        Int memory_gb = 16
        Int disk_gb = 100
        String method = "fastuniq"  # fastuniq, clumpify, seqkit
    }

    command <<<
        set -euo pipefail
        
        echo "🔄 开始去重复序列处理 (双端)..."
        echo "输入文件R1: ~{fastq_r1}"
        echo "输入文件R2: ~{fastq_r2}"
        echo "样本名称: ~{sample_name}"
        echo "去重复方法: ~{method}"
        echo "线程数: ~{threads}"
        
        mkdir -p dedup_output
        
        # 统计原始reads数量
        original_reads_r1=$(zcat ~{fastq_r1} | wc -l | awk '{print $1/4}')
        original_reads_r2=$(zcat ~{fastq_r2} | wc -l | awk '{print $1/4}')
        
        echo "原始reads数量 R1: $original_reads_r1"
        echo "原始reads数量 R2: $original_reads_r2"
        
        if [[ "~{method}" == "fastuniq" ]]; then
            echo "使用FastUniq进行去重复..."
            
            # 创建输入文件列表
            echo "~{fastq_r1}" > input_list.txt
            echo "~{fastq_r2}" >> input_list.txt
            
            # 运行FastUniq
            fastuniq \
                -i input_list.txt \
                -o dedup_output/~{sample_name}_R1_dedup.fastq \
                -p dedup_output/~{sample_name}_R2_dedup.fastq
            
            # 压缩输出文件
            gzip dedup_output/~{sample_name}_R1_dedup.fastq
            gzip dedup_output/~{sample_name}_R2_dedup.fastq
            
        elif [[ "~{method}" == "clumpify" ]]; then
            echo "使用Clumpify进行去重复..."
            
            clumpify.sh \
                in1=~{fastq_r1} \
                in2=~{fastq_r2} \
                out1=dedup_output/~{sample_name}_R1_dedup.fastq.gz \
                out2=dedup_output/~{sample_name}_R2_dedup.fastq.gz \
                dedupe=t \
                optical=t \
                threads=~{threads}
                
        elif [[ "~{method}" == "seqkit" ]]; then
            echo "使用SeqKit进行去重复..."
            
            # 对R1和R2分别去重复，然后同步
            seqkit rmdup \
                --by-seq \
                --threads ~{threads} \
                ~{fastq_r1} \
                -o dedup_output/temp_R1_dedup.fastq.gz
            
            seqkit rmdup \
                --by-seq \
                --threads ~{threads} \
                ~{fastq_r2} \
                -o dedup_output/temp_R2_dedup.fastq.gz
            
            # 同步配对reads (这里简化处理，实际应该更复杂)
            mv dedup_output/temp_R1_dedup.fastq.gz dedup_output/~{sample_name}_R1_dedup.fastq.gz
            mv dedup_output/temp_R2_dedup.fastq.gz dedup_output/~{sample_name}_R2_dedup.fastq.gz
        fi
        
        # 统计去重复后的reads数量
        dedup_reads_r1=$(zcat dedup_output/~{sample_name}_R1_dedup.fastq.gz | wc -l | awk '{print $1/4}')
        dedup_reads_r2=$(zcat dedup_output/~{sample_name}_R2_dedup.fastq.gz | wc -l | awk '{print $1/4}')
        
        # 计算去重复率
        duplicate_rate_r1=$(echo "scale=4; ($original_reads_r1 - $dedup_reads_r1) / $original_reads_r1 * 100" | bc)
        duplicate_rate_r2=$(echo "scale=4; ($original_reads_r2 - $dedup_reads_r2) / $original_reads_r2 * 100" | bc)
        
        # 生成详细统计报告
        echo "去重复处理完成" > dedup_output/~{sample_name}_dedup_summary.txt
        echo "处理时间: $(date)" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "去重复方法: ~{method}" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "统计信息:" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "  R1 原始reads: $original_reads_r1" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "  R1 去重复后: $dedup_reads_r1" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "  R1 重复率: ${duplicate_rate_r1}%" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "  R2 原始reads: $original_reads_r2" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "  R2 去重复后: $dedup_reads_r2" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "  R2 重复率: ${duplicate_rate_r2}%" >> dedup_output/~{sample_name}_dedup_summary.txt
        
        # 文件大小统计
        echo "" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "文件大小:" >> dedup_output/~{sample_name}_dedup_summary.txt
        for file in dedup_output/*.fastq.gz; do
            if [[ -f "$file" ]]; then
                size=$(ls -lh "$file" | awk '{print $5}')
                echo "  $(basename $file): $size" >> dedup_output/~{sample_name}_dedup_summary.txt
            fi
        done
        
        echo "✅ 去重复处理完成"
        echo "重复率 R1: ${duplicate_rate_r1}%, R2: ${duplicate_rate_r2}%"
        ls -la dedup_output/
    >>>

    runtime {
        docker: "quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File dedup_r1 = "dedup_output/${sample_name}_R1_dedup.fastq.gz"
        File dedup_r2 = "dedup_output/${sample_name}_R2_dedup.fastq.gz"
        File log_file = "dedup_output/${sample_name}_dedup_summary.txt"
    }
}

task RemoveDuplicatesSE {
    meta {
        description: "使用SeqKit去除单端测序数据中的重复序列"
        version: "2.3.1"
    }

    input {
        File fastq
        String sample_name
        Int threads = 8
        Int memory_gb = 8
        Int disk_gb = 50
        String method = "seqkit"
    }

    command <<<
        set -euo pipefail
        
        echo "🔄 开始去重复序列处理 (单端)..."
        echo "输入文件: ~{fastq}"
        echo "样本名称: ~{sample_name}"
        echo "方法: ~{method}"
        
        mkdir -p dedup_output
        
        # 统计原始reads数量
        original_reads=$(zcat ~{fastq} | wc -l | awk '{print $1/4}')
        echo "原始reads数量: $original_reads"
        
        if [[ "~{method}" == "seqkit" ]]; then
            # 使用SeqKit去重复
            seqkit rmdup \
                --by-seq \
                --threads ~{threads} \
                ~{fastq} \
                -o dedup_output/~{sample_name}_dedup.fastq.gz
                
        elif [[ "~{method}" == "clumpify" ]]; then
            # 使用Clumpify去重复
            clumpify.sh \
                in=~{fastq} \
                out=dedup_output/~{sample_name}_dedup.fastq.gz \
                dedupe=t \
                threads=~{threads}
        fi
        
        # 统计去重复后的reads数量
        dedup_reads=$(zcat dedup_output/~{sample_name}_dedup.fastq.gz | wc -l | awk '{print $1/4}')
        duplicate_rate=$(echo "scale=4; ($original_reads - $dedup_reads) / $original_reads * 100" | bc)
        
        # 生成统计报告
        echo "去重复处理完成" > dedup_output/~{sample_name}_dedup_summary.txt
        echo "原始reads: $original_reads" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "去重复后: $dedup_reads" >> dedup_output/~{sample_name}_dedup_summary.txt
        echo "重复率: ${duplicate_rate}%" >> dedup_output/~{sample_name}_dedup_summary.txt
        
        echo "✅ 去重复处理完成，重复率: ${duplicate_rate}%"
    >>>

    runtime {
        docker: "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File dedup = "dedup_output/${sample_name}_dedup.fastq.gz"
        File log_file = "dedup_output/${sample_name}_dedup_summary.txt"
    }
}
