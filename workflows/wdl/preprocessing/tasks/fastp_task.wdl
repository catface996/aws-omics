version 1.0

## FastP 高级质量过滤和长度过滤任务
## 提供比Trimmomatic更精细的过滤控制

task RunFastPPE {
    meta {
        description: "使用FastP进行双端测序数据的高级质量过滤和长度过滤"
        version: "0.23.4"
    }

    parameter_meta {
        fastq_r1: "正向测序文件 (R1)"
        fastq_r2: "反向测序文件 (R2)"
        sample_name: "样本名称"
        min_length: "最小读长阈值"
        min_quality: "最小质量分数"
        threads: "并行线程数"
    }

    input {
        File fastq_r1
        File fastq_r2
        String sample_name
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Int memory_gb = 16
        Int disk_gb = 100
        
        # 高级过滤参数
        Int max_length = 500
        Float length_required_percent = 0.8
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
    }

    command <<<
        set -euo pipefail
        
        echo "🚀 开始FastP双端数据高级过滤..."
        echo "输入文件R1: ~{fastq_r1}"
        echo "输入文件R2: ~{fastq_r2}"
        echo "样本名称: ~{sample_name}"
        echo "最小长度: ~{min_length}"
        echo "最大长度: ~{max_length}"
        echo "最小质量: ~{min_quality}"
        echo "线程数: ~{threads}"
        
        mkdir -p fastp_output
        
        # 构建FastP命令
        fastp_cmd="fastp \
            --in1 ~{fastq_r1} \
            --in2 ~{fastq_r2} \
            --out1 fastp_output/~{sample_name}_R1_filtered.fastq.gz \
            --out2 fastp_output/~{sample_name}_R2_filtered.fastq.gz \
            --html fastp_output/~{sample_name}_fastp.html \
            --json fastp_output/~{sample_name}_fastp.json \
            --thread ~{threads} \
            --length_required ~{min_length} \
            --length_limit ~{max_length} \
            --qualified_quality_phred ~{min_quality} \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            --low_complexity_filter \
            --complexity_threshold ~{complexity_threshold} \
            --overrepresentation_analysis"
        
        # 添加poly-G和poly-X修剪选项
        if [[ "~{enable_polyg_trimming}" == "true" ]]; then
            fastp_cmd="$fastp_cmd --trim_poly_g"
        fi
        
        if [[ "~{enable_polyx_trimming}" == "true" ]]; then
            fastp_cmd="$fastp_cmd --trim_poly_x"
        fi
        
        # 执行FastP
        echo "执行命令: $fastp_cmd"
        eval $fastp_cmd
        
        # 生成详细统计报告
        echo "FastP高级过滤完成" > fastp_output/~{sample_name}_fastp_summary.txt
        echo "处理时间: $(date)" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "过滤参数:" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "  最小长度: ~{min_length} bp" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "  最大长度: ~{max_length} bp" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "  最小质量: ~{min_quality}" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "  复杂度阈值: ~{complexity_threshold}" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "  Poly-G修剪: ~{enable_polyg_trimming}" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "  Poly-X修剪: ~{enable_polyx_trimming}" >> fastp_output/~{sample_name}_fastp_summary.txt
        
        # 统计输出文件
        echo "" >> fastp_output/~{sample_name}_fastp_summary.txt
        echo "输出文件统计:" >> fastp_output/~{sample_name}_fastp_summary.txt
        
        for file in fastp_output/*_filtered.fastq.gz; do
            if [[ -f "$file" ]]; then
                reads=$(zcat "$file" | wc -l | awk '{print $1/4}')
                size=$(ls -lh "$file" | awk '{print $5}')
                echo "  $(basename $file): $reads reads, $size" >> fastp_output/~{sample_name}_fastp_summary.txt
            fi
        done
        
        # 从JSON文件提取关键统计信息
        if [[ -f "fastp_output/~{sample_name}_fastp.json" ]]; then
            echo "" >> fastp_output/~{sample_name}_fastp_summary.txt
            echo "质量统计 (从JSON提取):" >> fastp_output/~{sample_name}_fastp_summary.txt
            
            # 使用python提取JSON统计信息
            python3 << 'EOF' >> fastp_output/~{sample_name}_fastp_summary.txt
import json
import sys

try:
    with open('fastp_output/~{sample_name}_fastp.json', 'r') as f:
        data = json.load(f)
    
    # 提取前后统计信息
    before = data.get('summary', {}).get('before_filtering', {})
    after = data.get('summary', {}).get('after_filtering', {})
    
    print(f"  处理前总reads: {before.get('total_reads', 'N/A'):,}")
    print(f"  处理后总reads: {after.get('total_reads', 'N/A'):,}")
    print(f"  reads保留率: {(after.get('total_reads', 0) / before.get('total_reads', 1) * 100):.2f}%")
    print(f"  处理前平均长度: {before.get('read1_mean_length', 'N/A')} bp")
    print(f"  处理后平均长度: {after.get('read1_mean_length', 'N/A')} bp")
    print(f"  处理前Q30率: {before.get('q30_rate', 'N/A'):.4f}")
    print(f"  处理后Q30率: {after.get('q30_rate', 'N/A'):.4f}")
    
except Exception as e:
    print(f"  JSON解析错误: {e}")
EOF
        fi
        
        echo "✅ FastP高级过滤完成"
        ls -la fastp_output/
    >>>

    runtime {
        docker: "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File filtered_r1 = "fastp_output/${sample_name}_R1_filtered.fastq.gz"
        File filtered_r2 = "fastp_output/${sample_name}_R2_filtered.fastq.gz"
        File fastp_html = "fastp_output/${sample_name}_fastp.html"
        File fastp_json = "fastp_output/${sample_name}_fastp.json"
        File summary = "fastp_output/${sample_name}_fastp_summary.txt"
    }
}

task RunFastPSE {
    meta {
        description: "使用FastP进行单端测序数据的高级质量过滤和长度过滤"
        version: "0.23.4"
    }

    input {
        File fastq
        String sample_name
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Int memory_gb = 8
        Int disk_gb = 50
        Int max_length = 500
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
    }

    command <<<
        set -euo pipefail
        
        echo "🚀 开始FastP单端数据高级过滤..."
        echo "输入文件: ~{fastq}"
        echo "样本名称: ~{sample_name}"
        
        mkdir -p fastp_output
        
        # 构建并执行FastP命令
        fastp \
            --in1 ~{fastq} \
            --out1 fastp_output/~{sample_name}_filtered.fastq.gz \
            --html fastp_output/~{sample_name}_fastp.html \
            --json fastp_output/~{sample_name}_fastp.json \
            --thread ~{threads} \
            --length_required ~{min_length} \
            --length_limit ~{max_length} \
            --qualified_quality_phred ~{min_quality} \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            --low_complexity_filter \
            --complexity_threshold ~{complexity_threshold} \
            --trim_poly_g \
            --trim_poly_x \
            --overrepresentation_analysis
        
        # 生成统计报告
        reads=$(zcat fastp_output/~{sample_name}_filtered.fastq.gz | wc -l | awk '{print $1/4}')
        echo "FastP SE高级过滤完成，输出 $reads reads" > fastp_output/~{sample_name}_fastp_summary.txt
        
        echo "✅ FastP SE高级过滤完成"
    >>>

    runtime {
        docker: "quay.io/biocontainers/fastp:0.23.4--h5f740d0_0"
        memory: "${memory_gb} GB"
        cpu: threads
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File filtered = "fastp_output/${sample_name}_filtered.fastq.gz"
        File fastp_html = "fastp_output/${sample_name}_fastp.html"
        File fastp_json = "fastp_output/${sample_name}_fastp.json"
        File summary = "fastp_output/${sample_name}_fastp_summary.txt"
    }
}
