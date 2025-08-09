version 1.0

## MultiQC 综合质量报告任务
## 整合所有质量控制工具的输出，生成统一的HTML报告

task RunMultiQC {
    meta {
        description: "使用MultiQC生成综合质量控制报告"
        version: "1.15"
    }

    parameter_meta {
        input_files: "所有质量控制工具的输出文件"
        sample_name: "样本名称"
        output_name: "输出报告名称"
    }

    input {
        Array[File] input_files
        String sample_name
        String output_name
        Int memory_gb = 8
        Int disk_gb = 50
    }

    command <<<
        set -euo pipefail
        
        echo "📊 开始生成MultiQC综合报告..."
        echo "样本名称: ~{sample_name}"
        echo "输出名称: ~{output_name}"
        echo "输入文件数量: ~{length(input_files)}"
        
        # 创建工作目录
        mkdir -p multiqc_input
        mkdir -p multiqc_output
        
        # 复制所有输入文件到工作目录
        file_count=0
        for file in ~{sep=' ' input_files}; do
            if [[ -f "$file" ]]; then
                cp "$file" multiqc_input/
                file_count=$((file_count + 1))
                echo "复制文件: $(basename $file)"
            fi
        done
        
        echo "成功复制 $file_count 个文件"
        
        # 创建MultiQC配置文件
        cat > multiqc_config.yaml << 'EOF'
title: "~{sample_name} 数据预处理质量报告"
subtitle: "奶牛基因组测序数据预处理流水线"
intro_text: "本报告展示了测序数据预处理各个步骤的质量控制结果，包括原始数据质量、接头去除、质量过滤、长度过滤和去重复处理的效果。"

report_comment: "生成时间: $(date)"

# 自定义样本名称
sample_names_rename:
    - ["_R1_initial", " (R1 原始)"]
    - ["_R2_initial", " (R2 原始)"]
    - ["_R1_final", " (R1 最终)"]
    - ["_R2_final", " (R2 最终)"]
    - ["_trimmed", " (修剪后)"]
    - ["_filtered", " (过滤后)"]
    - ["_dedup", " (去重复后)"]

# 模块顺序
module_order:
    - fastqc:
        name: "FastQC (原始数据)"
        path_filters:
            - "*initial*"
    - trimmomatic
    - fastp
    - fastqc:
        name: "FastQC (最终数据)"
        path_filters:
            - "*final*"

# 自定义颜色方案
colors:
    - "#1f77b4"  # 蓝色
    - "#ff7f0e"  # 橙色
    - "#2ca02c"  # 绿色
    - "#d62728"  # 红色
    - "#9467bd"  # 紫色
    - "#8c564b"  # 棕色

# 表格配置
table_columns_visible:
    FastQC:
        percent_duplicates: True
        percent_gc: True
        avg_sequence_length: True
        percent_fails: True
        total_sequences: True

# 删除不需要的部分
remove_sections:
    - fastqc_sequence_counts
    - fastqc_sequence_duplication_levels

# 文件名清理
fn_clean_exts:
    - ".fastq"
    - ".fastq.gz"
    - ".fq"
    - ".fq.gz"
    - "_fastqc"
    - "_trimmed"
    - "_filtered"
    - "_dedup"

# 数据目录
data_dir: True
EOF
        
        # 运行MultiQC
        multiqc \
            multiqc_input/ \
            --config multiqc_config.yaml \
            --outdir multiqc_output \
            --filename ~{output_name} \
            --title "~{sample_name} 预处理质量报告" \
            --comment "奶牛基因组测序数据预处理流水线质量控制报告" \
            --force \
            --verbose
        
        # 生成处理摘要
        echo "MultiQC报告生成完成" > multiqc_output/~{output_name}_summary.txt
        echo "生成时间: $(date)" >> multiqc_output/~{output_name}_summary.txt
        echo "样本名称: ~{sample_name}" >> multiqc_output/~{output_name}_summary.txt
        echo "处理文件数: $file_count" >> multiqc_output/~{output_name}_summary.txt
        echo "" >> multiqc_output/~{output_name}_summary.txt
        echo "包含的分析模块:" >> multiqc_output/~{output_name}_summary.txt
        
        # 检查生成的文件中包含哪些模块
        if ls multiqc_input/*fastqc* >/dev/null 2>&1; then
            echo "  ✅ FastQC - 测序质量评估" >> multiqc_output/~{output_name}_summary.txt
        fi
        
        if ls multiqc_input/*trimmomatic* >/dev/null 2>&1; then
            echo "  ✅ Trimmomatic - 接头去除和质量过滤" >> multiqc_output/~{output_name}_summary.txt
        fi
        
        if ls multiqc_input/*fastp* >/dev/null 2>&1; then
            echo "  ✅ FastP - 高级质量过滤和长度过滤" >> multiqc_output/~{output_name}_summary.txt
        fi
        
        if ls multiqc_input/*dedup* >/dev/null 2>&1; then
            echo "  ✅ 去重复处理统计" >> multiqc_output/~{output_name}_summary.txt
        fi
        
        # 检查输出文件
        echo "" >> multiqc_output/~{output_name}_summary.txt
        echo "生成的文件:" >> multiqc_output/~{output_name}_summary.txt
        for file in multiqc_output/*; do
            if [[ -f "$file" ]]; then
                size=$(ls -lh "$file" | awk '{print $5}')
                echo "  $(basename $file): $size" >> multiqc_output/~{output_name}_summary.txt
            fi
        done
        
        echo "✅ MultiQC报告生成完成"
        echo "📊 报告文件: multiqc_output/~{output_name}.html"
        ls -la multiqc_output/
    >>>

    runtime {
        docker: "quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0"
        memory: "${memory_gb} GB"
        cpu: 2
        disks: "local-disk ${disk_gb} SSD"
        preemptible: 2
        maxRetries: 1
    }

    output {
        File multiqc_html = "multiqc_output/${output_name}.html"
        File multiqc_data = "multiqc_output/${output_name}_data.zip"
        File summary = "multiqc_output/${output_name}_summary.txt"
    }
}
