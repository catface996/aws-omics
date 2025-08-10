version 1.0

workflow PreprocessingWorkflowReadSetV10 {
    meta {
        description: "修复版测序数据预处理工作流v10 - 使用私有ECR FastQC镜像"
        version: "10.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File input_fastq
        Int fastqc_memory_gb = 8
    }

    # 步骤1: 原始数据质量评估
    call RunFastQC as InitialQC {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name + "_qc",
            memory_gb = fastqc_memory_gb
    }

    output {
        File qc_report = InitialQC.fastqc_report
        File qc_zip = InitialQC.fastqc_zip
        File original_fastq = input_fastq
    }
}

task RunFastQC {
    input {
        File fastq_file
        String sample_name
        Int memory_gb
    }

    command <<<
        set -e
        
        echo "开始FastQC质量评估..."
        echo "输入文件: ~{fastq_file}"
        echo "样本名称: ~{sample_name}"
        
        # 检查FastQC版本和可用性
        fastqc --version
        
        # 创建输出目录
        mkdir -p fastqc_output
        
        # 运行FastQC
        fastqc \
            --outdir fastqc_output \
            --threads 8 \
            --format fastq \
            --extract \
            ~{fastq_file}
        
        echo "FastQC分析完成，处理输出文件..."
        
        # 列出输出文件
        ls -la fastqc_output/
        
        # 重命名输出文件
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
        
        echo "输出文件准备完成"
        ls -la ~{sample_name}_fastqc.*
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/fastqc:0.12.1"
        memory: "${memory_gb} GB"
        cpu: 8
        disks: "local-disk 20 SSD"
    }
}
