version 1.0

workflow PreprocessingWorkflowReadSetV8 {
    meta {
        description: "简化版测序数据预处理工作流v8 - 仅FastQC质量评估"
        version: "8.0"
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
        
        # 创建输出目录
        mkdir -p fastqc_output
        
        # 运行FastQC（已预装在镜像中）
        fastqc \
            --outdir fastqc_output \
            --threads 8 \
            --format fastq \
            ~{fastq_file}
        
        # 重命名输出文件
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/ubuntu-bioinformatics:latest"
        memory: "${memory_gb} GB"
        cpu: 8
        disks: "local-disk 20 SSD"
    }
}
