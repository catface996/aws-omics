version 1.0

workflow PreprocessingWorkflowReadSetV7 {
    meta {
        description: "修复版测序数据预处理工作流v7 - 使用私有ECR镜像"
        version: "7.0"
        author: "AWS Omics Demo Project"
    }

    input {
        String sample_name
        File input_fastq
        File reference_genome
        Int min_length = 50
        Int min_quality = 20
        Int threads = 8
        Int max_length = 500
        Int complexity_threshold = 30
        Boolean enable_polyg_trimming = true
        Boolean enable_polyx_trimming = true
        String dedup_method = "simple"
        
        Int fastqc_memory_gb = 8
        Int fastp_memory_gb = 16
        Int dedup_memory_gb = 16
    }

    call RunFastQC as InitialQC {
        input:
            fastq_file = input_fastq,
            sample_name = sample_name + "_initial",
            memory_gb = fastqc_memory_gb
    }

    call RunFastP {
        input:
            fastq = input_fastq,
            sample_name = sample_name,
            min_length = min_length,
            min_quality = min_quality,
            max_length = max_length,
            complexity_threshold = complexity_threshold,
            enable_polyg_trimming = enable_polyg_trimming,
            enable_polyx_trimming = enable_polyx_trimming,
            threads = threads,
            memory_gb = fastp_memory_gb
    }

    call RemoveDuplicates {
        input:
            fastq = RunFastP.cleaned_fastq,
            sample_name = sample_name,
            method = dedup_method,
            memory_gb = dedup_memory_gb
    }

    call RunFastQC as FinalQC {
        input:
            fastq_file = RemoveDuplicates.dedup_fastq,
            sample_name = sample_name + "_final",
            memory_gb = fastqc_memory_gb
    }

    output {
        File processed_fastq = RemoveDuplicates.dedup_fastq
        File initial_qc_report = InitialQC.fastqc_report
        File final_qc_report = FinalQC.fastqc_report
        File fastp_report = RunFastP.fastp_report
        File dedup_stats = RemoveDuplicates.dedup_stats
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
        apt-get update
        apt-get install -y openjdk-8-jre-headless wget unzip
        
        cd /tmp
        wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
        unzip -q fastqc_v0.11.9.zip
        chmod +x FastQC/fastqc
        
        mkdir -p fastqc_output
        /tmp/FastQC/fastqc \
            --outdir fastqc_output \
            --threads 8 \
            --format fastq \
            ~{fastq_file}
        
        mv fastqc_output/*.html ~{sample_name}_fastqc.html
        mv fastqc_output/*.zip ~{sample_name}_fastqc.zip
    >>>

    output {
        File fastqc_report = "${sample_name}_fastqc.html"
        File fastqc_zip = "${sample_name}_fastqc.zip"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/ubuntu-bioinformatics:20.04"
        memory: "${memory_gb} GB"
        cpu: 8
        disks: "local-disk 20 SSD"
    }
}

task RunFastP {
    input {
        File fastq
        String sample_name
        Int min_length
        Int min_quality
        Int max_length
        Int complexity_threshold
        Boolean enable_polyg_trimming
        Boolean enable_polyx_trimming
        Int threads
        Int memory_gb
    }

    command <<<
        set -e
        export DEBIAN_FRONTEND=noninteractive
        
        apt-get update
        apt-get install -y build-essential zlib1g-dev wget
        
        cd /tmp
        wget -q https://github.com/OpenGene/fastp/archive/v0.23.2.tar.gz
        tar -xzf v0.23.2.tar.gz
        cd fastp-0.23.2
        make -j~{threads}
        
        ./fastp \
            --in1 ~{fastq} \
            --out1 ~{sample_name}_cleaned.fastq.gz \
            --json ~{sample_name}_fastp.json \
            --html ~{sample_name}_fastp.html \
            --thread ~{threads} \
            --length_required ~{min_length} \
            --qualified_quality_phred ~{min_quality} \
            --length_limit ~{max_length} \
            --complexity_threshold ~{complexity_threshold} \
            ~{if enable_polyg_trimming then "--trim_poly_g" else ""} \
            ~{if enable_polyx_trimming then "--trim_poly_x" else ""} \
            --detect_adapter_for_pe \
            --correction \
            --cut_front \
            --cut_tail \
            --cut_window_size 4 \
            --cut_mean_quality 20
        
        mv ~{sample_name}_cleaned.fastq.gz ./
        mv ~{sample_name}_fastp.html ./
        mv ~{sample_name}_fastp.json ./
    >>>

    output {
        File cleaned_fastq = "${sample_name}_cleaned.fastq.gz"
        File fastp_report = "${sample_name}_fastp.html"
        File fastp_json = "${sample_name}_fastp.json"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/ubuntu-bioinformatics:20.04"
        memory: "${memory_gb} GB"
        cpu: "${threads}"
        disks: "local-disk 50 SSD"
    }
}

task RemoveDuplicates {
    input {
        File fastq
        String sample_name
        String method
        Int memory_gb
    }

    command <<<
        set -e
        apt-get update
        apt-get install -y gzip bc
        
        zcat ~{fastq} | \
        awk 'BEGIN{RS="@"; ORS=""} NR>1 {
            getline seq; getline plus; getline qual;
            if (!seen[seq]) {
                print "@" $0 "\n" seq "\n" plus "\n" qual "\n";
                seen[seq] = 1;
            }
        }' | gzip > ~{sample_name}_dedup.fastq.gz
        
        echo "Method: Simple deduplication" > ~{sample_name}_dedup_stats.txt
        echo "Input reads: $(zcat ~{fastq} | wc -l | awk '{print $1/4}')" >> ~{sample_name}_dedup_stats.txt
        echo "Output reads: $(zcat ~{sample_name}_dedup.fastq.gz | wc -l | awk '{print $1/4}')" >> ~{sample_name}_dedup_stats.txt
        
        INPUT_READS=$(zcat ~{fastq} | wc -l | awk '{print $1/4}')
        OUTPUT_READS=$(zcat ~{sample_name}_dedup.fastq.gz | wc -l | awk '{print $1/4}')
        DEDUP_RATE=$(echo "scale=2; (1 - $OUTPUT_READS / $INPUT_READS) * 100" | bc -l)
        echo "Deduplication rate: ${DEDUP_RATE}%" >> ~{sample_name}_dedup_stats.txt
    >>>

    output {
        File dedup_fastq = "${sample_name}_dedup.fastq.gz"
        File dedup_stats = "${sample_name}_dedup_stats.txt"
    }

    runtime {
        docker: "864899854573.dkr.ecr.us-east-1.amazonaws.com/omics/ubuntu-bioinformatics:20.04"
        memory: "${memory_gb} GB"
        cpu: 8
        disks: "local-disk 50 SSD"
    }
}
