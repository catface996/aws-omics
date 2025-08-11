# 使用Ubuntu 20.04作为基础镜像，明确指定AMD64架构
FROM --platform=linux/amd64 ubuntu:20.04

# 设置非交互式安装
ENV DEBIAN_FRONTEND=noninteractive

# 更新包列表并安装必要的工具
RUN apt-get update && apt-get install -y \
    wget \
    openjdk-8-jdk \
    python3 \
    python3-pip \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# 设置JAVA_HOME
ENV JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64

# 安装GATK 4.4.0.0
WORKDIR /opt
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip && \
    unzip gatk-4.4.0.0.zip && \
    rm gatk-4.4.0.0.zip && \
    ln -s /opt/gatk-4.4.0.0/gatk /usr/local/bin/gatk

# 安装samtools
WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.17*

# 安装bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.17*

# 设置工作目录
WORKDIR /data

# 验证安装
RUN gatk --version && samtools --version && bcftools --version

# 设置默认命令
CMD ["/bin/bash"]
