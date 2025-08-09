FROM ubuntu:20.04

# 设置非交互模式
ENV DEBIAN_FRONTEND=noninteractive

# 更新包管理器并安装所有需要的软件
RUN apt-get update && apt-get install -y \
    openjdk-8-jre-headless \
    wget \
    unzip \
    gzip \
    bc \
    && rm -rf /var/lib/apt/lists/*

# 复制本地FastQC文件并安装
COPY FastQC-0.11.9.zip /tmp/
RUN cd /opt && \
    unzip -q /tmp/FastQC-0.11.9.zip && \
    chmod +x FastQC-0.11.9/fastqc && \
    ln -s /opt/FastQC-0.11.9/fastqc /usr/local/bin/fastqc && \
    rm /tmp/FastQC-0.11.9.zip

# 设置工作目录
WORKDIR /work

# 默认命令
CMD ["/bin/bash"]
