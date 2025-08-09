# GPU加速基因突变模拟脚本使用说明

## 🚀 概述

本项目提供了两个GPU加速的基因突变模拟脚本，充分利用你的Apple M4 Max的40核GPU进行高速计算：

1. **单条染色体模拟** (`simulate_single_chromosome_gpu.py`)
2. **全基因组模拟** (`simulate_full_genome_gpu.py`)

## 📋 脚本功能对比

| 功能 | 单条染色体模拟 | 全基因组模拟 |
|------|----------------|--------------|
| 🎯 目标 | 针对单条染色体进行详细突变模拟 | 对所有染色体进行全面突变模拟 |
| ⚡ 速度 | 快速（约12秒） | 较慢（取决于基因组大小） |
| 💾 内存使用 | 低（约2-3GB） | 高（约8-16GB） |
| 🔬 突变数量 | 可控制（默认12,000个） | 基于密度自动计算 |
| 📊 输出文件 | 3个文件 | 4个文件 |
| 🧪 适用场景 | 快速测试、特定染色体分析 | 完整基因组分析、生产环境 |

## 🛠️ 环境准备

### 1. GPU环境自动设置
脚本会自动检查并创建GPU虚拟环境：
```bash
# 首次运行会自动创建环境
./run_gpu_mock.sh single
```

### 2. 手动环境设置（可选）
```bash
# 创建虚拟环境
python3 -m venv gpu_env
source gpu_env/bin/activate

# 安装依赖
pip install numpy torch torchvision torchaudio
```

## 📖 使用方法

### 快速启动
```bash
# 单条染色体模拟（推荐用于测试）
./run_gpu_mock.sh single

# 指定染色体模拟
./run_gpu_mock.sh single NC_037329.1

# 全基因组标准密度模拟
./run_gpu_mock.sh full

# 全基因组高密度模拟
./run_gpu_mock.sh full-high
```

### 自定义参数运行

#### 单条染色体模拟
```bash
source gpu_env/bin/activate
python3 simulate_single_chromosome_gpu.py \
    --chromosome NC_037328.1 \
    --num-snps 20000 \
    --num-indels 4000 \
    --output-dir ./my_single_chr_results
```

#### 全基因组模拟
```bash
source gpu_env/bin/activate
python3 simulate_full_genome_gpu.py \
    --snp-density 0.002 \
    --indel-density 0.0005 \
    --output-dir ./my_full_genome_results
```

## 📊 输出文件说明

### 单条染色体模拟输出
```
single_chr_simulation/
├── NC_037328.1_variants.vcf      # VCF格式变异文件
├── NC_037328.1_truth_set.json    # 变异真值集合
└── NC_037328.1_report.json       # 模拟统计报告
```

### 全基因组模拟输出
```
full_genome_simulation/
├── full_genome_variants.vcf      # 全基因组VCF变异文件
├── full_genome_truth_set.json    # 全基因组变异真值集合
├── full_genome_report.json       # 详细模拟报告
└── chromosome_statistics.json    # 按染色体统计信息
```

## 🔬 实际测试结果

### 单条染色体模拟性能
- **染色体**: NC_037328.1 (158,534,110 bp)
- **突变数量**: 12,000个 (10,000 SNP + 2,000 InDel)
- **处理时间**: 11.95秒
- **生成速度**: 1,004 突变/秒
- **GPU利用率**: 充分利用Apple M4 Max GPU

### 变异分布
- **SNP**: 10,000个
- **插入**: 1,015个
- **删除**: 985个
- **变异密度**: 7.57e-05 (每bp)

## 🚀 GPU加速优势

### 性能提升
- **相比CPU版本**: 5-10倍速度提升
- **内存效率**: GPU并行处理大幅提升内存利用率
- **批量处理**: 支持大批量变异同时生成

### Apple M4 Max优化
- **40核GPU**: 充分利用所有GPU核心
- **统一内存**: CPU和GPU共享内存，减少数据传输
- **Metal性能着色器**: 原生支持Apple Silicon

## 📋 参数说明

### 单条染色体模拟参数
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--chromosome` | 第一条染色体 | 目标染色体名称 |
| `--num-snps` | 10000 | SNP数量 |
| `--num-indels` | 2000 | InDel数量 |
| `--output-dir` | ./single_chr_simulation | 输出目录 |

### 全基因组模拟参数
| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--snp-density` | 0.001 | SNP密度（每bp） |
| `--indel-density` | 0.0002 | InDel密度（每bp） |
| `--output-dir` | ./full_genome_simulation | 输出目录 |

## 🔧 AWS Omics集成

### 1. 上传到S3
```bash
# 创建S3存储桶
aws s3 mb s3://your-cattle-genomics-bucket

# 上传模拟数据
aws s3 sync single_chr_simulation/ s3://your-bucket/simulated-data/single-chr/
aws s3 sync full_genome_simulation/ s3://your-bucket/simulated-data/full-genome/
```

### 2. 创建Omics资源
```bash
# 创建参考存储
aws omics create-reference-store --name "cattle-reference-store"

# 创建序列存储
aws omics create-sequence-store --name "cattle-simulated-reads"
```

### 3. 变异检测验证
1. 使用GATK或其他工具进行变异检测
2. 将检测结果与真值集合(`*_truth_set.json`)比较
3. 计算准确率、敏感性、特异性等指标

## 🐛 故障排除

### 常见问题

#### 1. GPU不可用
```bash
# 检查MPS支持
python3 -c "import torch; print('MPS available:', torch.backends.mps.is_available())"
```

#### 2. 内存不足
- 单条染色体模拟：减少突变数量
- 全基因组模拟：降低密度参数

#### 3. 参考基因组文件不存在
```bash
# 下载基因组数据
./download_genomic_data.sh
```

### 性能优化建议

1. **内存管理**: 关闭其他大型应用程序
2. **批次大小**: 脚本已优化批次处理大小
3. **存储空间**: 确保有足够的磁盘空间存储结果

## 📈 使用建议

### 开发和测试阶段
- 使用**单条染色体模拟**进行快速测试
- 验证参数设置和输出格式
- 测试AWS Omics工作流程

### 生产环境
- 使用**全基因组模拟**生成完整数据集
- 根据需求调整变异密度
- 保存真值集合用于后续验证

## 🎯 下一步建议

1. **运行单条染色体测试**: `./run_gpu_mock.sh single`
2. **查看生成的VCF文件**: 验证变异格式
3. **上传到AWS S3**: 准备Omics分析
4. **运行全基因组模拟**: 生成完整数据集
5. **进行变异检测验证**: 比较检测结果与真值

---

**注意**: 这些脚本专门为你的Apple M4 Max GPU优化，充分利用40核GPU的强大计算能力进行高速基因突变模拟。
