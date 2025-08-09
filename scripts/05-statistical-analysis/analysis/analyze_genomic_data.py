#!/usr/bin/env python3
"""
奶牛基因组数据分析入门脚本
用于分析下载的基因组数据并生成基本统计信息
"""

import os
import gzip
import sys
from pathlib import Path
from collections import defaultdict, Counter
import argparse

def analyze_fasta_file(file_path):
    """分析FASTA文件的基本统计信息"""
    print(f"\n分析文件: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"错误: 文件 {file_path} 不存在")
        return
    
    # 判断是否为压缩文件
    is_gzipped = file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    try:
        with open_func(file_path, mode) as f:
            sequences = {}
            current_seq_id = None
            current_seq = []
            
            line_count = 0
            for line in f:
                line_count += 1
                line = line.strip()
                
                if line.startswith('>'):
                    # 保存前一个序列
                    if current_seq_id:
                        sequences[current_seq_id] = ''.join(current_seq)
                    
                    # 开始新序列
                    current_seq_id = line[1:].split()[0]  # 取第一个空格前的部分作为ID
                    current_seq = []
                else:
                    current_seq.append(line)
                
                # 限制读取行数以避免内存问题
                if line_count > 100000:
                    print("注意: 文件较大，只分析前100,000行")
                    break
            
            # 保存最后一个序列
            if current_seq_id:
                sequences[current_seq_id] = ''.join(current_seq)
        
        # 生成统计信息
        print(f"序列数量: {len(sequences)}")
        
        if sequences:
            seq_lengths = [len(seq) for seq in sequences.values()]
            total_length = sum(seq_lengths)
            
            print(f"总长度: {total_length:,} bp")
            print(f"平均长度: {total_length/len(sequences):,.1f} bp")
            print(f"最长序列: {max(seq_lengths):,} bp")
            print(f"最短序列: {min(seq_lengths):,} bp")
            
            # 显示前几个序列的信息
            print("\n前5个序列:")
            for i, (seq_id, seq) in enumerate(list(sequences.items())[:5]):
                print(f"  {seq_id}: {len(seq):,} bp")
            
            # 分析碱基组成（仅分析第一个序列的前1000个碱基）
            if sequences:
                first_seq = list(sequences.values())[0][:1000]
                base_count = Counter(first_seq.upper())
                
                print(f"\n碱基组成分析 (前1000个碱基):")
                for base in ['A', 'T', 'G', 'C', 'N']:
                    count = base_count.get(base, 0)
                    percentage = (count / len(first_seq)) * 100 if first_seq else 0
                    print(f"  {base}: {count} ({percentage:.1f}%)")
    
    except Exception as e:
        print(f"分析文件时出错: {e}")

def analyze_gtf_file(file_path):
    """分析GTF注释文件的基本统计信息"""
    print(f"\n分析GTF注释文件: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"错误: 文件 {file_path} 不存在")
        return
    
    is_gzipped = file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'
    
    try:
        feature_counts = defaultdict(int)
        chromosomes = set()
        gene_count = 0
        
        with open_func(file_path, mode) as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chromosome = parts[0]
                    feature_type = parts[2]
                    
                    chromosomes.add(chromosome)
                    feature_counts[feature_type] += 1
                    
                    if feature_type == 'gene':
                        gene_count += 1
                
                # 限制读取行数
                if line_num > 50000:
                    print("注意: 文件较大，只分析前50,000行")
                    break
        
        print(f"染色体/序列数量: {len(chromosomes)}")
        print(f"基因数量: {gene_count}")
        
        print("\n特征类型统计:")
        for feature_type, count in sorted(feature_counts.items(), key=lambda x: x[1], reverse=True):
            print(f"  {feature_type}: {count}")
        
        print("\n前10个染色体/序列:")
        for chrom in sorted(list(chromosomes))[:10]:
            print(f"  {chrom}")
    
    except Exception as e:
        print(f"分析GTF文件时出错: {e}")

def main():
    parser = argparse.ArgumentParser(description='分析奶牛基因组数据')
    parser.add_argument('--data-dir', default='./genomic_data', 
                       help='基因组数据目录路径 (默认: ./genomic_data)')
    parser.add_argument('--quick', action='store_true',
                       help='快速分析模式，只分析关键文件')
    
    args = parser.parse_args()
    
    data_dir = Path(args.data_dir)
    
    if not data_dir.exists():
        print(f"错误: 数据目录 {data_dir} 不存在")
        print("请先运行 ./download_genomic_data.sh 下载数据")
        sys.exit(1)
    
    print("=" * 60)
    print("奶牛基因组数据分析报告")
    print("=" * 60)
    
    # 分析参考基因组
    reference_file = data_dir / "reference" / "GCF_002263795.1_ARS-UCD1.2_genomic.fna.gz"
    if reference_file.exists():
        analyze_fasta_file(str(reference_file))
    
    # 分析GTF注释文件
    gtf_file = data_dir / "annotation" / "GCF_002263795.1_ARS-UCD1.2_genomic.gtf.gz"
    if gtf_file.exists():
        analyze_gtf_file(str(gtf_file))
    
    if not args.quick:
        # 分析蛋白质序列
        protein_file = data_dir / "protein" / "GCF_002263795.1_ARS-UCD1.2_protein.faa.gz"
        if protein_file.exists():
            analyze_fasta_file(str(protein_file))
        
        # 分析编码序列
        cds_file = data_dir / "annotation" / "GCF_002263795.1_ARS-UCD1.2_cds_from_genomic.fna.gz"
        if cds_file.exists():
            analyze_fasta_file(str(cds_file))
    
    # 显示文件大小信息
    print(f"\n文件大小统计:")
    for file_path in data_dir.rglob("*.gz"):
        size_mb = file_path.stat().st_size / (1024 * 1024)
        print(f"  {file_path.name}: {size_mb:.1f} MB")
    
    print("\n" + "=" * 60)
    print("分析完成！")
    print("\n使用建议:")
    print("1. 解压文件: gunzip genomic_data/*/*.gz")
    print("2. 使用samtools为参考基因组建立索引:")
    print("   samtools faidx GCF_002263795.1_ARS-UCD1.2_genomic.fna")
    print("3. 使用BWA建立比对索引:")
    print("   bwa index GCF_002263795.1_ARS-UCD1.2_genomic.fna")
    print("4. 查看组装报告了解更多信息:")
    print("   less genomic_data/assembly/GCF_002263795.1_ARS-UCD1.2_assembly_report.txt")

if __name__ == "__main__":
    main()
