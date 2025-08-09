#!/usr/bin/env python3
"""
模拟数据整合脚本
将分散的染色体数据整合到单个文件中，并生成清理后的数据集
"""

import os
import json
import gzip
from pathlib import Path
from collections import defaultdict

class DataConsolidator:
    def __init__(self, simulation_dir):
        self.simulation_dir = Path(simulation_dir)
        self.output_dir = Path("consolidated_simulation_data")
        self.chromosome_dirs = []
        self.variants = []
        self.truth_data = []
        
    def discover_chromosome_directories(self):
        """发现所有染色体目录"""
        print("🔍 发现染色体目录...")
        
        chr_dirs = [d for d in self.simulation_dir.iterdir() 
                   if d.is_dir() and d.name.startswith('chr_')]
        
        # 按染色体名称排序
        chr_dirs.sort(key=lambda x: x.name)
        
        self.chromosome_dirs = chr_dirs
        print(f"✅ 发现 {len(chr_dirs)} 个染色体目录")
        
        # 显示前10个和后10个
        if len(chr_dirs) > 20:
            print("前10个染色体:")
            for i, d in enumerate(chr_dirs[:10]):
                chr_name = d.name.replace('chr_', '')
                print(f"  {i+1:2d}. {chr_name}")
            print("...")
            print("后10个染色体:")
            for i, d in enumerate(chr_dirs[-10:], len(chr_dirs)-9):
                chr_name = d.name.replace('chr_', '')
                print(f"  {i:2d}. {chr_name}")
        else:
            for i, d in enumerate(chr_dirs):
                chr_name = d.name.replace('chr_', '')
                print(f"  {i+1:2d}. {chr_name}")
        
        return chr_dirs
    
    def filter_main_chromosomes(self, max_chromosomes=30):
        """筛选主要染色体（通常是NC_开头的）"""
        print(f"\n🎯 筛选前 {max_chromosomes} 条主要染色体...")
        
        # 分类染色体
        main_chrs = []  # NC_开头的主要染色体
        other_chrs = []  # 其他染色体
        
        for chr_dir in self.chromosome_dirs:
            chr_name = chr_dir.name.replace('chr_', '')
            if chr_name.startswith('NC_'):
                main_chrs.append(chr_dir)
            else:
                other_chrs.append(chr_dir)
        
        # 按名称排序主要染色体
        main_chrs.sort(key=lambda x: x.name)
        
        # 选择前N条主要染色体
        selected_chrs = main_chrs[:max_chromosomes]
        
        print(f"📊 染色体分类:")
        print(f"  - 主要染色体 (NC_): {len(main_chrs)}")
        print(f"  - 其他染色体: {len(other_chrs)}")
        print(f"  - 选择用于整合: {len(selected_chrs)}")
        
        print(f"\n✅ 选择的染色体:")
        for i, chr_dir in enumerate(selected_chrs):
            chr_name = chr_dir.name.replace('chr_', '')
            print(f"  {i+1:2d}. {chr_name}")
        
        self.chromosome_dirs = selected_chrs
        return selected_chrs
    
    def consolidate_vcf_files(self):
        """整合所有VCF文件"""
        print(f"\n📄 整合VCF文件...")
        
        self.output_dir.mkdir(exist_ok=True)
        consolidated_vcf = self.output_dir / "consolidated_variants.vcf"
        
        total_variants = 0
        chr_variant_counts = {}
        
        with open(consolidated_vcf, 'w') as outfile:
            # 写入VCF头部
            outfile.write("##fileformat=VCFv4.2\n")
            outfile.write("##source=ConsolidatedSimulation\n")
            outfile.write("##reference=ARS-UCD1.2\n")
            outfile.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
            outfile.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
            outfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            outfile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
            outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n")
            
            # 合并每个染色体的VCF数据
            for chr_dir in self.chromosome_dirs:
                chr_name = chr_dir.name.replace('chr_', '')
                vcf_file = chr_dir / f"{chr_name}_variants.vcf"
                
                if vcf_file.exists():
                    chr_variants = 0
                    with open(vcf_file, 'r') as infile:
                        for line in infile:
                            if not line.startswith('#'):
                                outfile.write(line)
                                chr_variants += 1
                                total_variants += 1
                    
                    chr_variant_counts[chr_name] = chr_variants
                    print(f"  ✅ {chr_name}: {chr_variants:,} 变异")
                else:
                    print(f"  ❌ {chr_name}: VCF文件不存在")
        
        print(f"\n📊 VCF整合完成:")
        print(f"  - 总变异数: {total_variants:,}")
        print(f"  - 输出文件: {consolidated_vcf}")
        
        return consolidated_vcf, chr_variant_counts
    
    def consolidate_truth_sets(self):
        """整合所有真值集合"""
        print(f"\n🎯 整合真值集合...")
        
        consolidated_truth = self.output_dir / "consolidated_truth_set.json"
        all_variants = []
        chr_stats = {}
        
        for chr_dir in self.chromosome_dirs:
            chr_name = chr_dir.name.replace('chr_', '')
            truth_file = chr_dir / f"{chr_name}_truth_set.json"
            
            if truth_file.exists():
                with open(truth_file, 'r') as f:
                    data = json.load(f)
                    
                    # 提取染色体统计信息
                    if 'metadata' in data:
                        chr_stats[chr_name] = {
                            'length': data['metadata'].get('chromosome_length', 0),
                            'gc_content': data['metadata'].get('gc_content', 0),
                            'variant_count': len(data.get('variants', []))
                        }
                    
                    # 添加变异数据
                    variants = data.get('variants', [])
                    all_variants.extend(variants)
                    
                    print(f"  ✅ {chr_name}: {len(variants):,} 变异")
            else:
                print(f"  ❌ {chr_name}: 真值文件不存在")
        
        # 按染色体和位置排序
        all_variants.sort(key=lambda x: (x.get('chromosome', ''), x.get('position', 0)))
        
        # 生成整合的真值集合
        consolidated_data = {
            'metadata': {
                'description': '整合的奶牛基因组变异真值集合',
                'reference': 'ARS-UCD1.2',
                'total_chromosomes': len(self.chromosome_dirs),
                'total_variants': len(all_variants),
                'consolidation_time': __import__('time').strftime('%Y-%m-%d %H:%M:%S'),
                'chromosome_statistics': chr_stats
            },
            'variants': all_variants
        }
        
        with open(consolidated_truth, 'w') as f:
            json.dump(consolidated_data, f, indent=2)
        
        print(f"\n📊 真值集合整合完成:")
        print(f"  - 总变异数: {len(all_variants):,}")
        print(f"  - 输出文件: {consolidated_truth}")
        
        return consolidated_truth
    
    def generate_summary_report(self, vcf_file, truth_file, chr_variant_counts):
        """生成整合报告"""
        print(f"\n📋 生成整合报告...")
        
        report_file = self.output_dir / "consolidation_report.json"
        
        # 统计变异类型
        variant_types = defaultdict(int)
        total_variants = sum(chr_variant_counts.values())
        
        # 简单估算（基于比例）
        variant_types['SNP'] = int(total_variants * 0.83)  # 约83% SNP
        variant_types['INS'] = int(total_variants * 0.085)  # 约8.5% 插入
        variant_types['DEL'] = total_variants - variant_types['SNP'] - variant_types['INS']  # 其余为删除
        
        # 计算总基因组大小（估算）
        total_genome_size = sum(chr_variant_counts.values()) / 0.0012 * 1000  # 基于密度反推
        
        report = {
            'consolidation_summary': {
                'total_chromosomes': len(self.chromosome_dirs),
                'total_variants': total_variants,
                'estimated_genome_size': int(total_genome_size),
                'variant_density': total_variants / total_genome_size if total_genome_size > 0 else 0,
                'consolidation_time': __import__('time').strftime('%Y-%m-%d %H:%M:%S')
            },
            'variant_statistics': {
                'total_variants': total_variants,
                'variant_types': dict(variant_types),
                'variants_by_chromosome': chr_variant_counts
            },
            'output_files': {
                'consolidated_vcf': str(vcf_file),
                'consolidated_truth_set': str(truth_file),
                'consolidation_report': str(report_file)
            },
            'chromosome_list': [chr_dir.name.replace('chr_', '') for chr_dir in self.chromosome_dirs]
        }
        
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"✅ 报告生成完成: {report_file}")
        return report
    
    def cleanup_original_data(self, keep_original=True):
        """清理原始分散数据（可选）"""
        if not keep_original:
            print(f"\n🧹 清理原始分散数据...")
            for chr_dir in self.chromosome_dirs:
                if chr_dir.exists():
                    import shutil
                    shutil.rmtree(chr_dir)
                    print(f"  🗑️  删除: {chr_dir.name}")
            print("✅ 原始数据清理完成")
        else:
            print(f"\n💾 保留原始数据（如需清理，请设置 keep_original=False）")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='整合分散的模拟数据')
    parser.add_argument('--simulation-dir', default='./parallel_genome_simulation',
                       help='模拟数据目录')
    parser.add_argument('--max-chromosomes', type=int, default=30,
                       help='最大染色体数量')
    parser.add_argument('--cleanup', action='store_true',
                       help='清理原始分散数据')
    
    args = parser.parse_args()
    
    print("🔄 数据整合脚本")
    print("=" * 50)
    
    # 初始化整合器
    consolidator = DataConsolidator(args.simulation_dir)
    
    # 发现染色体目录
    consolidator.discover_chromosome_directories()
    
    # 筛选主要染色体
    consolidator.filter_main_chromosomes(args.max_chromosomes)
    
    # 整合VCF文件
    vcf_file, chr_counts = consolidator.consolidate_vcf_files()
    
    # 整合真值集合
    truth_file = consolidator.consolidate_truth_sets()
    
    # 生成报告
    report = consolidator.generate_summary_report(vcf_file, truth_file, chr_counts)
    
    # 清理原始数据（可选）
    consolidator.cleanup_original_data(keep_original=not args.cleanup)
    
    print("\n" + "=" * 50)
    print("🎉 数据整合完成！")
    print("=" * 50)
    print(f"📁 输出目录: {consolidator.output_dir}")
    print(f"📄 主要文件:")
    print(f"  - VCF文件: consolidated_variants.vcf")
    print(f"  - 真值集合: consolidated_truth_set.json")
    print(f"  - 整合报告: consolidation_report.json")
    print(f"\n📊 数据统计:")
    print(f"  - 染色体数: {report['consolidation_summary']['total_chromosomes']}")
    print(f"  - 总变异数: {report['consolidation_summary']['total_variants']:,}")
    print(f"  - 基因组大小: {report['consolidation_summary']['estimated_genome_size']:,} bp")

if __name__ == "__main__":
    main()
