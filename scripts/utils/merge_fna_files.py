#!/usr/bin/env python3
"""
FNA文件合并工具
将分散的染色体FNA文件合并成一个完整的基因组文件

特性:
- 按染色体编号排序合并
- 保持FASTA格式规范
- 生成合并统计报告
- 支持大文件高效处理
"""

import os
import sys
import time
import json
from pathlib import Path
import argparse
from typing import List, Dict
import re

class FNAMerger:
    def __init__(self, input_dir: str, output_file: str):
        self.input_dir = Path(input_dir)
        self.output_file = Path(output_file)
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # 统计信息
        self.stats = {
            'files_processed': 0,
            'total_sequences': 0,
            'total_bases': 0,
            'chromosomes': [],
            'processing_time': 0
        }

    def find_fna_files(self) -> List[Path]:
        """查找并排序FNA文件"""
        print("🔍 查找FNA文件...")
        
        fna_files = list(self.input_dir.glob("*.fna"))
        
        if not fna_files:
            print(f"❌ 在 {self.input_dir} 中未找到FNA文件")
            return []
        
        # 按染色体编号排序
        def extract_chromosome_number(filename):
            # 提取NC_037328.1格式中的数字
            match = re.search(r'NC_0373(\d+)\.1', filename.name)
            if match:
                return int(match.group(1))
            return 999  # 未匹配的文件排在最后
        
        fna_files.sort(key=extract_chromosome_number)
        
        print(f"✅ 找到 {len(fna_files)} 个FNA文件")
        for i, file in enumerate(fna_files[:5]):  # 显示前5个
            print(f"  {i+1}. {file.name}")
        if len(fna_files) > 5:
            print(f"  ... 还有 {len(fna_files)-5} 个文件")
        
        return fna_files

    def get_file_info(self, fna_file: Path) -> Dict:
        """获取FNA文件信息"""
        info = {
            'file': fna_file.name,
            'size': fna_file.stat().st_size,
            'sequences': 0,
            'bases': 0,
            'headers': []
        }
        
        try:
            with open(fna_file, 'r') as f:
                current_seq_length = 0
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_seq_length > 0:
                            info['bases'] += current_seq_length
                            current_seq_length = 0
                        info['sequences'] += 1
                        info['headers'].append(line)
                    else:
                        current_seq_length += len(line)
                
                # 添加最后一个序列的长度
                if current_seq_length > 0:
                    info['bases'] += current_seq_length
        
        except Exception as e:
            print(f"⚠️  读取文件 {fna_file.name} 时出错: {e}")
        
        return info

    def merge_files(self, fna_files: List[Path]) -> None:
        """合并FNA文件"""
        print(f"\n🔗 开始合并 {len(fna_files)} 个FNA文件...")
        print(f"📁 输出文件: {self.output_file}")
        
        start_time = time.time()
        
        with open(self.output_file, 'w') as outf:
            for i, fna_file in enumerate(fna_files):
                print(f"  📄 处理文件 {i+1}/{len(fna_files)}: {fna_file.name}")
                
                # 获取文件信息
                file_info = self.get_file_info(fna_file)
                self.stats['chromosomes'].append(file_info)
                
                # 复制文件内容
                try:
                    with open(fna_file, 'r') as inf:
                        for line in inf:
                            outf.write(line)
                    
                    # 确保文件间有换行分隔
                    if i < len(fna_files) - 1:  # 不是最后一个文件
                        outf.write('\n')
                    
                    self.stats['files_processed'] += 1
                    self.stats['total_sequences'] += file_info['sequences']
                    self.stats['total_bases'] += file_info['bases']
                    
                    print(f"    ✅ 完成: {file_info['sequences']} 序列, {file_info['bases']:,} bp")
                
                except Exception as e:
                    print(f"    ❌ 处理文件 {fna_file.name} 时出错: {e}")
        
        self.stats['processing_time'] = time.time() - start_time
        
        print(f"\n🎉 合并完成!")
        print(f"📊 处理文件: {self.stats['files_processed']} 个")
        print(f"🧬 总序列数: {self.stats['total_sequences']} 个")
        print(f"📏 总碱基数: {self.stats['total_bases']:,} bp")
        print(f"⏱️  处理时间: {self.stats['processing_time']:.2f}s")
        print(f"💾 输出大小: {self.output_file.stat().st_size / (1024**3):.2f} GB")

    def generate_report(self) -> None:
        """生成合并报告"""
        report_file = self.output_file.parent / f"{self.output_file.stem}_merge_report.json"
        
        report = {
            'merge_summary': {
                'input_directory': str(self.input_dir),
                'output_file': str(self.output_file),
                'files_processed': self.stats['files_processed'],
                'total_sequences': self.stats['total_sequences'],
                'total_bases': self.stats['total_bases'],
                'processing_time': self.stats['processing_time'],
                'output_size_bytes': self.output_file.stat().st_size,
                'output_size_gb': self.output_file.stat().st_size / (1024**3)
            },
            'chromosome_details': self.stats['chromosomes']
        }
        
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"📋 详细报告: {report_file}")

    def validate_merged_file(self) -> bool:
        """验证合并后的文件"""
        print(f"\n🔍 验证合并后的文件...")
        
        try:
            sequence_count = 0
            base_count = 0
            
            with open(self.output_file, 'r') as f:
                current_seq_length = 0
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_seq_length > 0:
                            base_count += current_seq_length
                            current_seq_length = 0
                        sequence_count += 1
                    else:
                        current_seq_length += len(line)
                
                # 添加最后一个序列
                if current_seq_length > 0:
                    base_count += current_seq_length
            
            print(f"✅ 验证通过:")
            print(f"  📊 序列数: {sequence_count} (预期: {self.stats['total_sequences']})")
            print(f"  📏 碱基数: {base_count:,} (预期: {self.stats['total_bases']:,})")
            
            if sequence_count == self.stats['total_sequences'] and base_count == self.stats['total_bases']:
                print(f"🎯 数据完整性验证成功!")
                return True
            else:
                print(f"⚠️  数据不匹配，请检查合并过程")
                return False
        
        except Exception as e:
            print(f"❌ 验证过程出错: {e}")
            return False

    def merge_fna_files(self) -> bool:
        """主合并流程"""
        print(f"🚀 开始FNA文件合并")
        print(f"📁 输入目录: {self.input_dir}")
        print(f"📄 输出文件: {self.output_file}")
        
        # 查找文件
        fna_files = self.find_fna_files()
        if not fna_files:
            return False
        
        # 合并文件
        self.merge_files(fna_files)
        
        # 验证结果
        validation_success = self.validate_merged_file()
        
        # 生成报告
        self.generate_report()
        
        return validation_success


def main():
    parser = argparse.ArgumentParser(description="FNA文件合并工具")
    parser.add_argument("--input-dir", "-i", required=True, 
                       help="包含FNA文件的输入目录")
    parser.add_argument("--output", "-o", required=True,
                       help="合并后的输出文件路径")
    parser.add_argument("--validate", action="store_true",
                       help="验证合并后的文件完整性")
    
    args = parser.parse_args()
    
    # 检查输入目录
    if not os.path.exists(args.input_dir):
        print(f"❌ 输入目录不存在: {args.input_dir}")
        sys.exit(1)
    
    # 创建合并器并执行
    merger = FNAMerger(args.input_dir, args.output)
    success = merger.merge_fna_files()
    
    if success:
        print(f"\n🎉 FNA文件合并成功完成!")
        sys.exit(0)
    else:
        print(f"\n❌ FNA文件合并失败!")
        sys.exit(1)


if __name__ == "__main__":
    main()
