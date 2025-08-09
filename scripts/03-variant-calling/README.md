# 03-变异检测 (Variant Calling)

## 目录说明
此目录包含基因组变异检测相关的脚本和工具。

## 主要功能
- SNP (单核苷酸多态性) 检测
- InDel (插入缺失) 检测
- 结构变异检测
- 变异质量过滤和筛选
- 变异格式标准化

## 预期脚本类型
- GATK HaplotypeCaller 脚本
- FreeBayes 变异检测脚本
- VarScan 变异检测脚本
- 变异质量过滤脚本
- VCF 文件处理脚本

## 输入
来自 02-alignment 的 BAM 文件

## 输出
- VCF (Variant Call Format) 文件
- 变异统计报告
- 质量过滤报告

## 执行顺序
1. 变异检测预处理
2. 变异调用
3. 变异质量评估
4. 变异过滤
5. VCF 文件标准化

## 下一步
输出的 VCF 文件将用于 04-variant-annotation 步骤。
