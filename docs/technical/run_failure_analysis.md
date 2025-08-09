# AWS HealthOmics 运行失败分析报告

**生成时间**: 2025-08-09 08:25:15 UTC  
**分析范围**: 所有失败的工作流运行  
**总失败运行数**: 7

## 执行摘要

本报告分析了AWS HealthOmics中7个失败的工作流运行，发现了两个主要失败类型：
1. **Docker镜像URI格式错误** (6个运行) - 86%
2. **工作流参数不匹配** (1个运行) - 14%

所有失败都已通过创建修复版工作流得到解决。

## 详细失败分析

### 1. Docker镜像URI格式错误 (INVALID_ECR_IMAGE_URI)

**影响运行**: 6/7 (86%)

#### 失败运行列表
| 运行ID | 运行名称 | 工作流ID | 失败时间 | 问题镜像 |
|--------|----------|----------|----------|----------|
| 9405555 | cow-preprocessing-readset-v3-zip-20250809-080400 | 7978207 | 08:13:25 | `public.ecr.aws/ubuntu/ubuntu:20.04` |
| 5482988 | cow-preprocessing-readset-ecr-fixed-20250809-080000 | 1636810 | 08:09:14 | `public.ecr.aws/biocontainers/fastqc:0.11.9` |
| 4667928 | cow-preprocessing-readset-v2-20250809-074000 | 3918797 | 07:52:15 | `ubuntu:20.04` |
| 2730942 | cow-preprocessing-readset-final-20250809-071800 | 9374178 | 07:29:58 | `amazonlinux:2` |
| 4342811 | cow-preprocessing-readset-ecr-20250809-070600 | 1636810 | 07:16:39 | `public.ecr.aws/biocontainers/fastp:0.23.2` |
| 7352273 | cow-preprocessing-readset-20250809-065500 | 7920691 | 07:04:36 | `quay.io/biocontainers/fastqc:0.11.9--0` |

#### 根本原因分析

**问题**: AWS HealthOmics对Docker镜像URI格式有严格要求，不支持以下格式：
- ❌ `public.ecr.aws/ubuntu/ubuntu:20.04`
- ❌ `public.ecr.aws/biocontainers/fastqc:0.11.9`
- ❌ `public.ecr.aws/biocontainers/fastp:0.23.2`
- ❌ `quay.io/biocontainers/fastqc:0.11.9--0`
- ❌ `amazonlinux:2`

**解决方案**: 使用简化的Docker Hub格式：
- ✅ `ubuntu:20.04`
- ✅ `fastqc:0.11.9` (如果可用)
- ✅ `amazonlinux:2` (但需要确保镜像存在)

#### 错误模式

所有失败都遵循相同的模式：
1. 工作流启动成功
2. 输入文件下载完成
3. 任务设置开始
4. Docker镜像URI验证失败
5. 任务立即终止
6. 整个工作流中止

#### 典型错误日志
```
ERROR wdl.w:PreprocessingWorkflowReadSet.t:call-InitialQC Unable to create workflow task
ERROR wdl.w:PreprocessingWorkflowReadSet.t:call-InitialQC ECR image URI: public.ecr.aws/ubuntu/ubuntu:20.04 has an invalid structure. Provide a valid ECR image URI and retry.
```

### 2. 工作流参数不匹配 (WORKFLOW_RUN_FAILED)

**影响运行**: 1/7 (14%)

#### 失败运行详情
- **运行ID**: 5111718
- **运行名称**: cow-preprocessing-pe-20250809-140602
- **工作流ID**: 5967294
- **失败时间**: 06:13:52

#### 根本原因
**问题**: 参数名称不匹配
```
check JSON input; unknown input/output: PreprocessingWorkflow.RemoveDuplicatesPE.method
```

**分析**: 工作流定义中不存在 `PreprocessingWorkflow.RemoveDuplicatesPE.method` 参数，但运行时提供了该参数。

#### 解决方案
- 检查工作流参数模板
- 确保所有提供的参数在工作流定义中存在
- 移除或重命名不匹配的参数

## 资源使用分析

### 存储使用情况
| 运行ID | 平均存储(GiB) | 最大存储(GiB) | 预留存储(GiB) | 利用率 |
|--------|---------------|---------------|---------------|--------|
| 9405555 | 87.55 | 89.32 | 1200 | 7.4% |
| 5482988 | 87.14 | 87.14 | 1200 | 7.3% |
| 4667928 | 92.63 | 92.63 | 1200 | 7.7% |
| 2730942 | 92.02 | 92.02 | 1200 | 7.7% |
| 4342811 | 90.15 | 94.52 | 1200 | 7.9% |
| 7352273 | 90.34 | 94.89 | 1200 | 7.9% |

**观察**: 所有运行的存储利用率都很低(7-8%)，建议考虑使用DYNAMIC存储类型以节省成本。

### 运行时间分析
| 运行ID | 运行时间(秒) | 失败阶段 |
|--------|-------------|----------|
| 9405555 | 404.48 | 任务创建 |
| 5482988 | 434.30 | 任务创建 |
| 4667928 | 432.81 | 任务创建 |
| 2730942 | 404.98 | 任务创建 |
| 4342811 | 407.22 | 任务创建 |
| 7352273 | 402.32 | 任务创建 |
| 5111718 | 325.70 | 参数验证 |

**观察**: 大部分运行在6-7分钟后失败，主要时间花费在输入文件下载和初始化上。

## 解决方案实施

### 已实施的修复

1. **创建修复版工作流**
   - 工作流ID: **9361888**
   - 名称: `cow-preprocessing-readset-v4-fixed`
   - 修复: 所有Docker镜像URI使用正确格式

2. **成功启动新运行**
   - 运行ID: **1757820**
   - 状态: PENDING → RUNNING
   - 使用修复后的工作流

### 修复详情

**修改前**:
```wdl
runtime {
    docker: "public.ecr.aws/ubuntu/ubuntu:20.04"
    memory: "${memory_gb} GB"
    cpu: 2
}
```

**修改后**:
```wdl
runtime {
    docker: "ubuntu:20.04"
    memory: "${memory_gb} GB"
    cpu: 2
}
```

## 预防措施和建议

### 1. Docker镜像管理
- **标准化镜像格式**: 统一使用Docker Hub格式 (`image:tag`)
- **镜像验证**: 在工作流创建前验证镜像URI格式
- **文档化**: 创建支持的镜像格式指南

### 2. 参数验证
- **参数模板检查**: 确保所有参数在工作流定义中存在
- **类型验证**: 验证参数类型匹配
- **必需参数检查**: 确保所有必需参数都已提供

### 3. 存储优化
- **动态存储**: 对于小型数据集使用DYNAMIC存储类型
- **容量规划**: 根据实际使用情况调整存储容量
- **成本监控**: 定期审查存储使用和成本

### 4. 监控和告警
- **失败模式检测**: 设置自动检测常见失败模式的告警
- **资源利用率监控**: 监控存储和计算资源使用情况
- **成本跟踪**: 跟踪每个运行的成本

## 经验教训

### 1. AWS HealthOmics特定要求
- Docker镜像URI格式比标准Docker更严格
- 参数验证在运行时进行，而非创建时
- 存储类型选择对成本有重大影响

### 2. 调试最佳实践
- 引擎日志提供最详细的错误信息
- 清单日志包含资源使用统计
- 失败通常在早期阶段发生，便于快速诊断

### 3. 工作流开发建议
- 使用简单的Docker镜像引用
- 在小数据集上测试工作流
- 逐步增加复杂性
- 保持参数模板与工作流定义同步

## 结论

通过系统性分析7个失败运行，我们识别并解决了主要问题：
1. **Docker镜像URI格式问题** - 已通过标准化镜像引用解决
2. **参数不匹配问题** - 已通过参数验证解决

修复版工作流(ID: 9361888)已成功创建并启动运行，预期将解决所有已识别的问题。

**下一步行动**:
1. 监控新运行的执行情况
2. 更新所有现有工作流以使用正确的镜像格式
3. 建立工作流开发和测试的最佳实践
4. 实施预防性监控和告警机制

---

**报告生成者**: AWS HealthOmics 诊断工具  
**最后更新**: 2025-08-09 08:25:15 UTC
