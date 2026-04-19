# CNVkit WES CNV Pipeline（阶段一）

本项目当前实现 **阶段一：构建 pooled normal 的 CNVkit reference（CNVKIT_REFERENCE）**，使用 Slurm 调度并支持并行计算 normal coverage。

## 已实现阶段

- `--stage build-reference`：已实现。
- `--stage run-cnv`：仅预留接口，暂未实现（阶段二）。

## 输入格式

`samples.tsv` 必须为 **TAB 分隔**，且表头必须严格为：

```tsv
sample_id	tumor_bam	normal_bam
```

说明：
- 阶段一仅使用 `sample_id` 与 `normal_bam` 构建 reference。
- 但会严格校验 `tumor_bam` 与 `normal_bam` 均存在且可读。

## 快速开始

```bash
bash 01_submit_cnvkit_pipeline.sh \
  --samples /path/to/samples.tsv \
  --stage build-reference \
  --mode wes \
  --threads 8 \
  --max-parallel 8
```

可选参数：
- `--reference-out`：指定 reference 输出路径。
- `--workdir`：指定中间目录与日志目录。
- `--overwrite-reference`：如果 reference 已存在则强制重建。
- `--partition` / `SLURM_PARTITION`：仅允许 `cpu1` 或 `cpu2`，默认 `cpu1`。

## 阶段一执行顺序

1. `01_submit_cnvkit_pipeline.sh`
   - 校验 `samples.tsv`
   - 去重 `normal_bam`
   - 提交 Slurm array（并行 coverage）
   - 提交依赖 array 的 pooled reference 汇总任务
2. `run_cnvkit_reference_array.sbatch`
   - 根据数组索引取 1 个 unique normal BAM
   - 调用单任务脚本
3. `02_run_cnvkit_reference_task.sh`
   - 生成 antitarget BED（存在则跳过）
   - 计算 target/antitarget coverage（存在则跳过）
4. `02_build_cnvkit_reference.sh`
   - 汇总所有 coverage `.cnn` 文件生成 pooled reference

## 关键 CNVkit 子命令

- `cnvkit.py antitarget`
- `cnvkit.py coverage`
- `cnvkit.py reference`

## 目录结构

- `config/00_config.sh`：统一管理 hg38 资源与默认参数
- `lib/validate_samples_tsv.sh`：严格校验输入 TSV
- `lib/steps_cnvkit.sh`：封装 CNVkit 阶段一核心步骤
- `01_submit_cnvkit_pipeline.sh`：阶段入口与 Slurm 提交
- `run_cnvkit_reference_array.sbatch`：数组任务包装
- `02_run_cnvkit_reference_task.sh`：单个 normal coverage 任务
- `02_build_cnvkit_reference.sh`：pooled reference 汇总任务

## 日志与中间文件

默认位于 `./work/cnvkit/`：
- `meta/`：校验后样本表、去重 normal 列表、antitarget BED
- `coverage/`：每个 normal 的 `.targetcoverage.cnn` 与 `.antitargetcoverage.cnn`
- `logs/`：Slurm array 与汇总任务日志
