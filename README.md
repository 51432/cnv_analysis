# CNVkit WES Pipeline（Slurm）

当前项目包含两个阶段：

1. **阶段一（已完成）**：构建 pooled normal `wes_reference.cnn`
2. **阶段二（本次新增）**：对每个 tumor 样本进行 CNVkit 后续分析（coverage → fix → segment → 可选 filtering/call/reports）

---

## 输入格式（固定）

`samples.tsv` 必须是 TAB 分隔，且表头严格为：

```tsv
sample_id	tumor_bam	normal_bam
```

校验规则：
- `sample_id` 唯一
- `tumor_bam` / `normal_bam` 必须存在且可读
- `tumor_bam` 与 `normal_bam` 不能相同

错误格式统一：

```text
[ERROR] line=<行号> sample_id=<ID> message=<原因>
```

---

## 阶段一：build reference（完整指令）

入口脚本：`01_submit_cnvkit_pipeline.sh`

如果你还没有 reference，可运行：

```bash
bash 01_submit_cnvkit_pipeline.sh   --samples samples.tsv   --stage build-reference   --mode wes   --threads 8   --max-parallel 8   --reference-out /data/person/wup/public/liusy_files/reference_genomes/hg38/resources/cnvkit/wes_reference.cnn
```

常用可选参数：
- `--workdir`：阶段一中间目录（默认 `./work/cnvkit`）
- `--overwrite-reference`：强制重建 reference（会清理 reference/coverage/BED 缓存）
- `--partition`：`cpu1` 或 `cpu2`
- 环境变量可覆盖 Slurm 资源：`SLURM_TIME`、`SLURM_MEM`

如果 reference 已存在并且你想直接复用，可使用：

> `/data/person/wup/public/liusy_files/reference_genomes/hg38/resources/cnvkit/wes_reference.cnn`

---

## 阶段二：tumor CNV 分析（新增）

> 说明：pipeline 现在直接使用带基因注释的 target BED（4列），不再生成或校验原始 3 列 BED。


入口脚本：`03_submit_cnvkit_tumor_array.sh`

### 主链步骤

每个样本（一个 Slurm array task）按需执行：

1. tumor `targetcoverage`
2. tumor `antitargetcoverage`
3. `fix`（使用 `wes_reference.cnn`）
4. `segment`

### 可选增强步骤（开关控制）

- `segmetrics --ci`
- `call --filter ci`
- `genemetrics`
- `export seg`
- `export bed`
- `export vcf`
- `scatter`
- `diagram`

默认开启 filtering/genemetrics/export/plots，可通过 CLI 关闭。

---

## 阶段控制

`--stage` 支持：
- `coverage`
- `fix`
- `segment`
- `call`
- `reports`
- `all`

说明：
- `fix` 会自动依赖并复用/生成 coverage
- `segment` 会自动依赖 fix
- `call` 会自动依赖 segment
- `reports` 会自动依赖 call

---

## 输出结构（每个样本）

根目录：`<output_root>/samples/<sample_id>/`

- `coverage/`
  - `<sample>.targetcoverage.cnn`
  - `<sample>.antitargetcoverage.cnn`
- `cnr/`
  - `<sample>.cnr`
- `cns/`
  - `<sample>.seg.cns`
- `call/`
  - `<sample>.call.cns`（启用 filtering 时）
- `metrics/`
  - `<sample>.segmetrics.cns`
  - `<sample>.genemetrics.tsv`
- `export/`
  - `<sample>.seg`
  - `<sample>.bed`
  - `<sample>.vcf`
- `plots/`
  - `<sample>.scatter.pdf`
  - `<sample>.diagram.pdf`

---

## 最小运行示例

### 只跑主链到 segment

```bash
bash 03_submit_cnvkit_tumor_array.sh \
  --samples samples.tsv \
  --stage segment \
  --reference-cnn /data/person/wup/public/liusy_files/reference_genomes/hg38/resources/cnvkit/wes_reference.cnn \
  --max-parallel 8
```

### 全流程（含 filtering + reports）

```bash
bash 03_submit_cnvkit_tumor_array.sh \
  --samples samples.tsv \
  --stage all \
  --reference-cnn /data/person/wup/public/liusy_files/reference_genomes/hg38/resources/cnvkit/wes_reference.cnn \
  --max-parallel 8
```

---

## 推荐执行方案

1. 先跑：`--stage segment`（先得到稳定主链结果）
2. 再跑：`--stage all` 或 `--stage reports`（启用 filtering/call/reports）
3. 对比：
   - 原始分段：`cns/<sample>.seg.cns`
   - 过滤后推荐结果：`call/<sample>.call.cns`

---

## 关键配置（config/00_config.sh）

统一管理：
- `CNVKIT_REFERENCE_DEFAULT`
- `CNVKIT_ANNOTATED_TARGET_BED`
- `CNVKIT_ACCESS_BED`
- `REFERENCE`
- `CNVKIT_CMD`
- 输出/日志/tmp 默认目录
- 可选步骤开关默认值

---

## 阶段三：分组统计与自动出图（新增 downstream）

### 文件结构

```text
05_run_cnvkit_group_downstream.sh
scripts/
  cnvkit_group_downstream.py
```

### 每个脚本职责

- `05_run_cnvkit_group_downstream.sh`
  - 命令行入口；校验输入目录与 `group.tsv`
  - 统一组装参数并调用 Python 主脚本
  - 约定默认输出目录：`<cnvkit-root>/downstream/group_stats`
- `scripts/cnvkit_group_downstream.py`
  - 自动识别 `group.tsv` 是 2 组还是 3 组
  - 汇总 sample-level CNV burden 指标
  - 构建 gene-level CNV matrix（loss=-1, neutral=0, gain=1）
  - 进行 burden / gene-level / arm-level 频率统计比较（含 BH 校正 q 值）
  - 自动输出固定命名的 PDF 图与 TSV 结果表

### 输入说明

1. `--cnvkit-root`
   - CNVkit 主流程输出根目录，内部需包含：`samples/<sample_id>/...`
   - 脚本优先读取 `call/<sample_id>.call.cns`，若不存在则回退到 `cns/<sample_id>.seg.cns`
2. `--group-tsv`
   - 必须为 TAB 分隔，表头固定：
   - `sample_id<TAB>group`
   - 仅支持两组或三组

### 输出说明（固定文件名）

默认输出到：`<cnvkit-root>/downstream/group_stats/`

表格：
- `burden_metrics.tsv`
- `burden_group_comparison.tsv`
- `gene_cnv_matrix.tsv`
- `gene_frequency_by_group.tsv`
- `gene_frequency_comparison.tsv`
- `arm_level_cnv_matrix.tsv`
- `arm_level_frequency_comparison.tsv`
- `run_summary.tsv`

图片：
- `burden_boxplot.pdf`
- `gain_loss_burden_boxplot.pdf`
- `gene_frequency_barplot_topN.pdf`
- `gene_cnv_heatmap.pdf`
- `arm_level_heatmap.pdf`（arm 由每条染色体的中点启发式划分 p/q）

### 统计方法

- 两组：
  - burden：Wilcoxon rank-sum（`mannwhitneyu`）
  - gene/arm 频率：Fisher's exact
- 三组：
  - burden：Kruskal-Wallis
  - gene/arm 频率：卡方检验（`chi2_contingency`）
- 所有比较均同时输出原始 `p_value` 与 BH 校正 `q_value`

### 最小运行示例

```bash
bash 05_run_cnvkit_group_downstream.sh \
  --cnvkit-root /path/to/work/cnvkit_tumor \
  --group-tsv /path/to/group.tsv
```

可选参数示例：

```bash
bash 05_run_cnvkit_group_downstream.sh \
  --cnvkit-root /path/to/work/cnvkit_tumor \
  --group-tsv /path/to/group.tsv \
  --output-dir /path/to/group_stats \
  --top-n 40 \
  --gain-threshold 0.25 \
  --loss-threshold -0.25
```
