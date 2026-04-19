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

## 阶段 2 输出结果说明（给新手）

下面这节专门解释阶段 2 每个样本目录里“文件是干什么的、先看什么、后看什么”。

示例目录：

```text
work/cnvkit_tumor/samples/TSDE014/
  call/
  cnr/
  cns/
  coverage/
  export/
  metrics/
  plots/
```

### 1）每个子目录的作用

- `coverage/`
  - 覆盖度中间结果（target + antitarget），主要用于 `fix` 生成 `.cnr`。
  - 一般不是人工解读 CNV 的第一选择，但用于排查流程是否跑通很有用。
- `cnr/`
  - bin-level（窗口级）log2 信号文件，分辨率高但噪声也更高。
  - 适合深度排查和自定义重分段，不是新手首选结果文件。
- `cns/`
  - 分段结果（segment-level），是“原始分段”主结果之一。
  - 适合看样本整体 CNV 结构，也常用于 burden / arm-level 统计。
- `call/`
  - 在分段基础上做了 calling/filtering 后的结果（如 `--filter ci`）。
  - 更适合做“保守版”事件列表展示，但可能比 `cns/` 更粗。
- `metrics/`
  - 统计摘要结果目录，常见有 `segmetrics` 和 `genemetrics`。
  - 尤其 `genemetrics.tsv` 很适合 gene-level 汇总分析。
- `export/`
  - 为其他软件或下游模块导出的交换格式（SEG/BED/VCF）。
  - 常用于跨工具对接，不是最先看的主结果。
- `plots/`
  - CNVkit 自动生成的可视化（scatter/diagram）。
  - 最适合人工质控和快速浏览样本 CNV 全貌。

### 2）典型文件说明（按常见使用场景）

- `coverage/*.targetcoverage.cnn`
  - 靶区覆盖度；`fix` 的输入之一。
- `coverage/*.antitargetcoverage.cnn`
  - 非靶区（背景）覆盖度；`fix` 的输入之一。
- `cnr/*.cnr`
  - 每个 bin 的归一化 log2 值，后续 `segment` 的直接输入。
- `cns/*.seg.cns`
  - 分段输出（未 calling 的主分段）；适合查看原始分段结构。
- `call/*.call.cns`
  - calling 后分段（可能更保守/更粗）；适合看“过滤后高置信事件”。
- `metrics/*.genemetrics.tsv`
  - gene 维度统计；做 gene-level CNV 统计/比较时优先参考。
- `export/*.seg`
  - SEG 交换格式，常给 IGV/GISTIC 等工具或下游流程使用。
- `export/*.bed`
  - BED 区间格式，便于区间类工具处理。
- `export/*.vcf`
  - VCF 表达的 CNV 事件，便于与变异流程或注释流程衔接。
- `plots/*.scatter.pdf`
  - 全基因组散点 + 分段图；最常用的人工 QC 图之一。
- `plots/*.diagram.pdf`
  - 染色体示意图（CNV diagram）；适合快速看大片段扩增/缺失模式。

### 3）“看哪类问题用哪类文件”快速对照

- 想看**原始分段结果**：
  - 优先 `cns/*.seg.cns`（分段细节更完整）。
- 想做 **gene-level 统计**：
  - 优先 `metrics/*.genemetrics.tsv`；
  - 无该文件时可用 `call/*.call.cns` 或 `cns/*.seg.cns` 映射到基因。
- 想做 **burden / arm-level 分析**：
  - 通常优先 `cns/*.seg.cns`（分段更细，长度统计更稳定）；
  - 若只关心高置信事件，也可用 `call/*.call.cns` 做更保守分析。
- 需要给**其他软件/下游工具**：
  - 使用 `export/*.seg`、`export/*.bed`、`export/*.vcf`。
- 做**人工质控（QC）**：
  - 优先看 `plots/*.scatter.pdf` 和 `plots/*.diagram.pdf`。

### 4）新手如何看结果（推荐顺序）

1. **先看 `plots/`**  
   先用 `scatter.pdf` / `diagram.pdf` 快速判断样本是否有明显 CNV 模式、是否存在异常噪声。
2. **再看 `cns/`**  
   读取 `*.seg.cns` 看原始分段是否合理（是否出现连续大片段 gain/loss、是否过度碎片化）。
3. **然后看 `call/` 和 `metrics/`**  
   用 `*.call.cns` 看过滤后的高置信事件；用 `*.genemetrics.tsv` 做基因层面解释与汇总。
4. **最后再看 `export/`**  
   `export/` 更偏“对接其他工具”的交换格式，不建议作为第一主结果入口。

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
   - burden 统计优先读取 `cns/<sample_id>.seg.cns`（分段更细）；若不存在则回退 `call/<sample_id>.call.cns`
   - gene matrix 优先读取 `metrics/<sample_id>.genemetrics.tsv`；若缺失则回退 `call/<sample_id>.call.cns`，再回退 `cns/<sample_id>.seg.cns`
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
- `data_source_by_sample.tsv`

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
