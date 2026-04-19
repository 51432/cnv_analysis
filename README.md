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
