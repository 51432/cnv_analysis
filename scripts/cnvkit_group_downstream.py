#!/usr/bin/env python3
import argparse
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, kruskal, mannwhitneyu

BURDEN_METRICS = [
    "segment_count",
    "gain_segment_count",
    "loss_segment_count",
    "gain_total_length",
    "loss_total_length",
    "altered_genome_fraction",
]


def parse_args():
    p = argparse.ArgumentParser(description="CNVkit group downstream analysis + auto report")
    p.add_argument("--cnvkit-root", required=True)
    p.add_argument("--group-tsv", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--top-n", type=int, default=30)
    p.add_argument("--gain-threshold", type=float, default=0.2)
    p.add_argument("--loss-threshold", type=float, default=-0.2)
    p.add_argument("--pathway-geneset", default="", help="Optional pathway geneset TSV: pathway<TAB>gene")
    return p.parse_args()


def bh_adjust(pvalues):
    if len(pvalues) == 0:
        return np.array([])
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adj = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = min(prev, ranked[i] * n / rank)
        prev = val
        adj[i] = val
    out = np.empty(n)
    out[order] = adj
    return out


def read_group_table(group_tsv):
    g = pd.read_csv(group_tsv, sep="\t")
    if list(g.columns) != ["sample_id", "group"]:
        raise ValueError("group.tsv header must be exactly: sample_id<TAB>group")
    if g["sample_id"].duplicated().any():
        raise ValueError("Duplicated sample_id in group.tsv")
    g["group"] = g["group"].astype(str)
    groups = sorted(g["group"].unique().tolist())
    if len(groups) not in (2, 3):
        raise ValueError(f"Only 2 or 3 groups are supported, got {len(groups)}")
    return g, groups


def classify_state_from_row(row, gain_thr, loss_thr):
    cn = row.get("cn", np.nan)
    if not pd.isna(cn):
        if cn > 2:
            return 1
        if cn < 2:
            return -1
        return 0
    log2 = row.get("log2", 0.0)
    if log2 >= gain_thr:
        return 1
    if log2 <= loss_thr:
        return -1
    return 0


def split_genes(g):
    if pd.isna(g):
        return []
    out = []
    for x in str(g).replace(";", ",").split(","):
        x = x.strip()
        if x and x != ".":
            out.append(x)
    return out


def locate_sample_files(sample_dir: Path, sample_id: str):
    files = {
        "seg_cns": sample_dir / "cns" / f"{sample_id}.seg.cns",
        "call_cns": sample_dir / "call" / f"{sample_id}.call.cns",
        "genemetrics": sample_dir / "metrics" / f"{sample_id}.genemetrics.tsv",
        "cnr": sample_dir / "cnr" / f"{sample_id}.cnr",
        "scatter": sample_dir / "plots" / f"{sample_id}.scatter.pdf",
        "diagram": sample_dir / "plots" / f"{sample_id}.diagram.pdf",
    }
    return {k: v if v.exists() and v.stat().st_size > 0 else None for k, v in files.items()}


def read_cns(path: Path, gain_thr, loss_thr):
    df = pd.read_csv(path, sep="\t", comment="#")
    req = {"chromosome", "start", "end", "gene", "log2"}
    if not req.issubset(df.columns):
        raise ValueError(f"{path} missing required columns: {req}")
    df = df.copy()
    df["state"] = df.apply(classify_state_from_row, axis=1, gain_thr=gain_thr, loss_thr=loss_thr)
    return df


def build_burden(df):
    lengths = (df["end"] - df["start"]).astype(float)
    gain = df["state"] == 1
    loss = df["state"] == -1
    alt = gain | loss
    total = lengths.sum()
    return {
        "segment_count": int(len(df)),
        "gain_segment_count": int(gain.sum()),
        "loss_segment_count": int(loss.sum()),
        "gain_total_length": float(lengths[gain].sum()),
        "loss_total_length": float(lengths[loss].sum()),
        "altered_genome_fraction": float(lengths[alt].sum() / total) if total > 0 else math.nan,
    }


def infer_arm(row, chr_max):
    chrom = str(row["chromosome"])
    mid = (float(row["start"]) + float(row["end"])) / 2.0
    arm = "p" if mid <= chr_max.get(chrom, float(row["end"])) / 2.0 else "q"
    return f"{chrom}{arm}"


def gene_state_from_genemetrics(path: Path, gain_thr, loss_thr):
    gm = pd.read_csv(path, sep="\t", comment="#")
    if "gene" not in gm.columns:
        raise ValueError(f"{path} missing gene column")
    if "cn" not in gm.columns and "log2" not in gm.columns:
        raise ValueError(f"{path} requires cn or log2")
    out = {}
    for _, row in gm.iterrows():
        gene = str(row.get("gene", "")).strip()
        if not gene or gene == ".":
            continue
        state = classify_state_from_row(row, gain_thr, loss_thr)
        if state == 0:
            continue
        prev = out.get(gene, 0)
        if abs(state) >= abs(prev):
            out[gene] = state
    return out


def gene_state_from_cns(cns_df):
    out = {}
    for _, row in cns_df.iterrows():
        state = int(row["state"])
        if state == 0:
            continue
        for gene in split_genes(row["gene"]):
            prev = out.get(gene, 0)
            if abs(state) >= abs(prev):
                out[gene] = state
    return out


def compare_metric(burden_df, groups):
    rows = []
    for m in BURDEN_METRICS:
        vals = [burden_df.loc[burden_df.group == g, m].dropna().astype(float).values for g in groups]
        if len(groups) == 2:
            _, p = mannwhitneyu(vals[0], vals[1], alternative="two-sided")
            test = "Wilcoxon_rank_sum"
        else:
            _, p = kruskal(*vals)
            test = "Kruskal_Wallis"
        rows.append({"metric": m, "test": test, "p_value": p})
    out = pd.DataFrame(rows)
    out["q_value"] = bh_adjust(out["p_value"].values)
    return out


def build_state_matrix(state_by_sample, sample_order, feature_name):
    feats = sorted(set(f for d in state_by_sample.values() for f in d))
    mat = pd.DataFrame({feature_name: feats})
    for s in sample_order:
        d = state_by_sample.get(s, {})
        mat[s] = [int(d.get(f, 0)) for f in feats]
    return mat


def freq_and_stats_direction(matrix_df, sample_to_group, groups, feature_col, direction, feature_label):
    sample_cols = [c for c in matrix_df.columns if c != feature_col]
    freq_rows, stat_rows = [], []
    for _, row in matrix_df.iterrows():
        feat = row[feature_col]
        cont = []
        for g in groups:
            members = [s for s in sample_cols if sample_to_group.get(s) == g]
            vals = row[members].astype(int) if members else pd.Series(dtype=int)
            pos = int((vals == direction).sum())
            neg = int((vals != direction).sum())
            freq = float(pos / len(members)) if members else math.nan
            freq_rows.append({feature_label: feat, "direction": "gain" if direction == 1 else "loss", "group": g, "freq": freq, "n": len(members)})
            cont.append([pos, neg])
        cont = np.asarray(cont)
        if len(groups) == 2:
            _, p = fisher_exact(cont)
            test = "Fishers_exact"
        else:
            _, p, _, _ = chi2_contingency(cont)
            test = "Chi_square"
        stat_rows.append({feature_label: feat, "direction": "gain" if direction == 1 else "loss", "test": test, "p_value": p})
    freq_df = pd.DataFrame(freq_rows)
    stat_df = pd.DataFrame(stat_rows)
    if not stat_df.empty:
        stat_df["q_value"] = bh_adjust(stat_df["p_value"].values)
    return freq_df, stat_df


def safe_import_matplotlib():
    try:
        import matplotlib.pyplot as plt

        return plt
    except Exception as e:
        raise RuntimeError(f"Plotting requires matplotlib. Please install dependencies. Detail: {e}")


def plot_qc_overview(qc_df, out_path):
    plt = safe_import_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].hist(qc_df["segment_count"].dropna(), bins=20)
    axes[0].set_title("Segment count distribution")
    axes[1].hist(qc_df["altered_genome_fraction"].dropna(), bins=20)
    axes[1].set_title("Altered genome fraction distribution")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_boxes(df, metrics, out_path, title):
    plt = safe_import_matplotlib()
    groups = sorted(df["group"].unique())
    n = len(metrics)
    r = int(math.ceil(n / 2))
    fig, axes = plt.subplots(r, 2, figsize=(12, 4 * r))
    axes = np.array(axes).reshape(-1)
    for i, m in enumerate(metrics):
        vals = [df.loc[df.group == g, m].dropna().values for g in groups]
        axes[i].boxplot(vals, labels=groups)
        axes[i].set_title(m)
        axes[i].tick_params(axis="x", labelrotation=20)
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_heatmap(matrix_df, feature_col, sample_order, out_path, title, top_n=40):
    plt = safe_import_matplotlib()
    score = (matrix_df.drop(columns=[feature_col]) != 0).sum(axis=1)
    keep = matrix_df.loc[score.sort_values(ascending=False).index[:top_n], feature_col].tolist()
    if not keep:
        return
    arr = matrix_df.set_index(feature_col).loc[keep, sample_order]
    fig, ax = plt.subplots(figsize=(max(10, len(sample_order) * 0.35), max(5, len(keep) * 0.23)))
    im = ax.imshow(arr.values, aspect="auto", cmap="bwr", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(keep)))
    ax.set_yticklabels(keep, fontsize=8)
    ax.set_title(title)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("state")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(fig)


def plot_top_gene_freq(freq_df, direction, out_path, top_n=20):
    plt = safe_import_matplotlib()
    sub = freq_df[freq_df.direction == direction].copy()
    if sub.empty:
        return
    overall = sub.groupby("gene", as_index=False)["freq"].mean().sort_values("freq", ascending=False)
    top = overall.head(top_n)["gene"].tolist()
    plot_df = sub[sub.gene.isin(top)].copy()
    plot_df["gene"] = pd.Categorical(plot_df["gene"], categories=top, ordered=True)
    piv = plot_df.pivot(index="gene", columns="group", values="freq").fillna(0.0)
    ax = piv.plot(kind="bar", figsize=(12, 6))
    ax.set_ylabel(f"{direction} frequency")
    ax.set_title(f"Top {direction} genes by group")
    plt.xticks(rotation=75, ha="right")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def summarize_recurrent(gene_freq_df, top_n=20):
    rows = []
    for direction in ("gain", "loss"):
        sub = gene_freq_df[gene_freq_df.direction == direction]
        if sub.empty:
            continue
        agg = sub.groupby("gene", as_index=False)["freq"].mean().sort_values("freq", ascending=False).head(top_n)
        for _, r in agg.iterrows():
            rows.append({"gene": r["gene"], "direction": direction, "mean_group_frequency": r["freq"]})
    return pd.DataFrame(rows)


def build_report(outdir: Path, groups, group_counts, warnings_df, burden_stats, gene_stats, arm_stats, pathway_enabled):
    lines = []
    lines.append("# CNVkit Group Downstream Report")
    lines.append("")
    lines.append("## 分析概览")
    lines.append(f"- 分析日期：{pd.Timestamp.utcnow().strftime('%Y-%m-%d %H:%M UTC')}")
    lines.append(f"- 分组数：{len(groups)}（{', '.join(groups)}）")
    lines.append("")
    lines.append("## 输入与分组信息")
    lines.append("- 输入为现有 CNVkit 结果目录 + group.tsv。")
    lines.append("- 各组样本数：")
    for g in groups:
        lines.append(f"  - {g}: {group_counts.get(g, 0)}")
    lines.append("")
    lines.append("## QC 概览")
    lines.append("本节检查样本结果完整性与异常样本信号。")
    lines.append("![fig1](report_assets/fig1_qc_overview.svg)")
    if len(warnings_df) > 0:
        lines.append(f"- 共记录 warning {len(warnings_df)} 条，详见 warning_summary.tsv。")
    else:
        lines.append("- 未检测到 warning。")
    lines.append("")
    lines.append("## CNV burden 分组比较")
    lines.append("本节比较 segment 数量与基因组改变比例等 burden 指标。")
    lines.append("![fig2](report_assets/fig2_burden_boxplot.svg)")
    lines.append("![fig3](report_assets/fig3_gain_loss_burden_boxplot.svg)")
    if not burden_stats.empty:
        best = burden_stats.sort_values("p_value").iloc[0]
        lines.append(f"- 最显著指标：{best['metric']}（p={best['p_value']:.3g}, q={best['q_value']:.3g}）。")
    lines.append("")
    lines.append("## gene-level CNV 结果")
    lines.append("本节展示基因层面 CNV 矩阵和按组频率比较。")
    lines.append("![fig4](report_assets/fig4_gene_heatmap.svg)")
    lines.append("![fig5](report_assets/fig5_gene_frequency_topN_gain.svg)")
    lines.append("![fig6](report_assets/fig6_gene_frequency_topN_loss.svg)")
    if not gene_stats.empty:
        best = gene_stats.sort_values("p_value").iloc[0]
        lines.append(f"- 最显著基因/方向：{best['gene']} {best['direction']}（p={best['p_value']:.3g}, q={best['q_value']:.3g}）。")
    lines.append("")
    lines.append("## arm-level / broad CNA 结果")
    lines.append("本节概览染色体臂层面的 gain/loss 频率。")
    lines.append("![fig7](report_assets/fig7_arm_level_heatmap.svg)")
    if not arm_stats.empty:
        best = arm_stats.sort_values("p_value").iloc[0]
        lines.append(f"- 最显著染色体臂/方向：{best['arm']} {best['direction']}（p={best['p_value']:.3g}, q={best['q_value']:.3g}）。")
    lines.append("")
    lines.append("## pathway-level CNV 结果")
    if pathway_enabled:
        lines.append("已启用 pathway 基因集分析。")
        lines.append("![fig8](report_assets/fig8_pathway_summary.svg)")
    else:
        lines.append("未提供 pathway_geneset.tsv，本节跳过。")
    lines.append("")
    lines.append("## 结果文件说明")
    lines.append("- 主要表格：sample_qc_summary.tsv, sample_burden.tsv, burden_group_stats.tsv, gene_cnv_matrix.tsv, gene_frequency_by_group.tsv, gene_group_comparison.tsv, arm_level_frequency_by_group.tsv, broad_cna_summary.tsv, recurrent_cnv_summary.tsv。")
    lines.append("- 审计与告警：data_source_by_sample.tsv, warning_summary.tsv。")
    lines.append("- 可视化：report_assets/fig1-fig8.svg。")
    (outdir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def maybe_pathway_analysis(pathway_file, gene_matrix, group_df, outdir):
    if not pathway_file:
        return pd.DataFrame(), False
    p = Path(pathway_file)
    if not (p.exists() and p.stat().st_size > 0):
        return pd.DataFrame(), False
    pg = pd.read_csv(p, sep="\t")
    if list(pg.columns) != ["pathway", "gene"]:
        raise ValueError("pathway_geneset.tsv header must be pathway<TAB>gene")
    sample_order = group_df["sample_id"].tolist()
    gene_idx = gene_matrix.set_index("gene")
    rows = []
    for pathway, sub in pg.groupby("pathway"):
        genes = sorted(set(sub["gene"]))
        genes = [g for g in genes if g in gene_idx.index]
        if not genes:
            continue
        mat = gene_idx.loc[genes, sample_order]
        for _, gr in group_df.groupby("group"):
            sids = gr["sample_id"].tolist()
            subm = mat[sids]
            gain_frac = float((subm == 1).sum().sum() / subm.size) if subm.size else math.nan
            loss_frac = float((subm == -1).sum().sum() / subm.size) if subm.size else math.nan
            altered_frac = float((subm != 0).sum().sum() / subm.size) if subm.size else math.nan
            rows.append({"pathway": pathway, "group": gr["group"].iloc[0], "gain_fraction": gain_frac, "loss_fraction": loss_frac, "altered_fraction": altered_frac, "n_genes": len(genes), "n_samples": len(sids)})
    out = pd.DataFrame(rows)
    out.to_csv(outdir / "pathway_cnv_summary.tsv", sep="\t", index=False)

    if not out.empty:
        plt = safe_import_matplotlib()
        piv = out.pivot(index="pathway", columns="group", values="altered_fraction").fillna(0.0)
        ax = piv.plot(kind="bar", figsize=(12, 6))
        ax.set_ylabel("Altered fraction")
        ax.set_title("Pathway CNV altered fraction by group")
        plt.xticks(rotation=60, ha="right")
        plt.tight_layout()
        plt.savefig(outdir / "fig8_pathway_summary.svg")
        plt.close()
    return out, True


def main():
    args = parse_args()
    cnv_root = Path(args.cnvkit_root).resolve()
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    group_df, groups = read_group_table(args.group_tsv)
    sample_to_group = dict(zip(group_df["sample_id"], group_df["group"]))

    warning_rows = []
    for g, sub in group_df.groupby("group"):
        if len(sub) < 3:
            warning_rows.append({"sample_id": "*group*", "group": g, "warning_type": "small_group", "warning": f"Group {g} has only {len(sub)} samples"})

    burden_rows, qc_rows, source_rows = [], [], []
    gene_state_by_sample, arm_state_by_sample = {}, {}
    chr_lengths = {}

    for sample_id, group in group_df[["sample_id", "group"]].itertuples(index=False):
        sample_dir = cnv_root / "samples" / sample_id
        files = locate_sample_files(sample_dir, sample_id)

        burden_source = ""
        gene_source = ""
        qc_source = ""
        sample_warn = []

        burden_file = files["seg_cns"] or files["call_cns"]
        if burden_file is None:
            sample_warn.append("missing_seg_and_call")
            warning_rows.append({"sample_id": sample_id, "group": group, "warning_type": "missing_result", "warning": "No seg.cns or call.cns"})
            source_rows.append({"sample_id": sample_id, "group": group, "burden_source": "NA", "gene_source": "NA", "qc_source": "NA", "warnings": ";".join(sample_warn)})
            qc_rows.append({"sample_id": sample_id, "group": group, "has_cnv_result": 0, "has_cnr": int(files["cnr"] is not None), "has_plots": int(files["scatter"] is not None or files["diagram"] is not None), "segment_count": np.nan, "altered_genome_fraction": np.nan, "is_outlier_segment_count": np.nan})
            continue

        burden_source = "seg_cns" if files["seg_cns"] is not None else "call_cns"
        cns = read_cns(burden_file, args.gain_threshold, args.loss_threshold)
        burden = build_burden(cns)
        burden["sample_id"] = sample_id
        burden["group"] = group
        burden_rows.append(burden)

        chr_max = cns.groupby("chromosome")["end"].max().to_dict()
        for k, v in chr_max.items():
            chr_lengths[k] = max(chr_lengths.get(k, 0), float(v))

        gene_state = {}
        if files["genemetrics"] is not None:
            gene_state = gene_state_from_genemetrics(files["genemetrics"], args.gain_threshold, args.loss_threshold)
            gene_source = "genemetrics"
        elif files["call_cns"] is not None:
            call_df = read_cns(files["call_cns"], args.gain_threshold, args.loss_threshold)
            gene_state = gene_state_from_cns(call_df)
            gene_source = "call_cns"
            sample_warn.append("genemetrics_missing_fallback_call")
        else:
            gene_state = gene_state_from_cns(cns)
            gene_source = "seg_cns"
            sample_warn.append("genemetrics_call_missing_fallback_seg")
        gene_state_by_sample[sample_id] = gene_state

        arm_state = {}
        for _, row in cns.iterrows():
            s = int(row["state"])
            if s == 0:
                continue
            arm = infer_arm(row, chr_max)
            prev = arm_state.get(arm, 0)
            if abs(s) >= abs(prev):
                arm_state[arm] = s
        arm_state_by_sample[sample_id] = arm_state

        has_plots = int(files["scatter"] is not None or files["diagram"] is not None)
        qc_source = "cnr+plots" if files["cnr"] is not None and has_plots else "partial_qc"
        if files["cnr"] is None:
            sample_warn.append("cnr_missing")
        if not has_plots:
            sample_warn.append("plots_missing")

        source_rows.append({"sample_id": sample_id, "group": group, "burden_source": burden_source, "gene_source": gene_source, "qc_source": qc_source, "warnings": ";".join(sample_warn)})
        qc_rows.append({"sample_id": sample_id, "group": group, "has_cnv_result": 1, "has_cnr": int(files["cnr"] is not None), "has_plots": has_plots, "segment_count": burden["segment_count"], "altered_genome_fraction": burden["altered_genome_fraction"], "is_outlier_segment_count": 0})

    source_df = pd.DataFrame(source_rows).sort_values(["group", "sample_id"])
    source_df.to_csv(outdir / "data_source_by_sample.tsv", sep="\t", index=False)

    qc_df = pd.DataFrame(qc_rows).sort_values(["group", "sample_id"])
    if not qc_df["segment_count"].dropna().empty:
        s = qc_df["segment_count"].dropna()
        q1, q3 = s.quantile(0.25), s.quantile(0.75)
        iqr = q3 - q1
        low, high = max(5, q1 - 3 * iqr), q3 + 3 * iqr
        outlier_idx = qc_df[(qc_df["segment_count"] < low) | (qc_df["segment_count"] > high)].index
        qc_df.loc[outlier_idx, "is_outlier_segment_count"] = 1
        for idx in outlier_idx:
            row = qc_df.loc[idx]
            warning_rows.append({"sample_id": row["sample_id"], "group": row["group"], "warning_type": "segment_count_outlier", "warning": f"segment_count={row['segment_count']} outside [{low:.1f}, {high:.1f}]"})
    qc_df.to_csv(outdir / "sample_qc_summary.tsv", sep="\t", index=False)

    burden_df = pd.DataFrame(burden_rows).sort_values(["group", "sample_id"])
    burden_df.to_csv(outdir / "sample_burden.tsv", sep="\t", index=False)

    if burden_df.empty:
        warning_rows.append({"sample_id": "*global*", "group": "NA", "warning_type": "fatal", "warning": "No analyzable samples with burden source"})
        pd.DataFrame(warning_rows).to_csv(outdir / "warning_summary.tsv", sep="\t", index=False)
        return

    burden_stats = compare_metric(burden_df, groups)
    burden_stats.to_csv(outdir / "burden_group_stats.tsv", sep="\t", index=False)

    sample_order = group_df[group_df["sample_id"].isin(burden_df["sample_id"])].sort_values(["group", "sample_id"])["sample_id"].tolist()
    gene_matrix = build_state_matrix(gene_state_by_sample, sample_order, "gene")
    gene_matrix.to_csv(outdir / "gene_cnv_matrix.tsv", sep="\t", index=False)

    gene_freq_gain, gene_stat_gain = freq_and_stats_direction(gene_matrix, sample_to_group, groups, "gene", 1, "gene")
    gene_freq_loss, gene_stat_loss = freq_and_stats_direction(gene_matrix, sample_to_group, groups, "gene", -1, "gene")
    gene_freq = pd.concat([gene_freq_gain, gene_freq_loss], ignore_index=True)
    gene_stats = pd.concat([gene_stat_gain, gene_stat_loss], ignore_index=True)
    if not gene_stats.empty:
        gene_stats["q_value"] = bh_adjust(gene_stats["p_value"].values)
    gene_freq.to_csv(outdir / "gene_frequency_by_group.tsv", sep="\t", index=False)
    gene_stats.to_csv(outdir / "gene_group_comparison.tsv", sep="\t", index=False)

    arm_matrix = build_state_matrix(arm_state_by_sample, sample_order, "arm")
    arm_freq_gain, arm_stat_gain = freq_and_stats_direction(arm_matrix, sample_to_group, groups, "arm", 1, "arm")
    arm_freq_loss, arm_stat_loss = freq_and_stats_direction(arm_matrix, sample_to_group, groups, "arm", -1, "arm")
    arm_freq = pd.concat([arm_freq_gain, arm_freq_loss], ignore_index=True)
    arm_stats = pd.concat([arm_stat_gain, arm_stat_loss], ignore_index=True)
    if not arm_stats.empty:
        arm_stats["q_value"] = bh_adjust(arm_stats["p_value"].values)
    arm_out = arm_freq.merge(arm_stats, on=["arm", "direction"], how="left")
    arm_out.to_csv(outdir / "arm_level_frequency_by_group.tsv", sep="\t", index=False)

    broad_rows = []
    for _, row in burden_df.iterrows():
        sid = row["sample_id"]
        group = row["group"]
        sample_dir = cnv_root / "samples" / sid
        files = locate_sample_files(sample_dir, sid)
        cns_path = files["seg_cns"] or files["call_cns"]
        cns = read_cns(cns_path, args.gain_threshold, args.loss_threshold)
        gain_chr = set()
        loss_chr = set()
        for _, r in cns.iterrows():
            chrom = str(r["chromosome"])
            length = float(r["end"] - r["start"])
            chr_len = max(chr_lengths.get(chrom, 0.0), 1.0)
            if length / chr_len >= 0.5:
                if int(r["state"]) == 1:
                    gain_chr.add(chrom)
                elif int(r["state"]) == -1:
                    loss_chr.add(chrom)
        broad_rows.append({"sample_id": sid, "group": group, "broad_gain_chr_count": len(gain_chr), "broad_loss_chr_count": len(loss_chr), "broad_total_chr_count": len(gain_chr | loss_chr)})
    broad_df = pd.DataFrame(broad_rows)
    broad_df.to_csv(outdir / "broad_cna_summary.tsv", sep="\t", index=False)

    recurrent_df = summarize_recurrent(gene_freq, top_n=args.top_n)
    recurrent_df.to_csv(outdir / "recurrent_cnv_summary.tsv", sep="\t", index=False)

    pathway_df, pathway_enabled = maybe_pathway_analysis(args.pathway_geneset, gene_matrix, group_df[group_df["sample_id"].isin(sample_order)], outdir)

    # plotting
    plot_qc_overview(qc_df, outdir / "fig1_qc_overview.svg")
    plot_boxes(burden_df, BURDEN_METRICS, outdir / "fig2_burden_boxplot.svg", "CNV burden metrics")
    plot_boxes(burden_df, ["gain_segment_count", "loss_segment_count", "gain_total_length", "loss_total_length"], outdir / "fig3_gain_loss_burden_boxplot.svg", "Gain/Loss burden")
    plot_heatmap(gene_matrix, "gene", sample_order, outdir / "fig4_gene_heatmap.svg", "Gene-level CNV heatmap", top_n=args.top_n)
    plot_top_gene_freq(gene_freq, "gain", outdir / "fig5_gene_frequency_topN_gain.svg", top_n=args.top_n)
    plot_top_gene_freq(gene_freq, "loss", outdir / "fig6_gene_frequency_topN_loss.svg", top_n=args.top_n)
    plot_heatmap(arm_matrix, "arm", sample_order, outdir / "fig7_arm_level_heatmap.svg", "Arm-level CNV heatmap", top_n=args.top_n)

    warn_df = pd.DataFrame(warning_rows)
    if warn_df.empty:
        warn_df = pd.DataFrame(columns=["sample_id", "group", "warning_type", "warning"])
    warn_df.to_csv(outdir / "warning_summary.tsv", sep="\t", index=False)

    group_counts = group_df.groupby("group").size().to_dict()
    build_report(outdir, groups, group_counts, warn_df, burden_stats, gene_stats, arm_stats, pathway_enabled)

    assets = outdir / "report_assets"
    assets.mkdir(exist_ok=True)
    for fig in [
        "fig1_qc_overview.svg",
        "fig2_burden_boxplot.svg",
        "fig3_gain_loss_burden_boxplot.svg",
        "fig4_gene_heatmap.svg",
        "fig5_gene_frequency_topN_gain.svg",
        "fig6_gene_frequency_topN_loss.svg",
        "fig7_arm_level_heatmap.svg",
        "fig8_pathway_summary.svg",
    ]:
        src = outdir / fig
        if src.exists():
            shutil.copy2(src, assets / fig)


if __name__ == "__main__":
    main()
