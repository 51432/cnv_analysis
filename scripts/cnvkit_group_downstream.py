#!/usr/bin/env python3
import argparse
import math
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
    parser = argparse.ArgumentParser(description="Group-level downstream analysis for CNVkit WES results")
    parser.add_argument("--cnvkit-root", required=True)
    parser.add_argument("--group-tsv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--top-n", type=int, default=30)
    parser.add_argument("--gain-threshold", type=float, default=0.2)
    parser.add_argument("--loss-threshold", type=float, default=-0.2)
    return parser.parse_args()


def bh_adjust(pvalues):
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
    result = np.empty(n)
    result[order] = adj
    return result


def read_group_table(group_tsv):
    group_df = pd.read_csv(group_tsv, sep="\t")
    if list(group_df.columns) != ["sample_id", "group"]:
        raise ValueError("group.tsv header must be exactly: sample_id<TAB>group")
    if group_df["sample_id"].duplicated().any():
        raise ValueError("Duplicated sample_id in group.tsv")
    group_df["group"] = group_df["group"].astype(str)
    groups = sorted(group_df["group"].unique().tolist())
    if len(groups) not in (2, 3):
        raise ValueError(f"Only 2 or 3 groups are supported, got {len(groups)} groups")
    return group_df, groups


def locate_sample_files(cnvkit_root: Path, sample_id: str):
    sample_dir = cnvkit_root / "samples" / sample_id
    call_file = sample_dir / "call" / f"{sample_id}.call.cns"
    seg_file = sample_dir / "cns" / f"{sample_id}.seg.cns"
    gene_file = sample_dir / "metrics" / f"{sample_id}.genemetrics.tsv"
    return {
        "call_cns": call_file if call_file.exists() and call_file.stat().st_size > 0 else None,
        "seg_cns": seg_file if seg_file.exists() and seg_file.stat().st_size > 0 else None,
        "genemetrics": gene_file if gene_file.exists() and gene_file.stat().st_size > 0 else None,
    }


def classify_state(row, gain_thr, loss_thr):
    if "cn" in row and not pd.isna(row["cn"]):
        cn = row["cn"]
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


def split_genes(gene_field):
    if pd.isna(gene_field):
        return []
    genes = []
    for item in str(gene_field).replace(";", ",").split(","):
        g = item.strip()
        if g and g != ".":
            genes.append(g)
    return genes


def classify_state_from_dict(row_dict, gain_thr, loss_thr):
    if "cn" in row_dict and not pd.isna(row_dict["cn"]):
        cn = row_dict["cn"]
        if cn > 2:
            return 1
        if cn < 2:
            return -1
        return 0
    log2 = row_dict.get("log2", 0.0)
    if log2 >= gain_thr:
        return 1
    if log2 <= loss_thr:
        return -1
    return 0


def build_sample_metrics(cnv_df):
    lengths = (cnv_df["end"] - cnv_df["start"]).astype(float)
    gain = cnv_df["state"] == 1
    loss = cnv_df["state"] == -1
    altered_len = lengths[gain | loss].sum()
    total_len = lengths.sum()
    return {
        "segment_count": int(len(cnv_df)),
        "gain_segment_count": int(gain.sum()),
        "loss_segment_count": int(loss.sum()),
        "gain_total_length": float(lengths[gain].sum()),
        "loss_total_length": float(lengths[loss].sum()),
        "altered_genome_fraction": float(altered_len / total_len) if total_len > 0 else math.nan,
    }


def read_cns_with_state(cns_file: Path, gain_thr: float, loss_thr: float):
    cns = pd.read_csv(cns_file, sep="\t", comment="#")
    required_cols = {"chromosome", "start", "end", "gene", "log2"}
    if not required_cols.issubset(cns.columns):
        raise ValueError(f"{cns_file} missing required columns: {required_cols}")
    cns = cns.copy()
    cns["state"] = cns.apply(classify_state, axis=1, gain_thr=gain_thr, loss_thr=loss_thr)
    return cns


def gene_state_from_genemetrics(gene_file: Path, gain_thr: float, loss_thr: float):
    gm = pd.read_csv(gene_file, sep="\t", comment="#")
    if "gene" not in gm.columns:
        raise ValueError(f"{gene_file} missing required column: gene")
    if "cn" not in gm.columns and "log2" not in gm.columns:
        raise ValueError(f"{gene_file} must contain cn or log2 column")
    gene_state = {}
    for _, row in gm.iterrows():
        gene_name = row.get("gene")
        if pd.isna(gene_name) or str(gene_name).strip() in ("", "."):
            continue
        state = classify_state_from_dict(row.to_dict(), gain_thr=gain_thr, loss_thr=loss_thr)
        if state == 0:
            continue
        prev = gene_state.get(str(gene_name), 0)
        if abs(state) >= abs(prev):
            gene_state[str(gene_name)] = state
    return gene_state


def gene_state_from_cns(cns_df):
    gene_state = {}
    for _, row in cns_df.iterrows():
        state = int(row["state"])
        if state == 0:
            continue
        for gene in split_genes(row["gene"]):
            prev = gene_state.get(gene, 0)
            if abs(state) >= abs(prev):
                gene_state[gene] = state
    return gene_state


def infer_arm(row, chr_max_end):
    chrom = str(row["chromosome"])
    midpoint = (float(row["start"]) + float(row["end"])) / 2.0
    chrom_max = chr_max_end.get(chrom, float(row["end"]))
    arm = "p" if midpoint <= chrom_max / 2.0 else "q"
    return f"{chrom}{arm}"


def summarize_burden_stats(burden_df, groups):
    rows = []
    for metric in BURDEN_METRICS:
        series_by_group = [
            burden_df.loc[burden_df["group"] == g, metric].dropna().astype(float).values
            for g in groups
        ]
        if len(groups) == 2:
            _, pval = mannwhitneyu(series_by_group[0], series_by_group[1], alternative="two-sided")
            test_name = "Wilcoxon_rank_sum"
        else:
            _, pval = kruskal(*series_by_group)
            test_name = "Kruskal_Wallis"
        rows.append({"metric": metric, "test": test_name, "p_value": pval})
    out = pd.DataFrame(rows)
    out["q_value"] = bh_adjust(out["p_value"].values)
    return out.sort_values("p_value")


def compare_binary_feature(matrix_df, group_map, groups, feature_col="feature"):
    rows = []
    sample_cols = [c for c in matrix_df.columns if c not in [feature_col]]
    for _, rec in matrix_df.iterrows():
        feature = rec[feature_col]
        data = []
        for g in groups:
            members = [s for s in sample_cols if group_map.get(s) == g]
            vals = rec[members].fillna(0).astype(int)
            altered = int((vals != 0).sum())
            unaltered = int((vals == 0).sum())
            data.append([altered, unaltered])
        data_arr = np.asarray(data)
        if len(groups) == 2:
            _, pval = fisher_exact(data_arr)
            test_name = "Fishers_exact"
        else:
            _, pval, _, _ = chi2_contingency(data_arr)
            test_name = "Chi_square"
        rows.append({"feature": feature, "test": test_name, "p_value": pval})
    out = pd.DataFrame(rows)
    out["q_value"] = bh_adjust(out["p_value"].values)
    return out.sort_values("p_value")


def build_frequency_table(matrix_df, group_map, groups, feature_col="feature"):
    records = []
    sample_cols = [c for c in matrix_df.columns if c != feature_col]
    for _, row in matrix_df.iterrows():
        for g in groups:
            members = [s for s in sample_cols if group_map.get(s) == g]
            values = row[members].fillna(0).astype(int)
            freq = float((values != 0).mean()) if len(values) > 0 else math.nan
            records.append({feature_col: row[feature_col], "group": g, "frequency": freq, "n": len(values)})
    return pd.DataFrame(records)


def draw_burden_boxplots(burden_df, outdir: Path):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    axes = axes.ravel()
    for i, metric in enumerate(BURDEN_METRICS):
        data = [
            burden_df.loc[burden_df["group"] == g, metric].dropna().values
            for g in sorted(burden_df["group"].unique())
        ]
        labels = sorted(burden_df["group"].unique())
        axes[i].boxplot(data, labels=labels)
        axes[i].set_title(metric)
        axes[i].tick_params(axis="x", labelrotation=30)
    fig.tight_layout()
    fig.savefig(outdir / "burden_boxplot.pdf")
    plt.close(fig)

    gl_metrics = ["gain_segment_count", "loss_segment_count", "gain_total_length", "loss_total_length"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.ravel()
    for i, metric in enumerate(gl_metrics):
        data = [
            burden_df.loc[burden_df["group"] == g, metric].dropna().values
            for g in sorted(burden_df["group"].unique())
        ]
        labels = sorted(burden_df["group"].unique())
        axes[i].boxplot(data, labels=labels)
        axes[i].set_title(metric)
        axes[i].tick_params(axis="x", labelrotation=30)
    fig.tight_layout()
    fig.savefig(outdir / "gain_loss_burden_boxplot.pdf")
    plt.close(fig)


def draw_gene_frequency_barplot(freq_df, top_genes, outdir: Path):
    import matplotlib.pyplot as plt

    plot_df = freq_df[freq_df["gene"].isin(top_genes)].copy()
    plot_df["gene"] = pd.Categorical(plot_df["gene"], categories=top_genes, ordered=True)
    pivot = plot_df.pivot(index="gene", columns="group", values="frequency").fillna(0.0)
    ax = pivot.plot(kind="bar", figsize=(12, 6))
    ax.set_ylabel("Altered frequency")
    ax.set_xlabel("Gene")
    ax.set_title("Top altered genes by group")
    plt.xticks(rotation=75, ha="right")
    plt.tight_layout()
    plt.savefig(outdir / "gene_frequency_barplot_topN.pdf")
    plt.close()


def draw_gene_heatmap(gene_matrix, sample_order, top_genes, outdir: Path):
    import matplotlib.pyplot as plt

    sub = gene_matrix.set_index("gene").loc[top_genes, sample_order]
    fig, ax = plt.subplots(figsize=(max(10, len(sample_order) * 0.35), max(6, len(top_genes) * 0.22)))
    im = ax.imshow(sub.values, aspect="auto", cmap="bwr", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(top_genes)))
    ax.set_yticklabels(top_genes, fontsize=8)
    ax.set_title("Gene-level CNV matrix (loss=-1, neutral=0, gain=1)")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("CNV state")
    plt.tight_layout()
    plt.savefig(outdir / "gene_cnv_heatmap.pdf")
    plt.close(fig)


def draw_arm_heatmap(arm_matrix, sample_order, outdir: Path, top_n=40):
    import matplotlib.pyplot as plt

    arm_score = (arm_matrix.drop(columns=["arm"]) != 0).sum(axis=1)
    keep_arms = arm_matrix.loc[arm_score.sort_values(ascending=False).index[:top_n], "arm"].tolist()
    if not keep_arms:
        return
    sub = arm_matrix.set_index("arm").loc[keep_arms, sample_order]
    fig, ax = plt.subplots(figsize=(max(10, len(sample_order) * 0.35), max(6, len(keep_arms) * 0.22)))
    im = ax.imshow(sub.values, aspect="auto", cmap="bwr", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=90, fontsize=8)
    ax.set_yticks(np.arange(len(keep_arms)))
    ax.set_yticklabels(keep_arms, fontsize=8)
    ax.set_title("Arm-level CNV matrix (heuristic p/q split)")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("CNV state")
    plt.tight_layout()
    plt.savefig(outdir / "arm_level_heatmap.pdf")
    plt.close(fig)


def main():
    args = parse_args()
    cnvkit_root = Path(args.cnvkit_root).resolve()
    outdir = Path(args.output_dir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    group_df, groups = read_group_table(args.group_tsv)
    group_map = dict(zip(group_df["sample_id"], group_df["group"]))

    burden_rows = []
    gene_state_by_sample = {}
    arm_state_by_sample = {}
    missing = []
    source_rows = []

    for sample_id in group_df["sample_id"]:
        sample_files = locate_sample_files(cnvkit_root, sample_id)

        burden_cns_file = sample_files["seg_cns"] or sample_files["call_cns"]
        if burden_cns_file is None:
            missing.append(sample_id)
            continue
        cns = read_cns_with_state(burden_cns_file, gain_thr=args.gain_threshold, loss_thr=args.loss_threshold)

        metrics = build_sample_metrics(cns)
        metrics["sample_id"] = sample_id
        metrics["group"] = group_map[sample_id]
        burden_rows.append(metrics)

        gene_source = ""
        if sample_files["genemetrics"] is not None:
            gene_state = gene_state_from_genemetrics(
                sample_files["genemetrics"], gain_thr=args.gain_threshold, loss_thr=args.loss_threshold
            )
            gene_source = "genemetrics"
        elif sample_files["call_cns"] is not None:
            call_cns = read_cns_with_state(sample_files["call_cns"], gain_thr=args.gain_threshold, loss_thr=args.loss_threshold)
            gene_state = gene_state_from_cns(call_cns)
            gene_source = "call_cns"
        else:
            gene_state = gene_state_from_cns(cns)
            gene_source = "seg_cns"
        gene_state_by_sample[sample_id] = gene_state

        chr_max_end = cns.groupby("chromosome")["end"].max().to_dict()
        arm_state = {}
        for _, row in cns.iterrows():
            state = int(row["state"])
            if state == 0:
                continue
            arm = infer_arm(row, chr_max_end)
            prev = arm_state.get(arm, 0)
            if abs(state) >= abs(prev):
                arm_state[arm] = state
        arm_state_by_sample[sample_id] = arm_state
        source_rows.append(
            {
                "sample_id": sample_id,
                "group": group_map[sample_id],
                "burden_source": "seg_cns" if sample_files["seg_cns"] is not None else "call_cns",
                "gene_source": gene_source,
            }
        )

    if missing:
        raise FileNotFoundError(f"Missing CNV files for samples: {', '.join(missing)}")

    burden_df = pd.DataFrame(burden_rows).sort_values(["group", "sample_id"])
    burden_df.to_csv(outdir / "burden_metrics.tsv", sep="\t", index=False)

    burden_cmp = summarize_burden_stats(burden_df, groups)
    burden_cmp.to_csv(outdir / "burden_group_comparison.tsv", sep="\t", index=False)

    all_genes = sorted(set(g for sample_data in gene_state_by_sample.values() for g in sample_data))
    gene_matrix = pd.DataFrame({"gene": all_genes})
    for sid in group_df["sample_id"]:
        state_map = gene_state_by_sample.get(sid, {})
        gene_matrix[sid] = [state_map.get(g, 0) for g in all_genes]
    gene_matrix.to_csv(outdir / "gene_cnv_matrix.tsv", sep="\t", index=False)

    gene_freq = build_frequency_table(gene_matrix.rename(columns={"gene": "feature"}), group_map, groups, feature_col="feature")
    gene_freq = gene_freq.rename(columns={"feature": "gene"})
    gene_freq.to_csv(outdir / "gene_frequency_by_group.tsv", sep="\t", index=False)

    gene_cmp = compare_binary_feature(gene_matrix.rename(columns={"gene": "feature"}), group_map, groups, feature_col="feature")
    gene_cmp = gene_cmp.rename(columns={"feature": "gene"})
    gene_cmp.to_csv(outdir / "gene_frequency_comparison.tsv", sep="\t", index=False)

    all_arms = sorted(set(a for sample_data in arm_state_by_sample.values() for a in sample_data))
    arm_matrix = pd.DataFrame({"arm": all_arms})
    for sid in group_df["sample_id"]:
        state_map = arm_state_by_sample.get(sid, {})
        arm_matrix[sid] = [state_map.get(a, 0) for a in all_arms]
    arm_matrix.to_csv(outdir / "arm_level_cnv_matrix.tsv", sep="\t", index=False)

    arm_cmp = compare_binary_feature(arm_matrix.rename(columns={"arm": "feature"}), group_map, groups, feature_col="feature")
    arm_cmp = arm_cmp.rename(columns={"feature": "arm"})
    arm_cmp.to_csv(outdir / "arm_level_frequency_comparison.tsv", sep="\t", index=False)

    sample_order = group_df.sort_values(["group", "sample_id"])["sample_id"].tolist()

    draw_burden_boxplots(burden_df, outdir)

    gene_total_freq = (gene_matrix.drop(columns=["gene"]) != 0).mean(axis=1)
    top_genes = gene_matrix.loc[gene_total_freq.sort_values(ascending=False).index[: args.top_n], "gene"].tolist()
    if top_genes:
        draw_gene_frequency_barplot(gene_freq, top_genes, outdir)
        draw_gene_heatmap(gene_matrix, sample_order, top_genes, outdir)

    if len(all_arms) > 0:
        draw_arm_heatmap(arm_matrix, sample_order, outdir, top_n=min(args.top_n, len(all_arms)))

    source_df = pd.DataFrame(source_rows).sort_values(["group", "sample_id"])
    source_df.to_csv(outdir / "data_source_by_sample.tsv", sep="\t", index=False)

    summary = {
        "group_count": len(groups),
        "groups": ",".join(groups),
        "sample_count": len(group_df),
        "top_n": args.top_n,
    }
    pd.DataFrame([summary]).to_csv(outdir / "run_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
