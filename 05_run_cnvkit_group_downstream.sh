#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<USAGE
Usage: $0 --cnvkit-root <dir> --group-tsv <group.tsv> [options]

Required:
  --cnvkit-root <dir>          CNVkit tumor output root (contains samples/<sample_id>/...)
  --group-tsv <path>           Group table, TAB-delimited with header: sample_id<tab>group

Optional:
  --output-dir <dir>           Output directory (default: <cnvkit-root>/downstream/group_stats)
  --top-n <int>                Top N genes to visualize (default: 30)
  --gain-threshold <float>     Gain threshold on log2 when cn column is unavailable (default: 0.2)
  --loss-threshold <float>     Loss threshold on log2 when cn column is unavailable (default: -0.2)
  --pathway-geneset <path>     Optional pathway geneset TSV: pathway<tab>gene
  -h, --help                   Show help
USAGE
}

cnvkit_root=""
group_tsv=""
output_dir=""
top_n="30"
gain_threshold="0.2"
loss_threshold="-0.2"
pathway_geneset=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cnvkit-root) cnvkit_root="$2"; shift 2 ;;
    --group-tsv) group_tsv="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --top-n) top_n="$2"; shift 2 ;;
    --gain-threshold) gain_threshold="$2"; shift 2 ;;
    --loss-threshold) loss_threshold="$2"; shift 2 ;;
    --pathway-geneset) pathway_geneset="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$cnvkit_root" ]] || { echo "--cnvkit-root is required" >&2; exit 1; }
[[ -n "$group_tsv" ]] || { echo "--group-tsv is required" >&2; exit 1; }
[[ -d "$cnvkit_root" ]] || { echo "CNVkit root dir not found: $cnvkit_root" >&2; exit 1; }
[[ -s "$group_tsv" ]] || { echo "group.tsv not found or empty: $group_tsv" >&2; exit 1; }

if [[ -z "$output_dir" ]]; then
  output_dir="${cnvkit_root%/}/downstream/group_stats"
fi
mkdir -p "$output_dir"

python3 "$script_dir/scripts/cnvkit_group_downstream.py" \
  --cnvkit-root "$cnvkit_root" \
  --group-tsv "$group_tsv" \
  --output-dir "$output_dir" \
  --top-n "$top_n" \
  --gain-threshold "$gain_threshold" \
  --loss-threshold "$loss_threshold" \
  --pathway-geneset "$pathway_geneset"

echo "[INFO] Done. Output: $output_dir"
