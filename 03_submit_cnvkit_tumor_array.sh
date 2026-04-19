#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/00_config.sh"

usage() {
  cat <<USAGE
Usage: $0 --samples <samples.tsv> [options]

Options:
  --samples <path>             samples.tsv path (required)
  --stage <name>               coverage|fix|segment|call|reports|all (default: all)
  --mode <name>                wes only (default: ${CNVKIT_MODE_DEFAULT})
  --reference-cnn <path>       pooled normal reference cnn (default: ${CNVKIT_REFERENCE_DEFAULT})
  --output-root <path>         output root (default: ${CNVKIT_OUTPUT_ROOT_DEFAULT})
  --log-dir <path>             log directory (default: ${CNVKIT_LOG_DIR_DEFAULT})
  --tmp-dir <path>             temp directory (default: ${CNVKIT_TMP_DIR_DEFAULT})
  --max-parallel <int>         Slurm array concurrency cap (default: ${CNVKIT_MAX_PARALLEL_DEFAULT})
  --threads <int>              reserved thread parameter (default: ${CNVKIT_THREADS_DEFAULT})
  --partition <name>           cpu1 or cpu2 (default: ${SLURM_PARTITION_DEFAULT})
  --disable-filtering          skip segmetrics --ci and call --filter ci
  --disable-genemetrics        skip genemetrics
  --disable-export             skip export seg/bed/vcf
  --disable-plots              skip scatter/diagram
  -h, --help                   show help
USAGE
}

samples_tsv=""
stage="all"
mode="${CNVKIT_MODE_DEFAULT}"
reference_cnn="${CNVKIT_REFERENCE_DEFAULT}"
output_root="${CNVKIT_OUTPUT_ROOT_DEFAULT}"
log_dir="${CNVKIT_LOG_DIR_DEFAULT}"
tmp_dir="${CNVKIT_TMP_DIR_DEFAULT}"
max_parallel="${CNVKIT_MAX_PARALLEL_DEFAULT}"
threads="${CNVKIT_THREADS_DEFAULT}"
slurm_partition="${SLURM_PARTITION:-${SLURM_PARTITION_DEFAULT}}"
slurm_time="${SLURM_TIME:-${SLURM_TIME_DEFAULT}}"
slurm_mem="${SLURM_MEM:-${SLURM_MEM_DEFAULT}}"
enable_filtering="${CNVKIT_ENABLE_FILTERING_DEFAULT}"
enable_genemetrics="${CNVKIT_ENABLE_GENEMETRICS_DEFAULT}"
enable_export="${CNVKIT_ENABLE_EXPORT_DEFAULT}"
enable_plots="${CNVKIT_ENABLE_PLOTS_DEFAULT}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) samples_tsv="$2"; shift 2 ;;
    --stage) stage="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --reference-cnn) reference_cnn="$2"; shift 2 ;;
    --output-root) output_root="$2"; shift 2 ;;
    --log-dir) log_dir="$2"; shift 2 ;;
    --tmp-dir) tmp_dir="$2"; shift 2 ;;
    --max-parallel) max_parallel="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --partition) slurm_partition="$2"; shift 2 ;;
    --disable-filtering) enable_filtering=0; shift ;;
    --disable-genemetrics) enable_genemetrics=0; shift ;;
    --disable-export) enable_export=0; shift ;;
    --disable-plots) enable_plots=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$samples_tsv" ]] || { echo "--samples is required" >&2; exit 1; }
[[ "$mode" == "wes" ]] || { echo "Only --mode wes is supported" >&2; exit 1; }
[[ -s "$reference_cnn" ]] || { echo "Missing reference.cnn: $reference_cnn" >&2; exit 1; }
[[ " ${SLURM_PARTITION_ALLOWED} " == *" ${slurm_partition} "* ]] || { echo "Unsupported partition: ${slurm_partition}" >&2; exit 1; }

case "$stage" in
  coverage|fix|segment|call|reports|all) ;;
  *) echo "Unsupported --stage: $stage" >&2; exit 1 ;;
esac

mkdir -p "$output_root" "$log_dir" "$tmp_dir" "$output_root/meta"
validated_samples="${output_root}/meta/samples.validated.tsv"
antitarget_bed="${output_root}/meta/antitarget.hg38.bed"
annotated_target_bed="${CNVKIT_ANNOTATED_TARGET_BED}"

[[ -s "$annotated_target_bed" ]] || { echo "Missing annotated target BED: $annotated_target_bed" >&2; exit 1; }

"${script_dir}/lib/validate_samples_tsv.sh" "$samples_tsv" "$validated_samples"
sample_count=$(( $(wc -l < "$validated_samples") - 1 ))
[[ "$sample_count" -gt 0 ]] || { echo "No samples in validated table" >&2; exit 1; }

array_range="0-$((sample_count - 1))%${max_parallel}"
common_sbatch_args=(
  --parsable
  --partition="$slurm_partition"
  --time="$slurm_time"
  --mem="$slurm_mem"
)

array_job_id=$(sbatch \
  "${common_sbatch_args[@]}" \
  --cpus-per-task="$threads" \
  --job-name=cnvkit_tumor \
  --output="${log_dir}/tumor_%A_%a.out" \
  --error="${log_dir}/tumor_%A_%a.err" \
  --array="$array_range" \
  --export=ALL,SCRIPT_DIR="$script_dir",VALIDATED_SAMPLES="$validated_samples",OUTPUT_ROOT="$output_root",TMP_DIR="$tmp_dir",REFERENCE_CNN="$reference_cnn",STAGE="$stage",ANNOTATED_TARGET_BED="$annotated_target_bed",ANTITARGET_BED="$antitarget_bed",THREADS="$threads",ENABLE_FILTERING="$enable_filtering",ENABLE_GENEMETRICS="$enable_genemetrics",ENABLE_EXPORT="$enable_export",ENABLE_PLOTS="$enable_plots" \
  "${script_dir}/run_cnvkit_tumor_array.sbatch")

echo "[INFO] Submitted tumor array job: ${array_job_id}"
echo "[INFO] stage=${stage}, output_root=${output_root}"
