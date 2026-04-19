#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/00_config.sh"

usage() {
  cat <<USAGE
Usage: $0 --samples <samples.tsv> [options]

Options:
  --samples <path>             samples.tsv path (required, TAB-delimited)
  --stage <name>               build-reference | run-cnv (default: ${CNVKIT_STAGE_DEFAULT})
  --mode <name>                wes only (default: ${CNVKIT_MODE_DEFAULT})
  --workdir <path>             work directory (default: ./work/cnvkit)
  --reference-out <path>       pooled reference output path (default: ${CNVKIT_REFERENCE_DEFAULT})
  --threads <int>              CNVkit threads for each array task (default: ${CNVKIT_THREADS_DEFAULT})
  --max-parallel <int>         Slurm array concurrency cap (default: ${CNVKIT_MAX_PARALLEL_DEFAULT})
  --partition <name>           Slurm partition, only cpu1 or cpu2 (default: ${SLURM_PARTITION_DEFAULT})
  --overwrite-reference        overwrite existing pooled reference

Environment overrides:
  SLURM_PARTITION              optional override (allowed: cpu1/cpu2)
  SLURM_TIME / SLURM_MEM
  -h, --help                   show help
USAGE
}

samples_tsv=""
stage="${CNVKIT_STAGE_DEFAULT}"
mode="${CNVKIT_MODE_DEFAULT}"
workdir="${PWD}/work/cnvkit"
reference_out="${CNVKIT_REFERENCE_DEFAULT}"
threads="${CNVKIT_THREADS_DEFAULT}"
max_parallel="${CNVKIT_MAX_PARALLEL_DEFAULT}"
overwrite_reference=0
slurm_partition="${SLURM_PARTITION:-${SLURM_PARTITION_DEFAULT}}"
slurm_time="${SLURM_TIME:-${SLURM_TIME_DEFAULT}}"
slurm_mem="${SLURM_MEM:-${SLURM_MEM_DEFAULT}}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) samples_tsv="$2"; shift 2 ;;
    --stage) stage="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --workdir) workdir="$2"; shift 2 ;;
    --reference-out) reference_out="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --max-parallel) max_parallel="$2"; shift 2 ;;
    --partition) slurm_partition="$2"; shift 2 ;;
    --overwrite-reference) overwrite_reference=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$samples_tsv" ]]; then
  echo "--samples is required" >&2
  usage
  exit 1
fi

if [[ "$mode" != "wes" ]]; then
  echo "Only --mode wes is supported now." >&2
  exit 1
fi

if [[ " ${SLURM_PARTITION_ALLOWED} " != *" ${slurm_partition} "* ]]; then
  echo "Unsupported partition: ${slurm_partition}. Allowed: ${SLURM_PARTITION_ALLOWED}" >&2
  exit 1
fi

case "$stage" in
  build-reference)
    ;;
  run-cnv)
    echo "Stage run-cnv is reserved for phase 2 and not implemented yet." >&2
    exit 2
    ;;
  *)
    echo "Unsupported --stage: $stage" >&2
    exit 1
    ;;
esac

mkdir -p "$workdir" "$workdir/logs" "$workdir/meta" "$workdir/coverage"
validated_samples="${workdir}/meta/samples.validated.tsv"
normal_list="${workdir}/meta/normal_bams.unique.list"
antitarget_bed="${workdir}/meta/antitarget.hg38.bed"

"${script_dir}/lib/validate_samples_tsv.sh" "$samples_tsv" "$validated_samples"

awk -F'\t' 'NR>1{if(!seen[$3]++){print $3}}' "$validated_samples" > "$normal_list"
normal_count="$(wc -l < "$normal_list")"
if [[ "$normal_count" -eq 0 ]]; then
  echo "No unique normal_bam found after validation." >&2
  exit 1
fi

if [[ -s "$reference_out" && "$overwrite_reference" -eq 0 ]]; then
  echo "[INFO] Reference already exists, skip submission: $reference_out"
  echo "[INFO] Use --overwrite-reference to rebuild."
  exit 0
fi

if [[ "$overwrite_reference" -eq 1 ]]; then
  rm -f "$reference_out"
fi

array_range="0-$((normal_count - 1))%${max_parallel}"

common_sbatch_args=(
  --parsable
  --partition="$slurm_partition"
  --time="$slurm_time"
  --mem="$slurm_mem"
)

array_job_id=$(sbatch \
  "${common_sbatch_args[@]}" \
  --cpus-per-task="$threads" \
  --job-name=cnvkit_ref_cov \
  --output="${workdir}/logs/cov_%A_%a.out" \
  --error="${workdir}/logs/cov_%A_%a.err" \
  --array="$array_range" \
  --export=ALL,SCRIPT_DIR="$script_dir",NORMAL_LIST="$normal_list",COVERAGE_DIR="${workdir}/coverage",ANTITARGET_BED="$antitarget_bed",THREADS="$threads" \
  "${script_dir}/run_cnvkit_reference_array.sbatch")

echo "[INFO] Submitted coverage array job: ${array_job_id}"

build_job_id=$(sbatch \
  "${common_sbatch_args[@]}" \
  --cpus-per-task=1 \
  --dependency="afterok:${array_job_id}" \
  --kill-on-invalid-dep=yes \
  --job-name=cnvkit_ref_build \
  --output="${workdir}/logs/build_%j.out" \
  --error="${workdir}/logs/build_%j.err" \
  --export=ALL,SCRIPT_DIR="$script_dir",COVERAGE_DIR="${workdir}/coverage",REFERENCE_OUT="$reference_out",ANTITARGET_BED="$antitarget_bed",OVERWRITE_REFERENCE="$overwrite_reference" \
  "${script_dir}/02_build_cnvkit_reference.sh")

echo "[INFO] Submitted reference build job: ${build_job_id}"
echo "[INFO] Reference target path: ${reference_out}"
echo "[INFO] If build job shows DependencyNeverSatisfied, check array logs: ${workdir}/logs/cov_<array_jobid>_<taskid>.err"
echo "[INFO] Suggested check: sacct -j ${array_job_id} --format=JobID,State,ExitCode"
