#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/00_config.sh"
source "${script_dir}/lib/steps_cnvkit.sh"

normal_bam=""
coverage_dir=""
antitarget_bed=""
threads="${CNVKIT_THREADS_DEFAULT}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --normal-bam) normal_bam="$2"; shift 2 ;;
    --coverage-dir) coverage_dir="$2"; shift 2 ;;
    --antitarget-bed) antitarget_bed="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

: "${normal_bam:?--normal-bam is required}"
: "${coverage_dir:?--coverage-dir is required}"
: "${antitarget_bed:?--antitarget-bed is required}"

cnvkit_build_antitarget_bed "$CNVKIT_TARGET_BED" "$CNVKIT_ACCESS_BED" "$antitarget_bed"
cnvkit_run_normal_coverage "$normal_bam" "$CNVKIT_TARGET_BED" "$antitarget_bed" "$coverage_dir" "$threads"
