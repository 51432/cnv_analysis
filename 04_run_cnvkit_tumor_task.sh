#!/usr/bin/env bash
set -euo pipefail

script_dir="${SCRIPT_DIR:-}"
if [[ -z "$script_dir" ]]; then
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${script_dir}/config/00_config.sh"
source "${script_dir}/lib/steps_cnvkit.sh"

sample_id=""
tumor_bam=""
stage="all"
output_root=""
reference_cnn=""
target_bed3=""
antitarget_bed=""
antitarget_bed3=""
threads="1"
enable_filtering=1
enable_genemetrics=1
enable_export=1
enable_plots=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) sample_id="$2"; shift 2 ;;
    --tumor-bam) tumor_bam="$2"; shift 2 ;;
    --stage) stage="$2"; shift 2 ;;
    --output-root) output_root="$2"; shift 2 ;;
    --reference-cnn) reference_cnn="$2"; shift 2 ;;
    --target-bed3) target_bed3="$2"; shift 2 ;;
    --antitarget-bed) antitarget_bed="$2"; shift 2 ;;
    --antitarget-bed3) antitarget_bed3="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --enable-filtering) enable_filtering="$2"; shift 2 ;;
    --enable-genemetrics) enable_genemetrics="$2"; shift 2 ;;
    --enable-export) enable_export="$2"; shift 2 ;;
    --enable-plots) enable_plots="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

[[ -n "$sample_id" && -n "$tumor_bam" ]] || { echo "Missing sample info" >&2; exit 1; }

sample_root="${output_root}/samples/${sample_id}"
coverage_dir="${sample_root}/coverage"
cnr_dir="${sample_root}/cnr"
cns_dir="${sample_root}/cns"
call_dir="${sample_root}/call"
metrics_dir="${sample_root}/metrics"
export_dir="${sample_root}/export"
plots_dir="${sample_root}/plots"
mkdir -p "$coverage_dir" "$cnr_dir" "$cns_dir" "$call_dir" "$metrics_dir" "$export_dir" "$plots_dir"

cnvkit_prepare_bed3 "$CNVKIT_TARGET_BED" "$target_bed3"
cnvkit_build_antitarget_bed "$target_bed3" "$CNVKIT_ACCESS_BED" "$antitarget_bed"
cnvkit_prepare_bed3 "$antitarget_bed" "$antitarget_bed3"

need_fix=0
need_segment=0
need_call=0
need_reports=0
case "$stage" in
  coverage) ;;
  fix) need_fix=1 ;;
  segment) need_fix=1; need_segment=1 ;;
  call) need_fix=1; need_segment=1; need_call=1 ;;
  reports) need_fix=1; need_segment=1; need_call=1; need_reports=1 ;;
  all) need_fix=1; need_segment=1; need_call=1; need_reports=1 ;;
  *) echo "Unsupported stage: $stage" >&2; exit 1 ;;
esac

cnvkit_run_tumor_coverage "$sample_id" "$tumor_bam" "$target_bed3" "$antitarget_bed3" "$coverage_dir"

if [[ "$need_fix" -eq 1 ]]; then
  cnvkit_run_fix "$sample_id" "$coverage_dir" "$reference_cnn" "$cnr_dir"
fi

if [[ "$need_segment" -eq 1 ]]; then
  cnvkit_run_segment "$sample_id" "$cnr_dir" "$cns_dir"
fi

report_cns="${cns_dir}/${sample_id}.seg.cns"
if [[ "$need_call" -eq 1 ]]; then
  if [[ "$enable_filtering" -eq 1 ]]; then
    cnvkit_run_segmetrics_ci "$sample_id" "$cnr_dir" "$cns_dir" "$metrics_dir"
    cnvkit_run_call_ci "$sample_id" "${metrics_dir}/${sample_id}.segmetrics.cns" "$call_dir"
    report_cns="${call_dir}/${sample_id}.call.cns"
  fi
fi

if [[ "$need_reports" -eq 1 ]]; then
  if [[ "$enable_genemetrics" -eq 1 ]]; then
    cnvkit_run_genemetrics "$sample_id" "$cnr_dir" "$report_cns" "$metrics_dir"
  fi
  if [[ "$enable_export" -eq 1 ]]; then
    cnvkit_run_exports "$sample_id" "$report_cns" "$export_dir"
  fi
  if [[ "$enable_plots" -eq 1 ]]; then
    cnvkit_run_plots "$sample_id" "$cnr_dir" "$report_cns" "$plots_dir"
  fi
fi
