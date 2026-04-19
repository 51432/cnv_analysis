#!/usr/bin/env bash
set -euo pipefail

cnvkit_prepare_bed3() {
  local input_bed="$1"
  local output_bed="$2"

  mkdir -p "$(dirname "$output_bed")"
  if [[ -s "$output_bed" ]]; then
    echo "[INFO] BED3 exists, skip: $output_bed"
    return 0
  fi

  awk 'BEGIN{FS=OFS="\t"} $0 !~ /^#/ && $0 !~ /^track/ && $0 !~ /^browser/ {if (NF<3) {exit 2}; print $1,$2,$3}' "$input_bed" > "$output_bed"
  echo "[INFO] Prepared BED3: $output_bed"
}

cnvkit_build_antitarget_bed() {
  local target_bed="$1"
  local access_bed="$2"
  local antitarget_bed="$3"

  if [[ -s "$antitarget_bed" ]]; then
    echo "[INFO] Antitarget BED exists, skip: $antitarget_bed"
    return 0
  fi

  mkdir -p "$(dirname "$antitarget_bed")"
  echo "[INFO] Building antitarget BED: $antitarget_bed"
  "${CNVKIT_CMD:-cnvkit.py}" antitarget "$target_bed" -g "$access_bed" -o "$antitarget_bed"
}

cnvkit_run_normal_coverage() {
  local normal_bam="$1"
  local target_bed="$2"
  local antitarget_bed="$3"
  local out_dir="$4"
  local threads="$5"

  local bam_base
  bam_base="$(basename "$normal_bam")"
  bam_base="${bam_base%.bam}"

  local target_cnn="${out_dir}/${bam_base}.targetcoverage.cnn"
  local antitarget_cnn="${out_dir}/${bam_base}.antitargetcoverage.cnn"

  mkdir -p "$out_dir"

  if [[ -s "$target_cnn" && -s "$antitarget_cnn" ]]; then
    echo "[INFO] Coverage exists, skip normal: $normal_bam"
    return 0
  fi

  if [[ "$threads" -ne 1 ]]; then
    echo "[WARN] For stability, cnvkit.py coverage is forced to -p 1 (array-level parallelism is still used)."
  fi

  local -a coverage_opts
  coverage_opts=(-p 1)
  if [[ "${CNVKIT_COVERAGE_USE_COUNT:-0}" -eq 1 ]]; then
    coverage_opts+=(--count)
    echo "[INFO] coverage mode: --count"
  fi

  echo "[INFO] Running target coverage for: $normal_bam"
  "${CNVKIT_CMD:-cnvkit.py}" coverage "$normal_bam" "$target_bed" "${coverage_opts[@]}" -o "$target_cnn"

  echo "[INFO] Running antitarget coverage for: $normal_bam"
  "${CNVKIT_CMD:-cnvkit.py}" coverage "$normal_bam" "$antitarget_bed" "${coverage_opts[@]}" -o "$antitarget_cnn"
}

cnvkit_run_tumor_coverage() {
  local sample_id="$1"
  local tumor_bam="$2"
  local target_bed="$3"
  local antitarget_bed="$4"
  local coverage_dir="$5"

  local target_cnn="${coverage_dir}/${sample_id}.targetcoverage.cnn"
  local antitarget_cnn="${coverage_dir}/${sample_id}.antitargetcoverage.cnn"
  mkdir -p "$coverage_dir"

  if [[ -s "$target_cnn" && -s "$antitarget_cnn" ]]; then
    echo "[INFO] Tumor coverage exists, skip: $sample_id"
    return 0
  fi

  local -a coverage_opts
  coverage_opts=(-p 1)
  if [[ "${CNVKIT_COVERAGE_USE_COUNT:-0}" -eq 1 ]]; then
    coverage_opts+=(--count)
  fi

  "${CNVKIT_CMD:-cnvkit.py}" coverage "$tumor_bam" "$target_bed" "${coverage_opts[@]}" -o "$target_cnn"
  "${CNVKIT_CMD:-cnvkit.py}" coverage "$tumor_bam" "$antitarget_bed" "${coverage_opts[@]}" -o "$antitarget_cnn"
}

cnvkit_run_fix() {
  local sample_id="$1"
  local coverage_dir="$2"
  local reference_cnn="$3"
  local cnr_dir="$4"

  local target_cnn="${coverage_dir}/${sample_id}.targetcoverage.cnn"
  local antitarget_cnn="${coverage_dir}/${sample_id}.antitargetcoverage.cnn"
  local cnr_out="${cnr_dir}/${sample_id}.cnr"

  [[ -s "$target_cnn" ]] || { echo "[ERROR] Missing target coverage: $target_cnn" >&2; return 1; }
  [[ -s "$antitarget_cnn" ]] || { echo "[ERROR] Missing antitarget coverage: $antitarget_cnn" >&2; return 1; }

  mkdir -p "$cnr_dir"
  if [[ -s "$cnr_out" ]]; then
    echo "[INFO] CNR exists, skip fix: $sample_id"
    return 0
  fi

  "${CNVKIT_CMD:-cnvkit.py}" fix "$target_cnn" "$antitarget_cnn" "$reference_cnn" -o "$cnr_out"
}

cnvkit_run_segment() {
  local sample_id="$1"
  local cnr_dir="$2"
  local cns_dir="$3"

  local cnr_in="${cnr_dir}/${sample_id}.cnr"
  local seg_out="${cns_dir}/${sample_id}.seg.cns"
  [[ -s "$cnr_in" ]] || { echo "[ERROR] Missing CNR: $cnr_in" >&2; return 1; }

  mkdir -p "$cns_dir"
  if [[ -s "$seg_out" ]]; then
    echo "[INFO] Segment CNS exists, skip: $sample_id"
    return 0
  fi

  "${CNVKIT_CMD:-cnvkit.py}" segment "$cnr_in" -o "$seg_out"
}

cnvkit_run_segmetrics_ci() {
  local sample_id="$1"
  local cnr_dir="$2"
  local cns_dir="$3"
  local metrics_dir="$4"

  local cnr_in="${cnr_dir}/${sample_id}.cnr"
  local seg_in="${cns_dir}/${sample_id}.seg.cns"
  local segmetrics_out="${metrics_dir}/${sample_id}.segmetrics.cns"
  mkdir -p "$metrics_dir"

  if [[ -s "$segmetrics_out" ]]; then
    echo "[INFO] Segmetrics exists, skip: $sample_id"
    return 0
  fi

  "${CNVKIT_CMD:-cnvkit.py}" segmetrics "$seg_in" -s "$cnr_in" --ci -o "$segmetrics_out"
}

cnvkit_run_call_ci() {
  local sample_id="$1"
  local input_cns="$2"
  local call_dir="$3"

  local call_out="${call_dir}/${sample_id}.call.cns"
  mkdir -p "$call_dir"
  if [[ -s "$call_out" ]]; then
    echo "[INFO] Call CNS exists, skip: $sample_id"
    return 0
  fi

  "${CNVKIT_CMD:-cnvkit.py}" call "$input_cns" --filter ci -o "$call_out"
}

cnvkit_run_genemetrics() {
  local sample_id="$1"
  local cnr_dir="$2"
  local cns_for_reports="$3"
  local metrics_dir="$4"

  local cnr_in="${cnr_dir}/${sample_id}.cnr"
  local out_tsv="${metrics_dir}/${sample_id}.genemetrics.tsv"
  mkdir -p "$metrics_dir"

  [[ -s "$cnr_in" ]] || { echo "[ERROR] Missing CNR: $cnr_in" >&2; return 1; }
  [[ -s "$cns_for_reports" ]] || { echo "[ERROR] Missing CNS for genemetrics: $cns_for_reports" >&2; return 1; }

  if [[ -s "$out_tsv" ]]; then
    echo "[INFO] Genemetrics exists, skip: $sample_id"
    return 0
  fi

  "${CNVKIT_CMD:-cnvkit.py}" genemetrics "$cnr_in" -s "$cns_for_reports" -o "$out_tsv"
}

cnvkit_run_exports() {
  local sample_id="$1"
  local cns_for_reports="$2"
  local export_dir="$3"

  mkdir -p "$export_dir"
  local seg_out="${export_dir}/${sample_id}.seg"
  local bed_out="${export_dir}/${sample_id}.bed"
  local vcf_out="${export_dir}/${sample_id}.vcf"

  [[ -s "$seg_out" ]] || "${CNVKIT_CMD:-cnvkit.py}" export seg "$cns_for_reports" -o "$seg_out"
  [[ -s "$bed_out" ]] || "${CNVKIT_CMD:-cnvkit.py}" export bed "$cns_for_reports" -o "$bed_out"
  [[ -s "$vcf_out" ]] || "${CNVKIT_CMD:-cnvkit.py}" export vcf "$cns_for_reports" -o "$vcf_out"
}

cnvkit_run_plots() {
  local sample_id="$1"
  local cnr_dir="$2"
  local cns_for_reports="$3"
  local plots_dir="$4"

  local cnr_in="${cnr_dir}/${sample_id}.cnr"
  local scatter_pdf="${plots_dir}/${sample_id}.scatter.pdf"
  local diagram_pdf="${plots_dir}/${sample_id}.diagram.pdf"

  mkdir -p "$plots_dir"
  [[ -s "$scatter_pdf" ]] || "${CNVKIT_CMD:-cnvkit.py}" scatter "$cnr_in" -s "$cns_for_reports" -o "$scatter_pdf"
  [[ -s "$diagram_pdf" ]] || "${CNVKIT_CMD:-cnvkit.py}" diagram "$cnr_in" -s "$cns_for_reports" -o "$diagram_pdf"
}

cnvkit_build_reference() {
  local coverage_dir="$1"
  local reference_fasta="$2"
  local output_reference="$3"

  local -a target_cov
  local -a antitarget_cov

  mapfile -t target_cov < <(find "$coverage_dir" -maxdepth 1 -name '*.targetcoverage.cnn' | sort)
  mapfile -t antitarget_cov < <(find "$coverage_dir" -maxdepth 1 -name '*.antitargetcoverage.cnn' | sort)

  if [[ ${#target_cov[@]} -eq 0 || ${#antitarget_cov[@]} -eq 0 ]]; then
    echo "[ERROR] No coverage CNN files found in ${coverage_dir}" >&2
    return 1
  fi

  if [[ ${#target_cov[@]} -ne ${#antitarget_cov[@]} ]]; then
    echo "[ERROR] targetcoverage and antitargetcoverage count mismatch in ${coverage_dir}" >&2
    return 1
  fi

  mkdir -p "$(dirname "$output_reference")"
  echo "[INFO] Building pooled normal CNVkit reference: $output_reference"
  "${CNVKIT_CMD:-cnvkit.py}" reference "${target_cov[@]}" "${antitarget_cov[@]}" -f "$reference_fasta" -o "$output_reference"
}
