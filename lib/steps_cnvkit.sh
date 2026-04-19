#!/usr/bin/env bash
set -euo pipefail

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
  cnvkit.py antitarget "$target_bed" -g "$access_bed" -o "$antitarget_bed"
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

  echo "[INFO] Running target coverage for: $normal_bam"
  cnvkit.py coverage "$normal_bam" "$target_bed" -p "$threads" -o "$target_cnn"

  echo "[INFO] Running antitarget coverage for: $normal_bam"
  cnvkit.py coverage "$normal_bam" "$antitarget_bed" -p "$threads" -o "$antitarget_cnn"
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
  cnvkit.py reference "${target_cov[@]}" "${antitarget_cov[@]}" -f "$reference_fasta" -o "$output_reference"
}
