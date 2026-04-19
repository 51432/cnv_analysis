#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <samples.tsv> <validated_output.tsv>" >&2
  exit 1
fi

samples_tsv="$1"
validated_out="$2"

if [[ ! -f "$samples_tsv" ]]; then
  echo "[ERROR] line=1 sample_id=NA message=samples.tsv not found: ${samples_tsv}" >&2
  exit 1
fi
if [[ ! -r "$samples_tsv" ]]; then
  echo "[ERROR] line=1 sample_id=NA message=samples.tsv not readable: ${samples_tsv}" >&2
  exit 1
fi

expected_header=$'sample_id\ttumor_bam\tnormal_bam'
header_line="$(head -n 1 "$samples_tsv" || true)"
if [[ "$header_line" != "$expected_header" ]]; then
  echo "[ERROR] line=1 sample_id=NA message=header must be exactly: sample_id\ttumor_bam\tnormal_bam" >&2
  exit 1
fi

tmp_out="${validated_out}.tmp"
mkdir -p "$(dirname "$validated_out")"
echo -e "$expected_header" > "$tmp_out"

declare -A seen_sample_ids=()
error_count=0
line_no=1

while IFS= read -r raw_line || [[ -n "$raw_line" ]]; do
  ((line_no += 1))

  tab_count="$(awk -F'\t' '{print NF-1}' <<< "$raw_line")"
  if [[ "$tab_count" -ne 2 ]]; then
    sid="$(cut -f1 <<< "$raw_line")"
    sid="${sid:-NA}"
    echo "[ERROR] line=${line_no} sample_id=${sid} message=row must be TAB-delimited with exactly 3 columns" >&2
    ((error_count += 1))
    continue
  fi

  IFS=$'\t' read -r sample_id tumor_bam normal_bam <<< "$raw_line"

  if [[ -z "${sample_id:-}" || -z "${tumor_bam:-}" || -z "${normal_bam:-}" ]]; then
    sid="${sample_id:-NA}"
    echo "[ERROR] line=${line_no} sample_id=${sid} message=each row must contain non-empty sample_id/tumor_bam/normal_bam" >&2
    ((error_count += 1))
    continue
  fi

  if [[ -n "${seen_sample_ids[$sample_id]:-}" ]]; then
    echo "[ERROR] line=${line_no} sample_id=${sample_id} message=duplicate sample_id" >&2
    ((error_count += 1))
    continue
  fi
  seen_sample_ids["$sample_id"]=1

  if [[ ! -f "$tumor_bam" || ! -r "$tumor_bam" ]]; then
    echo "[ERROR] line=${line_no} sample_id=${sample_id} message=tumor_bam not found/readable: ${tumor_bam}" >&2
    ((error_count += 1))
  fi

  if [[ ! -f "$normal_bam" || ! -r "$normal_bam" ]]; then
    echo "[ERROR] line=${line_no} sample_id=${sample_id} message=normal_bam not found/readable: ${normal_bam}" >&2
    ((error_count += 1))
  fi

  echo -e "${sample_id}\t${tumor_bam}\t${normal_bam}" >> "$tmp_out"
done < <(tail -n +2 "$samples_tsv")

if [[ "$line_no" -eq 1 ]]; then
  echo "[ERROR] line=2 sample_id=NA message=no sample rows found" >&2
  ((error_count += 1))
fi

if (( error_count > 0 )); then
  rm -f "$tmp_out"
  exit 1
fi

mv "$tmp_out" "$validated_out"
echo "Validated samples TSV: ${validated_out}"
