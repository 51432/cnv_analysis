#!/usr/bin/env bash
set -euo pipefail

: "${SCRIPT_DIR:?missing SCRIPT_DIR}"
: "${VALIDATED_SAMPLES:?missing VALIDATED_SAMPLES}"

line_no="$((SLURM_ARRAY_TASK_ID + 2))"
line="$(sed -n "${line_no}p" "$VALIDATED_SAMPLES")"
[[ -n "$line" ]] || { echo "Missing line for array index ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

IFS=$'\t' read -r sample_id tumor_bam normal_bam <<< "$line"

"${SCRIPT_DIR}/04_run_cnvkit_tumor_task.sh" \
  --sample-id "$sample_id" \
  --tumor-bam "$tumor_bam" \
  --stage "$STAGE" \
  --output-root "$OUTPUT_ROOT" \
  --reference-cnn "$REFERENCE_CNN" \
  --target-bed3 "$TARGET_BED3" \
  --antitarget-bed "$ANTITARGET_BED" \
  --antitarget-bed3 "$ANTITARGET_BED3" \
  --threads "$THREADS" \
  --enable-filtering "$ENABLE_FILTERING" \
  --enable-genemetrics "$ENABLE_GENEMETRICS" \
  --enable-export "$ENABLE_EXPORT" \
  --enable-plots "$ENABLE_PLOTS"
