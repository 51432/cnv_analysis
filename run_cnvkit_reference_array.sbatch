#!/usr/bin/env bash
set -euo pipefail

: "${SCRIPT_DIR:?missing SCRIPT_DIR}"
: "${NORMAL_LIST:?missing NORMAL_LIST}"
: "${COVERAGE_DIR:?missing COVERAGE_DIR}"
: "${ANTITARGET_BED:?missing ANTITARGET_BED}"
: "${TARGET_BED3:?missing TARGET_BED3}"
: "${ANTITARGET_BED3:?missing ANTITARGET_BED3}"
: "${THREADS:?missing THREADS}"

normal_bam="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$NORMAL_LIST")"
if [[ -z "$normal_bam" ]]; then
  echo "No normal_bam for array index ${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

"${SCRIPT_DIR}/02_run_cnvkit_reference_task.sh" \
  --normal-bam "$normal_bam" \
  --coverage-dir "$COVERAGE_DIR" \
  --antitarget-bed "$ANTITARGET_BED" \
  --target-bed3 "$TARGET_BED3" \
  --antitarget-bed3 "$ANTITARGET_BED3" \
  --threads "$THREADS"
