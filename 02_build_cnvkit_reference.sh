#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${script_dir}/config/00_config.sh"
source "${script_dir}/lib/steps_cnvkit.sh"

: "${COVERAGE_DIR:?missing COVERAGE_DIR}"
: "${REFERENCE_OUT:?missing REFERENCE_OUT}"
: "${OVERWRITE_REFERENCE:=0}"

if [[ -s "$REFERENCE_OUT" && "$OVERWRITE_REFERENCE" -ne 1 ]]; then
  echo "[INFO] Reference exists, skip: $REFERENCE_OUT"
  exit 0
fi

cnvkit_build_reference "$COVERAGE_DIR" "$REFERENCE" "$REFERENCE_OUT"
echo "[INFO] Reference generated: $REFERENCE_OUT"
