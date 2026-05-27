#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$ROOT"

OUTPUT_ROOT="${REPORT_REPRODUCTION_OUTPUT_ROOT:-${1:-}}"

if [[ -z "$OUTPUT_ROOT" ]]; then
  ./reproduce_project_data.sh all
  exit 0
fi

OUTPUT_ROOT="$(mkdir -p "$OUTPUT_ROOT" && cd "$OUTPUT_ROOT" && pwd)"
mkdir -p "$OUTPUT_ROOT/results" "$OUTPUT_ROOT/project_outputs"

if [[ -e "$ROOT/results" || -e "$ROOT/project_outputs" ]]; then
  echo "Refusing to create output links because results/ or project_outputs/ already exists in $ROOT" >&2
  echo "Move or remove those generated directories before running with an external output root." >&2
  exit 2
fi

cleanup() {
  [[ -L "$ROOT/results" ]] && rm "$ROOT/results"
  [[ -L "$ROOT/project_outputs" ]] && rm "$ROOT/project_outputs"
}
trap cleanup EXIT

ln -s "$OUTPUT_ROOT/results" "$ROOT/results"
ln -s "$OUTPUT_ROOT/project_outputs" "$ROOT/project_outputs"

./reproduce_project_data.sh all
