#!/usr/bin/env bash
##############################################################################
# refilter_all.sh
#   Re-run bin/filter.py over every <DONOR>/ subfolder of a results
#   directory, producing confident_allelic_table / final_allelic_table /
#   metrics inside each donor folder.
#
# Usage:
#   refilter_all.sh <results_dir> <ref.fa> [filter.py args...]
#
# Example:
#   refilter_all.sh /lustre/.../results/mitchell_2025 \
#                   /lustre/.../resources/GRCh38.d1.vd1.fa
##############################################################################
set -euo pipefail

if (( $# < 2 )); then
  echo "Usage: $0 <results_dir> <ref.fa> [extra filter.py args...]" >&2
  exit 2
fi

RESULTS_DIR=$1
REF=$2
shift 2
EXTRA=( "$@" )

FILTER_PY="/lustre/scratch126/cellgen/behjati/ac87/MT_chemo/CoMiCall/bin/filter.py"

[[ -d "$RESULTS_DIR" ]] || { echo "ERROR: not a directory: $RESULTS_DIR" >&2; exit 2; }
[[ -r "$REF" ]]         || { echo "ERROR: cannot read ref: $REF" >&2; exit 2; }
[[ -r "$FILTER_PY" ]]   || { echo "ERROR: cannot find $FILTER_PY" >&2; exit 2; }

ok=0; fail=0
for donor_dir in "$RESULTS_DIR"/*/; do
  donor=$(basename "$donor_dir")
  allelic="$donor_dir/${donor}.allelic_table.tsv.gz"
  coverage="$donor_dir/${donor}.coverage_table.tsv.gz"

  if [[ ! -r "$allelic" || ! -r "$coverage" ]]; then
    echo "skip $donor — missing allelic/coverage table"
    continue
  fi

  echo "==> $donor"
  if ( cd "$donor_dir" && python "$FILTER_PY" \
         --allelic_table  "${donor}.allelic_table.tsv.gz" \
         --coverage_table "${donor}.coverage_table.tsv.gz" \
         --ref            "$REF" \
         --donor          "$donor" \
         "${EXTRA[@]}" ); then
    ok=$(( ok + 1 ))
  else
    echo "FAILED $donor"
    fail=$(( fail + 1 ))
  fi
done

echo "----"
echo "ok:     $ok"
echo "failed: $fail"
