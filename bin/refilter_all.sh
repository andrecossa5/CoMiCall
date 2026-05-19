#!/usr/bin/env bash
##############################################################################
# refilter_all.sh
#   Re-run bin/filter.py (then bin/genotype.py) on every <DONOR>/ subfolder
#   of a results directory. For each donor, writes annotated_allelic_table /
#   confident_allelic_table / metrics / genotyped_allelic_table inside that
#   donor's folder (alongside the existing tables).
#
# Usage:
#   refilter_all.sh <results_dir> <ref.fa> [extra filter.py args...]
#
# Note: extra args are forwarded only to filter.py. genotype.py uses defaults.
##############################################################################
set -euo pipefail

# --- Argument parsing -------------------------------------------------------
if (( $# < 2 )); then
  echo "Usage: $0 <results_dir> <ref.fa> [extra filter.py args...]" >&2
  exit 2
fi

RESULTS_DIR=$1   # top-level dir containing per-donor subdirs
REF=$2           # faidx-indexed FASTA used by filter.py (3nt context)
shift 2
EXTRA=( "$@" )   # forwarded verbatim to filter.py (threshold overrides etc.)

# Hard-coded paths to the scripts — keep aligned with pipeline layout
FILTER_PY="/lustre/scratch126/cellgen/behjati/ac87/MT_chemo/CoMiCall/bin/filter.py"
GENOTYPE_PY="/lustre/scratch126/cellgen/behjati/ac87/MT_chemo/CoMiCall/bin/genotype.py"

# --- Sanity checks ----------------------------------------------------------
[[ -d "$RESULTS_DIR" ]] || { echo "ERROR: not a directory: $RESULTS_DIR" >&2; exit 2; }
[[ -r "$REF" ]]         || { echo "ERROR: cannot read ref: $REF" >&2; exit 2; }
[[ -r "$FILTER_PY" ]]   || { echo "ERROR: cannot find $FILTER_PY" >&2; exit 2; }
[[ -r "$GENOTYPE_PY" ]] || { echo "ERROR: cannot find $GENOTYPE_PY" >&2; exit 2; }

# --- Build donor list (skip non-dirs and donors missing inputs) -------------
donor_dirs=()
for d in "$RESULTS_DIR"/*/; do
  donor=$(basename "$d")
  if [[ -r "$d/${donor}.allelic_table.tsv.gz" \
     && -r "$d/${donor}.coverage_table.tsv.gz" ]]; then
    donor_dirs+=( "$d" )
  else
    echo "skip $donor — missing allelic/coverage table"
  fi
done

total=${#donor_dirs[@]}
if (( total == 0 )); then
  echo "no donor folders with required inputs under $RESULTS_DIR" >&2
  exit 1
fi
echo "found $total donor folder(s) to process"

# --- Main loop --------------------------------------------------------------
ok=0; fail=0; i=0
for donor_dir in "${donor_dirs[@]}"; do
  i=$(( i + 1 ))
  donor=$(basename "$donor_dir")
  printf '[%d/%d] %s\n' "$i" "$total" "$donor"

  # cd into the donor dir so filter.py's relative output paths land beside
  # the inputs; run in a subshell so the cd does not leak between donors.
  if ( cd "$donor_dir" && python "$FILTER_PY" \
         --allelic_table  "${donor}.allelic_table.tsv.gz" \
         --coverage_table "${donor}.coverage_table.tsv.gz" \
         --ref            "$REF" \
         --donor          "$donor" \
         "${EXTRA[@]}" \
       && python "$GENOTYPE_PY" \
         --annotated_allelic_table "${donor}.annotated_allelic_table.tsv.gz" \
         --confident_allelic_table "${donor}.confident_allelic_table.tsv.gz" \
         --donor                   "$donor" ); then
    ok=$(( ok + 1 ))
  else
    echo "  FAILED $donor"
    fail=$(( fail + 1 ))
  fi
done

# --- Summary ----------------------------------------------------------------
echo "----"
printf 'ok:     %d / %d\n' "$ok" "$total"
printf 'failed: %d / %d\n' "$fail" "$total"
(( fail == 0 ))
