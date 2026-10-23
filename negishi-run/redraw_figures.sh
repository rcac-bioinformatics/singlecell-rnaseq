#!/bin/bash
# redraw_figures.sh: redraw the committed figures listed in figures_redraw.tsv from a kit run's
# handoff RDS files, at their committed pixel size, with the episode's own chunk code.
# Usage: bash redraw_figures.sh --run RUN_ID
# The kit's own capture (collect.sh figures/episodes/) is 1400 x 1000 px, smaller than the
# committed files; use this when a code change means a committed figure must be replaced.
# Runs R in the configured environment (KIT_R_MODE; image mode uses OOD_SIF and the host
# library) on the login node, about a minute. Writes $RUN/records/figures/redraw/, which
# collect.sh brings back under figures/raw/redraw/; copy those files into episodes/fig/.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""
while [ $# -gt 0 ]; do
  case $1 in
    --run) RUN_ID=$2; shift 2;;
    -h|--help) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) kit_die "unknown option $1";;
  esac
done
[ -n "$RUN_ID" ] || kit_die "--run RUN_ID is required"
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN; kit_load_run
python3 "$KIT_DIR/build_kit.py" --check >&2 || kit_die "generated/ is stale or was edited; rebuild with python3 $KIT_DIR/build_kit.py and re-send"

TAB="$KIT_DIR/figures_redraw.tsv"
while read -r s; do
  kit_is_done "$s" || kit_die "step $s has not passed in $RUN_ID; its handoff RDS is not trustworthy"
done < <(awk -F'\t' '!/^#/ && $1 != "fig" { split($4, a, "/"); print a[1] }' "$TAB" | sort -u)
OUT="$RUN/records/figures/redraw"; LOG="$RUN/logs/redraw_figures.$$.log"
kit_log "redrawing $(grep -vc '^#\|^fig' "$TAB") figure(s) into $OUT ($KIT_R_MODE mode)"
kit_r "$KIT_R_MODE" redraw-figures "$KIT_USER_OVERRIDE" "$LOG" "$KIT_DIR/R/redraw_figures.R" \
  "$TAB" "$RUN/handoff" "$KIT_GEN" "$OUT" || kit_die "R failed; see $LOG"
grep -P '\t\d+x\d+\t' "$LOG" | while IFS= read -r l; do kit_log "  $l"; done
kit_log "done; run collect.sh --run $RUN_ID and copy figures/raw/redraw/*.png into episodes/fig/"
