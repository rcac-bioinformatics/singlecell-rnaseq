#!/bin/bash
# collect.sh: gather a kit run into one tarball to bring back (login node, 1 to 3 minutes).
# Usage: bash collect.sh --run RUN_ID
# Safe to rerun at any time (also while jobs are still running): sacct is queried again and a new
# timestamped tarball is written, so an earlier one is never overwritten:
#   $KIT_SCRATCH_BASE/scrna-kit-<RUN_ID>-<YYYYmmdd-HHMM>.tar.gz
# Contents: SUMMARY.md and tables/, sacct records, module versions, preflight, R console per chunk,
# events, anchors, sessionInfo, egress, episode 2 metrics, figures (raw per chunk, and copies named
# as in episodes/fig/), listings, markers and all logs. Handoff RDS files are not included (large;
# publish_rds.sh reads them on Negishi).
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""
while [ $# -gt 0 ]; do
  case $1 in --run) RUN_ID=$2; shift 2;; *) kit_die "usage: collect.sh --run RUN_ID";; esac
done
[ -n "$RUN_ID" ] || kit_die "usage: collect.sh --run RUN_ID"
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN
kit_load_run
NAME="scrna-kit-$RUN_ID-$(date +%Y%m%d-%H%M)"
C="$RUN/collect/$NAME"
mkdir -p "$C"/{figures,records,logs,checks}
kit_log "collecting $RUN into $C"

# ---- SLURM accounting for every job this run submitted (array tasks and job steps included)
ids=$(awk -F'\t' 'NR > 1 && $2 ~ /^[0-9]+$/ { printf "%s%s", sep, $2; sep = "," }' "$RUN/jobs.tsv" 2>/dev/null || true)
if [ -n "$ids" ] && command -v sacct >/dev/null 2>&1; then
  sacct -j "$ids" -P -n --format=JobID,JobName%40,Elapsed,MaxRSS,ReqMem,State,Submit,Start,End,AllocCPUS,Timelimit,ExitCode,TotalCPU \
    > "$RUN/records/sacct.txt" 2>/dev/null || kit_log "sacct failed; runtimes will be missing"
  python3 "$KIT_DIR/py/sacct_summary.py" "$RUN/records/sacct.txt" "$RUN/jobs.tsv" > "$RUN/records/sacct_summary.tsv" || true
fi

# ---- versions and preflight state
mkdir -p "$RUN/records/preflight"
latest=$(find "$KIT_SCRATCH_BASE/state" -maxdepth 1 -name 'preflight-*.txt' 2>/dev/null | sort | tail -1)
[ -n "$latest" ] && cp "$latest" "$RUN/records/preflight/"
[ -f "$KIT_SCRATCH_BASE/state/pkgcheck.txt" ] && cp "$KIT_SCRATCH_BASE/state/pkgcheck.txt" "$RUN/records/preflight/"
if command -v module >/dev/null 2>&1 || [ -n "${LMOD_PKG:-}" ]; then
  for mv in $EXPECTED_MODULES; do
    tool=${mv%%/*}
    printf '%s\t%s\n' "$tool" "$(bash -lc "module load biocontainers >/dev/null 2>&1; module load $tool >/dev/null 2>&1; module -t list 2>&1 | grep -E '^$tool/' | head -1" 2>/dev/null || echo '?')"
  done > "$RUN/records/module_versions.tsv"
fi
(cd "$KIT_DIR/.." && sha256sum -c "$KIT_GEN/episodes.sha256" > "$C/checks/episode_checksums.txt" 2>&1) || true

# ---- summary
python3 "$KIT_DIR/py/summarize.py" "$RUN" "$KIT_GEN" "$C"

# ---- figures: raw per chunk, and copies named as the committed episodes/fig/ files
cp -r "$RUN"/records/figures "$C/figures/raw" 2>/dev/null || true
while IFS=$'\t' read -r ep label k fig; do
  [ "$ep" = episode ] && continue
  for pass in r c; do
    src="$RUN/records/figures/$pass$ep/${ep}_${label}_$(printf '%02d' "$k").png"
    dst="$C/figures/episodes/$fig"; [ $pass = c ] && dst="$C/figures/episodes_cold/$fig"
    if [ -f "$src" ]; then mkdir -p "$(dirname "$dst")"; cp "$src" "$dst"; fi
  done
done < "$KIT_GEN/figmap.tsv"

# ---- records and logs (handoff RDS stay on Negishi)
for d in console events anchors sessioninfo egress env metrics listings preflight; do
  [ -d "$RUN/records/$d" ] && cp -r "$RUN/records/$d" "$C/records/"
done
cp "$RUN"/records/*.tsv "$RUN"/records/*.txt "$C/records/" 2>/dev/null || true
cp "$RUN/kit.env" "$RUN/config.sh.used" "$RUN"/jobs*.tsv "$C/records/" 2>/dev/null || true
cp -r "$RUN/markers" "$C/records/" 2>/dev/null || true
cp "$RUN"/logs/* "$C/logs/" 2>/dev/null || true
cp "$W"/*.out "$W"/*.err "$C/logs/" 2>/dev/null || true   # the learner scripts' own SLURM logs
cp "$KIT_GEN/steps.tsv" "$KIT_GEN/figmap.tsv" "$KIT_DIR/anchors.tsv" "$C/records/"
(cd "$RUN" && find handoff -type f -printf '%p\t%s\n' 2>/dev/null) > "$C/records/handoff-files.tsv" || true
date '+%F %T' > "$C/COLLECTED.txt"

TAR="$KIT_SCRATCH_BASE/$NAME.tar.gz"
tar -czf "$TAR" -C "$RUN/collect" "$NAME"
kit_log "wrote $TAR ($(du -h "$TAR" | cut -f1)); summary: $C/SUMMARY.md"
kit_log "bring it back (laptop): rsync -avP $(id -un)@negishi.rcac.purdue.edu:$TAR ~/svn/singlecell-rnaseq/negishi-run/records/ && tar -xzf ~/svn/singlecell-rnaseq/negishi-run/records/$NAME.tar.gz -C ~/svn/singlecell-rnaseq/negishi-run/records/"
