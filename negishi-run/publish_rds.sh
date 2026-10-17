#!/bin/bash
# publish_rds.sh: replace the five staged RDS files ($STAGED/rds/) with the ones a kit run produced.
# One of the two kit scripts allowed to write to $STAGED (the other is publish_staged.sh).
# Usage: bash publish_rds.sh --run RUN_ID [--dry-run]
# Sources are the handoff copies the continuous pass kept right after each save chunk
# ($RUN/handoff/r03/pbmc_filtered.rds ... r06/pbmc_annotated.rds, r07/ifnb_annotated.rds), so
# episode 8 overwriting ifnb_annotated.rds cannot leak into the staged copy. Steps r03 to r07 must
# have passed. Before touching anything it shows, per file, staged and new side by side (cells,
# features, clusters, idents, assays, layers, reductions, metadata columns, size, md5), using R in
# the configured environment, then asks you to type the destination directory; anything else
# cancels. --dry-run stops after the comparison.
# Replaced files are first copied to $KIT_SCRATCH_BASE/publish-backups/<time>/rds/. Nothing is deleted.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""; DRY=0
while [ $# -gt 0 ]; do
  case $1 in
    --run) RUN_ID=$2; shift 2;;
    --dry-run) DRY=1; shift;;
    -h|--help) sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) kit_die "unknown option $1";;
  esac
done
[ -n "$RUN_ID" ] || kit_die "--run RUN_ID is required"
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN; kit_load_run
DEST="$STAGED/rds"
declare -A FROM=( [pbmc_filtered.rds]=r03 [pbmc_normalized.rds]=r04 [pbmc_clustered.rds]=r05 [pbmc_annotated.rds]=r06 [ifnb_annotated.rds]=r07 )
NAMES="pbmc_filtered.rds pbmc_normalized.rds pbmc_clustered.rds pbmc_annotated.rds ifnb_annotated.rds"
new=(); old=()
for f in $NAMES; do
  s=${FROM[$f]}
  kit_is_done "$s" || kit_die "step $s has not passed in $RUN_ID (its marker: $(cat "$RUN/markers/$s.done" 2>/dev/null || echo none))"
  [ -s "$RUN/handoff/$s/$f" ] || kit_die "$RUN/handoff/$s/$f missing"
  new+=("$RUN/handoff/$s/$f"); old+=("$DEST/$f")
done

T="$RUN/tmp/publish_rds.$$"; mkdir -p "$T"
kit_log "reading the staged and new RDS files in R ($KIT_R_MODE mode; a few minutes)"
kit_r "$KIT_R_MODE" rds-summary-staged "$KIT_USER_OVERRIDE" "$T/staged.tsv" "$KIT_DIR/R/rds_summary.R" staged "${old[@]}" ||
  kit_die "R failed reading the staged RDS files; see $T/staged.tsv"
kit_r "$KIT_R_MODE" rds-summary-new "$KIT_USER_OVERRIDE" "$T/new.tsv" "$KIT_DIR/R/rds_summary.R" new "${new[@]}" ||
  kit_die "R failed reading the new RDS files; see $T/new.tsv"
python3 - "$T/staged.tsv" "$T/new.tsv" <<'PY' | tee "$T/side-by-side.txt"
import sys
def rows(p):
    lines = [l.rstrip("\n").split("\t") for l in open(p) if l.startswith(("label\t", "staged\t", "new\t"))]
    hdr = lines[0]
    return {r[1]: dict(zip(hdr, r)) for r in lines[1:]}
a, b = rows(sys.argv[1]), rows(sys.argv[2])
fields = ["cells", "features", "clusters", "idents", "assays", "layers", "reductions", "class", "bytes", "md5"]
for f in sorted(set(a) | set(b), key=lambda x: ["pbmc_filtered.rds", "pbmc_normalized.rds", "pbmc_clustered.rds", "pbmc_annotated.rds", "ifnb_annotated.rds"].index(x) if x.endswith(".rds") else 9):
    print("\n== %s" % f)
    print("  %-11s %-40s %-40s" % ("field", "staged", "new"))
    for k in fields:
        x, y = a.get(f, {}).get(k, "MISSING"), b.get(f, {}).get(k, "MISSING")
        print("  %-11s %-40s %-40s%s" % (k, x[:40], y[:40], "" if x == y else "   <- differs"))
    ma = set(a.get(f, {}).get("meta_columns", "").split(",")) - {""}
    mb = set(b.get(f, {}).get("meta_columns", "").split(",")) - {""}
    print("  meta columns in both: %s" % ", ".join(sorted(ma & mb)))
    if ma - mb: print("  only staged: %s" % ", ".join(sorted(ma - mb)))
    if mb - ma: print("  only new:    %s" % ", ".join(sorted(mb - ma)))
PY
TS=$(date +%Y%m%d-%H%M%S); BK="$KIT_SCRATCH_BASE/publish-backups/$TS/rds"
if [ "$STAGED_PERM" = group ]; then PF=g+r; PD=g+rx; else PF=o+r; PD=o+rx; fi
printf '\nCommands:\n'
printf '  mkdir -p "%s" "%s"\n' "$BK" "$DEST"
for i in "${!new[@]}"; do
  [ -f "${old[$i]}" ] && printf '  cp -p "%s" "%s/"   # backup\n' "${old[$i]}" "$BK"
  printf '  cp -p "%s" "%s" && chmod %s "%s"\n' "${new[$i]}" "${old[$i]}" "$PF" "${old[$i]}"
done
printf '  chmod %s "%s"\n' "$PD" "$DEST"
printf 'Side-by-side saved to %s\n' "$T/side-by-side.txt"
if [ $DRY = 1 ]; then printf '\ndry run: nothing written\n'; exit 0; fi
printf '\nTo publish, type the destination directory exactly (%s): ' "$DEST"
read -r ans
[ "$ans" = "$DEST" ] || kit_die "cancelled (typed '$ans')"
mkdir -p "$BK" "$DEST"
for i in "${!new[@]}"; do
  [ -f "${old[$i]}" ] && cp -p "${old[$i]}" "$BK/"
  cp -p "${new[$i]}" "${old[$i]}"
  chmod "$PF" "${old[$i]}"
done
chmod "$PD" "$DEST"
printf '%s\tpublish_rds\t%s\t%s\tbackup %s\n' "$(kit_ts)" "$RUN" "$DEST" "$BK" >> "$KIT_SCRATCH_BASE/publish.log"
kit_log "published 5 RDS files from run $RUN_ID; previous copies in $BK"
bash "$KIT_DIR/check_staged.sh" --quick || true
