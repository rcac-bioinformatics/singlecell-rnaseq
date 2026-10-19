#!/bin/bash
# publish_staged.sh: copy checked files into the staged learner copy on Depot ($STAGED).
# One of the two kit scripts allowed to write there (the other is publish_rds.sh).
# Usage:
#   bash publish_staged.sh --from STAGING_DIR --part fastq|reference [--dry-run]
#   bash publish_staged.sh --part readme [--dry-run]
#   bash publish_staged.sh --matrix-from-run RUN_ID [--dry-run]
#   bash publish_staged.sh --matrix-from-path DIR [--dry-run]
# --from             a staging tree written by stage_fetch.sh
# --part readme      staged-top/README.txt from this kit into the top of $STAGED (refused while it
#                    still contains FILL-IN)
# --matrix-from-run  the filtered_feature_bc_matrix/ that run's 02-cellranger produced (the step and
#                    02-compare must have passed); the 10x public matrix is not a substitute
# --matrix-from-path an existing filtered_feature_bc_matrix/ outside $STAGED (an earlier run of the
#                    episode 2 script); its gene and barcode counts must equal the episode 3
#                    output block (anchors.tsv ep3_read_genes, ep3_read_cells). 02-compare later
#                    checks a fresh Cell Ranger run against it.
# Always prints, before touching anything: every file as add / replace / same (md5 compared when
# sizes agree), and the exact commands it will run. Then it asks you to type the destination
# directory; anything else cancels. --dry-run stops after the listing.
# Files it replaces are first copied to $KIT_SCRATCH_BASE/publish-backups/<time>/. It never deletes
# files from $STAGED. Afterwards it makes the new files readable (STAGED_PERM) and runs
# check_staged.sh --quick.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

FROM=""; PART=""; RUN_ID=""; MPATH=""; DRY=0
while [ $# -gt 0 ]; do
  case $1 in
    --from) FROM=$2; shift 2;;
    --part) PART=$2; shift 2;;
    --matrix-from-run) RUN_ID=$2; shift 2;;
    --matrix-from-path) MPATH=$2; shift 2;;
    --dry-run) DRY=1; shift;;
    -h|--help) sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) kit_die "unknown option $1";;
  esac
done
if [ -n "$RUN_ID" ]; then
  case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
  export RUN; kit_load_run
  kit_step_passed 02-cellranger || kit_die "02-cellranger has not passed in $RUN_ID"
  kit_is_done 02-compare || kit_die "02-compare has not passed in $RUN_ID; check its SUMMARY.md (episode 2 anchors) first"
  SRC="$W/cellranger_output/pbmc10k/outs/filtered_feature_bc_matrix"; DEST="$STAGED/filtered_feature_bc_matrix"
elif [ -n "$MPATH" ]; then
  SRC=${MPATH%/}; DEST="$STAGED/filtered_feature_bc_matrix"
  for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do [ -s "$SRC/$f" ] || kit_die "$SRC/$f missing"; done
  want_g=$(awk -F'\t' '$3 == "ep3_read_genes" {print $4}' "$KIT_DIR/anchors.tsv")
  want_c=$(awk -F'\t' '$3 == "ep3_read_cells" {print $4}' "$KIT_DIR/anchors.tsv")
  got_g=$(zcat "$SRC/features.tsv.gz" | wc -l); got_c=$(zcat "$SRC/barcodes.tsv.gz" | wc -l)
  printf 'matrix %s: %s genes (episode 3 says %s), %s barcodes (episode 3 says %s)\n' "$SRC" "$got_g" "$want_g" "$got_c" "$want_c"
  { [ "$got_g" = "$want_g" ] && [ "$got_c" = "$want_c" ]; } || kit_die "counts differ from the episode 3 output block; not publishing"
elif [ "$PART" = readme ]; then
  SRC="$KIT_DIR/staged-top"; DEST="$STAGED"
  ! grep -n 'FILL-IN' "$SRC/README.txt" || kit_die "$SRC/README.txt still has FILL-IN lines (above); fill them first"
else
  { [ -n "$FROM" ] && [ -n "$PART" ]; } || kit_die "give --from DIR --part fastq|reference, --matrix-from-run RUN_ID or --matrix-from-path DIR"
  case $PART in
    fastq) SRC="$FROM/fastq"; DEST="$STAGED/fastq";;
    reference) SRC="$FROM/reference/refdata-gex-GRCh38-2024-A"; DEST="$STAGED/reference/refdata-gex-GRCh38-2024-A";;
    *) kit_die "--part must be fastq or reference";;
  esac
fi
[ -d "$SRC" ] || kit_die "source $SRC not found"
case "$(readlink -m "$SRC")/" in "$(readlink -m "$STAGED")"/*) kit_die "source $SRC is inside \$STAGED";; esac
case "$(readlink -m "$DEST")/" in "$(readlink -m "$STAGED")"/*) ;; *) kit_die "destination $DEST is not inside \$STAGED";; esac
[ -d "$STAGED" ] || kit_die "\$STAGED=$STAGED does not exist; create it by hand first (this script does not create the top level)"

PLAN=$(mktemp)
python3 "$KIT_DIR/py/staged.py" plan "$SRC" "$DEST" > "$PLAN"
printf '\nPublish %s\n   into %s\n\n' "$SRC" "$DEST"
column -t -s$'\t' "$PLAN" | sed 's/^/  /'
nadd=$(awk -F'\t' '$1 == "add"' "$PLAN" | wc -l); nrep=$(awk -F'\t' '$1 == "replace"' "$PLAN" | wc -l)
printf '\n%d to add, %d to replace, %d unchanged. Nothing is deleted.\n' "$nadd" "$nrep" "$(awk -F'\t' '$1 == "same"' "$PLAN" | wc -l)"
TS=$(date +%Y%m%d-%H%M%S); BK="$KIT_SCRATCH_BASE/publish-backups/$TS"
if [ "$STAGED_PERM" = group ]; then PF=g+r; PD=g+rx; else PF=o+r; PD=o+rx; fi
printf '\nCommands:\n'
while IFS=$'\t' read -r act rel _; do
  case $act in action|same) continue;; esac
  [ "$act" = replace ] && printf '  mkdir -p "%s" && cp -p "%s" "%s"   # backup\n' "$(dirname "$BK/$rel")" "$DEST/$rel" "$BK/$rel"
  printf '  mkdir -p "%s" && cp -p "%s" "%s" && chmod %s "%s"\n' "$(dirname "$DEST/$rel")" "$SRC/$rel" "$DEST/$rel" "$PF" "$DEST/$rel"
done < "$PLAN"
DMAX=(); [ "$PART" = readme ] && DMAX=(-maxdepth 0)   # the README touches only the top directory
printf '  find "%s" %s-type d -exec chmod %s {} +\n' "$DEST" "${DMAX[*]:+${DMAX[*]} }" "$PD"
if [ $DRY = 1 ]; then printf '\ndry run: nothing written\n'; rm -f "$PLAN"; exit 0; fi
if [ "$nadd" = 0 ] && [ "$nrep" = 0 ]; then printf '\nnothing to publish\n'; rm -f "$PLAN"; exit 0; fi
printf '\nTo publish, type the destination directory exactly (%s): ' "$DEST"
read -r ans
[ "$ans" = "$DEST" ] || { rm -f "$PLAN"; kit_die "cancelled (typed '$ans')"; }

while IFS=$'\t' read -r act rel _; do
  case $act in action|same) continue;; esac
  if [ "$act" = replace ]; then mkdir -p "$(dirname "$BK/$rel")"; cp -p "$DEST/$rel" "$BK/$rel"; fi
  mkdir -p "$(dirname "$DEST/$rel")"
  cp -p "$SRC/$rel" "$DEST/$rel"
  chmod "$PF" "$DEST/$rel"
done < "$PLAN"
find "$DEST" "${DMAX[@]}" -type d -exec chmod "$PD" {} +
printf '%s\tpublish_staged\t%s\t%s\t%d added\t%d replaced\tbackup %s\n' "$(kit_ts)" "$SRC" "$DEST" "$nadd" "$nrep" "$BK" >> "$KIT_SCRATCH_BASE/publish.log"
rm -f "$PLAN"
kit_log "published; log $KIT_SCRATCH_BASE/publish.log; replaced files backed up in $BK"
bash "$KIT_DIR/check_staged.sh" --quick || true
