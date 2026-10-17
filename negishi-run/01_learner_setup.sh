#!/bin/bash
# 01_learner_setup.sh: start a kit run by doing the learners/setup.md data copy exactly as
# written, into a fresh test directory, then the setup.md verification commands.
# Login node (this is where learners run it). Copies the whole staged tree (FASTQs about 30 GB
# plus the reference), so run it inside tmux or screen.
# Usage: bash 01_learner_setup.sh [RUN_ID]     (default RUN_ID: YYYY-MM-DD-HHMM)
# Creates $KIT_SCRATCH_BASE/runs/<RUN_ID>/ with scratch/scrna_workshop as the learner copy.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=${1:-$(date +%Y-%m-%d-%H%M)}
RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID"
if [ -e "$RUN/markers/01.done" ] && RUN=$RUN kit_is_done 01 2>/dev/null; then
  kit_log "run $RUN_ID already set up ($(cat "$RUN/markers/01.done")); nothing to do"
  exit 0
fi
[ -d "$STAGED" ] || kit_die "staged data $STAGED not found or not readable; run 00_preflight.sh"
grep -qF "rsync -avP $STAGED/ " "$KIT_GEN/setup/01-copy.sh" ||
  kit_die "learners/setup.md copies from '$(grep -o 'rsync -avP [^ ]*' "$KIT_GEN/setup/01-copy.sh")', but STAGED is $STAGED"
kit_init_run "$RUN"
export RUN
kit_load_run   # exports SCRATCH and RCAC_SCRATCH = $RUN/scratch (the test learner's scratch)
kit_marker 01 RUNNING
kit_log "run $RUN_ID: learner scratch is $RCAC_SCRATCH (staged source $STAGED is only read)"

# learner session: the setup.md blocks, verbatim
S="$RUN/tmp/01-setup-session.sh"
{
  echo "trap 'echo \"### KIT-ERR rc=\$? line=\$LINENO cmd=\$BASH_COMMAND\" >&2' ERR"
  # shellcheck disable=SC2016  # written into the session file, expanded there
  echo 'cd "$RCAC_SCRATCH"'
  for f in 01-copy 02-verify; do
    echo "echo \"### KIT-BLOCK setup/$f begin \$(date +%s)\""
    cat "$KIT_GEN/setup/$f.sh"
    echo "echo \"### KIT-BLOCK setup/$f end \$(date +%s)\""
  done
} > "$S"
kit_run_session 01 "$S" || { kit_marker 01 FAIL; kit_die "setup.md commands failed; see $RUN/logs/session-01.log"; }

# checks on what the learner now has
nfq=$(find "$W/fastq" -maxdepth 1 -name '*.fastq.gz' | wc -l)
nr=$(find "$W/fastq" -maxdepth 1 -name '*_R[12]_001.fastq.gz' | wc -l)
[ "$nr" = 4 ] || { kit_marker 01 FAIL; kit_die "learner copy has $nr R1/R2 FASTQ files, expected 4"; }
for f in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
  [ -s "$W/filtered_feature_bc_matrix/$f" ] || { kit_marker 01 FAIL; kit_die "filtered_feature_bc_matrix/$f missing in the learner copy"; }
done
for f in pbmc_filtered pbmc_normalized pbmc_clustered pbmc_annotated ifnb_annotated; do
  [ -s "$W/rds/$f.rds" ] || kit_log "WARNING: rds/$f.rds missing in the learner copy (the cold pass needs it)"
done
[ -d "$W/reference/refdata-gex-GRCh38-2024-A/star" ] || { kit_marker 01 FAIL; kit_die "reference/refdata-gex-GRCh38-2024-A/star missing"; }
# learners start without episode 2 outputs; a staged copy that has them would let finished results stand in
pre=$( (cd "$W" && find cellranger_output starsolo_output reference -maxdepth 1 \( -name pbmc10k -o -name 'starsolo_pbmc10k_*' -o -name star_index \) 2>/dev/null) | head -5) || true
if [ -n "$pre" ] && [ "${KIT_ALLOW_PRECOMPUTED:-0}" != 1 ]; then
  kit_marker 01 FAIL
  kit_die "the staged copy $STAGED contains episode 2 outputs, so this run would start with finished results:
$pre
Remove them from the staged directory, then start a new run ID. KIT_ALLOW_PRECOMPUTED=1 overrides."
fi
if [ -e "$W/reference/3M-february-2018.txt" ]; then
  kit_log "WARNING: the staged copy has reference/3M-february-2018.txt; episode 2 extracts it (it is not in the manifest)"
fi
LOG="$RUN/logs/session-01.log"
secs=$(awk '/KIT-BLOCK setup\/01-copy begin/{b=$NF} /KIT-BLOCK setup\/01-copy end/{e=$NF} END{print e-b}' "$LOG")
printf 'setup rsync seconds\t%s\nfastq files\t%s\nlearner copy size\t%s\n' "$secs" "$nfq" "$(du -sh "$W" | cut -f1)" > "$RUN/records/listings/01-setup.tsv"
(cd "$W" && find . -maxdepth 3 | sort) > "$RUN/records/listings/01-tree.txt"
kit_job_end 01 PASS
kit_log "setup done in ${secs}s ($nfq FASTQ files). Next: bash $KIT_DIR/submit_all.sh --run $RUN_ID --dry-run"
