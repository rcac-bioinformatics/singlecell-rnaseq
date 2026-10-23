#!/bin/bash
# run_R_chain.sh: run the R episodes (3 to 8) of one kit run in THIS shell, without sbatch.
# For an srun shell (after module load of the R module) or an Open OnDemand RStudio terminal,
# where the R the learners get is the R on PATH: use --r current. submit_all.sh runs the same
# episodes as batch jobs instead.
# Usage:
#   bash run_R_chain.sh --run RUN_ID --pass continuous|cold [--episodes 03,04,...|env]
#                       [--r image|module|current]
# --pass continuous  the RDS chain as a learner runs it, in the run's learner copy; 03 reads
#                    the matrix its set-paths names (generated/ep3_matrix_path.txt). A failed
#                    episode skips the rest of its chain (03 to 06, or 07 to 08).
# --pass cold        each episode in its own fresh directory, seeded from $STAGED (rds/ and the
#                    filtered matrix); all run even if one fails.
# --episodes         default: all of 03 04 05 06 07 08 (cold: 03 04 05 06 08); "env" runs the
#                    environment and egress check.
# --r                default KIT_R_MODE from config.sh (image).
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""; PASS=""; EPS=""; MODE=$KIT_R_MODE
usage() { sed -n '2,19p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [ $# -gt 0 ]; do
  case $1 in
    --run) RUN_ID=$2; shift 2;;
    --pass) PASS=$2; shift 2;;
    --episodes) EPS=${2//,/ }; shift 2;;
    --r) MODE=$2; shift 2;;
    -h|--help) usage 0;;
    *) echo "unknown option $1" >&2; usage 1;;
  esac
done
[ -n "$RUN_ID" ] || kit_die "--run RUN_ID is required"
case $PASS in continuous|cold) ;; *) [ "$EPS" = env ] || kit_die "--pass continuous|cold is required";; esac
case $MODE in image|module|current) ;; *) kit_die "--r must be image, module or current";; esac
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN
kit_load_run
export TMPDIR="$RUN/tmp/chain.$$"; mkdir -p "$TMPDIR"
if [ -z "$EPS" ]; then
  if [ "$PASS" = cold ]; then EPS="03 04 05 06 08"; else EPS=$KIT_R_EPS; fi
fi
kit_log "R chain: run $RUN_ID, pass ${PASS:-n/a}, episodes $EPS, R mode $MODE"
nfail=0; broke_pbmc=0; broke_ifnb=0
for ep in $EPS; do
  case $ep in
    env)
      log="$RUN/logs/R-env-egress.$$.log"
      kit_marker env-egress RUNNING
      if kit_r "$MODE" env-egress "$KIT_USER_OVERRIDE" "$log" "$KIT_DIR/R/env_egress.R"; then kit_job_end env-egress PASS
      else kit_job_end env-egress WARN; fi
      cp "$log" "$RUN/logs/R-env-egress.log"; cp "$log" "$RUN/records/console/R-env-egress.log"
      tail -12 "$log" >&2;;
    0[3-8])
      if [ "$PASS" = continuous ] && [ "$ep" = 03 ] && grep -q '^cellranger_output/' "$KIT_GEN/ep3_matrix_path.txt" && ! kit_step_passed 02-cellranger; then
        kit_log "warning: 02-cellranger is not done in this run; episode 03 as written will not find its matrix"
      fi
      if [ "$PASS" = continuous ]; then name="r$ep"; else name="c$ep"; fi
      # continuous: 03 to 06 are one RDS chain, 07 to 08 another; a failure only stops its own chain
      if [ "$PASS" = continuous ] && { { [ "$ep" -le 6 ] && [ "$broke_pbmc" = 1 ]; } || { [ "$ep" -ge 7 ] && [ "$broke_ifnb" = 1 ]; }; }; then
        kit_log "skip $name: an earlier episode of its RDS chain failed"; continue
      fi
      kit_marker "$name" RUNNING
      t0=$(date +%s)
      if kit_r_episode "$MODE" "$PASS" "$ep"; then
        kit_log "$name done in $(( $(date +%s) - t0 )) s ($(cut -d' ' -f1 "$RUN/markers/$name.done"))"
      else
        kit_marker "$name" FAIL
        nfail=$((nfail + 1))
        if [ "$PASS" = continuous ]; then
          if [ "$ep" -le 6 ]; then broke_pbmc=1; else broke_ifnb=1; fi
          kit_log "$name failed; later episodes of the same RDS chain are skipped (03 to 06, or 07 to 08)"
        fi
      fi;;
    *) kit_die "unknown episode $ep";;
  esac
done
kit_log "R chain finished: $nfail failed step(s). Records: $RUN/records; collect with: bash $KIT_DIR/collect.sh --run $RUN_ID"
[ $nfail -eq 0 ]
