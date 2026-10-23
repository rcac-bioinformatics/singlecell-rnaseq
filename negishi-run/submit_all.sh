#!/bin/bash
# submit_all.sh: submit the whole step graph for one kit run (login node, a few minutes because
# step 02-prep runs the episode 2 login-node commands first, as learners do).
# Usage:
#   bash submit_all.sh --run RUN_ID [--only all|shell|r|cold] [--from STEP] [--force] [--dry-run]
# Steps and dependencies come from generated/steps.tsv (built from the episodes). A step that is
# already done in this run (marker, or outputs for the learner's own scripts) is skipped unless
# --force; --from STEP skips everything before STEP in the table. Job IDs go to $RUN/jobs.tsv.
# --only shell: episode 2; r: env-egress and the continuous R chain; cold: the cold pass.
# --dry-run submits nothing: it prints each sbatch command, checks every script with bash -n, and
# asks SLURM to validate the job headers with sbatch --test-only.
# Kit jobs use KIT_ACCOUNT; R jobs use OOD_ACCOUNT, OOD_CORES and OOD_TIME (config.sh), as the
# OOD form does. The learner scripts are submitted as written (their own headers).
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

RUN_ID=""; ONLY=all; FROM=""; FORCE=0; DRY=0
usage() { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [ $# -gt 0 ]; do
  case $1 in
    --run) RUN_ID=$2; shift 2;;
    --only) ONLY=$2; shift 2;;
    --from) FROM=$2; shift 2;;
    --force) FORCE=1; shift;;
    --dry-run) DRY=1; shift;;
    -h|--help) usage 0;;
    *) echo "unknown option $1" >&2; usage 1;;
  esac
done
# A bind variable in this shell changes what every container sees; learners do not have one.
for v in SINGULARITY_BIND APPTAINER_BIND; do
  [ -z "${!v:-}" ] || kit_die "FAIL: $v is set in this shell ($v=${!v}); unset it (learners do not have it) and resubmit"
done
case $ONLY in all|shell|r|cold) ;; *) kit_die "--only must be all, shell, r or cold";; esac
[ -n "$RUN_ID" ] || kit_die "--run RUN_ID is required (the ID given to 01_learner_setup.sh)"
python3 "$KIT_DIR/build_kit.py" --check >&2 || kit_die "generated/ is stale or was edited; rebuild with python3 $KIT_DIR/build_kit.py and re-send"
if [ -n "$FROM" ] && ! awk -F'\t' -v s="$FROM" '$1 == s { f = 1 } END { exit !f }' "$KIT_GEN/steps.tsv"; then
  kit_die "--from $FROM: no such step; steps are: $(awk -F'\t' 'NR > 1 { printf "%s ", $1 }' "$KIT_GEN/steps.tsv")"
fi
case "$RUN_ID" in /*) RUN=$RUN_ID;; *) RUN="$KIT_SCRATCH_BASE/runs/$RUN_ID";; esac
export RUN
if [ $DRY = 1 ] && [ ! -f "$RUN/kit.env" ]; then
  kit_log "dry run: $RUN has no kit.env; using a temporary run directory (nothing is written to $RUN)"
  RUN=$(mktemp -d "${TMPDIR:-/tmp}/kit-dryrun.XXXXXX")
  KIT_ALLOW_ANY_ROOT=1 kit_init_run "$RUN"
  # shellcheck source=/dev/null
  source "$RUN/kit.env"
  export SCRATCH="$KIT_LEARNER_SCRATCH" RCAC_SCRATCH="$KIT_LEARNER_SCRATCH" W
else
  kit_load_run
  kit_is_done 01 || [ $DRY = 1 ] || kit_die "setup not done for $RUN_ID; run 01_learner_setup.sh $RUN_ID first"
fi
case "$OOD_ACCOUNT$R_MODULE$OOD_SIF" in *FILL-IN*)
  if [ "$ONLY" != shell ]; then
    [ "$KIT_R_MODE" = image ] && case $OOD_SIF in *FILL-IN*) kit_die "OOD_SIF is still a placeholder in config.sh";; esac
    [ "$KIT_R_MODE" = module ] && case $R_MODULE in *FILL-IN*) kit_die "R_MODULE is still a placeholder in config.sh";; esac
    case $OOD_ACCOUNT in *FILL-IN*) kit_die "OOD_ACCOUNT is still a placeholder in config.sh";; esac
  fi;;
esac

step_done() { [ $FORCE = 1 ] && return 1; kit_step_passed "$1"; }
in_scope() {
  case $ONLY in
    all) return 0;;
    shell) case $1 in 02-*) return 0;; esac;;
    r) case $1 in env-egress|r0*) return 0;; esac;;
    cold) case $1 in c0*) return 0;; esac;;
  esac
  return 1
}

# refuse while an earlier submission for this run is still queued or running: two copies of a
# step would work in the same learner directory at the same time
if [ $DRY = 0 ] && [ -f "$RUN/jobs.tsv" ]; then
  prev=$(awk -F'\t' 'NR > 1 && $2 ~ /^[0-9]+$/ { printf "%s%s", s, $2; s = "," }' "$RUN/jobs.tsv")
  if [ -n "$prev" ]; then
    live=$(squeue -h -j "$prev" -o '%i %j %T' 2>/dev/null || true)
    [ -z "$live" ] || kit_die "jobs from an earlier submission of this run are still queued or running; wait for them or scancel them first:
$live"
  fi
fi

declare -A JOB
SUBMITTED=0; FAKE=9000000; DRYFAIL=0
JOBS_TSV="$RUN/jobs.tsv"
[ $DRY = 1 ] && JOBS_TSV="$RUN/jobs-dryrun.tsv"
[ -f "$JOBS_TSV" ] || printf 'step\tjobid\tsubmitted\tdeps\tflags\n' > "$JOBS_TSV"

sbatch_or_dry() {  # workdir args... ; prints the job id; returns 1 if submission (or test-only) failed
  local wd=$1; shift
  if [ $DRY = 1 ]; then
    local args=() a r
    for a in "$@"; do case $a in --dependency=*|--parsable|--kill-on-invalid-dep=*) ;; *) args+=("$a");; esac; done
    printf '  (cd %s && sbatch %s)\n' "$wd" "$*" >&2
    r=$(cd "$wd" 2>/dev/null || cd "$RUN"; sbatch --test-only "${args[@]}" 2>&1 | tail -1) || true
    echo "$((FAKE + SUBMITTED + 1))"
    case "$r" in *"to start at"*) printf '  test-only: %s\n' "$r" >&2;;
      *) printf '  test-only FAILED: %s\n' "$r" >&2; return 1;; esac
  else
    (cd "$wd" && sbatch "$@")
  fi
}

started=0; [ -z "$FROM" ] && started=1
while IFS=$'\t' read -r step kind deps dep_type _ _ script desc; do
  [ "$step" = step ] && continue
  if [ $started = 0 ]; then
    if [ "$step" = "$FROM" ]; then started=1; else continue; fi
  fi
  in_scope "$step" || continue
  if [ "$step" != 02-prep ] && step_done "$step"; then kit_log "skip $step (already done)"; continue; fi

  # dependencies: only on steps submitted in this invocation; others must be done already
  dl=""
  for d in ${deps//,/ }; do
    [ "$d" = - ] && continue
    if [ -n "${JOB[$d]:-}" ]; then dl="$dl:${JOB[$d]}"
    elif kit_step_passed "$d" || [ "$dep_type" = afterany ]; then :
    elif [ $DRY = 1 ] || ! in_scope "$d"; then kit_log "note: $step depends on $d, which is not done (submit it, or run --only all)"
    else kit_die "$step depends on $d, which is neither done nor submitted; use --from $d"; fi
  done
  common=(--parsable --kill-on-invalid-dep=yes)
  [ -n "$KIT_QOS" ] && common+=("--qos=$KIT_QOS")
  [ -n "$dl" ] && common+=("--dependency=$dep_type$dl")

  case $kind in
    login)
      bash -n "$KIT_GEN/sessions/$step.sh"
      if [ $DRY = 1 ]; then kit_log "dry run: would run $step on this login node (generated/sessions/$step.sh)"; continue; fi
      if kit_is_done "$step" && [ $FORCE = 0 ]; then kit_log "skip $step (already done)"; continue; fi
      kit_log "run $step on the login node ($desc)"
      KIT_STEP=$step; kit_marker "$step" RUNNING
      if kit_run_session "$step" "$KIT_GEN/sessions/$step.sh"; then kit_job_end "$step" PASS
      else kit_marker "$step" FAIL; kit_die "$step failed; see $RUN/logs/session-$step.log"; fi
      continue;;
    learner)
      bash -n "$KIT_GEN/learner-scripts/$script"
      if [ $DRY = 0 ]; then
        cmp -s "$W/$script" "$KIT_GEN/learner-scripts/$script" ||
          kit_die "$W/$script is missing or differs from the episode text; run 02-prep (submit_all.sh --run $RUN_ID --from 02-prep)"
        wd=$W
      else
        wd="$KIT_GEN/learner-scripts"
      fi
      kit_log "submit $step ($desc)${dl:+ after${dl//:/ }}"
      rc=0; jid=$(sbatch_or_dry "$wd" "${common[@]}" "$script") || rc=$?;;
    kitjob)
      bash -n "$KIT_GEN/jobs/$step.sh"; bash -n "$KIT_GEN/sessions/$step.sh"
      kit_log "submit $step ($desc)${dl:+ after${dl//:/ }}"
      rc=0; jid=$(sbatch_or_dry "$RUN" "${common[@]}" "--account=$KIT_ACCOUNT" -o "$RUN/logs/$step.%j.out" -e "$RUN/logs/$step.%j.out" \
            "--export=ALL,RUN=$RUN,KIT_DIR=$KIT_DIR,KIT_GEN=$KIT_GEN" "$KIT_GEN/jobs/$step.sh") || rc=$?;;
    rjob)
      bash -n "$KIT_GEN/jobs/$step.sh"
      kit_log "submit $step ($desc)${dl:+ after${dl//:/ }}"
      rc=0; jid=$(sbatch_or_dry "$RUN" "${common[@]}" "--account=$OOD_ACCOUNT" "--ntasks=$OOD_CORES" "--time=$OOD_TIME" \
            -o "$RUN/logs/$step.%j.out" -e "$RUN/logs/$step.%j.out" \
            "--export=ALL,RUN=$RUN,KIT_DIR=$KIT_DIR,KIT_GEN=$KIT_GEN" "$KIT_GEN/jobs/$step.sh") || rc=$?;;
    *) kit_die "unknown step kind $kind";;
  esac
  if [ "$rc" -ne 0 ]; then
    if [ $DRY = 1 ]; then DRYFAIL=$((DRYFAIL + 1)); else kit_die "sbatch failed for $step"; fi
  fi
  jid=${jid%%;*}
  JOB[$step]=$jid
  SUBMITTED=$((SUBMITTED + 1))
  printf '%s\t%s\t%s\t%s\t%s\n' "$step" "$jid" "$(kit_ts)" "${dl#:}" "${common[*]}" >> "$JOBS_TSV"
done < "$KIT_GEN/steps.tsv"

if [ $DRY = 1 ]; then
  # the array example is a template (placeholder samples); SLURM only validates its header
  r=$(cd "$RUN" && sbatch --test-only "$KIT_GEN/learner-scripts/cellranger_array_template.sh" 2>&1 | tail -1) || true
  kit_log "array template (episode 2 callout), test-only: $r"
  kit_log "dry run: $SUBMITTED steps checked, $DRYFAIL sbatch --test-only failure(s); nothing submitted"
  [ $DRYFAIL = 0 ]
else
  kit_log "submitted $SUBMITTED step(s); job IDs in $JOBS_TSV"
  kit_log "watch: squeue -u $(id -un) -o '%.10i %.24j %.8T %.10M %R'   then: bash $KIT_DIR/collect.sh --run $(basename "$RUN")"
fi
