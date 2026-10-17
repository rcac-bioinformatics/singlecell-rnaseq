# negishi-run/lib.sh: functions shared by the kit scripts and the generated jobs.
# Sourced by scripts that already run under `set -euo pipefail`.
# shellcheck shell=bash

KIT_DIR=${KIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
export KIT_DIR
# shellcheck source=negishi-run/config.sh
source "$KIT_DIR/config.sh"
KIT_GEN=${KIT_GEN:-$KIT_DIR/$KIT_GEN_DIRNAME}
export KIT_GEN

# The R episodes in the order they run (used by run_R_chain.sh)
export KIT_R_EPS="03 04 05 06 07 08"

kit_ts() { date '+%Y-%m-%dT%H:%M:%S'; }
kit_log() { printf '[%s] %s\n' "$(kit_ts)" "$*" >&2; }
kit_die() { kit_log "ERROR: $*"; exit 1; }

kit_apptainer() {
  if command -v apptainer >/dev/null 2>&1; then command -v apptainer
  elif command -v singularity >/dev/null 2>&1; then command -v singularity
  else return 1; fi
}

# Refuse to run if a learner working directory could be the staged source.
kit_guard() {
  local w=${1:-${W:?W not set}}
  # KIT_ALLOW_ANY_ROOT=1 is only for testing the kit off-cluster
  case "$w" in /scratch/negishi/*) ;; *) [ "${KIT_ALLOW_ANY_ROOT:-0}" = 1 ] || kit_die "W=$w is not under /scratch/negishi/";; esac
  case "$(readlink -m "$w")/" in
    "$(readlink -m "$STAGED")"/*) kit_die "W=$w resolves inside the staged data $STAGED; refusing to run";;
  esac
  [ "$(readlink -m "$w")" != "$(readlink -m "$STAGED")" ] || kit_die "W is the staged directory"
  case "$(readlink -m "$w")" in "$(readlink -m "$KIT_SCRATCH_BASE")"/*) ;;
    *) kit_die "W=$w is not under KIT_SCRATCH_BASE=$KIT_SCRATCH_BASE";; esac
}

# Per-run settings frozen by 01_learner_setup.sh in $RUN/kit.env
kit_load_run() {
  : "${RUN:?RUN not set}"
  [ -f "$RUN/kit.env" ] || kit_die "$RUN/kit.env missing; start with 01_learner_setup.sh"
  # shellcheck source=/dev/null
  source "$RUN/kit.env"
  export RUN RUN_ID W KIT_REAL_USER KIT_LEARNER_SCRATCH KIT_USER_OVERRIDE STAGED
  export SCRATCH="$KIT_LEARNER_SCRATCH" RCAC_SCRATCH="$KIT_LEARNER_SCRATCH"
  kit_guard "$W"
}

# Write $RUN/kit.env for a new run directory
kit_init_run() {
  local run=$1
  case "$run" in /scratch/negishi/*) ;; *) [ "${KIT_ALLOW_ANY_ROOT:-0}" = 1 ] || kit_die "run directory $run must be under /scratch/negishi/";; esac
  mkdir -p "$run"/{scratch,cold,logs,markers,tmp,ood} \
           "$run"/records/{listings,console,events,figures,sessioninfo,anchors,metrics,egress,env}
  {
    printf 'RUN=%q\n' "$run"
    printf 'RUN_ID=%q\n' "$(basename "$run")"
    printf 'W=%q\n' "$run/scratch/scrna_workshop"
    printf 'KIT_REAL_USER=%q\n' "$(id -un)"
    printf 'KIT_LEARNER_SCRATCH=%q\n' "$run/scratch"
    printf 'KIT_USER_OVERRIDE=%q\n' "${run#/scratch/negishi/}/scratch"
    printf 'KIT_DIR=%q\n' "$KIT_DIR"
    printf 'KIT_GEN=%q\n' "$KIT_GEN"
    printf 'KIT_BUILT_FROM=%q\n' "$(sha256sum "$KIT_GEN/episodes.sha256" | cut -c1-16)"
  } > "$run/kit.env"
  cp "$KIT_DIR/config.sh" "$run/config.sh.used"
}

kit_marker() { printf '%s %s\n' "$2" "$(kit_ts)" > "$RUN/markers/$1.done"; }
# done = marker says PASS or WARN (a step that is running or failed says RUNNING or FAIL)
kit_is_done() { [ -f "$RUN/markers/$1.done" ] && grep -qE '^(PASS|WARN) ' "$RUN/markers/$1.done"; }

# Done predicates for the learner's own SLURM scripts, which write no marker.
# Only outputs written after this run started count (copied files keep older mtimes).
fresh() { [ -e "$1" ] && [ "$1" -nt "$RUN/kit.env" ]; }
outputs_done() {
  case $1 in
    02-cellranger) fresh "$W/cellranger_output/pbmc10k/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" &&
                   fresh "$W/cellranger_output/pbmc10k/outs/metrics_summary.csv";;
    02-starindex) fresh "$W/reference/star_index/SA" && fresh "$W/reference/star_index/Genome";;
    02-starsolo) fresh "$W/starsolo_output/starsolo_pbmc10k_Solo.out/Gene/filtered/barcodes.tsv";;
    *) return 1;;
  esac
}
# A step passed: its marker says PASS or WARN, or (learner steps) its outputs are fresh
kit_step_passed() { kit_is_done "$1" || outputs_done "$1"; }

kit_on_exit() {
  local rc=$?
  if [ "$rc" -ne 0 ] && [ -n "${KIT_STEP:-}" ] && [ -n "${RUN:-}" ]; then
    kit_marker "$KIT_STEP" FAIL
    printf '%s\t%s\t%s\t%s\n' "$KIT_STEP" "${SLURM_JOB_ID:-none}" "$(kit_ts)" "FAIL" >> "$RUN/records/step-end.tsv"
  fi
}

kit_job_begin() {
  local name=$1
  kit_load_run
  export KIT_STEP=$name
  # kit-only jobs: Negishi $TMPDIR is node-wide /tmp, so namespace it by job (learner scripts are untouched)
  export TMPDIR="$RUN/tmp/${SLURM_JOB_ID:-$$}"; mkdir -p "$TMPDIR"
  trap kit_on_exit EXIT
  kit_log "step $name begin: job ${SLURM_JOB_ID:-none} on $(hostname) cpus=${SLURM_CPUS_ON_NODE:-?}"
  kit_marker "$name" RUNNING
  printf '%s\t%s\t%s\n' "$name" "${SLURM_JOB_ID:-none}" "$(kit_ts)" >> "$RUN/records/step-start.tsv"
}

kit_job_end() {
  local name=$1 status=${2:-PASS}
  kit_marker "$name" "$status"
  printf '%s\t%s\t%s\t%s\n' "$name" "${SLURM_JOB_ID:-none}" "$(kit_ts)" "$status" >> "$RUN/records/step-end.tsv"
  kit_log "step $name end: $status"
}

# Run a generated learner session in a child bash without -e/-u/pipefail (as a learner's
# shell behaves). Any command that fails is logged as "### KIT-ERR" by the session's ERR
# trap; the step fails if there is at least one. No set -E: with it, the trap also fires
# inside Lmod's shell functions (SIGPIPE rc 141 seen in the variant-calling kit).
kit_run_session() {
  # one log per attempt (session-<step>.<jobid>.log); session-<step>.log is a copy of the latest
  local name=$1 file=$2 log="$RUN/logs/session-$1.${SLURM_JOB_ID:-$$}.log" rc=0
  kit_log "session $name: $file -> $log"
  ( cd "$KIT_LEARNER_SCRATCH" && bash --noprofile --norc "$file" ) > "$log" 2>&1 || rc=$?
  cp "$log" "$RUN/logs/session-$name.log"
  cp "$log" "$RUN/records/console/session-$name.log"
  local nerr
  nerr=$(grep -c '^### KIT-ERR' "$log" || true)
  if [ "$rc" -ne 0 ] || [ "$nerr" -gt 0 ]; then
    kit_log "session $name: exit $rc, $nerr failed command(s); first:"
    grep -m1 '^### KIT-ERR' "$log" >&2 || true
    return 1
  fi
}

# ---- R episodes ------------------------------------------------------------------------
# Working directory for one R pass: continuous = the learner copy; cold = one fresh directory
# per episode (cold episodes run in parallel and must not share files)
kit_r_scratch() {  # pass ep -> prints the scratch directory (the parent of scrna_workshop)
  if [ "$1" = continuous ]; then printf '%s\n' "$KIT_LEARNER_SCRATCH"
  else printf '%s\n' "$RUN/cold/$2/scratch"; fi
}

# Cold start: what a learner who skipped the previous episode needs, copied from $STAGED.
# Documented in docs/ADAPTATIONS.md; the episodes do not say how (updates.md, patch 05).
kit_cold_seed() {  # ep scratch_dir
  local ep=$1 d=$2/scrna_workshop
  mkdir -p "$d"
  case $ep in
    03) # the directory episode 3's set-paths reads (generated/ep3_matrix_path.txt, from the episode)
        local mp; mp=$(cat "$KIT_GEN/ep3_matrix_path.txt"); mp=${mp%/}
        mkdir -p "$d/$(dirname "$mp")"
        rm -rf "${d:?}/$mp"
        cp -r "$STAGED/filtered_feature_bc_matrix" "$d/$mp";;
    04) cp "$STAGED/rds/pbmc_filtered.rds" "$d/";;
    05) cp "$STAGED/rds/pbmc_normalized.rds" "$d/";;
    06) cp "$STAGED/rds/pbmc_clustered.rds" "$d/";;
    07) ;;
    08) cp "$STAGED/rds/ifnb_annotated.rds" "$d/";;
    *) kit_die "kit_cold_seed: unknown episode $ep";;
  esac
  kit_log "cold seed for $ep in $d: $(find "$d" -maxdepth 4 -type f | sed "s#^$d/##" | tr '\n' ' ')"
}

# Run one R script the way the selected R environment runs it. Prints nothing; returns R's status.
#   mode: image | module | current;  user_path: value for USER (test root below /scratch/negishi/)
kit_r() {  # mode name user_path log script [args...]
  local mode=$1 name=$2 user_path=$3 log=$4; shift 4
  local tag="$name.${SLURM_JOB_ID:-$$}" home rlib rc=0
  home="$RUN/ood/home/$tag"; rlib="$RUN/ood/rlib/$tag"
  mkdir -p "$home" "$rlib"
  # image and module mode get a fresh R_LIBS_USER (a new learner); current mode keeps the caller's
  local start envs
  start=$(stat -c %Y "$RUN/kit.env")
  envs="USER=$user_path,KIT_RECORDS=$RUN/records,KIT_STEP=$name,KIT_GEN=$KIT_GEN,KIT_DIR=$KIT_DIR,KIT_RUN_START=$start"
  [ "$mode" = current ] || envs="$envs,R_LIBS_USER=$rlib"
  case $mode in
    image)
      local app; app=$(kit_apptainer) || kit_die "apptainer/singularity not found"
      [ -r "$OOD_SIF" ] || kit_die "OOD_SIF $OOD_SIF not readable (config.sh)"
      local binds=(--bind "$KIT_SCRATCH_BASE" --bind "$KIT_DIR")
      [ -n "$OOD_HOST_LIB" ] && binds+=(--bind "$OOD_HOST_LIB:$OOD_HOST_LIB_MOUNT:ro")
      [ -d "$STAGED" ] && binds+=(--bind "$STAGED:$STAGED:ro")
      local x; for x in ${KIT_EXTRA_BINDS:-}; do binds+=(--bind "$x"); done   # off-cluster testing only
      if [ -f "$KIT_SCRATCH_BASE/state/ood_env.sh" ]; then
        # shellcheck source=/dev/null
        source "$KIT_SCRATCH_BASE/state/ood_env.sh"   # sets OOD_R_LIBS_SITE if preflight found it necessary
      fi
      [ -n "${OOD_R_LIBS_SITE:-}" ] && envs="$envs,R_LIBS_SITE=$OOD_R_LIBS_SITE"
      kit_log "R $name: $app exec --cleanenv $OOD_SIF (USER=$user_path)"
      "$app" exec --cleanenv --home "$home" "${binds[@]}" --env "$envs" "$OOD_SIF" Rscript "$@" > "$log" 2>&1 || rc=$?
      ;;
    module)
      kit_log "R $name: module load $R_MODULE (USER=$user_path)"
      # shellcheck source=/dev/null
      ( type module >/dev/null 2>&1 || source "${LMOD_PKG:-/opt/lmod/lmod}/init/bash" 2>/dev/null || source /etc/profile.d/lmod.sh
        module load biocontainers >/dev/null 2>&1 || true
        # shellcheck disable=SC2086  # R_MODULE may name more than one module
        module load $R_MODULE || exit 97
        export HOME="$home"
        IFS=, read -r -a kv <<< "$envs"; for x in "${kv[@]}"; do export "${x?}"; done
        Rscript "$@" ) > "$log" 2>&1 || rc=$?
      ;;
    current)
      command -v Rscript >/dev/null 2>&1 || kit_die "--r current: no Rscript on PATH (load the R module, or use an OOD RStudio terminal)"
      kit_log "R $name: $(command -v Rscript) from the current shell (USER=$user_path)"
      ( IFS=, read -r -a kv <<< "$envs"; for x in "${kv[@]}"; do export "${x?}"; done
        Rscript "$@" ) > "$log" 2>&1 || rc=$?
      ;;
    *) kit_die "R mode must be image, module or current (got $mode)";;
  esac
  return $rc
}

# One R episode pass: seed (cold), run the episode's plan, record status. Used by the
# generated jobs and by run_R_chain.sh.
kit_r_episode() {  # mode pass ep
  local mode=$1 pass=$2 ep=$3 name sd rc=0
  if [ "$pass" = continuous ]; then name="r$ep"; else name="c$ep"; fi
  sd=$(kit_r_scratch "$pass" "$ep")
  kit_guard "$sd/scrna_workshop"
  [ "$pass" = cold ] && kit_cold_seed "$ep" "$sd"
  local plan="$KIT_GEN/rplans/$ep.tsv"
  [ -f "$plan" ] || kit_die "no R plan $plan"
  local log="$RUN/logs/R-$name.${SLURM_JOB_ID:-$$}.log"
  kit_r "$mode" "$name" "${sd#/scratch/negishi/}" "$log" "$KIT_DIR/R/run_blocks.R" "$plan" || rc=$?
  cp "$log" "$RUN/logs/R-$name.log"; cp "$log" "$RUN/records/console/R-$name.log"
  case $rc in
    0) kit_job_end "$name" PASS ;;
    3) kit_job_end "$name" WARN ;;   # only optional blocks (challenge, solution, spoiler) failed
    *) kit_log "R step $name failed (exit $rc); see $RUN/logs/R-$name.log"; return 1 ;;
  esac
}

# ---- helpers used by "KIT EXTRA" lines in the generated sessions (exported below) ----------
# Episode 2 scripts are written by the learner's heredoc blocks; check they are the episode text
kit_check_written() {  # file expected
  if cmp -s "$1" "$2"; then echo "KIT: $1 matches the episode text"
  else echo "### KIT-ERR $1 differs from the episode heredoc ($2)" >&2; diff "$2" "$1" | head -20 >&2; fi
}

kit_ep2_metrics() {  # copy Cell Ranger and STARsolo metrics, compare barcodes with the staged matrix
  local m="$RUN/records/metrics/02" cr="$W/cellranger_output/pbmc10k/outs" so="$W/starsolo_output"
  mkdir -p "$m"
  local f
  for f in "$cr/metrics_summary.csv" "$cr/web_summary.html" "$so/starsolo_pbmc10k_Log.final.out" \
           "$so/starsolo_pbmc10k_Solo.out/Gene/Summary.csv" "$so/starsolo_pbmc10k_Solo.out/Barcodes.stats" \
           "$so/starsolo_pbmc10k_Solo.out/Gene/Features.stats"; do
    [ -f "$f" ] && cp "$f" "$m/" || echo "### KIT-ERR missing $f" >&2
  done
  local crb="$cr/filtered_feature_bc_matrix/barcodes.tsv.gz" sob="$so/starsolo_pbmc10k_Solo.out/Gene/filtered/barcodes.tsv"
  local stb="$STAGED/filtered_feature_bc_matrix/barcodes.tsv.gz"
  {
    printf 'item\tvalue\n'
    printf 'cellranger_cells\t%s\n' "$(zcat "$crb" 2>/dev/null | wc -l)"
    printf 'cellranger_genes\t%s\n' "$(zcat "$cr/filtered_feature_bc_matrix/features.tsv.gz" 2>/dev/null | wc -l)"
    printf 'starsolo_cells\t%s\n' "$(wc -l < "$sob" 2>/dev/null || echo NA)"
    printf 'starsolo_genes\t%s\n' "$(wc -l < "$so/starsolo_pbmc10k_Solo.out/Gene/filtered/features.tsv" 2>/dev/null || echo NA)"
    printf 'staged_cells\t%s\n' "$(zcat "$stb" 2>/dev/null | wc -l)"
    printf 'shared_barcodes_cr_starsolo\t%s\n' "$(comm -12 <(zcat "$crb" 2>/dev/null | sed 's/-1$//' | sort) <(sed 's/-1$//' "$sob" 2>/dev/null | sort) | wc -l)"
    printf 'fresh_barcodes_md5\t%s\n' "$(md5sum < "$crb" 2>/dev/null | cut -c1-32)"
    printf 'staged_barcodes_md5\t%s\n' "$(md5sum < "$stb" 2>/dev/null | cut -c1-32)"
    printf 'fresh_barcodes_content_md5\t%s\n' "$(zcat "$crb" 2>/dev/null | md5sum | cut -c1-32)"
    printf 'staged_barcodes_content_md5\t%s\n' "$(zcat "$stb" 2>/dev/null | md5sum | cut -c1-32)"
    for f in features.tsv.gz matrix.mtx.gz; do
      printf 'fresh_vs_staged_%s_content_identical\t%s\n' "$f" \
        "$(cmp -s <(zcat "$cr/filtered_feature_bc_matrix/$f" 2>/dev/null) <(zcat "$STAGED/filtered_feature_bc_matrix/$f" 2>/dev/null) && echo yes || echo NO)"
    done
  } > "$m/barcodes.tsv"
  cat "$m/barcodes.tsv"
}

export -f kit_ts kit_log kit_check_written kit_ep2_metrics
