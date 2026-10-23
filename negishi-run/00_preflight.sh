#!/bin/bash
# 00_preflight.sh: check Negishi before a kit run. Login node, about 3 to 10 minutes (most of it
# loading R in the configured environment and waiting for a 5-minute network test job).
# Usage: bash 00_preflight.sh [--no-srun]
# Writes a PASS/WARN/FAIL/INFO table to the screen and to $KIT_SCRATCH_BASE/state/preflight-<time>.txt,
# plus state/ood_env.sh (R library setting the R jobs use) and state/pkgcheck.txt.
# Exit status 1 if any check FAILs. Read only outside $KIT_SCRATCH_BASE/state.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

DO_SRUN=1
[ "${1:-}" = "--no-srun" ] && DO_SRUN=0
STATE="$KIT_SCRATCH_BASE/state"
mkdir -p "$STATE"
OUT="$STATE/preflight-$(date +%Y%m%d-%H%M%S).txt"
NFAIL=0; NWARN=0
row() {  # status check detail
  local s=$1; shift
  [ "$s" = FAIL ] && NFAIL=$((NFAIL + 1))
  [ "$s" = WARN ] && NWARN=$((NWARN + 1))
  printf '%-5s %-38s %s\n' "$s" "$1" "${2:-}" | tee -a "$OUT"
}
section() { printf '\n== %s ==\n' "$1" | tee -a "$OUT"; }
T0=$(date +%s)
printf 'scRNA-seq kit preflight, %s, user %s, host %s\n' "$(kit_ts)" "$(id -un)" "$(hostname)" | tee "$OUT"

# ---------------------------------------------------------------- kit and tools
section "Kit"
if (cd "$KIT_DIR/.." && sha256sum --quiet -c "$KIT_GEN/episodes.sha256" >/dev/null 2>&1); then
  row PASS "episodes match the kit build" "$(wc -l < "$KIT_GEN/episodes.sha256") files"
else
  row WARN "episodes differ from the kit build" "rerun: python3 $KIT_DIR/build_kit.py (episodes changed after the kit was generated)"
fi
if python3 "$KIT_DIR/build_kit.py" --check >/dev/null 2>&1; then row PASS "build_kit.py --check" "generated/ matches the episodes"
else row FAIL "build_kit.py --check" "generated/ is stale; python3 $KIT_DIR/build_kit.py"; fi
for t in python3 sha256sum sbatch sacct sacctmgr squeue sinfo srun rsync curl; do
  if command -v "$t" >/dev/null 2>&1; then row PASS "command $t" "$(command -v "$t")"; else row FAIL "command $t" "not found"; fi
done
if APP=$(kit_apptainer); then row PASS "apptainer/singularity" "$APP: $("$APP" --version 2>&1 | head -1)"; else row FAIL "apptainer/singularity" "not found"; APP=; fi
if bash --noprofile --norc -c 'type module' >/dev/null 2>&1; then
  row PASS "module function exported" "available in a non-login child shell"
elif [ -n "${LMOD_PKG:-}" ] && [ -f "$LMOD_PKG/init/bash" ]; then
  row PASS "module function" "not exported; sessions source \$LMOD_PKG/init/bash"
else
  row WARN "module function" "not exported and LMOD_PKG unset; sessions fall back to /etc/profile.d/lmod.sh"
fi
for v in SCRATCH RCAC_SCRATCH; do
  if [ "${!v:-}" = "/scratch/negishi/$(id -un)" ]; then row PASS "\$$v" "${!v}"
  else row WARN "\$$v" "'${!v:-unset}' (the R episodes assume /scratch/negishi/\$USER, episode 2 \$RCAC_SCRATCH)"; fi
done
for v in SINGULARITY_BIND APPTAINER_BIND; do
  if [ -z "${!v:-}" ]; then row PASS "\$$v unset" ""; else row FAIL "\$$v is set" "${!v} (learners do not have it; unset before submitting)"; fi
done
case "$KIT_SCRATCH_BASE" in /scratch/negishi/*) row PASS "KIT_SCRATCH_BASE" "$KIT_SCRATCH_BASE";;
  *) row FAIL "KIT_SCRATCH_BASE" "$KIT_SCRATCH_BASE must be under /scratch/negishi/";; esac
row INFO "KIT_R_MODE" "$KIT_R_MODE (image: OOD_SIF; module: R_MODULE)"
case "$OOD_ACCOUNT" in *FILL-IN*) row FAIL "OOD_ACCOUNT" "placeholder; set it in config.sh (the account on the OOD form)";; *) row PASS "OOD_ACCOUNT" "$OOD_ACCOUNT";; esac
if [ "$KIT_R_MODE" = image ]; then case "$OOD_SIF" in *FILL-IN*) row FAIL "OOD_SIF" "placeholder; set it in config.sh";; esac; fi
if [ "$KIT_R_MODE" = module ]; then case "$R_MODULE" in *FILL-IN*) row FAIL "R_MODULE" "placeholder; set it in config.sh";; esac; fi

# ---------------------------------------------------------------- modules and the episode 2 container
section "Modules (ml biocontainers, then the module string episode 2 loads; EXPECTED_MODULES)"
declare -A VCMD=( [cellranger]="cellranger --version" [star]="STAR --version" )
for mv in $EXPECTED_MODULES; do
  tool=${mv%%/*}; want=""; case $mv in */*) want=${mv#*/};; esac
  res=$(timeout 180 bash -c "module load biocontainers >/dev/null 2>&1; module load $mv >/dev/null 2>&1 || exit 3; \
        module -t list 2>&1 | grep -E '^$tool/' | head -1; ${VCMD[$tool]:-true} 2>&1 | grep -v -i warning | head -2 | tr '\n' ' '" 2>&1) || {
    row FAIL "module $tool" "did not load"; continue; }
  loaded=$(printf '%s\n' "$res" | head -1)
  if [ -z "$want" ] || [ "$loaded" = "$mv" ]; then row PASS "module $tool" "$loaded; $(printf '%s\n' "$res" | sed -n 2p)"
  else row FAIL "module $tool" "loaded '$loaded', expected $mv"; fi
done
bid=$(bash -c 'module load biocontainers >/dev/null 2>&1; printf %s "${BIOC_IMAGE_DIR:-}"' 2>/dev/null || true)
if [ -z "$bid" ]; then row FAIL "BIOC_IMAGE_DIR" "unset after module load biocontainers (episode 2 whitelist step uses it)"
elif [ -r "$bid/cumulusprod_cellranger:10.0.0.sif" ]; then row PASS "whitelist container" "$bid/cumulusprod_cellranger:10.0.0.sif"
else row FAIL "whitelist container" "$bid/cumulusprod_cellranger:10.0.0.sif not readable"; fi

# ---------------------------------------------------------------- SLURM
section "SLURM accounts, partition, queue"
for acct in $(printf '%s\n' "$KIT_ACCOUNT" "$OOD_ACCOUNT" | sort -u); do
  case $acct in *FILL-IN*) continue;; esac
  assoc=$(sacctmgr -nP show assoc where user="$(id -un)" account="$acct" format=account,partition,qos 2>/dev/null | head -3 | tr '\n' ' ')
  if [ -n "$assoc" ]; then row PASS "member of $acct" "$assoc"; else row FAIL "member of $acct" "no association for $(id -un)"; fi
done
part=$(scontrol show partition "$KIT_PARTITION" 2>/dev/null | tr ' ' '\n' | grep -E '^(DefMemPerCPU|MaxMemPerCPU|DefMemPerNode|MaxTime|TotalCPUs|TotalNodes|State)=' | tr '\n' ' ')
row INFO "partition $KIT_PARTITION" "${part:-scontrol failed}"
node=$(sinfo -h -p "$KIT_PARTITION" -o '%c %m' 2>/dev/null | sort -u | head -3 | tr '\n' ';')
row INFO "cpu node cores and memory (MB)" "${node:-?} (episode 2 array: 16 cores with --localmem=64)"
pend=$(squeue -h -p "$KIT_PARTITION" -t PD 2>/dev/null | wc -l); runn=$(squeue -h -p "$KIT_PARTITION" -t R 2>/dev/null | wc -l)
idle=$(sinfo -h -p "$KIT_PARTITION" -t idle -o '%D' 2>/dev/null | paste -sd+ | bc 2>/dev/null || echo '?')
row INFO "queue depth $KIT_PARTITION" "$pend pending, $runn running, $idle idle nodes at $(date +%H:%M)"
TMPJ=$(mktemp -d "${TMPDIR:-/tmp}/kitpf.XXXXXX")
for f in "$KIT_GEN"/learner-scripts/*.sh; do
  r=$(cd "$TMPJ" && sbatch --test-only "$f" 2>&1 | tail -1) || true
  case "$r" in *"to start at"*) row PASS "test-only $(basename "$f")" "${r#sbatch: }";;
    *) row FAIL "test-only $(basename "$f")" "$r";; esac
done
for f in "$KIT_GEN"/jobs/*.sh; do
  acct=$KIT_ACCOUNT; extra=()
  case $(basename "$f") in 02-*) ;; *) acct=$OOD_ACCOUNT; extra=("--ntasks=$OOD_CORES" "--time=$OOD_TIME");; esac
  case $acct in *FILL-IN*) row WARN "test-only $(basename "$f")" "skipped: account placeholder"; continue;; esac
  r=$(cd "$TMPJ" && sbatch --test-only "--account=$acct" "${extra[@]}" "$f" 2>&1 | tail -1) || true
  case "$r" in *"to start at"*) row PASS "test-only $(basename "$f")" "${r#sbatch: }";;
    *) row FAIL "test-only $(basename "$f")" "$r";; esac
done
rmdir "$TMPJ" 2>/dev/null || true

# ---------------------------------------------------------------- scratch space
section "Scratch space"
staged_gb=$(du -s --block-size=1G "$STAGED" 2>/dev/null | cut -f1 || echo 0)
need=$(( staged_gb + 150 ))
avail=$(df --output=avail -B1G "$(dirname "$KIT_SCRATCH_BASE")" 2>/dev/null | tail -1 | tr -d ' ' || echo 0)
detail="staged ${staged_gb} GB; one run needs about ${need} GB (copy, Cell Ranger outs with BAM, STAR index, STARsolo BAM, cold copies)"
command -v myquota >/dev/null 2>&1 && row INFO "myquota" "$(myquota 2>/dev/null | grep -i scratch | head -1 | tr -s ' ')"
if [ "${avail:-0}" -ge "$need" ]; then row PASS "filesystem free space" "${avail} GB free; $detail"
else row WARN "filesystem free space" "${avail} GB free; $detail (check your quota with myquota)"; fi

# ---------------------------------------------------------------- staged data
section "Staged data ($STAGED)"
if bash "$KIT_DIR/check_staged.sh" --quick > "$STATE/check_staged-quick.txt" 2>&1; then
  row PASS "check_staged.sh --quick" "$(tail -1 "$STATE/check_staged-quick.txt")"
else
  row FAIL "check_staged.sh --quick" "$(tail -1 "$STATE/check_staged-quick.txt"); details: $STATE/check_staged-quick.txt"
  grep -E '^(FAIL|WARN)' "$STATE/check_staged-quick.txt" | head -12 | sed 's/^/      /' | tee -a "$OUT"
fi

# ---------------------------------------------------------------- R environment
section "R for episodes 3 to 8 ($KIT_R_MODE mode)"
: > "$STATE/ood_env.sh"
# shellcheck disable=SC2016  # R code, expanded by R
RCHK='pk <- readLines(Sys.getenv("KIT_PK")); ok <- vapply(pk, requireNamespace, logical(1), quietly = TRUE)
cat("R", as.character(getRversion()), "Bioc", tryCatch(as.character(packageVersion("BiocVersion")), error = function(e) "?"), "\n")
cat("HOSTLIB", any(grepl("host-site-library", .libPaths())), "\n")
cat("LOADED", sum(ok), "of", length(ok), "\n"); cat("MISSING", paste(pk[!ok], collapse = ","), "\n")
for (p in pk[ok]) cat("VERSION", p, as.character(packageVersion(p)), "\n")'
r1=""
if [ "$KIT_R_MODE" = image ]; then
  if [ -r "$OOD_SIF" ] && [ -n "$APP" ]; then
    row PASS "SIF readable" "$OOD_SIF ($(du -h "$OOD_SIF" | cut -f1))"
    binds=(--bind "$KIT_DIR"); [ -n "$OOD_HOST_LIB" ] && binds+=(--bind "$OOD_HOST_LIB:$OOD_HOST_LIB_MOUNT:ro")
    mkdir -p "$STATE/ood-home" "$STATE/empty-rlib"
    r1=$(timeout 600 "$APP" exec --cleanenv --home "$STATE/ood-home" "${binds[@]}" \
         --env "KIT_PK=$KIT_GEN/r_packages.txt,R_LIBS_USER=$STATE/empty-rlib" "$OOD_SIF" Rscript -e "$RCHK" 2>/dev/null || true)
    if [ -n "$OOD_HOST_LIB" ] && ! printf '%s' "$r1" | grep -q "^HOSTLIB TRUE"; then
      r2=$(timeout 600 "$APP" exec --cleanenv --home "$STATE/ood-home" "${binds[@]}" \
           --env "KIT_PK=$KIT_GEN/r_packages.txt,R_LIBS_USER=$STATE/empty-rlib,R_LIBS_SITE=$OOD_HOST_LIB_MOUNT" "$OOD_SIF" Rscript -e "$RCHK" 2>/dev/null || true)
      if printf '%s' "$r2" | grep -q "^HOSTLIB TRUE"; then
        echo "OOD_R_LIBS_SITE=$OOD_HOST_LIB_MOUNT" > "$STATE/ood_env.sh"
        row WARN "host library on .libPaths()" "only with R_LIBS_SITE=$OOD_HOST_LIB_MOUNT; R jobs will set it (check the OOD app does the same)"
        r1=$r2
      else row FAIL "host library on .libPaths()" "$OOD_HOST_LIB_MOUNT not used by R in the image"; fi
    fi
  else row FAIL "SIF readable" "$OOD_SIF"; fi
else
  r1=$(timeout 600 bash -c "module load biocontainers >/dev/null 2>&1; module load $R_MODULE >/dev/null 2>&1 || exit 3
       KIT_PK=$KIT_GEN/r_packages.txt R_LIBS_USER=$STATE/empty-rlib Rscript -e '$RCHK'" 2>/dev/null || true)
  [ -n "$r1" ] || row FAIL "module $R_MODULE" "did not load or R failed"
fi
if [ -n "$r1" ]; then
  rv=$(printf '%s\n' "$r1" | grep '^R ' || true)
  if [ -z "$EXPECTED_R_VERSION" ]; then row INFO "R and Bioconductor" "$rv (record only; learners/setup.md says 4.4.0-bioconductor)"
  else case "$rv" in "R $EXPECTED_R_VERSION Bioc $EXPECTED_BIOC_VERSION"*) row PASS "R and Bioconductor" "$rv";;
    *) row WARN "R and Bioconductor" "'$rv', expected R $EXPECTED_R_VERSION Bioc $EXPECTED_BIOC_VERSION";; esac; fi
  ld=$(printf '%s\n' "$r1" | grep '^LOADED' || echo "LOADED ?"); ms=$(printf '%s\n' "$r1" | grep '^MISSING' | cut -c9- || true)
  row INFO "episode packages" "${ld#LOADED } load; missing: ${ms:- none}"
  # one row per package (generated/r_packages.tsv = episode code + package_requirements.tsv)
  while IFS=$'\t' read -r pkg eps blocking vle reason; do
    [ "$pkg" = package ] && continue
    ver=$(printf '%s\n' "$r1" | awk -v p="$pkg" '$1 == "VERSION" && $2 == p { print $3 }')
    if [ -z "$ver" ]; then
      if [ "$blocking" = yes ]; then row FAIL "package $pkg" "MISSING; blocks episode $eps: $reason"
      else row WARN "package $pkg" "MISSING; episode $eps (not blocking): $reason"; fi
    elif [ "$vle" != - ] && [ "$(printf '%s\n%s\n' "$ver" "$vle" | sort -V | tail -1)" = "$vle" ]; then
      row WARN "package $pkg" "$ver (<= $vle): $reason"
    else
      row PASS "package $pkg" "$ver (episode $eps${blocking:+, blocking=$blocking})"
    fi
  done < "$KIT_GEN/r_packages.tsv"
  printf '%s: %s; R/Bioc: %s; missing: %s\n' "$(date +%F)" "${ld#LOADED }" "${rv#R }" "${ms:-none}" > "$STATE/pkgcheck.txt"
fi

# ---------------------------------------------------------------- compute-node egress (plain HTTPS; per-call R checks run in step env-egress)
section "Network from a compute node (5-minute job; waits at most 4 minutes to start)"
if [ $DO_SRUN = 1 ]; then
  # shellcheck disable=SC2016  # expanded by the bash on the compute node, not here
  EG='env | grep -i _proxy | sed "s/^/PROXY /"; for u in "$@"; do printf "%s %s\n" "$(curl -s -o /dev/null -r 0-0 -L --max-time 30 -w "%{http_code}" "$u" || echo 000)" "$u"; done'
  eg=$(timeout 240 srun -A "$KIT_ACCOUNT" -p "$KIT_PARTITION" -N1 -n1 -t 00:05:00 --job-name=kit-egress \
       bash -c "$EG" _ https://experimenthub.bioconductor.org/ https://bioconductor.org/config.yaml \
       https://seurat.nygenome.org/src/contrib/PACKAGES https://github.com/satijalab/seurat-data \
       https://rest.kegg.jp/info/kegg https://cloud.r-project.org/ https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz 2>&1) || true
  if [ -z "$eg" ] || ! printf '%s' "$eg" | grep -qE '^[0-9]{3} '; then
    row WARN "compute-node network test" "job did not start or failed within 4 min: $(printf '%s' "$eg" | tail -1)"
  else
    printf '%s\n' "$eg" | grep '^PROXY' | sed 's/^PROXY /      proxy: /' | tee -a "$OUT" || true
    while read -r code url; do
      case "$code" in 2??|3??) row PASS "egress $(printf '%s' "$url" | cut -d/ -f3)" "HTTP $code";;
        *) row FAIL "egress $(printf '%s' "$url" | cut -d/ -f3)" "HTTP $code $url";; esac
    done < <(printf '%s\n' "$eg" | grep -E '^[0-9]{3} ')
  fi
else
  row WARN "compute-node network test" "skipped (--no-srun)"
fi

section "Result"
printf 'FAIL %d, WARN %d, %d s. Saved to %s\n' "$NFAIL" "$NWARN" "$(( $(date +%s) - T0 ))" "$OUT" | tee -a "$OUT"
[ "$NFAIL" -eq 0 ]
