#!/bin/bash
# submit_guard.sh: submit_all.sh must stop before creating or submitting anything when
# SINGULARITY_BIND or APPTAINER_BIND is set, and stay silent when both are unset. Called by run_all.sh.
set -u
: "${REPO:?}" "${WORK:?}"
G=$WORK/guard; mkdir -p "$G/bin"
# shellcheck disable=SC2016  # $0 belongs to the stub
for c in sbatch squeue sacct; do printf '#!/bin/sh\necho "stub $0: not called in tests" >&2\nexit 1\n' > "$G/bin/$c"; chmod +x "$G/bin/$c"; done
export PATH=$G/bin:$PATH
ok=0; n=0
for v in SINGULARITY_BIND APPTAINER_BIND; do
  n=$((n + 1))
  out=$(env -u SINGULARITY_BIND -u APPTAINER_BIND "$v=/var/opt" KIT_SCRATCH_BASE="$G/base" \
        bash "$REPO/negishi-run/submit_all.sh" --run guardtest --dry-run 2>&1); rc=$?
  if [ $rc -ne 0 ] && grep -q "FAIL: $v is set" <<< "$out" && [ ! -e "$G/base" ]; then ok=$((ok + 1)); echo "  ok: $v set -> exit $rc, nothing created"
  else echo "  FAILED: $v set -> exit $rc: $(tail -1 <<< "$out")"; fi
done
n=$((n + 1))
out=$(env -u SINGULARITY_BIND -u APPTAINER_BIND KIT_SCRATCH_BASE="$G/base" bash "$REPO/negishi-run/submit_all.sh" --run guardtest --dry-run 2>&1)
if ! grep -q "is set in this shell" <<< "$out"; then ok=$((ok + 1)); echo "  ok: both unset -> no guard message"; else echo "  FAILED: guard fired with both unset"; fi
# --force reruns the steps from --from on, but a dependency that already passed still counts as done
n=$((n + 1)); FB=$G/forcebase; FR=$FB/runs/forcetest
out=$(export KIT_DIR=$REPO/negishi-run KIT_GEN=$REPO/negishi-run/generated KIT_ALLOW_ANY_ROOT=1 KIT_SCRATCH_BASE=$FB
      # shellcheck source=negishi-run/lib.sh disable=SC1091
      source "$REPO/negishi-run/lib.sh" && kit_init_run "$FR"
      for s in r03 r04 r05 r06 r07; do echo "PASS 2026-10-22T00:00:00" > "$FR/markers/$s.done"; done
      env -u SINGULARITY_BIND -u APPTAINER_BIND bash "$REPO/negishi-run/submit_all.sh" --run forcetest --only r --from r06 --force --dry-run 2>&1)
if grep -q "submit r06" <<< "$out" && ! grep -q "r06 depends on r05" <<< "$out"; then ok=$((ok + 1)); echo "  ok: --from r06 --force accepts the passed r05"
else echo "  FAILED: --from r06 --force: $(grep -m1 'depends on\|ERROR' <<< "$out")"; fi
echo "submit guard: $ok of $n"
[ $ok -eq $n ]
