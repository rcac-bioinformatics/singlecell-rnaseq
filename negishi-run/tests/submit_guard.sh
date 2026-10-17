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
echo "submit guard: $ok of $n"
[ $ok -eq $n ]
