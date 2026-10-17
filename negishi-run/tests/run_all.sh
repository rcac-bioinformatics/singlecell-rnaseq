#!/bin/bash
# run_all.sh: the kit's laptop self-tests (about 1 minute; needs bash, python3, rsync; no SLURM, no data).
# Usage: bash negishi-run/tests/run_all.sh
#   1. mutate_kit.sh     edits to the lesson and the kit that build_kit.py (--check) must catch
#   2. submit_guard.sh   submit_all.sh refuses a shell with SINGULARITY_BIND or APPTAINER_BIND
#   3. summarize.sh      py/summarize.py on synthetic records: anchors judged, errors and warnings listed
#   4. rparse.sh         every generated R block parses (skipped without Rscript)
# Work goes to a new temporary directory, removed at the end (kept with KEEP_WORK=1).
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
WORK=$(mktemp -d "${TMPDIR:-/tmp}/scrna-kit-tests.XXXXXX") || exit 1
export REPO WORK
cleanup() { if [ "${KEEP_WORK:-0}" = 1 ]; then echo "work directory kept: $WORK"; else rm -rf "$WORK"; fi; }
trap cleanup EXIT
D=$REPO/negishi-run/tests
fail=0; n=0
for t in mutate_kit.sh submit_guard.sh summarize.sh rparse.sh; do
  n=$((n + 1)); echo "=== $t"; start=$SECONDS
  if bash "$D/$t" > "$WORK/$t.log" 2>&1; then r=PASS; else r=FAIL; fail=$((fail + 1)); fi
  if [ $r = FAIL ] || [ "${VERBOSE:-0}" = 1 ]; then cat "$WORK/$t.log"; else grep -E 'MISSED|FAILED|NOT APPLIED|SKIPPED' "$WORK/$t.log"; tail -1 "$WORK/$t.log"; fi
  echo "--- $t: $r ($((SECONDS - start)) s)"
done
if [ $fail -eq 0 ]; then echo "run_all: all $n suites pass"; else echo "run_all: $fail of $n suites FAIL"; fi
[ $fail -eq 0 ]
