#!/bin/bash
# check_staged.sh: check a staged learner copy against staged_manifest.tsv. Read only.
# Usage: bash check_staged.sh [DIR] [--quick] [--perm group|other]
#   DIR      copy to check (default: STAGED from config.sh); also works on a staging tree
#   --quick  md5 only for files under 50 MB (sizes are always checked once frozen)
#   --perm   who must be able to read it: group (Depot, default STAGED_PERM) or other
# Checks: every manifest row present with the frozen size and md5 (rows with "-" are reported as
# not frozen yet), readable, directories traversable on the way to it, and nothing extra (episode 2
# outputs, job logs, the 3M whitelist). Exit status 1 if anything fails.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"
DIR=$STAGED; args=(); PERM=$STAGED_PERM; next=0
for a in "$@"; do
  if [ $next = 1 ]; then PERM=$a; next=0; continue; fi
  case $a in
    --quick) args+=(--quick);;
    --perm) next=1;;
    -h|--help) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) DIR=$a;;
  esac
done
[ -d "$DIR" ] || kit_die "$DIR is not a directory (or not readable)"
python3 "$KIT_DIR/py/staged.py" check "$DIR" "$KIT_DIR/staged_manifest.tsv" --perm "$PERM" "${args[@]}"
