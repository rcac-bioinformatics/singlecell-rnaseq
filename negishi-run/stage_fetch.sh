#!/bin/bash
# stage_fetch.sh: fetch from 10x only what the staged copy is missing, into a staging tree on scratch.
# Never writes to $STAGED; publish_staged.sh copies a checked staging tree there.
# Usage: bash stage_fetch.sh [--to DIR] [--only fastqs|refdata] [--force] [--submit] [--dry-run]
#   --to      staging tree (default $KIT_SCRATCH_BASE/staging/scrna_workshop)
#   --only    one source from sources.tsv (default: every 10x source)
#   --force   fetch even if $STAGED already has everything the source provides
#   --submit  run this as a batch job (1 node, 4 tasks, 12 h, KIT_ACCOUNT) instead of here; the
#             FASTQ archive is 52 GB, so do not run it in the foreground on a login node
#   --dry-run say what would be fetched, fetch nothing
# Each archive is verified against sources.tsv (bytes and S3 multipart ETag), then extracted;
# staging/provenance.tsv records URL, bytes, ETag, archive md5 and the md5 of every extracted file.
# The filtered matrix and the five RDS files are not on 10x: see publish_staged.sh --matrix-from-run
# and publish_rds.sh.
set -euo pipefail
KIT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export KIT_DIR
# shellcheck source=negishi-run/lib.sh
source "$KIT_DIR/lib.sh"

TO="$KIT_SCRATCH_BASE/staging/scrna_workshop"; ONLY=""; FORCE=0; SUBMIT=0; DRY=0; ARGS=("$@")
while [ $# -gt 0 ]; do
  case $1 in
    --to) TO=$2; shift 2;;
    --only) ONLY=$2; shift 2;;
    --force) FORCE=1; shift;;
    --submit) SUBMIT=1; shift;;
    --dry-run) DRY=1; shift;;
    -h|--help) sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) kit_die "unknown option $1";;
  esac
done
case "$(readlink -m "$TO")/" in "$(readlink -m "$STAGED")"/*) kit_die "--to $TO is inside \$STAGED; fetch into a staging tree, then publish_staged.sh";; esac
if [ $SUBMIT = 1 ]; then
  mkdir -p "$KIT_SCRATCH_BASE/staging"
  a=(); for x in "${ARGS[@]}"; do [ "$x" = --submit ] || a+=("$x"); done
  sbatch --account="$KIT_ACCOUNT" --partition="$KIT_PARTITION" --nodes=1 --ntasks=4 --time=12:00:00 --job-name=kit-stage-fetch \
    -o "$KIT_SCRATCH_BASE/staging/stage_fetch.%j.out" --export=ALL "$0" "${a[@]}"
  exit 0
fi
DL="$(dirname "$TO")/downloads"
mkdir -p "$TO" "$DL"
PROV="$(dirname "$TO")/provenance.tsv"
[ -f "$PROV" ] || printf 'date\tsource\tpath\tbytes\tmd5\tdetail\n' > "$PROV"
[ -n "${SLURM_JOB_ID:-}" ] && export TMPDIR="$DL/tmp.$SLURM_JOB_ID" && mkdir -p "$TMPDIR"

while IFS=$'\t' read -r id url bytes etag part provides _; do
  case "$id" in id|\#*) continue;; esac
  case "$url" in kit:*) kit_log "$id: not on 10x ($url); skipped"; continue;; esac
  [ -n "$ONLY" ] && [ "$ONLY" != "$id" ] && continue
  # does $STAGED already have every required manifest row this source provides?
  missing=$(awk -F'\t' -v s="10x:$id" -v st="$STAGED" 'index($6, s) == 1 && $5 == "required" { print $1 }' "$KIT_DIR/staged_manifest.tsv" |
            while read -r p; do [ -e "$STAGED/$p" ] || echo "$p"; done)
  if [ -z "$missing" ] && [ $FORCE = 0 ]; then kit_log "$id: \$STAGED has everything it provides ($provides); nothing to fetch"; continue; fi
  kit_log "$id: missing in \$STAGED: $(printf '%s' "$missing" | tr '\n' ' ')"
  f="$DL/$(basename "$url")"
  if [ $DRY = 1 ]; then kit_log "dry run: would fetch $url ($bytes bytes) to $f and extract into $TO/$provides"; continue; fi
  if [ ! -f "$f" ] || [ "$(stat -c %s "$f")" != "$bytes" ]; then
    kit_log "$id: downloading $url"
    curl -fL --retry 5 --retry-delay 30 -C - -o "$f" "$url"
  fi
  [ "$(stat -c %s "$f")" = "$bytes" ] || kit_die "$id: $f is $(stat -c %s "$f") bytes, sources.tsv says $bytes"
  read -r got md5a < <(python3 "$KIT_DIR/py/s3etag.py" "$f" "$part")
  [ "$got" = "$etag" ] || kit_die "$id: ETag $got differs from sources.tsv $etag (corrupt or changed upstream)"
  printf '%s\t%s\t%s\t%s\t%s\tETag %s verified\n' "$(kit_ts)" "$id" "$(basename "$f")" "$bytes" "$md5a" "$etag" >> "$PROV"
  case $id in
    fastqs)
      mkdir -p "$TO/fastq"
      tar -tf "$f" | grep -E '\.fastq\.gz$' > "$DL/fastqs.list"
      tar -xf "$f" -C "$TO/fastq" --strip-components=1 --wildcards '*.fastq.gz';;
    refdata)
      mkdir -p "$TO/reference"
      tar -xzf "$f" -C "$TO/reference";;
    *) kit_die "$id: no extraction rule";;
  esac
  (cd "$TO" && find "${provides%/}" -type f | sort | while read -r p; do
     printf '%s\t%s\t%s\t%s\t%s\textracted from %s\n' "$(kit_ts)" "$id" "$p" "$(stat -c %s "$p")" "$(md5sum < "$p" | cut -c1-32)" "$(basename "$f")"
   done) >> "$PROV"
  kit_log "$id: extracted into $TO/$provides"
done < "$KIT_DIR/sources.tsv"
kit_log "staging tree $TO; provenance $PROV"
kit_log "check it with: bash $KIT_DIR/check_staged.sh $TO --quick   (the matrix and rds/ rows are expected to be missing there)"
kit_log "publish with:  bash $KIT_DIR/publish_staged.sh --from $TO --part fastq|reference --dry-run"
