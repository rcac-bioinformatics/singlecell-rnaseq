#!/bin/bash
# mutate_kit.sh: each case copies the repository, applies one edit, and runs build_kit.py --check
# (FAIL cases must fail with the expected text). Called by run_all.sh (REPO, WORK).
set -u
: "${REPO:?}" "${WORK:?}"
M=$WORK/mut; mkdir -p "$M"; n=0; caught=0
copy() { rm -rf "$M/r"; rsync -a --exclude site --exclude .git --exclude negishi-run/records --exclude maintenance "$REPO/" "$M/r/"; }
run() {  # name file old new expected-text
  local name=$1 file=$2 old=$3 new=$4 want=$5 count=${6:-1} out rc
  copy
  if ! python3 - "$M/r/$file" "$old" "$new" "$count" <<'PY'
import sys
p, a, b, c = sys.argv[1:]
s = open(p).read()
if a not in s:
    sys.exit(1)
open(p, "w").write(s.replace(a, b) if c == "all" else s.replace(a, b, 1))
PY
  then n=$((n + 1)); echo "NOT APPLIED (pattern not found): $name"; return; fi
  out=$(cd "$M/r" && python3 negishi-run/build_kit.py --check 2>&1); rc=$?
  n=$((n + 1))
  if [ $rc -ne 0 ] && grep -qF -- "$want" <<< "$out"; then caught=$((caught + 1)); echo "caught: $name -> $(grep -m1 -F -- "$want" <<< "$out" | cut -c1-140)"
  else echo "MISSED: $name (exit $rc)"; echo "$out" | tail -3; fi
}
E=episodes
run "learner script gets --mem" $E/raw-data-processing.Rmd "#SBATCH --job-name=cellranger
#SBATCH --account=rcac-rnaseq
#SBATCH --partition=cpu" "#SBATCH --job-name=cellranger
#SBATCH --account=rcac-rnaseq
#SBATCH --mem=100G
#SBATCH --partition=cpu" "sets --mem"
run "learner script loses its partition" $E/raw-data-processing.Rmd "#SBATCH --job-name=star_index
#SBATCH --account=rcac-rnaseq
#SBATCH --partition=cpu" "#SBATCH --job-name=star_index
#SBATCH --account=rcac-rnaseq" "--partition=MISSING"
run "cellranger module unpinned" $E/raw-data-processing.Rmd "ml cellranger/10.0.0" "ml cellranger" "EXPECTED_MODULES"
run "learner script back on the workshop account" $E/raw-data-processing.Rmd "#SBATCH --job-name=star_index
#SBATCH --account=rcac-rnaseq" "#SBATCH --job-name=star_index
#SBATCH --account=workshop" "--account=workshop (expected rcac-rnaseq)"
run "episode 8 reads an RDS no episode saves" $E/differential-expression.Rmd 'readRDS(paste0(work_dir, "pbmc_annotated.rds"))' 'readRDS(paste0(work_dir, "pbmc_annotate.rds"))' "which no earlier episode saves"
run "heredoc renamed" $E/raw-data-processing.Rmd "> \${RCAC_SCRATCH}/scrna_workshop/run_starsolo.sh" "> \${RCAC_SCRATCH}/scrna_workshop/starsolo.sh" "matched 0"
run "staged path changed in setup.md" learners/setup.md "rsync -avP /depot/workshop/data/scrna_workshop/" "rsync -avP /depot/workshop/data/scrna-workshop/" "matched 0"
run "chunk with a kit extra renamed" $E/quality-control.Rmd '```{r filter, eval=FALSE}' '```{r filter-cells, eval=FALSE}' "R chunk 'filter'"
run "??? removed from the skipped template" $E/quality-control.Rmd "???" "250" "no longer contains ???" all
run "episode order changed in config.yaml" config.yaml "- quality-control.Rmd" "- quality-control.Rmd
- quality-control-extra.Rmd" "episode order changed"
run "an R chunk edited (generated/ stale)" $E/normalization-and-features.Rmd "pbmc <- ScaleData(pbmc)" "pbmc <- ScaleData(pbmc, verbose = FALSE)" "differs from a fresh build"
run "generated file edited by hand" negishi-run/generated/rplans/04.tsv "run" "skip:hand edit" "differs from a fresh build"
echo "mutations caught: $caught of $n"
[ $caught -eq $n ]
