#!/bin/bash
# summarize.sh: py/summarize.py on synthetic records. Checks that anchors are judged both ways,
# that an R error and a warning reach SUMMARY.md verbatim, and that missing steps read "not run".
set -u
: "${REPO:?}" "${WORK:?}"
R=$WORK/sumrun; O=$WORK/sumout
mkdir -p "$R"/records/{anchors,events,metrics/02,console/r03} "$R/markers"
printf 'RUN_ID=synthetic\nKIT_BUILT_FROM=test\n' > "$R/kit.env"
printf 'PASS 2026-10-16T10:00:00\n' > "$R/markers/r03.done"
printf 'FAIL 2026-10-16T10:00:00\n' > "$R/markers/r06.done"
printf 'anchor\tvalue\nep3_cells_after_qc\t11000\nep3_read_genes\t33538\nep3_median_percent_mt_prefilter\t5.8\n' > "$R/records/anchors/r03.tsv"
printf 'anchor\tvalue\nep6_marker_MS4A1\tCD8 T\n' > "$R/records/anchors/c06.tsv"
printf 'episode\tlabel\tline\tcontext\ttype\tmessage\n06\ttop-markers\t181\t-\tERROR\tcould not find function "%%>%%"\n06\tsave\t531\t-\tWARNING\tsome deprecation\n06\tALL\tALL\t\tTIME_S\t600\n' > "$R/records/events/r06.tsv"
printf 'item\tvalue\ncellranger_cells\t11809\nstarsolo_cells\t11000\n' > "$R/records/metrics/02/barcodes.tsv"
printf 'step\tjobid\tsubmitted\tdeps\tflags\n02-cellranger\t11\t-\t\t-\n02-starsolo\t12\t-\t\t-\n' > "$R/jobs.tsv"
printf 'jobid\tstep\tname\tstate\tcpus\twait_s\telapsed_s\tcpu_efficiency\tmaxrss_gib\ttimelimit\texitcode\n11\t02-cellranger\tcellranger\tCOMPLETED\t64\t60\t19447\t0.06\t75.18\t1-00:00:00\t0:0\n12\t02-starsolo\tstarsolo\tFAILED\t64\t60\t16\t0\t0\t1-00:00:00\t1:0\n' > "$R/records/sacct_summary.tsv"
KIT_DIR=$REPO/negishi-run python3 "$REPO/negishi-run/py/summarize.py" "$R" "$REPO/negishi-run/generated" "$O" || exit 1
S=$O/SUMMARY.md; ok=0; n=0
chk() { n=$((n + 1)); if grep -qF -- "$2" "$S"; then ok=$((ok + 1)); echo "  ok: $1"; else echo "  FAILED: $1 (no '$2')"; fi; }
chk "cells after QC within 10 percent passes" "| 03 | cells after QC | 11310 | rel:0.10 | 11000 | PASS"
chk "wrong gene count fails" "| 03 | genes read | 38606 | exact | 33538 | FAIL"
chk "cold marker check fails" "| B | in:B | - | not run | CD8 T | FAIL"
chk "CR vs STARsolo difference over 3 percent fails" "| 6.85 | FAIL"
chk "R error verbatim" 'could not find function "%>%"'
chk "R warning verbatim" "some deprecation"
chk "unrun step reads not run" "| env-egress | not run"
chk "completed learner job reads PASS from sacct" "| 02-cellranger | PASS | COMPLETED"
chk "failed learner job reads FAIL from sacct" "| 02-starsolo | FAIL | FAILED"
chk "learner job without a job reads not run" "| 02-starindex | not run"
echo "summarize checks: $ok of $n"
[ $ok -eq $n ]
