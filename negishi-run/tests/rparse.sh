#!/bin/bash
# rparse.sh: every generated R block and extra parses in R (hard wraps or a stray character would not).
set -u
: "${REPO:?}"
command -v Rscript >/dev/null 2>&1 || { echo "SKIPPED: no Rscript"; exit 0; }
Rscript -e 'f <- list.files(commandArgs(TRUE)[1], pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
bad <- 0; f <- f[!vapply(f, function(x) any(grepl("???", readLines(x), fixed = TRUE)), TRUE)]  # fill-in templates, skipped by the plans
for (x in f) tryCatch(parse(x), error = function(e) { bad <<- bad + 1; cat("FAILED:", x, conditionMessage(e), "\n") })
cat("R files parsed:", length(f) - bad, "of", length(f), "\n"); quit(status = bad > 0)' "$REPO/negishi-run/generated"
