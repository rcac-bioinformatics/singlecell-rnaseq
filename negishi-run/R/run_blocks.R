# run_blocks.R <plan.tsv>
# Runs the R chunks of one episode, in order, in this one fresh R session, the way a learner
# pastes them into the RStudio console. Per chunk it records the console output (commands echoed,
# printed values, messages and warnings where they occur, errors) to
# $KIT_RECORDS/console/<step>/<episode>_<label>.txt (and to stdout, so R-<step>.log is the whole
# console), plots to $KIT_RECORDS/figures/<step>/<episode>_<label>_NN.png, and events (ERROR,
# WARNING, MESSAGE, TIME_S, VMHWM_KB, VMRSS_KB, SKIPPED, KITERR) to $KIT_RECORDS/events/<step>.tsv.
# An error stops the rest of that chunk (as source() does) and the run continues with the next
# chunk, so later failures are visible too. Rows with action "extra" are kit code (anchors,
# handoff copies): evaluated silently in a child environment, never on the learner console.
# Exit status: 0 all chunks ran, 3 only optional chunks (challenge, solution, spoiler, callout)
# failed, 1 otherwise.
args <- commandArgs(trailingOnly = TRUE)
plan <- read.delim(args[1], colClasses = "character")
rec <- Sys.getenv("KIT_RECORDS"); step <- Sys.getenv("KIT_STEP"); gen <- Sys.getenv("KIT_GEN")
kitdir <- Sys.getenv("KIT_DIR")
condir <- file.path(rec, "console", step); figdir <- file.path(rec, "figures", step)
for (d in c(condir, figdir, file.path(rec, c("events", "sessioninfo", "anchors")))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
evf <- file.path(rec, "events", paste0(step, ".tsv"))
cat("episode\tlabel\tline\tcontext\ttype\tmessage\n", file = evf)
ev <- function(b, type, msg) {
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\n", b$episode, b$label, b$line, b$context, type,
              gsub("[\t\r\n]+", " | ", paste(msg, collapse = " "))), file = evf, append = TRUE)
}
proc_kb <- function(key) {
  x <- tryCatch(grep(paste0("^", key, ":"), readLines("/proc/self/status"), value = TRUE), error = function(e) character())
  if (length(x)) sub("\\s*kB$", "", sub(paste0("^", key, ":\\s*"), "", x)) else "NA"
}

# kit helpers (kit_anchor, kit_handoff, ...) on the search path, not in the learner's globalenv
kit <- new.env()
sys.source(file.path(kitdir, "R", "kit_helpers.R"), envir = kit)
attach(kit, name = "kit", warn.conflicts = FALSE)

cat("R:", R.version.string, "\n")
cat("Bioconductor:", tryCatch(as.character(packageVersion("BiocVersion")), error = function(e) "BiocVersion not installed"), "\n")
cat(".libPaths():\n"); print(.libPaths())
cat("USER =", Sys.getenv("USER"), " HOME =", Sys.getenv("HOME"), " TMPDIR =", Sys.getenv("TMPDIR"), "\n")
cat("cwd =", getwd(), "\n")

t0 <- Sys.time(); n_err_required <- 0L; n_err_optional <- 0L
for (i in seq_len(nrow(plan))) {
  b <- as.list(plan[i, ])
  tag <- sprintf("%s_%s", b$episode, b$label)
  if (startsWith(b$action, "skip")) {
    cat(sprintf("\n######## %s chunk %s SKIPPED by the kit: %s\n", b$episode, b$label, sub("^skip:", "", b$action)))
    ev(b, "SKIPPED", sub("^skip:", "", b$action)); next
  }
  f <- file.path(gen, b$file)
  if (b$action == "extra") {
    r <- tryCatch({ eval(parse(f), envir = new.env(parent = globalenv())); "ok" },
                  error = function(e) conditionMessage(e), warning = function(w) conditionMessage(w))
    if (!identical(r, "ok")) ev(b, "KITERR", r)
    next
  }
  cat(sprintf("\n######## %s chunk %s (Rmd line %s, %s) ########\n", b$episode, b$label, b$line, b$context))
  con <- file(file.path(condir, paste0(tag, ".txt")), open = "wt")
  sink(con, split = TRUE)
  grDevices::png(file.path(figdir, paste0(tag, "_%02d.png")), width = 1400, height = 1000, res = 150)
  tb <- Sys.time()
  failed <- FALSE
  withCallingHandlers(
    tryCatch(
      source(f, echo = TRUE, print.eval = TRUE, max.deparse.length = Inf, local = globalenv(),
             spaced = FALSE, prompt.echo = "> ", continue.echo = "+ "),
      error = function(e) {
        failed <<- TRUE
        cat("Error:", conditionMessage(e), "\n")
        ev(b, "ERROR", conditionMessage(e))
      }),
    warning = function(w) {
      cat("Warning message:\n", conditionMessage(w), "\n", sep = "")
      ev(b, "WARNING", conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cat(conditionMessage(m))
      ev(b, "MESSAGE", conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  while (grDevices::dev.cur() > 1) grDevices::dev.off()
  sink(); close(con)
  if (failed) {
    if (b$optional == "1") n_err_optional <- n_err_optional + 1L else n_err_required <- n_err_required + 1L
  }
  ev(b, "TIME_S", sprintf("%.1f", as.numeric(difftime(Sys.time(), tb, units = "secs"))))
  ev(b, "VMHWM_KB", proc_kb("VmHWM"))
  ev(b, "VMRSS_KB", proc_kb("VmRSS"))
}
# devices opened for chunks that drew nothing leave no file or an empty one
for (p in list.files(figdir, pattern = "\\.png$", full.names = TRUE)) if (file.size(p) < 3000) file.remove(p)
tot <- list(episode = "ALL", label = "ALL", line = "", context = "")
ev(tot, "TIME_S", sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
ev(tot, "VMHWM_KB", proc_kb("VmHWM"))
writeLines(capture.output(sessionInfo()), file.path(rec, "sessioninfo", paste0(step, ".txt")))
cat(sprintf("\nrun_blocks: %s done, %d required and %d optional chunk(s) failed\n", step, n_err_required, n_err_optional))
quit(status = if (n_err_required > 0) 1L else if (n_err_optional > 0) 3L else 0L, save = "no")
