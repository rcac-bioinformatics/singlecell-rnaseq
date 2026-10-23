# redraw_figures.R <figures_redraw.tsv> <handoff dir> <generated dir> <out dir>
# Redraws committed figures from handoff RDS files with the episode's own chunk code: the episode's
# setup block, any `before` blocks, then the figure's block, printed to a PNG of the committed size;
# `part` picks which top-level expression of that block is the figure (- for the last). Prints one line per
# figure (file, pixel size, identities of the object it was drawn from).
args <- commandArgs(trailingOnly = TRUE)
tab <- utils::read.delim(args[1], comment.char = "#", stringsAsFactors = FALSE)
handoff <- args[2]; gen <- args[3]; out <- args[4]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
block <- function(ep, label) file.path(gen, "blocks", ep, paste0(label, ".R"))
loaded <- character(0); objects <- list()
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  ep <- sprintf("%02d", as.integer(r$episode))
  if (!ep %in% loaded) {
    env <- new.env()
    sys.source(block(ep, "setup"), envir = env)  # the episode's library() calls
    loaded <- c(loaded, ep)
  }
  rds <- file.path(handoff, r$handoff)
  if (is.null(objects[[rds]])) objects[[rds]] <- readRDS(rds)
  env <- new.env()
  assign(sub("_.*", "", basename(r$handoff)), objects[[rds]], envir = env)  # pbmc_annotated.rds -> pbmc
  if (r$before != "-") for (b in strsplit(r$before, ",")[[1]]) sys.source(block(ep, b), envir = env)
  set.seed(42)  # ggrepel label placement is random; a fixed seed makes the redraw byte-reproducible
  vals <- lapply(parse(block(ep, r$label)), eval, envir = env)  # in order, like the episode
  part <- if (is.null(r$part) || r$part %in% c("-", "")) length(vals) else as.integer(r$part)
  p <- vals[[part]]
  f <- file.path(out, r$fig)
  grDevices::png(f, width = r$width, height = r$height, res = r$res, type = "cairo")
  print(p)
  grDevices::dev.off()
  id <- SeuratObject::Idents(objects[[rds]])
  cat(sprintf("%s\t%dx%d\t%d identities: %s\n", r$fig, r$width, r$height, nlevels(id), paste(levels(id), collapse = ", ")))
}
