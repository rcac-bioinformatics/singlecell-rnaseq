# kit_helpers.R: functions for the kit's EXTRA rows in the R plans (sourced by run_blocks.R into
# an environment attached as "kit", so they never appear in the learner's globalenv). They read
# the learner's objects and write records only; they never change learner state or print.

.rec <- function() Sys.getenv("KIT_RECORDS")
.step <- function() Sys.getenv("KIT_STEP")

kit_anchor <- function(name, value) {
  f <- file.path(.rec(), "anchors", paste0(.step(), ".tsv"))
  if (!file.exists(f)) cat("anchor\tvalue\n", file = f)
  cat(sprintf("%s\t%s\n", name, paste(value, collapse = ",")), file = f, append = TRUE)
  invisible(value)
}

# keep a copy of a handoff RDS just saved by the episode (publish_rds.sh publishes from here)
kit_handoff <- function(file) {
  wd <- get("work_dir", envir = globalenv())
  src <- if (file.exists(file.path(wd, file))) file.path(wd, file) else file
  dst <- file.path(dirname(.rec()), "handoff", .step())
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, file.path(dst, file), overwrite = TRUE, copy.date = TRUE)
  kit_anchor(paste0("handoff_", file), if (ok) file.size(src) else "COPY-FAILED")
}

# which manual label has the highest mean expression of each canonical marker
kit_marker_check <- function(obj, group, prefix) {
  genes <- c("CD3D", "MS4A1", "LYZ", "GNLY", "PPBP", "FCER1A", "CD8A", "FCGR3A")
  genes <- genes[genes %in% rownames(obj)]
  m <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  g <- as.character(obj[[group, drop = TRUE]])
  for (gene in genes) {
    mu <- tapply(as.numeric(m[gene, ]), g, mean)
    kit_anchor(sprintf("%s_marker_%s", prefix, gene), names(mu)[which.max(mu)])
  }
  # labels per cluster, so round 3 can check the hardcoded mapping against the clusters
  tab <- table(cluster = obj$seurat_clusters, label = g)
  kit_table(tab, paste0(prefix, "_cluster_vs_label"))
}

kit_table <- function(tab, name) {
  f <- file.path(.rec(), "anchors", sprintf("%s_%s.tsv", .step(), name))
  utils::write.table(as.data.frame.matrix(tab), f, sep = "\t", quote = FALSE, col.names = NA)
}

kit_top_genes <- function(df, name, n = 20) {
  o <- order(df$p_val_adj, -df$avg_log2FC)
  kit_anchor(name, paste(head(rownames(df)[o], n), collapse = ","))
  up <- df[df$avg_log2FC > 0, , drop = FALSE]
  o <- order(-up$avg_log2FC)
  kit_anchor(paste0(name, "_by_log2fc"), paste(head(rownames(up)[o], n), collapse = ","))
}
