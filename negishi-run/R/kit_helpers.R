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
  genes <- c("CD3D", "MS4A1", "LYZ", "GNLY", "PPBP", "FCER1A", "CD8A", "FCGR3A", "LILRA4", "CD34")
  genes <- genes[genes %in% rownames(obj)]
  m <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  g <- as.character(obj[[group, drop = TRUE]])
  for (gene in genes) {
    mu <- tapply(as.numeric(m[gene, ]), g, mean)
    kit_anchor(sprintf("%s_marker_%s", prefix, gene), names(mu)[which.max(mu)])
  }
  # labels per cluster, so round 3 can check the hardcoded mapping against the clusters. The
  # episodes relabel the resolution 0.5 clusters; seurat_clusters holds the last FindClusters call
  # (resolution 1.2 after episode 5's multi-resolution chunk).
  cl <- if ("RNA_snn_res.0.5" %in% colnames(obj[[]])) obj$RNA_snn_res.0.5 else obj$seurat_clusters
  kit_table(table(cluster = cl, label = g), paste0(prefix, "_cluster_vs_label"))
  # mean log-normalised expression of a wider marker panel per cluster, to rebuild the mapping
  panel <- c("CD3D", "IL7R", "CCR7", "CD8A", "CD8B", "MS4A1", "CD79A", "LYZ", "CD14", "S100A9",
             "FCGR3A", "MS4A7", "GNLY", "NKG7", "KLRD1", "FCER1A", "CST3", "CLEC10A", "LILRA4",
             "IL3RA", "PPBP", "PF4", "CD34", "HBB")
  panel <- panel[panel %in% rownames(obj)]
  pm <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")[panel, , drop = FALSE]
  means <- t(sapply(split(seq_along(cl), cl), function(i) Matrix::rowMeans(pm[, i, drop = FALSE])))
  kit_table(round(means, 2), paste0(prefix, "_cluster_marker_means"))
}

# per identity: mean log-normalised expression and percent of cells expressing each gene
kit_label_profile <- function(obj, genes, name) {
  genes <- genes[genes %in% rownames(obj)]
  m <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  id <- Idents(obj)
  out <- do.call(rbind, lapply(levels(id), function(l) {
    x <- m[, id == l, drop = FALSE]
    c(n = ncol(x), setNames(round(Matrix::rowMeans(x), 2), paste0(genes, "_mean")),
      setNames(round(100 * Matrix::rowMeans(x > 0), 1), paste0(genes, "_pct")))
  }))
  rownames(out) <- levels(id)
  kit_table(out, name)
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
