# rds_summary.R <label> <file.rds>...: one row per RDS for publish_rds.sh's side-by-side listing.
# Columns: label file bytes md5 class cells features assays layers reductions idents clusters meta_columns
args <- commandArgs(trailingOnly = TRUE)
label <- args[1]; files <- args[-1]
cat("label\tfile\tbytes\tmd5\tclass\tcells\tfeatures\tassays\tlayers\treductions\tidents\tclusters\tmeta_columns\n")
for (f in files) {
  base <- basename(f)
  if (!file.exists(f)) { cat(label, base, "MISSING", rep("-", 10), sep = "\t"); cat("\n"); next }
  md5 <- unname(tools::md5sum(f))
  o <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(o, "error")) { cat(label, base, file.size(f), md5, paste("readRDS error:", conditionMessage(o)), rep("-", 8), sep = "\t"); cat("\n"); next }
  suppressPackageStartupMessages(requireNamespace("SeuratObject", quietly = TRUE))
  assays <- tryCatch(paste(SeuratObject::Assays(o), collapse = ","), error = function(e) "-")
  layers <- tryCatch(paste(SeuratObject::Layers(o[[SeuratObject::DefaultAssay(o)]]), collapse = ","), error = function(e) "-")
  reds <- tryCatch(paste(SeuratObject::Reductions(o), collapse = ","), error = function(e) "-")
  idents <- tryCatch(sprintf("%d (%s)", nlevels(SeuratObject::Idents(o)), paste(head(levels(SeuratObject::Idents(o)), 20), collapse = ",")), error = function(e) "-")
  cl <- if ("seurat_clusters" %in% colnames(o@meta.data)) nlevels(factor(o$seurat_clusters)) else "-"
  cat(label, base, file.size(f), md5, class(o)[1], ncol(o), nrow(o), assays, layers, reds, idents, cl,
      paste(colnames(o@meta.data), collapse = ","), sep = "\t")
  cat("\n")
}
