# env_egress.R: the R environment an episode session gets, and every network call the episodes
# make, tested one by one from this (compute-node) session with a fresh HOME and user library.
# Writes $KIT_RECORDS/env/r-env.tsv, env/r-sessioninfo.txt, egress/egress.tsv, egress/proxy-cache.tsv.
# Exit status 0 if every call passed, 2 otherwise (the job reports WARN; nothing depends on it).
rec <- Sys.getenv("KIT_RECORDS"); gen <- Sys.getenv("KIT_GEN")
dir.create(file.path(rec, "env"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(rec, "egress"), recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)

pk <- readLines(file.path(gen, "r_packages.txt"))
ver <- vapply(pk, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"), "")
loc <- vapply(pk, function(p) tryCatch(dirname(find.package(p)), error = function(e) "-"), "")
env <- rbind(
  data.frame(item = c("R", "Bioconductor", "platform", "libPaths", "nodename"),
             value = c(R.version.string,
                       tryCatch(as.character(BiocManager::version()), error = function(e) "?"),
                       R.version$platform, paste(.libPaths(), collapse = ":"), Sys.info()[["nodename"]]),
             where = "-"),
  data.frame(item = paste0("pkg:", pk), value = ver, where = loc))
write.table(env, file.path(rec, "env", "r-env.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(rec, "env", "r-sessioninfo.txt"))
cat("packages:", sum(ver != "MISSING"), "of", length(ver), "present; missing:", paste(pk[ver == "MISSING"], collapse = ", "), "\n")

# proxy variables and cache locations the downloads will use
pv <- c("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "no_proxy", "NO_PROXY", "ftp_proxy",
        "EXPERIMENT_HUB_CACHE", "ANNOTATION_HUB_CACHE", "R_USER_CACHE_DIR", "XDG_CACHE_HOME", "HOME", "R_LIBS_USER", "TMPDIR")
cache <- c(
  ExperimentHub = tryCatch(ExperimentHub::getExperimentHubOption("CACHE"), error = function(e) paste("?", conditionMessage(e))),
  AnnotationHub = tryCatch(AnnotationHub::getAnnotationHubOption("CACHE"), error = function(e) paste("?", conditionMessage(e))),
  BiocFileCache_default = tryCatch(tools::R_user_dir("BiocFileCache", which = "cache"), error = function(e) "?"),
  ExperimentHub_proxy = tryCatch(format(ExperimentHub::getExperimentHubOption("PROXY")), error = function(e) "?"))
write.table(rbind(data.frame(item = pv, value = Sys.getenv(pv)), data.frame(item = paste0("cache:", names(cache)), value = cache)),
            file.path(rec, "egress", "proxy-cache.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ef <- file.path(rec, "egress", "egress.tsv")
cat("call\tstatus\tseconds\tdetail\n", file = ef)
nfail <- 0L
probe <- function(name, expr) {
  t <- Sys.time()
  r <- tryCatch({ v <- force(expr); list("PASS", if (is.character(v)) v[1] else class(v)[1]) },
                error = function(e) list("FAIL", conditionMessage(e)))
  if (r[[1]] == "FAIL") nfail <<- nfail + 1L
  s <- sprintf("%.1f", as.numeric(difftime(Sys.time(), t, units = "secs")))
  cat(sprintf("%s\t%s\t%s\t%s\n", name, r[[1]], s, gsub("[\t\r\n]+", " | ", r[[2]])), file = ef, append = TRUE)
  cat(sprintf("%-45s %s %6ss %s\n", name, r[[1]], s, substr(gsub("[\r\n]+", " ", r[[2]]), 1, 120)))
}
url_ok <- function(u) {
  tf <- tempfile(); on.exit(unlink(tf))
  utils::download.file(u, tf, quiet = TRUE, mode = "wb")
  sprintf("%d bytes", file.size(tf))
}
# plain HTTPS to each host the episodes reach
for (u in c("https://experimenthub.bioconductor.org/metadata/experimenthub.sqlite3",
            "https://bioconductor.org/config.yaml",
            "https://seurat.nygenome.org/src/contrib/PACKAGES",
            "https://raw.githubusercontent.com/satijalab/seurat-data/master/DESCRIPTION",
            "https://rest.kegg.jp/info/kegg",
            "https://cloud.r-project.org/src/contrib/PACKAGES.rds"))
  probe(paste("url", sub("^https://([^/]+)/.*", "\\1", u)), url_ok(u))
# the episode calls themselves (fresh caches: HOME and R_LIBS_USER are new for this job)
probe("celldex::MonacoImmuneData()", { x <- celldex::MonacoImmuneData(); sprintf("%d x %d", nrow(x), ncol(x)) })
probe("celldex::HumanPrimaryCellAtlasData()", { x <- celldex::HumanPrimaryCellAtlasData(); sprintf("%d x %d", nrow(x), ncol(x)) })
probe("SeuratData::InstallData(\"ifnb\")", {
  suppressPackageStartupMessages(library(SeuratData))
  SeuratData::InstallData("ifnb")
  x <- SeuratData::LoadData("ifnb"); sprintf("%d cells", ncol(x)) })
isg <- c("9636", "3434", "3437", "2537", "4599", "91543", "4938", "3669", "3433", "10964", "6772", "3665")
probe("clusterProfiler::enrichGO(ISGs, BP)", {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  g <- clusterProfiler::enrichGO(isg, OrgDb = org.Hs.eg.db, ont = "BP")
  paste(head(as.data.frame(g)$Description, 3), collapse = "; ") })
probe("clusterProfiler::enrichKEGG(ISGs, hsa)", {
  k <- clusterProfiler::enrichKEGG(isg, organism = "hsa")
  paste(head(as.data.frame(k)$Description, 3), collapse = "; ") })
probe("BiocManager::available() (episode 4 install chunk)", { a <- BiocManager::available("^glmGamPoi$"); paste(a, collapse = ",") })
cat("egress:", nfail, "call(s) failed; details in", ef, "\n")
quit(status = if (nfail > 0) 2L else 0L, save = "no")
