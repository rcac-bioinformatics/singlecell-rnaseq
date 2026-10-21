# updates.md

Running log of testing findings and revision opportunities for this lesson. Severity tags: `[blocker]` breaks a learner or the live delivery, `[should-fix]` wrong or stale but survivable, `[nice]` polish, `[check]` needs verification before it becomes a finding. Append dated entries; move resolved items to Done with the fixing commit.

## Delivery constraint and condensation plan (2026-10-29)

Front matter totals 445 min; the day in `learners/setup.md` allows about 325 min of instruction. Everything stays published and everything gets tested; every episode is taught, but condensed.

| Episode | Front matter | Live on 2026-10-29 | What moves to self-paced |
|---|---|---|---|
| 1 Introduction | 30 | 20 min | Technology comparison detail |
| 2 Raw data processing | 60 | 45 min, submit-then-move-on | Waiting for Cell Ranger/STARsolo to finish; both challenges |
| 3 Quality control | 60 | 45 min | Challenge 2 (aggressive mt filtering) |
| 4 Normalization and features | 55 | 50 min | SCTransform comparison (challenge 1) |
| 5 Dim reduction and clustering | 60 | 45 min | Clustree walkthrough shortened; challenge 2 |
| 6 Cell type annotation | 60 | 45 min | HPCA vs Monaco comparison (challenge 2) |
| 7 Multi-sample integration | 60 | 40 min | Alternative integration methods; challenge 1 |
| 8 Differential expression | 60 | 35 min | sceasy export; GO challenge if short on time |

Adjust after the timed dry run.

## Open findings

### Static review, 2026-10-16

- `[should-fix]` `instructors/instructor-notes.md` describes a two-day format ("Suggested day 1 / day 2") and a 445 min total; the delivery is one day. Replace that paragraph with the condensation table above or a pointer to it.
- `[should-fix]` `instructors/instructor-notes.md` "STAR index building OOM" advises `--mem=40G`. RCAC allocates memory by core count and the lesson's SLURM scripts correctly omit `--mem`; rewrite the advice as "request enough cores (64 on a cpu node gives the full node's memory)" so the notes do not contradict the scripts.
- `[should-fix]` `learners/setup.md` OOD table: app `RStudio (bioconductor)`, R version `4.4.0-bioconductor`, account "(provided by instructor)". Confirm against the live Negishi OOD form for the Bioconductor release that will be deployed on Oct 29 (Goal 4b target is latest r-bioconductor; Bioconductor 3.22 shipped Oct 2025 on R 4.5) and fill in the account name. Re-take the two `learners/fig/ood_rstudio_*.png` screenshots if the form changed. 2026-10-19: the kit runs with `OOD_SIF=/apps/biocontainers/images/r-bioconductor_3.23-R-4.6.1.sif` and `OOD_ACCOUNT=rcac-rnaseq`, the same R 4.6.1 and Bioconductor 3.23 as the local pre-screen. If the OOD app runs this image on Oct 29, the table needs R `4.6.1` (in whatever form the app shows it) and account `rcac-rnaseq`; confirm with `docs/ood_manual_check.md`.
- `[should-fix]` `learners/setup.md` prerequisites say "Negishi cluster or Scholar (provided)". Episode paths, the OOD URL and the SLURM account are all Negishi-only. Remove Scholar or explain it.
- `[should-fix]` `learners/setup.md` local-install block installs `SeuratData` from GitHub via `remotes`; `episodes/multi-sample-integration.Rmd:92` then calls `InstallData("ifnb")`, which downloads from GitHub at run time. Both need egress from the compute node that RStudio runs on. Pre-install `ifnb.SeuratData` into the module library, or stage the `ifnb` object as an RDS and load that instead, and say so in the episode.
- `[should-fix]` `episodes/cell-type-annotation.Rmd:366,639,643`: `celldex::MonacoImmuneData()` and `HumanPrimaryCellAtlasData()` pull from ExperimentHub at run time. Same egress dependency. Pre-populate the ExperimentHub cache in the module (or set `EXPERIMENT_HUB_CACHE` to a shared Depot path) and verify from a compute node.
- `[should-fix]` `episodes/differential-expression.Rmd:595-600`: the `sceasy` export needs `reticulate` plus a Python `anndata` install. Either confirm it works in the OOD app or mark the section as optional/self-paced with a note that it needs a Python environment.
- `[check]` Episodes 3 to 8 `set-paths` chunks: confirm `work_dir` ends with a trailing slash in every episode, because `readRDS(paste0(work_dir, "pbmc_clustered.rds"))` concatenates without a separator while `saveRDS()` relies on `setwd(work_dir)` having run.
- `[check]` Seurat v5 API drift in the current r-bioconductor stack: `IntegrateLayers(method = CCAIntegration, ...)`, `JoinLayers()`, `AggregateExpression(return.seurat = TRUE)`, `FindMarkers(test.use = "DESeq2")`, `SCTransform(vst.flavor = "v2")`, `ScaleData` on split layers. Record `packageVersion("Seurat")` and `SeuratObject` and quote any deprecation message verbatim.
- `[check]` `instructors/instructor-notes.md` says "Request at least 16 GB of memory for R sessions"; `learners/setup.md` requests 16 cores and no memory. Confirm the OOD form has no memory field, and state the memory that 16 cores implies on a Negishi cpu node.
- `[check]` Depot staging `/depot/workshop/data/scrna_workshop/`: confirm all four subdirectories exist, `rds/` holds all five RDS files, and permissions are group-readable by `workshop`. Confirm the process and lead time for adding registrants to the `workshop` Unix group and SLURM account; on 2026-10-06 group provisioning lagged until the day before.
- `[nice]` `README.md` is one line. Add title, site link, license, delivery history.
- `[nice]` `instructors/instructor-notes.md` "Pre-Workshop Checklist" lists `sceasy` as optional for Episode 8 and `clustree` for Episode 5; cross-check against the actual `library()` calls (Seurat, ggplot2, patchwork, dplyr, clustree, SingleR, celldex, SeuratData, clusterProfiler, org.Hs.eg.db, sceasy) and add `enrichplot` if Episode 8 uses it.
- `[nice]` `config.yaml` `created: '2026-02-24'`; consider adding the delivery dates to `README.md` rather than the config.

### Execution testing (prompt-01-test-materials.md)

- `[blocker]` staged data (Negishi preflight 2026-10-19, `check_staged.sh --quick`, 50 problems): `/depot/workshop/data/scrna_workshop/` has no top-level `filtered_feature_bc_matrix/` and no `rds/`, but has `cellranger_output/pbmc10k/outs/filtered_feature_bc_matrix/`, which is not in the manifest. Effects: (a) since patch 01, episode 3 reads `scrna_workshop/filtered_feature_bc_matrix/`, which a learner copy will not have; (b) the cold-start copies in `rds/` do not exist, so no learner can start episodes 4 to 8 cold; (c) every learner's `rsync` brings in a `cellranger_output/pbmc10k/` directory, so the learner's own `cellranger count --id=pbmc10k` job starts in an existing directory that Cell Ranger did not create, which it is expected to refuse (confirm on Negishi). `01_learner_setup.sh:47` stops on (a) and its precomputed-output guard at line 57 would stop on (c). Proposed fix: publish the matrix at the top level and the five RDS files through the publish scripts, and take `cellranger_output/` out of the staged tree (a Depot removal, which Arun runs himself). Full listing (2026-10-19): the staged tree also has `results/`, an earlier complete run (Cell Ranger job 37306471, STAR index 37298434, STARsolo 37304588): its `filtered_feature_bc_matrix/`, all five RDS files, the whole `starsolo_output/` including `starsolo_pbmc10k_Aligned.sortedByCoord.out.bam`, the 3M whitelist and the job logs. No episode or setup.md mentions `results/`, so learners cannot use it, but every learner's `rsync -avP` copies it, BAM included. There is also an empty top-level `starsolo_output/`, and `reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf` (uncompressed; episode 2 makes it with `gunzip -c`). The `cellranger_output/` matrix has 38,606 genes and 11,809 barcodes, the episode 3 output block. Revised fix: move `results/` and `cellranger_output/` out of the learner tree (Arun), publish the top-level matrix with the new `publish_staged.sh --matrix-from-path` (counts gated on the episode 3 anchors), and publish `rds/` from this kit run with `publish_rds.sh`.
- `[should-fix]` kit (`stage_fetch.sh:16,38`, fixed 2026-10-19): `--submit` sbatched `$0`; SLURM runs a spool copy, so the job looked for `lib.sh` in `/var/spool/slurm/job<id>/` and exited at once (jobs 44486378 and 44486433). It now submits the kit path and passes `KIT_SELF_DIR`. The same submission also used `KIT_ACCOUNT=workshop`; that default is now `rcac-rnaseq`.
- `[check]` staged data, 2026-10-21: FASTQs and refdata fetched from 10x (both ETags verified) and published; the matrix published with `--matrix-from-path` from the backup copy (`scrna_workshop_backup-20261019`). Its `cellranger_output/` and `results/` copies are byte-identical (md5 cebbb56c..., 15e8e610..., c3913a6c...), but the job log (`results/cellranger_37306471.out`) does not state the Cell Ranger version (Martian v4.0.14, STAR 2.7.2a). `02-compare` checks it against the kit's fresh `cellranger/10.0.0` run; if they differ, republish with `--matrix-from-run`. Remaining: `rds/` (from the kit run) and `README.txt`.
- `[check]` staged data, 2026-10-19: Arun is moving the old `/depot/workshop/data/scrna_workshop` aside as a backup and the directory is being rebuilt from the manifest, in the layout `learners/setup.md` already describes (`fastq/`, `reference/refdata-gex-GRCh38-2024-A/`, `filtered_feature_bc_matrix/`, `rds/`) plus a top-level `README.txt` (`negishi-run/staged-top/README.txt`) that says which episode uses each folder and where it came from. Order in `negishi-run/README.md`, "Building a new staged directory". `learners/setup.md` "This will copy" could list the README too (not patched yet).
- `[nice]` kit (`py/staged.py`): `check_staged.sh` does not report `reference/refdata-gex-GRCh38-2024-A/genes/genes.gtf`, because the whole refdata directory is a `tree` entry and anything inside it passes. Proposed fix: once the manifest is frozen with `make_staged_manifest.py`, list the refdata files individually so extras inside it are reported.


Append here: episode, finding, severity, file:line, proposed fix.

#### Round 1, static checks on the laptop, 2026-10-16

Commands: `Rscript -e 'sandpaper::validate_lesson(); sandpaper::build_lesson()'` (sandpaper 0.17.1, R 4.5.0): clean, no warnings. `python3 negishi-run/tools/static_checks.py --links`: report in `negishi-run/records/static/static-report.md` (gitignored). Patches in `updates-patches/` (gitignored), each checked with `git apply --check`; with all six applied, `validate_lesson()` is still clean. The episode output blocks are the reference for every number below (decision 2026-10-16); the old CLAUDE.md anchors were wrong and are corrected in `updates-patches/07-claude-md.patch` (applied).

All episodes

- `[should-fix]` Cold start from the staged `rds/` does not work as described. `learners/setup.md:144` says `rds/` is "for jumping into later episodes", but `readRDS(paste0(work_dir, "<name>.rds"))` reads `scrna_workshop/<name>.rds` and the rsync puts the copies in `scrna_workshop/rds/`. No episode says to copy them up. Fix: `updates-patches/05-cold-start-callouts.patch` (a two-line callout before `load-data` in episodes 4, 5, 6, 8). The kit's cold pass copies `rds/*.rds` into the working directory and records this as an adaptation.
- `[nice]` Episodes 3 to 8 build `work_dir` from `"/scratch/negishi/", Sys.getenv("USER")` (`quality-control.Rmd:58`, `normalization-and-features.Rmd:55`, `dimensionality-reduction-clustering.Rmd:58`, `cell-type-annotation.Rmd:58`, `multi-sample-integration.Rmd:50`, `differential-expression.Rmd:60`) while episode 2 and `learners/setup.md:136` use `${RCAC_SCRATCH}`. Same path for a learner; two spellings to explain. Fix: `Sys.getenv("RCAC_SCRATCH")` in all six set-paths chunks, or one sentence saying they are the same directory. Every `work_dir` ends with `/` (checked).
- `[nice]` Output fences are spelled three ways: ```` ```output ```` (22), ```` ```{output} ```` (12, renders the same), and bare ```` ``` ```` (10, `multi-sample-integration.Rmd:101` onward and `raw-data-processing.Rmd:65,130`), which renders as plain code without the OUTPUT label. Fix: use ```` ```output ```` throughout.
- `[should-fix]` `episodes/fig/FIGURE_DESCRIPTIONS.md` is out of date. It has no entry for 15 referenced code-generated figures: `ep3_vlnplot-thresholds`, `ep4_challenge2`, `ep5_challenge1`, `ep5_challenge2`, `challenge2_plot1`, `ep6_load-data`, `ep6_challenge1`, `ep7_challenge1`, `ep7_challenge2`, and the 7 `ep8_*` files other than the volcano (episode 8 names figures `ep8_umap_split.png` and similar, not `ep8_<chunk-label>.png`). Every line number it cites is stale (for example `vlnplot` "line 240", now 271). It also describes a nonexistent `filename.png`. `cellranger-web-summary.png` and `seurat_object.png` are static figures it does not list. Fix: regenerate the file from the kit's `generated/figmap.tsv` after round 3. Rename the `ep8_*` files to `ep8_<chunk-label>.png` when they are regenerated.
- `[nice]` `episodes/fig/10x-chromium-workflow2.png` and `10x-chromium-workflow3.png` are not referenced. Remove them or move them to `old_plots/`.
- `[nice]` A local `sandpaper::build_lesson()` renders the gitignored root `CLAUDE.md` as `site/docs/CLAUDE.html`. CI does not see it, but it would be published if CLAUDE.md were ever committed. Fix: move it to `.claude/CLAUDE.md`, as in rnaseq-analysis.
- `[check]` `links.md:11,15`: the two 10x Genomics URLs return HTTP 429 to scripted requests; the other 23 external URLs return 200. Check them in a browser. `links.md:11` points at `cell-ranger/latest` while the lesson pins Cell Ranger 10.0.0.
- No smart quotes, typographic dashes, shell or R prompts, or broken line continuations in any R or bash block. Every bash block passes `bash -n`. Front matter matches the `instructors/instructor-notes.md` timing table for all eight episodes (325 + 120 = 445 min). The handoff names match on save and load (`pbmc_filtered`, `pbmc_normalized`, `pbmc_clustered`, `pbmc_annotated`, `ifnb_annotated`). Object and column names are consistent (`pbmc`, `ifnb`, `sce`, `bulk`, `mono`; `percent.mt`, `seurat_clusters`, `stim`, `celltype`). The lesson uses `celltype`, not `cell_type`, in episodes 7 and 8; `celltype.stim` is never created, so it is not checked.

Episode 2 (raw-data-processing)

- `[check]` The same array example requests `--ntasks=16` with `--localmem=64`. If Negishi gives about 2 GB per core, 16 cores is about 32 GB, below the 64 GB Cell Ranger is told it may use. Preflight records the memory per core of a cpu node; adjust either value if it is too low.
- `[check]` The whitelist extraction (`:83-86`) hardcodes `cumulusprod_cellranger:10.0.0.sif` and `/software/cellranger-10.0.0/` under `${BIOC_IMAGE_DIR}`, while the job scripts load the default `cellranger` module. Preflight checks that `BIOC_IMAGE_DIR` is set after `module load biocontainers`, that the SIF exists, and that the default cellranger is 10.0.0.
- `[nice]` The directory layout at `:60-71` shows `fastq/ reference/ cellranger_output/ starsolo_output/`, but the rsync at `:76` also brings `filtered_feature_bc_matrix/` and `rds/`. Add them to the tree.
- `[nice]` `learners/setup.md:141` and CLAUDE.md list the FASTQs as R1 and R2 only. Fix: `updates-patches/06-setup-i1-listing.patch` (CLAUDE.md already updated). The whitelist is dropped from the staged manifest; learners extract it in episode 2.
- Challenge 2 (`:771-800`) values 55.2 and 28.3 percent are a hypothetical colleague's metrics. The kit reports the measured Fraction Reads in Cells and Sequencing Saturation beside them as context, not as pass/fail.

Episode 3 (quality-control)


Episode 6 (cell-type-annotation)

- `[should-fix]` For review: the SingleR vs manual cross-table (`cell-type-annotation.Rmd:442-464`, discussed at `:475-505`) shows all 1,006 manual "CD8 T" cells (cluster 4, `:298`) called "B cells" by SingleR, and all 78 "Platelet" cells (cluster 15, `:309`) called "Dendritic cells". The prose says the CD8 T label "is wrong and should be reconsidered", yet the hardcoded `new.cluster.ids` (`:293-312`) keeps it, and the annotated UMAP and every later figure carry it. Check the cluster 4 markers (MS4A1, CD8A) in the round 2 run. Then either relabel cluster 4 and regenerate the table and figures, or rewrite the discussion. Not tolerated away in `anchors.tsv`.
- `[check]` The 18 hardcoded labels (`:293-312`) only work if resolution 0.5 gives exactly 18 clusters in the same order. `names(new.cluster.ids) <- levels(pbmc)` errors on any other count. `anchors.tsv` allows plus or minus 2 clusters for the cluster anchor, but episode 6 needs exactly 18; the kit reports both.

Episode 7 (multi-sample-integration)

- Covered by the earlier static review: `InstallData("ifnb")` egress (`:92`). The kit tests it inside a compute-node R session.

Episode 8 (differential-expression)

- `[should-fix]` `kegg` (`:511`) calls `enrichKEGG()`, which needs `rest.kegg.jp` from the compute node; this egress was not in the earlier list. Added to the kit's egress check.
- `[nice]` The `pseudo-count` code prints "(padj < 0.05)" (`:284`), but the expected output shows "(p_val_adj < 0.05)" (`:305`). Fix the output block when it is regenerated.
- The anchors are 765 Wilcoxon and 1,455 pseudobulk genes (`:151,305`), not "thousands and hundreds". The episode explains why pseudobulk finds more (pseudo-replicates, `:207-212,307-316`). No finding; CLAUDE.md corrected.

Episode 4 (normalization-and-features)

- `[check]` `install-glmgampoi` (`normalization-and-features.Rmd:540`) runs `BiocManager::install("glmGamPoi")` in the learner's session. That needs Bioconductor egress and a writable user library, and it is pointless if the OOD image already has glmGamPoi. Preflight reports whether the image has it; if so, drop the chunk or wrap it in `if (!requireNamespace(...))`.

#### Round 1, local R pre-screen in Docker, 2026-10-16

Run: `negishi-run/run_R_chain.sh --r current`, continuous pass, episodes 3 to 8 plus `env` and `sceasy`. It ran in Docker `arnstrm2/r-bioconductor:3.23-R-4.6.1` (R 4.6.1, Bioconductor 3.23), with the missing packages in a scratch library. Records are in `~/temp/scrna-kit-local/negishi/kit/scrna-kit/runs/local-2026-10-16/` (laptop only). Input: the 10x public `pbmc_10k_v3` filtered matrix (Cell Ranger 3.0.0, 33,538 genes). The staged matrix comes from Cell Ranger 10 (38,606 genes), so the episode 3 to 6 numbers differ from the lesson by design. Episodes 7 and 8 use the same ifnb data as the lesson. The decision rule: local errors are findings; local warnings are provisional until Negishi confirms them.

- `[should-fix]` sceasy cannot convert a Seurat v5 object. When the kit's `08-sceasy` step gives the chunk the `pbmc` object (`:597`), `sceasy::convertFormat()` 0.0.7 stops with: `The 'slot' argument of 'GetAssayData()' was deprecated in SeuratObject 5.0.0 and is now defunct. Please use the 'layer' argument instead.` Python and anndata were present (anndata 0.10.9), so the cause is sceasy itself. Fix options:
  1. convert the assay first (`pbmc[["RNA"]] <- as(pbmc[["RNA"]], "Assay")`, then `convertFormat`);
  2. replace sceasy with an exporter that supports Seurat v5;
  3. keep the section as reading only.

  Confirm on Negishi with the OOD image's sceasy version.
- `[check]` `FindConservedMarkers()` (`multi-sample-integration.Rmd:354`) stops unless the `metap` package is installed ("Please install the metap package to use FindConservedMarkers"). The local image lacked it, which made `r07` fail. `metap` and `multtest` are now in the kit's derived package list, so preflight checks the OOD image for them. If the image lacks them, `learners/setup.md` and the instructor-notes package checklist must add them.
- `[should-fix]` Episode 6 `find-markers` (`cell-type-annotation.Rmd:94`) took 708 s of the 776 s episode locally. Seurat printed "For a (much!) faster implementation of the Wilcoxon Rank Sum Test ... please install the presto package". presto is on GitHub only (immunogenomics/presto). Check whether the OOD image has it. Without it, `FindAllMarkers` takes about 12 minutes of a 45-minute live slot. Either way, the episode should warn that this step takes a while.
- `[should-fix]` (provisional) Warnings learners will see, verbatim:
  - Episode 4 `install-glmgampoi` (`:540`): "package(s) not installed when version(s) same as or greater than current; use `force = TRUE` to re-install: 'glmGamPoi'". The chunk is pointless when the image has glmGamPoi; see the earlier `[check]`.
  - Episode 5 `run-umap` (`:198`) and episode 7 `umap-before` (`:170`): "The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric ... This message will be shown once per session".
  - Episode 7 `load-ifnb` (`:91`): "Assay RNA changing from Assay to Assay5" (and "from Assay to Assay"), and, on rerun, "The following packages are already installed and will not be reinstalled: ifnb".
  - Episode 8 `convert-ids` (`:412`) and `challenge2` (`:729`): "5.54% of input gene IDs are fail to map...".

  Either mention them in the episodes ("you will see this warning; it is harmless") or suppress them where harmless.
- `[check]` (provisional) Episode 8 GO terms: locally the top 15 BP terms are virus terms plus "response to type II interferon", and no type I interferon term. The `go-dotplot` fig.alt (`differential-expression.Rmd:475`) says "response to type I interferon ... appearing lower in the list". Compare with the Negishi run before regenerating `ep8_go_dotplot.png`.
- `[check]` (provisional) Episode 6 marker check, local matrix: MS4A1 was highest in the cells the hardcoded mapping labels "CD8 T", PPBP in "CD4 T", and FCER1A in "pDC". The cluster numbering depends on the input matrix, so this only shows how fragile the hardcoded 18 labels are. It still matches the lesson's own SingleR table (all "CD8 T" called B cells). The Negishi run's `ep6_cluster_vs_label` table settles it.
- Episodes 7 and 8 reproduce the lesson:
  - ifnb 13,999 cells (6,548 CTRL, 7,451 STIM) and 14 integrated clusters;
  - the same top ISGs (IFIT1, IFIT3, TNFSF10, RSAD2, IFIT2, MX1, CXCL10);
  - Wilcoxon 759 genes (lesson 765) and pseudobulk 1,464 (lesson 1,455);
  - `AggregateExpression` returns a Seurat object;
  - the GO top terms are response to virus and defense response to virus.

  The avg_log2FC values drift in the second decimal (IFIT1 7.17 vs 7.12), so the output blocks need regenerating after the Negishi run, not before.
- Egress from the laptop (`env` step): every call passed. That covers ExperimentHub (Monaco 4.0 s, HPCA 0.8 s), `InstallData("ifnb")` (7.7 s), `enrichGO` (26 s), `enrichKEGG` (12 s, rest.kegg.jp) and BiocManager. Negishi compute nodes are the real test (`env-egress` step).

## Runtime record

Append per-step wall time and MaxRSS from testing (Cell Ranger count, STAR index, STARsolo), RStudio session memory peaks per episode, tool versions and `sessionInfo()`. Compare against the times claimed in the text and the schedule.

### Local pre-screen, 2026-10-16 (laptop Docker, 22 threads; not Negishi)

Versions: R 4.6.1, Bioconductor 3.23, Seurat 5.5.1, SeuratObject 5.4.0, SingleR 2.14.2, celldex 1.22.0, SeuratData 0.2.2.9002 (commit 3e51f44), clusterProfiler 4.20.2, DESeq2 1.52.0, glmGamPoi 1.24.0, clustree 0.5.1, sceasy 0.0.7, reticulate 1.46.0, metap 1.14, anndata 0.10.9 (Python). learners/setup.md still says R 4.4.0-bioconductor; Negishi records the OOD versions in `env-egress`.

| Episode | R compute (all chunks, challenges included) | Peak VmHWM | Slowest chunk |
|---|---|---|---|
| 3 | 0.4 min | 2.1 GiB | - |
| 4 | 2.7 min | 9.6 GiB | `scale-regress` 81 s, `challenge1` (SCTransform v2) 45 s and the 9.6 GiB peak |
| 5 | 0.8 min | 2.1 GiB | - |
| 6 | 12.9 min (with 3 failed chunks) | 4.6 GiB | `find-markers` 708 s (no presto) |
| 7 | 3.0 min | 2.7 GiB | `find-markers-int` 73 s |
| 8 | 1.2 min | 3.5 GiB | - |

## Revision opportunities

Non-blocking improvements worth doing before the next delivery cycle (the tentative 2026-11-12 repeat, and the Spring 2027 workshop in Goal 1b).

## Done

### 2026-10-21: patches 08 and 09 applied to episode 2 (branch prep-2026-10-29); applied, pending Negishi confirmation

- `[should-fix]` episode 2 (`raw-data-processing.Rmd:264,327`): `ml cellranger` loads the default, which on Negishi is now `cellranger/10.1.0` (preflight 2026-10-19). The whitelist step uses the `cumulusprod_cellranger:10.0.0` image and the staged matrix is documented as Cell Ranger 10.0.0, so learners would run a different version from the one in the reference results. `ml star` gives star/2.7.11b, which is fine for now but also unpinned. Proposed fix: `updates-patches/08-ep2-pin-cellranger.patch` (pins `cellranger/10.0.0` in the job script and the array template). Applied with `updates-patches/08-ep2-pin-cellranger.patch`, pending Negishi confirmation.
- `[check]` episode 2 (`raw-data-processing.Rmd:256,323,475,520`): the SLURM headers use `--account=workshop`, but the OOD account for this delivery is `rcac-rnaseq` (2026-10-19). If learners submit under `rcac-rnaseq`, all four headers and the CLAUDE.md rule need the change; if `workshop` is a separate SLURM account they can use, setup.md should name both. The kit submits episode 2 under `KIT_ACCOUNT=workshop` as written, so preflight `sbatch --test-only` and the `02-*` jobs will show whether `aseethar` can use it. 2026-10-19: Arun confirmed jobs run under `rcac-rnaseq`. The kit's own jobs now use it (`KIT_ACCOUNT`); the episode 2 scripts keep the account their headers name. Proposed fix: `updates-patches/09-ep2-account-rcac-rnaseq.patch` (all four headers); after it is applied, set `LESSON_ACCOUNT` in `build_kit.py` and the CLAUDE.md SLURM rule to `rcac-rnaseq`. Applied with `updates-patches/09-ep2-account-rcac-rnaseq.patch`, pending Negishi confirmation.

Resolved items with date and commit.

### 2026-10-16: patches 01 to 04 applied to the episodes (branch prep-2026-10-29); applied, pending Negishi confirmation

- `[blocker]` `quality-control.Rmd:58-62`: `data_dir` points at `cellranger_output/pbmc10k/outs/filtered_feature_bc_matrix/`, episode 2's own Cell Ranger output. On the day, learners' Cell Ranger jobs (hours) will not have finished when episode 3 starts at 10:20, so `Read10X()` fails for every learner. The episode 2 text (`:720-722`) already says episode 3 uses the pre-computed matrix. Fix: `updates-patches/01-ep3-load-staged-matrix.patch` loads `scrna_workshop/filtered_feature_bc_matrix/`, keeps the Cell Ranger path as an optional callout, and moves `setwd(work_dir)` from `load-data` (`:72`) into `set-paths` to match episodes 4 to 8. The kit's continuous pass runs episode 3 as written after `02-cellranger`. The cold pass seeds the Cell Ranger path from the staged matrix (documented in `negishi-run/docs/ADAPTATIONS.md`). Applied with `updates-patches/01-ep3-load-staged-matrix.patch`, pending Negishi confirmation.
- `[blocker]` `cell-type-annotation.Rmd:182-183,218-219` use `%>%` with `dplyr::group_by()` and `dplyr::slice_max()`, but the setup chunk (`:43-49`) loads only Seurat, ggplot2, patchwork, SingleR and celldex. Verified locally (Seurat 5.5.1): `exists("%>%")` is FALSE after those five `library()` calls, so `top-markers` and `heatmap` fail in a fresh RStudio session. Fix: `updates-patches/02-ep6-load-dplyr.patch`. Applied with `updates-patches/02-ep6-load-dplyr.patch`, pending Negishi confirmation.
- `[blocker]` Episode 6 confirmed. `top-markers` (`cell-type-annotation.Rmd:181`) and `heatmap` (`:217`) fail with `could not find function "%>%"`, and `dotplot` (`:237`) then fails with `object 'top3' not found`. So the heatmap and dot plot figures are not drawn. Fix: patch 02. Applied with `updates-patches/02-ep6-load-dplyr.patch`, pending Negishi confirmation.
- `[blocker]` `differential-expression.Rmd:572-603`: `export-csv`, `export-dense` and `convert-anndata` use `pbmc`, which episode 8 never creates (it loads only `ifnb`). `writeMM()` comes from Matrix, which is not attached (verified locally: `exists("writeMM")` is FALSE). Fix: `updates-patches/03-ep8-exports.patch` loads `pbmc_annotated.rds` first and calls `Matrix::writeMM()`. Applied with `updates-patches/03-ep8-exports.patch`, pending Negishi confirmation.
- `[should-fix]` `save-rds` (`:564`) overwrites `ifnb_annotated.rds`, the episode 7 to 8 handoff file and the name of the staged `rds/` copy. Fix in patch 03: save as `ifnb_de.rds`. Applied with `updates-patches/03-ep8-exports.patch`, pending Negishi confirmation.
- `[should-fix]` `export-dense` (`:587-590`) writes a dense 29,155 x 11,310 CSV (about 330 million values) despite the "for smaller datasets" prose. Fix in patch 03: limit it to 500 cells. Applied with `updates-patches/03-ep8-exports.patch`, pending Negishi confirmation.
- `[blocker]` Episode 8 confirmed. `export-csv` (`differential-expression.Rmd:572`): `object 'pbmc' not found`. `export-dense` (`:587`): same. Fix: patch 03. Applied with `updates-patches/03-ep8-exports.patch`, pending Negishi confirmation.
- `[check]` `episodes/raw-data-processing.Rmd:314-319`: the Cell Ranger array job header has `--account=workshop` but no `--partition=cpu`, unlike the other three scripts. Confirm the default partition is acceptable or add the line for consistency. Applied with `updates-patches/04-ep2-consistency.patch`, pending Negishi confirmation.
- `[should-fix]` `raw-data-processing.Rmd:139-140` says "four files total: two lanes times two reads"; the expected listing at `:100-105` shows six (I1 included). Fix in `updates-patches/04-ep2-consistency.patch`. Applied with `updates-patches/04-ep2-consistency.patch`, pending Negishi confirmation.
- `[should-fix]` Runtime claims contradict each other. `:288` says "30--60 minutes on 16 cores"; the script requests 64 cores (`:251`) and the callout at `:624` says about 7 hours on 64 cores. `:498` says the STAR index takes "about 30 minutes on 16 cores"; that script requests 64 too. Fix in patch 04. Round 3 replaces the numbers with the measured Elapsed times. Applied with `updates-patches/04-ep2-consistency.patch`, pending Negishi confirmation.
- `[should-fix]` `--localmem=100` in the script (`:276`) vs `64` in the parameter table (`:300`). Fix in patch 04 (the table now says 100 and that it must stay below the memory SLURM gives the job). Applied with `updates-patches/04-ep2-consistency.patch`, pending Negishi confirmation.
- `[should-fix]` Cell Ranger array example (`:312-337`): no `--partition=cpu` (the other three scripts have it), and the sample names WT_rep1 and so on do not exist in the data. Fix in patch 04: add the partition and label it an illustrative template. The kit checks it with `sbatch --test-only` only. Applied with `updates-patches/04-ep2-consistency.patch`, pending Negishi confirmation.
