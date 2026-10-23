---
title: Instructor Notes
---

## Workshop Timing

| Episode | Title | Teaching | Exercises | Total |
|:-------:|:------|:--------:|:---------:|:-----:|
| 1 | Introduction to Single-Cell RNA-Seq | 30 min | 0 min | 30 min |
| 2 | Raw Data Processing | 45 min | 15 min | 60 min |
| 3 | Quality Control | 40 min | 20 min | 60 min |
| 4 | Normalization and Feature Selection | 40 min | 15 min | 55 min |
| 5 | Dimensionality Reduction and Clustering | 45 min | 15 min | 60 min |
| 6 | Cell Type Annotation | 40 min | 20 min | 60 min |
| 7 | Multi-Sample Integration | 45 min | 15 min | 60 min |
| 8 | Differential Expression Analysis | 40 min | 20 min | 60 min |
| | **Total** | **325 min** | **120 min** | **445 min** |

The workshop is designed for a two-day schedule with approximately 3.5--4 hours
of instruction per day (excluding breaks).

**Suggested day 1:** Episodes 1--4 (205 min instruction + breaks)
**Suggested day 2:** Episodes 5--8 (240 min instruction + breaks)

Build in 10-minute breaks between episodes and a 45-minute lunch break. Expect
the actual pace to vary: Episode 2 may run long if SLURM queue times are slow,
and Episodes 6 and 8 often generate discussion.

## Delivery on Negishi (one day, 2026-10-29)

**RStudio sessions.** Learners launch RStudio (bioconductor) on Open OnDemand
with account `rcac-rnaseq`, 8 cores (about 16 GB), 8 hours, R 4.6.1 /
Bioconductor 3.23, and "Load extra RCAC site library" checked (setup.md,
Step 3). Seurat, SingleR, celldex, SeuratData, presto, metap and anndataR come
from that library; a learner with missing-package errors has most likely
unchecked it. 20 learners at 8 cores use 160 of the account's 192 cores.

**Restart R at the start of each episode** (Session > Restart R; setup.md,
Step 5). Each episode reloads its object from scratch, so nothing is lost. The
measured peak memory of a fresh session per episode (Negishi, 2026-10-23) is
3.5, 14.1, 2.5, 6.6, 2.7 and 7.3 GiB for episodes 3 to 8 (8 cores give about 14.9 GiB). A session that is
never restarted carries every earlier object along. The 14.1 GiB in episode 4
is the self-paced SCTransform challenge; its main path peaks at 3.0 GiB. That
challenge says to run it in a new session with 12 or more cores.

**Episode 2: check, do not run.** The learner scripts ask for 64 cores; Cell
Ranger runs about 5.5 hours and the account can run only a few such jobs at
once. In class, learners do the login-node steps, write the three job scripts
and check each with `sbatch --test-only`. Before the workshop, submit the same
three scripts under your own account and show the outputs from that run
(`web_summary.html`, `Log.final.out`, the cell and gene counts). At the end of
the day (16:00), learners delete their RStudio session and submit
`run_cellranger.sh` for real as follow-up (the STAR jobs are optional). With
the RStudio sessions closed, 192 cores run three 64-core jobs at once, so 20
Cell Ranger jobs take about a day and a half to clear; tell learners their job
may wait a day or more.

## Common Issues

### R package installation failures

BiocManager version mismatches are the most common source of installation
errors. If learners see messages like `version '3.19' is required but '3.18' is
installed`, have them run:

```r
install.packages("BiocManager")
BiocManager::install(version = "3.20")  # match the R version on the cluster
```

On RCAC clusters, R packages should be installed from a compute node or RStudio
Server session, not the login node. The login node may have an older R version
or restricted network access.

### Memory errors on login nodes

Seurat operations (especially `ScaleData`, `RunPCA`, `FindAllMarkers`, and
`SCTransform`) require substantial memory. Running these on a login node will
either fail with an out-of-memory error or get killed by the resource monitor.

**Solution:** Always use a compute node allocation or RStudio (bioconductor)
through Open OnDemand, with 8 or more cores (about 16 GB; see "Delivery on
Negishi" above).

### STAR index building OOM

Building a STAR genome index for the human genome requires approximately 32 GB
of RAM. If the SLURM job requests less memory, STAR will crash with a
segmentation fault or `std::bad_alloc` error.

**Solution:** Negishi allocates about 2 GB of memory per core, so the core count
sets the memory; do not add `--mem`. The episode's script asks for 64 cores
(about 128 GB), and the index build peaked at 34 GB on Negishi. A learner who
lowered `--ntasks` below about 18 cores will hit this.

### Cell Ranger license not found

Cell Ranger occasionally fails with a licensing error if the `--jobmode=local`
flag conflicts with SLURM environment variables, or if the binary path is
misconfigured in the module file.

**Solution:** Verify `which cellranger` resolves correctly after module load.
If licensing errors persist, STARsolo produces equivalent results and can be
used as a drop-in alternative.

### Seurat v4 vs v5 syntax confusion

Learners who have prior experience with Seurat v4 may try to use slot-based
syntax (`pbmc@assays$RNA@counts`) instead of the v5 layer syntax
(`pbmc[["RNA"]]$counts`). The old syntax will silently return incorrect results
or throw errors.

**Solution:** Emphasize the v5 layer syntax in Episode 3 (the callout box
explicitly shows both). If a learner has Seurat v4 installed locally, they must
upgrade or use the cluster installation.

## Pre-Workshop Checklist

- [ ] **Data staging:** Verify workshop data is accessible at
  `/depot/workshop/data/scrna_workshop/` and contains `fastq/`, `reference/`,
  and `filtered_feature_bc_matrix/` directories
- [ ] **Module availability:** Confirm that `ml biocontainers; ml cellranger/10.0.0`
  and `ml biocontainers; ml star` load without errors on Negishi (episode 2
  pins Cell Ranger 10.0.0; the Negishi default is newer)
- [ ] **R packages:** In an RStudio (bioconductor) session with the workshop
  settings (setup.md, Step 3), confirm these load. They come from the image
  or from the extra RCAC site library; `negishi-run/00_preflight.sh` checks
  them all:
  - `Seurat`, `ggplot2`, `patchwork`, `dplyr`
  - `clustree`, `glmGamPoi`
  - `SingleR`, `celldex`, `presto` (presto makes `FindAllMarkers` fast)
  - `SeuratData` (and run `InstallData("ifnb")`), `metap`, `multtest`
  - `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `anndataR` (Episode 8
    export; it replaced sceasy, which cannot convert Seurat v5 objects)
- [ ] **Pre-computed objects:** Generate and stage all intermediate `.rds` files
  so learners who fall behind can catch up:
  - `pbmc_filtered.rds` (end of Episode 3)
  - `pbmc_normalized.rds` (end of Episode 4)
  - `pbmc_clustered.rds` (end of Episode 5)
  - `pbmc_annotated.rds` (end of Episode 6)
  - `ifnb_annotated.rds` (end of Episode 7, used in Episode 8)
- [ ] **Build test:** Run `sandpaper::build_lesson()` and confirm it renders
  without errors
- [ ] **RStudio (bioconductor):** Verify the Open OnDemand app launches R 4.6.1 /
  Bioconductor 3.23 under `rcac-rnaseq` with 8 cores and the extra RCAC site
  library checked
- [ ] **SLURM account:** Confirm `rcac-rnaseq` is active, has its 192 cores, and
  that learners are members. `sbatch --test-only run_cellranger.sh` should
  print an estimated start, not an account error.
- [ ] **Pre-run for episode 2:** Submit `run_cellranger.sh`, `build_star_index.sh`
  and `run_starsolo.sh` under your own account at least a day ahead, so the
  outputs are ready to show
- [ ] **Network access:** Confirm compute nodes can reach Bioconductor and CRAN
  mirrors (needed for `celldex` reference downloads in Episode 6)

## Challenge Exercise Discussion Points

### Episode 2: Raw Data Processing

**Challenge 1 -- Modify STARsolo for v2 Chemistry:**
The key point is that v2 and v3 chemistry differ in UMI length (10 vs 12 bp)
and whitelist file (737K vs 3M barcodes). The cell barcode length (16 bp) stays
the same. Ask learners: "What would happen if you ran v3 parameters on v2
data?" Answer: very few valid barcodes would match because the wrong whitelist
is being searched, and the extra 2 bases of the "UMI" would actually be cDNA
sequence, producing incorrect UMI deduplication.

**Challenge 2 -- QC from Cell Ranger Summary:**
This exercises interpretation skills. The two concerning metrics are Fraction
Reads in Cells (55.2%, should be >70%) and Sequencing Saturation (28.3%,
suggesting under-sequencing). Discuss how ambient RNA contamination could be
addressed with tools like SoupX or CellBender. Note that low saturation does
not mean the experiment failed -- it means additional sequencing would be
beneficial.

### Episode 3: Quality Control

**Challenge 1 -- Propose Your Own Thresholds:**
There is no single correct answer. Encourage learners to justify their
choices based on the distributions they observed. Common proposals include
nFeature_RNA 300--4000 and percent.mt 10--12. Discuss the trade-off:
aggressive filtering removes more noise but may also remove real cells
(especially small cell types like platelets). Ask: "Which cell types are
most at risk from a low nFeature_RNA threshold?"

**Challenge 2 -- Aggressive Mitochondrial Filtering:**
With percent.mt < 5, only 3,276 of the 11,310 filtered cells remain, a loss of
71%, because the median mitochondrial percentage is about 5.8%. The loss is
uneven: about 45% of CD4 T and NK cells survive, but only 12% of CD8 T cells
and FCGR3A+ monocytes and 24% of CD14 monocytes. Use this to reinforce that
thresholds must be tuned to the specific tissue and dataset.

### Episode 4: Normalization and Feature Selection

**Challenge 1 -- LogNormalize vs. SCTransform Variable Features:**
Expect 7 of 20 shared (PPBP, PF4, GP1BB, IGKC, GNLY, LINGO2, CAVIN2). The
LogNormalize list leans to sparse, low-expressed genes (LINC01478, PTGDS,
SOX5), the SCTransform list to strongly expressed markers (S100A9, S100A8,
LYZ, NKG7). This challenge is self-paced and needs 12 or more cores. Emphasize that
downstream clustering is typically robust to these differences.

**Challenge 2 -- Housekeeping vs. Variable Gene Behavior:**
ACTB's distribution tightens after normalization (variation was mostly
technical), while LYZ retains its bimodal pattern (real biological variation
between monocytes and non-monocytes). This illustrates the purpose of both
normalization (remove technical variation) and variable feature selection
(keep biological variation).

### Episode 5: Dimensionality Reduction and Clustering

**Challenge 1 -- Choosing a Resolution with Clustree:**
Guide learners through reading the clustree plot. Stable splits (clean
branching) suggest real biology; unstable splits (cells scattering across
multiple clusters) suggest over-splitting. Resolutions 0.2, 0.5, 0.8 and 1.2
give 16, 18, 21 and 25 clusters. The workshop uses 0.5 (18 clusters, which
episode 6 merges into 10 cell types); 1.2 over-splits.

**Challenge 2 -- Mapping Markers to Clusters:**
This is a preview of Episode 6. Learners should identify that CD3D marks
multiple clusters (T cell subtypes), MS4A1 is specific to one cluster (B
cells), and LYZ marks one or two clusters (monocytes). Ask: "Why does CD3D
mark multiple clusters?" Answer: because CD3D is a pan-T cell marker expressed
in all T cell subtypes, which are resolved into separate clusters.

### Episode 6: Cell Type Annotation

**Challenge 1 -- Identifying a T Cell Subtype:**
CD3D and IL7R are expressed in almost every CD4 T cell, but SELL (69%) and
CCR7 (55%) only in part of the group: it mixes naive and central memory cells
(SELL+ CCR7+) with effector memory cells (SELL- CCR7-). CD69 is expressed in
about 90% of T, B and NK cells alike, so it does not mark activated cells here
(it is also induced during sample processing). Use this to discuss marker
combinations and checking that a marker is specific before trusting it.

**Challenge 2 -- Comparing Two SingleR References:**
Monaco Immune works better for PBMCs because it was built from sorted immune
populations. The two agree on B cells, monocytes, NK and T cells, but HPCA
cannot split T cells into CD4 and CD8, calls most Monaco dendritic cells
"Monocyte", and labels Monaco's progenitors mostly "Platelets", with a few
bone-marrow classes (CMP, GMP, MEP, HSC). Discuss the importance
of choosing a reference appropriate for your tissue type.

### Episode 7: Multi-Sample Integration

Since patch 24 the cell types come from the `seurat_annotations` column that
ships with ifnb (13 labels), not from a hardcoded vector; the
cluster-versus-label table is the integration check.

**Challenge 1 -- Quantifying Integration Quality:**
The proportion table shows 39--61% per condition within each cell type
(CD14 Mono 50.8/49.2), against a baseline of 47/53 because STIM has more
cells. If a cell type is heavily skewed (e.g., 90% from one condition), it
could be a condition-specific population or a sign of incomplete integration.
Ask: "What would it mean if a cluster was 95% STIM cells?" It could represent
a stimulation-induced cell state that does not exist in control.

**Challenge 2 -- Interferon Response Genes:**
The top 10 are led by IFN-inducible chemokines (CCL8, CXCL11, CXCL10) and
HESX1, followed by ISGs (IFIT1, IFIT2, TNFSF10, IFIT3) and CD38. Use this
to validate the integration: if the biological signal was lost during
integration, these genes would not be differentially expressed. This challenge
connects back to the biology -- IFN-beta activates the JAK-STAT pathway, which
induces ISG transcription.

**Challenge 3 -- Label the Integrated Clusters Yourself:**
The solution's 14 labels give each cluster the shipped label most of its
cells carry; the comparison table puts most cells on the diagonal. Discuss
the off-diagonal cells (activated B cells inside the B cluster, CD4 naive and
memory cells shared between clusters 1 and 2) and cluster 13, which has no
clear marker. Before patch 24 the episode hardcoded a vector like this one,
and 11 of its 14 labels were wrong (for example NK cells labelled CD14 Mono);
the guard and the kit's marker anchors exist because of that.

### Episode 8: Differential Expression Analysis

**Challenge 1 -- Wilcoxon vs. Pseudobulk:**
Wilcoxon finds 825 significant genes and pseudobulk 1,539 (patch 24 labels; Negishi, 2026-10-23). Pseudobulk
finds more here only because the episode builds artificial pseudo-replicates
from one sample per condition; the episode explains this. The key teaching
point is still pseudoreplication: per-cell tests treat each cell as an
independent sample.
Discuss when per-cell tests are acceptable (single-sample, exploratory) vs.
when pseudobulk is required (multi-sample, publication).

**Challenge 2 -- GO Enrichment:**
Top GO terms all relate to the antiviral response: "response to virus",
"defense response to virus", "regulation of viral process". No type I
interferon term reaches the top (patch 18).
Use this as the capstone: the full pipeline (from raw data through integration
and DE) correctly identifies the expected biology. Ask learners to think about
what GO terms they would expect for their own research questions.
