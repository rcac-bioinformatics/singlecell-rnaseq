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

**Solution:** Always use a compute node allocation or RStudio Server through
Open OnDemand. Request at least 16 GB of memory for R sessions.

### STAR index building OOM

Building a STAR genome index for the human genome requires approximately 32 GB
of RAM. If the SLURM job requests less memory, STAR will crash with a
segmentation fault or `std::bad_alloc` error.

**Solution:** Ensure the `build_star_index.sh` script requests `--mem=40G` or
higher. The pre-built index in the workshop data directory avoids this step
entirely.

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
- [ ] **Module availability:** Confirm that `ml biocontainers; ml cellranger`
  and `ml biocontainers; ml star` load without errors on Negishi
- [ ] **R packages:** Install all required R packages on training accounts and
  verify they load:
  - `Seurat`, `ggplot2`, `patchwork`, `dplyr`
  - `clustree`, `SingleR`, `celldex`
  - `SeuratData` (and run `InstallData("ifnb")`)
  - `clusterProfiler`, `org.Hs.eg.db`
  - `sceasy` (optional, for Episode 8 export)
- [ ] **Pre-computed objects:** Generate and stage all intermediate `.rds` files
  so learners who fall behind can catch up:
  - `pbmc_filtered.rds` (end of Episode 3)
  - `pbmc_normalized.rds` (end of Episode 4)
  - `pbmc_clustered.rds` (end of Episode 5)
  - `pbmc_annotated.rds` (end of Episode 6)
  - `ifnb_annotated.rds` (end of Episode 7, used in Episode 8)
- [ ] **Build test:** Run `sandpaper::build_lesson()` and confirm it renders
  without errors
- [ ] **RStudio Server:** Verify Open OnDemand RStudio Server app is available
  and launches R with the correct version and package library
- [ ] **SLURM queue:** Confirm the `workshop` account/queue is active and
  accepting jobs. Submit a test job to check wait times.
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
With percent.mt < 5, roughly 35--40% of cells are removed because the median
mitochondrial percentage is ~5.5%. This disproportionately removes metabolically
active immune cells (monocytes, activated T cells). Use this to reinforce that
thresholds must be tuned to the specific tissue and dataset.

### Episode 4: Normalization and Feature Selection

**Challenge 1 -- LogNormalize vs. SCTransform Variable Features:**
Expect ~60--70% overlap in the top 20 variable features. The shared genes
(S100A9, LYZ, GNLY, NKG7, IGKC) are strong cell-type markers that any
reasonable method will identify. The differences are in moderately variable
genes where the two variance estimation approaches disagree. Emphasize that
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
multiple clusters) suggest over-splitting. Resolution 0.5 typically gives
9--11 clusters matching the ~8 major PBMC types. Resolution 0.2 under-clusters
(merges T cell subtypes), while 1.2 over-splits.

**Challenge 2 -- Mapping Markers to Clusters:**
This is a preview of Episode 6. Learners should identify that CD3D marks
multiple clusters (T cell subtypes), MS4A1 is specific to one cluster (B
cells), and LYZ marks one or two clusters (monocytes). Ask: "Why does CD3D
mark multiple clusters?" Answer: because CD3D is a pan-T cell marker expressed
in all T cell subtypes, which are resolved into separate clusters.

### Episode 6: Cell Type Annotation

**Challenge 1 -- Identifying a T Cell Subtype:**
The pattern SELL-high, CCR7-high, CD69-low identifies naive T cells. Use this
to discuss the concept of marker combinations: no single gene is sufficient for
subtype identification. Ask learners what markers would distinguish effector
memory T cells (SELL-low, CCR7-low, possibly CD69-high).

**Challenge 2 -- Comparing Two SingleR References:**
Monaco Immune works better for PBMCs because it was built from sorted immune
populations. HPCA may assign non-immune labels (e.g., "Tissue_stem_cells")
because its reference includes non-immune cell types. Discuss the importance
of choosing a reference appropriate for your tissue type.

### Episode 7: Multi-Sample Integration

**Challenge 1 -- Quantifying Integration Quality:**
The proportion table should show ~45--55% per condition within each cell type.
If a cell type is heavily skewed (e.g., 90% from one condition), it could be a
condition-specific population or a sign of incomplete integration. Ask: "What
would it mean if a cluster was 95% STIM cells?" It could represent a
stimulation-induced cell state that does not exist in control.

**Challenge 2 -- Interferon Response Genes:**
All top 10 DE genes should be ISGs (IFIT1, ISG15, IFI6, MX1, etc.). Use this
to validate the integration: if the biological signal was lost during
integration, these genes would not be differentially expressed. This challenge
connects back to the biology -- IFN-beta activates the JAK-STAT pathway, which
induces ISG transcription.

### Episode 8: Differential Expression Analysis

**Challenge 1 -- Wilcoxon vs. Pseudobulk:**
Wilcoxon identifies thousands of significant genes; pseudobulk identifies
hundreds. The key teaching point is pseudoreplication: per-cell tests treat
each cell as independent (n = 4,000) when the true sample size is n = 2.
Discuss when per-cell tests are acceptable (single-sample, exploratory) vs.
when pseudobulk is required (multi-sample, publication).

**Challenge 2 -- GO Enrichment:**
Top GO terms should all relate to interferon/antiviral response: "defense
response to virus", "type I interferon signaling pathway", "response to virus".
Use this as the capstone: the full pipeline (from raw data through integration
and DE) correctly identifies the expected biology. Ask learners to think about
what GO terms they would expect for their own research questions.
