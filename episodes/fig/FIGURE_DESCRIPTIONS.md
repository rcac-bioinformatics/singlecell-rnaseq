---
title: Figure Descriptions and Specifications
---

This file describes every figure referenced across the workshop episodes.
Figures are divided into two categories:

- **Static figures** -- standalone image files that must be created manually
  and placed in `episodes/fig/`. These are referenced with markdown image
  syntax: `![caption](fig/filename.png)`.
- **Code-generated plots** -- produced by R/Seurat code chunks at build time
  (when `eval=TRUE`). These do not require separate image files but do need
  the analysis pipeline to run. They are listed here for completeness and to
  guide the creation of representative screenshots if the lesson is built
  with `eval=FALSE`.


## Static Figures (require manual creation)

### 1. 10x-chromium-workflow.png

| Field | Value |
|-------|-------|
| **Filename** | `fig/10x-chromium-workflow.png` |
| **Episode** | 1 -- Introduction to Single-Cell RNA-Seq |
| **Location** | `introduction.Rmd` line 94 |
| **Description** | Schematic of the 10x Chromium GEM generation process. Show a microfluidic chip with three input channels (cells, gel beads, partitioning oil) converging into a single channel that produces GEM droplets. Each GEM should contain one gel bead (shown as a sphere covered in oligonucleotide primers) and ideally one cell. Include labels for: cell suspension, gel beads, oil, GEM, cell barcode (on the bead), and UMI (on individual primers). Show a few representative GEMs: one with a cell + bead (correct), one empty (no cell), and optionally one doublet (two cells). |
| **Alt text** | Diagram showing the 10x Chromium microfluidic chip where cells, gel beads, and oil are combined to form gel bead-in-emulsion (GEM) droplets. Each GEM contains one gel bead and ideally one cell. |
| **Suggested dimensions** | 800 x 400 px (wide landscape) |
| **Suggested tool** | draw.io or BioRender |
| **Auto-generated?** | No -- requires manual illustration |
| **Style notes** | Use clean, flat vector style with labeled arrows. Match the color palette of the workshop theme. Keep text large enough to read at 50% zoom. |

### 2. fastq-read-structure.png

| Field | Value |
|-------|-------|
| **Filename** | `fig/fastq-read-structure.png` |
| **Episode** | 1 -- Introduction to Single-Cell RNA-Seq |
| **Location** | `introduction.Rmd` line 130 |
| **Description** | Diagram of the 10x Chromium library read structure. Show three horizontal bars representing Read 1, Read 2, and the I1 index read. Read 1 (28 bp) is divided into two color-coded segments: 16 bp cell barcode and 12 bp UMI. Read 2 (variable length, ~90-150 bp) is a single segment labeled "cDNA insert". I1 (8 bp) is labeled "sample index". Annotate each segment with its length in base pairs. Optionally include a small double-stranded DNA molecule above showing the full library construct with P5/P7 adapters at the ends, with arrows indicating where each read originates. |
| **Alt text** | Schematic showing three sequencing reads from a 10x Chromium library. Read 1 is 28 bp and contains the 16 bp cell barcode followed by the 12 bp UMI. Read 2 is variable length and contains the cDNA insert that maps to the transcriptome. The I1 index read contains the sample index for demultiplexing. |
| **Suggested dimensions** | 800 x 350 px (wide landscape) |
| **Suggested tool** | draw.io or BioRender |
| **Auto-generated?** | No -- requires manual illustration |
| **Style notes** | Use distinct colors for cell barcode (e.g., blue), UMI (e.g., orange), cDNA (e.g., green), and sample index (e.g., purple). Label bp lengths clearly. |

### 3. count-matrix-schematic.png

| Field | Value |
|-------|-------|
| **Filename** | `fig/count-matrix-schematic.png` |
| **Episode** | 1 -- Introduction to Single-Cell RNA-Seq |
| **Location** | `introduction.Rmd` line 179 |
| **Description** | Illustration of a sparse UMI count matrix. Show a grid/table with genes as rows (left axis, labeled with example gene names like ACTB, CD3D, LYZ, MT-CO1) and cell barcodes as columns (top axis, labeled with truncated barcodes like AACG..., TTCG...). Most cells in the grid contain "0" (greyed out or light-colored). A few scattered cells contain small integers (1, 2, 3, 5, 12) in a contrasting color. Include a brief annotation or callout indicating "~95% zeros". Optionally shade the non-zero entries to create a visual heat effect. |
| **Alt text** | Schematic of a count matrix with genes on rows and cells on columns. Most cells in the matrix contain zero, with occasional non-zero integer counts scattered throughout, illustrating the sparsity typical of single-cell RNA-seq data. |
| **Suggested dimensions** | 600 x 500 px (slightly tall) |
| **Suggested tool** | draw.io, Google Slides, or PowerPoint |
| **Auto-generated?** | No -- requires manual illustration. Could alternatively be generated with R (e.g., pheatmap on a toy matrix) but a hand-drawn schematic is clearer for teaching. |
| **Style notes** | Keep the matrix small (6-8 rows x 6-8 columns) for readability. Use grey for zeros and a warm color (red/orange) for non-zero values. |

### 4. workshop-pipeline-overview.png

| Field | Value |
|-------|-------|
| **Filename** | `fig/workshop-pipeline-overview.png` |
| **Episode** | 1 -- Introduction to Single-Cell RNA-Seq |
| **Location** | `introduction.Rmd` line 216 |
| **Description** | Flowchart of the complete scRNA-seq analysis pipeline covered in this workshop. Show 8 boxes arranged vertically (or in a slight zigzag), connected by downward arrows. Each box represents one episode: (1) Introduction, (2) Raw Data Processing, (3) Quality Control, (4) Normalization & Feature Selection, (5) Dimensionality Reduction & Clustering, (6) Cell Type Annotation, (7) Multi-Sample Integration, (8) Differential Expression. Each box should include the episode number and a one-line subtitle (e.g., "FASTQ to count matrix", "Filter low-quality cells"). Optionally color-code boxes by tool: grey for conceptual, blue for bash/HPC, green for R/Seurat. Include input/output annotations on the arrows (e.g., "FASTQ files" entering step 2, "filtered_feature_bc_matrix/" between steps 2 and 3, "Seurat object" flowing through steps 3-8). |
| **Alt text** | Flowchart of the scRNA-seq analysis pipeline covered in this workshop. Eight boxes arranged vertically show the progression: Raw Data Processing, Quality Control, Normalization and Feature Selection, Dimensionality Reduction and Clustering, Cell Type Annotation, Multi-Sample Integration, and Differential Expression. Arrows connect each step to the next. |
| **Suggested dimensions** | 500 x 900 px (tall portrait) |
| **Suggested tool** | draw.io or Mermaid (rendered to PNG) |
| **Auto-generated?** | No -- requires manual layout. A Mermaid diagram could be used as a starting point and exported to PNG. |
| **Style notes** | Use rounded rectangles. Color-code by phase if possible. Keep fonts large. Ensure the figure is readable at the width of the lesson page (~700 px rendered). |

### 5. cellranger-web-summary.png

| Field | Value |
|-------|-------|
| **Filename** | `fig/cellranger-web-summary.png` |
| **Episode** | 2 -- Raw Data Processing |
| **Location** | `raw-data-processing.Rmd` line 479 |
| **Description** | Screenshot of the Summary tab of `web_summary.html` from the pre-run Cell Ranger job (episode 2 runs Cell Ranger 10.0.0 on PBMC 10k). Retake it if the Cell Ranger version or the metrics in the episode change. |
| **Alt text** | Screenshot of the Cell Ranger web_summary.html report for the PBMC 10k dataset showing key QC metrics |
| **Auto-generated?** | No |

### 6. seurat_object.png

| Field | Value |
|-------|-------|
| **Filename** | `fig/seurat_object.png` |
| **Episode** | 3 -- Quality Control |
| **Location** | `quality-control.Rmd` line 169 |
| **Description** | Diagram of the Seurat v5 object: assays with their layers (counts, data, scale.data), cell metadata, reductions and graphs. Redraw it if the episode's description of the object changes. |
| **Alt text** | Diagram of the Seurat v5 object showing the RNA assay with its three layers (counts, data, scale.data) and the per-cell metadata columns |
| **Auto-generated?** | No |

### 7. ood_rstudio_dropdown.png (in `learners/fig/`)

| Field | Value |
|-------|-------|
| **Filename** | `learners/fig/ood_rstudio_dropdown.png` |
| **Page** | Setup |
| **Location** | `learners/setup.md` line 172 |
| **Description** | Screenshot of the Negishi Open OnDemand **Interactive Apps** menu, with RStudio (bioconductor) under Bioinformatics Apps marked as the app to select and RStudio Server under GUIs marked as the wrong one. Crop away the dashboard panels to the right of the menu (account and user names). |
| **Alt text** | Open OnDemand Interactive Apps menu with RStudio (bioconductor) under Bioinformatics Apps marked as the app to select, and RStudio Server under GUIs marked as the wrong one |
| **Auto-generated?** | No -- screenshot |

### 8. ood_rstudio_resources.png (in `learners/fig/`)

| Field | Value |
|-------|-------|
| **Filename** | `learners/fig/ood_rstudio_resources.png` |
| **Page** | Setup |
| **Location** | `learners/setup.md` line 196 |
| **Description** | Screenshot of the RStudio (bioconductor) form filled in with the values in setup.md Step 3 (partition `cpu`, account `rcac-rnaseq`, QoS `normal`, 8 hours, 8 cores, R 4.6.1 / Bioconductor 3.23, Load extra RCAC site library checked). Retake it for each delivery. |
| **Alt text** | Open OnDemand RStudio (bioconductor) form filled in with the workshop settings: partition cpu, account rcac-rnaseq, QoS normal, wall time 8 hours, 8 cores, R 4.6.1 / Bioconductor 3.23, and Load extra RCAC site library checked |
| **Auto-generated?** | No -- screenshot |


## Code-Generated Plots (produced by R chunks)

These figures are pre-rendered PNGs in `episodes/fig/`: the lesson builds with `eval=FALSE`, so
regenerate a figure by hand whenever its chunk or the data it plots changes. The list follows
`negishi-run/generated/figmap.tsv` (chunk to figure, built from the episodes). Each description is
the figure's alt text in the episode, which is the authoritative description; update both together.
Figures listed in `negishi-run/figures_redraw.tsv` are redrawn at their committed size with
`bash negishi-run/redraw_figures.sh --run RUN_ID`; the kit's own capture of every chunk is
1400 x 1000 px under `figures/raw/` in a kit bundle.

---

### Episode 3: Quality Control

#### 3.1 vlnplot (`ep3_vlnplot.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `vlnplot` (`quality-control.Rmd` line 286) |
| **Shown at** | line 295 |
| **Caption** | QC violin plots before filtering |
| **Description** | Violin plots of nFeature_RNA, nCount_RNA, and percent.mt showing the distribution of QC metrics across all cells before filtering |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 3.2 scatter (`ep3_scatter.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `scatter` (`quality-control.Rmd` line 309) |
| **Shown at** | line 317 |
| **Caption** | QC scatter plot colored by mitochondrial percentage |
| **Description** | Scatter plot of total UMI counts versus number of genes per cell colored by mitochondrial percentage, with high-mito cells highlighted in red |
| **Plot function** | `ggplot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 3.3 vlnplot-thresholds (`ep3_vlnplot-thresholds.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `vlnplot-thresholds` (`quality-control.Rmd` line 376) |
| **Shown at** | line 389 |
| **Caption** | Violin plots with filtering thresholds |
| **Description** | Violin plots of nFeature_RNA, nCount_RNA, and percent.mt with red dashed lines showing the filtering thresholds and individual cells visible as points |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 3.4 vlnplot-after (`ep3_vlnplot-after.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `vlnplot-after` (`quality-control.Rmd` line 430) |
| **Shown at** | line 437 |
| **Caption** | QC violin plots after filtering |
| **Description** | Violin plots of nFeature_RNA, nCount_RNA, and percent.mt after filtering, showing tighter distributions with extreme outliers removed |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 3.5 scatter-after (`ep3_scatter-after.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `scatter-after` (`quality-control.Rmd` line 445) |
| **Shown at** | line 453 |
| **Caption** | Scatter plot after filtering |
| **Description** | Scatter plot of total UMI counts versus number of genes after filtering, colored by mitochondrial percentage, showing a tighter cloud with outliers removed |
| **Plot function** | `ggplot()` |
| **Redraw** | By hand from a live session or a kit run |

---

### Episode 4: Normalization and Feature Selection

#### 4.1 norm-comparison (`ep4_norm-comparison.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `norm-comparison` (`normalization-and-features.Rmd` line 238) |
| **Shown at** | line 255 |
| **Caption** | LYZ raw vs. normalized expression |
| **Description** | Two histograms for LYZ. Raw counts show a zero-dominated right-skewed distribution. Log-normalized values reveal a bimodal pattern with peaks near 0.5 and 4.5, separating non-monocytes from LYZ-expressing monocytes |
| **Plot function** | `ggplot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 4.2 variable-plot (`ep4_variable-plot.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `variable-plot` (`normalization-and-features.Rmd` line 317) |
| **Shown at** | line 323 |
| **Caption** | Variable feature selection plot |
| **Description** | Variable feature plot showing mean expression versus standardized variance for all genes, with the top 2000 variable features highlighted in red and the top 10 labeled |
| **Plot function** | `LabelPoints()`, `VariableFeaturePlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 4.3 challenge2 (`ep4_challenge2.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge2` (`normalization-and-features.Rmd` line 578) |
| **Shown at** | line 604 |
| **Caption** | ACTB and LYZ raw vs. normalized expression |
| **Description** | Four histograms in a 2x2 grid. Top row: ACTB raw counts show a right-skewed distribution peaking near zero; ACTB normalized shows a unimodal bell shape centered around 3.5. Bottom row: LYZ raw counts show a massive zero spike with a flat tail to 500; LYZ normalized reveals a bimodal pattern with peaks near 0.5 and 4.5 |
| **Plot function** | `ggplot()` |
| **Redraw** | By hand from a live session or a kit run |

---

### Episode 5: Dimensionality Reduction and Clustering

#### 5.1 viz-loadings (`ep5_viz-loadings.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `viz-loadings` (`dimensionality-reduction-clustering.Rmd` line 122) |
| **Shown at** | line 125 |
| **Caption** | Dot plots showing the top 30 genes by loading magnitude for PC1 and PC2. PC1 shows myeloid markers such as LYZ, CST3, and S100A9. PC2 shows B cell markers such as BANK1, CD79A, and MS4A1. |
| **Description** | Dot plots of the top 30 genes by loading magnitude for PC1 and PC2, with PC1 dominated by myeloid markers and PC2 by B cell markers. |
| **Plot function** | `VizDimLoadings()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.2 dim-heatmap (`ep5_dim-heatmap.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `dim-heatmap` (`dimensionality-reduction-clustering.Rmd` line 137) |
| **Shown at** | line 152 |
| **Caption** | Grid of nine heatmaps showing gene expression patterns for PC1 through PC9. Each heatmap shows the top genes ordered by loading with cells ordered by PC score, revealing distinct expression programs. |
| **Description** | Grid of nine heatmaps for PC1 through PC9, each showing top genes ordered by loading and cells ordered by PC score to reveal distinct expression programs. |
| **Plot function** | `DimHeatmap()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.3 elbow (`ep5_elbow.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `elbow` (`dimensionality-reduction-clustering.Rmd` line 163) |
| **Shown at** | line 167 |
| **Caption** | Dot plot of standard deviation versus principal component number. The curve drops steeply for the first 7 PCs then gradually flattens, with an elbow around PC 8 to 10. |
| **Description** | Elbow plot of standard deviation versus PC number showing a steep decline through PC7 and a plateau after PC10, with the elbow around PC 8 to 10. |
| **Plot function** | `ElbowPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.4 umap-plot (`ep5_umap-plot.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `umap-plot` (`dimensionality-reduction-clustering.Rmd` line 217) |
| **Shown at** | line 221 |
| **Caption** | UMAP plot of all cells shown in a single color before clustering. Cells naturally organize into several spatially distinct groups corresponding to different cell types. |
| **Description** | UMAP scatter plot of 11310 cells in a single color before clustering, with cells forming several spatially distinct groups that correspond to different cell types. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.5 umap-clusters (`ep5_umap-clusters.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `umap-clusters` (`dimensionality-reduction-clustering.Rmd` line 374) |
| **Shown at** | line 379 |
| **Caption** | UMAP plot with 18 clusters labeled 0 through 17, each shown in a distinct color. Large clusters like 0 and 1 dominate the center and bottom, while smaller clusters appear at the periphery. |
| **Description** | UMAP plot colored by Louvain cluster identity at resolution 0.5, showing 18 clusters numbered 0 through 17 with large clusters in the center and small clusters at the periphery. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.6 clustree (`ep5_clustree.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `clustree` (`dimensionality-reduction-clustering.Rmd` line 399) |
| **Shown at** | line 402 |
| **Caption** | Clustree diagram showing how clusters split as resolution increases from 0.2 to 1.2. At low resolution a few large clusters exist, which progressively split into more subclusters at higher resolutions. Stable splits show clean branches while unstable splits show cells moving between multiple clusters. |
| **Description** | Clustree diagram tracing cluster membership from resolution 0.2 to 1.2, with clean branches indicating stable splits and tangled arrows indicating unstable over-splitting. |
| **Plot function** | `clustree()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.7 feature-plot (`ep5_feature-plot.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `feature-plot` (`dimensionality-reduction-clustering.Rmd` line 432) |
| **Shown at** | line 438 |
| **Caption** | Four UMAP plots showing expression of CD3D, MS4A1, LYZ, and GNLY. Each gene is highlighted in a different cluster region, indicating T cells, B cells, monocytes, and NK cells respectively. |
| **Description** | Four UMAP panels showing expression of CD3D, MS4A1, LYZ, and GNLY on a grey-to-purple gradient, each marking distinct cell populations for T cells, B cells, monocytes, and NK cells. |
| **Plot function** | `FeaturePlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.8 vln-markers (`ep5_vln-markers.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `vln-markers` (`dimensionality-reduction-clustering.Rmd` line 451) |
| **Shown at** | line 457 |
| **Caption** | Violin plots showing expression of CD3D, MS4A1, LYZ, and GNLY across all clusters. Each marker shows high expression in one or two specific clusters and low or zero expression in the others. |
| **Description** | Violin plots of CD3D, MS4A1, LYZ, and GNLY expression across 18 clusters, showing each marker with high expression in specific clusters and near-zero expression elsewhere. |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.9 challenge1 (`ep5_challenge1.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge1` (`dimensionality-reduction-clustering.Rmd` line 480) |
| **Shown at** | line 503 |
| **Caption** | Four UMAP plots comparing clustering at resolutions 0.2, 0.5, 0.8, and 1.2. Lower resolutions produce fewer large clusters while higher resolutions split cells into progressively more subclusters. |
| **Description** | Side-by-side UMAP plots at resolutions 0.2, 0.5, 0.8, and 1.2 showing progressively more clusters as resolution increases. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.10 challenge2 (`challenge2_plot1.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge2` (`dimensionality-reduction-clustering.Rmd` line 553) |
| **Shown at** | line 561 |
| **Caption** | FeaturePlot of CD3D, MS4A1, and LYZ on the UMAP, each highlighting distinct cell populations |
| **Description** | Three UMAP panels showing expression of CD3D, MS4A1, and LYZ. CD3D expression is concentrated in the upper-left and lower-left regions. MS4A1 lights up a compact group in the upper-right. LYZ is strongest in the lower-left clusters. |
| **Plot function** | `FeaturePlot()`, `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 5.11 challenge2 (`ep5_challenge2.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge2` (`dimensionality-reduction-clustering.Rmd` line 553) |
| **Shown at** | line 563 |
| **Caption** | Cluster-labeled UMAP for reference |
| **Description** | UMAP plot with 18 clusters labeled 0 through 17 for cross-referencing with the FeaturePlot above. |
| **Plot function** | `FeaturePlot()`, `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

---

### Episode 6: Cell Type Annotation

#### 6.1 load-data (`ep6_load-data.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `load-data` (`cell-type-annotation.Rmd` line 76) |
| **Shown at** | line 87 |
| **Caption** | Clustered UMAP from previous episode |
| **Description** | UMAP plot with 18 clusters labeled 0 through 17 at resolution 0.5, matching the clustering results from the previous episode. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 6.2 heatmap (`ep6_heatmap.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `heatmap` (`cell-type-annotation.Rmd` line 227) |
| **Shown at** | line 234 |
| **Caption** | Heatmap of top 3 markers per cluster |
| **Description** | Heatmap showing expression of the top 3 marker genes per cluster across all 18 clusters. Each cluster shows a distinct block of upregulated genes, confirming distinct transcriptional identities. |
| **Plot function** | `DoHeatmap()` |
| **Redraw** | By hand from a live session or a kit run |

#### 6.3 dotplot (`ep6_dotplot.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `dotplot` (`cell-type-annotation.Rmd` line 247) |
| **Shown at** | line 252 |
| **Caption** | Dot plot of top marker genes |
| **Description** | Dot plot showing expression of top 3 marker genes across all 18 clusters. Dot size indicates the percentage of cells expressing each gene; dot shading indicates average expression level. |
| **Plot function** | `DotPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 6.4 feature-markers (`ep6_feature-markers.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `feature-markers` (`cell-type-annotation.Rmd` line 275) |
| **Shown at** | line 282 |
| **Caption** | FeaturePlot of 8 canonical PBMC markers |
| **Description** | Grid of 8 UMAP panels showing expression of CD3D, IL7R, CD8A, MS4A1, LYZ, FCGR3A, GNLY, and FCER1A. Each marker highlights a distinct UMAP region: CD3D and IL7R in the upper-left T cell area, CD8A in a small cluster on the right, MS4A1 in B cells on the right, LYZ in monocytes at the lower-left, FCGR3A in a small cluster at the lower-right, GNLY in NK cells at the upper-center, and FCER1A in a small dendritic cell group. |
| **Plot function** | `FeaturePlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 6.5 vln-markers (`ep6_vln-markers.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `vln-markers` (`cell-type-annotation.Rmd` line 286) |
| **Shown at** | line 294 |
| **Caption** | Violin plots of 8 PBMC markers across clusters |
| **Description** | Violin plots of CD3D, IL7R, CD8A, MS4A1, LYZ, FCGR3A, GNLY, and PPBP across all 18 clusters. Each marker shows high expression in the clusters corresponding to its known cell type and low expression elsewhere. |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 6.6 annotated-umap (`ep6_annotated-umap.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `annotated-umap` (`cell-type-annotation.Rmd` line 374) |
| **Shown at** | line 378 |
| **Caption** | Annotated UMAP with cell type labels |
| **Description** | UMAP plot with cells labeled by manually assigned cell type names: CD4 T, CD14 Mono, B, CD8 T, NK, FCGR3A+ Mono, DC, Platelet, pDC, and HSPC. |
| **Plot function** | `DimPlot()` |
| **Redraw** | Yes -- `negishi-run/redraw_figures.sh` |

#### 6.7 singler-umap (`ep6_singler-umap.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `singler-umap` (`cell-type-annotation.Rmd` line 457) |
| **Shown at** | line 462 |
| **Caption** | SingleR annotated UMAP |
| **Description** | UMAP plot with cells colored by SingleR automated annotation labels from the Monaco Immune reference, showing CD4+ T cells, CD8+ T cells, T cells, Monocytes, B cells, NK cells, Dendritic cells, Progenitors, and Basophils. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 6.8 side-by-side (`ep6_side-by-side.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `side-by-side` (`cell-type-annotation.Rmd` line 543) |
| **Shown at** | line 552 |
| **Caption** | Manual vs SingleR annotation side by side |
| **Description** | Two UMAP plots side by side. Left panel shows manual annotation with CD4 T, CD14 Mono, B, CD8 T, NK, FCGR3A+ Mono, DC, Platelet, pDC, and HSPC labels. Right panel shows SingleR annotation with CD4+ T cells, CD8+ T cells, T cells, Monocytes, B cells, NK cells, Dendritic cells, Progenitors, and Basophils. Most cell groups receive consistent labels between the two methods. |
| **Plot function** | `DimPlot()` |
| **Redraw** | Yes -- `negishi-run/redraw_figures.sh` |

#### 6.9 challenge1 (`ep6_challenge1.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge1` (`cell-type-annotation.Rmd` line 619) |
| **Shown at** | line 627 |
| **Caption** | Violin plots of T cell subtype markers |
| **Description** | Violin plots of CD3D, SELL, CCR7, CD69, and IL7R across the ten annotated cell types. CD3D and IL7R are high in CD4 T and CD8 T cells and near zero elsewhere. SELL is expressed in most groups. CCR7 is expressed mainly in CD4 T cells. CD69 is high in CD4 T, CD8 T, B, NK and HSPC cells alike. |
| **Plot function** | `VlnPlot()` |
| **Redraw** | Yes -- `negishi-run/redraw_figures.sh` |

---

### Episode 7: Multi-Sample Integration

#### 7.1 plot-before (`ep7_plot-before.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `plot-before` (`multi-sample-integration.Rmd` line 176) |
| **Shown at** | line 181 |
| **Caption** | UMAP plot of IFNB dataset before integration. CTRL (salmon) and STIM (teal) cells form largely separate clusters. Several clusters are dominated by one condition, indicating batch-driven rather than cell-type-driven clustering. |
| **Description** | UMAP plot of the IFNB dataset before integration showing CTRL and STIM cells forming separate clusters driven by batch effects rather than cell type identity. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.2 plot-after (`ep7_plot-after.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `plot-after` (`multi-sample-integration.Rmd` line 230) |
| **Shown at** | line 233 |
| **Caption** | UMAP plot after CCA integration. CTRL and STIM cells are now intermingled within each cluster, with teal (STIM) dominating the overlay. Cell types co-cluster regardless of condition. |
| **Description** | UMAP plot of the IFNB dataset after CCA integration showing control and stimulated cells intermingled within each cluster, indicating successful batch correction. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.3 split-view (`ep7_split-view.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `split-view` (`multi-sample-integration.Rmd` line 240) |
| **Shown at** | line 243 |
| **Caption** | Two side-by-side UMAP panels split by condition (CTRL and STIM) after integration. Both panels show 14 clusters numbered 0 through 13, with the same spatial arrangement in each panel, confirming successful integration. |
| **Description** | Split UMAP view showing CTRL and STIM conditions side by side after integration, with both panels displaying the same 14 clusters in identical spatial arrangement. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.4 annotated-umap (`ep7_annotated-umap.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `annotated-umap` (`multi-sample-integration.Rmd` line 342) |
| **Shown at** | line 345 |
| **Caption** | UMAP plot with cells colored and labeled by annotated cell type. Nine cell types are visible: CD14 Mono (large cluster, right), CD4 T (large, upper left), CD8 T (upper right), B (center left), NK (left), FCGR3A Mono (lower right), DC (small, center), Mk (small, center right), and Eryth (small, bottom). |
| **Description** | Annotated UMAP plot of the integrated IFNB dataset with nine labeled cell types including CD14 Mono, CD4 T, CD8 T, B, NK, FCGR3A Mono, DC, Mk, and Eryth. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.5 isg-feature (`ep7_isg-feature.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `isg-feature` (`multi-sample-integration.Rmd` line 480) |
| **Shown at** | line 487 |
| **Caption** | Four-row FeaturePlot split by CTRL and STIM for ISG15, IFIT1, TNFSF10, and RSAD2. CTRL panels show mostly grey cells with sparse expression. STIM panels show strong red expression across monocyte and other clusters for all four genes. |
| **Description** | FeaturePlot of ISG genes split by condition showing near-absent expression in CTRL panels and strong red expression in STIM panels, especially in monocyte clusters. |
| **Plot function** | `FeaturePlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.6 isg-violin (`ep7_isg-violin.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `isg-violin` (`multi-sample-integration.Rmd` line 490) |
| **Shown at** | line 499 |
| **Caption** | Four violin plots for IFIT1, IFIT3, TNFSF10, and RSAD2 in CD14 monocytes. Each panel shows two violins: a thin line at zero for control and a broad teal shape peaking around expression level 3 to 4 for stimulated, confirming strong IFN-beta induction. |
| **Description** | Violin plots of IFIT1, IFIT3, TNFSF10, and RSAD2 in CD14 monocytes showing flat expression in control and broad high expression in stimulated cells. |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.7 challenge1 (`ep7_challenge1.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge1` (`multi-sample-integration.Rmd` line 516) |
| **Shown at** | line 542 |
| **Caption** | Split UMAP with CTRL and STIM panels, cells labeled by cell type. Both panels show identical cluster arrangement with CD4 T, CD8 T, B, NK, CD14 Mono, FCGR3A Mono, DC, Mk, and Eryth labels, confirming balanced condition mixing within each cell type. |
| **Description** | Split UMAP showing CTRL and STIM panels with cell type labels, demonstrating balanced mixing of conditions within each annotated cluster. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 7.8 challenge2 (`ep7_challenge2.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge2` (`multi-sample-integration.Rmd` line 582) |
| **Shown at** | line 599 |
| **Caption** | Four-row FeaturePlot split by CTRL and STIM showing CCL8, CXCL11, HESX1, and CXCL10. CTRL panels are almost entirely grey. STIM panels show strong red expression concentrated in the monocyte clusters for all four chemokine and ISG genes. |
| **Description** | FeaturePlot of CCL8, CXCL11, HESX1, and CXCL10 split by condition showing absent expression in control and strong monocyte-concentrated expression in stimulated cells. |
| **Plot function** | `FeaturePlot()` |
| **Redraw** | Yes -- `negishi-run/redraw_figures.sh` |

---

### Episode 8: Differential Expression Analysis

#### 8.1 load-data (`ep8_umap_split.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `load-data` (`differential-expression.Rmd` line 78) |
| **Shown at** | line 91 |
| **Caption** | Split UMAP showing CTRL and STIM panels with cells colored and labeled by cell type. |
| **Description** | UMAP split by CTRL and STIM conditions with nine cell types colored and labeled, showing identical spatial arrangement in both panels confirming successful integration. |
| **Plot function** | `DimPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.2 volcano (`ep8_volcano.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `volcano` (`differential-expression.Rmd` line 340) |
| **Shown at** | line 361 |
| **Caption** | Volcano plot for CD14 monocytes STIM vs CTRL. |
| **Description** | Volcano plot showing upregulated genes in red in the upper right, downregulated genes in blue in the upper left, and non-significant genes in grey near the bottom, with dashed cutoff lines for significance and fold change. |
| **Plot function** | `ggplot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.3 de-viz (`ep8_featureplot_ifit1.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `de-viz` (`differential-expression.Rmd` line 374) |
| **Shown at** | line 381 |
| **Caption** | FeaturePlot of IFIT1 split by CTRL and STIM. |
| **Description** | Feature plot of IFIT1 expression split by condition showing near-absent expression in control cells and strong expression across all cell types in stimulated cells, especially monocytes. |
| **Plot function** | `FeaturePlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.4 de-violin (`ep8_violin_de.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `de-violin` (`differential-expression.Rmd` line 383) |
| **Shown at** | line 391 |
| **Caption** | Violin plots for IFIT1, CXCL10, and MX1 across cell types, split by condition. |
| **Description** | Three violin plot panels showing IFIT1, CXCL10, and MX1 dramatically upregulated in stimulated cells across all cell types, with control violins near zero and stimulated violins broad and tall. |
| **Plot function** | `VlnPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.5 de-dotplot (`ep8_dotplot_de.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `de-dotplot` (`differential-expression.Rmd` line 398) |
| **Shown at** | line 408 |
| **Caption** | Dot plot of top 10 DE genes across cell types split by condition. |
| **Description** | Dot plot showing top 10 DE genes across all cell types split by CTRL and STIM, with dot size indicating percent expressed and color indicating average expression, showing consistently larger and redder dots in stimulated cells. |
| **Plot function** | `DotPlot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.6 go-dotplot (`ep8_go_dotplot.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `go-dotplot` (`differential-expression.Rmd` line 485) |
| **Shown at** | line 497 |
| **Caption** | GO Biological Process dot plot for STIM vs CTRL CD14 monocytes. |
| **Description** | GO enrichment dot plot showing 15 biological process terms with response to virus, defense response to virus, and viral process as the top three by gene ratio, with dot color from red to blue indicating adjusted p-value. |
| **Plot function** | `dotplot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.7 challenge1 (`ep8_pval_comparison.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge1` (`differential-expression.Rmd` line 689) |
| **Shown at** | line 722 |
| **Caption** | Scatter plot comparing p-values from Wilcoxon and pseudobulk DE tests. |
| **Description** | Scatter plot of negative log10 p-values from Wilcoxon versus pseudobulk for shared genes, with a red dashed diagonal showing equal significance and points broadly scattered indicating different significance rankings between methods. |
| **Plot function** | `ggplot()` |
| **Redraw** | By hand from a live session or a kit run |

#### 8.8 challenge2 (`ep8_go_top200.png`)

| Field | Value |
|-------|-------|
| **Chunk** | `challenge2` (`differential-expression.Rmd` line 767) |
| **Shown at** | line 800 |
| **Caption** | GO BP dot plot for the top 200 upregulated genes in STIM CD14 monocytes. |
| **Description** | GO enrichment dot plot for top 200 upregulated genes showing 10 biological process terms with response to virus and defense response to virus having the largest gene ratios and deepest red color indicating lowest adjusted p-values. |
| **Plot function** | `dotplot()` |
| **Redraw** | By hand from a live session or a kit run |


## Summary

| Category | Count | Action needed |
|----------|:-----:|---------------|
| Static figures and screenshots | 8 | Created by hand; screenshots retaken when the software or the OOD form changes |
| Code-generated plots | 44 | Pre-rendered PNGs; regenerate when the chunk or its data changes |
| **Total figures** | **52** | |
