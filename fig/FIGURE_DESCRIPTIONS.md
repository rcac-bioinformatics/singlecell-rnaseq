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
| **Location** | Line 93 |
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
| **Location** | Line 129 |
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
| **Location** | Line 178 |
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
| **Location** | Line 215 |
| **Description** | Flowchart of the complete scRNA-seq analysis pipeline covered in this workshop. Show 8 boxes arranged vertically (or in a slight zigzag), connected by downward arrows. Each box represents one episode: (1) Introduction, (2) Raw Data Processing, (3) Quality Control, (4) Normalization & Feature Selection, (5) Dimensionality Reduction & Clustering, (6) Cell Type Annotation, (7) Multi-Sample Integration, (8) Differential Expression. Each box should include the episode number and a one-line subtitle (e.g., "FASTQ to count matrix", "Filter low-quality cells"). Optionally color-code boxes by tool: grey for conceptual, blue for bash/HPC, green for R/Seurat. Include input/output annotations on the arrows (e.g., "FASTQ files" entering step 2, "filtered_feature_bc_matrix/" between steps 2 and 3, "Seurat object" flowing through steps 3-8). |
| **Alt text** | Flowchart of the scRNA-seq analysis pipeline covered in this workshop. Eight boxes arranged vertically show the progression: Raw Data Processing, Quality Control, Normalization and Feature Selection, Dimensionality Reduction and Clustering, Cell Type Annotation, Multi-Sample Integration, and Differential Expression. Arrows connect each step to the next. |
| **Suggested dimensions** | 500 x 900 px (tall portrait) |
| **Suggested tool** | draw.io or Mermaid (rendered to PNG) |
| **Auto-generated?** | No -- requires manual layout. A Mermaid diagram could be used as a starting point and exported to PNG. |
| **Style notes** | Use rounded rectangles. Color-code by phase if possible. Keep fonts large. Ensure the figure is readable at the width of the lesson page (~700 px rendered). |


## Code-Generated Plots (produced by R chunks)

These plots are generated by the R code in each episode when chunks run with
`eval=TRUE`. If the lesson is built with `eval=FALSE`, representative
screenshots of these plots should be captured from a live R session and placed
in `episodes/fig/` as fallback images.

The figures below are grouped by episode. Each entry lists the chunk name, a
description, and notes on what the figure should show.

---

### Episode 3: Quality Control

#### 3.1 vlnplot (line 240)

| Field | Value |
|-------|-------|
| **Chunk** | `vlnplot` |
| **Description** | Three violin plots showing distributions of nFeature_RNA, nCount_RNA, and percent.mt across all cells before filtering. Most cells cluster in a central distribution with outliers at the tails. |
| **Seurat function** | `VlnPlot()` |
| **Auto-generated?** | Yes -- `VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)` |

#### 3.2 scatter (line 258)

| Field | Value |
|-------|-------|
| **Chunk** | `scatter` |
| **Description** | Scatter plot of nCount_RNA vs nFeature_RNA, colored by percent.mt using viridis scale. Most cells form a dense diagonal cloud. High-mt cells (yellow) sit below the main diagonal. Outliers in the top-right are potential doublets. |
| **Seurat function** | `FeatureScatter()` + ggplot2 overlay |
| **Auto-generated?** | Yes |

#### 3.3 vlnplot-after (line 352)

| Field | Value |
|-------|-------|
| **Chunk** | `vlnplot-after` |
| **Description** | Same three violin plots as 3.1 but after QC filtering. Distributions are tighter with extreme outliers removed. percent.mt is capped below 15%. |
| **Seurat function** | `VlnPlot()` |
| **Auto-generated?** | Yes |

#### 3.4 scatter-after (line 364)

| Field | Value |
|-------|-------|
| **Chunk** | `scatter-after` |
| **Description** | Scatter plot of nCount_RNA vs nFeature_RNA after filtering. The cloud is bounded (no extreme outliers). |
| **Seurat function** | `FeatureScatter()` |
| **Auto-generated?** | Yes |

---

### Episode 4: Normalization and Feature Selection

#### 4.1 norm-comparison (line 137)

| Field | Value |
|-------|-------|
| **Chunk** | `norm-comparison` |
| **Description** | Two histograms side by side for the LYZ gene. Left: raw UMI counts (right-skewed, spike at zero). Right: log-normalized expression (more spread out, reduced skew for expressing cells). |
| **Seurat function** | Custom ggplot2 histograms using `pbmc[["RNA"]]$counts` and `pbmc[["RNA"]]$data` |
| **Auto-generated?** | Yes |

#### 4.2 variable-plot (line 213)

| Field | Value |
|-------|-------|
| **Chunk** | `variable-plot` |
| **Description** | Mean-variance scatter plot for all genes. X-axis: mean expression. Y-axis: standardized variance. 2,000 selected variable features in red, rest in black. Top 10 labeled with gene names (S100A9, LYZ, GNLY, etc.). |
| **Seurat function** | `VariableFeaturePlot()` + `LabelPoints()` |
| **Auto-generated?** | Yes |

---

### Episode 5: Dimensionality Reduction and Clustering

#### 5.1 viz-loadings (line 89)

| Field | Value |
|-------|-------|
| **Chunk** | `viz-loadings` |
| **Description** | Horizontal bar plots showing the top genes with highest positive and negative loadings for PC1 and PC2. PC1 separates myeloid (CST3, LYZ, S100A9) from lymphoid markers. PC2 captures other axes of variation (CD3D, IL32). |
| **Seurat function** | `VizDimLoadings()` |
| **Auto-generated?** | Yes |

#### 5.2 dim-heatmap (line 103)

| Field | Value |
|-------|-------|
| **Chunk** | `dim-heatmap` |
| **Description** | Grid of 9 heatmaps (PC1 through PC9). Each heatmap shows top genes (rows, ordered by loading) and cells (columns, ordered by PC score). Clear expression blocks indicate PCs capturing real biological structure. |
| **Seurat function** | `DimHeatmap()` |
| **Auto-generated?** | Yes |

#### 5.3 elbow (line 126)

| Field | Value |
|-------|-------|
| **Chunk** | `elbow` |
| **Description** | Line plot of standard deviation vs principal component number (1--30). Steep decline for PC1--10, then gradual flattening. Elbow around PC 12--15. |
| **Seurat function** | `ElbowPlot()` |
| **Auto-generated?** | Yes |

#### 5.4 umap-plot (line 177)

| Field | Value |
|-------|-------|
| **Chunk** | `umap-plot` |
| **Description** | UMAP of ~10,000 cells before clustering (all same color). Cells form distinct visual groups corresponding to future clusters. |
| **Seurat function** | `DimPlot()` |
| **Auto-generated?** | Yes |

#### 5.5 umap-clusters (line 318)

| Field | Value |
|-------|-------|
| **Chunk** | `umap-clusters` |
| **Description** | UMAP with cells colored by cluster identity (resolution 0.5). Approximately 9--11 colored groups with cluster numbers labeled at center of each group. No legend (labels on plot). |
| **Seurat function** | `DimPlot(label = TRUE) + NoLegend()` |
| **Auto-generated?** | Yes |

#### 5.6 clustree (line 341)

| Field | Value |
|-------|-------|
| **Chunk** | `clustree` |
| **Description** | Clustree diagram showing how clusters split across resolutions 0.2, 0.5, 0.8, and 1.2. Each row is a resolution. Nodes are clusters, arrows show cell flow. Stable splits appear as clean branches; unstable splits show cells scattering across multiple targets. |
| **Seurat function** | `clustree()` from the clustree package |
| **Auto-generated?** | Yes |

#### 5.7 feature-plot (line 373)

| Field | Value |
|-------|-------|
| **Chunk** | `feature-plot` |
| **Description** | 2x2 grid of UMAP feature plots for CD3D, MS4A1, LYZ, and GNLY. Each gene lights up a different cluster region: CD3D in T cell clusters, MS4A1 in B cells, LYZ in monocytes, GNLY in NK cells. Color scale from grey (no expression) to purple/blue (high expression). |
| **Seurat function** | `FeaturePlot()` |
| **Auto-generated?** | Yes |

#### 5.8 vln-markers (line 389)

| Field | Value |
|-------|-------|
| **Chunk** | `vln-markers` |
| **Description** | 2x2 grid of violin plots for CD3D, MS4A1, LYZ, and GNLY across all clusters. Each marker shows high expression in one or two clusters and near-zero in others. No individual points (pt.size = 0). |
| **Seurat function** | `VlnPlot()` |
| **Auto-generated?** | Yes |

---

### Episode 6: Cell Type Annotation

#### 6.1 heatmap (line 137)

| Field | Value |
|-------|-------|
| **Chunk** | `heatmap` |
| **Description** | Heatmap of top 3 marker genes per cluster. Cells (columns) grouped by cluster, genes (rows) grouped by cluster assignment. Each cluster shows a distinct block of upregulated genes (yellow/white) against a background of low expression (purple). |
| **Seurat function** | `DoHeatmap()` |
| **Auto-generated?** | Yes |

#### 6.2 dotplot (line 155)

| Field | Value |
|-------|-------|
| **Chunk** | `dotplot` |
| **Description** | Dot plot of top 3 marker genes across all clusters. Dot size = percent of cells expressing the gene. Dot color = average expression level. Genes on y-axis, clusters on x-axis (flipped with coord_flip). |
| **Seurat function** | `DotPlot() + coord_flip()` |
| **Auto-generated?** | Yes |

#### 6.3 feature-markers (line 179)

| Field | Value |
|-------|-------|
| **Chunk** | `feature-markers` |
| **Description** | 2x4 grid of UMAP feature plots for 8 canonical PBMC markers: CD3D, IL7R, CD8A, MS4A1, LYZ, FCGR3A, GNLY, FCER1A. Each marker lights up in its specific cluster region. |
| **Seurat function** | `FeaturePlot(ncol = 4)` |
| **Auto-generated?** | Yes |

#### 6.4 vln-markers (line 188)

| Field | Value |
|-------|-------|
| **Chunk** | `vln-markers` |
| **Description** | 2x4 grid of violin plots for 8 markers (CD3D, IL7R, CD8A, MS4A1, LYZ, FCGR3A, GNLY, PPBP) across all clusters. Each marker is high in its specific cluster(s). |
| **Seurat function** | `VlnPlot(ncol = 4)` |
| **Auto-generated?** | Yes |

#### 6.5 annotated-umap (line 241)

| Field | Value |
|-------|-------|
| **Chunk** | `annotated-umap` |
| **Description** | UMAP with cells colored and labeled by cell type name (CD4 T, CD8 T, B, CD14 Mono, FCGR3A Mono, NK, DC, Platelet). Labels placed with repel to avoid overlap. No legend. |
| **Seurat function** | `DimPlot(label = TRUE, repel = TRUE) + NoLegend()` |
| **Auto-generated?** | Yes |

#### 6.6 singler-umap (line 322)

| Field | Value |
|-------|-------|
| **Chunk** | `singler-umap` |
| **Description** | UMAP colored by SingleR automated annotation labels (Monaco Immune reference). Labels like "CD4+ T cells", "Classical monocytes", "B cells", "NK cells" should largely agree with manual annotation. |
| **Seurat function** | `DimPlot(group.by = "singler_labels")` |
| **Auto-generated?** | Yes |

#### 6.7 side-by-side (line 358)

| Field | Value |
|-------|-------|
| **Chunk** | `side-by-side` |
| **Description** | Two UMAP panels side by side: left = manual annotation, right = SingleR annotation. Most clusters receive consistent labels. Minor differences in naming conventions (e.g., "CD14 Mono" vs "Classical monocytes"). |
| **Seurat function** | `DimPlot()` + patchwork |
| **Auto-generated?** | Yes |

---

### Episode 7: Multi-Sample Integration

#### 7.1 plot-before (line 155)

| Field | Value |
|-------|-------|
| **Chunk** | `plot-before` |
| **Description** | UMAP of the IFNB dataset colored by condition (CTRL vs STIM) **before** integration. Control and stimulated cells form largely separate clusters, demonstrating the batch effect problem. |
| **Seurat function** | `DimPlot(group.by = "stim")` |
| **Auto-generated?** | Yes |

#### 7.2 plot-after (line 206)

| Field | Value |
|-------|-------|
| **Chunk** | `plot-after` |
| **Description** | UMAP colored by condition **after** CCA integration. Control and stimulated cells now intermingle within each cluster, showing successful batch correction. |
| **Seurat function** | `DimPlot(group.by = "stim")` |
| **Auto-generated?** | Yes |

#### 7.3 split-view (line 215)

| Field | Value |
|-------|-------|
| **Chunk** | `split-view` |
| **Description** | Two UMAP panels split by condition after integration. Both panels show the same cluster structure, confirming shared cell types across conditions. |
| **Seurat function** | `DimPlot(split.by = "stim")` |
| **Auto-generated?** | Yes |

#### 7.4 annotated-umap (line 316)

| Field | Value |
|-------|-------|
| **Chunk** | `annotated-umap` |
| **Description** | UMAP of the integrated IFNB dataset colored by annotated cell type (CD4 T, CD8 T, B, NK, CD14 Mono, FCGR3A Mono, DC, Mk, Eryth). Labels on plot, no legend. |
| **Seurat function** | `DimPlot(label = TRUE, repel = TRUE) + NoLegend()` |
| **Auto-generated?** | Yes |

#### 7.5 isg-feature (line 396)

| Field | Value |
|-------|-------|
| **Chunk** | `isg-feature` |
| **Description** | FeaturePlot of ISG15 and IFIT1, split by condition. Both genes are nearly absent in CTRL but strongly expressed across multiple cell types (especially monocytes) in STIM. Color scale: grey85 (low) to firebrick (high). |
| **Seurat function** | `FeaturePlot(split.by = "stim")` |
| **Auto-generated?** | Yes |

#### 7.6 isg-violin (line 403)

| Field | Value |
|-------|-------|
| **Chunk** | `isg-violin` |
| **Description** | Violin plot of ISG15 in CD14 monocytes, split by condition. STIM shows dramatically higher expression than CTRL. |
| **Seurat function** | `VlnPlot(split.by = "stim", idents = "CD14 Mono")` |
| **Auto-generated?** | Yes |

---

### Episode 8: Differential Expression Analysis

#### 8.1 volcano (line 239)

| Field | Value |
|-------|-------|
| **Chunk** | `volcano` |
| **Description** | Volcano plot for STIM vs CTRL in CD14+ monocytes. X-axis: log2 fold change. Y-axis: -log10(adjusted p-value). Upregulated genes (upper right) in firebrick, downregulated (upper left) in steelblue, non-significant (center) in grey. Dashed lines at |log2FC| = 0.5 and padj = 0.05. Top ISGs (ISG15, IFIT1, IFI6) appear in the upper right corner. |
| **Seurat function** | Custom ggplot2 |
| **Auto-generated?** | Yes |

#### 8.2 de-viz (line 271)

| Field | Value |
|-------|-------|
| **Chunk** | `de-viz` |
| **Description** | FeaturePlot of ISG15 split by condition. Nearly absent in CTRL, strongly expressed in STIM across multiple cell types (especially monocytes). Color: grey85 to firebrick. |
| **Seurat function** | `FeaturePlot(split.by = "stim")` |
| **Auto-generated?** | Yes |

#### 8.3 de-violin (line 278)

| Field | Value |
|-------|-------|
| **Chunk** | `de-violin` |
| **Description** | Three violin plots (ISG15, IFIT1, MX1) across cell types, split by condition. All three genes are dramatically upregulated in STIM across multiple cell types. |
| **Seurat function** | `VlnPlot(split.by = "stim")` |
| **Auto-generated?** | Yes |

#### 8.4 de-dotplot (line 291)

| Field | Value |
|-------|-------|
| **Chunk** | `de-dotplot` |
| **Description** | Dot plot of top 10 DE genes across all cell types, split by condition. ISGs are consistently upregulated in STIM across all cell types. Dot size = percent expressing, color = expression level (steelblue for CTRL, firebrick for STIM). |
| **Seurat function** | `DotPlot(split.by = "stim") + coord_flip()` |
| **Auto-generated?** | Yes |

#### 8.5 go-dotplot (line 358)

| Field | Value |
|-------|-------|
| **Chunk** | `go-dotplot` |
| **Description** | clusterProfiler dot plot of top 15 enriched GO Biological Process terms. Y-axis: GO term names. X-axis: gene ratio. Dot size = gene count. Dot color = adjusted p-value. Top terms include "defense response to virus", "type I interferon signaling pathway", "response to virus", "innate immune response". |
| **Seurat function** | `dotplot()` from clusterProfiler |
| **Auto-generated?** | Yes |


## Summary

| Category | Count | Action needed |
|----------|:-----:|---------------|
| Static figures (manual) | 4 | Must be created with draw.io, BioRender, or similar and saved as PNG in `episodes/fig/` |
| Code-generated plots | 31 | Produced automatically when R chunks run with `eval=TRUE`. For `eval=FALSE` builds, capture representative screenshots. |
| **Total figures** | **35** | |

### Priority for static figure creation

1. **workshop-pipeline-overview.png** -- highest impact, appears in the
   workshop overview section that every learner reads first
2. **10x-chromium-workflow.png** -- essential for understanding the core
   technology
3. **fastq-read-structure.png** -- directly supports Episode 2 interpretation
4. **count-matrix-schematic.png** -- illustrates the key data structure
