---
source: Rmd
title: 'Multi-Sample Integration'
teaching: 45
exercises: 15
---

:::::::::::::::::::::::::::::::::::::: questions

- Why do cells from different samples cluster separately even when they are the same cell type?
- How does CCA-based integration correct for batch effects?
- How do we evaluate whether integration was successful without removing real biology?
- How do we find genes that differ between conditions within a specific cell type?
- When should we NOT integrate datasets?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Explain the batch effect problem and why simple merging is insufficient
- Integrate two conditions using Seurat v5 CCA integration with IntegrateLayers
- Evaluate integration by checking condition mixing and preserved marker expression
- Identify conserved markers and condition-specific differentially expressed genes
- Recognize when integration is inappropriate for the experimental design

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::: prereq

## Prerequisites

This episode requires an RStudio session on the Negishi cluster. Launch
**RStudio (bioconductor)** under **Bioinformatics Apps** on
[Open OnDemand](https://gateway.negishi.rcac.purdue.edu) as described in the
[Setup](../learners/setup.md#starting-rstudio-on-open-ondemand) instructions.

:::::::::::::::::::::::::::::::::::::::

## Setup


``` r
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
```

Set up the working directory to match the previous episodes:


``` r
work_dir <- paste0(
    "/scratch/negishi/", Sys.getenv("USER"),
    "/scrna_workshop/"
)
setwd(work_dir)
```

## Why Integration Is Needed

In the previous episodes, we analyzed a single sample of PBMCs. But many
experimental designs involve **multiple samples** -- different patients,
different time points, or different conditions (treated vs. control). When we
combine these samples for joint analysis, a common problem emerges: cells
cluster by **which sample they came from** rather than by **what cell type they
are**.

This happens because of **batch effects**: technical differences between samples
caused by variation in cell handling, library preparation, sequencing depth,
or reagent lots. These technical differences can be large enough to dominate
over the biological differences between cell types, making it impossible to
compare the same cell type across conditions.

**Integration** algorithms correct for batch effects by finding shared
biological structure across samples. The goal is to align the same cell types
so they co-cluster, while preserving the real biological differences (like
changes in gene expression due to treatment).

Let's see this problem in practice using a dataset of PBMCs from control and
IFN-beta stimulated conditions.


## Loading the IFNB Dataset

The IFNB dataset from the `SeuratData` package contains PBMCs from two
conditions:

- **CTRL**: untreated control PBMCs
- **STIM**: PBMCs stimulated with interferon-beta (IFN-beta), a cytokine that
  activates antiviral immune responses


``` r
InstallData("ifnb")
ifnb <- LoadData("ifnb")
ifnb
```

:::::::::::::::::::::::::::::::::::::::::: spoiler

## Expected output

```
An object of class Seurat
14053 features across 13999 samples within 1 assay
Active assay: RNA (14053 features, 0 variable features)
 2 layers present: counts, data
```

The dataset has approximately 14,000 cells and 14,000 genes. The `stim` column
in the metadata records which condition each cell came from.

::::::::::::::::::::::::::::::::::::::::::::::::::


``` r
table(ifnb$stim)
```

:::::::::::::::::::::::::::::::::::::::::: spoiler

## Expected output

```
CTRL STIM
6548 7451
```

Roughly 6,500 control cells and 7,500 stimulated cells.

::::::::::::::::::::::::::::::::::::::::::::::::::


## Preparing Data

In Seurat v5, multi-sample data is handled by **splitting the assay into
layers** -- one layer per sample. This keeps everything in a single object
while allowing integration methods to process each sample separately.


``` r
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
```

:::::::::::::::::::::::::::::::::::::::::: spoiler

## Expected output

```
An object of class Seurat
14053 features across 13999 samples within 1 assay
Active assay: RNA (14053 features, 0 variable features)
 4 layers present: counts.CTRL, counts.STIM, data.CTRL, data.STIM
```

The counts are now stored in two separate layers: `counts.CTRL` and
`counts.STIM`.

::::::::::::::::::::::::::::::::::::::::::::::::::

Now run the standard preprocessing steps. When the assay is split,
`NormalizeData` and `FindVariableFeatures` process each layer independently:


``` r
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb, selection.method = "vst", nfeatures = 2000)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
```

Let's see what happens if we run UMAP and clustering **without** integration:


``` r
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 0.5)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca")
```


``` r
DimPlot(ifnb, reduction = "umap", group.by = "stim")
```


![UMAP plot of IFNB dataset before integration. CTRL (salmon) and STIM (teal) cells form largely separate clusters. Several clusters are dominated by one condition, indicating batch-driven rather than cell-type-driven clustering.](fig/ep7_plot-before.png){alt="UMAP plot of the IFNB dataset before integration showing CTRL and STIM cells forming separate clusters driven by batch effects rather than cell type identity."}

You should see a clear separation between control (CTRL) and stimulated (STIM)
cells. Instead of clustering by cell type (T cells with T cells, monocytes
with monocytes), the cells cluster by condition. This means we cannot directly
compare the same cell type between conditions -- which is exactly what we
want to do. This is the batch effect problem.


## Running Integration

Seurat v5 provides `IntegrateLayers()` as a unified interface for dataset
integration. We will use the **CCA (Canonical Correlation Analysis)** method,
which identifies shared sources of variation between the two conditions and
uses them to align the data.

CCA works by:

1. Finding **canonical correlation vectors** -- directions in gene expression
   space along which the two datasets are maximally correlated
2. Identifying **anchor pairs** -- cells from different conditions that are
   each other's mutual nearest neighbors in the shared CCA space
3. Using the anchors to compute a **correction vector** that aligns matching
   cell types across conditions


``` r
ifnb <- IntegrateLayers(object = ifnb,
                        method = CCAIntegration,
                        orig.reduction = "pca",
                        new.reduction = "integrated.cca")
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `method` | `CCAIntegration` | Use Canonical Correlation Analysis for integration. |
| `orig.reduction` | `"pca"` | The unintegrated PCA reduction to use as input. |
| `new.reduction` | `"integrated.cca"` | Name for the new integrated reduction that will be stored in the object. |

Now re-run the downstream steps using the **integrated** reduction instead
of the original PCA:


``` r
ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.5)
ifnb <- RunUMAP(ifnb, reduction = "integrated.cca", dims = 1:30)
```

Visualize the integrated result:


``` r
DimPlot(ifnb, reduction = "umap", group.by = "stim")
```
![UMAP plot after CCA integration. CTRL and STIM cells are now intermingled within each cluster, with teal (STIM) dominating the overlay. Cell types co-cluster regardless of condition.](fig/ep7_plot-after.png){alt="UMAP plot of the IFNB dataset after CCA integration showing control and stimulated cells intermingled within each cluster, indicating successful batch correction."}

The control and stimulated cells should now be intermingled within each
cluster. Cell types from both conditions co-cluster, which is what we want.

Let's verify with a split view:


``` r
DimPlot(ifnb, reduction = "umap", split.by = "stim")
```
![Two side-by-side UMAP panels split by condition (CTRL and STIM) after integration. Both panels show 14 clusters numbered 0 through 13, with the same spatial arrangement in each panel, confirming successful integration.](fig/ep7_split-view.png){alt="Split UMAP view showing CTRL and STIM conditions side by side after integration, with both panels displaying the same 14 clusters in identical spatial arrangement."}

Both panels should show the same overall structure, with the same clusters
present in both conditions. This confirms that integration successfully
aligned shared cell types.

Before proceeding to analysis, rejoin the layers so that expression data
from both conditions is accessible in a single matrix:


``` r
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
```

:::::::::::::::::::::::::::::::::::::::::: callout

## Other integration methods

Seurat v5 supports several integration methods through the same
`IntegrateLayers()` interface. You can swap `CCAIntegration` for any of these:

| Method | Function | Strengths |
|--------|:---------|:----------|
| **CCA** | `CCAIntegration` | Best for datasets with shared cell types. Robust to differences in cell type composition between batches. The default choice for most experiments. |
| **Harmony** | `HarmonyIntegration` | Very fast, scales well to large datasets (100k+ cells). Uses iterative soft clustering. Requires the `harmony` R package. Widely used and well-benchmarked. |
| **RPCA** | `RPCAIntegration` | Faster than CCA (uses reciprocal PCA instead of full CCA). Good when datasets are large or when CCA is too slow. May be less robust when cell type composition differs substantially between batches. |
| **scVI** | `scVIIntegration` | Deep learning approach. Excellent for complex batch structures (many batches, large composition differences). Requires Python and the `scvi-tools` package. Slower to run but often produces superior results on difficult datasets. |

For most standard experiments with two to five batches and shared cell types,
**CCA** or **Harmony** work well. Use **scVI** for more complex scenarios
(many batches, different tissues, cross-species integration).

:::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::: callout

## When NOT to integrate

Integration assumes that the same cell types exist in all batches and aligns
them. This is the wrong approach when:

- **All cells are expected to differ** between conditions. For example,
  comparing tumor cells to normal epithelial cells -- these are fundamentally
  different populations and should not be forced to co-cluster.
- **You want to find condition-specific cell types.** Integration pushes
  cells together, which can mask the appearance of a cell state that exists
  in only one condition.
- **The batch effect is minimal.** If cells already mix well by condition on
  the UMAP (as in some well-controlled experiments), integration is
  unnecessary and may introduce artifacts.

Always visualize the data **before integration** to assess whether batch
effects are actually present. If conditions already mix well, skip integration
and analyze directly.

:::::::::::::::::::::::::::::::::::::::::::::::::::


## Analysis on Integrated Data

Now that the data is integrated, we can annotate cell types and compare
between conditions.

### Annotating integrated clusters

Find markers for each integrated cluster:


``` r
ifnb.markers <- FindAllMarkers(ifnb,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25)
```

Annotate the clusters using known PBMC markers (the same markers from
Episode 6):


``` r
# Adjust cluster IDs to match YOUR results
new.cluster.ids <- c(
    "CD14 Mono",    # 0
    "CD4 T",        # 1
    "CD4 T",        # 2
    "CD4 T",        # 3
    "CD8 T",        # 4
    "NK",           # 5
    "B",            # 6
    "CD14 Mono",    # 7
    "FCGR3A Mono",  # 8
    "DC",           # 9
    "CD14 Mono",    # 10
    "Mk",           # 11
    "B",            # 12
    "Eryth"         # 13
)
names(new.cluster.ids) <- levels(ifnb)
ifnb <- RenameIdents(ifnb, new.cluster.ids)
ifnb$celltype <- Idents(ifnb)
```


``` r
DimPlot(ifnb, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
```
![UMAP plot with cells colored and labeled by annotated cell type. Nine cell types are visible: CD14 Mono (large cluster, right), CD4 T (large, upper left), CD8 T (upper right), B (center left), NK (left), FCGR3A Mono (lower right), DC (small, center), Mk (small, center right), and Eryth (small, bottom).](fig/ep7_annotated-umap.png){alt="Annotated UMAP plot of the integrated IFNB dataset with nine labeled cell types including CD14 Mono, CD4 T, CD8 T, B, NK, FCGR3A Mono, DC, Mk, and Eryth."}


### Conserved markers

**Conserved markers** are genes that are markers for a cell type in **both**
conditions. `FindConservedMarkers()` tests for differential expression within
each condition separately and then combines the results:


``` r
conserved <- FindConservedMarkers(ifnb,
                                  ident.1 = "CD14 Mono",
                                  grouping.var = "stim",
                                  only.pos = TRUE)


head(conserved, 10)
```

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ident.1` | `"CD14 Mono"` | The cell type to find markers for. |
| `grouping.var` | `"stim"` | The metadata column defining conditions. The test is run within each level of this variable. |
| `only.pos` | `TRUE` | Only return upregulated genes. |

The output has columns for each condition (e.g., `CTRL_avg_log2FC`,
`STIM_avg_log2FC`) showing that the marker is consistent across both. Genes
like LYZ, S100A9, and CD14 should appear as strong conserved markers for
monocytes.

:::::::::::::::::::::::::::::::::::::::::: spoiler

## Expected output

```
         CTRL_p_val CTRL_avg_log2FC CTRL_pct.1 CTRL_pct.2 CTRL_p_val_adj
TYROBP            0        2.357371      0.940      0.227              0
FCER1G            0        2.197611      0.930      0.233              0
S100A8            0        4.737642      0.742      0.071              0
C15orf48          0        2.719921      0.868      0.201              0
IL8               0        3.648031      0.798      0.144              0
CD63              0        2.739659      0.957      0.333              0
S100A9            0        4.597336      0.690      0.070              0
CTSB              0        3.108400      0.771      0.183              0
LGALS1            0        2.396676      0.904      0.322              0
TYMP              0        2.075884      0.859      0.283              0
            STIM_p_val STIM_avg_log2FC STIM_pct.1 STIM_pct.2 STIM_p_val_adj
TYROBP    0.000000e+00        2.344174      0.951      0.197   0.000000e+00
FCER1G    0.000000e+00        2.552006      0.910      0.196   0.000000e+00
S100A8    0.000000e+00        5.328514      0.444      0.020   0.000000e+00
C15orf48  0.000000e+00        2.758492      0.860      0.199   0.000000e+00
IL8      3.434056e-251        3.687341      0.294      0.028  4.825879e-247
CD63      0.000000e+00        2.896965      0.931      0.275   0.000000e+00
S100A9    0.000000e+00        5.195483      0.581      0.036   0.000000e+00
CTSB      0.000000e+00        3.421778      0.822      0.182   0.000000e+00
LGALS1    0.000000e+00        2.860965      0.866      0.216   0.000000e+00
TYMP      0.000000e+00        2.278326      0.927      0.437   0.000000e+00
              max_pval minimump_p_val
TYROBP    0.000000e+00              0
FCER1G    0.000000e+00              0
S100A8    0.000000e+00              0
C15orf48  0.000000e+00              0
IL8      3.434056e-251              0
CD63      0.000000e+00              0
S100A9    0.000000e+00              0
CTSB      0.000000e+00              0
LGALS1    0.000000e+00              0
TYMP      0.000000e+00              0
```

::::::::::::::::::::::::::::::::::::::::::::::::::


### Condition-specific differential expression

The most interesting biological question is: **which genes change between
control and stimulated cells within the same cell type?** To answer this, we
subset to a single cell type and compare conditions:


``` r
# Subset to CD14+ monocytes
mono <- subset(ifnb, idents = "CD14 Mono")

# Set the condition as the active identity
Idents(mono) <- "stim"

# Find DE genes between stimulated and control
mono.de <- FindMarkers(mono,
                       ident.1 = "STIM",
                       ident.2 = "CTRL",
                       min.pct = 0.25,
                       logfc.threshold = 0.25)
head(mono.de, 10)
```

:::::::::::::::::::::::::::::::::::::::::: spoiler

## Expected output

```
        p_val avg_log2FC pct.1 pct.2 p_val_adj
IFIT1       0   7.173024 0.965 0.033         0
IFIT3       0   6.798858 0.972 0.057         0
TNFSF10     0   6.370735 0.959 0.061         0
RSAD2       0   6.661976 0.920 0.041         0
IFIT2       0   6.905053 0.909 0.039         0
MX1         0   4.839203 0.953 0.099         0
CXCL10      0   8.007812 0.873 0.031         0
LY6E        0   4.263751 0.990 0.174         0
CXCL11      0   8.576537 0.801 0.011         0
OAS1        0   4.045792 0.915 0.131         0
```

The top DE genes are all **interferon-stimulated genes (ISGs)**: IFIT1, IFIT3,
TNFSF10, RSAD2, IFIT2. These are exactly the genes you would expect to be
upregulated by IFN-beta stimulation, confirming that our integration preserved
the biological signal while correcting the batch effect. The exact fold changes,
and the number of genes that pass a significance cutoff, vary slightly between
Seurat versions; the leading ISGs do not.

::::::::::::::::::::::::::::::::::::::::::::::::::



Save the annotated ifnb object for the next episode:


``` r
saveRDS(ifnb, file = "ifnb_annotated.rds")
```


### Visualizing condition-specific changes

Let's visualize a few ISGs to see the condition-specific response:


``` r
FeaturePlot(ifnb,
            features = c("IFIT3", "IFIT1"),
            split.by = "stim",
            cols = c("grey85", "firebrick"))
```

![Four-row FeaturePlot split by CTRL and STIM for ISG15, IFIT1, TNFSF10, and RSAD2. CTRL panels show mostly grey cells with sparse expression. STIM panels show strong red expression across monocyte and other clusters for all four genes.](fig/ep7_isg-feature.png){alt="FeaturePlot of ISG genes split by condition showing near-absent expression in CTRL panels and strong red expression in STIM panels, especially in monocyte clusters."}



``` r
VlnPlot(ifnb,
        features = c("IFIT1", "IFIT3", "TNFSF10", "RSAD2"),
        split.by = "stim",
        idents = "CD14 Mono",
        ncol = 4,
        pt.size = 0)
```

![Four violin plots for IFIT1, IFIT3, TNFSF10, and RSAD2 in CD14 monocytes. Each panel shows two violins: a thin line at zero for control and a broad teal shape peaking around expression level 3 to 4 for stimulated, confirming strong IFN-beta induction.](fig/ep7_isg-violin.png){alt="Violin plots of IFIT1, IFIT3, TNFSF10, and RSAD2 in CD14 monocytes showing flat expression in control and broad high expression in stimulated cells."}

IFIT1 and IFIT3 are strongly induced in the stimulated condition across
multiple cell types, with particularly strong expression in monocytes. This
demonstrates that integration corrected the batch effect (cell types co-cluster)
while preserving the biological effect of IFN-beta stimulation (ISGs are still
differentially expressed).


::::::::::::::::::::::::::::::::::::: challenge

## Challenge 1: Quantifying integration quality

After integration, check whether control and stimulated cells of the same type
are well mixed within each cluster. Use `DimPlot(split.by = "stim")` to
visualize, then compute the proportion of each condition per cluster.


``` r
# Split UMAP by condition
DimPlot(ifnb, reduction = "umap", split.by = "stim", label = TRUE, repel = TRUE) +
    NoLegend()

# Proportion table: condition per cell type
prop.table(table(Condition = ifnb$stim, CellType = Idents(ifnb)), margin = 2)
```

Are the proportions roughly balanced (close to 50/50) within each cell type?

:::::::::::::::::::::::: solution

## Solution

```
         CellType
Condition CD14 Mono     CD4 T     CD8 T        NK         B FCGR3A Mono
     CTRL 0.5033865 0.4326872 0.4775281 0.4287425 0.4812594   0.5378641
     STIM 0.4966135 0.5673128 0.5224719 0.5712575 0.5187406   0.4621359
         CellType
Condition        DC        Mk     Eryth
     CTRL 0.4473684 0.4561404 0.4583333
     STIM 0.5526316 0.5438596 0.5416667
```

![Split UMAP with CTRL and STIM panels, cells labeled by cell type. Both panels show identical cluster arrangement with CD4 T, CD8 T, B, NK, CD14 Mono, FCGR3A Mono, DC, Mk, and Eryth labels, confirming balanced condition mixing within each cell type.](fig/ep7_challenge1.png){alt="Split UMAP showing CTRL and STIM panels with cell type labels, demonstrating balanced mixing of conditions within each annotated cluster."}

Most cell types fall in the **0.43--0.57** range, confirming that integration
successfully mixed cells from both conditions. The split UMAP shows the same
cell types present in both panels with similar spatial arrangement -- each
cluster appears in both CTRL and STIM with comparable density.

Some observations worth discussing:

- **CD14 Mono** is the best-mixed cell type (50.3/49.7), essentially a perfect
  split. This is expected for the most abundant population.
- **CD4 T and NK** show the largest deviation (~43/57), with more cells in the
  stimulated condition. This could reflect a real biological shift -- IFN-beta
  stimulation may promote T cell and NK cell survival or proliferation in
  culture -- rather than a failure of integration.
- **Rare populations** (Mk, Eryth, DC) show reasonable balance despite small
  cell numbers, where stochastic variation has more impact.

The key diagnostic: if integration had failed, you would see a cell type that
is 90%+ from one condition, or you would see the same cell type forming
separate clusters in the two panels. Neither is the case here. The moderate
CD4 T skew (43/57) is well within the range of expected biological variation
between a control and stimulated sample.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::





::::::::::::::::::::::::::::::::::::: challenge

## Challenge 2: Interferon response genes in monocytes

Find the top 10 differentially expressed genes between stimulated and control
**CD14+ monocytes** (we computed this above as `mono.de`). Are these genes
consistent with an interferon stimulation response?


``` r
# Top 10 upregulated genes in stimulated vs control CD14+ monocytes
top10_de <- head(mono.de[order(mono.de$avg_log2FC, decreasing = TRUE), ], 10)
print(top10_de)

# Visualize the top 4 genes
FeaturePlot(ifnb,
            features = rownames(top10_de)[1:4],
            split.by = "stim",
            ncol = 4,
            cols = c("grey85", "firebrick"))
```

:::::::::::::::::::::::: solution

## Solution

![Four-row FeaturePlot split by CTRL and STIM showing CCL8, CXCL11, HESX1, and CXCL10. CTRL panels are almost entirely grey. STIM panels show strong red expression concentrated in the monocyte clusters for all four chemokine and ISG genes.](fig/ep7_challenge2.png){alt="FeaturePlot of CCL8, CXCL11, HESX1, and CXCL10 split by condition showing absent expression in control and strong monocyte-concentrated expression in stimulated cells."}


```
                 p_val avg_log2FC pct.1 pct.2     p_val_adj
CCL8      0.000000e+00   9.059454 0.798 0.015  0.000000e+00
CXCL11    0.000000e+00   8.576537 0.801 0.011  0.000000e+00
HESX1    3.530769e-218   8.283957 0.335 0.001 4.961789e-214
CXCL10    0.000000e+00   8.007812 0.873 0.031  0.000000e+00
IFIT1     0.000000e+00   7.173024 0.965 0.033  0.000000e+00
IFIT2     0.000000e+00   6.905053 0.909 0.039  0.000000e+00
IFIT3     0.000000e+00   6.798858 0.972 0.057  0.000000e+00
GMPR     2.606001e-242   6.750793 0.373 0.004 3.662213e-238
APOBEC3B 1.036577e-127   6.685079 0.272 0.034 1.456702e-123
RSAD2     0.000000e+00   6.661976 0.920 0.041  0.000000e+00
```

The top 10 genes fall into two functional categories, both consistent with
IFN-beta stimulation:

| Gene | log2FC | pct.1 / pct.2 | Function |
|------|:------:|:--------------|:---------|
| CCL8 | 9.06 | 0.80 / 0.02 | Chemokine; recruits monocytes and T cells to sites of inflammation |
| CXCL11 | 8.58 | 0.80 / 0.01 | IFN-inducible chemokine; attracts activated T cells via CXCR3 |
| HESX1 | 8.28 | 0.34 / 0.00 | Transcription factor; emerging role in IFN signaling |
| CXCL10 | 8.01 | 0.87 / 0.03 | IFN-inducible chemokine (IP-10); major recruiter of immune cells |
| IFIT1 | 7.17 | 0.97 / 0.03 | Directly induced by type I IFN; inhibits viral translation |
| IFIT2 | 6.91 | 0.91 / 0.04 | Antiviral effector; regulates translation during infection |
| IFIT3 | 6.80 | 0.97 / 0.06 | Antiviral effector; forms complex with IFIT1 and IFIT2 |
| GMPR | 6.75 | 0.37 / 0.00 | Guanosine monophosphate reductase; modulates purine metabolism |
| APOBEC3B | 6.69 | 0.27 / 0.03 | Cytidine deaminase; innate antiviral defense via RNA editing |
| RSAD2 | 6.66 | 0.92 / 0.04 | Viperin; broad-spectrum antiviral enzyme induced by IFN |

Two patterns stand out:

**Chemokines dominate the top.** CCL8, CXCL11, and CXCL10 are three of the
top four, with fold changes of 8.0-9.1, reflecting the monocyte's primary role as a sentinel
cell: upon IFN stimulation, monocytes broadcast chemokine signals to recruit
other immune cells. Note the pct.2 values (0.01-0.03) -- these chemokines are
essentially absent in control cells and massively induced by stimulation.

**Classical ISGs fill the rest.** IFIT1/2/3, RSAD2, and APOBEC3B are
canonical interferon-stimulated genes that directly inhibit viral replication.
The IFIT family shows the highest pct.1 values (0.91-0.97), meaning nearly
every stimulated monocyte expresses them -- a hallmark of a robust, uniform
interferon response.

The FeaturePlot confirms this visually: CCL8, CXCL11, HESX1, and CXCL10 are
virtually absent in the CTRL panels (gray) and strongly expressed in the STIM
panels (red), with the signal concentrated in the monocyte clusters (right
side of the UMAP). CXCL10 and CXCL11 also show some expression in the
stimulated CD4 T and NK clusters, consistent with these cell types also
responding to IFN-beta, though less strongly than monocytes.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::: keypoints

- Batch effects cause cells to cluster by sample rather than by cell type, preventing cross-condition comparisons
- CCA integration finds shared correlation structure between conditions and aligns matching cell types using anchor pairs
- In Seurat v5, multi-sample data is handled by splitting assay layers, integrating, then rejoining layers for downstream analysis
- Always verify integration by checking that conditions mix within clusters AND that known biological differences (like ISG expression) are preserved
- Integration is not always appropriate; skip it when conditions have no shared cell types or when batch effects are minimal

::::::::::::::::::::::::::::::::::::::::::::::::
