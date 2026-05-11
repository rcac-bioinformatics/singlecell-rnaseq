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

## Setup


``` r
library(Seurat)
```

``` error
Error in `library()`:
! there is no package called 'Seurat'
```

``` r
library(SeuratData)
```

``` error
Error in `library()`:
! there is no package called 'SeuratData'
```

``` r
library(ggplot2)
```

``` error
Error in `library()`:
! there is no package called 'ggplot2'
```

``` r
library(patchwork)
```

``` error
Error in `library()`:
! there is no package called 'patchwork'
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

```output
An object of class Seurat
14053 features across 13999 samples within 1 assay
Active assay: RNA (14053 features, 0 variable features)
 1 layer present: counts
```

The dataset has approximately 14,000 cells and 14,000 genes. The `stim` column
in the metadata records which condition each cell came from.

::::::::::::::::::::::::::::::::::::::::::::::::::


``` r
table(ifnb$stim)
```

:::::::::::::::::::::::::::::::::::::::::: spoiler

## Expected output

```output
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

```output
An object of class Seurat
14053 features across 13999 samples within 1 assay
Active assay: RNA (14053 features, 0 variable features)
 2 layers present: counts.CTRL, counts.STIM
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

The control and stimulated cells should now be intermingled within each
cluster. Cell types from both conditions co-cluster, which is what we want.

Let's verify with a split view:


``` r
DimPlot(ifnb, reduction = "umap", split.by = "stim")
```

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

```output
              p_val avg_log2FC pct.1 pct.2    p_val_adj
IFIT1  0.000e+00      4.528 0.992 0.056 0.000e+00
ISG15  0.000e+00      4.345 0.995 0.186 0.000e+00
IFI6   0.000e+00      4.123 0.968 0.042 0.000e+00
ISG20  0.000e+00      3.159 0.940 0.124 0.000e+00
MX1    0.000e+00      3.023 0.893 0.029 0.000e+00
IFIT3  0.000e+00      3.752 0.959 0.048 0.000e+00
IFI44L 0.000e+00      3.697 0.888 0.022 0.000e+00
IFIT2  0.000e+00      3.583 0.919 0.027 0.000e+00
LY6E   0.000e+00      3.243 0.919 0.077 0.000e+00
RSAD2  0.000e+00      3.443 0.828 0.013 0.000e+00
```

The top DE genes are all **interferon-stimulated genes (ISGs)**: IFIT1, ISG15,
IFI6, MX1, IFIT3. These are exactly the genes you would expect to be
upregulated by IFN-beta stimulation, confirming that our integration preserved
the biological signal while correcting the batch effect.

::::::::::::::::::::::::::::::::::::::::::::::::::

### Visualizing condition-specific changes

Let's visualize a few ISGs to see the condition-specific response:


``` r
FeaturePlot(ifnb,
            features = c("ISG15", "IFIT1"),
            split.by = "stim",
            cols = c("grey85", "firebrick"))
```


``` r
VlnPlot(ifnb,
        features = "ISG15",
        split.by = "stim",
        idents = "CD14 Mono",
        pt.size = 0)
```

ISG15 and IFIT1 are strongly induced in the stimulated condition across
multiple cell types, with particularly strong expression in monocytes. This
demonstrates that integration corrected the batch effect (cell types co-cluster)
while preserving the biological effect of IFN-beta stimulation (ISGs are still
differentially expressed).


::::::::::::::::::::::::::::::::::::: challenge

## Challenge 1: Quantifying Integration Quality

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

The proportion table should show values close to **0.45--0.55** for each
condition within each cell type, indicating good mixing. For example:

```output
          CellType
Condition  CD14 Mono CD4 T CD8 T   NK    B   FCGR3A Mono  DC    Mk
  CTRL         0.47  0.48  0.47 0.46 0.48         0.43 0.44  0.50
  STIM         0.53  0.52  0.53 0.54 0.52         0.57 0.56  0.50
```

The proportions are close to 50/50 within each cell type, confirming that
integration successfully mixed cells from both conditions. Small deviations
from 50/50 are normal and may reflect real differences in cell type abundance
between conditions (e.g., slightly more monocytes in the stimulated condition).

If you see a cell type that is heavily skewed (e.g., 90% from one condition),
that could indicate:

- A condition-specific cell population that should not have been integrated
- Poor integration for that particular cell type
- A real biological shift in cell type composition

The split UMAP should show the same clusters in both panels, with similar
density patterns.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: challenge

## Challenge 2: Interferon Response Genes in Monocytes

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

The top 10 upregulated genes should include several well-known
**interferon-stimulated genes (ISGs)**:

| Gene | Function | Evidence of IFN response |
|------|:---------|:------------------------|
| IFIT1 | Interferon-induced protein with tetratricopeptide repeats | Direct target of type I interferon signaling |
| ISG15 | Ubiquitin-like modifier | Conjugated to target proteins during antiviral response |
| IFI6 | Interferon alpha-inducible protein 6 | Inhibits apoptosis in response to viral infection |
| MX1 | Myxovirus resistance protein 1 | GTPase that blocks viral replication |
| IFIT3 | Interferon-induced protein with tetratricopeptide repeats 3 | Antiviral effector |
| IFI44L | Interferon-induced protein 44-like | Induced by type I interferons |
| IFIT2 | Interferon-induced protein with tetratricopeptide repeats 2 | Antiviral, regulates translation |
| ISG20 | Interferon-stimulated exonuclease gene 20 | Degrades viral RNA |
| LY6E | Lymphocyte antigen 6E | Modulates viral entry |
| RSAD2 (Viperin) | Radical SAM domain-containing 2 | Broad-spectrum antiviral enzyme |

Every gene in the top 10 is a known component of the **type I interferon
response pathway**. This makes perfect biological sense: IFN-beta is a type I
interferon, and stimulating PBMCs with it activates the JAK-STAT signaling
cascade, which induces transcription of hundreds of ISGs.

The FeaturePlot should show these genes as nearly absent in the CTRL panel
and strongly expressed across multiple cell types (especially monocytes) in
the STIM panel. This confirms that:

1. Integration successfully corrected technical batch effects
2. The biological response to IFN-beta stimulation was preserved
3. We can detect meaningful condition-specific changes after integration

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints

- Batch effects cause cells to cluster by sample rather than by cell type, preventing cross-condition comparisons
- CCA integration finds shared correlation structure between conditions and aligns matching cell types using anchor pairs
- In Seurat v5, multi-sample data is handled by splitting assay layers, integrating, then rejoining layers for downstream analysis
- Always verify integration by checking that conditions mix within clusters AND that known biological differences (like ISG expression) are preserved
- Integration is not always appropriate; skip it when conditions have no shared cell types or when batch effects are minimal

::::::::::::::::::::::::::::::::::::::::::::::::
