---
title: Reference
---

## Additional Reading

- **Hao, Y. et al. (2024).** Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nature Biotechnology*, 42, 293-304. [doi:10.1038/s41587-023-01767-y][seurat-v5-paper] - The Seurat v5 paper describing the dictionary learning framework and layer-based architecture used in this workshop.

- **Luecken, M.D. & Theis, F.J. (2019).** Current best practices in single-cell RNA-seq analysis: a tutorial. *Molecular Systems Biology*, 15(6), e8746. [doi:10.15252/msb.20188746][best-practices-paper] - A comprehensive overview of best practices for scRNA-seq analysis, covering QC, normalization, dimensionality reduction, clustering, and more.

- **Kaminow, B., Yunusov, D., & Dobin, A. (2021).** STARsolo: accurate, fast and versatile mapping/quantification of single-cell and single-nucleus RNA-seq data. *bioRxiv*. [doi:10.1101/2021.05.05.442755][starsolo-paper] - The STARsolo paper describing the open-source alternative to Cell Ranger used in Episode 2.

- **Aran, D. et al. (2019).** Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. *Nature Immunology*, 20, 163-172. [doi:10.1038/s41590-018-0276-y][singler-paper] - The SingleR paper describing the reference-based automated annotation method used in Episode 6.

## Glossary

**Barcode**
: A short DNA sequence (typically 16 nucleotides in 10x Genomics) that uniquely identifies a cell in a droplet-based scRNA-seq experiment.

**Batch effect**
: Systematic technical differences between samples processed at different times, on different lanes, or with different reagents that are unrelated to the biology of interest.

**CCA (Canonical Correlation Analysis)**
: A statistical method used by Seurat for data integration that identifies shared sources of variation between two datasets.

**Cell Ranger**
: A set of analysis pipelines from 10x Genomics for processing Chromium single-cell data, including alignment, barcode counting, and gene expression quantification.

**Cluster**
: A group of cells that share similar gene expression profiles, typically identified through graph-based community detection algorithms.

**Count matrix**
: A genes-by-cells matrix where each entry represents the number of UMI counts detected for a given gene in a given cell.

**Differential expression (DE)**
: Statistical testing to identify genes whose expression levels differ significantly between two or more groups of cells.

**Dimensionality reduction**
: Mathematical techniques (PCA, UMAP, t-SNE) that reduce high-dimensional gene expression data to a lower-dimensional representation for visualization and analysis.

**Doublet**
: An artifact where two cells are captured in the same droplet and assigned the same barcode, producing a hybrid expression profile.

**Feature**
: In scRNA-seq, typically refers to a gene. Feature selection identifies genes with the most variation across cells.

**GEM (Gel Bead-in-Emulsion)**
: A droplet in the 10x Chromium system containing a single gel bead (with barcode oligonucleotides) and ideally a single cell.

**Highly variable genes (HVGs)**
: Genes that show the most variation in expression across cells, selected as informative features for dimensionality reduction and clustering.

**Leiden algorithm**
: A community detection algorithm used for clustering cells in a shared nearest neighbor graph. An improvement over the Louvain algorithm.

**LogNormalize**
: A normalization method that divides gene counts by total counts per cell, multiplies by a scale factor (10,000), and log-transforms the result.

**Marker gene**
: A gene whose expression is specific to or enriched in a particular cell type, used for cell type annotation.

**PCA (Principal Component Analysis)**
: A linear dimensionality reduction technique that identifies the axes of greatest variance in the data.

**Pseudobulk**
: An approach that aggregates single-cell expression profiles by sample and cell type, creating pseudo-bulk samples for statistically valid differential expression analysis.

**QC (Quality Control)**
: The process of identifying and removing low-quality cells based on metrics such as total counts, number of genes detected, and mitochondrial gene percentage.

**SCTransform**
: A normalization method that uses regularized negative binomial regression to stabilize variance across genes with different expression levels.

**Seurat object**
: The primary data structure in Seurat that stores expression data, cell metadata, dimensionality reductions, and analysis results.

**SingleR**
: An R package for automated cell type annotation that assigns labels by comparing expression profiles to labeled reference datasets.

**SNN (Shared Nearest Neighbor) graph**
: A graph where cells are nodes and edges connect cells that share nearest neighbors in PCA space, used as input for clustering algorithms.

**STARsolo**
: An open-source tool integrated into the STAR aligner for processing droplet-based scRNA-seq data, serving as a fast alternative to Cell Ranger.

**UMAP (Uniform Manifold Approximation and Projection)**
: A nonlinear dimensionality reduction technique used primarily for 2D visualization of single-cell data.

**UMI (Unique Molecular Identifier)**
: A short random DNA sequence (typically 12 nucleotides in 10x Genomics) attached to each mRNA molecule before amplification, used to count unique transcripts and remove PCR duplicates.

[seurat-v5-paper]: https://doi.org/10.1038/s41587-023-01767-y
[best-practices-paper]: https://doi.org/10.15252/msb.20188746
[starsolo-paper]: https://doi.org/10.1101/2021.05.05.442755
[singler-paper]: https://doi.org/10.1038/s41590-018-0276-y
