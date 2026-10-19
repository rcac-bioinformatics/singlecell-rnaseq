Single-Cell RNA-Seq Analysis: Hands-on Training, workshop data (Negishi)
Lesson: https://rcac-bioinformatics.github.io/singlecell-rnaseq

This directory is read-only. Copy it to your scratch space and work there:

    mkdir -p ${RCAC_SCRATCH}/scrna_workshop
    rsync -avP /depot/workshop/data/scrna_workshop/ ${RCAC_SCRATCH}/scrna_workshop/

What is here, and which episode uses it
---------------------------------------

fastq/
    Episode 2 input. 10x Genomics PBMC 10k, 3' v3 chemistry, two lanes:
    pbmc_10k_v3_S1_L001_{I1,R1,R2}_001.fastq.gz and the same for L002.
    R1 holds the cell barcode and UMI, R2 the cDNA read. I1 is the sample
    index; Cell Ranger and STARsolo do not need it for this dataset.

reference/refdata-gex-GRCh38-2024-A/
    Episode 2 input. The 10x human reference GRCh38 2024-A: genome FASTA,
    gene annotation (genes/genes.gtf.gz) and the Cell Ranger STAR index.

filtered_feature_bc_matrix/
    Episode 3 input. The filtered count matrix that the episode 2 Cell Ranger
    script produces: 38,606 genes x 11,809 cell barcodes. Episode 3 reads
    this copy, so you do not have to wait for your own Cell Ranger job.

rds/
    Optional saved R objects, for starting an episode without finishing the
    one before it. Each episode reads its input from scrna_workshop/, so copy
    the file you need up one level first (in a terminal, from your
    scrna_workshop directory: cp rds/pbmc_filtered.rds .).

        pbmc_filtered.rds    saved by episode 3, read by episode 4
        pbmc_normalized.rds  saved by episode 4, read by episode 5
        pbmc_clustered.rds   saved by episode 5, read by episode 6
        pbmc_annotated.rds   saved by episode 6, read by episode 8 (export)
        ifnb_annotated.rds   saved by episode 7, read by episode 8

Not here on purpose
-------------------

- Episode 2 outputs (cellranger_output/, starsolo_output/,
  reference/star_index/, the 3M barcode whitelist): you create these
  yourself in episode 2.
- The ifnb dataset for episodes 7 and 8: SeuratData downloads it in
  episode 7.

Where the files come from
-------------------------

fastq/       https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar
reference/   https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
             (both fetched FILL-IN-DATE, archive size and S3 ETag checked)
filtered_feature_bc_matrix/
             FILL-IN-MATRIX-SOURCE
rds/         FILL-IN-RDS-SOURCE

Questions: ask an instructor, or rcac-help@purdue.edu.
