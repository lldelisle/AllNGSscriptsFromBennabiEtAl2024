# Single-cell RNA-seq analysis

All samples are described in [this table](./metadata.csv).

## FASTQ to matrices

The first steps of single-cell RNA-seq analysis was performed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

The workflow used has been exported [here](./scRNA-seq_preprocessing_10X_cellPlex_UPDUB.ga). They have been run with the following inputs:

- fastqPE collection GEX: fastqs are available on GEO/SRA.
- reference genome: `mm10` (from UCSC)
- gtf: downloaded from [zenodo](https://zenodo.org/records/10079673/files/mm10_custom102_allGastruloids_min10_extended.gtf.gz)
- cellranger_barcodes_3M-february-2018.txt: downloaded from [zenodo](https://zenodo.org/record/3457880/files/3M-february-2018.txt.gz).
- fastqPE collection CMO: fastqs are available on GEO/SRA
- cmo_10X_seq.txt: available [here](./cmo_10X_seq.txt)
- sample name and CMO sequence collection: see [CMO_samples](./CMO_samples) directory
- Number of expected cells (used by CITE-seq-Count): `24000`

The command lines used have been written [here](./scRNA-seq_preprocessing_CL.sh).

## Matrices to plots

### R configuration

The scripts have been launched on a RStudio server.

The session information is available [here](./sessionInfo.txt).

### Pipeline

All samples have been registered into a [csv file](./metadata.csv).

Once matrices have been generated in galaxy they were downloaded and organized like this:

```bash
.
├── CMO
│   └── Directory1
│       ├── barcodes.tsv
│       ├── genes.tsv
│       └── matrix.mtx
│   └── Directory2
│       ├── barcodes.tsv
│       ├── genes.tsv
│       └── matrix.mtx
└── GEX
    ├── Directory1
    │   ├── barcodes.tsv
    │   ├── genes.tsv
    │   └── matrix.mtx
    └── Directory2
        ├── barcodes.tsv
        ├── genes.tsv
        └── matrix.mtx
```

where Directory1 and Directory2 would be in the `Directory` column of the csv file.

[Step1](./Step1.Seurat.Demultiplexing.Analysis.R) has been run on all available samples. This R script generates a RDS for each sample. It demultiplexes CellPlex experiments.

[Step2](./Step2.Seurat.Analysis.and.Merging.R) has been run to generate a merged Seurat object with all samples.

[Step3 for WT](./Step3.qmd) contains all code used to generate figures.

All common functions useful for single-cell RNA-seq analysis to have been collected into a [single file](./scRNAseqFunctions.R).
