# BRBseq

## Samples

The samples are described in the [samples plan](./inputs/SamplePlan.txt).

## FASTQ to matrix

The preprocessing was computed in a local Galaxy server ([https://doi.org/10.1093/nar/gkae410](https://doi.org/10.1093/nar/gkae410)).

The workflow is [here](./BRB-seq.ga).

The parameters were:

- genome: mm10 (from UCSC)
- gtf (exteded for 3'UTR): [mm10_custom102_allGastruloids_min10_extended.gtf](https://zenodo.org/records/10079673/files/mm10_custom102_allGastruloids_min10_extended.gtf.gz)
- fastqs : will be available on GEO/SRA
- barcodes_96_V5A_brb.txt: [here](./inputs/barcodes_96_V5A_brb.txt)

The command lines and tool versions can be found [here](BRBseq/inputs/command_lines.sh).

The outputs are [here](./STAR_solo_results/).

Then, a [R script](./01_Mtx.generation.R) was used to generate the count matrices: [raw](./outputs/all_counts.txt) or [normalized](./outputs/all_counts_pm.txt):

```bash
Rscript BRBseq/01_Mtx.generation.R
```

## General QC

The PCA and clustering was performed with [this R script](./02_PCA_Clustering.R)

```bash
Rscript BRBseq/02_PCA_Clustering.R
```

The plots are [here](./outputs/general_plots/).

## Differential Expression Analysis

DEA was computed using DESeq2 using [this R script](./03_DESeq2.R)

```bash
Rscript BRBseq/03_DESeq2.R
```

All results are in [this directory](./outputs/DESeq2_pairwise/).

Then a deeper analysis was performed on genes that were deregulated in the 2 exterme conditions: 50 and 1800 cells compared to the 300 cells.

The R script is [here](./04_clustering.R) and the results are [here](./outputs/DESeq2_pairwise/).

```bash
Rscript BRBseq/04_clustering.R
```
