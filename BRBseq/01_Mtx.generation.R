# Load dependencies:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
# Load Seurat package
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("Seurat"))

source("BRBseq/config.R")

# Create output.dir
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

# Read samples plan, barcode file and merge both
samples.plan <- read.delim(samples.plan.file)
barcode <- read.delim(barcode.well.file)
samples.plan.with.barcode <- merge(samples.plan, barcode)
# Reorder samples plan
samples.plan.with.barcode <- samples.plan.with.barcode[match(samples.plan$sample, samples.plan.with.barcode$sample), ]
samples.plan <- samples.plan.with.barcode

# Read the matrix which is stored as sparse
mtx <- ReadMtx(
  list.files(input.matrix.dir, pattern = "mtx$", full.names = TRUE),
  list.files(input.matrix.dir, pattern = "Barcodes", full.names = TRUE),
  list.files(input.matrix.dir, pattern = "Genes", full.names = TRUE),
  cell.column = 1,
  # feature.column = 2,
  feature.column = 1,
  cell.sep = "\t",
  feature.sep = "\t",
  skip.cell = 0,
  skip.feature = 0,
  mtx.transpose = FALSE,
  unique.features = TRUE,
  strip.suffix = FALSE
)
# Check the number of UMI per barcode
colSums(mtx)
table(colSums(mtx) > 1000)
# FALSE  TRUE 
#    54    42

# Only keep columns in samples.plan
mtx <- mtx[, samples.plan$B1]

# Relabel column names
colnames(mtx) <- samples.plan$sample

# Check they all have above 1k UMIs
table(colSums(mtx) < 1000)
# FALSE 
#    21

head(samples.plan)

# Generate count.table to export
count.table <- cbind(data.frame(gene_id = rownames(mtx)), mtx)
write.table(count.table, file.path(output.dir, "all_counts.txt"), sep = "\t", quote = F, row.names = F)

# Normalize by million number of UMI per cell
big.mat.pm <- t(t(as.matrix(mtx)) / colSums(mtx) * 1e6)
rownames(big.mat.pm) <- rownames(mtx)
count.table.pm <- cbind(data.frame(gene_short_name = rownames(big.mat.pm)), big.mat.pm)
write.table(count.table.pm, norm.counts.file, sep = "\t", quote = F, row.names = F)

