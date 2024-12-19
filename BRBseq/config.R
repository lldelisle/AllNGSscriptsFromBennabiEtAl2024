# Define variables
# >>> from pylab import *
# >>> cmap = cm.get_cmap('gist_earth', 8)
# >>> for i in range(cmap.N):
# ...     rgba = cmap(i)
# ...     print(matplotlib.colors.rgb2hex(rgba))
# ...
# #000000
# #184279
# #30807f
# #429552
# #82a954
# #b9b35f
# #d4b093
# #fdfbfb
my.colors <- c("#000000", "#184279", "#30807f", "#429552", "#82a954", "#b9b35f", "#d4b093")
names(my.colors) <- c(50, 100, 300, 600, 1200, 1800, 3000)

my.cond.colors <- c("CTL" = "#00BFC4", "noGSKi" = "#D1D3D4")
# Not used:
# my.sign.colors <- c("down" = "#F8766D", "up" = "#00BFC4")

# Define paths:
samples.plan.file <- "BRBseq/inputs/SamplePlan.txt"
output.dir <- "output.files/BRBseq/"
barcode.well.file <- "BRBseq/inputs/barcodes_96_V5A_brb.txt"
input.matrix.dir <- "BRBseq/STAR_solo_results/"
norm.counts.file <- file.path(output.dir, "all_counts_pm.txt")
counts.file <- file.path(output.dir, "all_counts.txt")
pathForDESeq2 <- file.path(output.dir, "DESeq2_pairwise")

# Define thresholds
log2FC.threshold <- log2(1.5)
