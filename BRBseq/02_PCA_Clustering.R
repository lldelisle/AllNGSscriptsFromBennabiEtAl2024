# Load dependencies:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
# Load ggplot2 and pheatmap package
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("ggplot2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("scales"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("pheatmap"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("ggforce"))

source("BRBseq/config.R")

# Define paths:
output.dir <- file.path(output.dir, "general_plots")

# Define variables
ns.var.genes <- 1000

# Create output.dir
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

# Read files
samples.plan.df <- read.delim(samples.plan.file, check.names = FALSE)
expression.df <- read.delim(norm.counts.file, check.names = FALSE)

# Only keep colors in samples.plan.df
my.colors <- my.colors[intersect(names(my.colors), unique(samples.plan.df$NumberOfCells))]

# Add N_0=
samples.plan.df$NumberOfCellsFormatted <- paste0("N[0]~'= ", samples.plan.df$NumberOfCells, "'")
samples.plan.df$sampleFormatted <- paste0("N[0]~'= ", samples.plan.df$NumberOfCells, "_", samples.plan.df$Rep, "'")
samples.plan.df$sampleFormatted[samples.plan.df$Condition == "noGSKi"] <- paste0("N[0]~'= ", samples.plan.df$NumberOfCells[samples.plan.df$Condition == "noGSKi"], "_noGSKi_", samples.plan.df$Rep[samples.plan.df$Condition == "noGSKi"], "'")
my.colors.formatted <- my.colors
names(my.colors.formatted) <- samples.plan.df$NumberOfCellsFormatted[match(names(my.colors), samples.plan.df$NumberOfCells)]

# Transform NumberOfCells in factor:
samples.plan.df$NumberOfCells <- factor(samples.plan.df$NumberOfCells, levels = names(my.colors))
samples.plan.df$NumberOfCellsFormatted <- factor(samples.plan.df$NumberOfCellsFormatted, levels = names(my.colors.formatted))

# Remove non-expressed genes
data <- expression.df[, samples.plan.df$sample]
sumperline <- apply(data, 1, sum)
nonZdata <- data[sumperline != 0, ]
ldata <- log2(nonZdata + 1)

# Define sample order to follow metadata rather than well order
my.samples.order <- samples.plan.df[order(samples.plan.df$Condition, samples.plan.df$NumberOfCells, samples.plan.df$Rep), "sample"]

for (n.var.genes in ns.var.genes) {
  # Restrict to variable genes
  rldata <- ldata[order(apply(ldata, 1, var), decreasing = T)[1:min(nrow(ldata), n.var.genes )], ]
  # Compute PCA centered not scaled
  sample.pca <- prcomp(t(rldata), center = TRUE, scale. = FALSE)
  # Join with samples plan
  new.df <- data.frame(samples.plan.df, sample.pca$x[samples.plan.df$sample, ])
  # Get variable percentage
  var <- round((sample.pca$sdev)^2/sum(sample.pca$sdev^2) * 100)
  # Plot PCA 1vs2
  g <- ggplot(new.df, aes(PC1, PC2, color = NumberOfCellsFormatted)) +
    geom_point(aes(shape = Condition), size = 3) +
    theme_classic(base_size = 20) +
    xlab(paste0("PC1:", var[1], "% variance")) +
    ylab(paste0("PC2:", var[2], "% variance")) +
    scale_color_manual("Number of cells",
                       values = my.colors.formatted,
                       label = parse_format()) +
    ggtitle(paste0("PCA for ", n.var.genes, " variable genes"))
  ggsave(file.path(output.dir, paste0("PC1-PC2_", n.var.genes, "VarGenes.pdf")),
    g,
    width = 7, height = 5, units = "in"
  )
  correlation.mtx <- cor(rldata, method = "pearson")
  sample.dist <- as.dist(1 - correlation.mtx)
  clu <- hclust(sample.dist, method = "ward.D2")
  dd <- as.dendrogram(clu)
  dd2 <- reorder(dd, match(rownames(correlation.mtx), my.samples.order), agglo.FUN = min)
  annot <- subset(samples.plan.df, select = c("NumberOfCells", "Condition"))
  rownames(annot) <- samples.plan.df$sample
  pheatmap(
    correlation.mtx,
    cluster_rows = as.hclust(dd2),
    cluster_cols = as.hclust(dd2),
    cellwidth = 10,
    cellheight = 10,
    annotation = annot,
    annotation_colors = list(
    NumberOfCells = my.colors,
      Condition = my.cond.colors
    ),
    labels_row = parse(text = samples.plan.df$sampleFormatted[match(rownames(correlation.mtx), samples.plan.df$sample)]),
    show_colnames = FALSE,
    main = paste("spearmanCor - ward clustering -", n.var.genes, "genes"),
    clustering_method = "ward.D2",
    filename = file.path(output.dir, paste0("clustering_",n.var.genes,"VarGenes.pdf"))
  )
}
