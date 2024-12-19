# Load dependencies:
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)

suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("rtracklayer"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("DESeq2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("ggplot2"))
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("eulerr"))

source("BRBseq/config.R")

# Functions:
# Wrapper for DESeq2
deseqAnaWithCovariates <- function(count.table, factorForAna, covariates,
                                   pathOutput, samplesPlan, LRT=F, pvalT=0.05,
                                   lfcT=1.5, writeRLOG=F, gene_id = "gene_id", ...) {
  # Checking the conditions
  if (!(factorForAna %in% colnames(samplesPlan))) {
    stop("The factor is is not part of the column names.")
  }
  if (!is.null(covariates) & !(all(unlist(covariates) %in% colnames(samplesPlan)))) {
    stop("Not all covariates are part of the column names.")
  }
  if (length(levels(samplesPlan[, factorForAna])) == 1) {
    stop("The factor you chose have only 1 value. The analysis is not possible.")
  }
  if (length(levels(samplesPlan[, factorForAna])) > 2 & !LRT) {
    print("The factor you chose have more than 2 values. LRT will be applied.")
    LRT <- T
  }
  # Lauching DESeq2
  dds <- DESeqDataSetFromMatrix(countData = count.table[, match(rownames(samplesPlan), colnames(count.table))],
                                colData = samplesPlan,
                                design = as.formula(paste0("~", paste(c(unlist(covariates), factorForAna), collapse = " + "))))
  print("Design is:")
  print(design(dds))
  print("Genes that are never expressed are removed")
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  if (LRT) {
    # Here I am really not sure about the reduced
    reduced.formula <- as.formula("~1")
    if (!is.null(covariates)) {
      reduced.formula <- as.formula(paste0("~", paste(unlist(covariates), collapse = " + ")))
    }
    dds <- DESeq(dds, minReplicatesForReplace = Inf, test = "LRT", reduced = reduced.formula)
  } else {
    dds <- DESeq(dds, minReplicatesForReplace = Inf)
  }
  res <- results(dds, ...)
  resOrdered <- res[order(res$padj), ]
  # Subsetting the annotation file
  ann <- subset(count.table, select = intersect(colnames(count.table), c("gene_id", "gene_short_name", "locus")))
  annot.df <- data.frame(ann[match(rownames(resOrdered), ann[, gene_id]), ])
  colnames(annot.df) <- gene_id
  resToExport <- data.frame(annot.df, counts(dds, normalized = TRUE)[rownames(resOrdered), ], resOrdered)
  write.table(resToExport, file = paste(pathOutput, "DESeq2Results.txt", sep = ''), sep = '\t', row.names = F, quote = F)
  dfDiffExp <- subset(resToExport, resToExport$padj < pvalT & abs(resToExport$log2FoldChange) > lfcT)
  write.table(dfDiffExp, file = paste0(pathOutput, "DESeq2significant.txt"), sep = '\t', row.names = F, quote = F)
  rld <- rlog(dds)
  rlogdata <- assay(rld)
  if (writeRLOG) {
    resToExport2 <- data.frame(annot.df, rlogdata[rownames(resOrdered), ])
    write.table(resToExport2, file = paste0(pathOutput, "rlog.txt"), sep = '\t', row.names = F, quote = F)
  }
  return(invisible(dfDiffExp))
}
# Multiple intersect
mintersect <- function(list_l) {
  if (length(list_l) < 2) {
    return(list_l[[1]])
  }
  cur.intersect <- intersect(list_l[[1]], list_l[[2]])
  if (length(list_l) > 2) {
    for (i in 3:length(list_l)) {
      cur.intersect <- intersect(cur.intersect, list_l[[i]])
    }
  }
  return(cur.intersect)
}

# Fixed variables:
path <- "BRBseq/"
gtf.file <- file.path(path, "inputs", "mm10_custom102_allGastruloids_min10_extended.gtf.gz")
gtf.url <- "https://zenodo.org/records/10079673/files/mm10_custom102_allGastruloids_min10_extended.gtf.gz"
gene_id <- "gene_id"

# Prepare inputs
samples.plan.df <- read.delim(samples.plan.file, check.names = FALSE)
rownames(samples.plan.df) <- samples.plan.df$sample

# Only keep colors in samples.plan.df
my.colors <- my.colors[intersect(names(my.colors), unique(samples.plan.df$NumberOfCells))]

# Reorder factors
samples.plan.df$NumberOfCells <- factor(samples.plan.df$NumberOfCells, levels = c(300, unique(setdiff(unique(samples.plan.df$NumberOfCells), 300))))
samples.plan.df$Condition <- factor(samples.plan.df$Condition, levels = c("CTL", "noGSKi"))
all.analyses <- list("NumberOfCells" = list("NumberOfCells" = list("Condition" = "CTL")),
                     "Condition" = list("Condition" = list("NumberOfCells" = "300")))

count.table <- read.delim(counts.file, check.names = FALSE)
colnames(count.table)[1] <- gene_id

big.table.fn.long <- "summary_long.txt"

if (!file.exists(file.path(pathForDESeq2, "summary.txt"))) {
  # Prepare a big table with the results of all DESeq2
  if (! file.exists(gtf.file)) {
    download.file(gtf.url, gtf.file)
  }
  gtf <- readGFF(gtf.file)
  big.annot <- unique(gtf[, c("gene_id", "gene_name", "seqid", "gene_biotype")])
  colnames(big.annot) <- c("gene_id", "gene_short_name", "chr", "gene_biotype")
  ##### !  WARNING #### I am using only protein coding genes
  
  big.annot <- subset(big.annot, gene_biotype == "protein_coding")
  count.table <- merge(count.table, big.annot)
  rownames(count.table) <- count.table[, "gene_id"]
  rm(gtf)
  big.annot2 <- NULL
  if (!dir.exists(pathForDESeq2)) {
    dir.create(pathForDESeq2, recursive = TRUE)
  }
  # I will loop over the list all.analyses:
  # First is the factorToStudy
  for (factorToStudy in names(all.analyses)) {
    print("FACTOR TO STUDY")
    print(factorToStudy)
    # Second indicates the looping variable
    for (i in 1:length(all.analyses[[factorToStudy]])) {
      loopingVariable <- names(all.analyses[[factorToStudy]][i])
      print("LOOPING VARIABLE")
      print(loopingVariable)
      # Third indicates the subsetting
      pre.samples.plan <- samples.plan.df
      subsetting.list <- all.analyses[[factorToStudy]][i][[loopingVariable]]
      subsetting.name <- ""
      for (rn in names(subsetting.list)) {
        print("SUBSET")
        print(rn)
        pre.samples.plan <- pre.samples.plan[pre.samples.plan[, rn] %in% subsetting.list[[rn]], ] 
        subsetting.name <- paste0(subsetting.name, paste(subsetting.list[[rn]], collapse = "or"), "_")
      }
      looping.values <- unique(pre.samples.plan[, loopingVariable])
      ref.value <- NULL
      if (factorToStudy == loopingVariable) {
        ref.value <- levels(pre.samples.plan[, loopingVariable])[1]
        looping.values <- setdiff(looping.values, ref.value)
      }
      for (my.value in looping.values) {
        new.samples.plan <- pre.samples.plan[pre.samples.plan[, loopingVariable] %in% c(ref.value, my.value), ]
        # Drop levels for factorToStudy
        new.samples.plan[, factorToStudy] <- factor(new.samples.plan[, factorToStudy],
                                                    levels = intersect(levels(new.samples.plan[, factorToStudy]),
                                                                       unique(new.samples.plan[, factorToStudy])))
        base.filename <- paste0(factorToStudy, "_", subsetting.name, my.value)
        if (factorToStudy == loopingVariable) {
          base.filename <- paste0(base.filename, "vs", ref.value, "_")
        } else {
          base.filename <- paste0(base.filename, "_")
        }
        print(base.filename)
        # Run or read DESeq2 results with Wald test threshold of FC at 1.5
        if ( !file.exists(file.path(pathForDESeq2, paste0(base.filename, "DESeq2significant.txt")))) {
          print(new.samples.plan)
          deseqAnaWithCovariates(count.table, factorToStudy, NULL,
                                 file.path(pathForDESeq2, base.filename),
                                 new.samples.plan,
                                 LRT = F, lfcT = log2FC.threshold, writeRLOG = F,
                                 gene_id = gene_id)
          # theta = c(0.15, 0.99))
        } else {
          print("Exists")
        }
        # Add results to the dataframe
        all.res <- read.delim(file.path(pathForDESeq2, paste0(base.filename, "DESeq2Results.txt")))
        rownames(all.res) <- all.res[, gene_id]
        # Add results to the dataframe
        big.annot[, paste0(base.filename, "l2fc")] <- all.res$log2FoldChange[match(big.annot[, gene_id], all.res[, gene_id])]
        big.annot[, paste0(base.filename, "padj")] <- all.res$padj[match(big.annot[, gene_id], all.res[, gene_id])]
        big.annot[, paste0(base.filename, "signif")] <-
          with(all.res[match(big.annot[, gene_id], all.res[, gene_id]), ], !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > log2FC.threshold)
        tail(big.annot)
        all.res.fmt <- all.res[, c("gene_id", "baseMean", "log2FoldChange", "padj")]
        all.res.fmt$factor <- factorToStudy
        all.res.fmt$subsetting <- subsetting.name
        all.res.fmt$value <- my.value
        all.res.fmt$ref.value <- ref.value
        big.annot2 <- rbind(big.annot2, all.res.fmt)
        tail(big.annot2)
      }
    }
  }
  write.table(big.annot, file.path(pathForDESeq2, "summary.txt"), sep = "\t", quote = F, row.names = F)
  write.table(big.annot2, file.path(pathForDESeq2, big.table.fn.long), sep = "\t", quote = F, row.names = F)
} else {
  big.annot <- read.delim(file.path(pathForDESeq2, "summary.txt"))
  colnames(big.annot) <- gsub("^X", "", colnames(big.annot))
  colnames(big.annot) <- gsub("\\.", "-", colnames(big.annot))
  big.annot2 <- read.delim(file.path(pathForDESeq2, big.table.fn.long))
}

#### General stats ####
signif <- subset(big.annot2, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > log2FC.threshold)
signif$sign <- "up"
signif$sign[signif$log2FoldChange < 0] <- "down"
signif$value <- factor(signif$value, levels = c(
  sort(unique(as.numeric(as.character(samples.plan.df$NumberOfCells)))),
  setdiff(unique(signif$value), samples.plan.df$NumberOfCells)
))
signif$comparison.name <- paste0(signif$value, "VS", signif$ref.value)
signif$comparison.name <- factor(signif$comparison.name, levels = unique(signif$comparison.name[order(signif$value)]))
signif$factor <- factor(signif$factor, levels = c("NumberOfCells", "Condition"))

signif$valueFormatted <- ifelse(
  signif$comparison.name == "noGSKiVSCTL",
  "-CHI",
  as.character(signif$value)
)
signif$valueFormatted <- factor(signif$valueFormatted, levels = unique(signif$valueFormatted[order(signif$value)]))

my.colors.mixed <- c(my.colors, my.cond.colors)
names(my.colors.mixed) <- signif$valueFormatted[match(names(my.colors.mixed), signif$value)]

g <- ggplot(signif, aes(x = valueFormatted, fill = valueFormatted)) +
  geom_bar() +
  geom_text(
    stat = "count", aes(label = after_stat(count)),
    vjust = -1
  ) +
  facet_grid(. ~ subsetting,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(subsetting = c(
      "300_" = "300 cells",
      "CTL_" = "CTL condition"
    ))
  ) +
  ggtitle("Number of DE genes") +
  xlab(expression(N[0])) +
  ylab("Gene count") +
  scale_y_continuous(limits = c(0, 980), expand = c(0, NA)) +
  theme_classic(base_size = 20) +
  # theme(axis.text.x = element_text(size = 10)) +
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)
  ) +
  scale_fill_manual(expression(N[0]), values = my.colors.mixed)

ggsave(file.path(pathForDESeq2, "DESeq2_nb.pdf"), width = 7, height = 8)

#### Venn diagram ####

venn.list <- list()
venn.list[["Above 300 cells"]] <- paste0("NumberOfCells_CTL_", 600 * 1:3, "vs300_")
names(venn.list[["Above 300 cells"]]) <- paste(600 * 1:3, "vs300cells")
venn.list[["Below 300 cells"]] <- paste0("NumberOfCells_CTL_", 50 * 1:2, "vs300_")
names(venn.list[["Below 300 cells"]]) <- paste(50 * 1:2, "vs300cells")

# Intersection all below 300

test <- "Below 300 cells"
# Split up and down to get euler:
list.for.euler.up <- lapply(venn.list[[test]], function(base.filename) {
  big.annot$gene_short_name[big.annot[, paste0(base.filename, "signif")] & big.annot[, paste0(base.filename, "l2fc")] > 0]
})
list.for.euler.down <- lapply(venn.list[[test]], function(base.filename) {
  big.annot$gene_short_name[big.annot[, paste0(base.filename, "signif")] & big.annot[, paste0(base.filename, "l2fc")] < 0]
})
names(list.for.euler.up) <- paste0(gsub(" vs300cells", "", names(venn.list[[test]])), "_up")
names(list.for.euler.down) <- paste0(gsub(" vs300cells", "", names(venn.list[[test]])), "_down")
list.for.euler <- c(list.for.euler.up, list.for.euler.down)
e <- euler(list.for.euler)
my.colors.euler <- c(
  my.colors[gsub(" vs300cells", "", names(venn.list[[test]]))],
  adjustcolor(
    my.colors[gsub(" vs300cells", "", names(venn.list[[test]]))],
    alpha.f = 0.5
  )
)
pdf(file.path(pathForDESeq2, "below_300_euler.pdf"))
print(plot(e,
  main = test,
  fills = my.colors.euler,
  labels = list(col = "grey70"),
  quantities = list(col = "grey70")
))
dev.off()

test <- "Above 300 cells"
# Split up and down to get euler:
list.for.euler.up <- lapply(venn.list[[test]], function(base.filename) {
  big.annot$gene_short_name[big.annot[, paste0(base.filename, "signif")] & big.annot[, paste0(base.filename, "l2fc")] > 0]
})
list.for.euler.down <- lapply(venn.list[[test]], function(base.filename) {
  big.annot$gene_short_name[big.annot[, paste0(base.filename, "signif")] & big.annot[, paste0(base.filename, "l2fc")] < 0]
})
names(list.for.euler.up) <- paste0(gsub(" vs300cells", "", names(venn.list[[test]])), "_up")
names(list.for.euler.down) <- paste0(gsub(" vs300cells", "", names(venn.list[[test]])), "_down")
list.for.euler <- c(list.for.euler.up, list.for.euler.down)
e <- euler(list.for.euler)
my.colors.euler <- c(
  my.colors[gsub(" vs300cells", "", names(venn.list[[test]]))],
  adjustcolor(
    my.colors[gsub(" vs300cells", "", names(venn.list[[test]]))],
    alpha.f = 0.5
  )
)
pdf(file.path(pathForDESeq2, "above_300_euler.pdf"))
print(plot(e,
  quantities = T,
  main = test,
  fills = my.colors.euler
))
dev.off()

writeLines(capture.output(sessionInfo()), file.path(path, "sessionInfo_DESeq2.txt"))
