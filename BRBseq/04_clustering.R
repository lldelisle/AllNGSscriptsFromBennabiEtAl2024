# Load dependencies:
if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
suppressPackageStartupMessages(safelyLoadAPackageInCRANorBioconductor("pheatmap"))
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("biomaRt")
safelyLoadAPackageInCRANorBioconductor("goseq")
safelyLoadAPackageInCRANorBioconductor("TxDb.Mmusculus.UCSC.mm10.ensGene")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("viridis")

source("BRBseq/config.R")
source("BRBseq/add.flag.function.R")

nmodules <- 6

# Prepare inputs
samples.plan.df <- read.delim(samples.plan.file, check.names = FALSE)
rownames(samples.plan.df) <- samples.plan.df$sample
norm.exp.df <- read.delim(norm.counts.file, check.names = FALSE, row.names = 1)
big.annot <- read.delim(file.path(pathForDESeq2, "summary.txt"))
genes.to.highlight <- read.delim(file.path("BRBseq", "inputs", "genes_to_highlight.txt"), header = FALSE)$V1

my.samples.order <- samples.plan.df[order(samples.plan.df$Condition, as.numeric(as.character(samples.plan.df$NumberOfCells)), samples.plan.df$Rep), "sample"]
samples.plan.df$NumberOfCells <- factor(samples.plan.df$NumberOfCells, levels = c(300, unique(setdiff(unique(samples.plan.df$NumberOfCells), 300))))
annot <- subset(samples.plan.df, select = c("NumberOfCells"))
rownames(annot) <- samples.plan.df$sample

# Only keep colors in samples.plan.df
my.colors <- my.colors[intersect(names(my.colors), unique(samples.plan.df$NumberOfCells))]

# 50vs300 or 1800vs300
my.genes <- big.annot$gene_short_name[big.annot$NumberOfCells_CTL_50vs300_signif | big.annot$NumberOfCells_CTL_1800vs300_signif]
my.genes.id <- big.annot$gene_id[match(my.genes, big.annot$gene_short_name)]
names(my.genes) <- my.genes.id

scaled.df <- t(apply(norm.exp.df[my.genes.id, grep("noGSKi", my.samples.order, invert = TRUE, value = TRUE)], 1, scale))
colnames(scaled.df) <- grep("noGSKi", my.samples.order, invert = TRUE, value = TRUE)

avg.300 <- rowMeans(scaled.df[, grep(300, colnames(scaled.df))])

new.scaled.df <- scaled.df - matrix(rep(avg.300, ncol(scaled.df)), ncol = ncol(scaled.df))

annotation_colors <- list(NumberOfCells = my.colors)

d <- as.dist(1 - cor(t(new.scaled.df)))
hc <- hclust(d, method = "ward.D2")

annot.genes.module <- data.frame(gene_id = my.genes.id)
rownames(annot.genes.module) <- my.genes.id

annot.genes.module[, paste0("k", nmodules)] <-
    LETTERS[1:nmodules][cutree(hc, nmodules)][match(my.genes.id, hc$labels)]
annotation_colors[[paste0("k", nmodules)]] <- brewer.pal(nmodules, "Set3")
names(annotation_colors[[paste0("k", nmodules)]]) <- LETTERS[1:nmodules]
annot.genes.module$gene_id <- NULL

# From https://www.biostars.org/p/400381/
genes.formatted <- lapply(
    my.genes[hc$labels],
    function(x) bquote(italic(.(x)))
)

pdf(file.path(pathForDESeq2, "signif_extremes_pheatmap.pdf"))
pheatmap(
    new.scaled.df,
    cluster_cols = FALSE,
    cluster_rows = hc,
    annotation_col = annot,
    annotation_row = annot.genes.module,
    annotation_colors = annotation_colors,
    show_colnames = FALSE,
    labels_row = as.expression(genes.formatted),
    color = colorRampPalette(c(
        "#053061",
        "#6bacd0",
        "#f7f7f7",
        "#e58268",
        "#67001f"
    ))(100),
    breaks = c(min(new.scaled.df), seq(-2.5, 2.5, length.out = 98), max(new.scaled.df))
)
dev.off()

new.scaled.df2 <- new.scaled.df
rownames(new.scaled.df2) <- my.genes[rownames(new.scaled.df)]
annot.genes.module2 <- subset(annot.genes.module, select = "k6")
rownames(annot.genes.module2) <- my.genes[rownames(annot.genes.module)]

pdf(file.path(pathForDESeq2, "signif_extremes_pheatmap_highlight.pdf"))
heat <- pheatmap(
    new.scaled.df2,
    cluster_cols = FALSE,
    cluster_rows = hc,
    annotation_col = annot,
    annotation_row = annot.genes.module2,
    annotation_colors = annotation_colors,
    show_colnames = FALSE,
    color = colorRampPalette(c(
        "#053061",
        "#6bacd0",
        "#f7f7f7",
        "#e58268",
        "#67001f"
    ))(100),
    breaks = c(min(new.scaled.df), seq(-2.5, 2.5, length.out = 98), max(new.scaled.df)),
    silent = TRUE
)
add.flag(heat, kept.labels = genes.to.highlight, repel.degree = 0)
dev.off()
annot.genes.module$gene_short_name <- my.genes[rownames(annot.genes.module)]

# Output the table in the same order as the heatmap:
annot.genes.module <- annot.genes.module[hc$labels[hc$order], ]

write.table(
    annot.genes.module,
    file.path(pathForDESeq2, "signif_extremes_clustering.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

# To compute GO
# We use the same version as gtf (102):
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://nov2020.archive.ensembl.org")
mart <- useDataset("mmusculus_gene_ensembl", mart)
genes <- getBM(
    attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "go_id", "namespace_1003"),
    mart = mart
)

tested.genes.id <- big.annot$gene_id[!is.na(big.annot$NumberOfCells_CTL_50vs300_l2fc) | !is.na(big.annot$NumberOfCells_CTL_1800vs300_l2fc)]
tested.genes <- subset(genes, ensembl_gene_id %in% tested.genes.id)
genes.length <- getlength(tested.genes.id, "mm10", "ensGene")
gene2cat.bp <- with(
    subset(tested.genes, namespace_1003 == "biological_process"),
    split(go_id, ensembl_gene_id)
)
for (my.module in LETTERS[1:nmodules]) {
    output.file <- file.path(pathForDESeq2, paste0("signif_extremes_clustering_k", nmodules, "_module", my.module, "_BP.txt"))
    if (!file.exists(output.file)) {
        print(paste("module", my.module, "in", nmodules))
        selected.genes <- rownames(annot.genes.module[annot.genes.module[, paste0("k", nmodules)] == my.module, ])
        print(paste("Found", length(selected.genes), "genes"))
        gene.vector <- as.integer(tested.genes.id %in% selected.genes)
        names(gene.vector) <- tested.genes.id
        pwf <- nullp(gene.vector, bias.data = genes.length)
        GO.BP <- goseq(pwf, gene2cat = gene2cat.bp)
        write.table(GO.BP, output.file, sep = "\t", quote = F, row.names = F)
        dev.off()
    }
}
go.files <- file.path(pathForDESeq2, paste0("signif_extremes_clustering_k", nmodules, "_module", LETTERS[1:nmodules], "_BP.txt"))
names(go.files) <- LETTERS[1:nmodules]

go.res <- do.call(rbind, lapply(names(go.files), function(ngf) {
    if (!file.exists(go.files[ngf])) {
        return(NULL)
    }
    df <- read.delim(go.files[ngf])
    df$module <- ngf
    return(df)
}))

go.res$term[go.res$category == "GO:0060394"] <- "negative regulation of pathway-restricted SMAD protein phosphorylation"
go.res$term[go.res$category == "GO:0016576"] <- "histone dephosphorylation"

selected.go <- unique(unlist(lapply(
    lapply(
        split(go.res, go.res$module),
        head,
        n = 5
    ),
    "[[", "category"
)))

selected.df <- subset(go.res, category %in% selected.go)

selected.df$term[is.na(selected.df$term)] <- AnnotationDbi::Term(GO.db::GOOBSOLETE)[selected.df$category[is.na(selected.df$term)]]

# min.pval.per.go <- sapply(
#     lapply(
#         split(selected.df, selected.df$term),
#         "[[",
#         "over_represented_pvalue"
#     ),
#     min
# )
# selected.df$term <- factor(selected.df$term, levels = rev(names(sort(min.pval.per.go))))
selected.df$term <- factor(selected.df$term, levels = rev(selected.df$term[match(selected.go, selected.df$category)]))

g <- ggplot(
    selected.df,
    aes(x = module, y = term)
) +
    geom_point(aes(size = -log10(over_represented_pvalue), color = numDEInCat)) +
    scale_color_gradientn(colours = viridis(100)) +
    theme_classic()

ggsave(
    file.path(pathForDESeq2, paste0("signif_extremes_clustering_k", nmodules, "_summary_GO.pdf")),
    g,
    width = 10, height = 7
)

writeLines(capture.output(sessionInfo()), "BRBseq/sessionInfo_clustering.txt")
