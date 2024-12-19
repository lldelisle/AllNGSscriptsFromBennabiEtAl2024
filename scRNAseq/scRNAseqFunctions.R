# Install required packages
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
if (!"usefulLDfunctions" %in% installed.packages()) {
  devtools::install_github("lldelisle/usefulLDfunctions", upgrade = "never")
}
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("dplyr")
safelyLoadAPackageInCRANorBioconductor("Seurat")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("velocyto.R")
safelyLoadAPackageInCRANorBioconductor("pheatmap")
if (!"SeuratWrappers" %in% installed.packages()) {
  remotes::install_github('satijalab/seurat-wrappers', ref = 'b6d519b69e8a6364f7275d116aa8b0c1ee783364', upgrade = "never")
}
safelyLoadAPackageInCRANorBioconductor("SeuratWrappers")

# now we define our custom made functions

savePngPdf <- function(plot, directory, fileName, width = 6, height = 4){
  # This functions will save a ggPlot object as a pdf and a png
  ggsave(filename = paste0(directory, "/", fileName, ".pdf"),
         plot = plot, width = width, height = height)
  ggsave(filename = paste0(directory, "/", fileName, ".png"),
         plot = plot, width = width, height = height)
}

featurePlotMayran <- function(SampleName, gene.to.display, params,
                              width = 6, height = 4, pt.size = 0.3,
                              order = T, split.by = NULL, by.rows = T, ...) {
  # This functions displays featurePlot with consistent settings, color palette and save them to Png and Pdf
  seurat.object <- eval(parse(text=SampleName))
  g <- FeaturePlot(seurat.object, features = gene.to.display,
                   order = order, pt.size = pt.size, ...) +
    NoAxes() +
    scale_color_gradientn(gene.to.display, colours = params$Col.featurePlot) +
    theme(panel.background = element_blank(), strip.background =  element_blank(), title = element_text(face = "bold.italic")) +
    ggtitle("")
  # This part allows you to use the split.by argument of featurePlot while retaining the same aesthetics in both
  if(!is.null(split.by)){
    if(order){
      expression.values <- FetchData(seurat.object, vars = gene.to.display)
      group.values <- seurat.object[[]][order(!is.na(x = expression.values[, gene.to.display]),
                                             expression.values[, gene.to.display]),
                                       split.by]
      if (by.rows){
        g <- g +
          facet_grid(rows = vars(as.factor(group.values)))
      } else {
        g <- g +
          facet_grid(cols = vars(as.factor(group.values)))
      }

      
    } else{
      group.values <- seurat.object[[]][, split.by]
      if (by.rows){
        g <- g +
          facet_grid(rows = vars(as.factor(group.values)))
      } else {
        g <- g +
          facet_grid(cols = vars(as.factor(group.values)))
      }
      
    }
  }
  savePngPdf(g, params$current.fig,
             paste0("featurePlot.",gene.to.display,".",params$nameRDS,".",SampleName),
             width = width, height = height)
  return(g)
}


DimPlotMayran <- function(SampleName, params,
                          group.by = NULL,
                          split.by = NULL,
                          # CellColours = NULL,
                          width = 6, height = 4, label =T, ...) {
  # This functions displays DimPlot with consistent settings, color palette and save them to Png and Pdf
  seurat.object <- eval(parse(text=SampleName))
  g <- DimPlot(seurat.object,
               label = label,
               pt.size = 0.3, 
               label.size = 4,
               group.by = group.by, split.by = split.by,
               # cols = alpha(my.fate.colors, 0.3),
               repel = T, shuffle = T, raster = F, ...) +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
    ) + NoAxes()
  savePngPdf(g, params$current.fig,
             paste0("DimPlot.group.",group.by,".split", split.by, params$nameRDS,".",SampleName),
             width = width, height = height)
  return(g)
}


SwitchFigure <- function(Fig = "1", Panel = "a") {
  current.fig <- paste0(output.directory, "Fig.", Fig,"/", Panel)
  dir.create(current.fig, showWarnings = FALSE, recursive = TRUE)
  return(current.fig)
}


add.level <- function(seurat.object, my.gene, thresholds = c(-0.1, 0, 1, 2, 8)) {
  # This functions allows you to take a seurat object and define metaData according
  # to chosen threshold, 4 levels are defined here: ND, low, medium, high)
  stopifnot(length(thresholds) == 5)
  my.labels <- paste0(my.gene, "_", c("ND", "low", "medium", "high"))
  return(AddMetaData(seurat.object, cut(FetchData(object = seurat.object, vars = paste0('rna_', my.gene))[, 1], breaks = thresholds, labels = my.labels), col.name = paste0(my.gene, "_level")))
}
