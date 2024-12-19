# We define the color palette used for featurePlots 
pal <- c("#FEB24C",	"#FD8D3C",	"#FC4E2A",	"#E31A1C",	"#BD0026",	"#800026")
Col.featurePlot <- alpha(c("#D3D3D3", pal),0.9)

# And the color palette for UMAP:

# The time:
my.time.colors <- grey(level = c(0.8,0.4))
names(my.time.colors) <- c("120h", "144h")
# The treatment
my.treatment.colors <- c("#184279", "#30807f", "#82a954", "#b9b35f", "#d4b093")
names(my.treatment.colors) <- c( "100", "300", "600", "1800", "5400")
# The batch
my.batch.colors <- c("#4BABE6", "#CC9933")

# And the clusters:
# Group new cluster names by categories so they will have colors along a colormap:
# list.Fate.level is a list with the categories:
# Order matter! It will be used in all plots!

list.Fate.level <- list("Neuronal" = c("NMPs","Neural tube","Neuron progen.", "Neuron prec." ),
                        "Mesoderm" = c("Post. PSM", "Ant. PSM", "Somitic mes.", "Dermomyotome", "Sclerotome", "Cardiac mes."),
                        "Endoderm" = c("Endoderm", "Visceral endoderm"),
                        "Other" = c("Endothelium", "Axial mes.", "Pluripotent"))
# list.color is a list with the same categories and colors that describe the heatmap. There must be at least 2 colors:
list.color  <- list("Neuronal" = c('#DAE3F3', '#002060'),
                    "Mesoderm" = c('#FBE5D6','#5F2C09'),
                    "Endoderm" = c('#5FE756', '#70AD47'),
                    "Other" = c('#FBBEDE', '#7030A0'))
# Check all fates have a category name
if ("" %in% names(list.Fate.level)) {
  stop("Some fates level has no name!\n")
}
# Check all category names in list.Fate.level are also in list.color
if (any(! names(list.Fate.level) %in% names(list.color))) {
  stop("The following Fates are not in colors: ", paste(names(list.Fate.level)[!names(list.Fate.level) %in% names(list.color)], collapse = ", "), "\n")
}
# Create the colormaps
my.fate.colors <- unlist(lapply(names(list.Fate.level), function(fate) {
  colorRampPalette(list.color[[fate]])(length(list.Fate.level[[fate]]))
}))
# Give them as names the new cluster names
names(my.fate.colors) <- unlist(list.Fate.level)