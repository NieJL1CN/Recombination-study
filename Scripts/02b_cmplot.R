## 02b_cmplot.R -- CMplot marker-density figures (green scale, Mb axis)
suppressPackageStartupMessages({library(data.table); library(CMplot)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); FIG <- file.path(ROOT, "Figures")
L <- readRDS(file.path(DER, "maps.rds")); M <- L$markers; CM <- L$consensus
pos <- unique(M[, .(marker, ppos)])
setwd(FIG)
for (P in unique(CM$panel)) {
  D <- merge(CM[panel == P, .(chr, marker, cM)], pos, by = "marker")
  if (nrow(D) < 100) next
  Pmap <- data.frame(SNP = D$marker,
                     Chromosome = factor(D$chr, levels = as.vector(t(outer(1:7, c("A","B","D"), paste0)))),
                     Position = D$ppos)
  CMplot(Pmap, plot.type = "d", bin.size = 5e6,
         col = c("#F4F9F4", "#2E8B57"),
         file = "png", dpi = 400, file.output = TRUE, verbose = FALSE, main = "",
         file.name = paste0("Fig1a_CMplot_density_", P, "."))
}
f <- list.files(FIG, pattern = "^Marker_Density\\.")
for (x in f) {
  n <- sub("^Marker_Density\\.", "", x)
  n <- sub("\\.\\.png$", ".png", n)
  file.rename(file.path(FIG, x), file.path(FIG, n))
}
cat("cmplot done\n")