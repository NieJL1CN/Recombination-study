## 28_pd.R -- parental polymorphism (PD) per family x chromosome / region
suppressPackageStartupMessages({library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE, dataset]
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
M <- M[dataset %in% COH]
pos <- unique(M[, .(marker, ppos)])
D <- merge(unique(M[, .(dataset, marker, chr)]), pos, by = "marker")
D <- merge(D, B, by = "chr")
D[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
      fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
D[, family := sub("_WATDE.*$", "", dataset)]
CH <- D[, .(pd = .N), by = .(family, chr)]
RG <- D[, .(pd = .N), by = .(family, region)]
fwrite(CH, file.path(TAB, "Table19_PD_by_chr.csv"))
fwrite(RG, file.path(TAB, "Table19b_PD_by_region.csv"))
cat("families", uniqueN(CH$family), " chr rows", nrow(CH), " region rows", nrow(RG), "\n")
print(RG[, .(mean_pd = round(mean(pd))), by = region][order(-mean_pd)])
