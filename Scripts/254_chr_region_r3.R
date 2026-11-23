## 254_chr_region_r3.R -- chromosome-level CO by region and R3-vs-length data
## Produces Tables/Table183_chr_region_CO.csv used by Fig. 1d.
## Same data path as 31_regression.R (FigS30): marker-based region boundaries,
## events from Table18c, mean per chromosome.
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
DER <- file.path(ROOT, "Source/Derived"); TAB <- file.path(ROOT, "Tables")
CHR <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))
B <- data.table(chr = CHR,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84)*1e6,
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200)*1e6,
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373)*1e6,
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552)*1e6)
M <- readRDS(file.path(DER, "maps_min1.rds"))$markers[!is.na(g_rank)]
POS <- unique(M[, .(marker, chr, ppos)])
LEN <- merge(POS, B, by = "chr")
LEN[, region := fifelse(ppos < b1, "R1", fifelse(ppos < b2, "R2a",
      fifelse(ppos < b3, "C", fifelse(ppos < b4, "R2b", "R3"))))]
CHL <- LEN[, .(len = (max(ppos) - min(ppos)) / 1e6), by = chr]
RGL <- LEN[, .(len = (max(ppos) - min(ppos)) / 1e6), by = .(chr, region)]
EV <- fread(file.path(TAB, "Table18c_CO_events.csv"))
EV[, family := sub("_WATDE.*$", "", dataset)]
NIND <- EV[, .(n_ind = uniqueN(id)), by = dataset]
COC <- EV[, .(k = .N), by = .(dataset, chr)]
COC <- merge(COC, NIND, by = "dataset")[, co := k / n_ind][, family := sub("_WATDE.*$", "", dataset)]
COC <- COC[, .(co = mean(co)), by = .(family, chr)]
COR <- EV[, .(k = .N), by = .(dataset, chr, region)]
COR <- merge(COR, NIND, by = "dataset")[, co := k / n_ind][, family := sub("_WATDE.*$", "", dataset)]
COR <- COR[, .(co = mean(co)), by = .(family, chr, region)]
COC_c <- COC[, .(co = mean(co)), by = chr]
COR_c <- COR[, .(co = mean(co)), by = .(chr, region)]
w <- dcast(COR_c, chr ~ region, value.var = "co", fun.aggregate = mean)
w <- merge(w, COC_c, by = "chr", suffixes = c("", ".total"))
setnames(w, "co", "co_total")
w <- merge(w, CHL, by = "chr")
w <- merge(w, RGL[region == "R3", .(chr, R3_len = len)], by = "chr")
fwrite(w, file.path(TAB, "Table183_chr_region_CO.csv"))
cat("Table183 written: chr =", nrow(w), " r(R3 vs R3_len) =", round(cor(w$R3_len, w$R3), 3), "\n")
