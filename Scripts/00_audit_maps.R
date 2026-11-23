## 00_audit_maps.R -- audit the consensus maps
suppressPackageStartupMessages(library(data.table))
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
L <- readRDS(file.path(ROOT, "Source/Derived/maps_min1.rds"))
CM <- L$consensus; M <- L$markers[!is.na(g_rank)]
pos <- unique(M[, .(marker, ppos)])
D <- merge(CM, pos, by = "marker")
A <- D[, .(n = .N,
           rho = cor(ppos, cM, method = "spearman"),
           ties = round(mean(diff(sort(cM)) == 0), 3),
           len_cM = round(max(cM), 1),
           len_Mb = round(max(ppos) / 1e6, 0),
           cM_per_Mb = round(max(cM) / (max(ppos) / 1e6), 2)),
        by = .(panel, chr)]
fwrite(A, file.path(ROOT, "Tables", "Table1c_map_audit.csv"))
cat("chromosomes:", nrow(A), "\n")
cat("rho < 0.99:", sum(A$rho < 0.99), " min rho", round(min(A$rho), 4), "\n")
cat("median ties:", round(median(A$ties), 3), " max", max(A$ties), "\n")
cat("median cM/Mb:", round(median(A$cM_per_Mb), 2), "\n")
print(A[order(cM_per_Mb)][1:4])
print(A[order(-cM_per_Mb)][1:4])
