## ============================================================================
## 04_phenotypes.R -- tidy public phenotype data (JIC field trials)
## ============================================================================
suppressPackageStartupMessages({library(readxl); library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
P <- file.path(ROOT, "Source/Phenotype"); TAB <- file.path(ROOT, "Tables")
DER <- file.path(ROOT, "Source/Derived")

f1 <- file.path(P, "DFW_ParWat_JIC_phenotype_data_grassroot_format.xlsx")
sh <- grep("grassroot_data$", suppressMessages(excel_sheets(f1)), value = TRUE)
D <- rbindlist(lapply(sh, function(s) {
  d <- as.data.table(suppressMessages(read_excel(f1, sheet = s)))
  d[, year := as.integer(sub(".*H(20[0-9]{2}).*", "\\1", s))]
  d
}), fill = TRUE)
D <- D[!is.na(AccessionName) & AccessionName != "empty"]
D[, family := sub("_.*$", "", AccessionName)]
D[, line := AccessionName]
tr <- c("PH_M_cm", "GY_Calc_tha", "GW_M_g1000grn", "GN_Calc_Grnm2", "GWid_M_mm", "GLng_M_mm")
tr <- intersect(tr, names(D))
D <- D[, c("year", "family", "line", tr), with = FALSE]
L <- D[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = .(family, line), .SDcols = tr]
L[, n_year := uniqueN(D$year), by = line]
L <- L[!is.na(PH_M_cm) | !is.na(GY_Calc_tha)]
fwrite(L, file.path(TAB, "Table2a_phenotypes_by_line.csv"))
S <- L[, lapply(.SD, function(x) sum(!is.na(x))), by = family, .SDcols = tr]
S <- merge(S, L[, .(n_lines = uniqueN(line)), by = family], by = "family")
S <- S[order(-n_lines)]
fwrite(S, file.path(TAB, "Table2b_phenotype_coverage.csv"))
cat("lines:", nrow(L), "families:", uniqueN(L$family), "\n")
print(head(S[order(-n_lines)], 12))
print(S[, .(families = .N, med_lines = median(n_lines))])