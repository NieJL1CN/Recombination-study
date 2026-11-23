## 52_ril_heading_date.R -- multiseason RIL phenotypes from JIC grassroot data
suppressPackageStartupMessages({library(readxl); library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
fp <- file.path(WD, "DFW_ParWat_JIC_phenotype_data_grassroot_format.xlsx")
sh <- excel_sheets(fp); sh <- sh[grepl("grassroot_data", sh)]
out <- list()
for (s in sh) {
  d <- as.data.table(suppressWarnings(read_excel(fp, sheet = s)))
  if (!"AccessionName" %in% names(d)) next
  keep <- intersect(c("Hd_dto_day","PH_M_cm","GY_Calc_tha","GW_M_g1000grn"), names(d))
  if (!length(keep)) next
  dd <- d[!is.na(AccessionName), c("AccessionName", keep), with = FALSE]
  setnames(dd, "AccessionName", "line")
  for (k in keep) dd[[k]] <- suppressWarnings(as.numeric(dd[[k]]))
  dd <- dd[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = line, .SDcols = keep]
  dd[, season := s]
  out[[s]] <- dd
  cat(s, ":", nrow(dd), "lines, traits:", paste(keep, collapse = ","), "\n")
}
H <- rbindlist(out, fill = TRUE)
traits <- intersect(c("Hd_dto_day","PH_M_cm","GY_Calc_tha","GW_M_g1000grn"), names(H))
A <- H[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = line, .SDcols = traits]
A[, n_season := uniqueN(H$season[!is.na(match(H$line, line))]), by = line]
fwrite(H, file.path(TAB, "Table42_RIL_phenotypes_by_season.csv"))
fwrite(A, file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
cat("\nlines:", nrow(A), " traits:", paste(traits, collapse = ", "), "\n")
print(A[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols = traits])
