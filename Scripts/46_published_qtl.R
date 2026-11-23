## 46_published_qtl.R -- published Paragon x Watkins QTL (Cheng et al. 2024) for our cohort
suppressPackageStartupMessages({library(readxl); library(data.table)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
WD <- file.path(ROOT, "Source/WebData"); TAB <- file.path(ROOT, "Tables")
COH <- fread(file.path(TAB, "Table0_cohort.csv"))[keep == TRUE]
fams <- sub("_WATDE.*$", "", COH$dataset)
d <- as.data.table(suppressWarnings(read_excel(file.path(WD, "41586_2024_7682_MOESM5_ESM.xlsx"), sheet = "S17", skip = 1)))
setnames(d, c("QTL","oldQTL","population","chr","peak_marker","start_marker","end_marker","CI_markers",
              "trait","trait_full","trait_mean","add_effect","increase_allele","LOD","perc_var",
              "pos_cM","CIstart","CIend","chr_refseq","bp_start","bp_end","bp_peak","parent",
              "ancestral_group","environment","NIL_MAS","NIL_included","NIL_overlap"))
d <- d[population %in% fams]
cat("QTL for cohort families:", nrow(d), " traits:", uniqueN(d$trait), "\n")
d[, chr2 := sub("_.*$", "", chr)]
fwrite(d, file.path(TAB, "Table35_published_QTL_cohort.csv"))
print(d[, .N, by = trait][order(-N)][1:12])
cat("\nheading-date related rows:\n")
print(d[trait %in% c("Ant_dto_day","Hd_dto_day")][1:5,
      .(population, chr2, trait, LOD, pos_cM, bp_peak, environment)])
