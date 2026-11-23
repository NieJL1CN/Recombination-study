suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 5.8, colour = "black"),
  axis.title = element_text(size = 7), legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.6),
  plot.margin = margin(2,3,2,2)))
pub <- fread(file.path(TAB, "Table39_CO_QTL_trait_scan.csv"))
lab <- c(AGBMNitr_M_pct="Above-ground N concentration", BMNUtE_Calc_kgkgN="Biomass N utilisation efficiency",
  GPD_Calc_pct="Grain protein deviation", BMNUE_Calc_kgkgN="Biomass N use efficiency", HI_Calc_pct="Harvest index",
  STRNUP_Calc_kgha="Stem N uptake", StrawYld_Calc_tha="Straw yield", BM_Calc_gm2="Biomass",
  NHI_Calc_pct="N harvest index", PH_M_cm="Plant height", Ant_dto_day="Anthesis date",
  BMNUP_Calc_kgha="Biomass N uptake", GFP_calc_day="Grain filling period", GNitr_NIRS_pct="Grain N concentration",
  GNUP_Calc_kgha="Grain N uptake", GY_Calc_tha="Grain yield", Mat_dto_day="Maturity date",
  NDVI_M_Idx="NDVI", NUE_Calc_kgkgN="N use efficiency", NupE_Calc_kgkgN="N uptake efficiency",
  NutE_Calc_kgkgN="N utilisation efficiency", STRNitr_NIRS_pct="Straw N concentration",
  GFR_Calc_gday="Grain filling rate", Hd_dto_day="Heading date")
P <- data.table(trait = lab[pub$trait], n = pub$observed, source = "Published")
mq <- fread(file.path(TAB, "Table2f_metaQTL.csv"))
co <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
O <- rbindlist(lapply(c("GW","GY"), function(tt) {
  cl <- mq[grepl(tt, traits)]
  obs <- sum(sapply(seq_len(nrow(co)), function(i) any(co$chr[i] == cl$chr & co$pos[i] >= cl$start - 20 & co$pos[i] <= cl$end + 20)))
  data.table(trait = ifelse(tt == "GW", "Thousand-grain weight", "Grain yield"), n = obs, source = "This study")
}))
A <- rbind(P, O)
A[, trait := factor(trait, levels = unique(c(A[source=="Published"][order(n), trait], "Thousand-grain weight")))]
pa <- ggplot(A, aes(n, trait, colour = source)) + geom_point(size = 1.4) +
  scale_colour_manual(values = setNames(c("#2E6FA3", "#C0392B"), c("Published", "This study")), name = NULL) +
  labs(x = "Colocalised CO QTL (count)", y = NULL) + theme(legend.position = "top")
ggsave(file.path(OUT, "panels", "Fig3a_integrated.png"), pa, width = 100, height = 78, units = "mm", dpi = 600)
cat("a done\n")
