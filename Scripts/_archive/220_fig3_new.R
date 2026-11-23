suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6, colour = "black"),
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
P <- data.table(trait = lab[pub$trait], fold = pub$fold, source = "Published")
mq <- fread(file.path(TAB, "Table2f_metaQTL.csv"))
co <- fread(file.path(TAB, "Table46_CO_QTL_allele_direction.csv"))[significant == TRUE]
ours <- rbindlist(lapply(c("GW","GY"), function(tt) {
  cl <- mq[grepl(tt, traits)]
  span <- sum(cl$end - cl$start)
  obs <- sum(sapply(seq_len(nrow(co)), function(i) any(co$chr[i] == cl$chr & co$pos[i] >= cl$start - 20 & co$pos[i] <= cl$end + 20)))
  exp <- nrow(co) * span / 4112
  data.table(trait = lab[tt], fold = ifelse(exp > 0, obs/exp, NA_real_), source = "This study")
}))
A <- rbind(P, ours)
A[, trait := factor(trait, levels = A[source=="Published"][order(fold), trait])]
pa <- ggplot(A, aes(fold, trait, colour = source)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_point(size = 1.3) +
  scale_colour_manual(values = setNames(c("#2E6FA3", "#C0392B"), c("Published", "This study")), name = NULL) +
  labs(x = "Fold enrichment of CO QTL", y = NULL) + theme(legend.position = "top")
ggsave(file.path(OUT, "panels", "Fig3a_integrated.png"), pa, width = 95, height = 72, units = "mm", dpi = 600)
co2 <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
co2[, distal := R1 + R3][, middle := R2a + R2b][, centromere := C][, whole := R1 + R2a + C + R2b + R3]
ph1 <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
ph2 <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))[, .(line, Hd_dto_day, PH_M_cm2 = PH_M_cm)]
ph1 <- merge(ph1, ph2[, .(line, PH_M_cm2)], by = "line", all.x = TRUE)
ph1[, PH_M_cm := fifelse(is.na(PH_M_cm), PH_M_cm2, PH_M_cm)]
m <- merge(co2, ph1, by.x = "id", by.y = "line", all.x = TRUE)
m <- merge(m, ph2, by.x = "id", by.y = "line", all.x = TRUE)
tr <- c("Hd_dto_day","PH_M_cm","GY_Calc_tha","GW_M_g1000grn","GN_Calc_Grnm2","GWid_M_mm","GLng_M_mm")
lab2 <- c(Hd_dto_day="Heading date", PH_M_cm="Plant height", GY_Calc_tha="Grain yield",
          GW_M_g1000grn="Grain weight", GN_Calc_Grnm2="Grains per m2", GWid_M_mm="Grain width", GLng_M_mm="Grain length")
cot <- c("distal","middle","centromere","whole")
res <- rbindlist(lapply(cot, function(cc) rbindlist(lapply(tr, function(t) {
  mm <- m[is.finite(get(cc)) & is.finite(get(t))]
  rs <- mm[, { x <- get(cc); y <- get(t); if (.N >= 10) .(r = cor(x, y)) else .(r = NA_real_) }, by = dataset]
  rs <- rs[is.finite(r)]
  data.table(co = cc, trait = t, r = if (nrow(rs) > 0) tanh(mean(atanh(pmax(pmin(rs$r,0.99),-0.99)))) else NA_real_)
}))))
res[, trait := factor(lab2[trait], levels = lab2[tr])]
res[, co := factor(co, levels = rev(cot), labels = c("Whole","Centromere","Middle","Distal"))]
res[, lab := ifelse(abs(r) >= 0.2, sprintf("%.2f", r), "")]
pb <- ggplot(res, aes(trait, co, fill = r)) + geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.1) +
  scale_fill_gradient2(low = "#2E6FA3", mid = "white", high = "#C0392B", limits = c(-0.25,0.25), name = "r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5))
ggsave(file.path(OUT, "panels", "Fig3b_correlation.png"), pb, width = 85, height = 50, units = "mm", dpi = 600)
cat("Fig3 panels done\n")
