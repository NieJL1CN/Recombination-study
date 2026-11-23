suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_pub <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 6, colour = "black"), axis.title = element_text(size = 7),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        strip.text = element_text(size = 6), strip.background = element_blank(),
        plot.tag = element_text(size = 9, face = "bold"),
        panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.25),
        plot.margin = margin(2,3,2,2))
theme_set(theme_pub)
lab <- c(AGBMNitr_M_pct="Above-ground N concentration", BMNUtE_Calc_kgkgN="Biomass N utilisation efficiency",
  GPD_Calc_pct="Grain protein deviation", BMNUE_Calc_kgkgN="Biomass N use efficiency", HI_Calc_pct="Harvest index",
  STRNUP_Calc_kgha="Stem N uptake", StrawYld_Calc_tha="Straw yield", BM_Calc_gm2="Biomass",
  NHI_Calc_pct="N harvest index", PH_M_cm="Plant height", Ant_dto_day="Anthesis date",
  BMNUP_Calc_kgha="Biomass N uptake", GFP_calc_day="Grain filling period", GNitr_NIRS_pct="Grain N concentration",
  GNUP_Calc_kgha="Grain N uptake", GY_Calc_tha="Grain yield", Mat_dto_day="Maturity date",
  NDVI_M_Idx="NDVI", NUE_Calc_kgkgN="N use efficiency", NupE_Calc_kgkgN="N uptake efficiency",
  NutE_Calc_kgkgN="N utilisation efficiency", STRNitr_NIRS_pct="Straw N concentration",
  GFR_Calc_gday="Grain filling rate", Hd_dto_day="Heading date")
pub <- fread(file.path(TAB, "Table39_CO_QTL_trait_scan.csv"))
pub[, trait2 := lab[trait]][, trait2 := factor(trait2, levels = pub[order(fold), trait2])]
pub[, zero := fold == 0]
p3a <- ggplot(pub, aes(fold, trait2)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_segment(aes(x = 0, xend = fold, y = trait2, yend = trait2), linewidth = 0.3, colour = "grey70") +
  geom_point(aes(size = n_pub_qtl, shape = zero), colour = "#2E6FA3") +
  scale_shape_manual(values = setNames(c(16,1), c("FALSE","TRUE")), guide = "none") +
  scale_size_continuous(range = c(0.6, 2.4), name = "Published QTL") +
  scale_x_continuous(limits = c(0, 2.0), breaks = c(0, 0.5, 1, 1.5, 2)) +
  labs(x = "Fold enrichment near published trait QTL", y = NULL)
P <- fread(file.path(TAB, "Table168_HD_profiles.csv"))
P[, panel := paste0(family, " (", chr, ")")]
P[, panel := factor(panel, levels = unique(panel))]
ann <- data.table(panel = levels(P$panel), lab = c("peak 3.08 to 3.11","peak 3.01 to 2.96","peak 3.45 to 3.38"))
p3b <- ggplot(P, aes(pos, lod, colour = model)) + geom_line(linewidth = 0.5) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = 1.05, vjust = 1.6, size = 2.0) +
  scale_colour_manual(values = c("no covariate" = "#1E9E9E", "Hd covariate" = "#E0726A"), name = NULL) +
  labs(x = "Chromosome position (cM)", y = "LOD") +
  theme(legend.position = "top",
        plot.margin = margin(1, 2, 1, 0),
        axis.title.y = element_text(size = 7, margin = margin(r = 1)))
co <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
co[, distal := R1 + R3][, middle := R2a + R2b][, centromere := C][, whole := R1 + R2a + C + R2b + R3]
p43 <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
p2a <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
m <- merge(co, p43[, .(line, Hd_dto_day, PH_M_cm, GY_Calc_tha, GW_M_g1000grn)], by.x = "id", by.y = "line", all.x = TRUE)
m <- merge(m, p2a[, .(line, GY2 = GY_Calc_tha, GWid_M_mm, GLng_M_mm, GN_Calc_Grnm2)], by.x = "id", by.y = "line", all.x = TRUE)
m[, GY_Calc_tha := fifelse(is.na(GY_Calc_tha), GY2, GY_Calc_tha)]
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
p3c <- ggplot(res, aes(trait, co, fill = r)) + geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.1) +
  scale_fill_gradient2(low = "#2E6FA3", mid = "white", high = "#C0392B", limits = c(-1, 1), name = "r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5))
fig3 <- p3a / (p3b | p3c) + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.4, 0.8), widths = c(1.65, 0.7))
ggsave(file.path(OUT, "Fig3.tiff"), fig3, width = 183, height = 185, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig3, width = 183, height = 185, units = "mm", dpi = 600)
cat("Fig3 done\n")
