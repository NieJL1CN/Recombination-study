suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.5, colour = "black"),
  axis.title = element_text(size = 8), plot.tag = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5), plot.margin = margin(2,3,2,2)))
tr <- fread(file.path(TAB, "Table39_CO_QTL_trait_scan.csv"))
lab <- c(AGBMNitr_M_pct="Above-ground N concentration", BMNUtE_Calc_kgkgN="Biomass N utilisation efficiency",
  GPD_Calc_pct="Grain protein deviation", BMNUE_Calc_kgkgN="Biomass N use efficiency", HI_Calc_pct="Harvest index",
  STRNUP_Calc_kgha="Stem N uptake", StrawYld_Calc_tha="Straw yield", BM_Calc_gm2="Biomass",
  NHI_Calc_pct="N harvest index", PH_M_cm="Plant height", Ant_dto_day="Anthesis date",
  BMNUP_Calc_kgha="Biomass N uptake", GFP_calc_day="Grain filling period", GNitr_NIRS_pct="Grain N concentration",
  GNUP_Calc_kgha="Grain N uptake", GY_Calc_tha="Grain yield", Mat_dto_day="Maturity date",
  NDVI_M_Idx="NDVI", NUE_Calc_kgkgN="N use efficiency", NupE_Calc_kgkgN="N uptake efficiency",
  NutE_Calc_kgkgN="N utilisation efficiency", STRNitr_NIRS_pct="Straw N concentration",
  GFR_Calc_gday="Grain filling rate", Hd_dto_day="Heading date")
tr[, trait2 := lab[trait]]
tr[, trait2 := factor(trait2, levels = tr[order(fold), trait2])]
tr[, zero := fold == 0]
p3a <- ggplot(tr, aes(fold, trait2)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_segment(aes(x = 0, xend = fold, y = trait2, yend = trait2), linewidth = 0.3, colour = "grey70") +
  geom_point(aes(size = n_pub_qtl, shape = zero), colour = "#2E6FA3") +
  scale_shape_manual(values = setNames(c(16, 1), c("FALSE", "TRUE")), guide = "none") +
  scale_size_continuous(range = c(0.6, 2.4), name = "Published QTL") +
  scale_x_continuous(limits = c(0, 2.0), breaks = c(0, 0.5, 1, 1.5, 2)) +
  labs(x = "Fold enrichment of CO-QTL near published trait QTL", y = NULL) +
  annotate("text", x = 1.02, y = 1.2, label = "expected", size = 2.0, hjust = 0)
p3b <- { img <- readPNG(file.path(OUT, "panels", "Fig3b_profiles.png")); g <- rasterGrob(img, interpolate = TRUE)
  ggplot() + annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_void() }
fig3 <- p3a / p3b + plot_annotation(tag_levels = "a") + plot_layout(heights = c(1.5, 1))
ggsave(file.path(OUT, "Fig3.tiff"), fig3, width = 183, height = 160, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig3.png"), fig3, width = 183, height = 160, units = "mm", dpi = 600)
cat("done\n")
