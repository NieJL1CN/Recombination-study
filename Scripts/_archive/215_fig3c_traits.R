suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(png); library(grid)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 6.2, colour = "black"),
  axis.title = element_text(size = 7.2), legend.text = element_text(size = 6), legend.title = element_text(size = 6),
  plot.margin = margin(2,3,2,2)))
co <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
co[, distal := R1 + R3][, whole := R1 + R2a + C + R2b + R3]
ph1 <- fread(file.path(TAB, "Table2a_phenotypes_by_line.csv"))
ph2 <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))[, .(line, Hd_dto_day)]
m <- merge(co, ph1, by.x = "id", by.y = "line", all.x = TRUE)
m <- merge(m, ph2, by.x = "id", by.y = "line", all.x = TRUE)
tr <- c("Hd_dto_day","PH_M_cm","GY_Calc_tha","GW_M_g1000grn","GN_Calc_Grnm2","GWid_M_mm","GLng_M_mm")
lab <- c(Hd_dto_day="Heading date", PH_M_cm="Plant height", GY_Calc_tha="Grain yield",
         GW_M_g1000grn="Grain weight", GN_Calc_Grnm2="Grains per m2", GWid_M_mm="Grain width", GLng_M_mm="Grain length")
res <- rbindlist(lapply(c("whole","distal"), function(cot) rbindlist(lapply(tr, function(t) {
  mm <- m[is.finite(get(cot)) & is.finite(get(t))]
  rs <- mm[, { x <- get(cot); y <- get(t); if (.N >= 10) .(r = cor(x, y)) else .(r = NA_real_) }, by = dataset]
  rs <- rs[is.finite(r)]
  data.table(co = cot, trait = t, r = tanh(mean(atanh(pmax(pmin(rs$r,0.99),-0.99)))), nfam = nrow(rs))
}))))
res[, trait := factor(lab[trait], levels = lab[tr])]
res[, co := factor(co, levels = c("distal","whole"), labels = c("Distal","Whole"))]
res[, lab := ifelse(abs(r) >= 0.2, sprintf("%.2f", r), "")]
pc <- ggplot(res, aes(trait, co, fill = r)) + geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab), size = 2.4) +
  scale_fill_gradient2(low = "#2E6FA3", mid = "white", high = "#C0392B", limits = c(-0.25,0.25), name = "mean r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5))
ggsave(file.path(OUT, "panels", "Fig3c_within_family.png"), pc, width = 90, height = 45, units = "mm", dpi = 600)
print(res[, .(trait, co, r = round(r,3), nfam)])
