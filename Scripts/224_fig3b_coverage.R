suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
theme_set(theme_classic(base_size = 8) + theme(axis.text = element_text(size = 5.6, colour = "black"),
  axis.title = element_text(size = 7), legend.text = element_text(size = 5.6), legend.title = element_text(size = 5.6),
  plot.margin = margin(2,3,2,2)))
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
  data.table(co = cc, trait = t, r = if (nrow(rs) > 0) tanh(mean(atanh(pmax(pmin(rs$r,0.99),-0.99)))) else NA_real_, nfam = nrow(rs))
}))))
res[, trait := factor(lab2[trait], levels = lab2[tr])]
res[, co := factor(co, levels = rev(cot), labels = c("Whole","Centromere","Middle","Distal"))]
res[, lab := ifelse(abs(r) >= 0.2, sprintf("%.2f", r), "")]
print(res[, .(trait, co, r = round(r,2), nfam)])
pb <- ggplot(res, aes(trait, co, fill = r)) + geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.1) +
  scale_fill_gradient2(low = "#2E6FA3", mid = "white", high = "#C0392B", limits = c(-0.25,0.25), name = "r") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 5))
ggsave(file.path(OUT, "panels", "Fig3b_correlation.png"), pb, width = 88, height = 50, units = "mm", dpi = 600)
cat("b done\n")
