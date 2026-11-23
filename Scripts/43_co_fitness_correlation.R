## 43_co_fitness_correlation.R -- individual CO vs multiseason fitness traits within families
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
CO <- fread(file.path(TAB, "Table22c_CO_rate_individual.csv"))
CO[, family := sub("_WATDE.*$", "", dataset)]
CO[, co := R1 + R2a + C + R2b + R3]
PH <- fread(file.path(TAB, "Table43_RIL_phenotypes_multiseason.csv"))
TRAITS <- c(Hd_dto_day = "heading date", PH_M_cm = "plant height",
            GY_Calc_tha = "grain yield", GW_M_g1000grn = "grain weight")
out <- list()
for (fm in unique(CO$family)) {
  d <- CO[family == fm]
  d <- merge(d[, .(line = id, co)], PH, by = "line")
  if (!nrow(d)) next
  for (tr in names(TRAITS)) {
    x <- d$co; y <- d[[tr]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 40) next
    r <- cor.test(x[ok], y[ok], method = "spearman")
    out[[length(out) + 1]] <- data.table(family = fm, trait = tr, n = sum(ok),
                                         rho = round(unname(r$estimate), 3), p = r$p.value)
  }
}
R <- rbindlist(out)
fwrite(R, file.path(TAB, "Table32_CO_fitness_correlation.csv"))
cat("tests:", nrow(R), "\n")
print(R[, .(n = .N, sig05 = sum(p < 0.05), med_rho = round(median(rho), 3)), by = trait])
fisher <- R[, { stat <- -2 * sum(log(p)); df <- 2 * .N
  list(n_families = .N, chi2 = round(stat, 2), p_meta = signif(pchisq(stat, df, lower.tail = FALSE), 3),
       mean_rho = round(mean(rho), 3)) }, by = trait]
print(fisher)
fwrite(fisher, file.path(TAB, "Table32b_CO_fitness_fisher.csv"))
R[, trait_label := TRAITS[trait]]
p <- ggplot(R, aes(x = reorder(family, rho), y = rho, colour = trait_label, shape = trait_label)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
  geom_point(size = 0.9) +
  scale_colour_manual(values = c("heading date" = "#6A4C93", "plant height" = "#2C6E9B",
                                 "grain yield" = "#2E8B57", "grain weight" = "#E08A2E"), name = NULL) +
  scale_shape_manual(values = c(16, 15, 17, 18), name = NULL) +
  coord_flip() + theme_classic(base_size = 7) +
  theme(legend.position = c(0.8, 0.15)) +
  labs(title = "Within-family phenotypic correlation", x = NULL,
       y = "Spearman rho (CO vs trait within family)") +
  theme(plot.title = element_text(size = 8, hjust = 0))
ggsave(file.path(FIG, "Fig4c_CO_fitness_correlation.png"), p, width = 4.4, height = 4.2, dpi = 600)
cat("figure done\n")