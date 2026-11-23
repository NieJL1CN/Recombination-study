## 94_meiosis_expression.R -- CS meiosis-stage expression of RECQL1 homeologs
suppressPackageStartupMessages({library(data.table); library(jsonlite); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); FIG <- file.path(ROOT, "Figures")
out <- list()
for (p in c("PRJEB5029_tbl","PRJNA755254_tpm_tbl")) {
  f <- file.path(TAB, paste0("expr_", p, ".json"))
  if (!file.exists(f)) next
  j <- fromJSON(f, simplifyVector = FALSE)
  for (g in j$data$results) {
    for (pt in g$points) {
      out[[length(out) + 1]] <- data.table(project = p, gene = g$gene_id,
        stage = pt$label, value = as.numeric(pt$value), std = as.numeric(pt$std))
    }
  }
}
E <- rbindlist(out)
E[, copy := sub("^TraesCS([0-9][ABD]).*$", "\\1", gene)]
fwrite(E, file.path(TAB, "Table87_RECQL1_meiosis_expression.csv"))
cat("stages:", paste(unique(E[project == "PRJEB5029_tbl"]$stage), collapse = " | "), "\n\n")
print(dcast(E[project == "PRJEB5029_tbl"], gene ~ stage, value.var = "value"))
cat("\nper-project stages:\n"); print(E[, .N, by = .(project)])
E[, stage_f := factor(stage, levels = c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI"))]
p <- ggplot(E[project == "PRJEB5029_tbl"], aes(stage_f, value, colour = copy, group = copy)) +
  geom_line(linewidth = 0.5) + geom_point(size = 1) +
  geom_errorbar(aes(ymin = value - std, ymax = value + std), width = 0.15, linewidth = 0.3) +
  scale_colour_manual(values = c("4A" = "#2C6E9B", "4B" = "#E08A2E", "4D" = "#2E8B57"), name = "subgenome") +
  theme_classic(base_size = 7) + labs(x = NULL, y = "log2 expression (PRJEB5029)") +
  scale_x_discrete(labels = c("latent_lepto" = "leptotene", "zygo_pachy" = "zygo/pachy",
                              "diplo_dia" = "diplo/diak", "metaphaseI" = "meta I")) +
  theme(axis.text.x = element_text(size = 6))
ggsave(file.path(FIG, "FigS_RECQL1_meiosis_expression.png"), p, width = 3.6, height = 2.4, dpi = 600)
cat("figure saved\n")