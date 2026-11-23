suppressPackageStartupMessages({library(data.table); library(ggplot2)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
W <- fread(file.path(ROOT,"Tables/Table158_watkins_co_by_line.csv"))[, .(group="landrace", pop=family, co)][order(pop)]
C <- fread(file.path(ROOT,"Tables/Table166_cimmyt_3k_byline.csv"))[, .(group, pop, co)]
D <- rbind(W, C)
fwrite(D, file.path(ROOT,"Tables/Table167_co_groups_final.csv"))
P <- D[, .(mean=mean(co)), by=.(group,pop)]
cat("means:", P[, .(n=.N, mean=round(mean(mean),1), sd=round(sd(mean),1)), by=group]$group, "\n")
print(P[, .(n=.N, mean=round(mean(mean),1), sd=round(sd(mean),1)), by=group])
cat("t:", signif(t.test(mean~group, data=P)$p.value,3), " var:", signif(var.test(mean~group, data=P)$p.value,3), "\n")
p <- ggplot(D, aes(x=group, y=co, fill=group)) +
  geom_violin(trim=FALSE, alpha=0.65, colour=NA) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white") +
  theme_classic(base_size=12) +
  labs(x=NULL, y="Crossovers per RIL") +
  scale_fill_manual(values=c(landrace="#4C72B0", cultivar="#DD8452")) +
  theme(legend.position="none")
ggsave(file.path(ROOT,"Figures/FigS18_CO_landrace_vs_cultivar.png"), p, width=4.2, height=4.2, dpi=600)
cat("saved\n")
