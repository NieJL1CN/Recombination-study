## 250_unified_main_figures.R -- rebuild Fig1-3 with one shared theme
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel)})
ROOT <- "D:/CAAS/Recombination/WatkinsRILsRecombinationStudy"
TAB <- file.path(ROOT, "Tables"); OUT <- file.path(ROOT, "Figures/main")
BASE <- 8
theme_pub <- theme_classic(base_size = BASE) +
  theme(axis.text = element_text(size = 6, colour = "black"),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        strip.text = element_text(size = 6), strip.background = element_blank(),
        plot.tag = element_text(size = 9, face = "bold"),
        panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.25),
        plot.margin = margin(2,3,2,2))
theme_set(theme_pub)
CHR <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
LEN <- c("1A"=593,"1B"=689,"1D"=495,"2A"=780,"2B"=796,"2D"=649,"3A"=749,"3B"=830,"3D"=614,
         "4A"=743,"4B"=673,"4D"=488,"5A"=708,"5B"=713,"5D"=565,"6A"=617,"6B"=721,"6D"=473,
         "7A"=736,"7B"=750,"7D"=636)
CHR_ot <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))   # boundary vectors are in this order
BD <- data.table(chr = CHR_ot,
  b1 = c(59,62,29,42,59,37,62,66,49,41,42,10,39,52,46,46,56,44,89,12,84),
  b2 = c(151,172,98,206,248,192,249,257,167,180,186,135,140,140,128,216,221,164,239,146,200),
  b3 = c(231,277,171,379,433,338,414,407,287,414,360,288,260,221,207,409,429,280,416,418,373),
  b4 = c(480,534,385,662,660,520,670,728,543,594,537,427,427,430,345,556,651,410,659,660,552))

## ---------------- Fig 1 ----------------
CHRo <- as.vector(t(outer(1:7, c("A","B","D"), paste0)))   # 1A,1B,1D,...,7D (capsule y order)
A <- fread(file.path(TAB, "Table179_CO_10Mb_marker.csv")); A[, binMb := bin/1e6]
A <- merge(A, BD, by = "chr")
A[, region := fifelse(binMb < b1, "R1", fifelse(binMb < b2, "R2a", fifelse(binMb < b3, "C", fifelse(binMb < b4, "R2b", "R3"))))]
A[, region := factor(region, levels = c("R1","R2a","C","R2b","R3"))]
A[, y := match(chr, CHRo)]
A[, h := co / max(co) * 0.8]
p1a <- ggplot(A) +
  geom_rect(aes(xmin = binMb, xmax = binMb + 10, ymin = y - 0.4, ymax = y - 0.4 + h, fill = region),
            colour = "grey35", linewidth = 0.03) +
  scale_fill_manual(values = c(R1 = "#2C6E9B", R2a = "#6A4C93", C = "#BDBDBD", R2b = "#2E8B57", R3 = "#E08A2E"), name = NULL) +
  scale_y_continuous(breaks = 1:21, labels = CHRo, expand = expansion(add = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        axis.ticks = element_line(linewidth = 0.2),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        legend.position = "bottom", legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6)) +
  labs(x = "Physical position (Mb)", y = NULL)
ld <- fread(file.path(TAB, "Table80b_LD_curves_binned.csv"))[group %in% c("carrier","non")][order(bin)]
p1b <- ggplot(ld, aes(bin, r2, colour = group)) + geom_line(linewidth = 0.5) + geom_point(size = 0.3) +
  geom_hline(yintercept = 0.6, linetype = 2, linewidth = 0.3) + geom_vline(xintercept = c(13,15), linetype = 3, linewidth = 0.3) +
  scale_colour_manual(values = c(carrier = "#4C72B0", non = "#DD8452"), labels = c("high-CO","low-CO")) +
  labs(x = "Physical distance (Mb)", y = expression(r^2), colour = NULL) + theme(legend.position = c(0.72,0.8))
regm <- fread(file.path(TAB, "Table18b_CO_by_region_family.csv"))[, .(co = mean(co)), by = region][order(-co)]
regm[, region := factor(region, levels = regm$region)]
p1c <- ggplot(regm, aes(region, co)) + geom_col(fill = "#4C72B0", width = 0.65) +
  geom_text(aes(label = round(co,1)), vjust = -0.3, size = 2.2) + labs(x = NULL, y = "Crossovers per RIL")
Fd <- fread(file.path(TAB, "Table183_chr_region_CO.csv"))
p1d <- ggplot(Fd, aes(R3_len, R3)) + geom_point(size = 1, colour = "#E08214") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "black") +
  labs(x = "R3 length (Mb)", y = "R3 crossovers per RIL") +
  annotate("text", x = 60, y = max(Fd$R3)*0.98, label = sprintf("r = %.2f", cor(Fd$R3_len, Fd$R3)), size = 2.4, hjust = 0)
viol <- fread(file.path(TAB, "Table167_co_groups_final.csv"))[group == "landrace"]
ord <- viol[, .(m = mean(co)), by = pop][order(-m), pop]; viol[, pop := factor(pop, levels = ord)]
p1e <- ggplot(viol, aes(pop, co)) + geom_violin(fill = "#4C72B0", alpha = 0.6, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = NULL, y = "Crossovers per RIL") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5))
fig1 <- p1a + p1b + p1c + p1d + p1e + plot_layout(design = "abb\nacd\naee", widths = c(0.72,1,1)) + plot_annotation(tag_levels = "a")
ggsave(file.path(OUT, "Fig1.tiff"), fig1, width = 183, height = 225, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig1.png"), fig1, width = 183, height = 225, units = "mm", dpi = 600)
cat("Fig1 done\n")

## ---------------- Fig 2 ----------------
chrlev <- c(paste0(1:7,"A"), paste0(1:7,"B"), paste0(1:7,"D"))
len <- data.table(chr = c("1A","1B","1D","2A","2B","2D","3A","3B","3D","4A","4B","4D","5A","5B","5D","6A","6B","6D","7A","7B","7D"),
                  Mb = c(594,689,495,780,801,651,717,830,615,744,673,509,709,713,566,618,720,473,736,750,638))
q <- data.table(qtl = c("ParW579-3B","ParW729-1B","ParW729-4A","ParW352-7D","ParW579-4B","ParW049-7B","ParW007-6A","ParW145-5D","RobW145-6B"),
                chr = c("3B","1B","4A","7D","4B","7B","6A","5D","6B"),
                Mb = c(582.0,330.0,530.9,393.6,167.3,3.7,568.5,544.6,459.8),
                trait = c("distal","distal","distal","distal","distal","centromere","whole","whole","whole"))
len[, chr := factor(chr, levels = chrlev)][, idx := as.numeric(chr)]
q[, chr := factor(chr, levels = chrlev)][, idx := as.numeric(chr)]
len[, cen := Mb * 0.5]
p2a <- ggplot() +
  geom_rect(data = len, aes(xmin = idx - 0.14, xmax = idx + 0.14, ymin = 0, ymax = Mb), fill = "grey85", colour = "grey35", linewidth = 0.25) +
  geom_rect(data = len, aes(xmin = idx - 0.14, xmax = idx + 0.14, ymin = cen - 4, ymax = cen + 4), fill = "white", colour = NA) +
  geom_point(data = q, aes(idx, Mb, colour = trait), size = 1.1) +
  geom_text_repel(data = q, aes(idx, Mb, label = qtl), size = 2.0, direction = "y", nudge_x = 0.45,
                  segment.size = 0.15, min.segment.length = 0) +
  scale_colour_manual(values = c(distal = "#2E6FA3", centromere = "grey30", whole = "#E08214")) +
  scale_x_continuous(breaks = seq_along(chrlev), labels = chrlev, expand = expansion(add = 0.6)) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = NULL, y = "Physical position (Mb)", colour = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
eff <- fread(file.path(TAB, "Table115_QTL_effect_sizes.csv")); eff[, qtl := factor(qtl, levels = eff[order(pct), qtl])]
p2b <- ggplot(eff, aes(pct, qtl, fill = direction)) + geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.2, size = 2.1) +
  scale_fill_manual(values = c(Watkins = "#4C72B0", Paragon = "#DD8452")) +
  labs(x = "Effect (% of family mean)", y = NULL, fill = NULL) + xlim(0, max(eff$pct)*1.18) +
  theme(legend.position = c(0.8, 0.3), legend.key.size = unit(2.5,"mm"))
ind <- fread(file.path(TAB, "Table96_all_parw_induction.csv"))[!is.na(induction)]
ind[, hit := induction >= 1]
p2c <- ggplot(ind, aes(induction, fill = hit)) + geom_histogram(binwidth = 0.25, colour = NA) +
  scale_fill_manual(values = setNames(c("grey70","#C0392B"), c("FALSE","TRUE")), guide = "none") +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.3) +
  annotate("text", x = 1.05, y = Inf, label = paste0(sum(ind$hit), " of ", nrow(ind), " genes >= 2-fold"), vjust = 1.5, hjust = 0, size = 2.2) +
  labs(x = expression(log[2]~"meiosis / vegetative"), y = "Genes")
ex <- fread(file.path(TAB, "Table171_gene_stage_expr.csv"))
ex <- ex[stage %in% c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI")]
ex[, stage := factor(stage, levels = c("latent_lepto","zygo_pachy","diplo_dia","metaphaseI"),
                     labels = c("Leptotene","Zygotene-Pachytene","Diplotene-Diakinesis","Metaphase I"))]
keep <- c("RBR2","CPL3","RECQL1","RHC1A"); ex <- ex[gene %in% keep]; ex[, gene := factor(gene, levels = keep)]
p2d <- ggplot(ex, aes(stage, value, fill = gene)) + geom_col(position = "dodge", width = 0.75) +
  scale_fill_manual(values = c(RBR2 = "#C0392B", CPL3 = "#E08214", RECQL1 = "#2E6FA3", RHC1A = "#2F8F5B")) +
  labs(x = NULL, y = expression(log[2]~"TPM")) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 25, hjust = 1, size = 5.5), legend.key.size = unit(2.5,"mm"))
ft <- fread(file.path(TAB, "Table169_gene_features.csv")); vr <- fread(file.path(TAB, "Table170_gene_variants.csv"))
ord5 <- c("RBR2","CPL3","RECQL1","RHC1A","ATG12"); ypos <- setNames(rev(seq_along(ord5)), rev(ord5))
ft[, y := ypos[gene]]; vr[, y := ypos[gene]]
body <- ft[feature != "promoter", .(lo = min(x1), hi = max(x2), y = y[1]), by = gene]
p2e <- ggplot() +
  geom_segment(data = ft[feature == "promoter"], aes(x = x1, xend = x2, y = y, yend = y), colour = "grey55", linewidth = 0.4) +
  geom_segment(data = body, aes(x = lo, xend = hi, y = y, yend = y), colour = "grey35", linewidth = 0.4) +
  geom_rect(data = ft[feature %in% c("five_prime_UTR","three_prime_UTR")], aes(xmin = x1, xmax = x2, ymin = y-0.07, ymax = y+0.07), fill = "grey55", colour = NA) +
  geom_rect(data = ft[feature == "exon"], aes(xmin = x1, xmax = x2, ymin = y-0.17, ymax = y+0.17), fill = "#2C3E50", colour = NA) +
  geom_point(data = vr, aes(x = x, y = y), colour = "#C0392B", size = 0.4, alpha = 0.6, position = position_jitter(height = 0.07, seed = 1)) +
  scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1), labels = c("0","50","100")) +
  scale_y_continuous(breaks = ypos, labels = names(ypos), limits = c(0.4, 5.6)) +
  labs(x = "Position along locus (%)", y = NULL) + theme(axis.text.y = element_text(size = 5.5))
fig2 <- p2a / p2b / (p2c | p2d) / p2e + plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(1.0, 0.5, 0.75, 0.5))
ggsave(file.path(OUT, "Fig2.tiff"), fig2, width = 183, height = 235, units = "mm", dpi = 600, compression = "lzw")
ggsave(file.path(OUT, "Fig2.png"), fig2, width = 183, height = 235, units = "mm", dpi = 600)
cat("Fig2 done\n")
