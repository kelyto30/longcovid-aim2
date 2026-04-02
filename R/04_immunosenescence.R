# =============================================================================
# Big Aim 2 — Phase 4: Sub-Aim 2.2 — Immunosenescence & T Cell Exhaustion
# Gene panels: KLRG1, CD57, TIGIT, LAG3, PD-1, CD28 loss;
#              TOX, NR4A family, exhaustion signatures
# Dataset: GSE224615 (Long COVID PBMC RNA-seq)
#
# Author: Kelechi Wisdom Elechi
# Institution: Barshop Institute for Longevity and Aging Studies
#              UT Health San Antonio
# Date: 2026
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA); library(DESeq2)
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(dplyr); library(tidyr); library(patchwork)
  library(RColorBrewer); library(viridis); library(ggbeeswarm); library(scales)
})

cat("=== Loading Phase 1-3 objects ===\n")
vst_mat      <- readRDS("data/vst_matrix_gse224615.rds")
col_data     <- readRDS("data/coldata_gse224615.rds")
sen_scores   <- readRDS("data/gsva_senescence_scores.rds")
sen_scores_df <- readRDS("data/gsva_senescence_scores_df.rds")

grp_cols <- c("Healthy"="#2196F3","COVID_Convalescent"="#FF9800",
              "LongCOVID_NoBrainFog"="#9C27B0","LongCOVID_BrainFog"="#F44336")

expressed <- rownames(vst_mat)

# ── 1. DEFINE IMMUNOSENESCENCE GENE SETS ─────────────────────────────────────
cat("\n=== Defining immunosenescence gene sets ===\n")

# Core immunosenescence markers — senescent T cell surface phenotype
immunosen_panel <- c(
  "KLRG1",           # Terminal differentiation marker — senescent T/NK cells
  "B3GAT1",          # CD57 — canonical immunosenescence marker
  "TIGIT",           # Co-inhibitory; elevated in exhausted/senescent T cells
  "LAG3",            # Exhaustion checkpoint; co-expressed with PD-1
  "PDCD1",           # PD-1 — checkpoint receptor; exhaustion/senescence
  "HAVCR2",          # TIM-3 — exhaustion marker
  "ENTPD1",          # CD39 — exhaustion-associated
  "FCGR3A",          # CD16 — adaptive NK cell marker; immunosenescence
  "NCR1",            # NKp46 — NK cell functionality marker
  "SELL",            # CD62L — naive T cell marker (LOSS = senescence)
  "CCR7",            # Lymph node homing (LOSS = terminally differentiated)
  "IL7R",            # CD127 — naive/memory T cell (LOSS = exhaustion)
  "TCF7"             # TCF1 — stem-like T cell; loss = exhaustion
)

# CD28 — critical senescence marker (loss of CD28 = canonical immunosenescence)
cd28_loss_genes <- c("CD28","CD27","CD127")   # Loss genes (downregulation = senescence)
cd28_senescence <- c("KLRG1","B3GAT1","TIGIT","LAG3","PDCD1")  # Gain genes

# T cell exhaustion — TOX/NR4A transcription factor signature
exhaustion_tfs <- c(
  "TOX",   # Master exhaustion regulator
  "TOX2",  # TOX family
  "NR4A1", # Nur77 — drives exhaustion program
  "NR4A2", # Nurr1 — co-regulates exhaustion
  "NR4A3", # NOR1
  "BATF",  # AP-1 family; exhaustion-associated
  "IRF4",  # Exhaustion co-regulator
  "PRDM1", # BLIMP1 — terminal differentiation
  "EOMES", # T-bet/Eomes axis — exhaustion
  "TBX21"  # T-bet — effector/exhaustion balance
)

# Broader T cell exhaustion gene signature (Wherry/Bhatt modules)
exhaustion_signature <- c(
  "PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","CD160","2B4",
  "TOX","TOX2","NR4A1","NR4A2","NR4A3",
  "BATF","IRF4","PRDM1","EOMES","TBX21",
  "ENTPD1","CXCR5","TCF7",
  "GZMB","GZMK","PRF1",                      # Effector/cytotoxicity
  "IFNG","TNF","IL2",                          # Cytokine production (reduced)
  "SELL","CCR7","IL7R","CD28","CD27",          # Naive/memory loss markers
  "KLRG1","B3GAT1","CX3CR1"                   # Terminal differentiation
)

# NK cell immunosenescence
nk_immunosen <- c("FCGR3A","NCR1","KLRG1","B3GAT1","CX3CR1",
                   "FCER1G","TYROBP","KIR2DL1","KIR2DL3","LILRB1")

# Compiled gene set list
immunosen_sets <- list(
  Immunosenescence_Panel = immunosen_panel,
  Exhaustion_TFs         = exhaustion_tfs,
  Exhaustion_Signature   = exhaustion_signature,
  NK_Immunosenescence    = nk_immunosen
)

# Filter to expressed genes
for (gs in names(immunosen_sets)) {
  n_before <- length(immunosen_sets[[gs]])
  immunosen_sets[[gs]] <- intersect(immunosen_sets[[gs]], expressed)
  cat(gs, ":", n_before, "->", length(immunosen_sets[[gs]]), "genes in dataset\n")
}

# ── 2. GSVA IMMUNOSENESCENCE SCORING ─────────────────────────────────────────
cat("\n=== Running GSVA immunosenescence scoring ===\n")

gsva_isen_param <- GSVA::gsvaParam(
  exprData = vst_mat,
  geneSets = immunosen_sets,
  minSize  = 5,
  maxSize  = 500
)

gsva_isen_scores <- GSVA::gsva(gsva_isen_param, verbose = FALSE)
gsva_isen_df     <- as.data.frame(t(gsva_isen_scores))
gsva_isen_df$sample_id <- rownames(gsva_isen_df)
gsva_isen_df <- dplyr::left_join(gsva_isen_df, col_data, by = "sample_id")

cat("Immunosenescence scoring complete.\n")

# ── 3. DIRECT EXPRESSION: CD28 LOSS AND IMMUNOSENESCENCE MARKERS ──────────────
cat("\n=== Direct expression of immunosenescence marker genes ===\n")

all_marker_genes <- unique(c(immunosen_panel, exhaustion_tfs,
                               cd28_loss_genes, "CD28","CD27"))
all_marker_genes <- intersect(all_marker_genes, expressed)
cat("Immunosenescence marker genes in dataset:", length(all_marker_genes), "\n")

# Heatmap of immunosenescence markers
isen_heatmap_mat <- vst_mat[all_marker_genes, , drop = FALSE]

ann_col_df <- data.frame(Group = col_data$group, row.names = col_data$sample_id)
ann_colors  <- list(Group = grp_cols)

pheatmap::pheatmap(
  isen_heatmap_mat,
  scale             = "row",
  annotation_col    = ann_col_df,
  annotation_colors = ann_colors,
  color             = colorRampPalette(c("#2196F3","white","#F44336"))(100),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 8,
  main              = "Immunosenescence Marker Gene Expression — GSE224615",
  filename          = "results/immunosenescence/fig9_immunosenescence_heatmap.pdf",
  width = 10, height = 8
)
cat("Figure 9 (immunosenescence heatmap) saved.\n")

# ── 4. IMMUNOSENESCENCE SCORE BOXPLOTS ───────────────────────────────────────
cat("\n=== Generating immunosenescence score plots ===\n")

plot_isen <- function(df, score_col, title) {
  ggplot(df, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, linewidth = 0.5) +
    ggbeeswarm::geom_beeswarm(size = 2.5, alpha = 0.8, cex = 1.5) +
    scale_fill_manual(values = grp_cols, guide = "none") +
    labs(title = title, x = NULL, y = "GSVA Score") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          plot.title = element_text(face = "bold", size = 9))
}

pi1 <- plot_isen(gsva_isen_df, "Immunosenescence_Panel",
                  "Immunosenescence Panel\n(KLRG1, CD57, TIGIT, LAG3, PD-1...)")
pi2 <- plot_isen(gsva_isen_df, "Exhaustion_TFs",
                  "Exhaustion Transcription Factors\n(TOX, NR4A1/2/3, BATF...)")
pi3 <- plot_isen(gsva_isen_df, "Exhaustion_Signature",
                  "Full Exhaustion Signature\n(Wherry modules)")
pi4 <- plot_isen(gsva_isen_df, "NK_Immunosenescence",
                  "NK Cell Immunosenescence\n(CD16, NKp46...)")

combined_isen <- (pi1 + pi2) / (pi3 + pi4) +
  patchwork::plot_annotation(
    title    = "Sub-Aim 2.2: Immunosenescence & T Cell Exhaustion Scores",
    subtitle = "GSE224615 Long COVID PBMC RNA-seq | GSVA scoring",
    theme    = theme(plot.title = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(colour = "gray40", size = 9))
  )

ggsave("results/immunosenescence/fig10_immunosenescence_scores.pdf",
       combined_isen, width = 12, height = 9, dpi = 300)
ggsave("results/immunosenescence/fig10_immunosenescence_scores.png",
       combined_isen, width = 12, height = 9, dpi = 300)
cat("Figure 10 saved.\n")

# ── 5. STATISTICAL TESTING ────────────────────────────────────────────────────
cat("\n=== Statistical testing: BrainFog vs. NoBrainFog ===\n")

lc_isen <- gsva_isen_df %>%
  dplyr::filter(group %in% c("LongCOVID_NoBrainFog","LongCOVID_BrainFog"))

isen_stats <- lapply(names(immunosen_sets), function(gs) {
  bf  <- lc_isen[[gs]][lc_isen$group == "LongCOVID_BrainFog"]
  nbf <- lc_isen[[gs]][lc_isen$group == "LongCOVID_NoBrainFog"]
  wt  <- wilcox.test(bf, nbf, exact = FALSE)
  data.frame(
    gene_set      = gs,
    n_brainfog    = length(bf),
    n_no_brainfog = length(nbf),
    median_bf     = round(median(bf, na.rm=TRUE), 4),
    median_nbf    = round(median(nbf, na.rm=TRUE), 4),
    delta         = round(median(bf, na.rm=TRUE) - median(nbf, na.rm=TRUE), 4),
    p_value       = wt$p.value
  )
}) %>% dplyr::bind_rows()

isen_stats$padj <- p.adjust(isen_stats$p_value, method = "BH")

cat("Immunosenescence score comparison — Brain Fog vs. No Brain Fog:\n")
print(isen_stats)

# ── 6. CO-ENRICHMENT TEST: IMMUNOSENESCENCE vs. SENESCENCE ───────────────────
cat("\n=== KEY TEST: Immunosenescence vs. Senescence co-enrichment ===\n")

coenrich_df <- data.frame(
  sample_id    = colnames(vst_mat),
  ImmunoSen    = as.numeric(gsva_isen_scores["Immunosenescence_Panel", ]),
  Exhaustion   = as.numeric(gsva_isen_scores["Exhaustion_Signature", ]),
  SenMayo      = as.numeric(sen_scores["SenMayo", ]),
  SASP_Core    = as.numeric(sen_scores["SASP_Core", ])
) %>% dplyr::left_join(col_data, by = "sample_id")

# Correlation: immunosenescence vs. senescence
cor_isen_sen <- cor.test(coenrich_df$ImmunoSen, coenrich_df$SenMayo,
                          method = "pearson")
cor_exh_sasp  <- cor.test(coenrich_df$Exhaustion, coenrich_df$SASP_Core,
                           method = "pearson")

cat("Immunosenescence Panel vs. SenMayo: r =",
    round(cor_isen_sen$estimate, 3), ", p =", round(cor_isen_sen$p.value, 4), "\n")
cat("Exhaustion Signature vs. SASP Core: r =",
    round(cor_exh_sasp$estimate, 3), ", p =", round(cor_exh_sasp$p.value, 4), "\n")

# Scatter plot — the mechanistic test
p_coenrich <- ggplot(coenrich_df, aes(x = ImmunoSen, y = SenMayo, colour = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, colour = "black",
              linetype = "dashed", linewidth = 0.8) +
  ggrepel::geom_text_repel(aes(label = sample_id), size = 2.5,
                             max.overlaps = 15, segment.size = 0.3) +
  scale_colour_manual(values = grp_cols, name = "Group") +
  annotate("text",
           x = min(coenrich_df$ImmunoSen) + 0.05,
           y = max(coenrich_df$SenMayo)   - 0.05,
           label = paste0("r = ", round(cor_isen_sen$estimate, 3),
                           "\np = ", format(cor_isen_sen$p.value, digits = 3)),
           size = 4, fontface = "bold", colour = "black", hjust = 0) +
  labs(
    title    = "Co-Enrichment Test: Immunosenescence vs. Cellular Senescence",
    subtitle = "Mechanistic hypothesis: Immune aging drives senescence accumulation in Long COVID",
    x        = "Immunosenescence Panel Score (GSVA)",
    y        = "SenMayo Senescence Score (GSVA)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "gray40", size = 8))

ggsave("results/immunosenescence/fig11_coenrichment_isen_vs_senescence.pdf",
       p_coenrich, width = 9, height = 6, dpi = 300)
ggsave("results/immunosenescence/fig11_coenrichment_isen_vs_senescence.png",
       p_coenrich, width = 9, height = 6, dpi = 300)
cat("Figure 11 (co-enrichment scatter) saved.\n")

# ── 7. SAVE RESULTS ───────────────────────────────────────────────────────────
saveRDS(gsva_isen_scores, "data/gsva_immunosenescence_scores.rds")
saveRDS(gsva_isen_df,     "data/gsva_immunosenescence_scores_df.rds")
saveRDS(coenrich_df,      "data/coenrichment_scores.rds")
write.csv(isen_stats,     "results/immunosenescence/immunosenescence_statistics.csv", row.names=FALSE)
write.csv(gsva_isen_df,   "results/immunosenescence/immunosenescence_scores_per_sample.csv", row.names=FALSE)
write.csv(coenrich_df,    "results/integration/coenrichment_scores.csv", row.names=FALSE)

cat("\n=================================================================\n")
cat("  PHASE 4 (Sub-Aim 2.2) COMPLETE — IMMUNOSENESCENCE SUMMARY\n")
cat("=================================================================\n")
cat("Gene sets scored: Immunosenescence Panel, Exhaustion TFs,\n")
cat("                  Full Exhaustion Signature, NK Immunosenescence\n")
print(isen_stats[, c("gene_set","delta","p_value","padj")])
cat("\nCo-enrichment test (Immunosenescence vs. SenMayo):\n")
cat("  r =", round(cor_isen_sen$estimate, 3),
    ", p =", round(cor_isen_sen$p.value, 4), "\n")
cat("\nNext: Run R/05_coenrichment.R\n")
cat("=================================================================\n")
