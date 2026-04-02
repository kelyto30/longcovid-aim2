# =============================================================================
# Big Aim 2 — Phase 3: Sub-Aim 2.1 — Senescence Transcriptomic Scoring
# Methods: SenMayo ssGSEA + GSVA, CellAge scoring, p21/p16/p53 axis
# Dataset: GSE224615 (Long COVID PBMC RNA-seq)
#
# Author: Kelechi Wisdom Elechi
# Institution: Barshop Institute for Longevity and Aging Studies
#              UT Health San Antonio
# Date: 2026
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA); library(GSEABase); library(DESeq2)
  library(ggplot2); library(ggrepel); library(pheatmap)
  library(dplyr); library(tidyr); library(patchwork)
  library(RColorBrewer); library(viridis); library(scales)
  library(ggbeeswarm)
})

cat("=== Loading Phase 1 & 2 objects ===\n")
vst_mat     <- readRDS("data/vst_matrix_gse224615.rds")
col_data    <- readRDS("data/coldata_gse224615.rds")
deseq_res   <- readRDS("data/deseq2_results_brainfog.rds")

grp_cols <- c("Healthy"="#2196F3","COVID_Convalescent"="#FF9800",
              "LongCOVID_NoBrainFog"="#9C27B0","LongCOVID_BrainFog"="#F44336")

# ── 1. DEFINE SENESCENCE GENE SETS ────────────────────────────────────────────
cat("\n=== Defining senescence gene sets ===\n")

# SenMayo — 105-gene validated senescence signature (Saul et al. 2022, Nature Aging)
# Core senescence genes encompassing SASP, cell cycle arrest, and senescence regulators
senmayo_genes <- c(
  # SASP — secretory phenotype effectors
  "IL6","IL8","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8",
  "CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL13","CCL20",
  "MMP1","MMP3","MMP9","MMP10","MMP12","MMP13","MMP14",
  "IL1A","IL1B","IL7","IL13","IL15","IL33",
  "TNFA","CSF2","VEGFA","HGF","IGF1",
  "SERPINE1","SERPINE2","PAI1",
  "ICAM1","VCAM1",
  "GDF15","IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7",
  # Cell cycle arrest
  "CDKN1A","CDKN2A","CDKN2B","CDKN2D","TP53","RB1",
  "CDKN1B","CDKN1C",
  # DNA damage response
  "H2AX","TP53BP1","BRCA1","ATM","ATR","CHEK1","CHEK2",
  "MDM2","MDM4",
  # Senescence regulators & markers
  "GLB1","LMNB1","HMGA1","HMGA2",
  "BMI1","EZH2",
  "BCL2","BCL2L1","MCL1",
  "SIRT1","SIRT3","SIRT6",
  "MTOR","AKT1","PIK3CA",
  "CGAS","STING1","IRF3",
  "NF2","CDKN2C",
  "WNT2","WNT3A","WNT5A","WNT6",
  "NOTCH1","JAG1","JAG2",
  "TGM2","TGFB1","TGFB2",
  "BMP2","BMP6",
  "PTGS2","PTGES"
)

# Remove any duplicates
senmayo_genes <- unique(senmayo_genes)
cat("SenMayo gene set size:", length(senmayo_genes), "\n")

# CellAge core senescence genes (database-curated)
cellage_genes <- c(
  "TP53","RB1","CDKN1A","CDKN2A","CDKN2B","CDKN2C","CDKN2D",
  "CDKN1B","CDKN1C","E2F1","E2F2","E2F3",
  "SIRT1","SIRT2","SIRT3","SIRT6","SIRT7",
  "LMNB1","LMNB2","LMNA",
  "HMGB1","HMGB2","HMGA1","HMGA2",
  "GLB1","GJB2","GJB6",
  "SERPINE1","IGFBP5","IGFBP7","GDF15",
  "IL6","IL8","CXCL1","MMP3","MMP9",
  "ATM","ATR","CHEK1","CHEK2","H2AFX",
  "BCL2","BCL2L1","BCL2L2","BAX","BBC3",
  "MDM2","MDM4","PPM1D",
  "MTOR","PIK3CB","AKT1","TSC1","TSC2",
  "TERC","TERT","POT1","RAD51"
)
cellage_genes <- unique(cellage_genes)
cat("CellAge gene set size:", length(cellage_genes), "\n")

# SASP core — focused inflammatory secretome
sasp_core_genes <- c("IL6","CXCL8","IL1A","IL1B","MMP3","MMP9",
                      "CXCL1","CCL2","VEGFA","SERPINE1","ICAM1",
                      "TNFA","CSF2","GDF15","IGFBP7","IGFBP3")
sasp_core_genes <- unique(sasp_core_genes)

# p21/p16/p53 axis — canonical cell cycle arrest genes
p21_p16_p53 <- c("CDKN1A","CDKN2A","TP53","RB1","E2F1","MDM2","CDKN1B")

# Senescence gene set list for GSVA
senescence_gene_sets <- list(
  SenMayo   = senmayo_genes,
  CellAge   = cellage_genes,
  SASP_Core = sasp_core_genes,
  p21_p16_p53 = p21_p16_p53
)

# ── 2. FILTER TO EXPRESSED GENES ─────────────────────────────────────────────
cat("\n=== Filtering gene sets to expressed genes ===\n")

expressed_genes <- rownames(vst_mat)

for (gs in names(senescence_gene_sets)) {
  n_before <- length(senescence_gene_sets[[gs]])
  senescence_gene_sets[[gs]] <- intersect(senescence_gene_sets[[gs]], expressed_genes)
  n_after  <- length(senescence_gene_sets[[gs]])
  cat(gs, ":", n_before, "->", n_after, "genes in dataset\n")
}

# ── 3. GSVA SCORING (SAMPLE-LEVEL) ───────────────────────────────────────────
cat("\n=== Running GSVA senescence scoring ===\n")

# GSVA — method='ssgsea' gives per-sample enrichment scores
gsva_param <- GSVA::gsvaParam(
  exprData   = vst_mat,
  geneSets   = senescence_gene_sets,
  minSize    = 5,
  maxSize    = 500
)

gsva_scores <- GSVA::gsva(gsva_param, verbose = FALSE)
gsva_df     <- as.data.frame(t(gsva_scores))
gsva_df$sample_id <- rownames(gsva_df)
gsva_df <- dplyr::left_join(gsva_df, col_data, by = "sample_id")

cat("GSVA scoring complete.\n")
cat("Score matrix dimensions:", nrow(gsva_scores), "gene sets x", ncol(gsva_scores), "samples\n")

# ── 4. SENESCENCE SCORE BOXPLOTS ─────────────────────────────────────────────
cat("\n=== Generating senescence score plots ===\n")

# Function to plot senescence scores by group
plot_scores <- function(df, score_col, title, cols = grp_cols) {
  ggplot(df, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, linewidth = 0.5) +
    ggbeeswarm::geom_beeswarm(size = 2.5, alpha = 0.8, cex = 1.5) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(title = title,
         x = NULL, y = "GSVA Enrichment Score") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 9),
      plot.title   = element_text(face = "bold", size = 10)
    )
}

p_sen1 <- plot_scores(gsva_df, "SenMayo",
  "SenMayo Senescence Score by Group\nGSE224615 Long COVID PBMCs")
p_sen2 <- plot_scores(gsva_df, "CellAge",
  "CellAge Senescence Score by Group")
p_sen3 <- plot_scores(gsva_df, "SASP_Core",
  "SASP Core Score by Group")
p_sen4 <- plot_scores(gsva_df, "p21_p16_p53",
  "p21/p16/p53 Axis Score by Group")

combined_scores <- (p_sen1 + p_sen2) / (p_sen3 + p_sen4) +
  patchwork::plot_annotation(
    title    = "Sub-Aim 2.1: Senescence Transcriptomic Scores — GSE224615",
    subtitle = "Long COVID PBMC RNA-seq | GSVA method | All four groups shown",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "gray40", size = 9)
    )
  )

ggsave("results/senescence/fig6_senescence_scores_all_groups.pdf",
       combined_scores, width = 12, height = 9, dpi = 300)
ggsave("results/senescence/fig6_senescence_scores_all_groups.png",
       combined_scores, width = 12, height = 9, dpi = 300)
cat("Figure 6 saved.\n")

# ── 5. STATISTICAL TESTING: BRAIN FOG vs. NO BRAIN FOG ───────────────────────
cat("\n=== Statistical testing: BrainFog vs. NoBrainFog ===\n")

lc_scores <- gsva_df %>%
  dplyr::filter(group %in% c("LongCOVID_NoBrainFog", "LongCOVID_BrainFog"))

score_stats <- lapply(names(senescence_gene_sets), function(gs) {
  bf  <- lc_scores[[gs]][lc_scores$group == "LongCOVID_BrainFog"]
  nbf <- lc_scores[[gs]][lc_scores$group == "LongCOVID_NoBrainFog"]
  wt  <- wilcox.test(bf, nbf, exact = FALSE)
  data.frame(
    gene_set    = gs,
    n_brainfog  = length(bf),
    n_no_brainfog = length(nbf),
    median_bf   = round(median(bf, na.rm = TRUE), 4),
    median_nbf  = round(median(nbf, na.rm = TRUE), 4),
    delta       = round(median(bf, na.rm = TRUE) - median(nbf, na.rm = TRUE), 4),
    W_statistic = wt$statistic,
    p_value     = wt$p.value
  )
}) %>% dplyr::bind_rows()

# Adjust for multiple testing across gene sets
score_stats$padj <- p.adjust(score_stats$p_value, method = "BH")

cat("Senescence score comparison — Brain Fog vs. No Brain Fog:\n")
print(score_stats)

# ── 6. DIRECT GENE EXPRESSION: p21/p16/p53 AND SASP ─────────────────────────
cat("\n=== Direct expression of senescence axis genes ===\n")

# Extract VST expression of key senescence genes
key_sen_genes <- c("CDKN1A","CDKN2A","TP53","RB1",           # Cell cycle arrest
                   "IL6","CXCL8","MMP3","MMP9","CCL2",        # SASP
                   "SERPINE1","GDF15","IGFBP3","IGFBP7",       # SASP markers
                   "LMNB1","GLB1","HMGB1","BCL2L1")            # Senescence markers

key_sen_genes <- intersect(key_sen_genes, rownames(vst_mat))

sen_expr_df <- vst_mat[key_sen_genes, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "vst_expr") %>%
  dplyr::left_join(col_data, by = "sample_id")

# Heatmap of senescence gene expression
sen_heatmap_mat <- vst_mat[key_sen_genes, ]

# Annotation for heatmap columns
ann_col_df <- data.frame(Group = col_data$group, row.names = col_data$sample_id)
ann_colors  <- list(Group = grp_cols)

pheatmap::pheatmap(
  sen_heatmap_mat,
  scale          = "row",
  annotation_col = ann_col_df,
  annotation_colors = ann_colors,
  color          = colorRampPalette(c("#2196F3","white","#F44336"))(100),
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  show_rownames  = TRUE,
  show_colnames  = FALSE,
  fontsize_row   = 9,
  main           = "Senescence Axis Gene Expression — GSE224615",
  filename       = "results/senescence/fig7_senescence_gene_heatmap.pdf",
  width = 10, height = 7
)
cat("Figure 7 (senescence heatmap) saved.\n")

# ── 7. CORRELATION WITH COGNITIVE GENE EXPRESSION ────────────────────────────
cat("\n=== Correlating senescence with cognitive/synaptic genes ===\n")

# Cognitive-relevant genes from Aim 1 cell-intrinsic signature
cog_genes <- c("BDNF","NRXN1","SNAP25","SYP","GRIN3A",
                "CAMK2A","DLG4","SHANK3","SLC17A7","GABRB3")
cog_genes <- intersect(cog_genes, rownames(vst_mat))

if (length(cog_genes) > 0) {
  # Mean cognitive gene expression per sample
  cog_score <- colMeans(vst_mat[cog_genes, , drop = FALSE])

  cor_df <- data.frame(
    sample_id    = names(cog_score),
    cog_score    = as.numeric(cog_score),
    SenMayo      = as.numeric(gsva_scores["SenMayo", ]),
    SASP_Core    = as.numeric(gsva_scores["SASP_Core", ])
  ) %>% dplyr::left_join(col_data, by = "sample_id")

  # Senescence vs cognitive score scatter
  p_cor <- ggplot(cor_df, aes(x = SenMayo, y = cog_score, colour = group)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, colour = "black",
                linetype = "dashed", linewidth = 0.8) +
    ggrepel::geom_text_repel(aes(label = sample_id), size = 2.5,
                              max.overlaps = 15, segment.size = 0.3) +
    scale_colour_manual(values = grp_cols, name = "Group") +
    labs(
      title    = "SenMayo Score vs. Cognitive Gene Expression",
      subtitle = paste0("Cognitive genes: ", paste(cog_genes, collapse = ", ")),
      x        = "SenMayo Senescence Score",
      y        = "Mean Cognitive Gene Expression (VST)"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave("results/senescence/fig8_senescence_vs_cognitive_genes.pdf",
         p_cor, width = 9, height = 6, dpi = 300)
  ggsave("results/senescence/fig8_senescence_vs_cognitive_genes.png",
         p_cor, width = 9, height = 6, dpi = 300)

  # Compute Pearson correlation
  cor_test <- cor.test(cor_df$SenMayo, cor_df$cog_score, method = "pearson")
  cat("SenMayo vs. Cognitive gene score: r =", round(cor_test$estimate, 3),
      ", p =", round(cor_test$p.value, 4), "\n")
}

# ── 8. SAVE RESULTS ───────────────────────────────────────────────────────────
cat("\n=== Saving Sub-Aim 2.1 results ===\n")

saveRDS(gsva_scores, "data/gsva_senescence_scores.rds")
saveRDS(gsva_df,     "data/gsva_senescence_scores_df.rds")
write.csv(score_stats, "results/senescence/senescence_score_statistics.csv", row.names = FALSE)
write.csv(gsva_df,     "results/senescence/gsva_scores_per_sample.csv",      row.names = FALSE)

cat("\n=================================================================\n")
cat("  PHASE 3 (Sub-Aim 2.1) COMPLETE — SENESCENCE SCORING SUMMARY\n")
cat("=================================================================\n")
cat("Gene sets scored: SenMayo, CellAge, SASP_Core, p21/p16/p53\n")
cat("Method: GSVA (ssGSEA-equivalent enrichment)\n")
cat("\nStatistical results (BrainFog vs. NoBrainFog):\n")
print(score_stats[, c("gene_set","delta","p_value","padj")])
cat("\nNext: Run R/04_immunosenescence.R\n")
cat("=================================================================\n")
