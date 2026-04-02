# =============================================================================
# Big Aim 2 — Phase 2: Differential Expression Analysis
# Comparison: LongCOVID_BrainFog vs. LongCOVID_NoBrainFog
# Dataset: GSE224615 (Long COVID PBMC RNA-seq)
#
# Author: Kelechi Wisdom Elechi
# Institution: Barshop Institute for Longevity and Aging Studies
#              UT Health San Antonio
# Date: 2026
# =============================================================================

cat("=== Loading packages ===\n")
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel)
  library(dplyr); library(RColorBrewer); library(pheatmap)
  library(clusterProfiler); library(org.Hs.eg.db); library(fgsea)
  library(patchwork); library(scales)
})

# ── 1. LOAD PHASE 1 OBJECTS ───────────────────────────────────────────────────
cat("=== Loading Phase 1 objects ===\n")

dds      <- readRDS("data/dds_gse224615.rds")
vst_mat  <- readRDS("data/vst_matrix_gse224615.rds")
col_data <- readRDS("data/coldata_gse224615.rds")

cat("Samples loaded:", ncol(dds), "\n")
cat("Groups:\n"); print(table(col_data$group))

# Colour palette (consistent across all scripts)
grp_cols <- c(
  "Healthy"               = "#2196F3",
  "COVID_Convalescent"    = "#FF9800",
  "LongCOVID_NoBrainFog"  = "#9C27B0",
  "LongCOVID_BrainFog"    = "#F44336"
)

# ── 2. PRIMARY COMPARISON: BRAIN FOG vs. NO BRAIN FOG ────────────────────────
cat("\n=== Primary DESeq2 comparison: BrainFog vs. NoBrainFog ===\n")

# Subset to Long COVID groups only
lc_samples <- col_data$sample_id[col_data$group %in%
                c("LongCOVID_NoBrainFog", "LongCOVID_BrainFog")]

dds_lc <- dds[, lc_samples]
dds_lc$group <- droplevels(dds_lc$group)
dds_lc$group <- relevel(dds_lc$group, ref = "LongCOVID_NoBrainFog")

cat("Samples in primary comparison:\n")
print(table(dds_lc$group))

# Run DESeq2
dds_lc <- DESeq(dds_lc)

# Extract results: BrainFog vs NoBrainFog
res_bf <- results(dds_lc,
                  contrast  = c("group", "LongCOVID_BrainFog", "LongCOVID_NoBrainFog"),
                  alpha     = 0.05,
                  pAdjustMethod = "BH")

res_bf <- lfcShrink(dds_lc,
                    contrast  = c("group", "LongCOVID_BrainFog", "LongCOVID_NoBrainFog"),
                    type      = "ashr",
                    res       = res_bf)

# Summary
cat("\nDESeq2 results summary (BrainFog vs. NoBrainFog):\n")
summary(res_bf)

res_bf_df <- as.data.frame(res_bf) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(padj, desc(abs(log2FoldChange)))

deg_sig <- res_bf_df %>% dplyr::filter(padj < 0.05 & !is.na(padj))
cat("Significant DEGs (padj < 0.05):", nrow(deg_sig), "\n")
cat("  Upregulated in BrainFog:", sum(deg_sig$log2FoldChange > 0), "\n")
cat("  Downregulated in BrainFog:", sum(deg_sig$log2FoldChange < 0), "\n")

# ── 3. SECONDARY COMPARISON: LONG COVID vs. HEALTHY ──────────────────────────
cat("\n=== Secondary comparison: LongCOVID vs. Healthy ===\n")

dds_full <- dds
dds_full$group <- relevel(dds_full$group, ref = "Healthy")
dds_full <- DESeq(dds_full)

res_lc_v_hc <- results(dds_full,
                        contrast      = c("group", "LongCOVID_BrainFog", "Healthy"),
                        alpha         = 0.05,
                        pAdjustMethod = "BH")

res_lc_v_hc_df <- as.data.frame(res_lc_v_hc) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(padj)

cat("LongCOVID BrainFog vs. Healthy DEGs (padj<0.05):",
    sum(!is.na(res_lc_v_hc_df$padj) & res_lc_v_hc_df$padj < 0.05), "\n")

# ── 4. VOLCANO PLOT — PRIMARY COMPARISON ─────────────────────────────────────
cat("\n=== Generating Volcano Plot ===\n")

volcano_df <- res_bf_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up in BrainFog",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in BrainFog",
      padj < 0.05                        ~ "Sig (|LFC|<1)",
      TRUE                               ~ "NS"
    ),
    neg_log10p = -log10(padj + 1e-300)
  )

# Label top genes
top_genes <- volcano_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 20)

vol_cols <- c("Up in BrainFog" = "#F44336", "Down in BrainFog" = "#2196F3",
              "Sig (|LFC|<1)"   = "#FF9800", "NS"              = "grey80")

p_vol <- ggplot(volcano_df, aes(x = log2FoldChange, y = neg_log10p,
                                 colour = sig_label)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_point(data = top_genes,
             aes(x = log2FoldChange, y = -log10(padj + 1e-300)),
             size = 2, shape = 21, fill = "white", stroke = 0.8) +
  ggrepel::geom_text_repel(data = top_genes,
                            aes(label = gene), size = 2.5,
                            max.overlaps = 30, segment.size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  scale_colour_manual(values = vol_cols, name = NULL) +
  labs(
    title    = "Differential Gene Expression: Long COVID Brain Fog vs. No Brain Fog",
    subtitle = "GSE224615 PBMC RNA-seq | DESeq2 + apeglm shrinkage | padj < 0.05",
    x        = "log2 Fold Change (Brain Fog / No Brain Fog)",
    y        = "-log10(adjusted p-value)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 11),
    plot.subtitle   = element_text(size = 8, colour = "gray40"),
    legend.position = "top"
  )

ggsave("results/deg/fig4_volcano_brainfog.pdf", p_vol, width = 9, height = 6, dpi = 300)
ggsave("results/deg/fig4_volcano_brainfog.png", p_vol, width = 9, height = 6, dpi = 300)
cat("Volcano plot saved.\n")

# ── 5. GSEA ON RANKED GENE LIST ───────────────────────────────────────────────
cat("\n=== Running GSEA ===\n")

# Ranked by Wald stat for GSEA
wald_stats <- res_bf_df %>%
  dplyr::filter(!is.na(stat)) %>%
  dplyr::arrange(desc(stat))

ranked_genes <- wald_stats$stat
names(ranked_genes) <- wald_stats$gene

# Load GO:BP gene sets
go_sets <- msigdbr::msigdbr(species = "Homo sapiens",
                              category = "C5", subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(.$gs_name) %>%
  lapply(function(x) x$gene_symbol)

# Run fgsea
set.seed(42)
gsea_res <- fgsea::fgsea(pathways   = go_sets,
                          stats      = ranked_genes,
                          minSize    = 10,
                          maxSize    = 500,
                          nPermSimple = 1000)

gsea_sig <- gsea_res %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(NES)

cat("Significant GSEA terms (padj<0.05):", nrow(gsea_sig), "\n")
cat("Top upregulated (NES > 0):\n")
print(head(gsea_sig %>% dplyr::filter(NES > 0) %>% dplyr::arrange(padj), 5))
cat("Top downregulated (NES < 0):\n")
print(head(gsea_sig %>% dplyr::filter(NES < 0) %>% dplyr::arrange(padj), 5))

# GSEA dot plot — top 20 terms
top_gsea <- gsea_sig %>%
  dplyr::group_by(NES > 0) %>%
  dplyr::slice_min(padj, n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(pathway = stringr::str_replace_all(pathway, "GOBP_", "") %>%
                  stringr::str_replace_all("_", " ") %>%
                  tolower() %>%
                  stringr::str_to_sentence())

p_gsea <- ggplot(top_gsea, aes(x = NES, y = reorder(pathway, NES),
                                colour = -log10(padj), size = size)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_viridis_c(name = "-log10(padj)", option = "plasma") +
  scale_size_continuous(name = "Gene set size", range = c(2, 8)) +
  labs(
    title    = "GSEA: Long COVID Brain Fog vs. No Brain Fog",
    subtitle = "Gene Ontology Biological Process | padj < 0.05",
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y     = element_text(size = 7),
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(colour = "gray40", size = 8)
  )

ggsave("results/deg/fig5_gsea_brainfog.pdf", p_gsea, width = 11, height = 7, dpi = 300)
ggsave("results/deg/fig5_gsea_brainfog.png", p_gsea, width = 11, height = 7, dpi = 300)
cat("GSEA plot saved.\n")

# ── 6. SAVE RESULTS ───────────────────────────────────────────────────────────
cat("\n=== Saving DEG results ===\n")

write.csv(res_bf_df,          "results/deg/deseq2_brainfog_vs_noBrainfog.csv", row.names = FALSE)
write.csv(deg_sig,            "results/deg/significant_DEGs_padj05.csv",        row.names = FALSE)
write.csv(res_lc_v_hc_df,     "results/deg/deseq2_brainfog_vs_healthy.csv",     row.names = FALSE)
write.csv(as.data.frame(gsea_sig), "results/deg/gsea_results_brainfog.csv",     row.names = FALSE)

saveRDS(dds_lc,     "data/dds_lc_subset.rds")
saveRDS(dds_full,   "data/dds_full_allgroups.rds")
saveRDS(res_bf_df,  "data/deseq2_results_brainfog.rds")
saveRDS(ranked_genes, "data/wald_ranked_genes.rds")

cat("\n=================================================================\n")
cat("  PHASE 2 COMPLETE — DIFFERENTIAL EXPRESSION SUMMARY\n")
cat("=================================================================\n")
cat("Primary comparison: LongCOVID BrainFog vs. NoBrainFog\n")
cat("Total DEGs (padj<0.05):", nrow(deg_sig), "\n")
cat("  Upregulated         :", sum(deg_sig$log2FoldChange > 0), "\n")
cat("  Downregulated       :", sum(deg_sig$log2FoldChange < 0), "\n")
cat("GSEA significant terms:", nrow(gsea_sig), "\n")
cat("\nNext: Run R/03_senescence_scoring.R\n")
cat("=================================================================\n")
