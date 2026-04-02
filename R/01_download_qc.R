# =============================================================================
# Big Aim 2 — Phase 1: GEO Download, Preprocessing & Quality Control
# Dataset: GSE224615 (Long COVID PBMC RNA-seq)
# Groups: Healthy | COVID Convalescent | Long COVID | Long COVID Brain Fog
#
# Author: Kelechi Wisdom Elechi
# Institution: Barshop Institute for Longevity and Aging Studies
#              UT Health San Antonio
# Date: 2026
# =============================================================================

# ── 0. PACKAGE LOADING ────────────────────────────────────────────────────────
cat("=== Loading packages ===\n")

required_pkgs <- c("GEOquery", "DESeq2", "ggplot2", "ggrepel", "pheatmap",
                   "RColorBrewer", "dplyr", "tidyr", "stringr", "patchwork",
                   "viridis", "scales", "Biobase")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% rownames(installed.packages())) {
      library(pkg, character.only = TRUE)
    } else {
      message("Installing: ", pkg)
      if (pkg %in% c("GEOquery","DESeq2","Biobase")) {
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos = "https://cran.r-project.org")
      }
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("All packages loaded.\n\n")

# ── 1. DIRECTORY SETUP ────────────────────────────────────────────────────────
cat("=== Setting up directories ===\n")

dirs <- c("data", "results/qc", "results/deg", "results/senescence",
          "results/immunosenescence", "results/integration", "results/figures")

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat("Directories ready.\n\n")

# ── 2. DOWNLOAD GSE224615 FROM GEO ────────────────────────────────────────────
cat("=== Downloading GSE224615 from GEO ===\n")
cat("This may take 5-15 minutes depending on connection speed...\n")

gse_id <- "GSE224615"

# Download the series — getGEO fetches the series matrix file with metadata
gse <- tryCatch(
  GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE),
  error = function(e) {
    cat("ERROR downloading:", conditionMessage(e), "\n")
    stop(e)
  }
)

cat("Download complete.\n")
cat("Number of ExpressionSet objects:", length(gse), "\n\n")

# GSE224615 may return as a list — take the first element
gse_obj <- gse[[1]]

# ── 3. EXTRACT AND INSPECT METADATA ───────────────────────────────────────────
cat("=== Extracting metadata ===\n")

meta <- Biobase::pData(gse_obj)

cat("Total samples:", nrow(meta), "\n")
cat("Metadata columns:\n")
print(colnames(meta))
cat("\n")

# Display all characteristics columns to find group/cognitive status variable
char_cols <- colnames(meta)[grepl("characteristics", colnames(meta), ignore.case = TRUE)]
cat("Characteristics columns:\n")
for (col in char_cols) {
  cat("\n---", col, "---\n")
  print(table(meta[[col]]))
}

# ── 4. PARSE GROUP VARIABLE ───────────────────────────────────────────────────
cat("\n=== Parsing group assignments ===\n")

# GSE224615 has 4 groups: Healthy, COVID Convalescent, Long COVID, Long COVID Brain Fog
# The group variable is typically in characteristics_ch1 or similar
# We detect it automatically by searching for known group names

group_col <- NULL
for (col in char_cols) {
  vals <- tolower(as.character(meta[[col]]))
  if (any(grepl("brain fog|long covid|convalescent|healthy", vals))) {
    group_col <- col
    break
  }
}

# Also check 'source_name_ch1' and 'title' as fallbacks
if (is.null(group_col)) {
  for (col in c("source_name_ch1", "title", "description")) {
    if (col %in% colnames(meta)) {
      vals <- tolower(as.character(meta[[col]]))
      if (any(grepl("brain fog|long covid|convalescent|healthy", vals))) {
        group_col <- col
        break
      }
    }
  }
}

if (is.null(group_col)) {
  cat("WARNING: Could not automatically detect group column.\n")
  cat("Please inspect metadata above and set group_col manually.\n")
  # Manual fallback — print all columns for inspection
  cat("\nAll metadata for first 3 samples:\n")
  print(t(meta[1:min(3, nrow(meta)), ]))
} else {
  cat("Group column detected:", group_col, "\n")
  cat("Raw group values:\n")
  print(table(meta[[group_col]]))
}

# ── 5. CLEAN AND STANDARDIZE GROUP LABELS ─────────────────────────────────────
cat("\n=== Standardizing group labels ===\n")

if (!is.null(group_col)) {
  raw_group <- as.character(meta[[group_col]])

  # Standardize to clean labels
  meta$group <- dplyr::case_when(
    grepl("brain fog", tolower(raw_group))        ~ "LongCOVID_BrainFog",
    grepl("long covid|longcovid|long-covid", tolower(raw_group)) &
      !grepl("brain fog", tolower(raw_group))     ~ "LongCOVID_NoBrainFog",
    grepl("convalescent", tolower(raw_group))     ~ "COVID_Convalescent",
    grepl("healthy|control|hc", tolower(raw_group)) ~ "Healthy",
    TRUE ~ raw_group   # keep original if no match
  )

  # Factor with meaningful order
  meta$group <- factor(meta$group,
                       levels = c("Healthy", "COVID_Convalescent",
                                  "LongCOVID_NoBrainFog", "LongCOVID_BrainFog"))
} else {
  # Manual assignment fallback — edit below based on metadata inspection
  cat("MANUAL GROUP ASSIGNMENT REQUIRED\n")
  cat("Inspect metadata output above, then edit this section.\n")
  meta$group <- factor(rep(NA, nrow(meta)))
}

cat("Standardized group distribution:\n")
print(table(meta$group, useNA = "always"))

# ── 6. EXTRACT AGE AND SEX COVARIATES ─────────────────────────────────────────
cat("\n=== Extracting demographic covariates ===\n")

# Age extraction — search all columns
meta$age <- NA_real_
for (col in colnames(meta)) {
  if (any(grepl("^age:|age =", tolower(as.character(meta[[col]]))))) {
    meta$age <- as.numeric(gsub(".*age[: =]+", "", tolower(as.character(meta[[col]])),
                                ignore.case = TRUE))
    cat("Age found in column:", col, "\n")
    break
  }
}

# Sex extraction
meta$sex <- NA_character_
for (col in colnames(meta)) {
  vals <- tolower(as.character(meta[[col]]))
  if (any(grepl("^sex:|^gender:|male|female", vals))) {
    meta$sex <- dplyr::case_when(
      grepl("female", vals) ~ "F",
      grepl("male", vals)   ~ "M",
      TRUE ~ NA_character_
    )
    cat("Sex found in column:", col, "\n")
    break
  }
}

cat("Age summary:\n"); print(summary(meta$age))
cat("Sex distribution:\n"); print(table(meta$sex, useNA = "always"))

# ── 7. DOWNLOAD SUPPLEMENTARY COUNT MATRIX ────────────────────────────────────
cat("\n=== Downloading supplementary count matrix ===\n")

# Check what supplementary files are available
supp_info <- tryCatch(
  GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE,
                             baseDir = "data", fetch_files = FALSE),
  error = function(e) {
    cat("Could not list supplementary files:", conditionMessage(e), "\n")
    NULL
  }
)

if (!is.null(supp_info)) {
  cat("Available supplementary files:\n")
  print(supp_info$fname)
}

# Download all supplementary files into data/
cat("Downloading supplementary files to data/ ...\n")
tryCatch({
  GEOquery::getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = "data")
  cat("Supplementary download complete.\n")
}, error = function(e) {
  cat("ERROR downloading supplementary files:", conditionMessage(e), "\n")
  cat("You may need to download manually from:\n")
  cat("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_id, "\n")
})

# List downloaded files
cat("\nFiles in data/:\n")
print(list.files("data/", recursive = TRUE))

# ── 8. LOAD COUNT MATRIX ──────────────────────────────────────────────────────
cat("\n=== Loading count matrix ===\n")

# Find count matrix file (handles various common naming conventions)
data_files <- list.files("data/", pattern = "\\.txt\\.gz$|\\.csv\\.gz$|\\.tsv\\.gz$|\\.count.*gz$",
                          full.names = TRUE, recursive = TRUE)

cat("Candidate count files found:\n")
print(data_files)

if (length(data_files) == 0) {
  cat("\nNo count matrix found automatically.\n")
  cat("Please download the count matrix manually from GEO and place in data/\n")
  cat("Then reload this script from Step 8.\n")
  stop("Count matrix not found — see message above.")
}

# Load the first matching file (adjust index if needed)
count_file <- data_files[1]
cat("Loading:", count_file, "\n")

counts_raw <- tryCatch({
  read.table(gzfile(count_file), header = TRUE, sep = "\t",
             row.names = 1, check.names = FALSE)
}, error = function(e) {
  # Try comma separator
  read.csv(gzfile(count_file), header = TRUE, row.names = 1, check.names = FALSE)
})

cat("Count matrix dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
cat("First 5 genes, first 3 samples:\n")
print(counts_raw[1:5, 1:min(3, ncol(counts_raw))])

# ── 9. ALIGN METADATA WITH COUNT MATRIX ───────────────────────────────────────
cat("\n=== Aligning samples between metadata and count matrix ===\n")

# GEO sample IDs should match between the count matrix columns and metadata
common_samples <- intersect(colnames(counts_raw), rownames(meta))

if (length(common_samples) == 0) {
  # Try matching by GSM accession from column names
  cat("Direct ID match failed — attempting GSM-based matching...\n")
  # Some count matrices use GSM IDs in column names
  gsm_ids <- meta$geo_accession
  col_names <- colnames(counts_raw)

  # Try partial matching
  matched <- match(gsm_ids, col_names)
  if (any(!is.na(matched))) {
    common_samples <- gsm_ids[!is.na(matched)]
    counts_raw <- counts_raw[, col_names[matched[!is.na(matched)]]]
    meta_aligned <- meta[!is.na(matched), ]
    cat("GSM-based match found:", length(common_samples), "samples\n")
  } else {
    cat("WARNING: Could not match samples between metadata and count matrix.\n")
    cat("Column names in count matrix:\n"); print(head(colnames(counts_raw)))
    cat("Row names in metadata:\n"); print(head(rownames(meta)))
    cat("Please inspect and manually align.\n")
    meta_aligned <- meta
  }
} else {
  counts_aligned <- counts_raw[, common_samples, drop = FALSE]
  meta_aligned   <- meta[common_samples, ]
  cat("Aligned:", length(common_samples), "samples matched between metadata and counts.\n")
}

# Final count matrix and metadata
counts <- counts_aligned
cat("Final count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")
cat("Final metadata:", nrow(meta_aligned), "samples\n")

# ── 10. CREATE DESEQ2 OBJECT ──────────────────────────────────────────────────
cat("\n=== Creating DESeq2 object ===\n")

# Ensure counts are integers
counts_int <- round(counts)
storage.mode(counts_int) <- "integer"

# Remove genes with zero counts across all samples
keep_genes <- rowSums(counts_int) > 0
counts_int <- counts_int[keep_genes, ]
cat("Genes after removing all-zero rows:", nrow(counts_int), "\n")

# Build colData
col_data <- data.frame(
  sample_id = colnames(counts_int),
  group     = meta_aligned$group,
  age       = meta_aligned$age,
  sex       = meta_aligned$sex,
  row.names = colnames(counts_int),
  stringsAsFactors = FALSE
)

col_data$group <- factor(col_data$group,
                         levels = c("Healthy", "COVID_Convalescent",
                                    "LongCOVID_NoBrainFog", "LongCOVID_BrainFog"))

cat("colData summary:\n")
print(col_data)

# Create DESeq2 object
# Design uses group as primary factor; age/sex added if available
has_age <- !all(is.na(col_data$age))
has_sex <- !all(is.na(col_data$sex))

if (has_age & has_sex) {
  design_formula <- ~ age + sex + group
  cat("Design formula: ~ age + sex + group\n")
} else if (has_age) {
  design_formula <- ~ age + group
  cat("Design formula: ~ age + group\n")
} else {
  design_formula <- ~ group
  cat("Design formula: ~ group (no age/sex covariates available)\n")
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData   = col_data,
  design    = design_formula
)

cat("DESeq2 object created successfully.\n")
cat("DESeq2 object summary:\n"); print(dds)

# ── 11. QUALITY CONTROL ───────────────────────────────────────────────────────
cat("\n=== Running Quality Control ===\n")

# Low-count gene filtering
cat("Before filtering:", nrow(dds), "genes\n")
keep <- rowSums(counts(dds) >= 10) >= 3   # at least 3 samples with >= 10 counts
dds_filt <- dds[keep, ]
cat("After filtering (>=10 counts in >=3 samples):", nrow(dds_filt), "genes\n")

# VST normalization for QC visualization
cat("Running variance-stabilizing transformation...\n")
vst_obj  <- vst(dds_filt, blind = TRUE)
vst_mat  <- assay(vst_obj)

# ── Group colour palette ──
grp_cols <- c(
  "Healthy"               = "#2196F3",
  "COVID_Convalescent"    = "#FF9800",
  "LongCOVID_NoBrainFog"  = "#9C27B0",
  "LongCOVID_BrainFog"    = "#F44336"
)

# ── FIGURE 1: Library Size Distribution ──
cat("Generating Figure 1: Library size distribution...\n")

lib_df <- data.frame(
  sample = colnames(counts_int),
  library_size = colSums(counts_int) / 1e6,
  group = col_data$group
) %>%
  dplyr::arrange(group, library_size)

lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample)

p_lib <- ggplot(lib_df, aes(x = sample, y = library_size, fill = group)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 5, linetype = "dashed", colour = "red",
             linewidth = 0.6) +
  scale_fill_manual(values = grp_cols, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Library Size Distribution — GSE224615 (Long COVID PBMC)",
    subtitle = "Dashed red line = 5 million read minimum threshold",
    x        = "Sample",
    y        = "Library Size (millions of reads)",
    caption  = "QC Step 1: Samples below threshold flagged for exclusion"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "gray40")
  )

ggsave("results/qc/fig1_library_sizes.pdf", p_lib,
       width = 10, height = 5, dpi = 300)
ggsave("results/qc/fig1_library_sizes.png", p_lib,
       width = 10, height = 5, dpi = 300)
cat("Figure 1 saved.\n")

# ── FIGURE 2: PCA Plot ──
cat("Generating Figure 2: PCA plot...\n")

pca_res  <- prcomp(t(vst_mat), scale. = FALSE)
pca_var  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- data.frame(
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  group  = col_data$group,
  sample = colnames(vst_mat)
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group, label = sample)) +
  geom_point(size = 4, alpha = 0.85, stroke = 0.5) +
  ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20,
                            segment.colour = "grey60", segment.size = 0.3) +
  stat_ellipse(aes(colour = group), level = 0.80,
               linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = grp_cols, name = "Group") +
  labs(
    title    = "PCA of VST-Normalized Gene Expression — GSE224615",
    subtitle = "Long COVID PBMC RNA-seq | 80% confidence ellipses",
    x        = paste0("PC1 (", pca_var[1], "% variance)"),
    y        = paste0("PC2 (", pca_var[2], "% variance)"),
    caption  = "QC Step 2: Samples >3 SD from group centroid flagged as outliers"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "gray40")
  )

ggsave("results/qc/fig2_pca.pdf", p_pca, width = 8, height = 6, dpi = 300)
ggsave("results/qc/fig2_pca.png", p_pca, width = 8, height = 6, dpi = 300)
cat("Figure 2 saved.\n")

# ── FIGURE 3: Sample Distance Heatmap ──
cat("Generating Figure 3: Sample distance heatmap...\n")

samp_dist <- dist(t(vst_mat))
samp_dist_mat <- as.matrix(samp_dist)

ann_df <- data.frame(Group = col_data$group, row.names = colnames(vst_mat))
ann_cols_list <- list(Group = grp_cols)

pheatmap::pheatmap(
  samp_dist_mat,
  clustering_distance_rows = samp_dist,
  clustering_distance_cols = samp_dist,
  annotation_col  = ann_df,
  annotation_row  = ann_df,
  annotation_colors = ann_cols_list,
  color           = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize        = 8,
  fontsize_row    = 7,
  fontsize_col    = 7,
  main            = "Sample-to-Sample Distance Heatmap — GSE224615",
  filename        = "results/qc/fig3_sample_distance_heatmap.pdf",
  width           = 8,
  height          = 7
)
cat("Figure 3 saved.\n")

# ── 12. OUTLIER DETECTION ─────────────────────────────────────────────────────
cat("\n=== Outlier Detection ===\n")

# Detect outliers > 3 SD from group centroid on PC1/PC2
outlier_results <- pca_df %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    pc1_sd   = sd(PC1, na.rm = TRUE),
    pc2_sd   = sd(PC2, na.rm = TRUE),
    pc1_mean = mean(PC1, na.rm = TRUE),
    pc2_mean = mean(PC2, na.rm = TRUE),
    pc1_dev  = abs(PC1 - pc1_mean) / pc1_sd,
    pc2_dev  = abs(PC2 - pc2_mean) / pc2_sd,
    outlier  = (pc1_dev > 3 | pc2_dev > 3)
  ) %>%
  dplyr::ungroup()

outliers_detected <- outlier_results %>%
  dplyr::filter(outlier == TRUE) %>%
  dplyr::select(sample, group, pc1_dev, pc2_dev)

if (nrow(outliers_detected) == 0) {
  cat("No outliers detected (all samples within 3 SD of group centroid).\n")
} else {
  cat("OUTLIERS DETECTED:\n")
  print(outliers_detected)
  cat("\nThese samples are flagged for review.\n")
}

# ── 13. SAVE PROCESSED OBJECTS ────────────────────────────────────────────────
cat("\n=== Saving processed objects ===\n")

# Save main objects for downstream scripts
saveRDS(dds_filt,    "data/dds_gse224615.rds")
saveRDS(vst_mat,     "data/vst_matrix_gse224615.rds")
saveRDS(meta_aligned,"data/metadata_gse224615.rds")
saveRDS(col_data,    "data/coldata_gse224615.rds")
saveRDS(outlier_results, "data/qc_pca_outliers.rds")

# Save metadata as CSV for easy inspection
write.csv(col_data, "results/qc/sample_metadata_cleaned.csv", row.names = TRUE)

# Save QC summary table
qc_summary <- data.frame(
  Step        = c("GEO download", "Raw samples",
                  "Low-count gene filter (>=10 in >=3 samples)",
                  "Outliers detected", "Final samples", "Final genes"),
  Value       = c(gse_id, nrow(meta), nrow(dds_filt),
                  nrow(outliers_detected), ncol(counts_int), nrow(dds_filt))
)
write.csv(qc_summary, "results/qc/qc_summary.csv", row.names = FALSE)

cat("\nObjects saved:\n")
cat("  data/dds_gse224615.rds\n")
cat("  data/vst_matrix_gse224615.rds\n")
cat("  data/metadata_gse224615.rds\n")
cat("  data/coldata_gse224615.rds\n")
cat("  data/qc_pca_outliers.rds\n")
cat("  results/qc/sample_metadata_cleaned.csv\n")
cat("  results/qc/qc_summary.csv\n")

# ── 14. PHASE 1 SUMMARY ───────────────────────────────────────────────────────
cat("\n")
cat("=================================================================\n")
cat("  PHASE 1 COMPLETE — GSE224615 QC SUMMARY\n")
cat("=================================================================\n")
cat("Dataset          :", gse_id, "\n")
cat("Total samples    :", nrow(meta_aligned), "\n")
cat("Group breakdown  :\n")
print(table(col_data$group))
cat("Genes after QC   :", nrow(dds_filt), "\n")
cat("Outliers flagged :", nrow(outliers_detected), "\n")
cat("Design formula   :", deparse(design_formula), "\n")
cat("\nOutputs:\n")
cat("  Figures → results/qc/\n")
cat("  RDS objects → data/\n")
cat("\nNext: Run R/02_differential_expression.R\n")
cat("=================================================================\n")
