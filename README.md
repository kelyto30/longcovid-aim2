# Long COVID Brain Aging — Big Aim 2
## Cellular Senescence and Immunosenescence as Transcriptomic Mediators of COVID-19 Cognitive Decline

**Author:** Kelechi Wisdom Elechi  
**Program:** PhD Candidate, Biology of Aging — Integrated Biomedical Sciences Program  
**Institution:** Barshop Institute for Longevity and Aging Studies, UT Health San Antonio  
**Year:** 2026

---

## Overview

This repository contains the complete reproducible analysis pipeline for **Big Aim 2** of the doctoral dissertation:

> *"COVID-19-Associated Brain Aging and Cognitive Decline: Transcriptional Regulatory Mechanisms, Cellular Composition, and Immune Aging"*

Big Aim 2 investigates whether Long COVID patients with cognitive symptoms (brain fog) carry a measurable transcriptomic burden of **cellular senescence** and **immunosenescence** in peripheral blood, and whether these aging signals converge with the postmortem brain transcriptomic signature established in Big Aim 1.

---

## Datasets

| Role | Accession | Description |
|------|-----------|-------------|
| Primary discovery | GSE224615 | Long COVID PBMC bulk RNA-seq (n=22); 4 groups: Healthy, COVID Convalescent, Long COVID, Long COVID Brain Fog |
| Cross-tissue validation | GSE188847 | Postmortem frontal cortex bulk RNA-seq (n=51 post-QC); analyzed in Big Aim 1 |

---

## Repository Structure

```
longcovid-aim2/
├── R/
│   ├── 01_download_qc.R          # Phase 1: GEO download + QC
│   ├── 02_differential_expression.R  # Phase 2: DESeq2 brain fog vs. no brain fog
│   ├── 03_senescence_scoring.R    # Phase 3: Sub-Aim 2.1 — SenMayo/CellAge GSVA
│   ├── 04_immunosenescence.R      # Phase 4: Sub-Aim 2.2 — Immunosenescence/exhaustion
│   ├── 05_coenrichment.R          # Phase 5: Co-enrichment correlation test
│   └── 06_brain_validation.R      # Phase 6: GSE188847 cross-tissue validation
├── data/
│   └── (downloaded GEO data — not tracked in git)
├── results/
│   ├── qc/                        # QC plots, sample summary tables
│   ├── deg/                       # DESeq2 results, volcano plots
│   ├── senescence/                # GSVA scores, scoring plots
│   ├── immunosenescence/          # ssGSEA scores, exhaustion figures
│   ├── integration/               # Co-enrichment, brain validation
│   └── figures/                   # Final publication-quality figures
└── docs/
    └── methods.md                 # Detailed methods documentation
```

---

## Sub-Aims

### Sub-Aim 2.1 — Senescence Transcriptomic Scoring
- **Gene sets:** SenMayo (105 genes; Saul et al. 2022, *Nature Aging*), CellAge database
- **Axis genes:** CDKN1A (p21), CDKN2A (p16), TP53; SASP: IL-6, CXCL8, MMP3, CXCL1
- **Method:** GSVA + ssGSEA (R/Bioconductor)
- **Comparison:** Long COVID Brain Fog vs. Long COVID (no brain fog)
- **Validation:** Top senescence genes tested in GSE188847 postmortem brain

### Sub-Aim 2.2 — Immunosenescence and T Cell Exhaustion
- **Immunosenescence panel:** KLRG1, B3GAT1 (CD57), TIGIT, LAG3, PDCD1 (PD-1), CD28 (loss)
- **Exhaustion TFs:** TOX, TOX2, NR4A1, NR4A2, NR4A3
- **Method:** ssGSEA, direct expression scoring
- **Co-enrichment test:** Immunosenescence score vs. Senescence score correlation
- **Validation:** Key transcripts tested in GSE188847

---

## How to Reproduce

```r
# Run scripts in order:
source("R/01_download_qc.R")           # ~15 min (downloads from GEO)
source("R/02_differential_expression.R")
source("R/03_senescence_scoring.R")
source("R/04_immunosenescence.R")
source("R/05_coenrichment.R")
source("R/06_brain_validation.R")
```

---

## Key R Packages Required

```r
# Data acquisition
BiocManager::install(c("GEOquery", "DESeq2", "tximport"))

# Scoring
BiocManager::install(c("GSVA", "GSEABase", "fgsea", "clusterProfiler"))

# Annotation
BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))

# Visualization
install.packages(c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
                   "patchwork", "viridis", "scales"))
```

---

## Connection to Big Aim 1

Big Aim 1 (GSE188847) established:
- **STAT1** (NES=+7.23) and **RELA** (NES=+5.44) as primary TF drivers of COVID-19 brain aging
- **86.7%** of the transcriptomic signal is composition-driven (neuronal loss + gliosis)
- **751-gene cell-intrinsic signature** in surviving neurons (MT2A, GRIN3A, KIF21B, STAT1)

Big Aim 2 tests whether these same aging mechanisms are detectable in the periphery of living Long COVID patients with cognitive symptoms — creating a blood-to-brain translational arc.

---

## Citation

*Elechi, K.W. (2026). COVID-19-Associated Brain Aging and Cognitive Decline. Doctoral Dissertation, UT Health San Antonio.*

---

## License

Academic use only. All GEO datasets are publicly available under NCBI/GEO terms of use.
