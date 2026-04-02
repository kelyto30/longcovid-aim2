# Methods Documentation — Big Aim 2

## Dataset
- **GSE224615**: Long COVID PBMC bulk RNA-seq, 22 samples, 4 groups
  - Healthy controls
  - COVID Convalescent
  - Long COVID (no cognitive symptoms)
  - Long COVID Brain Fog

## Analysis Pipeline

### Phase 1 (R/01_download_qc.R)
- GEO download via GEOquery
- VST normalization (DESeq2)
- PCA-based outlier detection (>3 SD from group centroid)

### Phase 2 (R/02_differential_expression.R)
- DESeq2: LongCOVID_BrainFog vs. LongCOVID_NoBrainFog
- Design: ~ age + sex + group
- LFC shrinkage: ashr
- GSEA: fgsea, GO:BP gene sets

### Phase 3 (R/03_senescence_scoring.R) — Sub-Aim 2.1
- SenMayo gene set (105 genes; Saul et al. 2022, Nature Aging)
- CellAge database gene sets
- GSVA scoring (ssGSEA-equivalent)
- Statistical comparison: Wilcoxon + BH correction

### Phase 4 (R/04_immunosenescence.R) — Sub-Aim 2.2
- Immunosenescence panel: KLRG1, B3GAT1(CD57), TIGIT, LAG3, PDCD1
- Exhaustion TFs: TOX, NR4A1/2/3, BATF, PRDM1
- Co-enrichment test: Pearson correlation of immunosenescence vs. senescence scores

### Phase 5-6 (R/05_06_integration_brain_validation.R)
- RELA/SASP molecular overlap between Aim 1 and Aim 2
- STAT1 blood-brain bridge analysis
- Cross-tissue validation in GSE188847 postmortem frontal cortex

## Key References
- Saul et al. (2022). SenMayo gene set. Nature Aging.
- Wherry & Kurachi (2015). T cell exhaustion. Nature Reviews Immunology.
- Garcia-Alonso et al. (2019). VIPER. Genome Research.
- McKenzie et al. (2018). BRETIGEA. Scientific Reports.
