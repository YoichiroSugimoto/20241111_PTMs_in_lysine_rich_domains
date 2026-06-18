p2-09 · Export stoichiometry tables
================
Yoichiro Sugimoto and Pallavi Kesavan
18 June, 2026

- [Overview](#overview)
- [Setup](#setup)
- [Export tables to Excel](#export-tables-to-excel)
- [Session information](#session-information)

# Overview

**Purpose:** Export processed BRD stoichiometry (MS_SS, MS_KR_1, PNAS)
to a formatted Excel workbook (supplementary data).

**Inputs:** stoichiometry tables from p2-01; reference proteome.

**Outputs:** Excel workbook under `results/`.

**Upstream:** p2-01. **Downstream:** none.

# Setup

``` r
## Resolve the repository root (via the .here sentinel) and load the shared
## setup: packages, helper functions, ggplot/knitr settings, and project paths.
repo_root <- local({
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  while (!file.exists(file.path(p, ".here")) && !identical(dirname(p), p)) p <- dirname(p)
  p
})
source(file.path(repo_root, "R", "functions", "_setup.R"))
library("janitor")
library("openxlsx")
```

# Export tables to Excel

Assembles the per-site stoichiometry tables (PNAS data-A/B/C, MS_KR_1,
MS_SS) and writes them to a formatted Excel workbook for the
supplementary data.

``` r
# A diagnostic-peak flag shared by all loaded tables
add_di_flag <- function(dt) dt[, is_diagnostic_peak := diagnostic_peak == "+"][]

mq_di_dir <- file.path(results.dir, "p2-analysis-setting", "MQ_DI")

#--- PNAS dataset A (MaxQuant with diagnostic ions) ---
mq_di_a_stoic_dt <- read_stoic_data(
  prefix = "DI_data-A_trp_m7_v7_def_", pre_prefix = "", post_fix = "_DI",
  dir_path = mq_di_dir
) %>% add_di_flag

#--- PNAS datasets B1 and B2 ---
mq_di_sample_dt <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_DI"
) %>% data.table
mq_di_sample_dt[, sample_id := 1:.N]

mq_di_b_stoic_dt <- lapply(
  mq_di_sample_dt[data != c("data-A", "data-C"), prefix],
  read_stoic_data,
  pre_prefix = "DI_", post_fix = "_DI", dir_path = mq_di_dir
) %>% rbindlist %>% add_di_flag

#--- PNAS dataset C ---
mq_di_c_stoic_dt <- read_stoic_data(
  prefix = "DI_data-C_trp_m7_v7_mCC_", pre_prefix = "", post_fix = "_DI",
  dir_path = mq_di_dir
) %>% add_di_flag

#--- PNAS manually curated sites ---
pnas2022_stoic_dt <- fread(file.path(data.dir, "processed_data_from_PNAS2022/long_K_stoichiometry_data.csv"))
setnames(
  pnas2022_stoic_dt,
  old = c("uniprot_id", "position", "residue", "oxK_ratio", "data_source", "total_n_feature_oxK", "total_n_feature_K"),
  new = c("protein_accession", "aa_pos", "aa", "stoichiometry", "sample_name", "sum_psm_mapped", "sum_psm_mapped_per_position")
)
pnas2022_stoic_dt <- merge(
  ref_protein_dt[, .(protein_accession, gene_name)],
  pnas2022_stoic_dt,
  by = "protein_accession"
)
pnas2022_stoic_dt[, `:=`(
  accession_position = paste0(protein_accession, "_", aa_pos),
  ptm = fifelse(curated_oxK_site == TRUE, "[Oxidation (K)]", " ")
)]
# Harmonise column names with the exported tables
setnames(
  pnas2022_stoic_dt,
  c("stoichiometry", "protein_accession", "gene_name", "aa_pos", "curated_oxK_site"),
  c("Stoichiometry_PNAS2022", "Protein_accession", "Gene_name", "Amino_acid_pos", "PNAS2022_curated_oxK_site")
)

#--- MS_KR_1 (dataset D) ---
ms_kr1_stoic_dt <- read_stoic_data(
  prefix = "MS_KR_1_", pre_prefix = "", post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MS_KR_1")
) %>% add_di_flag

#--- MS_SS (dataset E) ---
ms_ss_stoic_dt <- read_stoic_data(
  prefix = "MS_SS_", pre_prefix = "", post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MS_SS")
) %>% add_di_flag
```

``` r
# Output directory for the supplementary workbooks (created if missing)
export_out_dir <- file.path(results.dir, "p2_MS_SS_KR")
create.dir(export_out_dir)

# Working columns dropped from every stoichiometry table before export
export_drop_cols <- c(
  "sum_intensity_per_position", "sum_psm_mapped_per_position",
  "sum_the_number_of_peptide", "localization_prob", "score_for_localization",
  "best_localization_ms_ms_id", "best_localization_raw_file", "max_score",
  "diagnostic_peak", "condition"
)

# Publication-ready column renaming, applied identically to every sheet
export_old_names <- c(
  "sample_name", "protein_accession", "gene_name", "aa", "aa_pos", "ptm",
  "stoichiometry", "sum_peak_intensity", "sum_psm_mapped",
  "the_number_of_peptide", "is_diagnostic_peak"
)
export_new_names <- c(
  "Sample_Name", "Protein_accession", "Gene_name", "Amino_acid", "Amino_acid_pos",
  "PTM_type", "Stoichiometry", "Total_peak_intensity", "Total_PSM_mapped",
  "Number_of_peptide", "Diagnostic_peak_presence"
)

## Format a per-site stoichiometry table for the supplementary workbook and write
## it to <results.dir>/p2_MS_SS_KR/<file_name>. Subsets to lysines, drops working
## columns, reorders + renames to publication headers, flags PNAS2022 curated
## sites, and writes the xlsx. `extra_filter` optionally transforms the table
## (original column names still in place) before reordering — used for the MS_SS
## cross-protein clean-up.
export_stoic_table <- function(stoic_dt, file_name, curated_dt, extra_filter = NULL) {
  dt <- stoic_dt[aa == "K"]                       # subset to lysines (copy)
  dt[, (export_drop_cols) := NULL]                # drop working columns

  if (!is.null(extra_filter)) dt <- extra_filter(dt)

  # Move stoichiometry / diagnostic-peak flag to their publication positions
  setcolorder(dt, append(setdiff(names(dt), "stoichiometry"), "stoichiometry", after = 6))
  setcolorder(dt, append(setdiff(names(dt), "is_diagnostic_peak"), "is_diagnostic_peak", after = 13))

  setnames(dt, export_old_names, export_new_names)  # publication-ready headers

  # Flag sites manually curated as oxK in PNAS2022
  dt[, PNAS2022_curated_oxK_site :=
       paste0(Protein_accession, "_", Amino_acid_pos) %in%
       curated_dt[PNAS2022_curated_oxK_site == TRUE,
                  paste0(Protein_accession, "_", Amino_acid_pos)]]

  write.xlsx(dt, file.path(export_out_dir, file_name))
  invisible(dt)
}
```

``` r
# PNAS dataset A (diagnostic ions)
export_stoic_table(mq_di_a_stoic_dt, "MQ_DI_dtA_stoic_dt.xlsx", pnas2022_stoic_dt)
```

``` r
# PNAS datasets B1 and B2 (diagnostic ions)
export_stoic_table(mq_di_b_stoic_dt, "MQ_DI_dtB_stoic_dt.xlsx", pnas2022_stoic_dt)
```

``` r
# PNAS dataset C (diagnostic ions)
export_stoic_table(mq_di_c_stoic_dt, "MQ_DI_dtC_stoic_dt.xlsx", pnas2022_stoic_dt)
```

``` r
# MS_KR_1 (dataset D)
export_stoic_table(ms_kr1_stoic_dt, "MS_KR1_stoic_K.xlsx", pnas2022_stoic_dt)
```

``` r
# MS_SS (dataset E). Drops cross-protein false positives (BRD4 signal reported in
# BRD2/3 samples and vice versa, e.g. from shared / homologous peptides) before
# formatting.
export_stoic_table(
  ms_ss_stoic_dt, "MS_SS_stoic_K.xlsx", pnas2022_stoic_dt,
  extra_filter = function(dt) dt[
    !((grepl("_BRD23$", sample_name) & gene_name == "BRD4") |
        (grepl("_BRD4$", sample_name) & gene_name %in% c("BRD2", "BRD3")))
  ]
)
```

# Session information

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.4.3 (2025-02-28)
    ##  os       Ubuntu 24.04.2 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       Europe/Berlin
    ##  date     2026-06-18
    ##  pandoc   3.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.5.57 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  BiocGenerics      * 0.52.0     2024-10-29 [1] Bioconduc~
    ##  Biostrings        * 2.74.1     2024-12-16 [1] Bioconduc~
    ##  cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.4.3)
    ##  cli                 3.6.5      2025-04-23 [1] CRAN (R 4.4.3)
    ##  crayon              1.5.3      2024-06-20 [1] CRAN (R 4.4.3)
    ##  data.table        * 1.17.8     2025-07-10 [1] CRAN (R 4.4.3)
    ##  digest              0.6.37     2024-08-19 [1] CRAN (R 4.4.3)
    ##  dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.4.3)
    ##  evaluate            1.0.5      2025-08-27 [1] CRAN (R 4.4.3)
    ##  farver              2.1.2      2024-05-13 [1] CRAN (R 4.4.3)
    ##  fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.4.3)
    ##  generics            0.1.4      2025-05-09 [1] CRAN (R 4.4.3)
    ##  GenomeInfoDb      * 1.42.3     2025-01-27 [1] Bioconduc~
    ##  GenomeInfoDbData    1.2.13     2025-07-21 [1] Bioconductor
    ##  ggplot2           * 4.0.0      2025-09-11 [1] CRAN (R 4.4.3)
    ##  glue                1.8.0      2024-09-30 [1] CRAN (R 4.4.3)
    ##  gtable              0.3.6      2024-10-25 [1] CRAN (R 4.4.3)
    ##  htmltools           0.5.8.1    2024-04-04 [1] CRAN (R 4.4.3)
    ##  httr                1.4.7      2023-08-15 [1] CRAN (R 4.4.3)
    ##  IRanges           * 2.40.1     2024-12-05 [1] Bioconduc~
    ##  janitor           * 2.2.1      2024-12-22 [1] CRAN (R 4.4.3)
    ##  jsonlite            2.0.0      2025-03-27 [1] CRAN (R 4.4.3)
    ##  khroma            * 1.16.0     2025-02-25 [1] CRAN (R 4.4.3)
    ##  knitr             * 1.50       2025-03-16 [1] CRAN (R 4.4.3)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.4.3)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.4.3)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.4.3)
    ##  openxlsx          * 4.2.8.1    2025-10-31 [1] CRAN (R 4.4.3)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.4.3)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.3)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-12-13 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.4.3)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.4.3)
    ##  Rcpp                1.1.0      2025-07-02 [1] CRAN (R 4.4.3)
    ##  readxl            * 1.4.5      2025-03-07 [1] CRAN (R 4.4.3)
    ##  rlang               1.1.6      2025-04-11 [1] CRAN (R 4.4.3)
    ##  rmarkdown           2.29       2024-11-04 [1] CRAN (R 4.4.3)
    ##  rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.4.3)
    ##  S4Vectors         * 0.44.0     2024-10-29 [1] Bioconduc~
    ##  S7                  0.2.0      2024-11-07 [1] CRAN (R 4.4.3)
    ##  scales              1.4.0      2025-04-24 [1] CRAN (R 4.4.3)
    ##  sessioninfo         1.2.3      2025-02-05 [1] CRAN (R 4.4.3)
    ##  snakecase           0.11.1     2023-08-27 [1] CRAN (R 4.4.3)
    ##  stringi             1.8.7      2025-03-27 [1] CRAN (R 4.4.3)
    ##  stringr           * 1.5.2      2025-09-08 [1] CRAN (R 4.4.3)
    ##  tibble              3.3.0      2025-06-08 [1] CRAN (R 4.4.3)
    ##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.4.3)
    ##  timechange          0.3.0      2024-01-18 [1] CRAN (R 4.4.3)
    ##  UCSC.utils          1.2.0      2024-10-29 [1] Bioconduc~
    ##  vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.4.3)
    ##  withr               3.0.2      2024-10-28 [1] CRAN (R 4.4.3)
    ##  xfun                0.53       2025-08-19 [1] CRAN (R 4.4.3)
    ##  XVector           * 0.46.0     2024-10-29 [1] Bioconduc~
    ##  yaml                2.3.10     2024-07-26 [1] CRAN (R 4.4.3)
    ##  zip                 2.3.3      2025-05-13 [1] CRAN (R 4.4.3)
    ##  zlibbioc            1.52.0     2024-10-29 [1] Bioconduc~
    ## 
    ##  [1] /home/ysugimo/R/x86_64-pc-linux-gnu-library/4.4
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
