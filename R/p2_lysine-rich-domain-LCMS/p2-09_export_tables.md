p2-09 · Export stoichiometry tables
================
Yoichiro Sugimoto and Pallavi Kesavan
24 August, 2026

- <a href="#overview" id="toc-overview">Overview</a>
- <a href="#setup" id="toc-setup">Setup</a>
- <a href="#export-tables-to-excel" id="toc-export-tables-to-excel">Export
  tables to Excel</a>
- <a href="#session-information" id="toc-session-information">Session
  information</a>

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
```

    ## - The project is out-of-sync -- use `renv::status()` for details.

``` r
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
    ##  version  R version 4.5.0 (2025-04-11)
    ##  os       Red Hat Enterprise Linux 9.6 (Plow)
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Europe/Berlin
    ##  date     2026-08-24
    ##  pandoc   2.19.2 @ /gnu/store/sqwwnsp5xb8yd3z1a57lhldcsvx3z9gb-profile/bin/ (via rmarkdown)
    ##  quarto   NA
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package           * version    date (UTC) lib source
    ##  P AnnotationDbi     * 1.72.0     2025-10-29 [?] Bioconduc~
    ##  P beeswarm            0.4.0      2021-06-01 [?] CRAN (R 4.5.0)
    ##  P Biobase           * 2.70.0     2025-10-29 [?] Bioconduc~
    ##  P BiocGenerics      * 0.56.0     2025-10-29 [?] Bioconduc~
    ##  P Biostrings        * 2.78.0     2025-10-29 [?] Bioconduc~
    ##  P bit                 4.6.0      2025-03-06 [?] CRAN (R 4.5.0)
    ##  P bit64               4.8.2      2026-05-19 [?] CRAN (R 4.5.0)
    ##  P blob                1.3.0      2026-01-14 [?] CRAN (R 4.5.0)
    ##  P cachem              1.1.0      2024-05-16 [?] CRAN (R 4.5.0)
    ##  P cellranger          1.1.0      2016-07-27 [?] CRAN (R 4.5.0)
    ##  P cli                 3.6.6      2026-04-09 [?] CRAN (R 4.5.0)
    ##  P colourpicker        1.3.0      2023-08-21 [?] CRAN (R 4.5.0)
    ##  P crayon              1.5.3      2024-06-20 [?] CRAN (R 4.5.0)
    ##  P data.table        * 1.18.4     2026-05-06 [?] CRAN (R 4.5.0)
    ##  P DBI                 1.3.0      2026-02-25 [?] CRAN (R 4.5.0)
    ##  P digest              0.6.39     2025-11-19 [?] CRAN (R 4.5.0)
    ##  P dplyr             * 1.2.1      2026-04-03 [?] CRAN (R 4.5.0)
    ##  P eulerr            * 7.0.2      2024-03-28 [?] CRAN (R 4.5.0)
    ##  P evaluate            1.0.5      2025-08-27 [?] CRAN (R 4.5.0)
    ##  P farver              2.1.2      2024-05-13 [?] CRAN (R 4.5.0)
    ##  P fastmap             1.2.0      2024-05-15 [?] CRAN (R 4.5.0)
    ##  P formattable         0.2.1      2021-01-07 [?] CRAN (R 4.5.0)
    ##  P generics          * 0.1.4      2025-05-09 [?] CRAN (R 4.5.0)
    ##  P ggbeeswarm        * 0.7.3      2025-11-29 [?] CRAN (R 4.5.0)
    ##  P ggplot2           * 4.0.3      2026-04-22 [?] CRAN (R 4.5.0)
    ##  P glue                1.8.1      2026-04-17 [?] CRAN (R 4.5.0)
    ##  P gridExtra           2.3        2017-09-09 [?] CRAN (R 4.5.0)
    ##  P gtable              0.3.6      2024-10-25 [?] CRAN (R 4.5.0)
    ##  P htmltools           0.5.9      2025-12-04 [?] CRAN (R 4.5.0)
    ##  P htmlwidgets         1.6.4      2023-12-06 [?] CRAN (R 4.5.0)
    ##  P httpuv              1.6.17     2026-03-18 [?] CRAN (R 4.5.0)
    ##  P httr                1.4.8      2026-02-13 [?] CRAN (R 4.5.0)
    ##  P IRanges           * 2.44.0     2025-10-29 [?] Bioconduc~
    ##  P janitor           * 2.2.1      2024-12-22 [?] CRAN (R 4.5.0)
    ##  P jsonlite            2.0.0      2025-03-27 [?] CRAN (R 4.5.0)
    ##  P KEGGREST            1.50.0     2025-10-29 [?] Bioconduc~
    ##  P khroma            * 1.17.0     2025-09-29 [?] CRAN (R 4.5.0)
    ##  P knitr             * 1.51       2025-12-20 [?] CRAN (R 4.5.0)
    ##  P labeling            0.4.3      2023-08-29 [?] CRAN (R 4.5.0)
    ##  P later               1.4.8      2026-03-05 [?] CRAN (R 4.5.0)
    ##  P lattice             0.22-9     2026-02-09 [?] CRAN (R 4.5.0)
    ##  P lazyeval            0.2.3      2026-04-04 [?] CRAN (R 4.5.0)
    ##  P lifecycle           1.0.5      2026-01-08 [?] CRAN (R 4.5.0)
    ##  P lubridate           1.9.5      2026-02-04 [?] CRAN (R 4.5.0)
    ##  P magrittr          * 2.0.5      2026-04-04 [?] CRAN (R 4.5.0)
    ##  P Matrix              1.7-5      2026-03-21 [?] CRAN (R 4.5.0)
    ##  P memoise             2.0.1      2021-11-26 [?] CRAN (R 4.5.0)
    ##  P mgcv              * 1.9-4      2025-11-07 [?] CRAN (R 4.5.0)
    ##  P mime                0.13       2025-03-17 [?] CRAN (R 4.5.0)
    ##  P miniUI              0.1.2      2025-04-17 [?] CRAN (R 4.5.0)
    ##  P nlme              * 3.1-169    2026-03-27 [?] CRAN (R 4.5.0)
    ##  P openxlsx          * 4.2.8.1    2025-10-31 [?] CRAN (R 4.5.0)
    ##  P org.Hs.eg.db      * 3.22.0     2026-06-24 [?] Bioconductor
    ##  P otel                0.2.0      2025-08-29 [?] CRAN (R 4.5.0)
    ##  P patchwork         * 1.3.2      2025-08-25 [?] CRAN (R 4.5.0)
    ##  P pillar              1.11.1     2025-09-17 [?] CRAN (R 4.5.0)
    ##  P pkgconfig           2.0.3      2019-09-22 [?] CRAN (R 4.5.0)
    ##  P plotly              4.12.0     2026-01-24 [?] CRAN (R 4.5.0)
    ##  P plyr                1.8.9      2023-10-02 [?] CRAN (R 4.5.0)
    ##  P png                 0.1-9      2026-03-15 [?] CRAN (R 4.5.0)
    ##  P polyclip            1.10-7     2024-07-23 [?] CRAN (R 4.5.0)
    ##  P polylabelr          1.0.0      2026-01-19 [?] CRAN (R 4.5.0)
    ##  P promises            1.5.0      2025-11-01 [?] CRAN (R 4.5.0)
    ##    ptm.stoichiometry * 0.0.0.9000 2026-06-24 [1] local (/fast/AG_Sugimoto/home/users/yoichiro/projects/ptm.stoichiometry)
    ##  P purrr               1.2.2      2026-04-10 [?] CRAN (R 4.5.0)
    ##  P R6                  2.6.1      2025-02-15 [?] CRAN (R 4.5.0)
    ##  P RColorBrewer      * 1.1-3      2022-04-03 [?] CRAN (R 4.5.0)
    ##  P Rcpp                1.1.1-1.1  2026-04-24 [?] CRAN (R 4.5.0)
    ##  P readxl            * 1.5.0      2026-05-16 [?] CRAN (R 4.5.0)
    ##    renv                1.1.5      2025-07-24 [1] CRAN (R 4.5.0)
    ##  P rlang               1.2.0      2026-04-06 [?] CRAN (R 4.5.0)
    ##  P rmarkdown           2.31       2026-03-26 [?] CRAN (R 4.5.0)
    ##  P RSQLite             3.53.2     2026-06-17 [?] CRAN (R 4.5.0)
    ##  P S4Vectors         * 0.48.1     2026-04-05 [?] Bioconduc~
    ##  P S7                  0.2.2      2026-04-22 [?] CRAN (R 4.5.0)
    ##  P scales              1.4.0      2025-04-24 [?] CRAN (R 4.5.0)
    ##  P Seqinfo           * 1.0.0      2025-10-29 [?] Bioconduc~
    ##  P sessioninfo         1.2.4      2026-06-04 [?] CRAN (R 4.5.0)
    ##  P shiny               1.14.0     2026-06-21 [?] CRAN (R 4.5.0)
    ##  P shinythemes         1.2.0      2021-01-25 [?] CRAN (R 4.5.0)
    ##  P snakecase           0.11.1     2023-08-27 [?] CRAN (R 4.5.0)
    ##  P stringi             1.8.7      2025-03-27 [?] CRAN (R 4.5.0)
    ##  P stringr           * 1.6.0      2025-11-04 [?] CRAN (R 4.5.0)
    ##    subcellularvis    * 0.0.0.9000 2026-06-24 [1] local (/fast/AG_Sugimoto/home/users/yoichiro/software/R_packages/subcellularvis)
    ##  P tibble              3.3.1      2026-01-11 [?] CRAN (R 4.5.0)
    ##  P tidyr               1.3.2      2025-12-19 [?] CRAN (R 4.5.0)
    ##  P tidyselect          1.2.1      2024-03-11 [?] CRAN (R 4.5.0)
    ##  P timechange          0.4.0      2026-01-29 [?] CRAN (R 4.5.0)
    ##  P UpSetR              1.4.1      2026-05-25 [?] CRAN (R 4.5.0)
    ##  P vctrs               0.7.3      2026-04-11 [?] CRAN (R 4.5.0)
    ##  P vipor               0.4.7      2023-12-18 [?] CRAN (R 4.5.0)
    ##  P viridisLite         0.4.3      2026-02-04 [?] CRAN (R 4.5.0)
    ##  P withr               3.0.3      2026-06-19 [?] CRAN (R 4.5.0)
    ##  P xfun                0.59       2026-06-19 [?] CRAN (R 4.5.0)
    ##  P xtable              1.8-8      2026-02-22 [?] CRAN (R 4.5.0)
    ##  P XVector           * 0.50.0     2025-10-29 [?] Bioconduc~
    ##  P yaml                2.3.12     2025-12-10 [?] CRAN (R 4.5.0)
    ##  P zip                 3.0.0      2026-06-10 [?] CRAN (R 4.5.0)
    ## 
    ##  [1] /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/renv/library/linux-rhel-9.6/R-4.5/x86_64-unknown-linux-gnu
    ##  [2] /fast/home/y/ysugimo/.cache/R/renv/sandbox/linux-rhel-9.6/R-4.5/x86_64-unknown-linux-gnu/cb72a45c
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
