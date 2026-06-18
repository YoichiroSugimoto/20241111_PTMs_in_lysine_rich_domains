p2-01 · Compute PTM stoichiometry
================
Yoichiro Sugimoto and Pallavi Kesavan
18 June, 2026

- [Overview](#overview)
- [Setup](#setup)
- [Helper functions and shared
  settings](#helper-functions-and-shared-settings)
- [Stoichiometry without diagnostic
  ions](#stoichiometry-without-diagnostic-ions)
- [Stoichiometry with diagnostic ions —
  PNAS](#stoichiometry-with-diagnostic-ions--pnas)
- [Stoichiometry with diagnostic ions —
  MS_KR_1](#stoichiometry-with-diagnostic-ions--ms_kr_1)
- [Stoichiometry with diagnostic ions —
  MS_SS](#stoichiometry-with-diagnostic-ions--ms_ss)
- [Session information](#session-information)

# Overview

**Purpose:** Compute site-specific lysine-hydroxylation (PTM)
stoichiometry from MaxQuant outputs for every dataset in the study. This
is the single computation step of Part 2; all downstream scripts read
the tables written here.

**Inputs:**

- MaxQuant `evidence.txt` per dataset under `data/MQ_Std/` and
  `data/MQ_DI/`
- sample matrix `data/analysis_setting/PXD031221_sample_matrix.xlsx` and
  per-dataset `sample_info`
- PTM mapping
  `data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv`
- UniProt reference proteome (`reference_fasta`; see Setup)

**Outputs:** stoichiometry tables under `results/p2-analysis-setting/`
(`MQ_Std_MSMS_SECPEP/`, `MQ_Std_MSMS/`, `MQ_DI/`, `MS_KR_1/`, `MS_SS/`).

**Upstream:** none. **Downstream:** p2-02 … p2-09 read these tables.

# Setup

``` r
## Resolve the repository root (via the .here sentinel) and load the shared
## setup: packages, helper functions, ggplot/knitr settings, and project paths
## (project.dir, data.dir, results.dir, reference_fasta).
repo_root <- local({
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  while (!file.exists(file.path(p, ".here")) && !identical(dirname(p), p)) p <- dirname(p)
  p
})
source(file.path(repo_root, "R", "functions", "_setup.R"))
```

# Helper functions and shared settings

Defines the two stoichiometry wrappers used throughout this script and
the common PTM-mapping file. Every dataset section below routes through
one of these wrappers:

- `process_stoichiometry()` — standard MaxQuant search (no diagnostic
  ions); loops over a sample matrix and calls
  `calculate_stoichiometry()`.
- `process_stoichiometry2()` — diagnostic-ion search; processes a single
  dataset (locating its per-PTM site files) and calls
  `calculate_stoichiometry2()`.

``` r
# Common PTM mapping file, shared by all datasets
ptm_mapping_file <- file.path(
  project.dir,
  "data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv"
)

# Automate stoichiometry calculations across standard-search MaxQuant datasets
# (e.g. MULTI-MSMS + MULTI-SECPEP and MULTI-MSMS only). `sample_dt` is the sample
# matrix (one row per run); `data_dir` is the MaxQuant output root for the dataset.
process_stoichiometry <- function(sample_dt,
                                  data_dir,
                                  out_dir,
                                  ref_protein_dt,
                                  ptm_mapping_file,
                                  selected_type = NA,
                                  prefix_label = "") {

  # Create output directory (and parents) if it does not already exist
  create.dir(out_dir)

  # For each run in the sample matrix, calculate the stoichiometry
  for (i in seq_len(nrow(sample_dt))) {
    prefix     <- sample_dt$prefix[i]                          # run prefix
    sample     <- sample_dt$data[i]                            # dataset folder
    mq_setting <- gsub("_$", "", sub("^data-[^_]+_", "", prefix))  # search setting

    # Path to this run's MaxQuant evidence file
    mq_evidence_file <- file.path(
      data_dir, sample, mq_setting, paste0(prefix, "evidence.txt")
    )

    message("Processing: ", basename(mq_evidence_file))

    # Proceed only if the evidence file exists, otherwise skip with a message
    if (file.exists(mq_evidence_file)) {

      # Corresponding sample information (metadata) file
      sample_info_file <- file.path(
        data_dir, sample,
        paste0("MS_dataset_overview_PXD031221_", sample, ".csv")
      )

      # Run main stoichiometry calculation
      calculate_stoichiometry(
        mq_evidence_data = mq_evidence_file,
        sample_info_file = sample_info_file,
        ref_protein_dt   = ref_protein_dt,
        output_prefix    = file.path(out_dir, paste0(prefix_label, prefix)),
        ptm_mapping_file = ptm_mapping_file,
        K_only           = FALSE,
        selected_type    = selected_type
      )

    } else {
      message("File does not exist: ", mq_evidence_file)
    }
    gc()  # free memory between runs
  }
}

# Automate diagnostic-ion stoichiometry calculations for a single dataset.
# Counterpart to process_stoichiometry(), but wraps calculate_stoichiometry2()
# (which additionally needs the per-PTM site files in `ptm_dir`).
process_stoichiometry2 <- function(mq_evidence_file,
                                   sample_info_file,
                                   ptm_dir,
                                   output_prefix,
                                   ref_protein_dt,
                                   ptm_mapping_file,
                                   selected_type = "MULTI-MSMS",
                                   parse_protein_accession_function = NA) {

  message("Processing: ", basename(mq_evidence_file))

  # Proceed only if the evidence file exists, otherwise skip with a message
  if (file.exists(mq_evidence_file)) {

    # Locate all PTM site files within the ptm folder
    ptm_files <- list.files(ptm_dir, full.names = TRUE)

    # Generate PTM names from file names: add brackets, drop "Sites.txt"
    ptm_names <- paste0(
      "[", str_replace_all(basename(ptm_files), "Sites.txt", ""), "]"
    )

    # Run the stoichiometry calculation
    calculate_stoichiometry2(
      mq_evidence_data = mq_evidence_file,
      sample_info_file = sample_info_file,
      ref_protein_dt   = ref_protein_dt,
      ptm_files        = ptm_files,
      ptm_names        = ptm_names,
      output_prefix    = output_prefix,
      ptm_mapping_file = ptm_mapping_file,
      K_only           = FALSE,
      selected_type    = selected_type,
      parse_protein_accession_function = parse_protein_accession_function
    )
  } else {
    message("File does not exist: ", mq_evidence_file)
  }
  gc()  # free memory
}
```

# Stoichiometry without diagnostic ions

Computes stoichiometry for the PNAS2022 standard MaxQuant search — once
including secondary-peptide identifications (MULTI-MSMS + MULTI-SECPEP)
and once with MULTI-MSMS only.

``` r
# MaxQuant output root for the standard (no diagnostic ion) PNAS2022 search
mq_std_data_dir <- file.path(project.dir, "data/MQ_Std/PNAS2022")

# Sample matrix (one row per run)
mq_std_sample_dt <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_Std"
) %>% data.table

# Sequential sample id (1:N)
mq_std_sample_dt[, sample_id := 1:.N]
```

``` r
## Stoichiometry analysis using MULTI-MSMS with MULTI-SECPEP (secondary peptide)
## and, separately, MULTI-MSMS only.

## --- MULTI-MSMS + MULTI-SECPEP ---------------------------------------------

# Restrict the standard search to data-A and data-D for the MULTI-MSMS +
# MULTI-SECPEP comparison.
# Data-A and data-D — BRD proteins with or without derivatisation — are analysed to assess the effect of derivatisation.
mq_std_secpep_sample_dt <- mq_std_sample_dt[data %in% c("data-A", "data-D")]

mq_std_msms_secpep_dir <- file.path(results.dir, "p2-analysis-setting", "MQ_Std_MSMS_SECPEP")

process_stoichiometry(
  sample_dt        = mq_std_secpep_sample_dt,
  data_dir         = mq_std_data_dir,
  out_dir          = mq_std_msms_secpep_dir,
  ref_protein_dt   = ref_protein_dt,
  ptm_mapping_file = ptm_mapping_file,
  selected_type    = NA,
  prefix_label     = "including_SECPEP_"
)
```

    ## Processing: data-A_trp_m2_v2_def_evidence.txt

    ## Processing: data-A_trp_m3_v3_def_evidence.txt

    ## Processing: data-A_trp_m4_v4_def_evidence.txt

    ## Processing: data-A_trp_m5_v5_def_evidence.txt

    ## Processing: data-A_trp_m6_v6_def_evidence.txt

    ## Processing: data-A_trp_m7_v7_def_evidence.txt

    ## Processing: data-A_trp_m8_v8_def_evidence.txt

    ## Processing: data-A_argC_m2_v7_def_evidence.txt

    ## Processing: data-D_trp_m2_v2_def_evidence.txt

    ## Processing: data-D_trp_m5_v5_def_evidence.txt

    ## Processing: data-D_trp_m7_v7_def_evidence.txt

``` r
## --- MULTI-MSMS only -------------------------------------------------------

mq_std_msms_dir <- file.path(results.dir, "p2-analysis-setting", "MQ_Std_MSMS")

process_stoichiometry(
  sample_dt        = mq_std_sample_dt,
  data_dir         = mq_std_data_dir,
  out_dir          = mq_std_msms_dir,
  ref_protein_dt   = ref_protein_dt,
  ptm_mapping_file = ptm_mapping_file,
  selected_type    = c("MULTI-MSMS"),
  prefix_label     = ""
)
```

    ## Processing: data-A_trp_m2_v2_def_evidence.txt

    ## Processing: data-A_trp_m3_v3_def_evidence.txt

    ## Processing: data-A_trp_m4_v4_def_evidence.txt

    ## Processing: data-A_trp_m5_v5_def_evidence.txt

    ## Processing: data-A_trp_m6_v6_def_evidence.txt

    ## Processing: data-A_trp_m7_v7_def_evidence.txt

    ## Processing: data-A_trp_m8_v8_def_evidence.txt

    ## Processing: data-A_argC_m2_v7_def_evidence.txt

    ## Processing: data-B_trp_m2_v2_mCC_evidence.txt

    ## Processing: data-B_trp_m5_v5_mCC_evidence.txt

    ## Processing: data-B_trp_m7_v7_mCC_evidence.txt

    ## Processing: data-C_trp_m2_v2_mCC_evidence.txt

    ## Processing: data-C_trp_m5_v5_mCC_evidence.txt

    ## Processing: data-C_trp_m7_v7_mCC_evidence.txt

    ## Processing: data-D_trp_m2_v2_def_evidence.txt

    ## Processing: data-D_trp_m5_v5_def_evidence.txt

    ## Processing: data-D_trp_m7_v7_def_evidence.txt

# Stoichiometry with diagnostic ions — PNAS

Computes stoichiometry for the diagnostic-ion (no water-loss) MaxQuant
search of the PNAS2022 datasets (A, B1, B2, C).

``` r
# MaxQuant output root for the diagnostic-ion PNAS2022 search
mq_di_data_dir <- file.path(project.dir, "data/MQ_DI/PNAS2022")

# Output directory for results
mq_di_dir <- file.path(results.dir, "p2-analysis-setting", "MQ_DI")
create.dir(mq_di_dir)

# Sample matrix (one row per run)
mq_di_sample_dt <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_DI"
) %>% data.table

# Sequential sample id (1:N)
mq_di_sample_dt[, sample_id := 1:.N]
```

``` r
# For each run in the sample matrix, calculate the stoichiometry
for (i in seq_len(nrow(mq_di_sample_dt))) {

  prefix     <- mq_di_sample_dt$prefix[i]
  mq_setting <- gsub("_$", "", sub("^data-[^_]+_", "", prefix))  # search setting

  # Per-run directory holding the evidence file and the ptm/ folder
  run_dir <- file.path(mq_di_data_dir, mq_di_sample_dt[i, data], mq_setting)

  process_stoichiometry2(
    mq_evidence_file = file.path(run_dir, paste0(prefix, "evidence.txt")),
    sample_info_file = file.path(
      mq_di_data_dir,
      mq_di_sample_dt[i, data],
      paste0("MS_dataset_overview_PXD031221_", mq_di_sample_dt[i, data], ".csv")
    ),
    ptm_dir          = file.path(run_dir, "ptm"),
    output_prefix    = file.path(mq_di_dir, paste0("DI_", prefix)),
    ref_protein_dt   = ref_protein_dt,
    ptm_mapping_file = ptm_mapping_file
  )
}
```

    ## Processing: data-A_trp_m7_v7_def_evidence.txt

    ## Processing: data-B1_trp_m7_v7_mCC_evidence.txt

    ## Processing: data-B2_trp_m7_v7_mCC_evidence.txt

    ## Processing: data-C_trp_m7_v7_mCC_evidence.txt

# Stoichiometry with diagnostic ions — MS_KR_1

Computes stoichiometry for the MS_KR_1 dataset (dataset D: hypoxia vs
normoxia, with JMJD6 re-expression).

``` r
# MaxQuant output directory for the MS_KR_1 dataset
ms_kr1_data_dir <- file.path(project.dir, "data/MQ_DI/MS_KR_1")

# Output directory for results
ms_kr1_dir <- file.path(results.dir, "p2-analysis-setting", "MS_KR_1")
create.dir(ms_kr1_dir)
```

``` r
process_stoichiometry2(
  mq_evidence_file = file.path(ms_kr1_data_dir, "MS_KR_1_evidence.txt"),
  sample_info_file = file.path(ms_kr1_data_dir, "sample_info.csv"),
  ptm_dir          = file.path(ms_kr1_data_dir, "ptm"),
  output_prefix    = file.path(ms_kr1_dir, "MS_KR_1_"),
  ref_protein_dt   = ref_protein_dt,
  ptm_mapping_file = ptm_mapping_file
)
```

    ## Processing: MS_KR_1_evidence.txt

    ##            used  (Mb) gc trigger  (Mb)  max used   (Mb)
    ## Ncells  4347231 232.2   13789037 736.5  13789037  736.5
    ## Vcells 11849723  90.5   57369272 437.7 187872791 1433.4

# Stoichiometry with diagnostic ions — MS_SS

Computes stoichiometry for the MS_SS dataset (dataset E: O2 and
doxycycline titration).

``` r
# MaxQuant output directory for the MS_SS dataset
ms_ss_data_dir <- file.path(project.dir, "data/MQ_DI/MS_SS")

# Output directory for results
ms_ss_dir <- file.path(results.dir, "p2-analysis-setting", "MS_SS")
create.dir(ms_ss_dir)
```

``` r
process_stoichiometry2(
  mq_evidence_file = file.path(ms_ss_data_dir, "MS_SS_evidence.txt"),
  sample_info_file = file.path(ms_ss_data_dir, "sample_info.csv"),
  ptm_dir          = file.path(ms_ss_data_dir, "ptm"),
  output_prefix    = file.path(ms_ss_dir, "MS_SS_"),
  ref_protein_dt   = ref_protein_dt,
  ptm_mapping_file = ptm_mapping_file
)
```

    ## Processing: MS_SS_evidence.txt

    ##            used  (Mb) gc trigger  (Mb)  max used   (Mb)
    ## Ncells  4348040 232.3   13789037 736.5  13789037  736.5
    ## Vcells 11764388  89.8   45895418 350.2 187872791 1433.4

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
    ##  bit                 4.6.0      2025-03-06 [1] CRAN (R 4.4.3)
    ##  bit64               4.6.0-1    2025-01-16 [1] CRAN (R 4.4.3)
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
    ##  janitor             2.2.1      2024-12-22 [1] CRAN (R 4.4.3)
    ##  jsonlite            2.0.0      2025-03-27 [1] CRAN (R 4.4.3)
    ##  khroma            * 1.16.0     2025-02-25 [1] CRAN (R 4.4.3)
    ##  knitr             * 1.50       2025-03-16 [1] CRAN (R 4.4.3)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.4.3)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.4.3)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.4.3)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.4.3)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.3)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-12-13 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.4.3)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.4.3)
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
    ##  zlibbioc            1.52.0     2024-10-29 [1] Bioconduc~
    ## 
    ##  [1] /home/ysugimo/R/x86_64-pc-linux-gnu-library/4.4
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
