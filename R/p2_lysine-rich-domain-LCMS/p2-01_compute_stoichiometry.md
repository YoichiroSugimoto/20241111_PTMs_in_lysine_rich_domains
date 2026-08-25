p2-01 · Compute PTM stoichiometry
================
Yoichiro Sugimoto and Pallavi Kesavan
25 August, 2026

- <a href="#overview" id="toc-overview">Overview</a>
- <a href="#setup" id="toc-setup">Setup</a>
- <a href="#helper-functions-and-shared-settings"
  id="toc-helper-functions-and-shared-settings">Helper functions and
  shared settings</a>
- <a href="#stoichiometry-without-diagnostic-ions"
  id="toc-stoichiometry-without-diagnostic-ions">Stoichiometry without
  diagnostic ions</a>
- <a href="#stoichiometry-with-diagnostic-ions--pnas"
  id="toc-stoichiometry-with-diagnostic-ions--pnas">Stoichiometry with
  diagnostic ions — PNAS</a>
- <a href="#stoichiometry-with-diagnostic-ions--ms_kr_1"
  id="toc-stoichiometry-with-diagnostic-ions--ms_kr_1">Stoichiometry with
  diagnostic ions — MS_KR_1</a>
- <a href="#stoichiometry-with-diagnostic-ions--ms_ss"
  id="toc-stoichiometry-with-diagnostic-ions--ms_ss">Stoichiometry with
  diagnostic ions — MS_SS</a>
- <a href="#session-information" id="toc-session-information">Session
  information</a>

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

    ## - The project is out-of-sync -- use `renv::status()` for details.

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

    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  19127763 1021.6   33303302 1778.6  33303302 1778.6
    ## Vcells 210777540 1608.2  498108776 3800.3 498108773 3800.3

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

    ##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
    ## Ncells  19128537 1021.6   33303302 1778.6  33303302 1778.6
    ## Vcells 210778799 1608.2  498108776 3800.3 498108773 3800.3

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
    ##  date     2026-08-25
    ##  pandoc   2.19.2 @ /gnu/store/sqwwnsp5xb8yd3z1a57lhldcsvx3z9gb-profile/bin/ (via rmarkdown)
    ##  quarto   NA
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package           * version    date (UTC) lib source
    ##  P AnnotationDbi     * 1.72.0     2025-10-29 [?] Bioconduc~
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
    ##  P evaluate            1.0.5      2025-08-27 [?] CRAN (R 4.5.0)
    ##  P farver              2.1.2      2024-05-13 [?] CRAN (R 4.5.0)
    ##  P fastmap             1.2.0      2024-05-15 [?] CRAN (R 4.5.0)
    ##  P formattable         0.2.1      2021-01-07 [?] CRAN (R 4.5.0)
    ##  P generics          * 0.1.4      2025-05-09 [?] CRAN (R 4.5.0)
    ##  P ggplot2           * 4.0.3      2026-04-22 [?] CRAN (R 4.5.0)
    ##  P glue                1.8.1      2026-04-17 [?] CRAN (R 4.5.0)
    ##  P gridExtra           2.3        2017-09-09 [?] CRAN (R 4.5.0)
    ##  P gtable              0.3.6      2024-10-25 [?] CRAN (R 4.5.0)
    ##  P htmltools           0.5.9      2025-12-04 [?] CRAN (R 4.5.0)
    ##  P htmlwidgets         1.6.4      2023-12-06 [?] CRAN (R 4.5.0)
    ##  P httpuv              1.6.17     2026-03-18 [?] CRAN (R 4.5.0)
    ##  P httr                1.4.8      2026-02-13 [?] CRAN (R 4.5.0)
    ##  P IRanges           * 2.44.0     2025-10-29 [?] Bioconduc~
    ##  P janitor             2.2.1      2024-12-22 [?] CRAN (R 4.5.0)
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
    ##  P org.Hs.eg.db      * 3.22.0     2026-06-24 [?] Bioconductor
    ##  P otel                0.2.0      2025-08-29 [?] CRAN (R 4.5.0)
    ##  P patchwork         * 1.3.2      2025-08-25 [?] CRAN (R 4.5.0)
    ##  P pillar              1.11.1     2025-09-17 [?] CRAN (R 4.5.0)
    ##  P pkgconfig           2.0.3      2019-09-22 [?] CRAN (R 4.5.0)
    ##  P plotly              4.12.0     2026-01-24 [?] CRAN (R 4.5.0)
    ##  P plyr                1.8.9      2023-10-02 [?] CRAN (R 4.5.0)
    ##  P png                 0.1-9      2026-03-15 [?] CRAN (R 4.5.0)
    ##  P promises            1.5.0      2025-11-01 [?] CRAN (R 4.5.0)
    ##    ptm.stoichiometry * 0.0.0.9000 2026-06-24 [1] local (/fast/AG_Sugimoto/home/users/yoichiro/projects/ptm.stoichiometry)
    ##  P purrr               1.2.2      2026-04-10 [?] CRAN (R 4.5.0)
    ##  P R6                  2.6.1      2025-02-15 [?] CRAN (R 4.5.0)
    ##  P RColorBrewer        1.1-3      2022-04-03 [?] CRAN (R 4.5.0)
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
    ##  P viridisLite         0.4.3      2026-02-04 [?] CRAN (R 4.5.0)
    ##  P withr               3.0.3      2026-06-19 [?] CRAN (R 4.5.0)
    ##  P xfun                0.59       2026-06-19 [?] CRAN (R 4.5.0)
    ##  P xtable              1.8-8      2026-02-22 [?] CRAN (R 4.5.0)
    ##  P XVector           * 0.50.0     2025-10-29 [?] Bioconduc~
    ##  P yaml                2.3.12     2025-12-10 [?] CRAN (R 4.5.0)
    ## 
    ##  [1] /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/renv/library/linux-rhel-9.6/R-4.5/x86_64-unknown-linux-gnu
    ##  [2] /fast/home/y/ysugimo/.cache/R/renv/sandbox/linux-rhel-9.6/R-4.5/x86_64-unknown-linux-gnu/cb72a45c
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
