p2-07 · Plot site-specific stoichiometry
================
Yoichiro Sugimoto and Pallavi Kesavan
18 June, 2026

- [Overview](#overview)
- [Setup](#setup)
- [Helper functions](#helper-functions)
- [Load and preprocess data](#load-and-preprocess-data)
- [Normoxia time course](#normoxia-time-course)
- [Normoxia vs hypoxia](#normoxia-vs-hypoxia)
- [Session information](#session-information)

# Overview

**Purpose:** Visualise site-specific hydroxylation stoichiometry in
BRD2/3/4 across O2% and dox induction (via `plot_ptm_stoichiometry`).

**Inputs:** MS_KR_1 / MS_SS / PNAS stoichiometry tables from p2-01;
reference proteome.

**Outputs:** figures (rendered on knit).

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

# Helper functions

``` r
# Order the sample factor and draw the per-site stoichiometry plot for one accession
plot_brd_stoic <- function(subset_dt, accession, plot_range, sample_levels) {
  subset_dt[, sample_name := factor(sample_name, levels = sample_levels)]
  plot_ptm_stoichiometry(
    subset_dt,
    accession = accession, plot_range = plot_range,
    all.protein.bs, sample_colors = NA
  )
}
```

# Load and preprocess data

Loads and preprocesses the MS_KR_1 / MS_SS / PNAS stoichiometry tables
(via load_stoichiometry_datasets()).

``` r
# Load and preprocess the MS_KR_1, MS_SS and PNAS stoichiometry tables
# (see load_stoichiometry_datasets() in R/functions/2-useful_functions.R).
.stoic_data       <- load_stoichiometry_datasets(results.dir, data.dir)
ms_kr1_stoic_dt   <- .stoic_data$MS_KR1_stoic_dt
ms_ss_stoic_dt    <- .stoic_data$MS_SS_stoic_dt
pnas2022_stoic_dt <- .stoic_data$pnas2022_stoic_dt
```

``` r
# Remove BRD2/3 samples reporting BRD4 stoichiometry and vice versa: these are
# cross-protein false positives (e.g. from shared / homologous peptides).
ms_ss_stoic_dt <- ms_ss_stoic_dt[
  !((grepl("_BRD23$", sample_name) & gene_name == "BRD4") |
       (grepl("_BRD4$", sample_name) & gene_name %in% c("BRD2", "BRD3")))
]

# Combine the three datasets and keep BRD2/3/4 only
brd_stoic_dt <- rbindlist(
  list(ms_ss_stoic_dt, ms_kr1_stoic_dt, pnas2022_stoic_dt),
  use.names = TRUE
)
brd_stoic_dt <- brd_stoic_dt[gene_name %in% c("BRD2", "BRD3", "BRD4")]

# Simplified sample-group label shared by the plots below
brd_stoic_dt[, sample_group := case_when(
  condition == "MS_SS" & grepl("minusDox", sample_name) ~ "iJ6_0h_21pc_SS",
  condition == "MS_SS" ~ paste0("iJ6_", str_split_fixed(sample_name, "_", 3)[, 1], "_", str_split_fixed(sample_name, "_", 3)[, 2], "_SS"),
  sample_name == "HeLaiJMJD6_noDox_N_NA"       ~ "iJ6_0h_21pc_KR",
  sample_name == "HeLaiJMJD6_Dox_N_NA"         ~ "iJ6_24h_21pc_KR",
  sample_name == "HeLaWT_NA_N_NA"              ~ "WT_Inf_21pc_KR",
  sample_name == "HeLaiJMJD6_Dox_01O224h_NA"   ~ "iJ6_24h_01pc_KR",
  sample_name == "JQ1_HeLaWT_derivatised"      ~ "WT_Inf_21pc_PNAS",
  sample_name == "JQ1_HeLaJMJD6KO_derivatised" ~ "iJ6_Inf_21pc_PNAS"
)]

brd_stoic_dt[grepl("iJ6_0h", sample_group), .N, sample_group]
```

    ##      sample_group     N
    ##            <char> <int>
    ## 1: iJ6_0h_21pc_SS  1706
    ## 2: iJ6_0h_21pc_KR  1938

# Normoxia time course

Plots site-specific stoichiometry under normoxia across the JMJD6
re-expression time course for BRD2/3/4. The plotted amino-acid range
matches the region shown in the PNAS2022 paper.

``` r
# Select the appropriate minus-Dox (un-induced) control and exclude the JQ1/PNAS
# derivatised samples from these plots.
brd2_21pc_dt <- brd_stoic_dt[
  sample_group %in% c(
    "WT_Inf_21pc_KR", "iJ6_0h_21pc_KR", "iJ6_4h_21pc_SS", "iJ6_8h_21pc_SS", "iJ6_18h_21pc_SS", "iJ6_24h_21pc_KR"
    
  ) &
    gene_name == "BRD2"
]

plot_brd_stoic(
  brd2_21pc_dt, accession = "P25440", plot_range = c(540, 590),
  sample_levels = c("HeLaWT_NA_N_NA", "HeLaiJMJD6_noDox_N_NA",
                    "4h_21pc_BRD23", "8h_21pc_BRD23", "18h_21pc_BRD23",
                    "HeLaiJMJD6_Dox_N_NA")
)
```

![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd2_normoxia-1.png)<!-- -->![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd2_normoxia-2.png)<!-- -->

``` r
brd3_21pc_dt <- brd_stoic_dt[
  grepl("Inf|0h|4h|8h|18h|24h", sample_group) &
    grepl("21pc", sample_group) &
    !sample_name %in% c("minusDox_BRD23", "JQ1_HeLaWT_derivatised", "JQ1_HeLaJMJD6KO_derivatised") &
    gene_name == "BRD3"
]

plot_brd_stoic(
  brd3_21pc_dt, accession = "Q15059", plot_range = c(483, 533),
  sample_levels = c("HeLaWT_NA_N_NA", "HeLaiJMJD6_noDox_N_NA",
                    "4h_21pc_BRD23", "8h_21pc_BRD23", "18h_21pc_BRD23",
                    "HeLaiJMJD6_Dox_N_NA")
)
```

![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd3_normoxia-1.png)<!-- -->![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd3_normoxia-2.png)<!-- -->

``` r
brd4_21pc_dt <- brd_stoic_dt[
  grepl("Inf|0h|4h|8h|18h|24h", sample_group) &
    grepl("21pc", sample_group) &
    !sample_name %in% c("minusDox_BRD4", "JQ1_HeLaWT_derivatised", "JQ1_HeLaJMJD6KO_derivatised") &
    gene_name == "BRD4"
]

plot_brd_stoic(
  brd4_21pc_dt, accession = "O60885", plot_range = c(531, 581),
  sample_levels = c("HeLaWT_NA_N_NA", "HeLaiJMJD6_noDox_N_NA",
                    "4h_21pc_BRD4", "8h_21pc_BRD4", "18h_21pc_BRD4",
                    "HeLaiJMJD6_Dox_N_NA")
)
```

![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd4_normoxia-1.png)<!-- -->![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd4_normoxia-2.png)<!-- -->

# Normoxia vs hypoxia

Plots site-specific stoichiometry comparing normoxia and hypoxia for
BRD2/3/4.

Sample-group ↔ sample-name key:

    # HeLaWT_NA_N_NA            - WT_Inf_21pc_KR
    # JQ1_HeLaWT_derivatised    - WT_Inf_21pc_PNAS
    # JQ1_HeLaJMJD6KO_derivatised - iJ6_Inf_21pc_PNAS
    # HeLaiJMJD6_Dox_N_NA       - iJ6_24h_21pc_KR
    # 18h_4pc_BRD23             - iJ6_18h_4pc_SS
    # 18h_1pc_BRD23             - iJ6_18h_1pc_SS
    # HeLaiJMJD6_Dox_01O224h_NA - iJ6_24h_01pc_KR
    # minusDox_BRD23            - iJ6_0h_21pc_SS
    # HeLaiJMJD6_noDox_N_NA     - iJ6_0h_21pc_KR

``` r
brd2_o2_dt <- brd_stoic_dt[
  grepl("Inf|18h|24h|0h", sample_group) &
    !sample_name %in% c("minusDox_BRD23", "18h_21pc_BRD23", "HeLaiJMJD6_noDox_N_NA",
                        "JQ1_HeLaWT_derivatised", "JQ1_HeLaJMJD6KO_derivatised") &
    gene_name == "BRD2"
]

plot_brd_stoic(
  brd2_o2_dt, accession = "P25440", plot_range = c(540, 590),
  sample_levels = c("HeLaWT_NA_N_NA", "HeLaiJMJD6_Dox_N_NA",
                    "18h_4pc_BRD23", "18h_1pc_BRD23", "HeLaiJMJD6_Dox_01O224h_NA")
)
```

![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd2_hypoxia-1.png)<!-- -->![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd2_hypoxia-2.png)<!-- -->

``` r
brd3_o2_dt <- brd_stoic_dt[
  grepl("Inf|18h|24h|0h", sample_group) &
    !sample_name %in% c("minusDox_BRD23", "18h_21pc_BRD23", "HeLaiJMJD6_noDox_N_NA",
                        "JQ1_HeLaWT_derivatised", "JQ1_HeLaJMJD6KO_derivatised") &
    gene_name == "BRD3"
]

plot_brd_stoic(
  brd3_o2_dt, accession = "Q15059", plot_range = c(483, 533),
  sample_levels = c("HeLaWT_NA_N_NA", "HeLaiJMJD6_Dox_N_NA",
                    "18h_4pc_BRD23", "18h_1pc_BRD23", "HeLaiJMJD6_Dox_01O224h_NA")
)
```

![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd3_hypoxia-1.png)<!-- -->![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd3_hypoxia-2.png)<!-- -->

``` r
brd4_o2_dt <- brd_stoic_dt[
  grepl("Inf|18h|24h|0h", sample_group) &
    !sample_name %in% c("minusDox_BRD4", "18h_21pc_BRD4", "HeLaiJMJD6_noDox_N_NA",
                        "JQ1_HeLaWT_derivatised", "JQ1_HeLaJMJD6KO_derivatised") &
    gene_name == "BRD4"
]

plot_brd_stoic(
  brd4_o2_dt, accession = "O60885", plot_range = c(531, 581),
  sample_levels = c("HeLaWT_NA_N_NA", "HeLaiJMJD6_Dox_N_NA",
                    "18h_4pc_BRD4", "18h_1pc_BRD4", "HeLaiJMJD6_Dox_01O224h_NA")
)
```

![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd4_hypoxia-1.png)<!-- -->![](p2-07_plot_stoichiometry_files/figure-gfm/plot_brd4_hypoxia-2.png)<!-- -->

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
    ##  janitor             2.2.1      2024-12-22 [1] CRAN (R 4.4.3)
    ##  jsonlite            2.0.0      2025-03-27 [1] CRAN (R 4.4.3)
    ##  khroma            * 1.16.0     2025-02-25 [1] CRAN (R 4.4.3)
    ##  knitr             * 1.50       2025-03-16 [1] CRAN (R 4.4.3)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.4.3)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.4.3)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.4.3)
    ##  patchwork           1.3.2      2025-08-25 [1] CRAN (R 4.4.3)
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
