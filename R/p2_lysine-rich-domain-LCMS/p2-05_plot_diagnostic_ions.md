p2-05 · Plot diagnostic-ion analysis
================
Yoichiro Sugimoto and Pallavi Kesavan
18 June, 2026

- [Overview](#overview)
- [Setup](#setup)
- [Helper functions](#helper-functions)
- [Import data](#import-data)
- [Diagnostic ions and precision](#diagnostic-ions-and-precision)
- [WT vs KO site overlap](#wt-vs-ko-site-overlap)
- [Comparison with PNAS 2022](#comparison-with-pnas-2022)
- [Session information](#session-information)

# Overview

**Purpose:** Evaluate how diagnostic ions improve hydroxylation analysis
— precision, WT / KO overlap, and comparison with PNAS 2022.

**Inputs:** stoichiometry tables from p2-01; PNAS2022 reference data.

**Outputs:** figures (rendered on knit).

**Upstream:** p2-01, p2-04. **Downstream:** none.

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

## Script-specific packages
library("patchwork")
library("eulerr")
library("RColorBrewer")
```

# Helper functions

``` r
# Area-proportional Euler diagram of WT vs KO site overlap
plot_overlap_venn <- function(wt_ids, ko_ids, main) {
  venn_list <- list(WT_DI_site = wt_ids, KO_DI_site = ko_ids)
  plot(
    euler(venn_list),
    fills      = list(fill = brewer.pal(3, "Set2"), alpha = 0.4),
    legend     = list(side = "right"),
    quantities = TRUE,
    main       = main
  )
}
```

# Import data

Loads the diagnostic-ion sample matrix (data paths come from
`_setup.R`).

``` r
# Sample run info (MQ diagnostic-ion search)
mq_di_sample_dt <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_DI"
) %>% data.table

# Sequential sample id (1:N)
mq_di_sample_dt[, sample_id := 1:.N]
```

# Diagnostic ions and precision

Assesses how requiring a diagnostic ion changes the precision of
hydroxylation calls.

``` r
# Read stoichiometry data for MQ_DI (change from the older version: a potential
# diagnostic ion with water loss is no longer considered)
mq_di_stoic_dt <- lapply(
  mq_di_sample_dt[, prefix],
  read_stoic_data,
  pre_prefix = "DI_",
  post_fix   = "_DI",
  dir_path   = file.path(results.dir, "p2-analysis-setting", "MQ_DI")
) %>% rbindlist

# fwrite(mq_di_stoic_dt[ptm != ""], file = file.path(results.dir, "p2-analysis-setting", "all_hydroxylysine.csv"))

mq_di_stoic_dt <- mq_di_stoic_dt[aa == "K"]
mq_di_stoic_dt[, genotype :=
  str_split_fixed(sample_name, "_", n = 3)[, 2] %>% factor(levels = c("HeLaWT", "HeLaJMJD6KO"))
]

# Sites with a diagnostic ion, and a per-row flag for them
di_sites_dt <- mq_di_stoic_dt[diagnostic_peak == "+"][
  !duplicated(paste(protein_accession, aa_pos)), .(protein_accession, aa_pos)
]
mq_di_stoic_dt[, DI_site :=
  paste(protein_accession, aa_pos) %in% di_sites_dt[, paste(protein_accession, aa_pos)]
]

# Protein feature data (lysines only), with a methionine-proximity annotation
protein_feature_dt <- fread(
  file.path(data.dir, "processed_data_from_PNAS2022/all_protein_feature_per_position.csv")
)
protein_feature_dt <- protein_feature_dt[residue == "K"]
setnames(protein_feature_dt, old = c("uniprot_id", "position"), new = c("protein_accession", "aa_pos"))

protein_feature_dt[, met_within_2 := case_when(
  nchar(Window) == 11 & grepl("M", substr(Window, start = 4, stop = 8)) ~ "Yes",  # M in the middle
  nchar(Window) == 11 ~ "No",                                                     # no M in window
  TRUE ~ "edge"                                                                   # M present but off-centre
)]

# Attach feature annotation to the stoichiometry data
wt_vs_ko_stoic_dt <- merge(
  copy(mq_di_stoic_dt),
  protein_feature_dt,
  by = c("protein_accession", "aa_pos")
)

# Keep only sites covered (>= 1 PSM) in both WT and KO so presence/absence is comparable
wt_coverage_dt <- wt_vs_ko_stoic_dt[genotype == "HeLaWT"][
  , list(total_wt_sum_psm_mapped = sum(sum_psm_mapped)), by = list(protein_accession, aa_pos)]
ko_coverage_dt <- wt_vs_ko_stoic_dt[genotype == "HeLaJMJD6KO"][
  , list(total_ko_sum_psm_mapped = sum(sum_psm_mapped)), by = list(protein_accession, aa_pos)]

all_coverage_dt <- merge(wt_coverage_dt, ko_coverage_dt, by = c("protein_accession", "aa_pos"))

wt_vs_ko_stoic_dt <- wt_vs_ko_stoic_dt[
  paste(protein_accession, aa_pos) %in%
    all_coverage_dt[total_wt_sum_psm_mapped > 0 & total_ko_sum_psm_mapped > 0][
      , paste(protein_accession, aa_pos)]
]

# Number of sites analysed, per genotype
wt_vs_ko_stoic_dt[!duplicated(paste(genotype, protein_accession, aa_pos))][, .N, by = genotype]
```

    ##       genotype     N
    ##         <fctr> <int>
    ## 1:      HeLaWT 26074
    ## 2: HeLaJMJD6KO 26074

``` r
# Hydroxylated K rows only, one (DI-preferred) row per genotype/site
koh_wt_vs_ko_dt <- wt_vs_ko_stoic_dt[ptm == "[Oxidation (K)]"]
koh_wt_vs_ko_dt <- koh_wt_vs_ko_dt[order(-DI_site)][
  !duplicated(paste(genotype, protein_accession, aa_pos))
]

g1 <- ggplot(
  koh_wt_vs_ko_dt[met_within_2 != "edge"],
  aes(x = genotype, fill = met_within_2)
) +
  geom_bar() +
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g2 <- ggplot(
  koh_wt_vs_ko_dt[met_within_2 != "edge" & DI_site == TRUE],
  aes(x = genotype, fill = met_within_2)
) +
  geom_bar() +
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

# Number of sites for plotting g2
koh_wt_vs_ko_dt[met_within_2 != "edge" & DI_site == TRUE][, .N, by = list(met_within_2, genotype)]
```

    ##    met_within_2    genotype     N
    ##          <char>      <fctr> <int>
    ## 1:           No HeLaJMJD6KO    20
    ## 2:           No      HeLaWT   143
    ## 3:          Yes HeLaJMJD6KO     2
    ## 4:          Yes      HeLaWT     5

``` r
g1 + g2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/precision-1.png)<!-- -->

# WT vs KO site overlap

Compares hydroxylysine sites detected in JMJD6 WT vs KO cells (overlap
and methionine proximity).

``` r
koh_per_site_dt <- dcast(
  koh_wt_vs_ko_dt,
  protein_accession + aa_pos + DI_site + met_within_2 ~ genotype,
  fun.aggregate = length
)

koh_per_site_dt[, Hyl_found_in := case_when(
  HeLaWT == 1 & HeLaJMJD6KO == 1 ~ "both",
  HeLaWT == 1 ~ "WT",
  HeLaJMJD6KO == 1 ~ "JMJD6KO",
  TRUE ~ "not_found"
) %>% factor(levels = c("WT", "JMJD6KO", "both", "not_found"))]

koh_per_site_dt <- koh_per_site_dt[met_within_2 != "edge"]

koh_per_site_count_dt <- koh_per_site_dt[, .N, by = list(Hyl_found_in, met_within_2, DI_site)]

g1 <- ggplot(koh_per_site_count_dt, aes(x = Hyl_found_in, y = N, fill = met_within_2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g2 <- ggplot(koh_per_site_count_dt[DI_site == TRUE], aes(x = Hyl_found_in, y = N, fill = met_within_2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g1 + g2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/overlap_wt_ko-1.png)<!-- -->

``` r
merge(
  koh_wt_vs_ko_dt[!duplicated(paste(protein_accession, gene_name)), .(protein_accession, gene_name)],
  koh_per_site_dt[DI_site == TRUE & Hyl_found_in == "both"]
)
```

    ## Key: <protein_accession>
    ##     protein_accession gene_name aa_pos DI_site met_within_2 HeLaWT HeLaJMJD6KO
    ##                <char>    <char>  <int>  <lgcl>       <char>  <int>       <int>
    ##  1:            O60814    H2BC12     21    TRUE           No      1           1
    ##  2:            P13639      EEF2    159    TRUE          Yes      1           1
    ##  3:            P16403      H1-2    137    TRUE           No      1           1
    ##  4:            P16403      H1-2    148    TRUE           No      1           1
    ##  5:            P16403      H1-2    149    TRUE           No      1           1
    ##  6:            P16403      H1-2    152    TRUE           No      1           1
    ##  7:            P16403      H1-2    153    TRUE           No      1           1
    ##  8:            P19338       NCL    282    TRUE          Yes      1           1
    ##  9:            P20908    COL5A1    535    TRUE           No      1           1
    ## 10:            P29375     KDM5A   1495    TRUE           No      1           1
    ## 11:            P29375     KDM5A   1497    TRUE           No      1           1
    ## 12:            P46777      RPL5    264    TRUE           No      1           1
    ## 13:            P62753      RPS6    221    TRUE           No      1           1
    ## 14:            P62805     H4C16     13    TRUE           No      1           1
    ## 15:            P62805     H4C16     92    TRUE           No      1           1
    ## 16:            Q13428     TCOF1   1444    TRUE           No      1           1
    ## 17:            Q8TA86       RP9    195    TRUE           No      1           1
    ## 18:            Q9Y3S2    ZNF330     10    TRUE           No      1           1
    ## 19:            Q9Y3S2    ZNF330     11    TRUE           No      1           1
    ##     Hyl_found_in
    ##           <fctr>
    ##  1:         both
    ##  2:         both
    ##  3:         both
    ##  4:         both
    ##  5:         both
    ##  6:         both
    ##  7:         both
    ##  8:         both
    ##  9:         both
    ## 10:         both
    ## 11:         both
    ## 12:         both
    ## 13:         both
    ## 14:         both
    ## 15:         both
    ## 16:         both
    ## 17:         both
    ## 18:         both
    ## 19:         both

``` r
# Statistical significance of methionine enrichment (all sites vs hydroxylated)
all_k_met_count_dt <- rbind(
  copy(protein_feature_dt[met_within_2 != "edge", .N, by = met_within_2])[, data_type := "all_K"],
  copy(koh_per_site_count_dt[, list(N = sum(N)), by = met_within_2])[, data_type := "koh"]
) %>%
  dcast(met_within_2 ~ data_type, value.var = "N")

fisher.test(all_k_met_count_dt[, .(all_K, koh)], alternative = "two.sided")
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  all_k_met_count_dt[, .(all_K, koh)]
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  4.577065 5.630746
    ## sample estimates:
    ## odds ratio 
    ##   5.079513

``` r
# WT vs KO site overlap (all sites)
plot_overlap_venn(
  wt_ids = koh_wt_vs_ko_dt[genotype == "HeLaWT", paste(protein_accession, aa_pos)],
  ko_ids = koh_wt_vs_ko_dt[genotype == "HeLaJMJD6KO", paste(protein_accession, aa_pos)],
  main   = "Hydroxylated Sites Overlap (All)"
)
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/overlap_wt_ko-2.png)<!-- -->

``` r
# WT vs KO site overlap (diagnostic-ion sites only)
plot_overlap_venn(
  wt_ids = koh_wt_vs_ko_dt[genotype == "HeLaWT" & DI_site == TRUE, paste(protein_accession, aa_pos)],
  ko_ids = koh_wt_vs_ko_dt[genotype == "HeLaJMJD6KO" & DI_site == TRUE, paste(protein_accession, aa_pos)],
  main   = "Hydroxylated Sites Overlap (with DI)"
)
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/overlap_wt_ko-3.png)<!-- -->

# Comparison with PNAS 2022

Benchmarks the diagnostic-ion-based sites against the manually curated
sites from the PNAS 2022 study.

``` r
# PNAS2022 curated stoichiometry data
pnas2022_stoic_dt <- fread(
  file.path(data.dir, "processed_data_from_PNAS2022/long_K_stoichiometry_data.csv")
)
setnames(pnas2022_stoic_dt, old = c("uniprot_id", "position", "residue"), new = c("protein_accession", "aa_pos", "aa"))
pnas2022_stoic_dt[, accession_position := paste0(protein_accession, "_", aa_pos)]

# Flag MQ_DI sites that PNAS2022 curated as hydroxylated
mq_di_stoic_dt[, curated_oxK_site :=
  paste0(protein_accession, "_", aa_pos) %in%
  pnas2022_stoic_dt[curated_oxK_site == TRUE, paste0(protein_accession, "_", aa_pos)]
]

# One hydroxylated WT row per site (DI-preferred)
mq_di_hydroxyk_dt <- mq_di_stoic_dt[ptm == "[Oxidation (K)]" & genotype == "HeLaWT"][
  order(genotype, DI_site)
][!duplicated(paste(protein_accession, aa_pos))]

# Cross-tabulate DI sites against PNAS2022 curated sites
hyl_precision_dt <- mq_di_hydroxyk_dt[, table(DI_site, curated_oxK_site) %>% addmargins] %>% data.table
hyl_precision_dt <- hyl_precision_dt[DI_site %in% c("Sum", "TRUE") & curated_oxK_site != "Sum"]

hyl_precision_dt
```

    ##    DI_site curated_oxK_site     N
    ##     <char>           <char> <num>
    ## 1:    TRUE            FALSE   157
    ## 2:     Sum            FALSE  1782
    ## 3:    TRUE             TRUE    78
    ## 4:     Sum             TRUE   120

``` r
ggplot(
  hyl_precision_dt,
  aes(x = DI_site, y = N, fill = curated_oxK_site)
) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("TRUE" = "coral2", "FALSE" = "#BBBBBB")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(aspect.ratio = 3.5)
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/compare_pnas2022-1.png)<!-- -->

``` r
dcast(hyl_precision_dt, DI_site ~ curated_oxK_site, value.var = "N") %>%
  {fisher.test(.[, c("FALSE", "TRUE"), with = FALSE], alternative = "two.sided")}
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  .[, c("FALSE", "TRUE"), with = FALSE]
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##   5.224579 10.359296
    ## sample estimates:
    ## odds ratio 
    ##   7.365665

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
    ##  eulerr            * 7.0.4      2025-09-24 [1] CRAN (R 4.4.3)
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
    ##  labeling            0.4.3      2023-08-29 [1] CRAN (R 4.4.3)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.4.3)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.4.3)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.4.3)
    ##  patchwork         * 1.3.2      2025-08-25 [1] CRAN (R 4.4.3)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.4.3)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.3)
    ##  polyclip            1.10-7     2024-07-23 [1] CRAN (R 4.4.3)
    ##  polylabelr          0.3.0      2024-11-19 [1] CRAN (R 4.4.3)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-12-13 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.4.3)
    ##  RColorBrewer      * 1.1-3      2022-04-03 [1] CRAN (R 4.4.3)
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
    ##  zlibbioc            1.52.0     2024-10-29 [1] Bioconduc~
    ## 
    ##  [1] /home/ysugimo/R/x86_64-pc-linux-gnu-library/4.4
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
