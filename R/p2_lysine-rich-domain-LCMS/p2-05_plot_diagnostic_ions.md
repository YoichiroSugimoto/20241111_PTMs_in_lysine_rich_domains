p2-05 · Plot diagnostic-ion analysis
================
Yoichiro Sugimoto and Pallavi Kesavan
26 August, 2026

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

# Methionine within +/-2 residues, computed from the reference proteome rather
# than from the pre-computed 11-mer `Window` column. The window is clamped to
# the protein, so lysines near a terminus are classified from the residues that
# exist instead of being discarded: residues that do not exist cannot be
# methionine, so there is no indeterminate ("edge") case.
protein_feature_dt[, met_within_2 := met_within_n(
  protein_accession, aa_pos, all.protein.bs, n = 2L
)]

# `aa_pos` must index the same sequence as the reference proteome. If the FASTA
# used for the MaxQuant search differed from the one loaded by _setup.R, every
# sequence-derived annotation below would be silently wrong.
bad_pos_dt <- check_residue_at_position(
  protein_feature_dt[, protein_accession], protein_feature_dt[, aa_pos],
  all.protein.bs, expected = "K"
)
if (nrow(bad_pos_dt) > 0) {
  warning(nrow(bad_pos_dt), " lysine positions do not match the reference proteome")
  print(head(bad_pos_dt))
}
```

    ## Warning: 3063 lysine positions do not match the reference proteome

    ##    protein_accession aa_pos observed in_reference
    ##               <char>  <int>   <char>       <lgcl>
    ## 1:        A0A0U1RRI6     31        M         TRUE
    ## 2:        A0A0U1RRI6     67        A         TRUE
    ## 3:        A0A0U1RRI6     72        S         TRUE
    ## 4:        A0A0U1RRI6     76        L         TRUE
    ## 5:        A0A0U1RRI6    107        R         TRUE
    ## 6:        A0A0U1RRI6    123        D         TRUE

``` r
# met_within_n() returns NA where the accession is absent from the reference
# proteome. Dropped explicitly: an NA group would turn the 2 x 2 contingency
# table of the Fisher test below into 3 x 2.
n_missing_acc <- protein_feature_dt[is.na(met_within_2), .N]
if (n_missing_acc > 0) {
  warning(n_missing_acc, " lysines dropped: accession absent from the reference proteome")
  protein_feature_dt <- protein_feature_dt[!is.na(met_within_2)]
}
```

    ## Warning: 2283 lysines dropped: accession absent from the reference proteome

``` r
# Lysines at aa_pos <= 3 can be flagged "Yes" on the strength of the initiator
# methionine, which is frequently cleaved co-translationally. No special
# handling is applied: across the canonical human proteome these are ~0.3% of
# lysines, and including all terminal lysines (~1.8% of lysines) raises the
# all_K Yes:No odds by ~3%, which shifts the enrichment odds ratio below by <3%
# -- well within its confidence interval.
protein_feature_dt[aa_pos <= 3 & met_within_2 == "Yes", .N]
```

    ## [1] 1792

``` r
# Category sizes, for the record (there is no longer an "edge" category)
protein_feature_dt[, .N, by = met_within_2]
```

    ##    met_within_2      N
    ##          <char>  <int>
    ## 1:          Yes  56869
    ## 2:           No 588743

``` r
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
  koh_wt_vs_ko_dt,
  aes(x = genotype, fill = met_within_2)
) +
  geom_bar() +
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g2 <- ggplot(
  koh_wt_vs_ko_dt[DI_site == TRUE],
  aes(x = genotype, fill = met_within_2)
) +
  geom_bar() +
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

# Number of sites for plotting g2
koh_wt_vs_ko_dt[DI_site == TRUE][, .N, by = list(met_within_2, genotype)]
```

    ##    met_within_2    genotype     N
    ##          <char>      <fctr> <int>
    ## 1:           No HeLaJMJD6KO    20
    ## 2:           No      HeLaWT   145
    ## 3:          Yes HeLaJMJD6KO     3
    ## 4:          Yes      HeLaWT     6

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
    ## 18:            Q9Y3S2    ZNF330      3    TRUE          Yes      1           1
    ## 19:            Q9Y3S2    ZNF330     10    TRUE           No      1           1
    ## 20:            Q9Y3S2    ZNF330     11    TRUE           No      1           1
    ##     protein_accession gene_name aa_pos DI_site met_within_2 HeLaWT HeLaJMJD6KO
    ##                <char>    <char>  <int>  <lgcl>       <char>  <int>       <int>
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
    ## 20:         both
    ##     Hyl_found_in
    ##           <fctr>

``` r
# Counts underlying the enrichment test, plus the diagnostic-ion subset
met_count_long_dt <- rbind(
  copy(protein_feature_dt[, .N, by = met_within_2])[, data_type := "all_K"],
  copy(koh_wt_vs_ko_dt[genotype == "HeLaWT", .N, by = met_within_2])[, data_type := "koh"],
  copy(koh_wt_vs_ko_dt[genotype == "HeLaWT" & DI_site == TRUE, .N, by = met_within_2])[
    , data_type := "koh_DI"]
)

# Exact n per condition, with the percentage methionine-proximal
met_count_wide_dt <- dcast(met_count_long_dt, data_type ~ met_within_2, value.var = "N")
met_count_wide_dt <- met_count_wide_dt[, .(
  data_type,
  met_absent       = No,
  met_within_2     = Yes,
  total            = No + Yes,
  pct_met_within_2 = round(100 * Yes / (No + Yes), 2)
)]

met_count_wide_dt
```

    ## Key: <data_type>
    ##    data_type met_absent met_within_2  total pct_met_within_2
    ##       <char>      <int>        <int>  <int>            <num>
    ## 1:     all_K     588743        56869 645612             8.81
    ## 2:       koh        898          422   1320            31.97
    ## 3:    koh_DI        145            6    151             3.97

``` r
# Statistical significance of methionine enrichment (all sites vs hydroxylated)
all_k_met_count_dt <- dcast(
  met_count_long_dt[data_type %in% c("all_K", "koh")],
  met_within_2 ~ data_type, value.var = "N"
)

f1 <- fisher.test(all_k_met_count_dt[, .(all_K, koh)], alternative = "two.sided")
f1
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  all_k_met_count_dt[, .(all_K, koh)]
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  4.321882 5.468991
    ## sample estimates:
    ## odds ratio 
    ##    4.86505

``` r
# Statistical significance of DI sites (all sites vs hydroxylated with DI)
all_k_met_di_count_dt <- dcast(
  met_count_long_dt[data_type %in% c("all_K", "koh_DI")],
  met_within_2 ~ data_type, value.var = "N"
)

f2 <- fisher.test(all_k_met_di_count_dt[, .(all_K, koh_DI)], alternative = "two.sided")
f2
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  all_k_met_di_count_dt[, .(all_K, koh_DI)]
    ## p-value = 0.03103
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.1546476 0.9553334
    ## sample estimates:
    ## odds ratio 
    ##  0.4283845

``` r
setNames(
  p.adjust(c(f1$p.value, f2$p.value), method = "holm"),
  c("all_K vs koh", "all_K vs koh_DI")
)
```

    ##    all_K vs koh all_K vs koh_DI 
    ##   1.533732e-123    3.102705e-02

``` r
# Composition compared by the tests above. "Yes" is stacked from the axis so the
# bars are directly comparable; exact counts are in met_count_wide_dt above.
met_plot_dt <- met_count_long_dt[, .(
  data_type = factor(
    data_type,
    levels = c("all_K", "koh", "koh_DI"),
    labels = c("All lysines", "Hydroxylysines", "Hydroxylysines (DI)")
  ),
  met_within_2 = factor(met_within_2, levels = c("No", "Yes")),
  N
)]

ggplot(met_plot_dt, aes(x = data_type, y = N, fill = met_within_2)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_manual(
    values = c("Yes" = "#4477AA", "No" = "#BBBBBB"),
    name   = "Met within ±2 aa",
    breaks = c("Yes", "No")
  ) +
  labs(x = NULL, y = "Proportion of lysines") +
  theme(aspect.ratio = 2) +
  scale_x_discrete(guide = guide_axis(angle = 90))
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/overlap_wt_ko-2.png)<!-- -->

``` r
# WT vs KO site overlap (all sites)
plot_overlap_venn(
  wt_ids = koh_wt_vs_ko_dt[genotype == "HeLaWT", paste(protein_accession, aa_pos)],
  ko_ids = koh_wt_vs_ko_dt[genotype == "HeLaJMJD6KO", paste(protein_accession, aa_pos)],
  main   = "Hydroxylated Sites Overlap (All)"
)
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/overlap_wt_KO_venn-1.png)<!-- -->

``` r
# WT vs KO site overlap (diagnostic-ion sites only)
plot_overlap_venn(
  wt_ids = koh_wt_vs_ko_dt[genotype == "HeLaWT" & DI_site == TRUE, paste(protein_accession, aa_pos)],
  ko_ids = koh_wt_vs_ko_dt[genotype == "HeLaJMJD6KO" & DI_site == TRUE, paste(protein_accession, aa_pos)],
  main   = "Hydroxylated Sites Overlap (with DI)"
)
```

![](p2-05_plot_diagnostic_ions_files/figure-gfm/overlap_wt_KO_venn-2.png)<!-- -->

``` r
# Genotype totals quoted in the main text. These must equal the Euler diagram
# totals above; if they do not, a filter is being applied somewhere unintended.
# Reported for all sites regardless of methionine proximity -- met_within_2 is
# used for the stacked bars and the Fisher test above, not to subset counts.
koh_wt_vs_ko_dt[DI_site == TRUE, .N, by = genotype]   # diagnostic-ion Euler totals
```

    ##       genotype     N
    ##         <fctr> <int>
    ## 1: HeLaJMJD6KO    23
    ## 2:      HeLaWT   151

``` r
koh_wt_vs_ko_dt[, .N, by = genotype]                  # all-sites Euler totals
```

    ##       genotype     N
    ##         <fctr> <int>
    ## 1: HeLaJMJD6KO  1107
    ## 2:      HeLaWT  1320

# Comparison with PNAS 2022

Benchmarks the diagnostic-ion-based sites against the manually curated
sites from the PNAS 2022 study.

Counts here are the wild-type slice of the same site set used above —
sites with at least one PSM in both genotypes — so they match the
wild-type circle of the Euler diagrams.

``` r
# PNAS2022 curated stoichiometry data
pnas2022_stoic_dt <- fread(
  file.path(data.dir, "processed_data_from_PNAS2022/long_K_stoichiometry_data.csv")
)
setnames(pnas2022_stoic_dt, old = c("uniprot_id", "position", "residue"), new = c("protein_accession", "aa_pos", "aa"))
pnas2022_stoic_dt[, accession_position := paste0(protein_accession, "_", aa_pos)]

# One hydroxylated WT row per site (DI-preferred), taken from the same filtered
# site set used throughout this script -- sites with >=1 PSM in both genotypes --
# so every count reported in the manuscript refers to one defined population.
mq_di_hydroxyk_dt <- koh_wt_vs_ko_dt[genotype == "HeLaWT"]

# Flag sites that PNAS2022 curated as hydroxylated
mq_di_hydroxyk_dt[, curated_oxK_site :=
  paste0(protein_accession, "_", aa_pos) %in%
  pnas2022_stoic_dt[curated_oxK_site == TRUE, paste0(protein_accession, "_", aa_pos)]
]

# Numbers quoted in the main text: wild-type sites from the standard filtered
# set, so these must equal the WT totals printed in the previous chunk.
mq_di_hydroxyk_dt[, .N]                    # lysines assigned as hydroxylated
```

    ## [1] 1320

``` r
mq_di_hydroxyk_dt[DI_site == TRUE, .N]     # ... supported by a diagnostic ion
```

    ## [1] 151

``` r
# Cross-tabulate DI sites against PNAS2022 curated sites
hyl_precision_dt <- mq_di_hydroxyk_dt[, table(DI_site, curated_oxK_site) %>% addmargins] %>% data.table
hyl_precision_dt <- hyl_precision_dt[DI_site %in% c("Sum", "TRUE") & curated_oxK_site != "Sum"]

hyl_precision_dt
```

    ##    DI_site curated_oxK_site     N
    ##     <char>           <char> <num>
    ## 1:    TRUE            FALSE    87
    ## 2:     Sum            FALSE  1218
    ## 3:    TRUE             TRUE    64
    ## 4:     Sum             TRUE   102

``` r
ggplot(
  hyl_precision_dt,
  aes(x = DI_site, y = N, fill = curated_oxK_site)
) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, 0.5)) +
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
    ##   5.875848 13.055330
    ## sample estimates:
    ## odds ratio 
    ##   8.761173

# Session information

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.1 (2025-06-13)
    ##  os       Ubuntu 24.04.2 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       Europe/Berlin
    ##  date     2026-08-26
    ##  pandoc   3.4 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.6.42 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package           * version    date (UTC) lib source
    ##  P BiocGenerics      * 0.56.0     2025-10-29 [?] Bioconduc~
    ##  P BiocManager         1.30.27    2025-11-14 [?] CRAN (R 4.5.1)
    ##  P Biostrings        * 2.78.0     2025-10-29 [?] Bioconduc~
    ##  P cellranger          1.1.0      2016-07-27 [?] CRAN (R 4.5.1)
    ##  P cli                 3.6.6      2026-04-09 [?] CRAN (R 4.5.1)
    ##  P crayon              1.5.3      2024-06-20 [?] CRAN (R 4.5.1)
    ##  P data.table        * 1.18.4     2026-05-06 [?] CRAN (R 4.5.1)
    ##  P digest              0.6.39     2025-11-19 [?] CRAN (R 4.5.1)
    ##  P dplyr             * 1.2.1      2026-04-03 [?] CRAN (R 4.5.1)
    ##  P eulerr            * 7.0.2      2024-03-28 [?] CRAN (R 4.5.1)
    ##  P evaluate            1.0.5      2025-08-27 [?] CRAN (R 4.5.1)
    ##  P farver              2.1.2      2024-05-13 [?] CRAN (R 4.5.1)
    ##  P fastmap             1.2.0      2024-05-15 [?] CRAN (R 4.5.1)
    ##  P generics          * 0.1.4      2025-05-09 [?] CRAN (R 4.5.1)
    ##  P ggplot2           * 4.0.3      2026-04-22 [?] CRAN (R 4.5.1)
    ##  P glue                1.8.1      2026-04-17 [?] CRAN (R 4.5.1)
    ##  P gtable              0.3.6      2024-10-25 [?] CRAN (R 4.5.1)
    ##  P htmltools           0.5.9      2025-12-04 [?] CRAN (R 4.5.1)
    ##  P IRanges           * 2.44.0     2025-10-29 [?] Bioconduc~
    ##  P janitor             2.2.1      2024-12-22 [?] CRAN (R 4.5.1)
    ##  P khroma            * 1.17.0     2025-09-29 [?] CRAN (R 4.5.1)
    ##  P knitr             * 1.51       2025-12-20 [?] CRAN (R 4.5.1)
    ##  P labeling            0.4.3      2023-08-29 [?] CRAN (R 4.5.1)
    ##  P lifecycle           1.0.5      2026-01-08 [?] CRAN (R 4.5.1)
    ##  P lubridate           1.9.5      2026-02-04 [?] CRAN (R 4.5.1)
    ##  P magrittr          * 2.0.5      2026-04-04 [?] CRAN (R 4.5.1)
    ##  P otel                0.2.0      2025-08-29 [?] CRAN (R 4.5.1)
    ##  P patchwork         * 1.3.2      2025-08-25 [?] CRAN (R 4.5.1)
    ##  P pillar              1.11.1     2025-09-17 [?] CRAN (R 4.5.1)
    ##  P pkgconfig           2.0.3      2019-09-22 [?] CRAN (R 4.5.1)
    ##  P polyclip            1.10-7     2024-07-23 [?] CRAN (R 4.5.1)
    ##  P polylabelr          1.0.0      2026-01-19 [?] CRAN (R 4.5.1)
    ##    ptm.stoichiometry * 0.0.0.9000 2026-08-25 [1] local (/fast/AG_Sugimoto/home/users/yoichiro/projects/ptm.stoichiometry)
    ##  P R6                  2.6.1      2025-02-15 [?] CRAN (R 4.5.1)
    ##  P RColorBrewer      * 1.1-3      2022-04-03 [?] CRAN (R 4.5.1)
    ##  P Rcpp                1.1.1-1.1  2026-04-24 [?] CRAN (R 4.5.1)
    ##  P readxl            * 1.5.0      2026-05-16 [?] CRAN (R 4.5.1)
    ##    renv                1.1.5      2025-07-24 [1] CRAN (R 4.5.1)
    ##  P rlang               1.2.0      2026-04-06 [?] CRAN (R 4.5.1)
    ##    rmarkdown           2.31       2026-03-26 [1] CRAN (R 4.5.1)
    ##  P S4Vectors         * 0.48.1     2026-04-05 [?] Bioconduc~
    ##  P S7                  0.2.2      2026-04-22 [?] CRAN (R 4.5.1)
    ##  P scales              1.4.0      2025-04-24 [?] CRAN (R 4.5.1)
    ##  P Seqinfo           * 1.0.0      2025-10-29 [?] Bioconduc~
    ##  P sessioninfo         1.2.4      2026-06-04 [?] CRAN (R 4.5.1)
    ##  P snakecase           0.11.1     2023-08-27 [?] CRAN (R 4.5.1)
    ##  P stringi             1.8.7      2025-03-27 [?] CRAN (R 4.5.1)
    ##  P stringr           * 1.6.0      2025-11-04 [?] CRAN (R 4.5.1)
    ##  P tibble              3.3.1      2026-01-11 [?] CRAN (R 4.5.1)
    ##  P tidyselect          1.2.1      2024-03-11 [?] CRAN (R 4.5.1)
    ##  P timechange          0.4.0      2026-01-29 [?] CRAN (R 4.5.1)
    ##  P vctrs               0.7.3      2026-04-11 [?] CRAN (R 4.5.1)
    ##  P withr               3.0.3      2026-06-19 [?] CRAN (R 4.5.1)
    ##  P xfun                0.59       2026-06-19 [?] CRAN (R 4.5.1)
    ##  P XVector           * 0.50.0     2025-10-29 [?] Bioconduc~
    ##  P yaml                2.3.12     2025-12-10 [?] CRAN (R 4.5.1)
    ## 
    ##  [1] /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/renv/library/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu
    ##  [2] /home/ysugimo/.cache/R/renv/sandbox/linux-ubuntu-noble/R-4.5/x86_64-pc-linux-gnu/9a444a72
    ## 
    ##  * ── Packages attached to the search path.
    ##  P ── Loaded and on-disk path mismatch.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
