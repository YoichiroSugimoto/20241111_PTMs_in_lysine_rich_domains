2-4. Analysis of lysine hydroxylations using diagnostic ions
================
Yoichiro Sugimoto and Pallavi Kesavan
01 June, 2026

- [Overview](#overview)
- [Environment setup](#environment-setup)
- [2.4.1 Import basic data](#241-import-basic-data)
- [2.4.2 Definition of functions](#242-definition-of-functions)
- [2.4.3 The effect of diagnostic ion on
  precision](#243-the-effect-of-diagnostic-ion-on-precision)
- [2.4.4 Overlaps of Hyl sites in WT and JMJD6 KO
  cells](#244-overlaps-of-hyl-sites-in-wt-and-jmjd6-ko-cells)
- [2.4.5 Comparison with previous PNAS 2022
  paper](#245-comparison-with-previous-pnas-2022-paper)
- [Session information](#session-information)

# Overview

This script examines how the use of diagnostic ions improve the analysis
of lysine hydroxylations.

# Environment setup

``` r
# renv::init(
#           "/fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/R"
#       )

# Define project directory - contains R scripts, data and results folders
project.dir <-
  file.path(
    "/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains"
  )

# renv::snapshot(file.path(project.dir, "R"))

# renv::restore(file.path(project.dir, "R"))
```

``` r
# Load all R scripts from the 'functions' folder into the current session
P2_functions <-
  sapply(list.files(
    file.path(project.dir, "R/functions"),
    pattern = "*.R",
    full.names = TRUE
  ),
  source)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     first, second

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:data.table':
    ## 
    ##     shift

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
### Install private packages 
# Install ptm.stiochiometry package - package installed 16.12.2025
# install.packages("/fast/AG_Sugimoto/home/users/pallavi/projects/ptm.stoichiometry", repos = NULL, type = "source")


# Load Libraries - ptm.stiochiometry,readxl and janitor
library("readxl")
library("janitor")
```

    ## 
    ## Attaching package: 'janitor'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
library("ptm.stoichiometry")
```

# 2.4.1 Import basic data

``` r
# Define path to data directory
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

# Import human protein reference data from specified file path 
ref_protein_dt <- import_reference_fasta(
  file.path(
    "/fast/AG_Sugimoto/reference/uniprot/human",
    "UP000005640_9606.fasta"
  )
)

#----------------
# MQ_DI_noH2Oloss
#----------------

# Load sample run info data (MQ DI with no waterloss)
MQ_DI_noH2O_run_info <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_with_DI_noH2O" 
) %>% data.table

# Add a new column 'sample_id' that assigns a unique sequential number to each 
#row (1:N), where N is the total number of rows in the data.table
MQ_DI_noH2O_run_info[, sample_id := 1:.N]
```

# 2.4.2 Definition of functions

``` r
# Define function - 'read_stoic_data'
read_stoic_data <- function(prefix, pre_prefix, post_fix = "", dir_path){
  # Read csv file with defined file path,  
  dt <- fread(file.path(dir_path, paste0(pre_prefix, prefix, "PTM_stoichiometry", post_fix, ".csv"))) 
  dt[, condition := gsub("_$", "", prefix)] # Removes '_' at the end of the prefix and creates a new column called condition
  return(dt)
}
```

# 2.4.3 The effect of diagnostic ion on precision

``` r
# Read stoichiometry data for MQ_DI (change from the older version: a potential diagnostic ions with water loss is no longer considered)
MQ_DI_noH20_stoic_dt <- lapply(
  MQ_DI_noH2O_run_info[, prefix],
  read_stoic_data,
  pre_prefix = "DI_noH2O_", 
  post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MQ_DI_noH2O")
) %>% rbindlist

# fwrite(MQ_DI_noH20_stoic_dt[ptm != ""], file = file.path(results.dir, "p2-analysis-setting", "all_hydroxylysine.csv"))

MQ_DI_noH20_stoic_dt <- MQ_DI_noH20_stoic_dt[aa == "K"]
MQ_DI_noH20_stoic_dt[, genotype := str_split_fixed(sample_name, "_", n = 3)[, 2] %>% factor(levels = c("HeLaWT", "HeLaJMJD6KO"))]

# Extract the sites with diagnostic ions
DI_sites <- MQ_DI_noH20_stoic_dt[diagnostic_peak == "+"][!duplicated(paste(protein_accession, aa_pos)), .(protein_accession, aa_pos)]

MQ_DI_noH20_stoic_dt[, DI_site :=
  paste(protein_accession, aa_pos) %in% DI_sites[, paste(protein_accession, aa_pos)] 
]

# Load protein feature data 
protein.feature.dt <- fread(file.path(data.dir, "PNAS2022/all_protein_feature_per_position.csv"))
protein.feature.dt <- protein.feature.dt[residue == "K"]
setnames(protein.feature.dt, old = c("uniprot_id", "position"), new = c("protein_accession", "aa_pos")) # change column names 

protein.feature.dt[, `:=`( #create a new column 'met_within_2'
  met_within_2 = case_when(
    nchar(Window) == 11 & grepl("M", substr(Window, start = 4, stop = 8)) ~ "Yes", # if M is in the middle, Yes
    nchar(Window) == 11 ~ "No", # if M is not in the sequence, then No 
    TRUE ~ "edge" # if M is in the sequence but not in the middle, then edge
  )
)]

# Merge data tables by protein_accession and aa_pos (MQ_DI)
wt_vs_KO_MQ_DI_noH20_stoic_dt <- copy(MQ_DI_noH20_stoic_dt)

wt_vs_KO_MQ_DI_noH20_stoic_dt <- merge(
  wt_vs_KO_MQ_DI_noH20_stoic_dt,
  protein.feature.dt,
  by = c("protein_accession", "aa_pos")
)

# Only analyse the sites with a coverage with 2 PSMs in both JMJD6 WT and KO data
wt_data_coverage <- wt_vs_KO_MQ_DI_noH20_stoic_dt[genotype == "HeLaWT"][
  , list(total_wt_sum_psm_mapped = sum(sum_psm_mapped)), by = list(protein_accession, aa_pos)]
KO_data_coverage <- wt_vs_KO_MQ_DI_noH20_stoic_dt[genotype == "HeLaJMJD6KO"][
  , list(total_ko_sum_psm_mapped = sum(sum_psm_mapped)), by = list(protein_accession, aa_pos)]

all_data_coverage <- merge(
  wt_data_coverage,
  KO_data_coverage,
  by = c("protein_accession", "aa_pos")
)

wt_vs_KO_MQ_DI_noH20_stoic_dt <- wt_vs_KO_MQ_DI_noH20_stoic_dt[
  paste(protein_accession, aa_pos) %in% 
    all_data_coverage[total_wt_sum_psm_mapped > 0 & total_ko_sum_psm_mapped > 0][, paste(protein_accession, aa_pos)]
]

# The # of sites analysed
wt_vs_KO_MQ_DI_noH20_stoic_dt[!duplicated(paste(genotype, protein_accession, aa_pos))][, .N, by = genotype]
```

    ##       genotype     N
    ##         <fctr> <int>
    ## 1:      HeLaWT 26074
    ## 2: HeLaJMJD6KO 26074

``` r
koh_wt_vs_KO_MQ_DI_noH20_stoic_dt <- wt_vs_KO_MQ_DI_noH20_stoic_dt[ptm == "[Oxidation (K)]"]
koh_wt_vs_KO_MQ_DI_noH20_stoic_dt <- koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[order(-DI_site)][!duplicated(paste(genotype, protein_accession, aa_pos))]

g1 <- ggplot(
  data = koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[met_within_2 != "edge"],
  aes(
    x = genotype,
    fill = met_within_2
  )
) +
  geom_bar() +  
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g2 <- ggplot(
  data = koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[met_within_2 != "edge" & DI_site == TRUE],
  aes(
    x = genotype,
    fill = met_within_2
  )
) +
  geom_bar() +  
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

# Number of sites for plotting g2
koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[met_within_2 != "edge" & DI_site == TRUE][
  , .N, by = list(met_within_2, genotype)
]
```

    ##    met_within_2    genotype     N
    ##          <char>      <fctr> <int>
    ## 1:           No HeLaJMJD6KO    20
    ## 2:           No      HeLaWT   143
    ## 3:          Yes HeLaJMJD6KO     2
    ## 4:          Yes      HeLaWT     5

``` r
library("patchwork")

g1 + g2 + plot_layout(guides = "collect") & theme(legend.position = "bottom") 
```

![](p2-4-diagnostic-ions_files/figure-gfm/precision-1.png)<!-- -->

# 2.4.4 Overlaps of Hyl sites in WT and JMJD6 KO cells

``` r
koh_per_site <- dcast(
  koh_wt_vs_KO_MQ_DI_noH20_stoic_dt,
  protein_accession + aa_pos + DI_site + met_within_2 ~ genotype,
  fun.aggregate = length 
)

koh_per_site[, `:=`(Hyl_found_in = case_when(
  HeLaWT == 1 & HeLaJMJD6KO == 1 ~ "both",
  HeLaWT == 1 ~ "WT",
  HeLaJMJD6KO == 1 ~ "JMJD6KO",
  TRUE ~ "not_found"
) %>% factor(levels = c("WT", "JMJD6KO", "both", "not_found")))]

koh_per_site <- koh_per_site[met_within_2 != "edge"]

koh_per_site_count <- koh_per_site[, .N, by = list(Hyl_found_in, met_within_2, DI_site)]

g1 <- ggplot(koh_per_site_count, aes(x = Hyl_found_in, y = N, fill = met_within_2)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g2 <- ggplot(koh_per_site_count[DI_site == TRUE], aes(x = Hyl_found_in, y = N, fill = met_within_2)) +
  geom_bar(stat = "identity") +  
  scale_fill_manual(values = c("Yes" = "#4477AA", "No" = "#BBBBBB")) +
  theme(aspect.ratio = 3) +
  scale_x_discrete(guide = guide_axis(angle = 90))

g1 + g2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
```

![](p2-4-diagnostic-ions_files/figure-gfm/overlap_wt_ko-1.png)<!-- -->

``` r
merge(
  koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[!duplicated(paste(protein_accession, gene_name)), .(protein_accession, gene_name)],
  koh_per_site[DI_site == TRUE & Hyl_found_in == "both"]
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
# Statistical significance of methionine enrichment
# All sites, found in both
all.K.m.count <- rbind(
  copy(protein.feature.dt[met_within_2 != "edge", .N, by = met_within_2])[, data_type := "all_K"],
  copy(koh_per_site_count[, list(N = sum(N)), by = met_within_2])[, data_type := "koh"]
) %>%
  dcast(
    met_within_2 ~ data_type,
    value.var = "N"
  )

fisher.test(all.K.m.count[, .(all_K, koh)], alternative = "two.sided")
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  all.K.m.count[, .(all_K, koh)]
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  4.577065 5.630746
    ## sample estimates:
    ## odds ratio 
    ##   5.079513

``` r
library("eulerr")
library("RColorBrewer")

vennlist <-  list(
  WT_DI_site = koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[genotype == "HeLaWT", paste(protein_accession, aa_pos)],
  KO_DI_site = koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[genotype == "HeLaJMJD6KO", paste(protein_accession, aa_pos)]
)

# Changing the circle size based on the number in the circle
venn_size_based <- euler(vennlist)

cols <- brewer.pal(3, "Set2")

# Plot the diagram
plot(venn_size_based,
     fills = list(fill = cols, alpha = 0.4),
     legend = list(side = "right"),
     quantities = TRUE,
     main = "Hydroxylated Sites Overlap (All)")
```

![](p2-4-diagnostic-ions_files/figure-gfm/overlap_wt_ko-2.png)<!-- -->

``` r
vennlist <-  list(
  WT_DI_site = koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[genotype == "HeLaWT" & DI_site == TRUE, paste(protein_accession, aa_pos)],
  KO_DI_site = koh_wt_vs_KO_MQ_DI_noH20_stoic_dt[genotype == "HeLaJMJD6KO" & DI_site == TRUE, paste(protein_accession, aa_pos)]
)

# Changing the circle size based on the number in the circle
venn_size_based <- euler(vennlist)

cols <- brewer.pal(3, "Set2")

# Plot the diagram
plot(venn_size_based,
     fills = list(fill = cols, alpha = 0.4),
     legend = list(side = "right"),
     quantities = TRUE,
     main = "Hydroxylated Sites Overlap (with DI)")
```

![](p2-4-diagnostic-ions_files/figure-gfm/overlap_wt_ko-3.png)<!-- -->

# 2.4.5 Comparison with previous PNAS 2022 paper

``` r
# Define file path to load 'long_K_stiochiometry_data'  
pnas2022.stoic.dt <- fread(file.path(data.dir, "PNAS2022/long_K_stoichiometry_data.csv"))

#change the column names
setnames(pnas2022.stoic.dt, old = c("uniprot_id", "position", "residue"), new = c("protein_accession", "aa_pos", "aa")) 

# Create new column 'accession_position'
pnas2022.stoic.dt[, `:=`(
  accession_position = paste0(protein_accession, "_", aa_pos) # paste protein_accession and aa_pos together '_'
)]

# Check how many hydroxylation sites that were reported by PNAS2022 are identified by the new workflow
MQ_DI_noH20_stoic_dt[, `:=`( 
  curated_oxK_site = 
    paste0(protein_accession, "_", aa_pos) %in%
    pnas2022.stoic.dt[curated_oxK_site == TRUE, paste0(protein_accession, "_", aa_pos)] 
)]

MQ_DI_d.hydroxyK_dt <- MQ_DI_noH20_stoic_dt[ptm == "[Oxidation (K)]" & genotype == "HeLaWT"][order(genotype, DI_site)][!duplicated(paste(protein_accession, aa_pos))]

# Convert the data table into table and sum the number of diagnostic peak and curated oxK sites
hyl_precision <- MQ_DI_d.hydroxyK_dt[, table(DI_site, curated_oxK_site) %>% addmargins] %>% data.table
hyl_precision <- hyl_precision[DI_site %in% c("Sum", "TRUE") & curated_oxK_site != "Sum"]

hyl_precision
```

    ##    DI_site curated_oxK_site     N
    ##     <char>           <char> <num>
    ## 1:    TRUE            FALSE   157
    ## 2:     Sum            FALSE  1782
    ## 3:    TRUE             TRUE    78
    ## 4:     Sum             TRUE   120

``` r
ggplot(
  hyl_precision,
  aes(
    x = DI_site,
    y = N,
    fill = curated_oxK_site
  )
) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("TRUE" = "coral2", "FALSE" = "#BBBBBB")) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(aspect.ratio = 3.5)
```

![](p2-4-diagnostic-ions_files/figure-gfm/pnas2022_comparison_MQ_DI_noH20loss-1.png)<!-- -->

``` r
dcast(
  hyl_precision,
  DI_site ~ curated_oxK_site,
  value.var = "N"
) %>%
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
    ##  version  R version 4.5.1 (2025-06-13)
    ##  os       Ubuntu 24.04.2 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       Europe/Berlin
    ##  date     2026-06-01
    ##  pandoc   3.4 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.6.42 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  BiocGenerics      * 0.54.0     2025-04-15 [1] Bioconduc~
    ##  Biostrings        * 2.76.0     2025-04-15 [1] Bioconduc~
    ##  cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.5.1)
    ##  cli                 3.6.5      2025-04-23 [1] CRAN (R 4.5.1)
    ##  crayon              1.5.3      2024-06-20 [1] CRAN (R 4.5.1)
    ##  data.table        * 1.17.8     2025-07-10 [1] CRAN (R 4.5.1)
    ##  digest              0.6.37     2024-08-19 [1] CRAN (R 4.5.1)
    ##  dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.5.1)
    ##  eulerr            * 7.0.4      2025-09-24 [1] CRAN (R 4.5.1)
    ##  evaluate            1.0.5      2025-08-27 [1] CRAN (R 4.5.1)
    ##  farver              2.1.2      2024-05-13 [1] CRAN (R 4.5.1)
    ##  fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.5.1)
    ##  generics          * 0.1.4      2025-05-09 [1] CRAN (R 4.5.1)
    ##  GenomeInfoDb      * 1.44.3     2025-09-21 [1] Bioconduc~
    ##  GenomeInfoDbData    1.2.14     2025-09-24 [1] Bioconductor
    ##  ggplot2           * 4.0.1      2025-11-14 [1] CRAN (R 4.5.1)
    ##  glue                1.8.0      2024-09-30 [1] CRAN (R 4.5.1)
    ##  gtable              0.3.6      2024-10-25 [1] CRAN (R 4.5.1)
    ##  htmltools           0.5.8.1    2024-04-04 [1] CRAN (R 4.5.1)
    ##  httr                1.4.7      2023-08-15 [1] CRAN (R 4.5.1)
    ##  IRanges           * 2.42.0     2025-04-15 [1] Bioconduc~
    ##  janitor           * 2.2.1      2024-12-22 [1] CRAN (R 4.5.1)
    ##  jsonlite            2.0.0      2025-03-27 [1] CRAN (R 4.5.1)
    ##  khroma            * 1.16.0     2025-02-25 [1] CRAN (R 4.5.1)
    ##  knitr             * 1.50       2025-03-16 [1] CRAN (R 4.5.1)
    ##  labeling            0.4.3      2023-08-29 [1] CRAN (R 4.5.1)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.5.1)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.5.1)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.5.1)
    ##  patchwork         * 1.3.2      2025-08-25 [1] CRAN (R 4.5.1)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.5.1)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.5.1)
    ##  polyclip            1.10-7     2024-07-23 [1] CRAN (R 4.5.1)
    ##  polylabelr          0.3.0      2024-11-19 [1] CRAN (R 4.5.1)
    ##  ptm.stoichiometry * 0.0.0.9000 2026-05-15 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.5.1)
    ##  RColorBrewer      * 1.1-3      2022-04-03 [1] CRAN (R 4.5.1)
    ##  Rcpp                1.1.0      2025-07-02 [1] CRAN (R 4.5.1)
    ##  readxl            * 1.4.5      2025-03-07 [1] CRAN (R 4.5.1)
    ##  rlang               1.1.6      2025-04-11 [1] CRAN (R 4.5.1)
    ##  rmarkdown           2.29       2024-11-04 [1] CRAN (R 4.5.1)
    ##  rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.5.1)
    ##  S4Vectors         * 0.46.0     2025-04-15 [1] Bioconduc~
    ##  S7                  0.2.0      2024-11-07 [1] CRAN (R 4.5.1)
    ##  scales              1.4.0      2025-04-24 [1] CRAN (R 4.5.1)
    ##  sessioninfo         1.2.3      2025-02-05 [1] CRAN (R 4.5.1)
    ##  snakecase           0.11.1     2023-08-27 [1] CRAN (R 4.5.1)
    ##  stringi             1.8.7      2025-03-27 [1] CRAN (R 4.5.1)
    ##  stringr           * 1.5.2      2025-09-08 [1] CRAN (R 4.5.1)
    ##  tibble              3.3.0      2025-06-08 [1] CRAN (R 4.5.1)
    ##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.5.1)
    ##  timechange          0.3.0      2024-01-18 [1] CRAN (R 4.5.1)
    ##  UCSC.utils          1.4.0      2025-04-15 [1] Bioconduc~
    ##  vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.5.1)
    ##  withr               3.0.2      2024-10-28 [1] CRAN (R 4.5.1)
    ##  xfun                0.53       2025-08-19 [1] CRAN (R 4.5.1)
    ##  XVector           * 0.48.0     2025-04-15 [1] Bioconduc~
    ##  yaml                2.3.10     2024-07-26 [1] CRAN (R 4.5.1)
    ## 
    ##  [1] /home/pkesava/R/x86_64-pc-linux-gnu-library/4.5
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
