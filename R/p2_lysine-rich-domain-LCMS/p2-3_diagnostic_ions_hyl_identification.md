2-3. Identify diagnostic ions that can mark lysine hydroxylations
================
Yoichiro Sugimoto and Pallavi Kesavan
03 April, 2026

- [Overview](#overview)
- [Environment setup](#environment-setup)
- [2.3.1 Import basic data](#231-import-basic-data)
- [2.3.2 Identify useful diagnostic ions to identify lysine
  hydroxylations](#232-identify-useful-diagnostic-ions-to-identify-lysine-hydroxylations)
- [Session information](#session-information)

# Overview

This script identify diagnostic ions that can aid to determine
hydroxylation sites.

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

#renv::restore(file.path(project.dir, "R"))
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

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, saveRDS, setdiff,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

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
# Install ptm.stiochiometry package - package installed 
#install.packages("/fast/AG_Sugimoto/home/users/pallavi/projects/ptm.stoichiometry", repos = NULL, type = "source")

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
# library("ggpubr")
```

# 2.3.1 Import basic data

``` r
# Define path to data directory
data.dir <- file.path(project.dir, "data")

# Load diagnostic ion data
diagnostic_ion_data <- file.path(
  data.dir, "FP_diagnostic_ion_search/fragpipe_dataset-A/dataset01.diagnosticIons.tsv" 
) %>% 
  fread %>%
  clean_names
```

    ## Warning in fread(.): Discarded single-line footer: <<COMPLETE>>

``` r
# Load Fragpipe psm data 
fragpipe_psm <- file.path(
  data.dir, "FP_diagnostic_ion_search/fragpipe_dataset-A/psm.tsv"
) %>%
  fread %>%
  clean_names

# Load MQ standard data
MQ_Std_data <- file.path(
  data.dir,
  "MQ_standard/PNAS2022" 
)

# Load Data A sample info from data directory 
dataA_sample_info <- fread(file.path(
  MQ_Std_data, "sample_info/MS_dataset_overview_PXD031221_data-A.csv"
))

# Import human protein reference data
ref_protein_dt <- import_reference_fasta(
  file.path(
    "/fast/AG_Sugimoto/reference/uniprot/human",
    "UP000005640_9606.fasta"
  )
)
```

# 2.3.2 Identify useful diagnostic ions to identify lysine hydroxylations

``` r
# Merge 'diagnostic_ion_data' and 'fragpipe_psm' data by column 'spectrum'
diagnostic_ion_data <- merge(
  fragpipe_psm[, .(spectrum, protein_id, gene, protein_start, protein_end)], 
  diagnostic_ion_data, 
  by = "spectrum"
)

# Create new column 'file_name' from column 'spectrum'  
diagnostic_ion_data[, `:=`(
  file_name = str_split_fixed(spectrum, "\\.", n = 2)[, 1] # split 'spectrum' by the first '.', take the part before it 
)] 

# Merge 'dataA_sample_info' and 'diagnostic_ion_data' by column 'file_name'
diagnostic_ion_data <- merge(
  dataA_sample_info,
  diagnostic_ion_data,
  by = "file_name"
)

# Retrieve hydroxylation info for BRD4
brd4_hydroxylation_site_di <- diagnostic_ion_data[
  gene == "BRD4" & (555 >= protein_start & 535 <= protein_end)
]

# Reshape data table from wide to long format 
m.brd4_hydroxylation_site_di <- melt(
  brd4_hydroxylation_site_di,
  measure.vars = grep("intensity$",colnames(brd4_hydroxylation_site_di), value = TRUE), # get column names whose names end with "intensity" 
  value.name = "intensity" # new column - "intensity" 
)

# Add columns 'monoisotopic_mass', 'diagnostic_ion' and 'genotype'
m.brd4_hydroxylation_site_di[, `:=`(
  monoisotopic_mass = variable %>% # take values from column 'variable'
    str_replace_all("ox_", "") %>% # remove the prefix "ox_"
    str_replace_all("_intensity", "") %>% # remove the suffix "_intensity"
    str_replace_all("_", ".") %>% # replace remaining underscores with dots
    factor(levels = c(
      "101.1079", #Lysine   K           immonium ion    C5 H13 N2 +
      "100.0762", #Hydroxylation (K)    K   O   15.9949 diagnostic ion  C5 H10 N O+
      "82.0657", #Hydroxylation (K) K   O   15.9949 diagnostic ion (water loss) C5 H8 N+
      "117.1028", #Hydroxylation (K)    K   O   15.9949 immonium ion    C5 H13 O N2 + calculated
      "117.0658", #Hydroxylation (K)    K   O   15.9949 immonium ion    C5 H13 O N2 + (from citation)
      "138.0919", #Hydroxylation-Propionylation (K) K   C3 H4 O2    72.02112937 diagnostic ion (water loss) C8 H12 N1 O+
      "145.0977", #Hydroxylation (K)    K   O   15.9949 (intact, water loss)    C6H13N2O2+
      "156.1025", #Hydroxylation-Propionylation (K) K   C3 H4 O2    72.02112937 diagnostic ion  C8 H14 N1 O2+
      "173.1290" #Hydroxylation-Propionylation (K)  K   C3 H4 O2    72.02112937 immonium ion    C8 H17 O2 N2 +
    )),
  diagnostic_ion = intensity > 0, 
  genotype = factor(genotype, levels = c("WT", "JMJD6KO"))
)]


# Load libraries
library("ggbeeswarm")
library("khroma")

# Plot - Quasirandom plot of monoisotopic mass versus intensity in WT and JMJD6KO genotype
ggplot(
  data = m.brd4_hydroxylation_site_di,
  aes(
    x = monoisotopic_mass,
    y = intensity,
    color = genotype
  )
) +
  geom_quasirandom(dodge.width=0.5, size = 0.2) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_color_bright()
```

![](p2-3-diagnostic-ions-lysine-hydorxylation-stoic_files/figure-gfm/identify_diagnostic_ions-1.png)<!-- -->

``` r
# Plot - Violin plot of monoisotopic mass versus intensity in WT and JMJD6KO genotype
ggplot(
  data = m.brd4_hydroxylation_site_di,
  aes(
    x = monoisotopic_mass,
    y = intensity,
    fill = genotype
  )
) +
  geom_violin(scale='width') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_bright()
```

![](p2-3-diagnostic-ions-lysine-hydorxylation-stoic_files/figure-gfm/identify_diagnostic_ions-2.png)<!-- -->

``` r
selected.ms <- c("100.0762", "156.1025")

# Count the number of PSMs for each selected monoisotopic mass 
psm.count.by.di <- m.brd4_hydroxylation_site_di[monoisotopic_mass %in% c(selected.ms)][, .N, by = list(
  monoisotopic_mass, genotype, diagnostic_ion
)]

print("Ratio of PSM with diagnostic ion: 100.0762")
```

    ## [1] "Ratio of PSM with diagnostic ion: 100.0762"

``` r
psm.count.by.di[monoisotopic_mass == "100.0762" & genotype == "WT" & diagnostic_ion == TRUE, N] /
  psm.count.by.di[monoisotopic_mass == "100.0762" & genotype == "WT", sum(N)]
```

    ## [1] 0.04388298

``` r
print("Ratio of PSM with diagnostic ion: 156.1025")
```

    ## [1] "Ratio of PSM with diagnostic ion: 156.1025"

``` r
psm.count.by.di[monoisotopic_mass == "156.1025" & genotype == "WT" & diagnostic_ion == TRUE, N] /
  psm.count.by.di[monoisotopic_mass == "156.1025" & genotype == "WT", sum(N)]
```

    ## [1] 0.2087766

``` r
# Perform Fisher's test for each selected monoisotopic mass 
for(mi.ms in selected.ms){
  print(mi.ms)
  
  psm.count.by.di[monoisotopic_mass == mi.ms] %>% # subsets psm counts for current mi.ms 
    dcast(genotype ~ diagnostic_ion, value.var = "N") %>% # convert long to wide format data table
    setnafill(cols = c("FALSE", "TRUE"), fill = 0)  %>% # Fill NA values as 0 
    {as.matrix(.[, c("FALSE", "TRUE"), with = FALSE])} %>% # convert into data matrix 
    fisher.test %>% print # perform Fisher's test and print results 
}
```

    ## [1] "100.0762"
    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  .
    ## p-value = 4.92e-09
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.0000000 0.1460815
    ## sample estimates:
    ## odds ratio 
    ##          0 
    ## 
    ## [1] "156.1025"
    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  .
    ## p-value < 2.2e-16
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.0001688554 0.0369678626
    ## sample estimates:
    ##  odds ratio 
    ## 0.006489756

``` r
# Plot - Violin plot of selected monoisotopic mass versus intensity in WT and JMJD6KO genotype
ggplot(
  data = m.brd4_hydroxylation_site_di[monoisotopic_mass %in% c(selected.ms)],
  aes(
    x = monoisotopic_mass,
    y = intensity,
    fill = genotype
  )
) +
  geom_violin(scale='width') +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_bright()
```

![](p2-3-diagnostic-ions-lysine-hydorxylation-stoic_files/figure-gfm/plot_with_only_diagnostic_ions-1.png)<!-- -->

``` r
m.brd4_hydroxylation_site_di[monoisotopic_mass %in% c(selected.ms)][, .N, by = list(monoisotopic_mass, genotype)]
```

    ##    monoisotopic_mass genotype     N
    ##               <fctr>   <fctr> <int>
    ## 1:          100.0762       WT   752
    ## 2:          100.0762  JMJD6KO   586
    ## 3:          156.1025       WT   752
    ## 4:          156.1025  JMJD6KO   586

``` r
# Plot - Bar chart of proportion of PSMs with diagnostic ion 
ggplot(
  data = m.brd4_hydroxylation_site_di[monoisotopic_mass %in% c(selected.ms)],
  aes(
    x = genotype,
    fill = diagnostic_ion
  )
) +
  geom_bar(position = "fill") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(~ monoisotopic_mass) +
  scale_fill_manual(values = c("TRUE" = "#A50026", "FALSE" = "#DDDDDD")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ylab("Proportion of PSMs with diagnostic ions")
```

![](p2-3-diagnostic-ions-lysine-hydorxylation-stoic_files/figure-gfm/plot_with_only_diagnostic_ions-2.png)<!-- -->

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
    ##  date     2026-04-03
    ##  pandoc   3.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.5.57 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  beeswarm            0.4.0      2021-06-01 [1] CRAN (R 4.4.3)
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
    ##  ggbeeswarm        * 0.7.2      2023-04-29 [1] CRAN (R 4.4.3)
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
    ##  labeling            0.4.3      2023-08-29 [1] CRAN (R 4.4.3)
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
    ##  vipor               0.4.7      2023-12-18 [1] CRAN (R 4.4.3)
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
