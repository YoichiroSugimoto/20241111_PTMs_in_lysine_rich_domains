2-3. Identify diagnostic ions that can mark lysine hydroxylations
================
Yoichiro Sugimoto
25 November, 2025

- [Environment setup](#environment-setup)
- [Import basic data](#import-basic-data)
- [Identify useful diagnostic ions to identify lysine
  hydroxylations](#identify-useful-diagnostic-ions-to-identify-lysine-hydroxylations)
- [Session information](#session-information)

This script identify diagnostic ions to identify confident hydroxylation
sites.

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
library("ggpubr")
```

# Import basic data

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

# Load PNAS2022 MQ standard data
pnas2022_data <- file.path(
  data.dir,
  "MQ_standard/PNAS2022" 
)

# Load Data A sample info from data directory 
dataA_sample_info <- fread(file.path(
  pnas2022_data, "sample_info/MS_dataset_overview_PXD031221_data-A.csv"
))

# Import human protein reference data
ref_protein_dt <- import_reference_fasta(
  file.path(
    "/fast/AG_Sugimoto/reference/uniprot/human",
    "UP000005640_9606.fasta"
  )
)
```

# Identify useful diagnostic ions to identify lysine hydroxylations

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
m.brd4_hydroxylation_site_di[monoisotopic_mass %in% c(selected.ms)]
```

    ##                   file_name                 sample_name   cell genotype
    ##                      <char>                      <char> <char>   <fctr>
    ##    1:          201021_MC264      JQ1_HeLaWT_derivatised   HeLa       WT
    ##    2:          201021_MC264      JQ1_HeLaWT_derivatised   HeLa       WT
    ##    3:          201021_MC264      JQ1_HeLaWT_derivatised   HeLa       WT
    ##    4:          201021_MC264      JQ1_HeLaWT_derivatised   HeLa       WT
    ##    5:          201021_MC264      JQ1_HeLaWT_derivatised   HeLa       WT
    ##   ---                                                                  
    ## 2672: 20201119_GV2048_MC278 JQ1_HeLaJMJD6KO_derivatised   HeLa  JMJD6KO
    ## 2673: 20201119_GV2048_MC278 JQ1_HeLaJMJD6KO_derivatised   HeLa  JMJD6KO
    ## 2674: 20201119_GV2048_MC278 JQ1_HeLaJMJD6KO_derivatised   HeLa  JMJD6KO
    ## 2675: 20201119_GV2048_MC278 JQ1_HeLaJMJD6KO_derivatised   HeLa  JMJD6KO
    ## 2676: 20201119_GV2048_MC278 JQ1_HeLaJMJD6KO_derivatised   HeLa  JMJD6KO
    ##       purification derivitisation replicate                            spectrum
    ##             <char>         <char>     <int>                              <char>
    ##    1:          JQ1            yes         1          201021_MC264.03134.03134.3
    ##    2:          JQ1            yes         1          201021_MC264.03140.03140.2
    ##    3:          JQ1            yes         1          201021_MC264.03206.03206.3
    ##    4:          JQ1            yes         1          201021_MC264.03215.03215.3
    ##    5:          JQ1            yes         1          201021_MC264.03218.03218.2
    ##   ---                                                                          
    ## 2672:          JQ1            yes         5 20201119_GV2048_MC278.62287.62287.3
    ## 2673:          JQ1            yes         5 20201119_GV2048_MC278.62308.62308.3
    ## 2674:          JQ1            yes         5 20201119_GV2048_MC278.62374.62374.3
    ## 2675:          JQ1            yes         5 20201119_GV2048_MC278.62780.62780.4
    ## 2676:          JQ1            yes         5 20201119_GV2048_MC278.62813.62813.4
    ##       protein_id   gene protein_start protein_end
    ##           <char> <char>         <int>       <int>
    ##    1:     O60885   BRD4           553         562
    ##    2:     O60885   BRD4           553         562
    ##    3:     O60885   BRD4           538         546
    ##    4:     O60885   BRD4           538         546
    ##    5:     O60885   BRD4           538         546
    ##   ---                                            
    ## 2672:     O60885   BRD4           511         539
    ## 2673:     O60885   BRD4           511         539
    ## 2674:     O60885   BRD4           511         539
    ## 2675:     O60885   BRD4           511         541
    ## 2676:     O60885   BRD4           511         541
    ##                               peptide
    ##                                <char>
    ##    1:                      RKEEVEENKK
    ##    2:                      RKEEVEENKK
    ##    3:                       KKEKDKKEK
    ##    4:                       KKEKDKKEK
    ##    5:                       KKEKDKKEK
    ##   ---                                
    ## 2672:   LAELQEQLKAVHEQLAALSQPQQNKPKKK
    ## 2673:   LAELQEQLKAVHEQLAALSQPQQNKPKKK
    ## 2674:   LAELQEQLKAVHEQLAALSQPQQNKPKKK
    ## 2675: LAELQEQLKAVHEQLAALSQPQQNKPKKKEK
    ## 2676: LAELQEQLKAVHEQLAALSQPQQNKPKKKEK
    ##                                                                                    mods
    ##                                                                                  <char>
    ##    1:                                                                                  
    ##    2:                                                                                  
    ##    3:                                1K(72.0211), 2K(15.9949), 4K(72.0211), 7K(72.0211)
    ##    4:                                             1K(72.0211), 4K(72.0211), 7K(72.0211)
    ##    5:                   1K(15.9949), 2K(56.0262), 4K(56.0262), 6K(72.0211), 7K(72.0211)
    ##   ---                                                                                  
    ## 2672:               25K(56.0262), 27K(56.0262), 28K(56.0262), 29K(56.0262), 9K(56.0262)
    ## 2673:               25K(56.0262), 27K(56.0262), 28K(56.0262), 29K(56.0262), 9K(56.0262)
    ## 2674:               25K(56.0262), 27K(56.0262), 28K(56.0262), 29K(56.0262), 9K(56.0262)
    ## 2675: 25K(56.0262), 27K(56.0262), 28K(56.0262), 29K(56.0262), 31K(56.0262), 9K(56.0262)
    ## 2676: 25K(56.0262), 27K(56.0262), 28K(56.0262), 29K(56.0262), 31K(56.0262), 9K(56.0262)
    ##       pep_mass mass_shift              variable intensity monoisotopic_mass
    ##          <num>      <num>                <fctr>     <num>            <fctr>
    ##    1: 1287.678    -0.0003 ox_100_0762_intensity         0          100.0762
    ##    2: 1287.678     0.0001 ox_100_0762_intensity         0          100.0762
    ##    3: 1391.751     0.0004 ox_100_0762_intensity         0          100.0762
    ##    4: 1375.756     0.0001 ox_100_0762_intensity         0          100.0762
    ##    5: 1375.756     0.0003 ox_100_0762_intensity         0          100.0762
    ##   ---                                                                      
    ## 2672: 3518.946    -0.0036 ox_156_1025_intensity         0          156.1025
    ## 2673: 3518.946     1.0012 ox_156_1025_intensity         0          156.1025
    ## 2674: 3518.946     0.0017 ox_156_1025_intensity         0          156.1025
    ## 2675: 3832.110     3.0114 ox_156_1025_intensity         0          156.1025
    ## 2676: 3832.110     0.9897 ox_156_1025_intensity         0          156.1025
    ##       diagnostic_ion
    ##               <lgcl>
    ##    1:          FALSE
    ##    2:          FALSE
    ##    3:          FALSE
    ##    4:          FALSE
    ##    5:          FALSE
    ##   ---               
    ## 2672:          FALSE
    ## 2673:          FALSE
    ## 2674:          FALSE
    ## 2675:          FALSE
    ## 2676:          FALSE

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
    ##  version  R version 4.5.1 (2025-06-13)
    ##  os       Ubuntu 24.04.2 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       Europe/Berlin
    ##  date     2025-11-25
    ##  pandoc   3.4 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.6.42 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  abind               1.4-8      2024-09-12 [1] CRAN (R 4.5.1)
    ##  backports           1.5.0      2024-05-23 [1] CRAN (R 4.5.1)
    ##  beeswarm            0.4.0      2021-06-01 [1] CRAN (R 4.5.1)
    ##  BiocGenerics      * 0.54.0     2025-04-15 [1] Bioconduc~
    ##  Biostrings        * 2.76.0     2025-04-15 [1] Bioconduc~
    ##  broom               1.0.10     2025-09-13 [1] CRAN (R 4.5.1)
    ##  car                 3.1-3      2024-09-27 [1] CRAN (R 4.5.1)
    ##  carData             3.0-5      2022-01-06 [1] CRAN (R 4.5.1)
    ##  cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.5.1)
    ##  cli                 3.6.5      2025-04-23 [1] CRAN (R 4.5.1)
    ##  crayon              1.5.3      2024-06-20 [1] CRAN (R 4.5.1)
    ##  data.table        * 1.17.8     2025-07-10 [1] CRAN (R 4.5.1)
    ##  digest              0.6.37     2024-08-19 [1] CRAN (R 4.5.1)
    ##  dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.5.1)
    ##  evaluate            1.0.5      2025-08-27 [1] CRAN (R 4.5.1)
    ##  farver              2.1.2      2024-05-13 [1] CRAN (R 4.5.1)
    ##  fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.5.1)
    ##  Formula             1.2-5      2023-02-24 [1] CRAN (R 4.5.1)
    ##  generics          * 0.1.4      2025-05-09 [1] CRAN (R 4.5.1)
    ##  GenomeInfoDb      * 1.44.3     2025-09-21 [1] Bioconduc~
    ##  GenomeInfoDbData    1.2.14     2025-09-24 [1] Bioconductor
    ##  ggbeeswarm        * 0.7.2      2023-04-29 [1] CRAN (R 4.5.1)
    ##  ggplot2           * 4.0.0      2025-09-11 [1] CRAN (R 4.5.1)
    ##  ggpubr            * 0.6.2      2025-10-17 [1] CRAN (R 4.5.1)
    ##  ggsignif            0.6.4      2022-10-13 [1] CRAN (R 4.5.1)
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
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.5.1)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.5.1)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-11-07 [1] local
    ##  purrr               1.1.0      2025-07-10 [1] CRAN (R 4.5.1)
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.5.1)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.5.1)
    ##  readxl            * 1.4.5      2025-03-07 [1] CRAN (R 4.5.1)
    ##  rlang               1.1.6      2025-04-11 [1] CRAN (R 4.5.1)
    ##  rmarkdown           2.29       2024-11-04 [1] CRAN (R 4.5.1)
    ##  rstatix             0.7.3      2025-10-18 [1] CRAN (R 4.5.1)
    ##  S4Vectors         * 0.46.0     2025-04-15 [1] Bioconduc~
    ##  S7                  0.2.0      2024-11-07 [1] CRAN (R 4.5.1)
    ##  scales              1.4.0      2025-04-24 [1] CRAN (R 4.5.1)
    ##  sessioninfo         1.2.3      2025-02-05 [1] CRAN (R 4.5.1)
    ##  snakecase           0.11.1     2023-08-27 [1] CRAN (R 4.5.1)
    ##  stringi             1.8.7      2025-03-27 [1] CRAN (R 4.5.1)
    ##  stringr           * 1.5.2      2025-09-08 [1] CRAN (R 4.5.1)
    ##  tibble              3.3.0      2025-06-08 [1] CRAN (R 4.5.1)
    ##  tidyr               1.3.1      2024-01-24 [1] CRAN (R 4.5.1)
    ##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.5.1)
    ##  timechange          0.3.0      2024-01-18 [1] CRAN (R 4.5.1)
    ##  UCSC.utils          1.4.0      2025-04-15 [1] Bioconduc~
    ##  vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.5.1)
    ##  vipor               0.4.7      2023-12-18 [1] CRAN (R 4.5.1)
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
