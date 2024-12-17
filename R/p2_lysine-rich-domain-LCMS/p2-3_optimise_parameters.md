2-3. Optimise MQ parameters
================
Yoichiro Sugimoto
17 December, 2024

- [Environment setup](#environment-setup)
- [Import basic data](#import-basic-data)
- [MS/MS count by different setting](#msms-count-by-different-setting)
- [QC by the position of propionylated
  lysines](#qc-by-the-position-of-propionylated-lysines)
- [The effect of MQ setting by K
  score](#the-effect-of-mq-setting-by-k-score)
- [Comparison of the results of PNAS paper (hydoxylation
  stoichometry)](#comparison-of-the-results-of-pnas-paper-hydoxylation-stoichometry)
- [Session information](#session-information)

This script calculates the stoichiometry of PTMs.

# Environment setup

``` r
# renv::init(
#           "/fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/R"
#       )

project.dir <-
  file.path(
    "/fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains"
  )

# renv::snapshot(file.path(project.dir, "R"))

renv::restore(file.path(project.dir, "R"))
```

    ## The following package(s) will be updated:
    ## 
    ## # https://bioconductor.org/packages/3.18/bioc --------------------------------
    ## - AnnotationDbi      [1.64.1 -> 1.64.1]
    ## - Biobase            [2.62.0 -> 2.62.0]
    ## - BiocGenerics       [0.48.1 -> 0.48.1]
    ## - BiocVersion        [3.18.1 -> 3.18.1]
    ## - Biostrings         [2.70.3 -> 2.70.3]
    ## - GenomeInfoDb       [1.38.8 -> 1.38.8]
    ## - IRanges            [2.36.0 -> 2.36.0]
    ## - KEGGREST           [1.42.0 -> 1.42.0]
    ## - S4Vectors          [0.40.2 -> 0.40.2]
    ## - XVector            [0.42.0 -> 0.42.0]
    ## - zlibbioc           [1.48.2 -> 1.48.2]
    ## 
    ## # https://bioconductor.org/packages/3.18/data/annotation ---------------------
    ## - GenomeInfoDbData   [1.2.11 -> 1.2.11]
    ## - org.Hs.eg.db       [3.18.0 -> 3.18.0]
    ## 
    ## # Installing packages --------------------------------------------------------
    ## - Installing BiocVersion ...                    OK [copied from cache]
    ## - Installing zlibbioc ...                       OK [copied from cache]
    ## - Installing BiocGenerics ...                   OK [copied from cache]
    ## - Installing GenomeInfoDbData ...               OK [copied from cache]
    ## - Installing S4Vectors ...                      OK [copied from cache]
    ## - Installing IRanges ...                        OK [copied from cache]
    ## - Installing GenomeInfoDb ...                   OK [copied from cache]
    ## - Installing XVector ...                        OK [copied from cache]
    ## - Installing Biobase ...                        OK [copied from cache]
    ## - Installing Biostrings ...                     OK [copied from cache]
    ## - Installing KEGGREST ...                       OK [copied from cache]
    ## - Installing AnnotationDbi ...                  OK [copied from cache]
    ## - Installing org.Hs.eg.db ...                   OK [copied from cache in 0.49s]

``` r
temp <-
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
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
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

    ## 
    ## Attaching package: 'janitor'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     chisq.test, fisher.test

``` r
# install.packages("/fast/AG_Sugimoto/home/users/yoichiro/projects/ptm.stoichiometry", repos = NULL, type = "source")
library("readxl")

p2.res.dir <- file.path(project.dir,
                        "results",
                        "p2-analysis-setting")
```

# Import basic data

``` r
pnas2022_data <- file.path(
  project.dir,
  "data/MQ_output/PNAS2022" 
)

all_sample_run_info <- read_excel(
  file.path(pnas2022_data, "PXD031221_sample_matrix.xlsx"),
  sheet = "run_setting"
) %>% data.table

data.dir <- file.path(project.dir, "data")
protein.feature.dt <- fread(file.path(data.dir, "PNAS2022/all_protein_feature_per_position.csv"))

setnames(protein.feature.dt, old = c("uniprot_id", "position"), new = c("protein_accession", "aa_pos"))
```

# MS/MS count by different setting

``` r
read_evidence <- function(prefix, dir_path){
  input.file <- file.path(
    dir_path, 
    paste0("including_SECPEP_", prefix, "all_processed_evidence_data.csv")
  )
  if(file.exists(input.file)){
    dt <- fread(input.file)
    dt[, condition := gsub("_$", "", prefix)]
  } else {
    dt <- data.table()
  }
  return(dt)
}

evidence.dt <- lapply(
  all_sample_run_info[data %in% c("data-A", "data-D") & !is.na(prefix), prefix],
  read_evidence,
  dir_path = p2.res.dir
  ) %>% rbindlist

data.count.dt <- evidence.dt[, list(
  total_peptide_count = .N,
  total_msms_count = sum(psm_mapped),
  total_intensity = sum(peak_intensity) 
), by = list(file_name, sample_name, condition, type)]

data.count.dt[, `:=`(
  propionylation = case_when(
    grepl("^data-D", condition) ~ "no_propionylation",
    TRUE ~ "propionylation"
  ),
  MQ_setting = str_extract(condition, "m\\d+_v\\d+"),
  protease_setting = str_split_fixed(condition, "_", n = 3)[, 2]
)]

data_col = c("total_peptide_count", "total_msms_count", "total_intensity")

m2_v2.dt <- data.count.dt[
  MQ_setting == "m2_v2" & type == "MULTI-MSMS"
  ][, c("file_name", data_col), with = FALSE]
setnames(m2_v2.dt, old = data_col, new = paste0("m2_v2_", data_col))

rel_count.dt <- data.count.dt

rel_count.dt <- merge(
  rel_count.dt,
  m2_v2.dt,
  by = c("file_name")
)

rel_count.dt[, `:=`(
  rel_total_peptide_count = total_peptide_count / m2_v2_total_peptide_count,
  rel_msms_count = total_msms_count / m2_v2_total_msms_count,
  rel_intensity = total_intensity / m2_v2_total_intensity 
)]

rel_count_stat.dt <- rel_count.dt[, list(
  meab_rel_total_peptide_count = mean(rel_total_peptide_count),
  sd_rel_total_peptide_count = sd(rel_total_peptide_count),
  mean_rel_msms_count = mean(rel_msms_count),
  sd_rel_msms_count = sd(rel_msms_count),
  mean_rel_intensity = mean(rel_intensity),
  sd_rel_intensity = sd(rel_intensity)
), by = list(propionylation, MQ_setting, protease_setting, type)]

rel_count_stat.dt[, `:=`(
  protease_setting = factor(protease_setting, levels = c("trp", "argC")),
  type = factor(type, levels = c("MULTI-MSMS", "MULTI-SECPEP")),
  stacked_mean_rel_msms_count = mean_rel_msms_count
)]

rel_count_stat.dt[type == "MULTI-SECPEP", `:=`(
  stacked_mean_rel_msms_count = 
    with(rel_count_stat.dt, mean_rel_msms_count[type == "MULTI-MSMS"]) +
    with(rel_count_stat.dt, mean_rel_msms_count[type == "MULTI-SECPEP"])
)]

ggplot(
  rel_count_stat.dt,
  aes(
    x = MQ_setting,
    y = mean_rel_msms_count,
    fill = type
  )
) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(
    ymin = mean_rel_msms_count - sd_rel_msms_count,
    ymax = mean_rel_msms_count + sd_rel_msms_count
  ), width = 0.5, position = position_dodge(width = 0.9)) +
  facet_grid(~ propionylation + protease_setting, scales = "free", space = "free") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  ylab("Relative MS/MS count") +
  scale_fill_manual(values = c("MULTI-MSMS" = "#4477AA", "MULTI-SECPEP" = "#66CCEE")) #+
```

![](p2-3_optimise_parameters_files/figure-gfm/MS_MS_count-1.png)<!-- -->

``` r
  # scale_color_manual(values = c("MULTI-MSMS" = "#4477AA", "MULTI-SECPEP" = "#66CCEE"))
```

The results indicate that the setting with 7 miscleavage and 7 variable
modifications identify the highest MS/MS count.

# QC by the position of propionylated lysines

``` r
read_all_per_pos_data_with_SECPEP <- function(prefix, dir_path){
  dt <- fread(file.path(dir_path, paste0("including_SECPEP_", prefix, "all_per_pos_data.csv")))
  dt[, condition := gsub("_$", "", prefix)]
  return(dt)
}

all_per_pos_with_SECPEP.dt <- lapply(
  all_sample_run_info[data %in% c("data-A") & grepl("m7_v7", prefix) & !is.na(prefix), prefix],
  read_all_per_pos_data_with_SECPEP,
  dir_path = p2.res.dir
  ) %>% rbindlist

all_per_pos_with_SECPEP.dt[, `:=`(
  propionylation = case_when(
    grepl("^data-D", condition) ~ "no_propionylation",
    TRUE ~ "propionylation"
  ),
  MQ_setting = str_extract(condition, "m\\d+_v\\d+"),
  protease_setting = str_split_fixed(condition, "_", n = 3)[, 2],
  is_peptide_C_term = peptide_end == aa_pos,
  is_propionylated = ptm %in% c("[Propionylation]", "[Oxidised Propionylation]")
)]

all_per_K_pos_with_SECPEP.dt <- all_per_pos_with_SECPEP.dt[aa == "K"]

propionylated_K_count.dt <- all_per_K_pos_with_SECPEP.dt[
  , list(total_MS_MS_count = sum(psm_mapped)),
  by = list(is_propionylated, is_peptide_C_term, type)
]

propionylated_K_count.dt
```

    ##    is_propionylated is_peptide_C_term         type total_MS_MS_count
    ##              <lgcl>            <lgcl>       <char>             <int>
    ## 1:            FALSE             FALSE MULTI-SECPEP               272
    ## 2:             TRUE             FALSE MULTI-SECPEP              1652
    ## 3:            FALSE              TRUE MULTI-SECPEP               273
    ## 4:             TRUE             FALSE   MULTI-MSMS            180191
    ## 5:            FALSE              TRUE   MULTI-MSMS             33088
    ## 6:            FALSE             FALSE   MULTI-MSMS             18788
    ## 7:             TRUE              TRUE   MULTI-MSMS              1124
    ## 8:             TRUE              TRUE MULTI-SECPEP                35

``` r
ggplot(
  propionylated_K_count.dt,
  aes(
    x = is_peptide_C_term,
    y = total_MS_MS_count,
    fill = is_propionylated
  )
) +
  geom_bar(stat = "identity", position="fill") +
  scale_fill_manual(values = c("TRUE" = "#4477AA", "FALSE" = "#BBBBBB")) +
  facet_grid(~ type) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
```

![](p2-3_optimise_parameters_files/figure-gfm/QC_by_propionylated_lysine_position-1.png)<!-- -->

Based on the analyses here, we decided to focus on using MULTI-MSMS
data.

# The effect of MQ setting by K score

``` r
read_stoic_data <- function(prefix, dir_path){
  dt <- fread(file.path(dir_path, paste0(prefix, "PTM_stoichiometry.csv")))
  dt[, condition := gsub("_$", "", prefix)]
  return(dt)
}

all_stoic_pos.dt <- lapply(
  all_sample_run_info[
    data != "data-D" & 
      grepl("_trp_", prefix) &
      !grepl("including_SECPEP", prefix) &
      grepl("m[257]_v[257]_(def|mCC)", prefix) & 
      !is.na(prefix), prefix
  ],
  read_stoic_data,
  dir_path = p2.res.dir
) %>% rbindlist

all_stoic_pos.dt[, `:=`(
  MQ_setting = str_extract(condition, "m\\d+_v\\d+")
)]

ms_ms_count.dt <- all_stoic_pos.dt[
  , list(total_ms_ms_count = sum(sum_psm_mapped)), by = list(protein_accession, aa_pos, MQ_setting)
]

ms_ms_count.dt <- merge(
  ms_ms_count.dt,
  protein.feature.dt[, .(protein_accession, aa_pos, K_ratio, K_ratio_score)],
  by = c("protein_accession", "aa_pos")
)

d.ms_ms_count.dt <- dcast(
  ms_ms_count.dt,
  protein_accession + aa_pos + K_ratio + K_ratio_score ~ MQ_setting,
  value.var = "total_ms_ms_count", fill = 0
  
)

library("RColorBrewer")
library("patchwork")

d.ms_ms_count.dt[, table(K_ratio)]
```

    ## K_ratio
    ##      0    0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9 
    ## 281712 188153  76269  24156   7014   2074    604    221     60     10

``` r
p1 <- ggplot(
  d.ms_ms_count.dt[order(K_ratio)],
  aes(
    x = m2_v2,
    y = m5_v5
  )
) +
  geom_abline(slope = 1, intercept = 0, color = "gray60") +
  geom_point() +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1000)) +
  theme(aspect.ratio = 1)

p2 <- ggplot(
  d.ms_ms_count.dt[order(K_ratio)],
  aes(
    x = m5_v5,
    y = m7_v7
  )
) +
  geom_abline(slope = 1, intercept = 0, color = "gray60") +
  geom_point() +
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 1000)) +
  theme(aspect.ratio = 1)

p1 + p2
```

![](p2-3_optimise_parameters_files/figure-gfm/high_confident_PTM_sites-1.png)<!-- -->

``` r
p3 <- ggplot(
  d.ms_ms_count.dt,
  aes(
    x = factor(K_ratio),
    y = m5_v5 - m2_v2
  )
) +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_boxplot(outlier.shape = NA, fill = "steelblue") +
  coord_cartesian(ylim = c(-15, 25)) +
  theme(aspect.ratio = 1)

p4 <- ggplot(
  d.ms_ms_count.dt,
  aes(
    x = factor(K_ratio),
    y = m7_v7 - m5_v5
  )
) +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_boxplot(outlier.shape = NA, fill = "steelblue") +
  coord_cartesian(ylim = c(-15, 25))

p3 + p4
```

![](p2-3_optimise_parameters_files/figure-gfm/high_confident_PTM_sites-2.png)<!-- -->

# Comparison of the results of PNAS paper (hydoxylation stoichometry)

``` r
all_stoic_pos.dt[, table(sample_name)]
```

    ## sample_name
    ##              FLAG_HeLaWT_derivatised JMJD6peptide_HeLaJMJD6KO_derivatised 
    ##                               438940                              1285413 
    ##      JMJD6peptide_HeLaWT_derivatised          JQ1_HeLaJMJD6KO_derivatised 
    ##                              1014522                               325618 
    ##               JQ1_HeLaWT_derivatised 
    ##                               318286

``` r
wt_ko_comparison.dt <- all_stoic_pos.dt[
  sample_name %in% c(
    "JMJD6peptide_HeLaWT_derivatised",
    "JMJD6peptide_HeLaJMJD6KO_derivatised",
    "JQ1_HeLaWT_derivatised",
    "JQ1_HeLaJMJD6KO_derivatised"
  )
]

wt_ko_comparison.dt[, `:=`(
  JMJD6_GT = str_extract(sample_name, "(?<=HeLa)(WT|JMJD6KO)")
)]


zero.stoic.dt <- wt_ko_comparison.dt[aa == "K" & grepl("[Oxidation (K)]", ptm, fixed = TRUE)]


zero.stoic.dt[, `:=`(
  sum_psm_mapped = 0,
  stoichiometry = 0,
  ptm = "[Oxidation (K)]"
)]

wt_ko_comparison.dt[aa == "K" & grepl("[Oxidation (K)]", ptm, fixed = TRUE)]
```

    ##                                sample_name protein_accession gene_name     aa
    ##                                     <char>            <char>    <char> <char>
    ##    1:               JQ1_HeLaWT_derivatised            A6NIE6    RRN3P2      K
    ##    2:               JQ1_HeLaWT_derivatised            A7E2F4   GOLGA8A      K
    ##    3:               JQ1_HeLaWT_derivatised            A7E2F4   GOLGA8A      K
    ##    4:          JQ1_HeLaJMJD6KO_derivatised            B1AJZ9     FHAD1      K
    ##    5:          JQ1_HeLaJMJD6KO_derivatised            B1AJZ9     FHAD1      K
    ##   ---                                                                        
    ## 9192: JMJD6peptide_HeLaJMJD6KO_derivatised            Q9Y5Q9    GTF3C3      K
    ## 9193: JMJD6peptide_HeLaJMJD6KO_derivatised            Q9Y608   LRRFIP2      K
    ## 9194: JMJD6peptide_HeLaJMJD6KO_derivatised            Q9Y608   LRRFIP2      K
    ## 9195: JMJD6peptide_HeLaJMJD6KO_derivatised            Q9Y6J9     TAF6L      K
    ## 9196:      JMJD6peptide_HeLaWT_derivatised            Q9Y6J9     TAF6L      K
    ##       aa_pos             ptm sum_peak_intensity sum_psm_mapped
    ##        <int>          <char>              <num>          <int>
    ##    1:     62 [Oxidation (K)]          385580000              2
    ##    2:    121 [Oxidation (K)]          778810000              1
    ##    3:    124 [Oxidation (K)]          778810000              1
    ##    4:    734 [Oxidation (K)]          635280000              1
    ##    5:    738 [Oxidation (K)]          635280000              1
    ##   ---                                                         
    ## 9192:    134 [Oxidation (K)]            9580000              1
    ## 9193:    692 [Oxidation (K)]           90984000              1
    ## 9194:    705 [Oxidation (K)]            3550000              1
    ## 9195:    484 [Oxidation (K)]           12774000              2
    ## 9196:    484 [Oxidation (K)]            5141600              1
    ##       the_number_of_peptide max_score max_ptm_probability
    ##                       <int>     <num>               <num>
    ##    1:                     1    70.622               0.969
    ##    2:                     1    95.502               0.667
    ##    3:                     1    95.502               0.667
    ##    4:                     1   109.860               1.000
    ##    5:                     1   109.860               1.000
    ##   ---                                                    
    ## 9192:                     1    47.288               0.291
    ## 9193:                     1    46.663               0.910
    ## 9194:                     1    56.121               0.711
    ## 9195:                     2    85.737               1.000
    ## 9196:                     1    62.362               1.000
    ##       sum_intensity_per_position sum_psm_mapped_per_position
    ##                            <num>                       <int>
    ##    1:                  385580000                           2
    ##    2:                  778810000                           1
    ##    3:                  778810000                           1
    ##    4:                  635280000                           1
    ##    5:                  635280000                           1
    ##   ---                                                       
    ## 9192:                   28228000                           4
    ## 9193:                  339284200                          17
    ## 9194:                  339284200                          17
    ## 9195:                   12774000                           2
    ## 9196:                    5141600                           1
    ##       sum_the_number_of_peptide stoichiometry            condition MQ_setting
    ##                           <int>         <num>               <char>     <char>
    ##    1:                         1     1.0000000 data-A_trp_m2_v2_def      m2_v2
    ##    2:                         1     1.0000000 data-A_trp_m2_v2_def      m2_v2
    ##    3:                         1     1.0000000 data-A_trp_m2_v2_def      m2_v2
    ##    4:                         1     1.0000000 data-A_trp_m2_v2_def      m2_v2
    ##    5:                         1     1.0000000 data-A_trp_m2_v2_def      m2_v2
    ##   ---                                                                        
    ## 9192:                         3     0.3393793 data-B_trp_m7_v7_mCC      m7_v7
    ## 9193:                        15     0.2681646 data-B_trp_m7_v7_mCC      m7_v7
    ## 9194:                        15     0.0104632 data-B_trp_m7_v7_mCC      m7_v7
    ## 9195:                         2     1.0000000 data-B_trp_m7_v7_mCC      m7_v7
    ## 9196:                         1     1.0000000 data-B_trp_m7_v7_mCC      m7_v7
    ##       JMJD6_GT
    ##         <char>
    ##    1:       WT
    ##    2:       WT
    ##    3:       WT
    ##    4:  JMJD6KO
    ##    5:  JMJD6KO
    ##   ---         
    ## 9192:  JMJD6KO
    ## 9193:  JMJD6KO
    ## 9194:  JMJD6KO
    ## 9195:  JMJD6KO
    ## 9196:       WT

``` r
## Hydroxylation stoichiometry comparisons
pnas2022.stoic.dt <- fread(file.path(data.dir, "PNAS2022/long_K_stoichiometry_data.csv"))
setnames(pnas2022.stoic.dt, old = c("uniprot_id", "position"), new = c("protein_accession", "aa_pos"))

pnas2022.stoic.dt[, `:=`(
  accession_position = paste0(protein_accession, "_", aa_pos)
)]

non.duplicated.pnas2022.stoic.dt <- pnas2022.stoic.dt[
  (data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep")) &
    (curated_oxK_site == TRUE) &
    total_n_feature_K > 2
][
  order(
    data_source %in% c("HeLa_WT_JQ1", "HeLa_WT_J6pep"),
    curated_oxK_site == TRUE,
    oxK_ratio,
    decreasing = TRUE
  )][
    !duplicated(accession_position)
]

selected_stoic_pos.dt <- all_stoic_pos.dt[
  sample_name %in% c("JQ1_HeLaWT_derivatised", "JMJD6peptide_HeLaWT_derivatised") &
    aa == "K" & grepl("[Oxidation (K)]", ptm, fixed = TRUE)
][order(stoichiometry, decreasing = TRUE)][!duplicated(paste0(protein_accession, "_", aa_pos))]

pnas2022.mq1.stoic.dt <- merge(
  non.duplicated.pnas2022.stoic.dt,
  selected_stoic_pos.dt[, .(protein_accession, aa_pos, stoichiometry)],
  by = c("protein_accession", "aa_pos")
)

pnas2022.mq1.stoic.dt[, cor.test(oxK_ratio, stoichiometry)]
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  oxK_ratio and stoichiometry
    ## t = 17.957, df = 48, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.8842532 0.9615994
    ## sample estimates:
    ##       cor 
    ## 0.9329693

``` r
ggplot(
  pnas2022.mq1.stoic.dt,
  aes(
    x = oxK_ratio,
    y = stoichiometry
  )
) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Stoichiometry in PNAS2022") +
  ylab("Stoichiometry in new workflow")
```

![](p2-3_optimise_parameters_files/figure-gfm/comparison_with_PNAS_paper-1.png)<!-- -->

# Session information

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.3.2 (2023-10-31)
    ##  os       Ubuntu 22.04.3 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       Europe/Berlin
    ##  date     2024-12-17
    ##  pandoc   3.1.1 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package          * version   date (UTC) lib source
    ##  BiocGenerics     * 0.48.1    2023-11-01 [1] Bioconductor
    ##  BiocManager        1.30.25   2024-08-28 [1] CRAN (R 4.3.2)
    ##  Biostrings       * 2.70.3    2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
    ##  bitops             1.0-9     2024-10-03 [1] CRAN (R 4.3.2)
    ##  cellranger         1.1.0     2016-07-27 [1] CRAN (R 4.3.2)
    ##  cli                3.6.3     2024-06-21 [1] CRAN (R 4.3.2)
    ##  colorspace         2.1-1     2024-07-26 [1] CRAN (R 4.3.2)
    ##  crayon             1.5.3     2024-06-20 [1] CRAN (R 4.3.2)
    ##  data.table       * 1.15.4    2024-03-30 [1] CRAN (R 4.3.2)
    ##  digest             0.6.36    2024-06-23 [1] CRAN (R 4.3.2)
    ##  dplyr            * 1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
    ##  evaluate           0.24.0    2024-06-10 [1] CRAN (R 4.3.2)
    ##  fansi              1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
    ##  farver             2.1.2     2024-05-13 [1] CRAN (R 4.3.2)
    ##  fastmap            1.2.0     2024-05-15 [1] CRAN (R 4.3.2)
    ##  generics           0.1.3     2022-07-05 [1] CRAN (R 4.3.2)
    ##  GenomeInfoDb     * 1.38.8    2024-03-15 [1] Bioconductor 3.18 (R 4.3.2)
    ##  GenomeInfoDbData   1.2.11    2024-11-18 [1] Bioconductor
    ##  ggplot2          * 3.5.1     2024-04-23 [1] CRAN (R 4.3.2)
    ##  glue               1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
    ##  gtable             0.3.5     2024-04-22 [1] CRAN (R 4.3.2)
    ##  highr              0.11      2024-05-26 [1] CRAN (R 4.3.2)
    ##  htmltools          0.5.8.1   2024-04-04 [1] CRAN (R 4.3.2)
    ##  IRanges          * 2.36.0    2023-10-24 [1] Bioconductor
    ##  janitor          * 2.2.0     2023-02-02 [1] CRAN (R 4.3.2)
    ##  khroma           * 1.14.0    2024-08-26 [1] CRAN (R 4.3.2)
    ##  knitr            * 1.48      2024-07-07 [1] CRAN (R 4.3.2)
    ##  labeling           0.4.3     2023-08-29 [1] CRAN (R 4.3.2)
    ##  lifecycle          1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
    ##  lubridate          1.9.3     2023-09-27 [1] CRAN (R 4.3.2)
    ##  magrittr         * 2.0.3     2022-03-30 [1] CRAN (R 4.3.2)
    ##  munsell            0.5.1     2024-04-01 [1] CRAN (R 4.3.2)
    ##  patchwork        * 1.3.0     2024-09-16 [1] CRAN (R 4.3.2)
    ##  pillar             1.9.0     2023-03-22 [1] CRAN (R 4.3.2)
    ##  pkgconfig          2.0.3     2019-09-22 [1] CRAN (R 4.3.2)
    ##  R6                 2.5.1     2021-08-19 [1] CRAN (R 4.3.2)
    ##  RColorBrewer     * 1.1-3     2022-04-03 [1] CRAN (R 4.3.2)
    ##  RCurl              1.98-1.16 2024-07-11 [1] CRAN (R 4.3.2)
    ##  readxl           * 1.4.3     2023-07-06 [1] CRAN (R 4.3.2)
    ##  renv               1.0.7     2024-04-11 [1] CRAN (R 4.3.2)
    ##  rlang              1.1.4     2024-06-04 [1] CRAN (R 4.3.2)
    ##  rmarkdown          2.27      2024-05-17 [1] CRAN (R 4.3.2)
    ##  rstudioapi         0.17.1    2024-10-22 [1] CRAN (R 4.3.2)
    ##  S4Vectors        * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
    ##  scales           * 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
    ##  sessioninfo        1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
    ##  snakecase          0.11.1    2023-08-27 [1] CRAN (R 4.3.2)
    ##  stringi            1.8.4     2024-05-06 [1] CRAN (R 4.3.2)
    ##  stringr          * 1.5.1     2023-11-14 [1] CRAN (R 4.3.2)
    ##  tibble             3.2.1     2023-03-20 [1] CRAN (R 4.3.2)
    ##  tidyselect         1.2.1     2024-03-11 [1] CRAN (R 4.3.2)
    ##  timechange         0.3.0     2024-01-18 [1] CRAN (R 4.3.2)
    ##  utf8               1.2.4     2023-10-22 [1] CRAN (R 4.3.2)
    ##  vctrs              0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
    ##  withr              3.0.1     2024-07-31 [1] CRAN (R 4.3.2)
    ##  xfun               0.46      2024-07-18 [1] CRAN (R 4.3.2)
    ##  XVector          * 0.42.0    2023-10-24 [1] Bioconductor
    ##  yaml               2.3.10    2024-07-26 [1] CRAN (R 4.3.2)
    ##  zlibbioc           1.48.2    2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
    ## 
    ##  [1] /home/ysugimo/R/x86_64-pc-linux-gnu-library/4.3
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
