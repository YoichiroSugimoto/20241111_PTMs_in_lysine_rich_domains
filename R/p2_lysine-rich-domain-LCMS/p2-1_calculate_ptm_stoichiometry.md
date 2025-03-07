2-1. Calculate PTM stoichiometry
================
Yoichiro Sugimoto
07 March, 2025

- [Environment setup](#environment-setup)
- [Import basic data](#import-basic-data)
- [Calculation of stoichiometry for data without diagnostic ion
  consideration (for optimisation of MQ run
  parameters)](#calculation-of-stoichiometry-for-data-without-diagnostic-ion-consideration-for-optimisation-of-mq-run-parameters)
- [Calculation of stoichiometry for data with diagnostic ion
  consideration](#calculation-of-stoichiometry-for-data-with-diagnostic-ion-consideration)
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
    ## - Installing BiocGenerics ...                   OK [copied from cache]
    ## - Installing Biobase ...                        OK [copied from cache]
    ## - Installing S4Vectors ...                      OK [copied from cache]
    ## - Installing IRanges ...                        OK [copied from cache]
    ## - Installing zlibbioc ...                       OK [copied from cache]
    ## - Installing XVector ...                        OK [copied from cache]
    ## - Installing GenomeInfoDbData ...               OK [copied from cache]
    ## - Installing GenomeInfoDb ...                   OK [copied from cache]
    ## - Installing Biostrings ...                     OK [copied from cache]
    ## - Installing KEGGREST ...                       OK [copied from cache]
    ## - Installing AnnotationDbi ...                  OK [copied from cache]
    ## - Installing org.Hs.eg.db ...                   OK [copied from cache in 0.57s]

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

``` r
# install.packages("/fast/AG_Sugimoto/home/users/yoichiro/projects/ptm.stoichiometry", repos = NULL, type = "source")
library("readxl")
library("ptm.stoichiometry")

p2.res.dir <- file.path(project.dir,
                        "results",
                        "p2-analysis-setting")

create.dirs(c(p2.res.dir))
```

# Import basic data

``` r
ref_protein_dt <- import_reference_fasta(
  file.path(
    "/fast/AG_Sugimoto/reference/uniprot/human",
    "UP000005640_9606.fasta"
  )
)
```

# Calculation of stoichiometry for data without diagnostic ion consideration (for optimisation of MQ run parameters)

``` r
pnas2022_data <- file.path(
  project.dir,
  "data/MQ_output/PNAS2022" 
)

all_sample_run_info <- read_excel(
  file.path(pnas2022_data, "PXD031221_sample_matrix.xlsx"),
  sheet = "run_setting"
) %>% data.table

all_sample_run_info[, sample_id := 1:.N]
```

``` r
## First, analyse with both MULTI-MSMS and MULTI-SECPEP
for(i in all_sample_run_info[data %in% c("data-A", "data-D"), sample_id]) {

  mq_evidence_data = file.path(pnas2022_data,
                               "evidence",
                               paste0(all_sample_run_info[sample_id == i, prefix], "evidence.txt"))
  
  print(paste0("Processing: ", basename(mq_evidence_data)))
  
  if (file.exists(mq_evidence_data)) {
    stoic.dt <- calculate_stoichiometry(
      mq_evidence_data = mq_evidence_data,
      sample_info_file = file.path(
        pnas2022_data,
        "sample_info",
        paste0(
          "MS_dataset_overview_PXD031221_",
          all_sample_run_info[sample_id == i, data],
          ".csv"
        )
      ),
      ref_protein_dt = ref_protein_dt,
      output_prefix = file.path(p2.res.dir, paste0(
        "including_SECPEP_", all_sample_run_info[sample_id == i, prefix] 
        )),
      ptm_mapping_file = file.path(
        project.dir, 
        "data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv"
      ),
      K_only = FALSE,
      selected_type = NA
    )
  } else {
    print(paste0("File does not exist: ", mq_evidence_data))
  }
  
  gc()
}
```

    ## [1] "Processing: data-A_trp_m2_v2_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m3_v3_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m4_v4_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m5_v5_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m6_v6_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m7_v7_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m8_v8_def_evidence.txt"
    ## [1] "Processing: data-A_argC_m2_v7_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-D_trp_m2_v2_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-D_trp_m5_v5_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-D_trp_m7_v7_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"

``` r
for(i in 1:nrow(all_sample_run_info)) {

  mq_evidence_data = file.path(pnas2022_data,
                               "evidence",
                               paste0(all_sample_run_info[i, prefix], "evidence.txt"))
  
  print(paste0("Processing: ", basename(mq_evidence_data)))
  
  if (file.exists(mq_evidence_data)) {
    stoic.dt <- calculate_stoichiometry(
      mq_evidence_data = mq_evidence_data,
      sample_info_file = file.path(
        pnas2022_data,
        "sample_info",
        paste0(
          "MS_dataset_overview_PXD031221_",
          all_sample_run_info[i, data],
          ".csv"
        )
      ),
      ref_protein_dt = ref_protein_dt,
      output_prefix = file.path(p2.res.dir, all_sample_run_info[i, prefix]),
      ptm_mapping_file = file.path(
        project.dir, 
        "data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv"
      ),
      K_only = FALSE,
      selected_type = c("MULTI-MSMS")
    )
  } else {
    print(paste0("File does not exist: ", mq_evidence_data))
  }
  
  gc()
}
```

    ## [1] "Processing: data-A_trp_m2_v2_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m3_v3_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m4_v4_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m5_v5_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m6_v6_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m7_v7_def_evidence.txt"
    ## [1] "Processing: data-A_trp_m8_v8_def_evidence.txt"
    ## [1] "Processing: data-A_argC_m2_v7_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-B_trp_m2_v2_mCC_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-B_trp_m5_v5_mCC_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-B_trp_m7_v7_mCC_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-C_trp_m2_v2_mCC_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-C_trp_m5_v5_mCC_evidence.txt"
    ## [1] "Processing: data-C_trp_m6_v6_mCC_evidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/data-C_trp_m6_v6_mCC_evidence.txt"
    ## [1] "Processing: data-C_trp_m7_v7_mCC_evidence.txt"
    ## [1] "Processing: data-C_trp_m8_v8_mCC_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-D_trp_m2_v2_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-D_trp_m5_v5_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: data-D_trp_m7_v7_def_evidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"
    ## [1] "Processing: NAevidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_output/PNAS2022/evidence/NAevidence.txt"

# Calculation of stoichiometry for data with diagnostic ion consideration

``` r
pnas2022_DI_data <- file.path(
  project.dir,
  "data/MQ_DI_output/PNAS2022" 
)
```

``` r
di_sample_run_info <- all_sample_run_info[grepl("data-[A-C]_trp_m7_v7_", prefix)]

for(i in 1:nrow(di_sample_run_info)) {

  mq_evidence_data = file.path(pnas2022_DI_data,
                               "evidence",
                               paste0(di_sample_run_info[i, prefix], "evidence.txt"))
  
  print(paste0("Processing: ", basename(mq_evidence_data)))
  
  if (file.exists(mq_evidence_data)) {
    
    ptm_files <- list.files(file.path(
      pnas2022_DI_data, "ptm", 
      paste0(di_sample_run_info[i, prefix], "ptm")
    ), full.names = TRUE)
    
    ptm_names <- paste0(
      "[", str_replace_all(basename(ptm_files), "Sites.txt", ""), "]"
    )
    
    stoic.dt <- calculate_stoichiometry2(
      mq_evidence_data = mq_evidence_data,
      sample_info_file = file.path(
        pnas2022_data,
        "sample_info",
        paste0(
          "MS_dataset_overview_PXD031221_",
          di_sample_run_info[i, data],
          ".csv"
        )
      ),
      ref_protein_dt = ref_protein_dt,
      ptm_files = ptm_files,
      ptm_names = ptm_names,
      output_prefix = file.path(p2.res.dir, paste0("DI_", di_sample_run_info[i, prefix])),
      ptm_mapping_file = file.path(
        project.dir, 
        "data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv"
      ),
      K_only = FALSE,
      selected_type = c("MULTI-MSMS")
    )
  } else {
    print(paste0("File does not exist: ", mq_evidence_data))
  }
  
  gc()
}
```

    ## [1] "Processing: data-A_trp_m7_v7_def_evidence.txt"
    ## [1] "Processing: data-B_trp_m7_v7_mCC_evidence.txt"
    ## [1] "File does not exist: /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/data/MQ_DI_output/PNAS2022/evidence/data-B_trp_m7_v7_mCC_evidence.txt"
    ## [1] "Processing: data-C_trp_m7_v7_mCC_evidence.txt"

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
    ##  date     2025-03-07
    ##  pandoc   3.1.1 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  BiocGenerics      * 0.48.1     2023-11-01 [1] Bioconductor
    ##  BiocManager         1.30.25    2024-08-28 [1] CRAN (R 4.3.2)
    ##  Biostrings        * 2.70.3     2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
    ##  bit                 4.5.0      2024-09-20 [1] CRAN (R 4.3.2)
    ##  bit64               4.5.2      2024-09-22 [1] CRAN (R 4.3.2)
    ##  bitops              1.0-9      2024-10-03 [1] CRAN (R 4.3.2)
    ##  cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.3.2)
    ##  cli                 3.6.3      2024-06-21 [1] CRAN (R 4.3.2)
    ##  colorspace          2.1-1      2024-07-26 [1] CRAN (R 4.3.2)
    ##  crayon              1.5.3      2024-06-20 [1] CRAN (R 4.3.2)
    ##  data.table        * 1.15.4     2024-03-30 [1] CRAN (R 4.3.2)
    ##  digest              0.6.36     2024-06-23 [1] CRAN (R 4.3.2)
    ##  dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
    ##  evaluate            0.24.0     2024-06-10 [1] CRAN (R 4.3.2)
    ##  fansi               1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
    ##  fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.3.2)
    ##  generics            0.1.3      2022-07-05 [1] CRAN (R 4.3.2)
    ##  GenomeInfoDb      * 1.38.8     2024-03-15 [1] Bioconductor 3.18 (R 4.3.2)
    ##  GenomeInfoDbData    1.2.11     2024-11-18 [1] Bioconductor
    ##  ggplot2           * 3.5.1      2024-04-23 [1] CRAN (R 4.3.2)
    ##  glue                1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
    ##  gtable              0.3.5      2024-04-22 [1] CRAN (R 4.3.2)
    ##  htmltools           0.5.8.1    2024-04-04 [1] CRAN (R 4.3.2)
    ##  IRanges           * 2.36.0     2023-10-24 [1] Bioconductor
    ##  janitor             2.2.0      2023-02-02 [1] CRAN (R 4.3.2)
    ##  khroma            * 1.14.0     2024-08-26 [1] CRAN (R 4.3.2)
    ##  knitr             * 1.48       2024-07-07 [1] CRAN (R 4.3.2)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
    ##  lubridate           1.9.3      2023-09-27 [1] CRAN (R 4.3.2)
    ##  magrittr          * 2.0.3      2022-03-30 [1] CRAN (R 4.3.2)
    ##  munsell             0.5.1      2024-04-01 [1] CRAN (R 4.3.2)
    ##  pillar              1.9.0      2023-03-22 [1] CRAN (R 4.3.2)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.3.2)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-03-03 [1] local
    ##  R6                  2.5.1      2021-08-19 [1] CRAN (R 4.3.2)
    ##  RCurl               1.98-1.16  2024-07-11 [1] CRAN (R 4.3.2)
    ##  readxl            * 1.4.3      2023-07-06 [1] CRAN (R 4.3.2)
    ##  renv                1.0.7      2024-04-11 [1] CRAN (R 4.3.2)
    ##  rlang               1.1.4      2024-06-04 [1] CRAN (R 4.3.2)
    ##  rmarkdown           2.27       2024-05-17 [1] CRAN (R 4.3.2)
    ##  rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.3.2)
    ##  S4Vectors         * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
    ##  scales              1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
    ##  sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.3.2)
    ##  snakecase           0.11.1     2023-08-27 [1] CRAN (R 4.3.2)
    ##  stringi             1.8.4      2024-05-06 [1] CRAN (R 4.3.2)
    ##  stringr           * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
    ##  tibble              3.2.1      2023-03-20 [1] CRAN (R 4.3.2)
    ##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.3.2)
    ##  timechange          0.3.0      2024-01-18 [1] CRAN (R 4.3.2)
    ##  utf8                1.2.4      2023-10-22 [1] CRAN (R 4.3.2)
    ##  vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
    ##  withr               3.0.1      2024-07-31 [1] CRAN (R 4.3.2)
    ##  xfun                0.46       2024-07-18 [1] CRAN (R 4.3.2)
    ##  XVector           * 0.42.0     2023-10-24 [1] Bioconductor
    ##  yaml                2.3.10     2024-07-26 [1] CRAN (R 4.3.2)
    ##  zlibbioc            1.48.2     2024-03-13 [1] Bioconductor 3.18 (R 4.3.2)
    ## 
    ##  [1] /home/ysugimo/R/x86_64-pc-linux-gnu-library/4.3
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
