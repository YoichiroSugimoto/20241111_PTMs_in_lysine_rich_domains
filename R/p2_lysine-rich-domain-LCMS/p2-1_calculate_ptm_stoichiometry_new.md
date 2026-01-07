2-1. Calculate PTM stoichiometry
================
Yoichiro Sugimoto and Pallavi Kesavan
07 January, 2026

- [Overview](#overview)
- [This script calculates the stiochiomtery using MaxQuant
  data](#this-script-calculates-the-stiochiomtery-using-maxquant-data)
- [Environment setup](#environment-setup)
- [2.1.1 Install,load essential functions and
  libraries](#211-installload-essential-functions-and-libraries)
- [2.1.2 Import human protein reference
  data](#212-import-human-protein-reference-data)
- [2.1.3 Calculation of Stoichiometry without diagnostic ion
  consideration](#213-calculation-of-stoichiometry-without-diagnostic-ion-consideration)
- [2.1.4 Calculation of Stoichiometry of MQ_DI without
  waterloss](#214-calculation-of-stoichiometry-of-mq_di-without-waterloss)
- [Session information](#session-information)

# Overview

# This script calculates the stiochiomtery using MaxQuant data

# Environment setup

``` r
## Initialize renv (first time only) - re-installed 24.09.2025
# Creates project specific library and renv.lock file. 
# Use 'renv::init(filepath)' to create project library
# renv::init(
#        "/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains/R"
#    )

# Define project directory - contains the R scripts, data and result folders
project.dir <-
  file.path(
    "/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains"
  )

# renv::snapshot(file.path(project.dir, "R"))

#renv::restore completed 24.09.2025
#renv::restore(file.path(project.dir, "R"))
```

# 2.1.1 Install,load essential functions and libraries

``` r
## Load all R scripts from the 'functions' folder into the current session
P2_functions <- sapply(list.files(file.path(project.dir, "R/functions"), pattern="*.R", full.names = TRUE), source)
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
## Install private package 
# Install ptm.stiochiometry package
#install.packages("/fast/AG_Sugimoto/home/users/pallavi/projects/ptm.stoichiometry", repos = NULL, type = "source")

# Load Libraries - ptm.stiochiometry and readxl
library(ptm.stoichiometry)
library("readxl")
```

``` r
# Define path to data directory
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

### To be deleted
results.dir <- file.path("/fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains/results")
```

# 2.1.2 Import human protein reference data

``` r
# Import human protein reference data from specified file path 
ref_protein_data <- import_reference_fasta(file.path
                                           ("/fast/AG_Sugimoto/reference/uniprot/human", 
                                             "UP000005640_9606.fasta")) 

all.protein.bs <- Biostrings::readAAStringSet(
  file.path(
    "/fast/AG_Sugimoto/reference/uniprot/human",
    "UP000005640_9606.fasta"
  )
)
```

# 2.1.3 Calculation of Stoichiometry without diagnostic ion consideration

``` r
# Load data into environment
pnas2022_data <- file.path(
  project.dir,
  "data/MQ_standard/PNAS2022" 
)

# Load sample info data 
MQ_standard_data <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_standard"
) %>% data.table

# Add a new column 'sample_id' that assigns a unique sequential number to each row (1:N), 
#where N is the total number of rows in the data.table
MQ_standard_data[, sample_id := 1:.N]
```

``` r
## Stoichiomtery analysis using with MULTI-MSMS with 
# MULTI-SECPEP (secondary peptide) and MULTI-MSMS only

# Automate stoichiometry calculations across MaxQuant datasets (MULTI-MSMS + 
# MULTI-SECPEP and MULTI-MSMS only)
process_stoichiometry <- function(MQ_standard_data,
                                  output_dir,
                                  ref_protein_data,
                                  ptm_mapping_file,
                                  selected_type = NA,
                                  prefix_label = "") {
  
  # Create output directory if it does not already exist 
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # For each sample_id in the MQ_standard data, calculate the stoichiometry
  for (i in seq_len(nrow(MQ_standard_data))) {
    prefix    <- MQ_standard_data$prefix[i] # Extract the prefix column
    MQ_standard_sample <- MQ_standard_data$data[i] # Extract the data column 
    
    # Define file path to MaxQuant evidence files
    mq_evidence_data <- file.path(pnas2022_data, 
                                  "evidence", 
                                  paste0(prefix, "evidence.txt"))
    
    # Print which sample id is being processed
    message("Processing: ", basename(mq_evidence_data)) 
    
    
    # Checks whether the sample ID has the corresponding evidence.txt file. If yes, then proceed
    if (file.exists(mq_evidence_data)) {
      
      # Define corresponding sample information (metadata) file
      sample_info_file <- file.path(
        pnas2022_data, "sample_info",
        paste0("MS_dataset_overview_PXD031221_", MQ_standard_sample, ".csv")
      )
      
      # Define output prefix for results
      output_prefix <- file.path(output_dir, paste0(prefix_label, prefix))
      
      # Run main stoichiometry calculation
      calculate_stoichiometry(
        mq_evidence_data = mq_evidence_data,
        sample_info_file = sample_info_file,
        ref_protein_dt   = ref_protein_data,
        output_prefix    = output_prefix,
        ptm_mapping_file = ptm_mapping_file,
        K_only           = FALSE,
        selected_type    = selected_type
      )
      
    } else {
      # If evidence.txt does not exist, skip the function and print below message
      message("File does not exist: ", mq_evidence_data)
    }
    # Garbage collection - to free up memory 
    gc()
  }
}


# Define file path common PTM mapping file
ptm_mapping_file <- file.path(
  project.dir,
  "data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv"
)

## MSMS + SECEP

# Subset data-A and data-D from MQ_standard_data
samples_to_process <- MQ_standard_data[data %in% c("data-A", "data-D")]

# Define output directory for MSMS + SECPEP results
MQ_Std_MSMS_SECEP_dir <- file.path(results.dir, "p2-analysis-setting","MQ_Std_MSMS_SECEP")

# Run stoichiometry processing for MSMS + SECPEP
process_stoichiometry(
  MQ_standard_data = samples_to_process,
  output_dir       = MQ_Std_MSMS_SECEP_dir,
  ref_protein_data = ref_protein_data,
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
##MSMS only

# Define output directory for MSMS-only results
MQ_Std_MSMS_dir <- file.path(results.dir,"p2-analysis-setting", "MQ_Std_MSMS")

# Run stoichiometry processing for MSMS-only
process_stoichiometry(
  MQ_standard_data = MQ_standard_data,
  output_dir       = MQ_Std_MSMS_dir,
  ref_protein_data = ref_protein_data,
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

    ## Processing: data-C_trp_m8_v8_mCC_evidence.txt

    ## Processing: data-D_trp_m2_v2_def_evidence.txt

    ## Processing: data-D_trp_m5_v5_def_evidence.txt

    ## Processing: data-D_trp_m7_v7_def_evidence.txt

# 2.1.4 Calculation of Stoichiometry of MQ_DI without waterloss

``` r
# Load data into environment
MQ_DI_noH2O_data <- file.path(
  project.dir,
  "data/MQ_with_DI_no_waterloss/PNAS2022" 
)

# Create file path for results
MQ_DI_noH2O_dir <- file.path(results.dir, "p2-analysis-setting", "MQ_DI_noH2O")
create.dir(MQ_DI_noH2O_dir)
```

``` r
# Load sample run info data 
MQ_DI_noH2O_run_info <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_with_DI_noH2O" 
) %>% data.table

# Add a new column 'sample_id' that assigns a unique sequential number to each 
#row (1:N), where N is the total number of rows in the data.table
MQ_DI_noH2O_run_info[, sample_id := 1:.N]
```

``` r
# For each sample ID calculate the stoichiomtery
for (i in seq_len(nrow(MQ_DI_noH2O_run_info))){ # Takes each row from the MQ_DI_noH2O_run_info
  
  prefix <- MQ_DI_noH2O_run_info$prefix[i]
  
  # Retrieves the evidence.txt file for each sample ID 
  mq_evidence_data = file.path(MQ_DI_noH2O_data,
                               "evidence",
                               paste0(prefix, "evidence.txt"))
  
  print(paste0("Processing: ", basename(mq_evidence_data)))
  
  # Checks whether the sample ID has the corresponding evidence.txt file. If yes, then proceed
  if (file.exists(mq_evidence_data)) { 
    
    # Locate and fetch all PTM site files within the ptm folder
    ptm_files <- list.files(file.path(
      MQ_DI_noH2O_data, "ptm",
      paste0(MQ_DI_noH2O_run_info[i, prefix], "ptm")), 
      full.names = TRUE)
    
    # Generate PTM names based on file names, add brackets and remove "Sites.txt"
    ptm_names <- paste0(
      "[", str_replace_all(basename(ptm_files), "Sites.txt", ""), "]"
    )
    
    # Run the stoichiometry calculation 
    stoic.dt <- calculate_stoichiometry2(
      mq_evidence_data = mq_evidence_data,
      sample_info_file = file.path( # Corresponding sample info data 
        MQ_DI_noH2O_data,
        "sample_info",
        paste0(
          "MS_dataset_overview_PXD031221_",
          MQ_DI_noH2O_run_info[i, data],
          ".csv"
        )
      ),
      ref_protein_dt = ref_protein_data,
      ptm_files = ptm_files,
      ptm_names = ptm_names,
      output_prefix = file.path(MQ_DI_noH2O_dir,  # Define output prefix for results
                                paste0("DI_", "noH2O_", MQ_DI_noH2O_run_info[i, prefix])),
      ptm_mapping_file = ptm_mapping_file,
      K_only = FALSE,
      selected_type = c("MULTI-MSMS"),
    )
    # If evidence.txt does not exist, skip the function and print below message
  } else {
    print(paste0("File does not exist: ", mq_evidence_data))
  }
  
  
  gc()
  
}
```

    ## [1] "Processing: data-A_trp_m7_v7_def_evidence.txt"
    ## [1] "Processing: data-B1_trp_m7_v7_mCC_evidence.txt"
    ## [1] "Processing: data-B2_trp_m7_v7_mCC_evidence.txt"
    ## [1] "Processing: data-C_trp_m7_v7_mCC_evidence.txt"

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
    ##  date     2026-01-07
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
    ##  evaluate            1.0.4      2025-06-18 [1] CRAN (R 4.4.3)
    ##  farver              2.1.2      2024-05-13 [1] CRAN (R 4.4.3)
    ##  fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.4.3)
    ##  generics            0.1.4      2025-05-09 [1] CRAN (R 4.4.3)
    ##  GenomeInfoDb      * 1.42.3     2025-01-27 [1] Bioconduc~
    ##  GenomeInfoDbData    1.2.13     2025-07-21 [1] Bioconductor
    ##  ggplot2           * 3.5.2      2025-04-09 [1] CRAN (R 4.4.3)
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
    ##  magrittr          * 2.0.3      2022-03-30 [1] CRAN (R 4.4.3)
    ##  pillar              1.11.0     2025-07-04 [1] CRAN (R 4.4.3)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.3)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-12-13 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.4.3)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.4.3)
    ##  readxl            * 1.4.5      2025-03-07 [1] CRAN (R 4.4.3)
    ##  rlang               1.1.6      2025-04-11 [1] CRAN (R 4.4.3)
    ##  rmarkdown           2.29       2024-11-04 [1] CRAN (R 4.4.3)
    ##  rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.4.3)
    ##  S4Vectors         * 0.44.0     2024-10-29 [1] Bioconduc~
    ##  scales              1.4.0      2025-04-24 [1] CRAN (R 4.4.3)
    ##  sessioninfo         1.2.3      2025-02-05 [1] CRAN (R 4.4.3)
    ##  snakecase           0.11.1     2023-08-27 [1] CRAN (R 4.4.3)
    ##  stringi             1.8.7      2025-03-27 [1] CRAN (R 4.4.3)
    ##  stringr           * 1.5.1      2023-11-14 [1] CRAN (R 4.4.3)
    ##  tibble              3.3.0      2025-06-08 [1] CRAN (R 4.4.3)
    ##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.4.3)
    ##  timechange          0.3.0      2024-01-18 [1] CRAN (R 4.4.3)
    ##  UCSC.utils          1.2.0      2024-10-29 [1] Bioconduc~
    ##  vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.4.3)
    ##  withr               3.0.2      2024-10-28 [1] CRAN (R 4.4.3)
    ##  xfun                0.52       2025-04-02 [1] CRAN (R 4.4.3)
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
