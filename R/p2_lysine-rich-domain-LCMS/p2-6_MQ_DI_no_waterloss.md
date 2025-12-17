2-6. Stoichiometry calculations without considering water loss
================
Yoichiro Sugimoto and Pallavi Kesavan
16 December, 2025

- [This script calculates the stoichiometry for MQ data with diagnostic
  ions and no water
  loss](#this-script-calculates-the-stoichiometry-for-mq-data-with-diagnostic-ions-and-no-water-loss)
- [Environment setup](#environment-setup)
- [2.6.1 Install,load essential functions and
  libraries](#261-installload-essential-functions-and-libraries)
- [2.6.2 Import human protein reference
  data](#262-import-human-protein-reference-data)
- [2.6.3 Calculation of Stoichiometry of MQ_DI without
  waterloss](#263-calculation-of-stoichiometry-of-mq_di-without-waterloss)
- [2.5.3 Calculation of Stoichiometry with diagnostic ion in hypoxia and
  normoxia](#253-calculation-of-stoichiometry-with-diagnostic-ion-in-hypoxia-and-normoxia)
- [Session information](#session-information)

# This script calculates the stoichiometry for MQ data with diagnostic ions and no water loss

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

# 2.6.1 Install,load essential functions and libraries

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
## Install private package 
# Install ptm.stiochiometry package
#install.packages("/fast/AG_Sugimoto/home/users/pallavi/projects/ptm.stoichiometry", repos = NULL, type = "source")

# Load Libraries
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

# 2.6.2 Import human protein reference data

``` r
# Import human protein reference data from specified file path 
ref_protein_dt <- import_reference_fasta(file.path
                                         ("/fast/AG_Sugimoto/reference/uniprot/human", 
                                           "UP000005640_9606.fasta")) 
```

# 2.6.3 Calculation of Stoichiometry of MQ_DI without waterloss

``` r
# Define path to data directory
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

# Define file path common PTM mapping file
ptm_mapping_file <- file.path(
  project.dir,
  "data/analysis_setting/ptm_replacement_de-propionylate_for_hydroxylysine.csv"
)

# Load data into environment
MQ_DI_noH2O_data <- file.path(
  project.dir,
  "data/MQ_with_DI_no_waterloss/PNAS2022" 
)

# Create file path for results
MQ_DI_noH2O_dir <- file.path(project.dir, "results", "p2-analysis-setting", "MQ_DI_noH2O")
# dir.create(MQ_DI_noH2O_dir, recursive = TRUE)
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
      ref_protein_dt = ref_protein_dt,
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

# 2.5.3 Calculation of Stoichiometry with diagnostic ion in hypoxia and normoxia

``` r
# Load data into environment
MS_KR1_no_H2O_data <- file.path(
  project.dir,
  "data/MQ_with_DI_no_waterloss/MS_KR_1" 
)

# Create file path for results
MS_KR_1_noH2O_dir <- file.path(project.dir, "results", "p2-analysis-setting", "MS_KR_1_noH2O")
#dir.create(MS_KR_1_noH2O_dir, recursive = TRUE)
```

``` r
# Define file path to MaxQuant evidence files
mq_evidence_data <- file.path(MS_KR1_no_H2O_data, 
                              "MS_KR_1_evidence.txt")

# Print which sample id is being processed
message("Processing: ", basename(mq_evidence_data)) 
```

    ## Processing: MS_KR_1_evidence.txt

``` r
# Checks whether the sample ID has the corresponding evidence.txt file. If yes, then proceed
if (file.exists(mq_evidence_data)) {
  
  # Locate and fetch all PTM site files within the ptm folder
  ptm_files <- list.files(file.path(
    MS_KR1_no_H2O_data, "ptm"), 
    full.names = TRUE)
  
  # Generate PTM names based on file names, add brackets and remove "Sites.txt"
  ptm_names <- paste0(
    "[", str_replace_all(basename(ptm_files), "Sites.txt", ""), "]"
  )
  
  # Run the stoichiometry calculation
  stoic.dt <- calculate_stoichiometry2(
    mq_evidence_data = mq_evidence_data,
    sample_info_file = file.path(
      MS_KR1_no_H2O_data,
      "sample_info.csv"
    ),
    ref_protein_dt = ref_protein_dt,
    ptm_files = ptm_files,
    ptm_names = ptm_names,
    output_prefix = file.path(MS_KR_1_noH2O_dir, paste0("MS_KR_1_", "noH2O_")),
    ptm_mapping_file = ptm_mapping_file,
    K_only = FALSE,
    selected_type = "MULTI-MSMS",
    parse_protein_accession_function = NA
  )
} else { # If evidence.txt does not exist, skip the function and print below message
  print(paste0("File does not exist: ", mq_evidence_data))
}

# Garbage collection - to free up memory
gc()
```

    ##            used  (Mb) gc trigger  (Mb)  max used  (Mb)
    ## Ncells  4246087 226.8    8206968 438.3   8206968 438.3
    ## Vcells 16971015 129.5   59388007 453.1 115666095 882.5

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
    ##  date     2025-12-16
    ##  pandoc   3.4 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.6.42 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  BiocGenerics      * 0.54.0     2025-04-15 [1] Bioconduc~
    ##  Biostrings        * 2.76.0     2025-04-15 [1] Bioconduc~
    ##  bit                 4.6.0      2025-03-06 [1] CRAN (R 4.5.1)
    ##  bit64               4.6.0-1    2025-01-16 [1] CRAN (R 4.5.1)
    ##  cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.5.1)
    ##  cli                 3.6.5      2025-04-23 [1] CRAN (R 4.5.1)
    ##  crayon              1.5.3      2024-06-20 [1] CRAN (R 4.5.1)
    ##  data.table        * 1.17.8     2025-07-10 [1] CRAN (R 4.5.1)
    ##  digest              0.6.37     2024-08-19 [1] CRAN (R 4.5.1)
    ##  dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.5.1)
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
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.5.1)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.5.1)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.5.1)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.5.1)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.5.1)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-12-16 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.5.1)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.5.1)
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
