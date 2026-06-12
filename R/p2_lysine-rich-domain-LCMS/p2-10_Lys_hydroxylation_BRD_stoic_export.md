2-10. Lysine hydroxylations in varying O2 pc and dox induction time -
excel data export
================
Yoichiro Sugimoto and Pallavi Kesavan
12 June, 2026

- [Overview](#overview)
- [Environment setup](#environment-setup)
- [2.10.1 Install,load essential functions and
  libraries](#2101-installload-essential-functions-and-libraries)
- [2.10.2 Import human protein reference
  data](#2102-import-human-protein-reference-data)
- [2.10.3 MS_SS, MS_KR, PNAS Data tables converted to
  excel](#2103-ms_ss-ms_kr-pnas-data-tables-converted-to-excel)
- [Session information](#session-information)

# Overview

This script exports processed lysine hydroxylation stoichiometry data
for BRD protein sites to a formatted Excel sheet.

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

# 2.10.1 Install,load essential functions and libraries

``` r
## Load all R scripts from the 'functions' folder into the current session
P2_functions <- sapply(list.files
                       (file.path(project.dir, "R/functions"), 
                         pattern="*.R", 
                         full.names = TRUE), 
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
## Install private package 
# Install ptm.stiochiometry package
#install.packages("/fast/AG_Sugimoto/home/users/pallavi/projects/ptm.stoichiometry", repos = NULL, type = "source")

# install package for writing data table to excel 
# install.packages('openxlsx')

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
library("openxlsx")
```

# 2.10.2 Import human protein reference data

``` r
# Import human protein reference data from specified file path 
ref_protein_dt <- import_reference_fasta(file.path
                                         ("/fast/AG_Sugimoto/reference/uniprot/human", 
                                           "UP000005640_9606.fasta")) 
```

``` r
# Define path to data directory
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

#Create p1 results directory
p2_MS_SS_KR <- file.path(results.dir, "p2_MS_SS_KR")
# create.dirs(c(results.dir, p2_MS_SS_KR))
```

# 2.10.3 MS_SS, MS_KR, PNAS Data tables converted to excel

``` r
# Read FASTA data 
all.protein.bs <- Biostrings::readAAStringSet(
  file.path(
    "/fast/AG_Sugimoto/reference/uniprot/human",
    "UP000005640_9606.fasta"
  )
)

# Define function - 'read_stoic_data'
read_stoic_data <- function(prefix, pre_prefix, post_fix = "", dir_path){
  # Read csv file with defined file path,  
  dt <- fread(file.path(dir_path, paste0(pre_prefix, prefix, "PTM_stoichiometry", post_fix, ".csv"))) 
  dt[, condition := gsub("_$", "", prefix)] # Removes '_' at the end of the prefix and creates a new column called condition
  return(dt)
}

#---------------------------------------------------------------------
# Load PNAS stoichiometry data analysed by MaxQuant with DI (dataset A)
#---------------------------------------------------------------------

 MQ_DI_dtA_stoic_dt <- read_stoic_data(
  prefix = "DI_noH2O_data-A_trp_m7_v7_def_",
  pre_prefix = "",
  post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MQ_DI_noH2O"))

# Create a logical column for diagnostic peak. For rows where diagnostic peak is 
# present "+", logical vector is "TRUE", or else "FALSE"
MQ_DI_dtA_stoic_dt[, `:=`(
  is_diagnostic_peak = diagnostic_peak == "+" 
)]

#--------------------------------------------------------------------------
# Load PNAS stoichiometry data analysed by MaxQuant with DI (dataset B1 and B2)
#--------------------------------------------------------------------------

# Load sample run info data (MQ DI with no waterloss)
MQ_DI_noH2O_run_info <- read_excel(
  file.path(project.dir, "data/analysis_setting/PXD031221_sample_matrix.xlsx"),
  sheet = "MQ_with_DI_noH2O" 
) %>% data.table

# Add a new column 'sample_id' that assigns a unique sequential number to each 
#row (1:N), where N is the total number of rows in the data.table
MQ_DI_noH2O_run_info[, sample_id := 1:.N]

MQ_DI_dtB_stoic_dt <- lapply(
  MQ_DI_noH2O_run_info[data != c("data-A", "data-C"), prefix],
  read_stoic_data,
  pre_prefix = "DI_noH2O_", 
  post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MQ_DI_noH2O")
) %>% rbindlist

# Create a logical column for diagnostic peak. For rows where diagnostic peak is 
# present "+", logical vector is "TRUE", or else "FALSE"
MQ_DI_dtB_stoic_dt[, `:=`(
  is_diagnostic_peak = diagnostic_peak == "+" 
)]

#---------------------------------------------------------------------
# Load PNAS stoichiometry data analysed by MaxQuant with DI (dataset C)
#---------------------------------------------------------------------

 MQ_DI_dtC_stoic_dt <- read_stoic_data(
  prefix = "DI_noH2O_data-C_trp_m7_v7_mCC_",
  pre_prefix = "",
  post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MQ_DI_noH2O"))

# Create a logical column for diagnostic peak. For rows where diagnostic peak is 
# present "+", logical vector is "TRUE", or else "FALSE"
MQ_DI_dtC_stoic_dt[, `:=`(
  is_diagnostic_peak = diagnostic_peak == "+" 
)]


#------------------------------------------------------
# Load PNAS stoichiometry data manually curated sites
#------------------------------------------------------

# Load PNAS stoichiometry data
pnas2022_stoic_dt <- fread(file.path(data.dir, "PNAS2022/long_K_stoichiometry_data.csv"))

# Replace old column names with new ones
setnames(
  pnas2022_stoic_dt, 
  old = c("uniprot_id", "position", "residue", "oxK_ratio", "data_source", "total_n_feature_oxK", "total_n_feature_K"), 
  new = c("protein_accession", "aa_pos", "aa", "stoichiometry", "sample_name", "sum_psm_mapped", "sum_psm_mapped_per_position")
)

# Merge the reference protein data with the PNAS stoichiometry data by protein accession ID
pnas2022_stoic_dt <- merge(
  ref_protein_dt[, .(protein_accession, gene_name)], # Take columns protein_accession and gene_name from the ref data 
  pnas2022_stoic_dt,
  by = "protein_accession" #combine both data.table by protein accession ID
)

## Add columns to pnas2022.stoic.dt metadata
# new columns: accession_position, ptm, sample_name 
pnas2022_stoic_dt[, `:=`(
  accession_position = paste0(protein_accession, "_", aa_pos), #this column combines data from protein_accession and aa_pos
    ptm = fifelse(curated_oxK_site == TRUE, "[Oxidation (K)]", " ")
    )
  ]

#Change col names to make it same across datasets
setnames(pnas2022_stoic_dt, c("stoichiometry", "protein_accession", "gene_name", "aa_pos", "curated_oxK_site"),  c("Stoichiometry_PNAS2022", "Protein_accession", "Gene_name", "Amino_acid_pos", "PNAS2022_curated_oxK_site"))

#-----------------------------
# Load MS_KR_1 stoichiometry data (Dataset D)
#-----------------------------

## Read stoichiometry data for MS_KR1 data (noH2O loss)
 MS_KR1_stoic_dt <- read_stoic_data(
  prefix = "MS_KR_1_noH2O_",
  pre_prefix = "",
  post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MS_KR_1_noH2O"))


# Create a logical column for diagnostic peak. For rows where diagnostic peak is 
# present "+", logical vector is "TRUE", or else "FALSE"
MS_KR1_stoic_dt[, `:=`(
  is_diagnostic_peak = diagnostic_peak == "+" 
)]

#-----------------------------
# Load MS_SS stoichiometry data (Dataset E)
#-----------------------------

## Read stoichiometry data for MS_SS data 
 MS_SS_stoic_dt <- read_stoic_data(
  prefix = "MS_SS_",
  pre_prefix = "",
  post_fix = "_DI",
  dir_path = file.path(results.dir, "p2-analysis-setting", "MS_SS_noH2O"))

 #[Oxidation (K)] 

# Create a logical column for diagnostic peak. For rows where diagnostic peak is 
# present "+", logical vector is "TRUE", or else "FALSE"
MS_SS_stoic_dt[, `:=`(
  is_diagnostic_peak = diagnostic_peak == "+" 
)]
```

``` r
# Subset the data to contain only Lysine
MQ_DI_dtA_stoic_dt <- MQ_DI_dtA_stoic_dt[aa == "K"]

# removing unused data columns from data table
MQ_DI_dtA_stoic_dt <-  MQ_DI_dtA_stoic_dt[, c("sum_intensity_per_position",
                                          "sum_psm_mapped_per_position",
                                          "sum_the_number_of_peptide",
                                          "localization_prob",
                                          "score_for_localization",
                                          "best_localization_ms_ms_id",
                                          "best_localization_raw_file",
                                          "max_score",
                                          "diagnostic_peak",
                                          "condition"
                                          
) := NULL]

# changing stoichiometry column position 
setcolorder(
  MQ_DI_dtA_stoic_dt,
  append(setdiff(names(MQ_DI_dtA_stoic_dt), "stoichiometry"), "stoichiometry", after = 6)
)
 
# changing sample_group column position 
setcolorder(
  MQ_DI_dtA_stoic_dt,
  append(setdiff(names(MQ_DI_dtA_stoic_dt), "is_diagnostic_peak"), "is_diagnostic_peak", after = 13)
)

 # changing column names - publication ready
setnames(MQ_DI_dtA_stoic_dt, c("sample_name",
                              "protein_accession",
                              "gene_name",
                              "aa",
                              "aa_pos",
                              "ptm", 
                              "stoichiometry",
                              "sum_peak_intensity",
                              "sum_psm_mapped",
                              "the_number_of_peptide",
                              "is_diagnostic_peak"), 
         c("Sample_Name",
           "Protein_accession",
           "Gene_name",
           "Amino_acid",
           "Amino_acid_pos",
           "PTM_type",
           "Stoichiometry",
           "Total_peak_intensity",
           "Total_PSM_mapped",
           "Number_of_peptide",
           "Diagnostic_peak_presence"))

# Merge MQ_DI stoic dataset A and pnas2022 manually curated sites by 'Protein_accession' and 'Amino_acid_pos'
 MQ_DI_dtA_stoic_dt <- MQ_DI_dtA_stoic_dt[, `:=`( 
  PNAS2022_curated_oxK_site = 
    paste0(Protein_accession, "_", Amino_acid_pos) %in%
    pnas2022_stoic_dt[PNAS2022_curated_oxK_site == TRUE, paste0(Protein_accession, "_", Amino_acid_pos)] 
)]

 # Write to csv table to identify "K" with DI 
write.xlsx(MQ_DI_dtA_stoic_dt, file.path(results.dir, "p2_MS_SS_KR", "MQ_DI_dtA_stoic_dt.xlsx"))
```

``` r
# Subset the data to contain only Lysine
MQ_DI_dtB_stoic_dt <- MQ_DI_dtB_stoic_dt[aa == "K"]

# removing unused data columns from data table
MQ_DI_dtB_stoic_dt <-  MQ_DI_dtB_stoic_dt[, c("sum_intensity_per_position",
                                          "sum_psm_mapped_per_position",
                                          "sum_the_number_of_peptide",
                                          "localization_prob",
                                          "score_for_localization",
                                          "best_localization_ms_ms_id",
                                          "best_localization_raw_file",
                                          "max_score",
                                          "diagnostic_peak",
                                          "condition"
                                          
) := NULL]

# changing stoichiometry column position 
setcolorder(
  MQ_DI_dtB_stoic_dt,
  append(setdiff(names(MQ_DI_dtB_stoic_dt), "stoichiometry"), "stoichiometry", after = 6)
)
 
# changing sample_group column position 
setcolorder(
  MQ_DI_dtB_stoic_dt,
  append(setdiff(names(MQ_DI_dtB_stoic_dt), "is_diagnostic_peak"), "is_diagnostic_peak", after = 13)
)

 # changing column names - publication ready
setnames(MQ_DI_dtB_stoic_dt, c("sample_name",
                              "protein_accession",
                              "gene_name",
                              "aa",
                              "aa_pos",
                              "ptm", 
                              "stoichiometry",
                              "sum_peak_intensity",
                              "sum_psm_mapped",
                              "the_number_of_peptide",
                              "is_diagnostic_peak"), 
         c("Sample_name",
           "Protein_accession",
           "Gene_name",
           "Amino_acid",
           "Amino_acid_pos",
           "PTM_type",
           "Stoichiometry",
           "Total_peak_intensity",
           "Total_PSM_mapped",
           "Total_number_of_peptide",
           "Diagnostic_peak_presence"))

# Merge MQ_DI stoic dataset B1 and B2 and pnas2022 manually curated sites by 'Protein_accession' and 'Amino_acid_pos'
 MQ_DI_dtB_stoic_dt <- MQ_DI_dtB_stoic_dt[, `:=`( 
  PNAS2022_curated_oxK_site = 
    paste0(Protein_accession, "_", Amino_acid_pos) %in%
    pnas2022_stoic_dt[PNAS2022_curated_oxK_site == TRUE, paste0(Protein_accession, "_", Amino_acid_pos)] 
)]

 # Write to csv table to identify "K" with DI 
write.xlsx(MQ_DI_dtB_stoic_dt, file.path(results.dir, "p2_MS_SS_KR", "MQ_DI_dtB_stoic_dt.xlsx"))
```

``` r
# Subset the data to contain only Lysine
MQ_DI_dtC_stoic_dt <- MQ_DI_dtC_stoic_dt[aa == "K"]

# removing unused data columns from data table
MQ_DI_dtC_stoic_dt <-  MQ_DI_dtC_stoic_dt[, c("sum_intensity_per_position",
                                          "sum_psm_mapped_per_position",
                                          "sum_the_number_of_peptide",
                                          "localization_prob",
                                          "score_for_localization",
                                          "best_localization_ms_ms_id",
                                          "best_localization_raw_file",
                                          "max_score",
                                          "diagnostic_peak",
                                          "condition"
                                          
) := NULL]

# changing stoichiometry column position 
setcolorder(
  MQ_DI_dtC_stoic_dt,
  append(setdiff(names(MQ_DI_dtC_stoic_dt), "stoichiometry"), "stoichiometry", after = 6)
)
 
# changing sample_group column position 
setcolorder(
  MQ_DI_dtC_stoic_dt,
  append(setdiff(names(MQ_DI_dtC_stoic_dt), "is_diagnostic_peak"), "is_diagnostic_peak", after = 13)
)

# changing column names - publication ready
setnames(MQ_DI_dtC_stoic_dt, c("sample_name",
                              "protein_accession",
                              "gene_name",
                              "aa",
                              "aa_pos",
                              "ptm", 
                              "stoichiometry",
                              "sum_peak_intensity",
                              "sum_psm_mapped",
                              "the_number_of_peptide",
                              "is_diagnostic_peak"), 
         c("Sample_Name",
           "Protein_accession",
           "Gene_name",
           "Amino_acid",
           "Amino_acid_pos",
           "PTM_type",
           "Stoichiometry",
           "Total_peak_intensity",
           "Total_PSM_mapped",
           "Number_of_peptide",
           "Diagnostic_peak_presence"))

# Merge MQ_DI stoic dataset C and pnas2022 manually curated sites by 'Protein_accession' and 'Amino_acid_pos'
 MQ_DI_dtC_stoic_dt <- MQ_DI_dtC_stoic_dt[, `:=`( 
  PNAS2022_curated_oxK_site = 
    paste0(Protein_accession, "_", Amino_acid_pos) %in%
    pnas2022_stoic_dt[PNAS2022_curated_oxK_site == TRUE, paste0(Protein_accession, "_", Amino_acid_pos)] 
)]
 
 # Write to csv table to identify "K" with DI 
write.xlsx(MQ_DI_dtC_stoic_dt, file.path(results.dir, "p2_MS_SS_KR", "MQ_DI_dtC_stoic_dt.xlsx"))
```

``` r
# Subset the data to contain only Lysine
MS_KR1_stoic_K_dt <- MS_KR1_stoic_dt[aa == "K"]

# removing unused data columns from data table
MS_KR1_stoic_K_dt <-  MS_KR1_stoic_K_dt[, c("sum_intensity_per_position",
                                          "sum_psm_mapped_per_position",
                                          "sum_the_number_of_peptide",
                                          "localization_prob",
                                          "score_for_localization",
                                          "best_localization_ms_ms_id",
                                          "best_localization_raw_file",
                                          "max_score",
                                          "diagnostic_peak",
                                          "condition") := NULL]

# changing stoichiometry column position 
setcolorder(
  MS_KR1_stoic_K_dt,
  append(setdiff(names(MS_KR1_stoic_K_dt), "stoichiometry"), "stoichiometry", after = 6)
)
 
# changing sample_group column position 
setcolorder(
  MS_KR1_stoic_K_dt,
  append(setdiff(names(MS_KR1_stoic_K_dt), "is_diagnostic_peak"), "is_diagnostic_peak", after = 13)
)

 # changing column names - publication ready
setnames(MS_KR1_stoic_K_dt, c("sample_name",
                              "protein_accession",
                              "gene_name",
                              "aa",
                              "aa_pos",
                              "ptm", 
                              "stoichiometry",
                              "sum_peak_intensity",
                              "sum_psm_mapped",
                              "the_number_of_peptide",
                              "is_diagnostic_peak"), 
         c("Sample_Name",
           "Protein_accession",
           "Gene_name",
           "Amino_acid",
           "Amino_acid_pos",
           "PTM_type",
           "Stoichiometry",
           "Total_peak_intensity",
           "Total_PSM_mapped",
           "Number_of_peptide",
           "Diagnostic_peak_presence"))


# Merge MS_KR_1 and pnas2022 manually curated sites by 'Protein_accession' and 'Amino_acid_pos'
 MS_KR1_stoic_K_dt <- MS_KR1_stoic_K_dt[, `:=`( 
  PNAS2022_curated_oxK_site = 
    paste0(Protein_accession, "_", Amino_acid_pos) %in%
    pnas2022_stoic_dt[PNAS2022_curated_oxK_site == TRUE, paste0(Protein_accession, "_", Amino_acid_pos)] 
)]

 # Write to Excel table
write.xlsx(MS_KR1_stoic_K_dt, file.path(results.dir, "p2_MS_SS_KR", "MS_KR1_stoic_K.xlsx"))
```

``` r
# Subset the data to contain only Lysine
MS_SS_stoic_K_dt <- MS_SS_stoic_dt[aa == "K"]

# removing unsed data columns from data table
MS_SS_stoic_K_dt <-  MS_SS_stoic_K_dt[, c("sum_intensity_per_position",
                                          "sum_psm_mapped_per_position",
                                          "sum_the_number_of_peptide",
                                          "localization_prob",
                                          "score_for_localization",
                                          "best_localization_ms_ms_id",
                                          "best_localization_raw_file",
                                          "max_score",
                                          "diagnostic_peak",
                                          "condition") := NULL]

# Remove BRD2/3 samples reporting BRD4 stoichiomtery and vice versa
# These samples are removed as it could be false positives reported by the system 
MS_SS_stoic_K_dt <- MS_SS_stoic_K_dt[
  !((grepl("_BRD23$", sample_name) & gene_name == "BRD4") |
       (grepl("_BRD4$", sample_name) & gene_name %in% c("BRD2", "BRD3")))
]

# changing stoichiometry column position 
setcolorder(
  MS_SS_stoic_K_dt,
  append(setdiff(names(MS_SS_stoic_K_dt), "stoichiometry"), "stoichiometry", after = 6)
)

# changing sample_group column position 
setcolorder(
  MS_SS_stoic_K_dt,
  append(setdiff(names(MS_SS_stoic_K_dt), "is_diagnostic_peak"), "is_diagnostic_peak", after = 13)
)

 # changing column names - publication ready
setnames(MS_SS_stoic_K_dt, c("sample_name",
                              "protein_accession",
                              "gene_name",
                              "aa",
                              "aa_pos",
                              "ptm", 
                              "stoichiometry",
                              "sum_peak_intensity",
                              "sum_psm_mapped",
                              "the_number_of_peptide",
                              "is_diagnostic_peak"), 
         c("Sample_Name",
           "Protein_accession",
           "Gene_name",
           "Amino_acid",
           "Amino_acid_pos",
           "PTM_type",
           "Stoichiometry",
           "Total_peak_intensity",
           "Total_PSM_mapped",
           "Number_of_peptide",
           "Diagnostic_peak_presence"))


# Merge MS_SS and pnas2022 manually curated sites by 'Protein_accession' and 'Amino_acid_pos'
 MS_SS_stoic_K_dt <- MS_SS_stoic_K_dt[, `:=`( 
  PNAS2022_curated_oxK_site = 
    paste0(Protein_accession, "_", Amino_acid_pos) %in%
    pnas2022_stoic_dt[PNAS2022_curated_oxK_site == TRUE, paste0(Protein_accession, "_", Amino_acid_pos)] 
)]

 # Write to Excel table
write.xlsx(MS_SS_stoic_K_dt, file.path(results.dir, "p2_MS_SS_KR", "MS_SS_stoic_K.xlsx"))
```

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
    ##  date     2026-06-12
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
    ##  openxlsx          * 4.2.8.1    2025-10-31 [1] CRAN (R 4.5.1)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.5.1)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.5.1)
    ##  ptm.stoichiometry * 0.0.0.9000 2026-05-15 [1] local
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.5.1)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.5.1)
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
    ##  zip                 2.3.3      2025-05-13 [1] CRAN (R 4.5.1)
    ## 
    ##  [1] /home/pkesava/R/x86_64-pc-linux-gnu-library/4.5
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
