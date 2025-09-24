2-1. Calculate PTM stoichiometry
================
Yoichiro Sugimoto and Pallavi Kesavan
24 September, 2025

- [Overview](#overview)
- [Environment setup](#environment-setup)
- [Install and load essentials
  functions](#install-and-load-essentials-functions)
- [Import data](#import-data)
- [Calculation of Stoichiometry without diagnostic ion
  consideration](#calculation-of-stoichiometry-without-diagnostic-ion-consideration)
- [Calculation of Stoichiometry with diagnostic ion
  consideration](#calculation-of-stoichiometry-with-diagnostic-ion-consideration)
- [Calculation of Stoichiometry with diagnostic ion and
  iterative](#calculation-of-stoichiometry-with-diagnostic-ion-and-iterative)
- [Calculation of Stoichiometry with diagnostic ion, iterative and
  additional
  PTMs](#calculation-of-stoichiometry-with-diagnostic-ion-iterative-and-additional-ptms)
- [Session information](#session-information)

# Overview

# Environment setup

``` r
## Initialize renv (first time only) - re-installed 24.09.2025
# Creates project specific library and renv.lock file. 
# Use 'renv::init(filepath)' to create project library
# renv::init(
#        "/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains/R"
#    )


project.dir <-
  file.path(
    "/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains"
  )

# renv::snapshot(file.path(project.dir, "R"))

#renv::restore completed 24.09.2025
#renv::restore(file.path(project.dir, "R"))
```

# Install and load essentials functions

``` r
## Load all R scripts from the 'functions' folder into the current session
temp <- sapply(list.files(file.path(project.dir, "R/functions"), pattern="*.R", full.names = TRUE), source)
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
## Install private packages 
# Install ptm.stiochiometry package
install.packages("/fast/AG_Sugimoto/home/users/pallavi/projects/ptm.stoichiometry", repos = NULL, type = "source")
```

    ## Installing package into '/home/pkesava/R/x86_64-pc-linux-gnu-library/4.5'
    ## (as 'lib' is unspecified)

``` r
# Load Libraries - ptm.stiochiometry
library(ptm.stoichiometry)
```

# Import data

``` r
# Import human protien reference data from specified file path 
ref_protein_data <- import_reference_fasta(file.path
("/fast/AG_Sugimoto/reference/uniprot/human", 
  "UP000005640_9606.fasta")) 
```

# Calculation of Stoichiometry without diagnostic ion consideration

# Calculation of Stoichiometry with diagnostic ion consideration

# Calculation of Stoichiometry with diagnostic ion and iterative

# Calculation of Stoichiometry with diagnostic ion, iterative and additional PTMs

# Session information

``` r
## sessioninfo::session_info()
```
