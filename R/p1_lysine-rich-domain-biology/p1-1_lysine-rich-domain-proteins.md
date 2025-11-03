1-1. Proteins with lysine-rich domains
================
Pallavi Kesavan and Yoichiro Sugimoto
03 November, 2025

- [Environment setup](#environment-setup)
- [Set up environment](#set-up-environment)
- [Import data](#import-data)
- [Maximum K score and Protein
  Length](#maximum-k-score-and-protein-length)
- [Data Visualization - K score, Maximum K score and Protein
  Length](#data-visualization---k-score-maximum-k-score-and-protein-length)
- [RNA-Binding proteins](#rna-binding-proteins)
- [Data Visualization - RNA-Binding Proteins (Max K score versus
  RBP)](#data-visualization---rna-binding-proteins-max-k-score-versus-rbp)
- [Subcellular localisation](#subcellular-localisation)
- [Data Visualization - Subcellular
  Localisation](#data-visualization---subcellular-localisation)
- [Session information](#session-information)

# Environment setup

``` r
## Initialize renv (first time only) - re-installed 24.09.2025
# Creates project specific library and renv.lock file. 
# Use 'renv::init(filepath)' to create project library
# renv::init(
#        "/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains/R"
#    )

## Define project directory 
project.dir <- file.path("/fast/AG_Sugimoto/home/users/pallavi/projects/20241111_PTMs_in_lysine_rich_domains")

## Use 'renv::snapshot', if there are updates in the current dependencies 
# renv::snapshot(file.path(project.dir, "R"))

## Restore packages
renv::restore(file.path(project.dir, "R"))
```

    ## The following required system packages are not installed:
    ## - pandoc  [required by knitr, rmarkdown]
    ## The R packages depending on these system packages may fail to install.
    ## 
    ## An administrator can install these packages with:
    ## - sudo apt install pandoc
    ## 
    ## - The library is already synchronized with the lockfile.

``` r
## Check package versions of current loaded library and lockfile 
renv::status
```

    ## function (project = NULL, ..., library = NULL, lockfile = NULL, 
    ##     sources = TRUE, cache = FALSE, dev = FALSE) 
    ## {
    ##     renv_scope_error_handler()
    ##     renv_dots_check(...)
    ##     renv_snapshot_auto_suppress_next()
    ##     renv_scope_options(renv.prompt.enabled = FALSE)
    ##     the$status_running <- TRUE
    ##     defer(the$status_running <- FALSE)
    ##     project <- renv_project_resolve(project)
    ##     renv_project_lock(project = project)
    ##     if (!renv_status_check_initialized(project, library, lockfile)) {
    ##         result <- list(library = list(Packages = named(list())), 
    ##             lockfile = list(Packages = named(list())), synchronized = FALSE)
    ##         return(invisible(result))
    ##     }
    ##     libpaths <- library %||% renv_libpaths_resolve()
    ##     lockpath <- lockfile %||% renv_paths_lockfile(project = project)
    ##     dependencies <- renv_snapshot_dependencies(project, dev = dev)
    ##     packages <- sort(union(dependencies, "renv"))
    ##     paths <- renv_package_dependencies(packages, libpaths = libpaths, 
    ##         project = project)
    ##     packages <- as.character(names(paths))
    ##     lockfile <- if (file.exists(lockpath)) 
    ##         renv_lockfile_read(lockpath)
    ##     else renv_lockfile_init(project = project)
    ##     library <- renv_lockfile_create(libpaths = libpaths, type = "all", 
    ##         prompt = FALSE, project = project)
    ##     ignored <- c(renv_project_ignored_packages(project), renv_packages_base(), 
    ##         if (renv_tests_running()) "renv")
    ##     packages <- setdiff(packages, ignored)
    ##     renv_lockfile_records(lockfile) <- omit(renv_lockfile_records(lockfile), 
    ##         ignored)
    ##     renv_lockfile_records(library) <- omit(renv_lockfile_records(library), 
    ##         ignored)
    ##     synchronized <- all(renv_status_check_consistent(lockfile, 
    ##         library, packages), renv_status_check_synchronized(lockfile, 
    ##         library), renv_status_check_version(lockfile))
    ##     if (sources) {
    ##         synchronized <- synchronized && renv_status_check_unknown_sources(project, 
    ##             lockfile)
    ##     }
    ##     if (cache) 
    ##         renv_status_check_cache(project)
    ##     if (synchronized) 
    ##         writef("No issues found -- the project is in a consistent state.")
    ##     else writef("See `?renv::status` for advice on resolving these issues.")
    ##     result <- list(library = library, lockfile = lockfile, synchronized = synchronized)
    ##     invisible(result)
    ## }
    ## <bytecode: 0x64d250aa5680>
    ## <environment: namespace:renv>

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

# Set up environment

``` r
## Define paths to directory
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

#Create p1 results directory
p1_K_rich_domains <- file.path(results.dir, "p1_K_rich_domains")
dir.create(p1_K_rich_domains, showWarnings = FALSE, recursive = TRUE)
```

# Import data

``` r
## Load data into environment 
protein.feature.dt <- fread(file.path(data.dir, "PNAS2022/all_protein_feature_per_position.csv"))

## Load library  
library("org.Hs.eg.db")
```

    ## Loading required package: AnnotationDbi

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

    ## 

``` r
## Convert gene id and names
# Map UniProt IDs to gene symbols and ENSEMBL IDs
gene.id.dt <- select(
    org.Hs.eg.db,
    keys = protein.feature.dt[, unique(uniprot_id)],
    columns = c("SYMBOL","ENSEMBL"),
    keytype = "UNIPROT"
) %>%
    data.table
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
# Rename columns 
setnames(
    gene.id.dt,
    old = c("UNIPROT", "SYMBOL","ENSEMBL"),
    new = c("uniprot_id", "gene_name", "gene_id")
)

## Clean and organize gene mapping table
# Sort gene_id in ascending order, remove duplicate UniProt id and remove rows missing gene_id 
gene.id.dt <- gene.id.dt[order(gene_id)][!duplicated(uniprot_id)][!is.na(gene_id)]

# Rename data file and save to directory 
fwrite(gene.id.dt, file.path(results.dir, "p1_K_rich_domains","uniprot-ensembl.csv"))
```

# Maximum K score and Protein Length

``` r
#-----------------------------------
# Maximum K score and Protein Length
#-----------------------------------

## Compute maximum K_ratio and protein length for each protein
max.k.score.dt <- protein.feature.dt[ 
  , list(
    max_k_ratio = max(K_ratio), 
    protein_len = max(position)
    ), 
  by = list(Accession, uniprot_id)] 

## Merge gene annotations with protein summary data (max K_ratio, protein length)
# Merge tables using the UniProt ID as the key
max.k.score.dt <- merge(
    gene.id.dt,
    max.k.score.dt,
    all.y = TRUE,
    by = "uniprot_id"
)

# Create dataset with max k score 
fwrite(max.k.score.dt, file.path(results.dir, "p1_K_rich_domains", "max-k-score-per-protein.csv")) 
```

# Data Visualization - K score, Maximum K score and Protein Length

``` r
#---------------------------------------------------------
# Data Visualization - K score, Maximum K score and Protein Length
#---------------------------------------------------------

## K_ratio distribution 
# Count how many proteins have each K_ratio value
protein.feature.dt[, table(K_ratio)] 
```

    ## K_ratio
    ##       0     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9 
    ## 6672182 3249314 1066633  255694   52006   11627    2951     979     322      54 
    ##       1 
    ##       2

``` r
# Create a new column `K_ratio_group` based on the K_ratio
protein.feature.dt[, K_ratio_group := case_when(
  K_ratio >= 0.9 ~ 0.9, #If K_ratio is 0.9 or higher, assign value 0.9
  TRUE ~ K_ratio #If K_ratio is lower than 0.9, then assign original K_ratio value
)]

# Plot - The number of proteins per K score
ggplot(
    protein.feature.dt,
    aes(
        x = K_ratio_group
    )
) +
    geom_bar() +
    scale_x_continuous(breaks=seq(0, 1, 0.1))
```

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-1.png)<!-- -->

``` r
# Plot - The number of proteins per max K score
ggplot(
    max.k.score.dt,
    aes(
        x = max_k_ratio
    )
) +
    geom_bar() +
    scale_x_continuous(breaks=seq(0, 1, 0.1)) + 
  labs(x = "max_k_ratio", y = "number of protein_count", title = "Number of proteins per max K score")
```

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-2.png)<!-- -->

``` r
# Plot - Max K ratio vs protein length
ggplot(
  max.k.score.dt,
  aes(
    y = protein_len,
    x = factor(max_k_ratio)
  )
) + geom_boxplot(fill = "steelblue", outlier.shape = NA) +
  labs(x = "max_K_ratio", y = "protein_length", title = "Max K ratio vs protein length") +
  coord_cartesian(ylim = c(0, 2000))
```

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-3.png)<!-- -->

``` r
# Plot - K_score of disordered regions of proteins 
# This plot indicates K score depending on the level of protein disorderliness

ggplot(
  protein.feature.dt,
  aes(
    x = factor(K_ratio_group),
    y = IUPRED2
  )
) + 
  geom_boxplot(fill = "steelblue", outlier.shape = NA)
```

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-4.png)<!-- -->

``` r
# Plot - Protein charge versus K_ratio_group
ggplot(
  protein.feature.dt,
  aes(
    x = factor(K_ratio_group),
    y = windowCharge
  )
) + 
  geom_hline(yintercept = 0, color = "gray60") +
  geom_boxplot(fill = "steelblue", outlier.shape = NA)
```

    ## Warning: Removed 211124 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-5.png)<!-- -->

# RNA-Binding proteins

``` r
#---------------------------
# RNA-Binding Proteins (RBP)
#---------------------------

# Load data into environment 
rbp.dt <- readxl::read_excel(file.path(data.dir, "20241113_RBPbase_Hs_DescriptiveID.xlsx")) %>%
  data.table()

# Remove line breaks and split into two columns
setnames(rbp.dt, old = colnames(rbp.dt), new = str_split_fixed(colnames(rbp.dt), "\n", n = 2)[, 1])
rbp.dt <- janitor::clean_names(rbp.dt, case = "none") #Rename duplicate column names
setnames(rbp.dt, old = "UnitProtSwissProtID_Hs", new = "uniprot_id") 


# Dataset containing only UniProt id and Hs_ columns
rbp.dt <- rbp.dt[
  , c("uniprot_id", grep("^Hs_", colnames(rbp.dt), value = TRUE)), with = FALSE
  ]

# Remove NA and duplicated values
rbp.dt <- rbp.dt[!is.na(uniprot_id) & !duplicated(uniprot_id)]

# Reshape dataset into long format for easy data management and visualization 
m.rbp.dt <- melt(
    rbp.dt,
    id.vars = c("uniprot_id"),
    variable.name = "dataset",
    value.name = "RBP"
)

# Dataset containing count of identified RBP and total RBP 
rbp.count.dt <- m.rbp.dt[!is.na(uniprot_id),
  list(
    identified_RBP = sum(RBP == "YES"),
    total = sum(RBP == "YES") + sum(RBP == "no") 
    ), by = uniprot_id
  ]

rbp.count.dt[, table(identified_RBP)]
```

    ## identified_RBP
    ##     0     1     2     3     4     5     6     7     8     9    10    11    12 
    ## 14104  1936   694   388   278   205   174   126   104    98    86    68    65 
    ##    13    14    15    16    17    18    19    20    21    22    23    24    25 
    ##    49    46    50    52    49    60    41    55    26    21    20    29    13 
    ##    26    27 
    ##     9     6

``` r
#-----------------------------------
## Comparing Max K score versus RBP
#-----------------------------------

# Create dataset containing max K score and RBP count
rbp.count.dt <- merge(
  max.k.score.dt, rbp.count.dt, by = "uniprot_id"
)

data.set.names <- c(
  "Others", 
  "RBP (1 dataset)", "RBP (2 dataset)", "RBP (> 2 dataset)"
  )

# Assign values to RBP dataset columns based on RBP count 
rbp.count.dt[, RBP := case_when(
  identified_RBP > 2 ~ "RBP (> 2 dataset)",
  identified_RBP > 1 ~ "RBP (2 dataset)",
  identified_RBP > 0 ~ "RBP (1 dataset)", 
  TRUE ~ "Others") %>% factor(levels = data.set.names)]
rbp.count.dt[, table(RBP)]
```

    ## RBP
    ##            Others   RBP (1 dataset)   RBP (2 dataset) RBP (> 2 dataset) 
    ##             14045              1928               694              2111

# Data Visualization - RNA-Binding Proteins (Max K score versus RBP)

``` r
# Plot - Max K score versus RBP
max_K_score_rbp <- ggplot(
  data = rbp.count.dt,
  aes(
    x = RBP,
    y = max_k_ratio
    )
  ) +
  geom_boxplot(fill = "steelblue", outlier.shape = NA) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme(aspect.ratio = 3)
```

# Subcellular localisation

``` r
#--------------------------------------
# Sub-cellular localisation analysis
#--------------------------------------

# Packaged installed 20th Oct 2025
# install.packages(file.path("/fast/AG_Sugimoto/home/users/pallavi/Software_packages/R/subcellularvis"), repos = NULL, type = "source")

## Load library
library("subcellularvis")

## Run sub-cellular compartment analysis with enrichment score according Gene Ontology
# The 'extractCompartment' function runs sub-cellular compartment analysis, retrieves 
extractCompartment <- function(gene.set, class.name){
  comp.out <- compartmentData(  # Run compartment enrichment analysis using the compartmentData() function
    gene.set,
    id_type = "UNIPROT",
    aspect = c("Whole cell"),
    organism = c("Human"),
    annotationSource =c("Gene Ontology") # Alternatively, Human Protein Atlas can also be used
  )
  
  # Convert enrichment results from comp.out into a data.table 
  comp.out.dt <- data.table(comp.out$enrichment)
  
  # Add metadata columns 
  comp.out.dt[, `:=`(
    Genes = NA,
    n_unmapped = length(str_split(comp.out$unmapped, ",")[[1]]), # Counting the number of unmapped IDs
    n_mapped = comp.out$nMapped, # adding mapped IDs
    class_name = class.name 
  )]
  
  #  Updated data table 
  return(comp.out.dt)
}

# Add column containing values from 0 - 0.9 only
max.k.score.dt[, max_k_ratio_group := case_when(
  max_k_ratio >= 0.9 ~ 0.9, # if ≥ 0.9 → set value to 0.9
  TRUE ~ max_k_ratio  #if < 0.9, then retain original max k ratio value
)]

# Apply the 'extractCompartment' function to each subset of protein grouped by
# max_K_ratio
k.score.compartment.dt <- mapply(
  extractCompartment,
  gene.set = lapply(
    seq(0, 9, by = 1) / 10, # create a sequence of 0 to 0.9 
    function(x){ # For each x, extract the 'uniport_id' belonging to that max_k_ratio_group
      max.k.score.dt[max_k_ratio_group == x, uniprot_id]}
    ),
  class.name = as.character(seq(0, 0.9, by = 0.1)), # convert and assign class names as strings
  SIMPLIFY = FALSE # Return list instead of matrix or array
) %>% rbindlist # combine all lists into a data table
  
# Calculate the negative log10 of FDR and create a new column called 'mlog10FDR'
k.score.compartment.dt[, mlog10FDR := -log10(FDR)]

# Convert data table into wide format 
# rows: 'Compartment'
# column: 'class_name'
# values: 'mlog10FDR'
d.k.score.compartment.dt <- dcast(
  k.score.compartment.dt,
  Compartment ~ factor(class_name),
  value.var = "mlog10FDR"
)

# convert data table into matrix 
d.k.score.compartment.mat <- as.matrix(d.k.score.compartment.dt[, 2:ncol(d.k.score.compartment.dt)])

# Assign Compartment names as rownames of the matrix 
rownames(d.k.score.compartment.mat) <- d.k.score.compartment.dt[, Compartment]
```

# Data Visualization - Subcellular Localisation

``` r
# Load libraries -  'RColorBrewer' and 'pheatmap'
library(RColorBrewer)
library("pheatmap")

mat_breaks <- seq(min(d.k.score.compartment.mat), max(d.k.score.compartment.mat), length.out = 40)

# Plot - heatmap of subcellular localisation
pheatmap(
  d.k.score.compartment.mat,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "steelblue"))(length(mat_breaks)),
  breaks = mat_breaks
)
```

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Subcellular_Localisation-1.png)<!-- -->

``` r
# Plot - max_k_score versus −log10(FDR) of nuclear localisation enrichment
ggplot(
  data = k.score.compartment.dt[Compartment == "Nucleus"],
  aes(
    x = factor(class_name),
    y = mlog10FDR
  )
) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), color = "gray60") +
  ylab("-log10(FDR) of nuclear localisation enrichment") +
  xlab("Maximum K score of protein")
```

![](p1-1_lysine-rich-domain-proteins_files/figure-gfm/Data_Visualization_Subcellular_Localisation-2.png)<!-- -->

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
    ##  date     2025-11-03
    ##  pandoc   3.4 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.6.42 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package          * version    date (UTC) lib source
    ##  AnnotationDbi    * 1.70.0     2025-04-15 [1] Bioconduc~
    ##  Biobase          * 2.68.0     2025-04-15 [1] Bioconduc~
    ##  BiocGenerics     * 0.54.0     2025-04-15 [1] Bioconduc~
    ##  Biostrings       * 2.76.0     2025-04-15 [1] Bioconduc~
    ##  bit                4.6.0      2025-03-06 [1] CRAN (R 4.5.1)
    ##  bit64              4.6.0-1    2025-01-16 [1] CRAN (R 4.5.1)
    ##  blob               1.2.4      2023-03-17 [1] CRAN (R 4.5.1)
    ##  cachem             1.1.0      2024-05-16 [1] CRAN (R 4.5.1)
    ##  cellranger         1.1.0      2016-07-27 [1] CRAN (R 4.5.1)
    ##  cli                3.6.5      2025-04-23 [1] CRAN (R 4.5.1)
    ##  colourpicker       1.3.0      2023-08-21 [1] CRAN (R 4.5.1)
    ##  crayon             1.5.3      2024-06-20 [1] CRAN (R 4.5.1)
    ##  data.table       * 1.17.8     2025-07-10 [1] CRAN (R 4.5.1)
    ##  DBI                1.2.3      2024-06-02 [1] CRAN (R 4.5.1)
    ##  digest             0.6.37     2024-08-19 [1] CRAN (R 4.5.1)
    ##  dplyr            * 1.1.4      2023-11-17 [1] CRAN (R 4.5.1)
    ##  evaluate           1.0.5      2025-08-27 [1] CRAN (R 4.5.1)
    ##  farver             2.1.2      2024-05-13 [1] CRAN (R 4.5.1)
    ##  fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.5.1)
    ##  formattable        0.2.1      2021-01-07 [1] CRAN (R 4.5.1)
    ##  generics         * 0.1.4      2025-05-09 [1] CRAN (R 4.5.1)
    ##  GenomeInfoDb     * 1.44.3     2025-09-21 [1] Bioconduc~
    ##  GenomeInfoDbData   1.2.14     2025-09-24 [1] Bioconductor
    ##  ggplot2          * 4.0.0      2025-09-11 [1] CRAN (R 4.5.1)
    ##  glue               1.8.0      2024-09-30 [1] CRAN (R 4.5.1)
    ##  gridExtra          2.3        2017-09-09 [1] CRAN (R 4.5.1)
    ##  gtable             0.3.6      2024-10-25 [1] CRAN (R 4.5.1)
    ##  htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.5.1)
    ##  htmlwidgets        1.6.4      2023-12-06 [1] CRAN (R 4.5.1)
    ##  httpuv             1.6.16     2025-04-16 [1] CRAN (R 4.5.1)
    ##  httr               1.4.7      2023-08-15 [1] CRAN (R 4.5.1)
    ##  IRanges          * 2.42.0     2025-04-15 [1] Bioconduc~
    ##  janitor            2.2.1      2024-12-22 [1] CRAN (R 4.5.1)
    ##  jsonlite           2.0.0      2025-03-27 [1] CRAN (R 4.5.1)
    ##  KEGGREST           1.48.1     2025-06-22 [1] Bioconduc~
    ##  khroma           * 1.16.0     2025-02-25 [1] CRAN (R 4.5.1)
    ##  knitr            * 1.50       2025-03-16 [1] CRAN (R 4.5.1)
    ##  labeling           0.4.3      2023-08-29 [1] CRAN (R 4.5.1)
    ##  later              1.4.4      2025-08-27 [1] CRAN (R 4.5.1)
    ##  lazyeval           0.2.2      2019-03-15 [1] CRAN (R 4.5.1)
    ##  lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.5.1)
    ##  lubridate          1.9.4      2024-12-08 [1] CRAN (R 4.5.1)
    ##  magrittr         * 2.0.4      2025-09-12 [1] CRAN (R 4.5.1)
    ##  memoise            2.0.1      2021-11-26 [1] CRAN (R 4.5.1)
    ##  mime               0.13       2025-03-17 [1] CRAN (R 4.5.1)
    ##  miniUI             0.1.2      2025-04-17 [1] CRAN (R 4.5.1)
    ##  org.Hs.eg.db     * 3.21.0     2025-10-17 [1] Bioconductor
    ##  pheatmap         * 1.0.13     2025-06-05 [1] CRAN (R 4.5.1)
    ##  pillar             1.11.1     2025-09-17 [1] CRAN (R 4.5.1)
    ##  pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.5.1)
    ##  plotly             4.11.0     2025-06-19 [1] CRAN (R 4.5.1)
    ##  plyr               1.8.9      2023-10-02 [1] CRAN (R 4.5.1)
    ##  png                0.1-8      2022-11-29 [1] CRAN (R 4.5.1)
    ##  promises           1.3.3      2025-05-29 [1] CRAN (R 4.5.1)
    ##  purrr              1.1.0      2025-07-10 [1] CRAN (R 4.5.1)
    ##  R6                 2.6.1      2025-02-15 [1] CRAN (R 4.5.1)
    ##  RColorBrewer     * 1.1-3      2022-04-03 [1] CRAN (R 4.5.1)
    ##  Rcpp               1.1.0      2025-07-02 [1] CRAN (R 4.5.1)
    ##  readxl             1.4.5      2025-03-07 [1] CRAN (R 4.5.1)
    ##  renv               1.1.5      2025-07-24 [1] CRAN (R 4.5.1)
    ##  rlang              1.1.6      2025-04-11 [1] CRAN (R 4.5.1)
    ##  rmarkdown          2.29       2024-11-04 [1] CRAN (R 4.5.1)
    ##  RSQLite            2.4.3      2025-08-20 [1] CRAN (R 4.5.1)
    ##  S4Vectors        * 0.46.0     2025-04-15 [1] Bioconduc~
    ##  S7                 0.2.0      2024-11-07 [1] CRAN (R 4.5.1)
    ##  scales             1.4.0      2025-04-24 [1] CRAN (R 4.5.1)
    ##  sessioninfo        1.2.3      2025-02-05 [1] CRAN (R 4.5.1)
    ##  shiny              1.11.1     2025-07-03 [1] CRAN (R 4.5.1)
    ##  shinythemes        1.2.0      2021-01-25 [1] CRAN (R 4.5.1)
    ##  snakecase          0.11.1     2023-08-27 [1] CRAN (R 4.5.1)
    ##  stringi            1.8.7      2025-03-27 [1] CRAN (R 4.5.1)
    ##  stringr          * 1.5.2      2025-09-08 [1] CRAN (R 4.5.1)
    ##  subcellularvis   * 0.0.0.9000 2025-10-20 [1] local
    ##  tibble             3.3.0      2025-06-08 [1] CRAN (R 4.5.1)
    ##  tidyr              1.3.1      2024-01-24 [1] CRAN (R 4.5.1)
    ##  tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.5.1)
    ##  timechange         0.3.0      2024-01-18 [1] CRAN (R 4.5.1)
    ##  UCSC.utils         1.4.0      2025-04-15 [1] Bioconduc~
    ##  UpSetR             1.4.0      2019-05-22 [1] CRAN (R 4.5.1)
    ##  vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.5.1)
    ##  viridisLite        0.4.2      2023-05-02 [1] CRAN (R 4.5.1)
    ##  withr              3.0.2      2024-10-28 [1] CRAN (R 4.5.1)
    ##  xfun               0.53       2025-08-19 [1] CRAN (R 4.5.1)
    ##  xtable             1.8-4      2019-04-21 [1] CRAN (R 4.5.1)
    ##  XVector          * 0.48.0     2025-04-15 [1] Bioconduc~
    ##  yaml               2.3.10     2024-07-26 [1] CRAN (R 4.5.1)
    ## 
    ##  [1] /home/pkesava/R/x86_64-pc-linux-gnu-library/4.5
    ##  [2] /usr/local/lib/R/site-library
    ##  [3] /usr/lib/R/site-library
    ##  [4] /usr/lib/R/library
    ##  * ── Packages attached to the search path.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
