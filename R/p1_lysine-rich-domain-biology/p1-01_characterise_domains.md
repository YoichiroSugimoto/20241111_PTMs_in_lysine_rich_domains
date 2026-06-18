p1-01 · Characterise lysine-rich domains
================
Yoichiro Sugimoto and Pallavi Kesavan
17 June, 2026

- [Overview](#overview)
- [Setup](#setup)
- [Import data](#import-data)
- [Maximum K score and protein
  length](#maximum-k-score-and-protein-length)
- [K score distributions](#k-score-distributions)
- [Subcellular localisation](#subcellular-localisation)
- [Functional class enrichment](#functional-class-enrichment)
- [K score in histones](#k-score-in-histones)
- [K score in disordered regions](#k-score-in-disordered-regions)
- [Session information](#session-information)

# Overview

**Purpose:** Characterise lysine-rich domains across the human proteome
— maximum K-score vs protein length, subcellular localisation,
functional-class enrichment, and K-scores in histones.

**Inputs:**
`data/processed_data_from_PNAS2022/all_protein_feature_per_position.csv`,
the cd-code annotation table under `data/public_data/`, and the UniProt
reference proteome; annotation packages `org.Hs.eg.db` and
`subcellularvis` (loaded where used).

**Outputs:** figures (rendered on knit).

**Upstream:** none. **Downstream:** none.

# Setup

``` r
## Resolve the repository root (via the .here sentinel) and load the shared
## setup: packages, helper functions, ggplot/knitr settings, and project paths.
repo_root <- local({
  p <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  while (!file.exists(file.path(p, ".here")) && !identical(dirname(p), p)) p <- dirname(p)
  p
})
source(file.path(repo_root, "R", "functions", "_setup.R"))
```

# Import data

Loads the per-position protein feature table, maps UniProt IDs to gene
symbols, and prepares the K-score data used throughout.

``` r
## Define paths to directory
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

#Create p1 results directory
p1.results.dir <- file.path(results.dir, "p1_characterise_domains")
create.dirs(c(results.dir, p1.results.dir))

## Load data into environment 
protein.feature.dt <- fread(file.path(data.dir, "processed_data_from_PNAS2022/all_protein_feature_per_position.csv"))

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
fwrite(gene.id.dt, file.path(p1.results.dir,"uniprot-ensembl.csv"))
```

# Maximum K score and protein length

Computes each protein’s maximum K score (lysine richness) and its
length.

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
fwrite(max.k.score.dt, file.path(p1.results.dir, "max-k-score-per-protein.csv")) 

print("The ratio of proteins with the max K ratio > 0.3")
```

    ## [1] "The ratio of proteins with the max K ratio > 0.3"

``` r
nrow(max.k.score.dt[max_k_ratio > 0.3]) / nrow(max.k.score.dt)
```

    ## [1] 0.3809687

``` r
# Median K score of all region in human proteome#
protein.feature.dt[, median(K_ratio)]
```

    ## [1] 0

# K score distributions

Visualises the distribution of K scores and the relationship between
maximum K score and protein length.

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

![](p1-01_characterise_domains_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-1.png)<!-- -->

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
  labs(x = "max_k_ratio", y = "number of protein_count")
```

![](p1-01_characterise_domains_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-2.png)<!-- -->

``` r
#Median - 24112025
median(max.k.score.dt[, max_k_ratio])
```

    ## [1] 0.3

``` r
# 0.3

# Plot - Max K ratio vs protein length
ggplot(
  max.k.score.dt,
  aes(
    y = protein_len,
    x = factor(max_k_ratio)
  )
) + geom_boxplot(fill = "steelblue", outlier.shape = NA) +
  labs(x = "max_K_ratio", y = "protein_length") +
  coord_cartesian(ylim = c(0, 2000))
```

![](p1-01_characterise_domains_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-3.png)<!-- -->

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

![](p1-01_characterise_domains_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-4.png)<!-- -->

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

![](p1-01_characterise_domains_files/figure-gfm/Data_Visualization_Max_K_Score_&_Protein_Length-5.png)<!-- -->

# Subcellular localisation

Annotates proteins with subcellular localisation to relate lysine
richness to cellular compartment.

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

# Functional class enrichment

Tests which functional protein classes are enriched among lysine-rich
proteins.

``` r
library("mgcv")
```

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     collapse

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     collapse

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## This is mgcv 1.9-1. For overview type 'help("mgcv-package")'.

``` r
library("subcellularvis")

comp.out <- compartmentData(
    max.k.score.dt[, uniprot_id],
    id_type = "UNIPROT",
    aspect = c("Whole cell"),
    organism = c("Human"),
    annotationSource =c("Gene Ontology") #c("Human Protein Atlas")
  )

localisation_coef = function(i){
  localised_genes <- comp.out$enrichment[i, "Genes"] %>% {strsplit(., split = ",", fixed = TRUE)[[1]]}
  max.k.score.dt[, loc := uniprot_id %in% localised_genes]
  len.gam <- gam(max_k_ratio ~ s(protein_len) + loc, data = max.k.score.dt, method = "REML")
  
  dt <- data.table(
    functional_type = "localisation", 
    subtype = comp.out$enrichment[, "Compartment"][i],
    coef = summary(len.gam)$p.coeff["locTRUE"],
    se = summary(len.gam)$se["locTRUE"],
    p_value = summary(len.gam)$p.pv["locTRUE"],
    n = length(localised_genes)
  )
  
  return(dt)
}

localisation_dt <- lapply(1:length(comp.out$enrichment[, "Compartment"]), localisation_coef) %>%
  rbindlist

# Molecular condensates
cd_code <- fread(file.path(data.dir, "public_data/20241113_cd-code.csv"))

condensates_coef <- function(i){  
  condensate_genes <- cd_code[i, Proteins] %>% {strsplit(., split = "\t", fixed = TRUE)[[1]]}
  max.k.score.dt[, loc := uniprot_id %in% condensate_genes]
  len.gam <- gam(max_k_ratio ~ s(protein_len) + loc, data = max.k.score.dt, method = "REML")
  
  dt <- data.table(
    functional_type = "condensates", 
    subtype = cd_code[i, Name],
    coef = summary(len.gam)$p.coeff["locTRUE"],
    se = summary(len.gam)$se["locTRUE"],
    p_value = summary(len.gam)$p.pv["locTRUE"],
    n = length(condensate_genes)
  )
  
  return(dt)
}

condensate_dt <- lapply(1:nrow(cd_code), condensates_coef) %>%
  rbindlist

condensate_dt <- condensate_dt[order(p_value)][1:5]

# histone
max.k.score.dt[, histone_protein := grepl("^H[1-4]", gene_name)]
len.gam <- gam(max_k_ratio ~ s(protein_len) + histone_protein, data = max.k.score.dt, method = "REML")
  
histone_dt <- data.table(
  functional_type = "localisation", 
  subtype = "Histone",
  coef = summary(len.gam)$p.coeff["histone_proteinTRUE"],
  se = summary(len.gam)$se["histone_proteinTRUE"],
  p_value = summary(len.gam)$p.pv["histone_proteinTRUE"],
  n = nrow(max.k.score.dt[histone_protein == TRUE])
)

# Merge all three data
functional_class_summary <- rbindlist(list(
  localisation_dt, histone_dt, condensate_dt
))

functional_class_summary <- functional_class_summary[order(-functional_type, -coef)]
functional_class_summary[, `:=`(
  functional_type = factor(functional_type, levels = c("localisation", "condensates")),
  subtype = factor(subtype, levels = subtype)
)]

functional_class_summary
```

    ##     functional_type               subtype         coef          se
    ##              <fctr>                <fctr>        <num>       <num>
    ##  1:    localisation               Histone  0.163573821 0.015153259
    ##  2:    localisation              Ribosome  0.094719538 0.007785273
    ##  3:    localisation               Nucleus  0.044130892 0.001641308
    ##  4:    localisation          Cytoskeleton  0.008931795 0.002546917
    ##  5:    localisation             Cytoplasm  0.004000254 0.001617667
    ##  6:    localisation         Mitochondrion -0.002150944 0.002939521
    ##  7:    localisation Endoplasmic reticulum -0.010155542 0.002687986
    ##  8:    localisation Intracellular vesicle -0.011941747 0.002451759
    ##  9:    localisation  Extracellular region -0.013043133 0.001990705
    ## 10:    localisation              Endosome -0.017430567 0.003699277
    ## 11:    localisation       Plasma membrane -0.018304547 0.001817210
    ## 12:    localisation       Golgi apparatus -0.018694085 0.002980457
    ## 13:    localisation            Peroxisome -0.019337922 0.009531725
    ## 14:    localisation               Vacuole -0.029113620 0.004065612
    ## 15:    localisation              Lysosome -0.031381737 0.004308968
    ## 16:     condensates            Cajal body  0.091028126 0.013163794
    ## 17:     condensates             Nucleolus  0.076694198 0.003220552
    ## 18:     condensates       Nuclear speckle  0.054212632 0.007903436
    ## 19:     condensates                P-body  0.048065116 0.005231583
    ## 20:     condensates        Stress granule  0.034788726 0.003997924
    ##     functional_type               subtype         coef          se
    ##           p_value     n
    ##             <num> <int>
    ##  1:  4.320651e-27    56
    ##  2:  6.148159e-34   214
    ##  3: 1.643062e-156  7329
    ##  4:  4.543272e-04  2276
    ##  5:  1.341198e-02 11635
    ##  6:  4.643401e-01  1618
    ##  7:  1.584649e-04  1965
    ##  8:  1.120462e-06  2421
    ##  9:  5.812286e-11  4033
    ## 10:  2.470576e-06   984
    ## 11:  8.278909e-24  5351
    ## 12:  3.630958e-10  1566
    ## 13:  4.249208e-02   142
    ## 14:  8.285065e-13   806
    ## 15:  3.386400e-13   714
    ## 16:  4.815930e-12    74
    ## 17: 1.169497e-123  1411
    ## 18:  7.114597e-12   207
    ## 19:  4.396641e-20   536
    ## 20:  3.513771e-18   905
    ##           p_value     n

``` r
ggplot(
  data = functional_class_summary,
  aes(
    x = subtype,
    y = coef
  )) +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_point() +
  geom_errorbar(aes(
    ymin = coef - se, ymax = coef + se, width = 0.3
  )) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(~ functional_type, scales = "free", space = "free")
```

![](p1-01_characterise_domains_files/figure-gfm/functional_classes-1.png)<!-- -->

``` r
ggplot(
  data = max.k.score.dt,
  aes(
    x = factor(max_k_ratio),
    y = protein_len,
    fill = histone_protein
  )
) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(preserve = "single")) +
  scale_fill_bright() +
  coord_cartesian(ylim = c(0, 2500))
```

![](p1-01_characterise_domains_files/figure-gfm/functional_classes-2.png)<!-- -->

# K score in histones

Examines histones — well-defined lysine-rich proteins — as a reference
point.

The definition of histone tail was obtained from “Histone
post-translational modifications — cause and consequence of genome
function (<https://doi.org/10.1038/s41576-022-00468-7>)”.

``` r
#H2A
ggplot(
  protein.feature.dt[grepl("Q6FI13", Accession)],
  aes(
    x = position,
    y = K_ratio_score
  )
) +
  coord_cartesian(ylim = c(0,0.5)) +
  annotate("rect", xmin = 2, xmax = 26, ymin = 0, ymax = 1, fill = "red", alpha = 0.2) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("Q6FI13", Accession)][1, Accession])
```

![](p1-01_characterise_domains_files/figure-gfm/histone_k_score_per_position-1.png)<!-- -->

``` r
#H2B
ggplot(
  protein.feature.dt[grepl("P33778", Accession)],
  aes(
    x = position,
    y = K_ratio_score
  )
) +
  coord_cartesian(ylim = c(0,0.5)) +
  annotate("rect", xmin = 2, xmax = 33, ymin = 0, ymax = 1, fill = "red", alpha = 0.2) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("P33778", Accession)][1, Accession])
```

![](p1-01_characterise_domains_files/figure-gfm/histone_k_score_per_position-2.png)<!-- -->

``` r
#H3
ggplot(
  protein.feature.dt[grepl("Q16695", Accession)],
  aes(
    x = position,
    y = K_ratio_score
  )
) +
  coord_cartesian(ylim = c(0,0.5)) +
  annotate("rect", xmin = 2, xmax = 44, ymin = 0, ymax = 1, fill = "red", alpha = 0.2) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("Q16695", Accession)][1, Accession])
```

![](p1-01_characterise_domains_files/figure-gfm/histone_k_score_per_position-3.png)<!-- -->

``` r
#H4
ggplot(
  protein.feature.dt[grepl("P62805", Accession)],
  aes(
    x = position,
    y = K_ratio_score
  )
) +
  coord_cartesian(ylim = c(0,0.5)) +
  annotate("rect", xmin = 2, xmax = 24, ymin = 0, ymax = 1, fill = "red", alpha = 0.2) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("P62805", Accession)][1, Accession])
```

![](p1-01_characterise_domains_files/figure-gfm/histone_k_score_per_position-4.png)<!-- -->

``` r
print("Median K score of histone tails")
```

    ## [1] "Median K score of histone tails"

``` r
rbindlist(list(
  protein.feature.dt[grepl("Q6FI13", Accession)][2:26],
  protein.feature.dt[grepl("P33778", Accession)][2:33],
  protein.feature.dt[grepl("Q16695", Accession)][2:44],
  protein.feature.dt[grepl("P62805", Accession)][2:24]
))[, median(K_ratio_score)]
```

    ## [1] 0.3

# K score in disordered regions

Relates lysine richness to protein disorder.

``` r
library("patchwork")

dkc1.1 <- ggplot(
  protein.feature.dt[grepl("O60832", Accession)],
  aes(
    x = position,
    y = K_ratio_score
  )
) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("O60832", Accession)][1, Accession]) +
  geom_point(x = 475, y = 0, color = "red") +
  geom_point(data = data.table(x = c(501:505), y = rep(0, 5)), aes(x = x, y = y), color = "blue")

dkc1.2 <- ggplot(
  protein.feature.dt[grepl("O60832", Accession)],
  aes(
    x = position,
    y = IUPRED2
  )
) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("O60832", Accession)][1, Accession])

dkc1.1 / dkc1.2
```

![](p1-01_characterise_domains_files/figure-gfm/kscore_disorderedness_other_proteins-1.png)<!-- -->

``` r
ARL6IP4.1 <- ggplot(
  protein.feature.dt[grepl("Q66PJ3", Accession)],
  aes(
    x = position,
    y = K_ratio_score
  )
) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("Q66PJ3", Accession)][1, Accession]) +
  geom_point(data = data.table(x = c(290, 292, 294), y = rep(0, 3)), aes(x = x, y = y), color = "red") +
  geom_point(data = data.table(x = c(290, 292, 294, 303, 304, 306, 307), y = rep(0, 7)), aes(x = x, y = y), color = "blue")

ARL6IP4.2 <- ggplot(
  protein.feature.dt[grepl("Q66PJ3", Accession)],
  aes(
    x = position,
    y = IUPRED2
  )
) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_area(fill = "gray40") +
  ggtitle(protein.feature.dt[grepl("Q66PJ3", Accession)][1, Accession])

ARL6IP4.1 / ARL6IP4.2
```

![](p1-01_characterise_domains_files/figure-gfm/kscore_disorderedness_other_proteins-2.png)<!-- -->

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
    ##  date     2026-06-17
    ##  pandoc   3.2 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64/ (via rmarkdown)
    ##  quarto   1.5.57 @ /usr/lib/rstudio-server/bin/quarto/bin/quarto
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package           * version    date (UTC) lib source
    ##  AnnotationDbi     * 1.68.0     2024-10-29 [1] Bioconduc~
    ##  Biobase           * 2.66.0     2024-10-29 [1] Bioconduc~
    ##  BiocGenerics      * 0.52.0     2024-10-29 [1] Bioconduc~
    ##  Biostrings        * 2.74.1     2024-12-16 [1] Bioconduc~
    ##  bit                 4.6.0      2025-03-06 [1] CRAN (R 4.4.3)
    ##  bit64               4.6.0-1    2025-01-16 [1] CRAN (R 4.4.3)
    ##  blob                1.2.4      2023-03-17 [1] CRAN (R 4.4.3)
    ##  cachem              1.1.0      2024-05-16 [1] CRAN (R 4.4.3)
    ##  cellranger          1.1.0      2016-07-27 [1] CRAN (R 4.4.3)
    ##  cli                 3.6.5      2025-04-23 [1] CRAN (R 4.4.3)
    ##  colourpicker        1.3.0      2023-08-21 [1] CRAN (R 4.4.3)
    ##  crayon              1.5.3      2024-06-20 [1] CRAN (R 4.4.3)
    ##  data.table        * 1.17.8     2025-07-10 [1] CRAN (R 4.4.3)
    ##  DBI                 1.2.3      2024-06-02 [1] CRAN (R 4.4.3)
    ##  digest              0.6.37     2024-08-19 [1] CRAN (R 4.4.3)
    ##  dplyr             * 1.1.4      2023-11-17 [1] CRAN (R 4.4.3)
    ##  evaluate            1.0.5      2025-08-27 [1] CRAN (R 4.4.3)
    ##  farver              2.1.2      2024-05-13 [1] CRAN (R 4.4.3)
    ##  fastmap             1.2.0      2024-05-15 [1] CRAN (R 4.4.3)
    ##  formattable         0.2.1      2021-01-07 [1] CRAN (R 4.4.3)
    ##  generics            0.1.4      2025-05-09 [1] CRAN (R 4.4.3)
    ##  GenomeInfoDb      * 1.42.3     2025-01-27 [1] Bioconduc~
    ##  GenomeInfoDbData    1.2.13     2025-07-21 [1] Bioconductor
    ##  ggplot2           * 4.0.0      2025-09-11 [1] CRAN (R 4.4.3)
    ##  glue                1.8.0      2024-09-30 [1] CRAN (R 4.4.3)
    ##  gridExtra           2.3        2017-09-09 [1] CRAN (R 4.4.3)
    ##  gtable              0.3.6      2024-10-25 [1] CRAN (R 4.4.3)
    ##  htmltools           0.5.8.1    2024-04-04 [1] CRAN (R 4.4.3)
    ##  htmlwidgets         1.6.4      2023-12-06 [1] CRAN (R 4.4.3)
    ##  httpuv              1.6.16     2025-04-16 [1] CRAN (R 4.4.3)
    ##  httr                1.4.7      2023-08-15 [1] CRAN (R 4.4.3)
    ##  IRanges           * 2.40.1     2024-12-05 [1] Bioconduc~
    ##  janitor             2.2.1      2024-12-22 [1] CRAN (R 4.4.3)
    ##  jsonlite            2.0.0      2025-03-27 [1] CRAN (R 4.4.3)
    ##  KEGGREST            1.46.0     2024-10-29 [1] Bioconduc~
    ##  khroma            * 1.16.0     2025-02-25 [1] CRAN (R 4.4.3)
    ##  knitr             * 1.50       2025-03-16 [1] CRAN (R 4.4.3)
    ##  labeling            0.4.3      2023-08-29 [1] CRAN (R 4.4.3)
    ##  later               1.4.2      2025-04-08 [1] CRAN (R 4.4.3)
    ##  lattice             0.22-5     2023-10-24 [4] CRAN (R 4.3.3)
    ##  lazyeval            0.2.2      2019-03-15 [1] CRAN (R 4.4.3)
    ##  lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.4.3)
    ##  lubridate           1.9.4      2024-12-08 [1] CRAN (R 4.4.3)
    ##  magrittr          * 2.0.4      2025-09-12 [1] CRAN (R 4.4.3)
    ##  Matrix              1.7-3      2025-03-11 [1] CRAN (R 4.4.3)
    ##  memoise             2.0.1      2021-11-26 [1] CRAN (R 4.4.3)
    ##  mgcv              * 1.9-1      2023-12-21 [1] CRAN (R 4.4.3)
    ##  mime                0.13       2025-03-17 [1] CRAN (R 4.4.3)
    ##  miniUI              0.1.2      2025-04-17 [1] CRAN (R 4.4.3)
    ##  nlme              * 3.1-168    2025-03-31 [4] CRAN (R 4.4.3)
    ##  org.Hs.eg.db      * 3.20.0     2025-10-31 [1] Bioconductor
    ##  patchwork         * 1.3.2      2025-08-25 [1] CRAN (R 4.4.3)
    ##  pillar              1.11.1     2025-09-17 [1] CRAN (R 4.4.3)
    ##  pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.4.3)
    ##  plotly              4.11.0     2025-06-19 [1] CRAN (R 4.4.3)
    ##  plyr                1.8.9      2023-10-02 [1] CRAN (R 4.4.3)
    ##  png                 0.1-8      2022-11-29 [1] CRAN (R 4.4.3)
    ##  promises            1.3.3      2025-05-29 [1] CRAN (R 4.4.3)
    ##  ptm.stoichiometry * 0.0.0.9000 2025-12-13 [1] local
    ##  purrr               1.1.0      2025-07-10 [1] CRAN (R 4.4.3)
    ##  R6                  2.6.1      2025-02-15 [1] CRAN (R 4.4.3)
    ##  RColorBrewer        1.1-3      2022-04-03 [1] CRAN (R 4.4.3)
    ##  Rcpp                1.1.0      2025-07-02 [1] CRAN (R 4.4.3)
    ##  readxl            * 1.4.5      2025-03-07 [1] CRAN (R 4.4.3)
    ##  rlang               1.1.6      2025-04-11 [1] CRAN (R 4.4.3)
    ##  rmarkdown           2.29       2024-11-04 [1] CRAN (R 4.4.3)
    ##  RSQLite             2.4.3      2025-08-20 [1] CRAN (R 4.4.3)
    ##  rstudioapi          0.17.1     2024-10-22 [1] CRAN (R 4.4.3)
    ##  S4Vectors         * 0.44.0     2024-10-29 [1] Bioconduc~
    ##  S7                  0.2.0      2024-11-07 [1] CRAN (R 4.4.3)
    ##  scales              1.4.0      2025-04-24 [1] CRAN (R 4.4.3)
    ##  sessioninfo         1.2.3      2025-02-05 [1] CRAN (R 4.4.3)
    ##  shiny               1.11.1     2025-07-03 [1] CRAN (R 4.4.3)
    ##  shinythemes         1.2.0      2021-01-25 [1] CRAN (R 4.4.3)
    ##  snakecase           0.11.1     2023-08-27 [1] CRAN (R 4.4.3)
    ##  stringi             1.8.7      2025-03-27 [1] CRAN (R 4.4.3)
    ##  stringr           * 1.5.2      2025-09-08 [1] CRAN (R 4.4.3)
    ##  subcellularvis    * 0.0.0.9000 2025-10-31 [1] local
    ##  tibble              3.3.0      2025-06-08 [1] CRAN (R 4.4.3)
    ##  tidyr               1.3.1      2024-01-24 [1] CRAN (R 4.4.3)
    ##  tidyselect          1.2.1      2024-03-11 [1] CRAN (R 4.4.3)
    ##  timechange          0.3.0      2024-01-18 [1] CRAN (R 4.4.3)
    ##  UCSC.utils          1.2.0      2024-10-29 [1] Bioconduc~
    ##  UpSetR              1.4.0      2019-05-22 [1] CRAN (R 4.4.3)
    ##  vctrs               0.6.5      2023-12-01 [1] CRAN (R 4.4.3)
    ##  viridisLite         0.4.2      2023-05-02 [1] CRAN (R 4.4.3)
    ##  withr               3.0.2      2024-10-28 [1] CRAN (R 4.4.3)
    ##  xfun                0.53       2025-08-19 [1] CRAN (R 4.4.3)
    ##  xtable              1.8-4      2019-04-21 [1] CRAN (R 4.4.3)
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
