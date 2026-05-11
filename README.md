# Introduction
**HylQuant** is a R-based pipeline for studying lysine hydroxylations from proteomics data. This pipeline uses the package "ptm.stoichiometry" to calculate the stoichiometry of lysine hydroxylations. 

---

**Abstract**

Hundreds of lysines in human proteins undergo hydroxylation in cells. Jumonji domain-containing protein 6 (JMJD6) hydroxylates multiple lysines within lysine-rich regions and is a major contributor to this modification across the proteome. As JMJD6 requires oxygen as a co-substrate, lysine hydroxylation is proposed to couple oxygen availability to cellular functions. However, the functional significance of this widespread modification remains incompletely understood, in part due to analytical challenges. Here, we systematically evaluated key steps in the analysis of mass spectrometry (MS) data from lysine-derivatised samples, and developed a workflow for comprehensive, accurate, and quantitative analysis of lysine hydroxylation in lysine-rich regions. Optimising database search strategies increases data coverage in lysine-rich regions. Leveraging immonium ion detection substantially improves confidence in hydroxylysine identification from MS data. We further demonstrated that stoichiometry derived from peptide MS intensity faithfully captures biologically relevant changes in lysine hydroxylation at amino acid resolution. We applied the workflow to bromodomain (BRD) proteins, epigenetic readers with lysine-rich regions that are extensively hydroxylated by JMJD6. Hydroxylation kinetics differ markedly between JMJD6 target lysines, with neighbouring sites displaying interdependence. Hypoxia suppresses hydroxylation in a site-dependent manner, with sites that exhibit slower apparent kinetics showing greater suppression. Overall, this work establishes a methodological and biological framework for understanding how widespread lysine hydroxylation links protein function to oxygen availability.

The repository contains:

-
-

---
# Repository Structure 

```
.
├── R.Rproj
├── functions
│   ├── 0-load_essential_packages.R
│   ├── 1-data_visualization_setting_and_functions.R
│   └── 2-useful_functions.R
├── p1_lysine-rich-domain-biology/
├── p2_lysine-rich-domain-LCMS/
│   ├── p2-1_calculate_ptm_stoichiometry_new.rmd
│   ├── p2-2_optimise_parameters.rmd
│   ├── p2-3-diagnostic-ions-lysine-hydorxylation-stoic.rmd
│   ├── p2-4-diagnostic-ions-v2.rmd
│   ├── p2-4-diagnostic-ions.rmd
│   ├── p2-5_MS_KR1.Rmd
│   ├── p2-6_PSM_PTM_comparisons.Rmd
│   ├── p2-7_MS_SS.Rmd
│   ├── p2-8_MS_SS_KR_plots_v2.Rmd
│   ├── p2-9_MS_SS_KR_plots_v3-1.Rmd
│   ├── p2-10_MS_SS_KR_plots_v3-2.Rmd
│   └── p2-11_MS_SS_KR_plots_v3-3.Rmd
├── renv/
│   └── settings.json
└── renv.lock
```


The rmd documents are stored in the R directory. The knitted documents are stored together, too.


