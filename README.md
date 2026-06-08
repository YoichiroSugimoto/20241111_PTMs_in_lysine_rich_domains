# Introduction
An R-based pipeline for studying lysine hydroxylations from proteomics data. This pipeline uses the package "ptm.stoichiometry" to calculate the stoichiometry of lysine hydroxylations. 

---

The repository contains:

- **Biological characterisation of lysine-rich domains**
  - Analysis of the domain biology of proteins with lysine-rich regions
- **LC-MS data processing and stoichiometry calculation**
- **Visualisation of site-specific hydroxylation stoichiometry**

---
# Repository Structure 

```
.
├── R.Rproj
├── functions/
├── p1_lysine-rich-domain-biology/   # Characterises lysine-rich domains of proteins
├── p2_lysine-rich-domain-LCMS/
│   ├── p2-1_calculate_ptm_stoichiometry_new.rmd # Calculates stoichiometry from MaxQuant data 
│   ├── p2-2_optimise_parameters.rmd # Optimal settings for the database search
│   ├── p2-3-diagnostic-ions-lysine-hydorxylation-stoic.rmd # Identify diagnostic ions(DI)
│   ├── p2-4-diagnostic-ions-v2.rmd # Using DI to identify Hyl
│   ├── p2-4-diagnostic-ions.rmd # Using DI to identify Hyl
│   ├── p2-5_MS_KR1.Rmd # Calculates stoichiometry of dataset D
│   ├── p2-6_PSM_PTM_comparisons.Rmd # PSm coverage in different database search settings
│   ├── p2-7_MS_SS.Rmd # Calculates stoichiometry of dataset E
│   ├── p2-8_MS_SS_KR_plots_v2.Rmd # Visualises site-specific stoichiometry of lysine hydroxylations in varying O2% and dox induction
│   ├── p2-9_MS_SS_KR_plots_v3-1.Rmd
│   └── p2-10_MS_SS_KR_plots_v3-2.Rmd # Script for Supplementar data table 
├── renv/
│   └── settings.json
└── renv.lock
```


The rmd documents are stored in the R directory. The knitted documents are stored together, too.


