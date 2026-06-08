# Introduction
An R-based pipeline for studying lysine hydroxylations from proteomics data. This pipeline uses the package "ptm.stoichiometry" to calculate the stoichiometry of lysine hydroxylations. 

---

The repository contains:

- **Biological characterisation of lysine-rich domains**
  - Analysis of the domain biology of proteins with lysine-rich regions
- **LC-MS data processing and stoichiometry calculation**
  - Workflows for calculating PTM stoichiometry from MaxQuant output, optimising database search parameters, identifying diagnostic ions for hydroxylysine (Hyl), and computing site-specific stoichiometry across multiple datasets
- **Visualisation of site-specific hydroxylation stoichiometry**
  - Plots examining lysine hydroxylation across varying O₂ conditions and doxycycline induction levels

---
# Directory Structure 
## R(markdown) scripts
- p1_lysine-rich-domain-biology
  - This characterises lysine-rich domains of proteins
- p2_lysine-rich-domain-LCMS
  - The scripts take MaxQuant output as input and progress through PTM stoichiometry calculation, parameter optimisation, diagnostic ion analysis, and visualisation of site-specific lysine hydroxylations.

```
.
├── R.Rproj
├── functions/
├── p1_lysine-rich-domain-biology/   # Characterises lysine-rich domains of proteins
├── p2_lysine-rich-domain-LCMS/
│   ├── p2-1_calculate_ptm_stoichiometry_new.rmd # Calculates stoichiometry from MaxQuant data 
│   ├── p2-2_optimise_parameters.rmd # Optimal settings for the database search
│   ├── p2-3_diagnostic_ions_hyl_identification.rmd # Identify diagnostic ions(DI)
│   ├── p2-4_diagnostic_ions_plots.rmd # Using DI to identify Hyl
│   ├── p2-5_MS_KR1.Rmd # Calculates stoichiometry of dataset D
│   ├── p2-6_PSM_PTM_comparisons.Rmd # PSm coverage in different database search settings
│   ├── p2-7_MS_SS.Rmd # Calculates stoichiometry of dataset E
│   ├── p2-8_ptm_stoichiometry_plots.Rmd # Visualises site-specific stoichiometry of lysine hydroxylations in varying O2% and dox induction
│   ├── p2-9_lysine_hydroxylation_O2_dox_visualisation.Rmd # Visualises site-specific stoichiometry of lysine hydroxylations in varying O2% and dox induction
│   └── p2-10_Lys_hydroxylation_BRD_stoic_export.Rmd # Script for Supplementary data table 
├── renv/
│   └── settings.json
└── renv.lock
```


The rmd documents are stored in the R directory. The knitted documents are stored together, too.


