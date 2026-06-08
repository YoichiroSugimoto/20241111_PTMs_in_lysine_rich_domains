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
  - This script characterises lysine-rich domains of proteins
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

## Data structure
```
.
├── 20241113_RBPbase_Hs_DescriptiveID.xlsx
├── 20241113_cd-code.csv
├── FP_diagnostic_ion_search
│   └── fragpipe_dataset-A
│       ├── dataset01.diagnosticIons.tsv
│       ├── fragger.params
│       ├── fragpipe.workflow
│       └── psm.tsv
├── MQ_standard
│   └── PNAS2022
│       ├── evidence/
│       ├── mqpar/
│       ├── runtime/
│       └── sample_info
│           ├── MS_dataset_overview_PXD031221_data-A.csv
│           ├── MS_dataset_overview_PXD031221_data-B.csv
│           ├── MS_dataset_overview_PXD031221_data-C.csv
│           └── MS_dataset_overview_PXD031221_data-D.csv
├── MQ_with_DI_no_waterloss
│   ├── MS_KR_1
│   │   ├── MS_KR_1_evidence.txt
│   │   ├── MS_KR_1_mqpar.xml
│   │   ├── MS_KR_1_proteinGroups.txt
│   │   ├── ptm
│   │   │   ├── Oxidation (K) DISites.txt
│   │   │   └── Oxidised Propionylation (K) DISites.txt
│   │   └── sample_info.csv
│   ├── MS_SS
│   │   ├── MS_SS_evidence.txt
│   │   ├── MS_SS_mqpar.xml
│   │   ├── MS_SS_proteinGroups.txt
│   │   ├── ptm
│   │   │   ├── Oxidation (K) DISites.txt
│   │   │   └── Oxidised Propionylation (K) DISites.txt
│   │   └── sample_info.csv
│   └── PNAS2022
│       ├── evidence
│       │   ├── data-A_trp_m7_v7_def_evidence.txt
│       │   ├── data-B1_trp_m7_v7_mCC_evidence.txt
│       │   ├── data-B2_trp_m7_v7_mCC_evidence.txt
│       │   └── data-C_trp_m7_v7_mCC_evidence.txt
│       ├── mqpar
│       │   ├── data-A_trp_m7_v7_def_mqpar.xml
│       │   ├── data-B1_trp_m7_v7_mCC_mqpar.xml
│       │   ├── data-B2_trp_m7_v7_mCC_mqpar.xml
│       │   └── data-C_trp_m7_v7_mCC_mqpar.xml
│       ├── proteingroups
│       │   ├── data-A_trp_m7_v7_def_proteinGroups.txt
│       │   ├── data-B1_trp_m7_v7_mCC_proteinGroups.txt
│       │   ├── data-B2_trp_m7_v7_mCC_proteinGroups.txt
│       │   └── data-C_trp_m7_v7_mCC_proteinGroups.txt
│       ├── ptm
│       │   ├── data-A_trp_m7_v7_def_ptm
│       │   │   ├── Oxidation (K) DISites.txt
│       │   │   └── Oxidised Propionylation (K) DISites.txt
│       │   ├── data-B1_trp_m7_v7_mCC_ptm
│       │   │   ├── Oxidation (K) DISites.txt
│       │   │   └── Oxidised Propionylation (K) DISites.txt
│       │   ├── data-B2_trp_m7_v7_mCC_ptm
│       │   │   ├── Oxidation (K) DISites.txt
│       │   │   └── Oxidised Propionylation (K) DISites.txt
│       │   └── data-C_trp_m7_v7_mCC_ptm
│       │       ├── Oxidation (K) DISites.txt
│       │       └── Oxidised Propionylation (K) DISites.txt
│       └── sample_info
│           ├── MS_dataset_overview_PXD031221_data-A.csv
│           ├── MS_dataset_overview_PXD031221_data-B1.csv
│           ├── MS_dataset_overview_PXD031221_data-B2.csv
│           └── MS_dataset_overview_PXD031221_data-C.csv
├── PNAS2022
│   ├── all_protein_feature_per_position.csv
│   └── long_K_stoichiometry_data.csv
├── analysis_setting
│   ├── PXD031221_sample_matrix.xlsx
│   ├── ptm_replacement_de-propionylate_for_hydroxylysine.csv
│   └── ptm_replacement_de-propionylate_for_hydroxylysine_fragpipe.csv
└── xic_MS_SS.csv
```

The rmd documents are stored in the R directory. The knitted documents are stored together, too.


