
# Introduction
An R-based pipeline for studying lysine hydroxylations from proteomics data, 
supporting the study:

 **Quantitative profiling of JMJD6-catalysed lysine hydroxylation reveals 
 residue-dependent oxygen sensitivity**

> *Preprint available on bioRxiv:* DOI: https://doi.org/10.64898/2026.06.08.730680

This pipeline uses the **`ptm.stoichiometry`** package to calculate the stoichiometry
of lysine hydroxylations from LC-MS/MS data processed by **`MaxQuant`**.

![R](https://img.shields.io/badge/R-4.5.1-blue)
![bioRxiv](https://img.shields.io/badge/bioRxiv-preprint-red)
![MaxQuant](https://img.shields.io/badge/MaxQuant-input-orange)

---

The repository contains:

- **Biological characterisation of lysine-rich domains**
- **LC-MS data processing and stoichiometry calculation**
  - Workflows for calculating PTM stoichiometry from MaxQuant output, optimising database search parameters, identifying diagnostic ions for hydroxylysine (Hyl), and computing site-specific stoichiometry across multiple datasets
- **Visualisation of site-specific hydroxylation stoichiometry**
  - Plots examining lysine hydroxylation across varying O₂ conditions and doxycycline induction levels

---

## Required Package Installation

Install the required [`ptm.stoichiometry`](https://github.com/YoichiroSugimoto/ptm.stoichiometry) package from GitHub:

```r
install.packages("devtools")
devtools::install_github("YoichiroSugimoto/ptm.stoichiometry")
```

---

# Directory Structure 
## R Markdown scripts
- **p1_lysine-rich-domain-biology**
  - This script characterises lysine-rich domains of proteins
- **p2_lysine-rich-domain-LCMS**
  - The scripts take MaxQuant output as input and progress through PTM stoichiometry calculation, parameter optimisation, diagnostic ion analysis, and visualisation of site-specific lysine hydroxylations.

> Figures for the manuscript were generated using the scripts below and are saved in the corresponding `_files` folder created when knitting the `.Rmd` file.

```
.
├── R.Rproj
├── functions
│   ├── 0-load_essential_packages.R
│   ├── 1-data_visualization_setting_and_functions.R
│   └── 2-useful_functions.R
├── p1_lysine-rich-domain-biology 
│   ├── p1-1_lysine-rich-domain-proteins.md
│   ├── p1-1_lysine-rich-domain-proteins.rmd
│   └── p1-1_lysine-rich-domain-proteins_files
│       └── ...
├── p2_lysine-rich-domain-LCMS
│   ├── p2-1_calculate_ptm_stoichiometry.md # Calculates stoichiometry of lysine hydroxylations (Hyl) from MaxQuant data 
│   ├── p2-1_calculate_ptm_stoichiometry.rmd
│   ├── p2-1_calculate_ptm_stoichiometry_files
│   │   └── ...
│   ├── p2-2_optimise_parameters.md # Optimal settings for the database search
│   ├── p2-2_optimise_parameters.rmd
│   ├── p2-2_optimise_parameters_files
│   │   └── ...
│   ├── p2-3_diagnostic_ions_hyl_identification.md # Identify diagnostic ions(DI)
│   ├── p2-3_diagnostic_ions_hyl_identification.rmd
│   ├── p2-3_diagnostic_ions_hyl_identification_files
│   │   └── ...
│   ├── p2-4_diagnostic_ions_plots.md # Using DI to identify Hyl
│   ├── p2-4_diagnostic_ions_plots.rmd
│   ├── p2-4_diagnostic_ions_plots_files
│   │   └── ...
│   ├── p2-5_MS_KR1.Rmd # Calculates stoichiometry of dataset D
│   ├── p2-5_MS_KR1.md
│   ├── p2-5_MS_KR1_files
│   │   └── ...
│   ├── p2-6_PSM_PTM_comparisons.Rmd # PSM coverage in different database search settings
│   ├── p2-6_PSM_PTM_comparisons.md
│   ├── p2-6_PSM_PTM_comparisons_files
│   │   └── ...
│   ├── p2-7_MS_SS.Rmd # Calculates stoichiometry of dataset E
│   ├── p2-7_MS_SS.md
│   ├── p2-8_ptm_stoichiometry_plots.Rmd # Visualises site-specific stoichiometry of Hyl in varying O2% and dox induction
│   ├── p2-8_ptm_stoichiometry_plots.md
│   ├── p2-8_ptm_stoichiometry_plots_files
│   │   └── ...
│   ├── p2-9_lysine_hydroxylation_O2_dox_visualisation.Rmd # Stoichiometry of Hyl in varying O2% and dox induction
│   ├── p2-9_lysine_hydroxylation_O2_dox_visualisation.md
│   ├── p2-9_lysine_hydroxylation_O2_dox_visualisation_files
│   │   └── ...
│   ├── p2-10_Lys_hydroxylation_BRD_stoic_export.Rmd # Script for Supplementary data table 
│   ├── p2-10_Lys_hydroxylation_BRD_stoic_export.md
├── renv
│   ├── activate.R
│   ├── library
│   │   └── ...
│   └── settings.json
└── renv.lock
```


## Data Availability

To run this pipeline, LC-MS/MS data must first be processed using `MaxQuant`. 
Data A–E can be retrieved from the PRoteomics IDEntification Database (PRIDE) 
using the following accession numbers:

| Dataset | Description | PRIDE Accession |
|---------|-------------|-----------------|
| A, B1, B2, C | PNAS 2022 | `PXD031221` |
| D | MS_KR_1 | `PXD079306` |
| E | MS_SS | `PXD031221` |

### Data Requirements
To calculate stoichiometry, the pipeline uses the following `MaxQuant` output files:

- `evidence.txt` — peptide-level search results including retention times, intensities, and modifications
- `sample_info.csv` — sample metadata and experimental design information
- `Oxidation (K) DISites.txt` — lysine hydroxylation sites identified via diagnostic ions
- `Oxidised Propionylation (K) DISites.txt` — propionylated lysine hydroxylation sites identified via diagnostic ions

### MaxQuant Parameters
The parameters of the MaxQuant search are available as mqpar.xml file, downloadable from the above-mentioned PRIDE Accession `PXD079306` and `PXD031221`.

### Data Structure

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
│       ├── evidence
│           ├── ...
│       ├── mqpar
│           ├── ...
│       ├── runtime
│           ├── ...
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
## Results 

The results directory contains the calculated stoichiometry files generated for each dataset.



