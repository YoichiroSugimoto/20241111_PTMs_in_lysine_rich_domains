
# Introduction
An R-based pipeline for the study:

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

# Required Package

Install the required [`ptm.stoichiometry`](https://github.com/YoichiroSugimoto/ptm.stoichiometry) package from GitHub:

```r
install.packages("devtools")
devtools::install_github("YoichiroSugimoto/ptm.stoichiometry")
```

All other packages are available from CRAN or Bioconductor.

---

# Directory Structure 
## R
Analyses were performed with Rmarkdown scripts stored in R directory.
- **p1_lysine-rich-domain-biology**
  - This script characterises lysine-rich domains of proteins
- **p2_lysine-rich-domain-LCMS**
  - The scripts take MaxQuant output as input and progress through PTM stoichiometry calculation, parameter optimisation, diagnostic ion analysis, and visualisation of site-specific lysine hydroxylations.

> Figures for the manuscript were generated using these scripts, and are saved in the corresponding `_files` folder created when knitting the `.Rmd` file.

## data

LC-MS/MS data must first be processed using `MaxQuant`. The pipeline uses the following `MaxQuant` output files:

- `evidence.txt` — peptide-level search results including retention times, intensities, and modifications
- `sample_info.csv` — sample metadata and experimental design information
- `Oxidation (K) DISites.txt` — lysine hydroxylation sites identified via diagnostic ions
- `Oxidised Propionylation (K) DISites.txt` — propionylated lysine hydroxylation sites identified via diagnostic ions

The output files can be retrieved from the PRoteomics IDEntification Database (PRIDE) using the following accession numbers:

| Dataset | Description | PRIDE Accession |
|---------|-------------|-----------------|
| A, B1, B2, C | PNAS 2022 | `PXD031221` |
| D | MS_KR_1 | `PXD079306` |
| E | MS_SS | `PXD031221` |


To run the pipeline, the `MaxQuant` output files and other raw data shuold be placed in the following locations:

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
│       │   ├── data-A_argC_m2_v7_def_evidence.txt
│       │   ├── data-A_trp_m2_v2_def_evidence.txt
│       │   ├── data-A_trp_m3_v3_def_evidence.txt
│       │   ├── data-A_trp_m4_v4_def_evidence.txt
│       │   ├── data-A_trp_m5_v5_AcP_evidence.txt
│       │   ├── data-A_trp_m5_v5_def_evidence.txt
│       │   ├── data-A_trp_m6_v6_def_evidence.txt
│       │   ├── data-A_trp_m7_v7_def_evidence.txt
│       │   ├── data-A_trp_m8_v8_def_evidence.txt
│       │   ├── data-B_trp_m2_v2_mCC_evidence.txt
│       │   ├── data-B_trp_m5_v5_mCC_evidence.txt
│       │   ├── data-B_trp_m7_v7_mCC_evidence.txt
│       │   ├── data-C_trp_m2_v2_mCC_evidence.txt
│       │   ├── data-C_trp_m5_v5_mCC_evidence.txt
│       │   ├── data-C_trp_m7_v7_mCC_evidence.txt
│       │   ├── data-C_trp_m8_v8_mCC_evidence.txt
│       │   ├── data-D_trp_m2_v2_def_evidence.txt
│       │   ├── data-D_trp_m5_v5_def_evidence.txt
│       │   └── data-D_trp_m7_v7_def_evidence.txt
│       ├── mqpar
│       │   ├── data-A_argC_m2_v7_def_mqpar.xml
│       │   ├── data-A_trp_m2_v2_def_mqpar.xml
│       │   ├── data-A_trp_m3_v3_def_mqpar.xml
│       │   ├── data-A_trp_m4_v4_def_mqpar.xml
│       │   ├── data-A_trp_m5_v5_def_mqpar.xml
│       │   ├── data-A_trp_m6_v6_def_mqpar.xml
│       │   ├── data-A_trp_m7_v7_def_mqpar.xml
│       │   ├── data-A_trp_m8_v8_def_mqpar.xml
│       │   ├── data-B_trp_m2_v2_mCC_mqpar.xml
│       │   ├── data-B_trp_m5_v5_mCC_mqpar.xml
│       │   ├── data-B_trp_m7_v7_mCC_mqpar.xml
│       │   ├── data-C_trp_m2_v2_mCC_mqpar.xml
│       │   ├── data-C_trp_m5_v5_mCC_mqpar.xml
│       │   ├── data-C_trp_m7_v7_mCC_mqpar.xml
│       │   ├── data-C_trp_m8_v8_mCC_mqpar.xml
│       │   ├── data-D_trp_m2_v2_def_mqpar.xml
│       │   ├── data-D_trp_m5_v5_def_mqpar.xml
│       │   └── data-D_trp_m7_v7_def_mqpar.xml
│       ├── runtime
│       │   ├── data-A_argC_m2_v7_def_runningTimes.txt
│       │   ├── data-A_trp_m2_v2_def_runningTimes.txt
│       │   ├── data-A_trp_m3_v3_def_runningTimes.txt
│       │   ├── data-A_trp_m4_v4_def_runningTimes.txt
│       │   ├── data-A_trp_m5_v5_def_runningTimes.txt
│       │   ├── data-A_trp_m6_v6_def_runningTimes.txt
│       │   ├── data-A_trp_m7_v7_def_runningTimes.txt
│       │   ├── data-A_trp_m8_v8_def_runningTimes.txt
│       │   ├── data-B_trp_m2_v2_mCC_runningTimes.txt
│       │   ├── data-B_trp_m5_v5_mCC_runningTimes.txt
│       │   ├── data-B_trp_m7_v7_mCC_runningTimes.txt
│       │   ├── data-C_trp_m2_v2_mCC_runningTimes.txt
│       │   ├── data-C_trp_m5_v5_mCC_runningTimes.txt
│       │   ├── data-C_trp_m7_v7_mCC_runningTimes.txt
│       │   ├── data-C_trp_m8_v8_mCC_runningTimes.txt
│       │   ├── data-D_trp_m2_v2_def_runningTimes.txt
│       │   ├── data-D_trp_m5_v5_def_runningTimes.txt
│       │   └── data-D_trp_m7_v7_def_runningTimes.txt
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





