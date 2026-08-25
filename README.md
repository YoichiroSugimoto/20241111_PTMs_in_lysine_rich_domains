
# Introduction
An R-based pipeline for studying lysine hydroxylations from proteomics data, 
supporting the study:

 **Quantitative profiling of JMJD6-catalysed lysine hydroxylation reveals a graded, residue-dependent readout of oxygen availability**

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

## Installation and setup

### 1. R packages
R packages are managed with [`renv`](https://rstudio.github.io/renv/). Restore the
project library from the lockfile:

```r
renv::restore(lockfile = "R/renv.lock", prompt = FALSE)
```

[`ptm.stoichiometry`](https://github.com/YoichiroSugimoto/ptm.stoichiometry) and [`subcellularvis`](https://github.com/JoWatson2011/subcellularvis) are recorded in the lockfile as local
sources, so they must be installed from a local copy:

```r
renv::install("/path/to/ptm.stoichiometry")
renv::install("/path/to/subcellularvis")
```

### 2. Path
Paths are resolved automatically from the repository root (marked by the `.here`
file).
The reference proteome `UP000005640_9606.fasta` data is not included in this repository due to its size. Thus, the data should be downloaded from UniProt, and the path should be specified in [`R/functions/_setup.R`](R/functions/_setup.R).

### 3. Input data
The MaxQuant inputs should be placed in `data/`, as documented in [`data/README.md`](data/README.md).


---

# Pipeline

Scripts are stored in `R/`, numbered in run order. Knit them in sequence, or run
`Rscript R/_run_all.R`. 
A full run takes approximately 20 minutes on a 4-CPU, 16 GB node.

| Script | Purpose |
|--------|---------|
| `p1-01_characterise_domains` | Proteome-wide characterisation of lysine-rich domains |
| `p2-01_compute_stoichiometry` | Compute PTM stoichiometry for **all** datasets (the engine) |
| `p2-02_optimise_search` | Optimise MaxQuant database-search settings |
| `p2-03_compare_miscleavages` | Compare PSM coverage across miscleavage settings |
| `p2-04_identify_diagnostic_ions` | Identify immonium / diagnostic ions marking hydroxylysine |
| `p2-05_plot_diagnostic_ions` | Evaluate diagnostic-ion gains (precision, WT/KO overlap) |
| `p2-06_analyse_MS_KR1` | Validate stoichiometry against hypoxia (MS_KR_1) |
| `p2-07_plot_stoichiometry` | Site-specific stoichiometry plots (BRD2/3/4) |
| `p2-08_analyse_kinetics` | Hydroxylation kinetics, O2 / dox sensitivity, site interaction |
| `p2-09_export_tables` | Export stoichiometry supplementary tables (Excel) |

### Repository layout

```
R/
├── R.Rproj
├── _run_all.R            # knit all scripts in run order
├── functions/
│   ├── _setup.R          # shared setup: paths, packages, helpers (sourced by every script)
│   ├── 0-load_essential_packages.R
│   ├── 1-data_visualization_setting_and_functions.R
│   └── 2-useful_functions.R
├── p1_lysine-rich-domain-biology/
│   └── p1-01_characterise_domains.Rmd
└── p2_lysine-rich-domain-LCMS/
    ├── p2-01_compute_stoichiometry.Rmd
    ├── p2-02_optimise_search.Rmd
    ├── p2-03_compare_miscleavages.Rmd
    ├── p2-04_identify_diagnostic_ions.Rmd
    ├── p2-05_plot_diagnostic_ions.Rmd
    ├── p2-06_analyse_MS_KR1.Rmd
    ├── p2-07_plot_stoichiometry.Rmd
    ├── p2-08_analyse_kinetics.Rmd
    └── p2-09_export_tables.Rmd
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

### Data structure

The full data-directory layout — which files are **bundled** in this repository
(PTM site lists, MaxQuant parameters, run times, sample sheets, settings) and
which must be **downloaded**  — is
documented in [`data/README.md`](data/README.md).
