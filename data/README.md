# Data directory

The pipeline uses the results of MaxQuant searches (with or without specifying diagnostic immonium ions) plus a few annotation files and public data.
**Only small files are committed to the repository** (PTM site lists, MaxQuant
parameters, run times, sample sheets, settings, small annotation tables). The
large mass-spectrometry inputs must be downloaded and placed at the paths
below.

- **`MQ_Std`** — MaxQuant search **without** diagnostic immonium ion consideration.
- **`MQ_DI`**  — MaxQuant search **with** diagnostic immonium ions.

The data folder structure is as follows:

```
data/
├── MQ_Std/                                    # standard MaxQuant search
│   └── PNAS2022/
│       └── data-A/   (and data-B/, data-C/, data-D/)        # one folder per source
│           ├── MS_dataset_overview_PXD031221_data-A.csv   [bundled]   (sample sheet)
│           └── trp_m7_v7_def/  (one folder per MaxQuant setting)
│               ├── data-A_trp_m7_v7_def_evidence.txt      [download]
│               ├── data-A_trp_m7_v7_def_mqpar.xml         [bundled]
│               └── data-A_trp_m7_v7_def_runningTimes.txt  [bundled]
├── MQ_DI/                                      # diagnostic-ion MaxQuant search
│   ├── PNAS2022/
│   │   └── data-A/   (and data-B1/, data-B2/, data-C/)
│   │       ├── MS_dataset_overview_PXD031221_data-A.csv   [bundled]
│   │       └── trp_m7_v7_def/
│   │           ├── data-A_trp_m7_v7_def_evidence.txt      [download]
│   │           ├── data-A_trp_m7_v7_def_mqpar.xml         [bundled]
│   │           └── ptm/
│   │               ├── Oxidation (K) DISites.txt             [bundled]
│   │               └── Oxidised Propionylation (K) DISites.txt  [bundled]
│   ├── MS_KR_1/                               # dataset D (single run)
│   │   ├── MS_KR_1_evidence.txt                           [download]
│   │   ├── ptm/  *DISites.txt                             [bundled]
│   │   └── sample_info.csv                                [bundled]
│   └── MS_SS/                                 # dataset E (single run)
│       ├── MS_SS_evidence.txt                             [download]
│       ├── ptm/  *DISites.txt                             [bundled]
│       └── sample_info.csv                                [bundled]
├── processed_data_from_PNAS2022/              # derived tables from the PNAS 2022 study
│   ├── all_protein_feature_per_position.csv              [download, ~1 GB]
│   └── long_K_stoichiometry_data.csv                     [download, ~16 MB]
├── XIC_MS_SS/
│   └── xic_MS_SS.csv                                     [bundled]
├── public_data/                               # public data
│   └── 20241113_cd-code.csv                              [bundled]
├── FP_diagnostic_ion_search/fragpipe_dataset-A/
│   ├── fragger.params / fragpipe.workflow                [bundled]
│   ├── dataset01.diagnosticIons.tsv                      
│   └── psm.tsv                                           [download, ~143 MB]
└── analysis_setting/                          [bundled]
    ├── PXD031221_sample_matrix.xlsx
    └── ptm_replacement_*.csv
```

## Notes

- The reference proteome FASTA is **not** stored under `data/`. Its path is configured as `reference_fasta` in `R/functions/_setup.R`.

