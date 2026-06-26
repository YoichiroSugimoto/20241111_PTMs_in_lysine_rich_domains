## _run_all.R — knit every analysis script in dependency order.
##
## Run from anywhere inside the repository, e.g.:
##   Rscript R/_run_all.R
##
## Prerequisites:
##   - the renv project library is restored (see README)
##   - input data is present under data/ (see README "Data Availability")
##
## Order reflects the pipeline dependencies:
##   p2-01 computes all stoichiometry tables that the later scripts read.

scripts <- c(
  ## Part 1 — biology of lysine-rich domains
  "p1_lysine-rich-domain-biology/p1-01_characterise_domains.Rmd",

  ## Part 2 — MS workflow: compute -> method evaluation -> validation -> biology
  "p2_lysine-rich-domain-LCMS/p2-01_compute_stoichiometry.Rmd",   # all stoichiometry tables
  "p2_lysine-rich-domain-LCMS/p2-02_optimise_search.Rmd",         # search optimisation
  "p2_lysine-rich-domain-LCMS/p2-03_compare_miscleavages.Rmd",    # search optimisation
  "p2_lysine-rich-domain-LCMS/p2-04_identify_diagnostic_ions.Rmd",# diagnostic-ion identification
  "p2_lysine-rich-domain-LCMS/p2-05_plot_diagnostic_ions.Rmd",    # diagnostic-ion evaluation
  "p2_lysine-rich-domain-LCMS/p2-06_analyse_MS_KR1.Rmd",          # validation (hypoxia vs normoxia)
  "p2_lysine-rich-domain-LCMS/p2-07_plot_stoichiometry.Rmd",      # integrated stoichiometry plots
  "p2_lysine-rich-domain-LCMS/p2-08_analyse_kinetics.Rmd",        # kinetics / O2 / dox biology
  "p2_lysine-rich-domain-LCMS/p2-09_export_tables.Rmd"            # supplementary data export
)

## Resolve the R/ directory that contains this script, regardless of caller wd.
.this_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
  error = function(e) file.path(getwd(), "R", "_run_all.R")
)
r.dir <- dirname(.this_file)

for (s in scripts) {
  rmd <- file.path(r.dir, s)
  message("\n=== Knitting: ", s, " ===")
  rmarkdown::render(rmd, envir = new.env(), quiet = TRUE)
}

message("\nAll scripts knitted.")
