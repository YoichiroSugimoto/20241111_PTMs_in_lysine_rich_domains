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

## --- Knit, timing each script ------------------------------------------------
fmt_hms <- function(secs) {
  sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, round(secs %% 60))
}

timings <- data.frame(
  script  = character(0),
  seconds = numeric(0),
  status  = character(0),
  stringsAsFactors = FALSE
)

run_start <- proc.time()[["elapsed"]]

for (s in scripts) {
  rmd <- file.path(r.dir, s)
  message("\n=== Knitting: ", s, " ===")

  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({
    rmarkdown::render(rmd, envir = new.env(), quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("!!! FAILED: ", conditionMessage(e))
    FALSE
  })
  elapsed <- proc.time()[["elapsed"]] - t0

  message(sprintf("    -> %s in %s", if (ok) "done" else "FAILED", fmt_hms(elapsed)))

  timings <- rbind(timings, data.frame(
    script  = basename(s),
    seconds = elapsed,
    status  = if (ok) "ok" else "FAILED",
    stringsAsFactors = FALSE
  ))

  ## p2-01 produces the stoichiometry tables every later script reads, so there
  ## is no point continuing if it failed.
  if (!ok && grepl("p2-01", s)) stop("p2-01 failed; downstream scripts cannot run.")
}

total <- proc.time()[["elapsed"]] - run_start

## --- Summary table -----------------------------------------------------------
message("\n===============================================================")
message(" Per-script run times")
message("===============================================================")
for (i in seq_len(nrow(timings))) {
  message(sprintf(" %-42s %10s  %8s  %s",
                  timings$script[i],
                  fmt_hms(timings$seconds[i]),
                  sprintf("%.1f min", timings$seconds[i] / 60),
                  timings$status[i]))
}
message(sprintf("%s\n %-42s %10s  %8s",
                strrep("-", 63), "TOTAL", fmt_hms(total),
                sprintf("%.1f min", total / 60)))
message(" R ", getRversion(), " | ", R.version$platform)
message("===============================================================")

## Machine-readable copy, so the README figures can be regenerated.
results.dir <- file.path(dirname(r.dir), "results")
dir.create(results.dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(
  timings, file.path(results.dir, "run_all_timings.csv"), row.names = FALSE
)

if (any(timings$status == "FAILED")) {
  quit(status = 1L)
}
message("\nAll scripts knitted.")
