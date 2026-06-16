## _setup.R — shared environment setup, sourced by every analysis .Rmd.
##
## Replaces the per-script "Environment setup" boilerplate. It resolves the
## project root portably, loads the helper functions, packages, ggplot/knitr
## settings, and defines the standard paths.
##
## NOTE: this does NOT manage renv. It assumes the project library has already
## been restored (see README). It only loads packages and sets up paths.

## --- Resolve the project root (portable, no extra dependencies) -------------
## Walks up from the working directory to the repo root, identified by the
## `.here` sentinel file shipped at the top level. This avoids the hard-coded
## absolute paths that previously appeared in every script.
.find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (!file.exists(file.path(path, ".here"))) {
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate project root: no '.here' file found above ", start)
    }
    path <- parent
  }
  path
}

project.dir <- .find_project_root()
data.dir    <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")

## --- External UniProt reference proteome ------------------------------------
## The reference FASTA is NOT stored in the repository (it is large and freely
## available from UniProt: proteome UP000005640, Homo sapiens).
##
## >>> EDIT THIS LINE <<< — set it to your local copy of the FASTA.
## This is the only path that needs to be configured per machine.
reference_fasta <- "/fast/AG_Sugimoto/reference/uniprot/human/UP000005640_9606.fasta"

## --- Shared helper functions and settings -----------------------------------
## Sources, in order:
##   functions/0-load_essential_packages.R           (core data/plot packages)
##   functions/1-data_visualization_setting_*.R       (ggplot theme + knitr opts)
##   functions/2-useful_functions.R                   (helpers + Biostrings)
## The "^[0-9]" pattern deliberately excludes this _setup.R file itself.
invisible(lapply(
  list.files(
    file.path(project.dir, "R", "functions"),
    pattern = "^[0-9].*\\.R$", full.names = TRUE
  ),
  source
))

## --- Analysis packages used across most scripts -----------------------------
## Script-specific packages (e.g. janitor, openxlsx, ggpubr) are loaded in the
## individual scripts that need them.
library("ptm.stoichiometry")
library("readxl")


## --- Reference proteome (loaded once; used across scripts) -------------------
ref_protein_dt <- import_reference_fasta(reference_fasta)
all.protein.bs <- Biostrings::readAAStringSet(reference_fasta)
