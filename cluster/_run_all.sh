#!/bin/bash
#SBATCH --job-name=K_ox_run_all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3-00:00:0
#SBATCH --mem=8G

## e.g. sbatch cluster/_run_all.sh

source /home/ysugimo/.bashrc

# Avoid conda(base) libraries clashing with Guix's
conda deactivate 2>/dev/null || true

# Fail the job if any command below errors (set after sourcing .bashrc/conda,
# which may reference unset vars under `set -u`).
set -eo pipefail

# Locale (fixes the "Setting LC_CTYPE failed, using C" warnings).
# glibc-locales is in manifest.scm, so GUIX_LOCPATH is set inside the shell.
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# Run from the project root (it must contain cluster/ and R/).
# Submit as:  sbatch cluster/_run_all.sh
cd /fast/AG_Sugimoto/home/users/yoichiro/projects/20241111_PTMs_in_lysine_rich_domains

RUN_START=$(date +%s)

report_runtime() {
  local status=$?
  local elapsed=$(( $(date +%s) - RUN_START ))
  echo
  echo " Pipeline run summary"
  printf ' exit status : %d\n' "$status"
  printf ' wall-clock  : %02d:%02d:%02d\n' \
    $(( elapsed / 3600 )) $(( (elapsed % 3600) / 60 )) $(( elapsed % 60 ))
  if [[ -n "${SLURM_JOB_ID:-}" ]] && command -v sacct >/dev/null 2>&1; then
    sacct -j "$SLURM_JOB_ID" --units=G \
      --format=JobID%-20,Elapsed,TotalCPU,MaxRSS,MaxVMSize,State 2>/dev/null
  fi
  return $status
}
trap report_runtime EXIT

# Pin Guix (cluster/channels.scm) + reproducible toolchain (cluster/manifest.scm),
# then knit. renv (activated in R/functions/_setup.R) restores the package library
# on first run; with the consistent GCC from manifest.scm, Rcpp et al. compile cleanly.
guix time-machine -C cluster/channels.scm -- \
  shell -m cluster/manifest.scm -- \
  Rscript R/_run_all.R
