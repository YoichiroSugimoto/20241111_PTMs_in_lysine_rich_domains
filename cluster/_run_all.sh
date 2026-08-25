#!/bin/bash
#SBATCH --job-name=K_ox_run_all
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3-00:00:0
#SBATCH --mem=16G

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

## --- Run-time reporting -----------------------------------------------------
## The README quotes the wall-clock time and peak memory of a full pipeline run.
## Those figures come from the block below, which is printed to the SLURM output
## file. `report_runtime` is registered on EXIT so the summary is printed whether
## the job succeeds or fails (needed because of `set -e` above).

RUN_START=$(date +%s)

report_runtime() {
  local status=$?
  local elapsed=$(( $(date +%s) - RUN_START ))

  echo
  echo "==============================================================="
  echo " Pipeline run summary"
  echo "==============================================================="
  printf ' exit status      : %d\n' "$status"
  printf ' wall-clock       : %02d:%02d:%02d (hh:mm:ss)\n' \
    $(( elapsed / 3600 )) $(( (elapsed % 3600) / 60 )) $(( elapsed % 60 ))
  printf ' host             : %s\n' "$(hostname)"
  printf ' cpus-per-task    : %s\n' "${SLURM_CPUS_PER_TASK:-NA}"
  printf ' mem requested    : %s\n' "${SLURM_MEM_PER_NODE:-NA} MB"

  # Peak memory, read from the cgroup the job runs in. Unlike sacct's MaxRSS,
  # this IS available from inside the job. memory.peak is cgroup v2; the
  # max_usage_in_bytes path is the cgroup v1 equivalent.
  for f in /sys/fs/cgroup/memory.peak \
           /sys/fs/cgroup/memory/memory.max_usage_in_bytes; do
    if [[ -r "$f" ]]; then
      printf ' peak memory      : %s MiB (from %s)\n' \
        "$(( $(cat "$f") / 1048576 ))" "$f"
      break
    fi
  done

  # sacct CANNOT report MaxRSS or TotalCPU while the job is still RUNNING, and
  # this trap runs inside the job -- so those fields would come back blank.
  # Print the command to run once the job has finished instead.
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    echo
    echo " For peak memory and CPU time from the accounting database, run this"
    echo " AFTER the job has completed:"
    echo "   sacct -j $SLURM_JOB_ID --units=M \\"
    echo "     --format=JobID%-20,Elapsed,TotalCPU,MaxRSS,State"
  fi
  echo "==============================================================="

  return $status
}
trap report_runtime EXIT

# Pin Guix (cluster/channels.scm) + reproducible toolchain (cluster/manifest.scm),
# then knit. renv (activated in R/functions/_setup.R) restores the package library
# on first run; with the consistent GCC from manifest.scm, Rcpp et al. compile cleanly.
guix time-machine -C cluster/channels.scm -- \
  shell -m cluster/manifest.scm -- \
  Rscript R/_run_all.R
