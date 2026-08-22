# Running the pipeline on the cluster (Guix + SLURM)

This documents how to knit the analysis on the MDC cluster, which uses **GNU Guix**
(not environment modules). The goal is a reproducible environment: a pinned Guix
revision (`channels.scm`) + a fixed package set (`manifest.scm`), with renv managing
the R packages.

These cluster files live in `cluster/`. **Run all commands from the project root**
(the directory that contains `cluster/` and `R/`), so the relative paths below resolve.

renv is activated automatically at R startup by the project-root `.Rprofile`
(`source("renv/activate.R")`). This must happen *before* `rmarkdown`/`knitr` load, so
do **not** run `renv::restore()` inside a knit — restore is a separate, clean-session
step (below). `R/functions/_setup.R` only re-activates renv as a fallback; it does
not restore.

## Files

| File | Purpose |
|------|---------|
| `cluster/channels.scm` | Pins the Guix revision (toolchain, R version, build deps). |
| `cluster/manifest.scm` | Declares R + pandoc + **one consistent GCC toolchain** + system libs. |
| `cluster/_run_all.sh`  | SLURM batch script: wraps `Rscript R/_run_all.R` in the Guix env. |
| `R/renv.lock`          | The R package lockfile (note: kept under `R/`, not the root). |

## Why a single GCC version (the Rcpp fix)

Compiling Rcpp (and other C++ packages) fails with errors like
`_GLIBCXX20_DEPRECATED expected identifier` / `__is_nothrow_new_constructible was
not declared` when **two different GCC versions** are mixed in one environment —
the g++ compiler and the libstdc++ headers no longer match.

Causes seen here:
- `gcc-toolchain` and `gfortran-toolchain` resolving to **different** GCC versions, and/or
- `gcc-toolchain` alone (it does **not** include gfortran) letting `gfortran` leak
  in from the host `/usr/bin` — a different GCC again.

`manifest.scm` fixes this by pinning **both** `gcc-toolchain` and `gfortran-toolchain`
to the **same** version (`@11`). The pipeline needs Fortran (`mgcv`, `nlme`, `survival`),
and in this channel `gfortran-toolchain` exists **only at 11.4.0** — so 11.4.0 is the
only version both toolchains share. gcc 11 fully supports C++17 (R 4.5's default) and
compiles the whole lockfile. With one consistent compiler + standard library, renv
compiles everything from source normally — no Rcpp special-casing.

## One-time setup (login node, has internet)

```bash
cd <project root>

# Keep conda from leaking incompatible libs into the build
conda deactivate 2>/dev/null || true

# 1. Pin the Guix revision you validate against (overwrites the template):
guix describe -f channels > cluster/channels.scm

# 2. Verify the toolchain: BOTH must be /gnu/store paths and the SAME 11.x version.
#    (bash + coreutils are added so the --pure shell can run bash/head.)
guix shell --pure bash coreutils gcc-toolchain@11 gfortran-toolchain@11 -- \
  bash -c 'command -v g++; command -v gfortran; g++ --version | head -1; gfortran --version | head -1'
#   - If a path is /usr/bin/... or versions differ -> fix cluster/manifest.scm first.
#   - List available toolchain versions: guix package -A '^(gcc|gfortran)-toolchain$'
#     (gfortran-toolchain exists only at 11.4.0 here, so both are pinned to @11.)

# 3. Build the renv library once (compiles with the consistent toolchain).
#    The lockfile records every package the pipeline uses, including the
#    Bioconductor ones, so a plain restore is enough. Only the two packages
#    recorded as local sources have to be installed separately.
#    renv is auto-activated by the root .Rprofile, so no source("renv/activate.R").
#    Note: eulerr is pinned to 7.0.2 in the lockfile because v8 requires
#    Rust >= 1.91 (only rustc 1.88 is available here).
guix time-machine -C cluster/channels.scm -- shell -m cluster/manifest.scm -- Rscript -e '
  renv::restore(lockfile = "R/renv.lock", prompt = FALSE)
  renv::install("/path/to/ptm.stoichiometry")   # local copy
  renv::install("/path/to/subcellularvis")      # local copy
'

# 4. (optional) Confirm the library loads everything before submitting:
guix time-machine -C cluster/channels.scm -- shell -m cluster/manifest.scm -- Rscript -e '
  for (p in c("Biostrings","org.Hs.eg.db","ptm.stoichiometry","subcellularvis",
              "eulerr","ggpubr","openxlsx","khroma","knitr","rmarkdown"))
    suppressMessages(library(p, character.only = TRUE)); cat("all packages load OK\n")'
```

Notes:
- Step 3 needs network access (CRAN/Bioconductor/GitHub) and warms the renv cache +
  Guix store, so compute nodes don't download/compile anything.
- The renv library lives on the shared `/fast` filesystem, so the batch job only
  *activates and loads* it (it does not restore) — see `R/functions/_setup.R`.
- Missing packages (not in the original lockfile) and their sources:
  `eulerr`, `ggpubr`, `openxlsx` (CRAN); `Biostrings`, `org.Hs.eg.db` (Bioconductor);
  `ptm.stoichiometry` (local: `…/projects/ptm.stoichiometry`); `subcellularvis`
  (GitHub: `JoWatson2011/subcellularvis`). `mgcv` is a recommended package shipped
  with R, so it needs no entry.
- After step 3, **commit the updated `R/renv.lock`** so the additions are permanent.
- The lockfile pins R 4.5.1 while Guix provides R 4.5.0 — only a harmless
  version-mismatch *warning* (patch releases are ABI-compatible).

## Submit the job

From the project root:

```bash
sbatch cluster/_run_all.sh
```

`cluster/_run_all.sh` activates the environment with:

```bash
guix time-machine -C cluster/channels.scm -- shell -m cluster/manifest.scm -- Rscript R/_run_all.R
```

and sets the locale (`LC_ALL=en_US.UTF-8`, with `glibc-locales` in the manifest) to
silence the `Setting LC_CTYPE failed, using "C"` warnings.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `_GLIBCXX...` / `__is_nothrow...` C++ errors | GCC mismatch — ensure `gcc-toolchain` and `gfortran-toolchain` are the same version (verify step 2). |
| `gfortran: command not found` or `/usr/bin/gfortran` used | `gfortran-toolchain` missing from `manifest.scm`, or `--pure` not used so host PATH leaked. |
| `there is no package called '<pkg>'` (e.g. khroma) | renv library not (fully) restored — run one-time setup step 3 to completion on a login node. |
| `Package 'knitr' ... cannot be unloaded: namespace ... imported by 'rmarkdown'` | renv activated mid-knit instead of at startup. Ensure the root `.Rprofile` exists, and never call `renv::restore()` during a render. |
| `cannot create file '.../results/...'` | output dir missing — created by the scripts now; otherwise `mkdir -p results`. |
| Compile fails for a missing system lib | add it to `manifest.scm` (e.g. `gmp`, `mpfr`, `libgit2`, `unixodbc`) and re-run step 3. |
| Locale warnings | ensure `glibc-locales` is in `manifest.scm` and `LC_ALL`/`LANG` are exported. |

## Optional: stricter isolation

For maximum hygiene (prevents host `/usr/bin` and conda from leaking), add `--pure`
to the shell, preserving only the locale vars:

```bash
guix time-machine -C cluster/channels.scm -- \
  shell --pure --preserve='^(LC_|LANG$|GUIX_LOCPATH$)' -m cluster/manifest.scm -- \
  Rscript R/_run_all.R
```

Note `--pure` also drops proxy variables — if the login node needs an HTTP proxy
to reach CRAN during the one-time restore, keep that step non-`--pure` (or also
preserve `^(http_proxy|https_proxy)$`).
