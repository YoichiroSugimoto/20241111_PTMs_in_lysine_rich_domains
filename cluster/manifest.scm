;; manifest.scm — build/runtime environment for knitting the R analysis pipeline.
;;
;; Use with:
;;   guix time-machine -C channels.scm -- shell -m manifest.scm -- Rscript R/_run_all.R
;;
;; Rcpp / compiled-package strategy
;; --------------------------------
;; The recurring Rcpp build failure ("_GLIBCXX20_DEPRECATED expected identifier",
;; "__is_nothrow_new_constructible was not declared") is caused by mixing TWO GCC
;; versions in one environment: previously `gcc-toolchain` and `gfortran-toolchain`
;; resolved to different GCC versions, so the g++ compiler and the libstdc++
;; headers did not match.
;;
;; The fix here: pin BOTH toolchains to a SINGLE version. gcc-toolchain provides
;; gcc/g++ (C/C++); gfortran-toolchain provides gfortran (Fortran) — gcc-toolchain
;; alone does NOT, so without it `gfortran` leaks in from the host /usr/bin (a
;; different GCC). The pipeline needs Fortran (mgcv, nlme, survival), and in this
;; channel gfortran-toolchain exists ONLY at 11.4.0 — while gcc-toolchain offers
;; many versions. The only version BOTH share is 11.4.0, so both are pinned @11.
;; (The earlier failure was gcc-toolchain default = 15.1.0 mixed with
;; gfortran-toolchain = 11.4.0 → two libstdc++ header trees → the Rcpp errors.)
;; gcc 11 fully supports C++17 (R 4.5's default) and compiles the whole lockfile.
;;
;; The ONLY hard rule: gcc-toolchain and gfortran-toolchain must be the SAME
;; version. To verify before a run (see cluster/README.md):
;;   guix shell --pure gcc-toolchain@11 gfortran-toolchain@11 -- \
;;     bash -c 'command -v g++ gfortran; g++ --version|head -1; gfortran --version|head -1'
;; Both must be /gnu/store paths reporting the same 11.x version.
;; (List available versions with: guix package -A '^(gcc|gfortran)-toolchain$')

(specifications->manifest
 (list
  ;; R + document rendering
  "r"
  "pandoc"

  ;; ONE consistent GCC toolchain (C / C++ / Fortran at the same version).
  ;; 11.4.0 is the only version both gcc-toolchain and gfortran-toolchain share.
  "gcc-toolchain@11"
  "gfortran-toolchain@11"
  "make"
  "cmake"        ; some packages (e.g. nloptr) build a bundled C library with cmake
  "pkg-config"

  ;; Core build utilities — R CMD INSTALL shells out to these during compilation.
  ;; Included so the environment is self-contained (works under `guix shell --pure`
  ;; and does not depend on the host /usr/bin).
  "bash" "coreutils" "sed" "grep" "gawk" "findutils" "tar" "gzip" "which"

  ;; Locales (fixes the "Setting LC_CTYPE failed, using C" warnings)
  "glibc-locales"

  ;; System libraries common CRAN packages link against (extend as needed)
  "libxml2" "openssl" "curl" "zlib"
  "fontconfig" "freetype" "harfbuzz" "fribidi"
  "cairo" "libpng" "libjpeg-turbo" "libtiff"
  "icu4c"        ; stringi (Unicode)
  "nlopt"        ; nloptr links this instead of building it via cmake
  "gmp" "mpfr")) ; common numeric-package deps (safe to include)
