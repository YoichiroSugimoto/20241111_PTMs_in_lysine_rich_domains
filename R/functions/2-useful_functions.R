library("Biostrings")

## Create a directory (and any missing parents) if it does not already exist.
## Idempotent and silent: no-op / no warning when the directory is already there.
create.dir <- function(dir.name){
    dir.create(dir.name, recursive = TRUE, showWarnings = FALSE)
}

create.dirs <- function(dirs){
    for(dir.name in dirs){
        create.dir(dir.name)
    }
}

system.cat <- function(cmd){
    stdout.text <- system(cmd, intern = TRUE)
    return(stdout.text)
}

setMethod("getSeq", "XStringSet",
    function(x, names)
    {
        stopifnot(is.character(names) || is(names, "GRanges") ||
                  is(names, "GRangesList"),
                  !is.null(names(x)))
        if (is.character(names)) {
            found <- names %in% names(x)
            regexNames <- unlist(lapply(names[!found], grep, names(x),
                                        value=TRUE))
            names <- c(names[found], regexNames)
            return(x[names])
        } else if (is(names, "GRangesList")) {
            gr <- unlist(names, use.names=FALSE)
        } else {
            gr <- names
        }
        ignoringStrand <- any(strand(gr) != "*") &&
            !hasMethod(reverseComplement, class(x))
        if (ignoringStrand) {
            warning("some strand(x) != '*' but ",
                    "strand has no meaning for ", class(x))
        }
        rl <- as(gr, "IntegerRangesList")
        ans <- unsplit(extractAt(x[names(rl)], unname(rl)), seqnames(gr))
        if (!ignoringStrand) {
            minus <- strand(gr) == "-"
            ans[minus] <- reverseComplement(ans[minus])
        }
        if (is(names, "GRangesList")) {
            ans <- relist(ans, names)
        }
        ans
    }
)


## Read a PTM-stoichiometry table written by ptm.stoichiometry.
## Reconstructs the path as <dir_path>/<pre_prefix><prefix>PTM_stoichiometry<post_fix>.csv
## (pre_prefix/post_fix default to "" so MQ_Std tables read with prefix + dir_path alone).
read_stoic_data <- function(prefix, pre_prefix = "", post_fix = "", dir_path) {
    dt <- fread(file.path(
        dir_path, paste0(pre_prefix, prefix, "PTM_stoichiometry", post_fix, ".csv")
    ))
    dt[, condition := gsub("_$", "", prefix)] # drop trailing "_" from prefix
    return(dt)
}


## Load and preprocess the MS_KR_1, MS_SS and PNAS stoichiometry tables that
## p2-07 and p2-08 share. Returns a named list of data.tables. Each diagnostic-ion
## table gains a logical `is_diagnostic_peak`; the curated PNAS table is column-renamed.
load_stoichiometry_datasets <- function(results.dir, data.dir) {
    read_di <- function(prefix, subdir) read_stoic_data(
        prefix = prefix, pre_prefix = "", post_fix = "_DI",
        dir_path = file.path(results.dir, "p2-analysis-setting", subdir)
    )

    MS_KR1_stoic_dt <- read_di("MS_KR_1_", "MS_KR_1")
    MS_KR1_stoic_dt[, is_diagnostic_peak := diagnostic_peak == "+"]

    MS_SS_stoic_dt <- read_di("MS_SS_", "MS_SS")
    MS_SS_stoic_dt[, is_diagnostic_peak := diagnostic_peak == "+"]

    pnas2022_stoic_dt <- read_di("DI_data-A_trp_m7_v7_def_", "MQ_DI")
    pnas2022_stoic_dt[, is_diagnostic_peak := diagnostic_peak == "+"]

    pnas2022_dt <- fread(file.path(data.dir, "processed_data_from_PNAS2022/long_K_stoichiometry_data.csv"))
    setnames(
        pnas2022_dt,
        old = c("uniprot_id", "position", "residue", "oxK_ratio", "data_source",
                "total_n_feature_oxK", "total_n_feature_K"),
        new = c("protein_accession", "aa_pos", "aa", "stoichiometry", "sample_name",
                "sum_psm_mapped", "sum_psm_mapped_per_position")
    )

    list(MS_KR1_stoic_dt = MS_KR1_stoic_dt, MS_SS_stoic_dt = MS_SS_stoic_dt,
         pnas2022_stoic_dt = pnas2022_stoic_dt, pnas2022_dt = pnas2022_dt)
}


## Zero-fill non-hydroxylated lysines and reshape hydroxylation stoichiometry to
## one column per group. Shared by p2-06 and p2-08, which differ only in which
## column defines the groups (`group_col`, used for the wide output) and which
## defines a sample for the per-sample zero-fill (`sample_col`, defaults to
## `group_col`).
##
## For every position with >=1 oxidation event anywhere, samples lacking an
## oxidation row get an explicit stoichiometry = 0 row, so the wide output holds
## a value (not NA) wherever the site was observed. Positions are kept only with
## sum_psm_mapped_per_position > 2, and the most-hydroxylated row per
## protein/gene/position/group is retained before casting wide.
##
## NOTE: modifies `all_stoic_dt` by reference (adds pos_id / sample_pos_id), as
## the original per-script functions did.
contrast_hydroxylation <- function(all_stoic_dt, group_col, sample_col = group_col) {
    dt <- all_stoic_dt

    # Position and per-sample-position identifiers used for the zero-fill
    dt[, `:=`(
        pos_id        = paste0(protein_accession, "_", aa_pos),
        sample_pos_id = paste0(get(sample_col), "_", protein_accession, "_", aa_pos)
    )]

    # Positions / sample-positions that carry at least one oxidation event
    oxidation_ids <- dt[
        grepl("[Oxidation (K)]", ptm, fixed = TRUE), unique(pos_id)
    ]
    oxidation_sample_ids <- dt[
        grepl("[Oxidation (K)]", ptm, fixed = TRUE), unique(sample_pos_id)
    ]

    # Non-hydroxylated K rows for hydroxylated positions that lack an oxidation
    # row in a given sample -> these become explicit stoichiometry = 0 entries
    no_hydroxyK_dt <- copy(dt[aa == "K"])
    no_hydroxyK_dt <- no_hydroxyK_dt[pos_id %in% oxidation_ids] %>%
        {.[!sample_pos_id %in% oxidation_sample_ids]}

    no_hydroxyK_dt[, `:=`(
        sum_psm_mapped = 0,
        stoichiometry  = 0,
        ptm            = "[Oxidation (K)]"
    )]

    # Combine measured oxidation rows with the zero-filled rows
    hydroxyK_dt <- rbind(
        dt[grepl("[Oxidation (K)]", ptm, fixed = TRUE)],
        no_hydroxyK_dt
    )

    # Drop sparsely covered positions (>2 PSMs required for a reliable estimate)
    hydroxyK_dt <- hydroxyK_dt[sum_psm_mapped_per_position > 2]

    # Keep the most-hydroxylated row per protein/gene/position/group
    hydroxyK_dt <- hydroxyK_dt[order(stoichiometry, decreasing = TRUE)][
        !duplicated(paste0(protein_accession, gene_name, aa_pos, get(group_col)))
    ]

    # Reshape long -> wide: one stoichiometry column per group
    dcast(
        hydroxyK_dt,
        stats::as.formula(paste("protein_accession + gene_name + aa_pos ~", group_col)),
        value.var = "stoichiometry"
    )
}
