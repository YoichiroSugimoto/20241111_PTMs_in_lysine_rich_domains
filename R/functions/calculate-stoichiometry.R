library("stringr")                                                                                                                                                                                                                               
library("data.table")                                                                                                                                                                                                                            
library("dplyr")                                                                                                                                                                                                                                 
library("ggplot2")                                                                                                                                                                                                                               
library("janitor")                                                                                                                                                                                                                               
library("scales")                                                                                                                                                                                                                                
library("knitr")                                                                                                                                                                                                                                 
library("Biostrings")

## Referece protein information
all.protein.bs <- readAAStringSet(
    file.path(
        "/fast/AG_Sugimoto/reference/uniprot/human",
        "UP000005640_9606.fasta")
)

ref.protein.dt <- data.table(
    protein_accession = str_split_fixed(names(all.protein.bs), "\\|", n = 3)[, 2],
    gene_name = str_extract(names(all.protein.bs), "(?<=GN=)\\S+"),
    protein_seq = as.character(all.protein.bs)
)



## read MS data
## Useful function
replace_outer_parentheses <- function(str){
    n <- nchar(str)
    brackets <- integer(n)
    cumsum <- 0
    for (i in seq_len(n)) {
        if (substr(str, i, i) == "(") {
            cumsum <- cumsum + 1
            if (cumsum == 1) brackets[i] <- 1
        } else if (substr(str, i, i) == ")") {
            if (cumsum == 1) brackets[i] <- 1
            cumsum <- cumsum - 1
        }
    }
    str <- strsplit(str, "")[[1]]
    str[brackets == 1] <- ifelse(str[brackets == 1] == "(", "[", "]")
    paste(str, collapse = "")
}



## Read MaxQuant Data
readMaxQuantData <- function(mq.ev.input, sample.info.file, mod.dt, output.prefix){

    mq.ev.df <- lapply(
        mq.ev.input,
        FUN = fread,
        header=TRUE,
        quote=""
    ) %>%
        rbindlist(fill = TRUE)

    mq.ev.df <- clean_names(mq.ev.df)

    mq.ev.df <-   mq.ev.df[
        potential_contaminant != "+" &
        reverse != "+"
    ]

    ## Data for sample QC
    qc.dt <- mq.ev.df[, list(
        total_peptide_count = .N, 
        total_MSMS_count = sum(ms_ms_count),
        total_intensity = sum(intensity, na.rm = TRUE)
    ),
    by = raw_file][order(raw_file)]

    fwrite(qc.dt, file = paste0(output.prefix, "QC.csv"))  
    
    selected.columns <- c(
        "file_name", "protein_accession", "base_sequence", "full_sequence", "peak_intensity", "PSMs_mapped", "retention_time"
    )

    setnames(
        mq.ev.df,
        old = c("raw_file", "leading_razor_protein", "sequence", "modified_sequence", "intensity", "ms_ms_count", "calibrated_retention_time"),
        new = selected.columns
    )

    mq.ev.df <- mq.ev.df[, selected.columns, with = FALSE]

    mq.ev.df[, `:=`(
        protein_accession =
            str_split_fixed(protein_accession, ";", n = 2)[, 1],
        full_sequence = gsub("_", "", full_sequence)
    )]

    mq.ev.df[
      , full_sequence := replace_outer_parentheses(full_sequence),
        by = seq_len(nrow(mq.ev.df))
    ]

    mod.vec <- setNames(mod.dt[, mod_name], nm = mod.dt[, mq_mod_name])
    
    mq.ev.df[, full_sequence := str_replace_all(full_sequence, mod.vec)]
    
    
    samples.dt <- fread(sample.info.file)
    samples.dt <- samples.dt[, .(file_name, sample_name)] 
        
    mq.ev.df <- merge(
        mq.ev.df, samples.dt, by = "file_name"
    )
    
    return(mq.ev.df)
}


preprocessPeakData <- function(qp.dt, ref.protein.dt, output.prefix){
    
    qp.ref.dt <- ref.protein.dt |>  
        merge(qp.dt, by = c("protein_accession")) |>
        mutate(the_number_of_match = str_count(protein_seq, base_sequence)) |>
        filter(the_number_of_match == 1,
               !peak_intensity == 0, !is.null(peak_intensity))

    ## add columns: peptide start and end positions
    qp.ref.dt <- qp.ref.dt |>
        mutate(
            peptide_start = str_locate(protein_seq, base_sequence)[,1],
            peptide_end = str_locate(protein_seq, base_sequence)[,2]
        )

    qp.ref.dt <- qp.ref.dt %>% mutate(protein_seq = NULL)
    ## Filtering & string replacements
    ## terminal modifications
    qp.ref.dt <- qp.ref.dt |>
        filter(
            !grepl("\\|", full_sequence)
        ) |>
        mutate(
            full_sequence = str_replace(full_sequence, "^\\[.*?\\]", "")
        )

    fwrite(qp.ref.dt, file = paste0(output.prefix, "all_processed_peak_data.csv"))
    
    return(qp.ref.dt)
}

extractPerPositionInfo <- function(qp.ref.dt, output.prefix, K.only = TRUE){
    mm.qp.seq.dt <- qp.ref.dt |>
        mutate(
            full_sequence_for_PTM_pos = str_replace_all(full_sequence, "[[:alpha:]]\\[.*?\\]", "x")
        )

    mm.qp.seq.filtered.dt <- mm.qp.seq.dt |>
        mutate(
            peak_id = row_number(),
            PTM_list = str_extract_all(full_sequence, "\\[.*?\\]"),
            PTM_pos = str_locate_all(full_sequence_for_PTM_pos, "x") %>%
                {lapply(., function(x) x[,2])}
            ## map(~ .x[,2]) # ?? 
        ) 

    mm.qp.seq.aapos.dt <- mm.qp.seq.filtered.dt[rep(1:nrow(mm.qp.seq.filtered.dt), times = peptide_end - peptide_start + 1)]

    mm.qp.seq.aapos.dt[, aa_rel_pos := 1:.N, by = peak_id]

    mm.qp.seq.aapos.dt[, `:=`(
        aa = substr(base_sequence, start = aa_rel_pos, stop = aa_rel_pos),
        aa_pos = peptide_start + aa_rel_pos - 1
    )]

    mm.qp.seq.long.dt <- mm.qp.seq.aapos.dt[, `:=`(
        ptm = unlist(PTM_list)[unlist(PTM_pos) %in% aa_rel_pos]
    ), by = seq_len(nrow(mm.qp.seq.aapos.dt))]

    if(K.only){
        mm.qp.seq.long.dt <- mm.qp.seq.long.dt[aa == "K"]
    } else {}

    fwrite(mm.qp.seq.long.dt, file = paste0(output.prefix, "all_per_pos_data.csv"))

    return(mm.qp.seq.long.dt)
}

dePropionylateData <- function(mm.qp.seq.long.dt){
    mm.qp.seq.long.dt <- copy(mm.qp.seq.long.dt)[
      , ptm := str_replace_all(
            ptm,
            c(
                "Less Common:Propionylation on K" = NA,
                "Custom:Oxidised Propionylation on K" = "Common Biological:Hydroxylation on K",
                "Custom:Hydroxylation-Propionylation on K" = "Common Biological:Hydroxylation on K",
                "Common Biological:Butyrylation on K" = "Common Biological:Methylation on K",
                "Common Biological:Hydroxybutyrylation on K" = "Custom:Oxidised methylation on K",
                "\\[|\\]" = "",
                "^.*:" = ""
            )
        )
    ]

    return(mm.qp.seq.long.dt)
}

calculateStoichiometryFromPerPosData <- function(mm.qp.seq.long.dt, output.prefix){

    mm.per.aa.ptm.stoichiomtery.dt <-  mm.qp.seq.long.dt[, list(
        sum_peak_intensity = sum(peak_intensity),
        sum_PSMs_mapped = sum(PSMs_mapped),
        the_number_of_peptide = length(unique(peak_id))
    ), by = list(sample_name, protein_accession, gene_name, aa, aa_pos, ptm)]

    mm.per.aa.ptm.stoichiomtery.dt[, `:=`(
        sum_intensity_per_position = sum(sum_peak_intensity),
        sum_PSMs_mapped_per_position = sum(sum_PSMs_mapped),
        sum_the_number_of_peptide = sum(the_number_of_peptide)
    ), by = list(sample_name, protein_accession, gene_name, aa, aa_pos)]

    mm.per.aa.ptm.stoichiomtery.dt[, `:=`(
        stoichiometry = sum_peak_intensity / sum_intensity_per_position
    )]

    fwrite(
        mm.per.aa.ptm.stoichiomtery.dt,
        paste0(output.prefix, "hydroxylation_stoichiomtery.csv")
    )

    return(mm.per.aa.ptm.stoichiomtery.dt)
}

calculateStoichiometry <- function(mq.ev.input, sample.info.file, ref.protein.dt, output.prefix, K.only = TRUE){
    qp.dt <- readMaxQuantData(
        mq.ev.input = mq.ev.input,
        sample.info.file = sample.info.file,
        mod.dt = data.table( # make a csv table for this; do need []
            mq_mod_name = c("\\[Propionylation\\]", "\\[Oxidised Propionylation\\]", "\\[Oxidation \\(K\\)\\]"),
            mod_name = c("[Less Common:Propionylation on K]", "[Custom:Hydroxylation-Propionylation on K]", "[Common Biological:Hydroxylation on K]")
        ),
        output.prefix = output.prefix
    )

    mq.filtered.dt <- preprocessPeakData(
        qp.dt = qp.dt,
        ref.protein.dt = ref.protein.dt,
        output.prefix = output.prefix
    )


    mm.qp.seq.long.dt <- extractPerPositionInfo(
        qp.ref.dt = mq.filtered.dt,
        output.prefix = output.prefix,
        K.only = K.only
    )

    mm.qp.seq.long.dt <- dePropionylateData(mm.qp.seq.long.dt)

    mm.per.aa.ptm.stoichiomtery.dt <- calculateStoichiometryFromPerPosData(
        mm.qp.seq.long.dt = mm.qp.seq.long.dt,
        output.prefix = output.prefix
    )
       
    return(mm.per.aa.ptm.stoichiomtery.dt)
}

addAnnotation <- function(mm.per.aa.ptm.stoichiomtery.dt, protein.feature.per.pos.dt, output.prefix){
    setnames(
        protein.feature.per.pos.dt,
        c("uniprot_id", "position"),
        c("protein_accession", "aa_pos")
    )

    protein.feature.per.pos.dt <- protein.feature.per.pos.dt[
    , .(protein_accession, aa_pos, IUPRED2, K_ratio_score, Window)
    ]

    protein.feature.per.pos.dt <-
        protein.feature.per.pos.dt[
            protein_accession %in% mm.per.aa.ptm.stoichiomtery.dt[, unique(protein_accession)]
        ]

    protein.feature.per.pos.dt[, `:=`(
        M1 = case_when(
            Window == "" ~ FALSE,
            grepl("^(.{4}M.{1}|.{6}M)", Window) ~ TRUE,
            TRUE ~ FALSE
        ),
        M2 = case_when(
            Window == "" ~ FALSE,
            grepl("^(.{3}M.{1}|.{7}M)", Window) ~ TRUE,
            TRUE ~ FALSE
        )
    )]

    mm.per.aa.ptm.stoichiomtery.dt <- mm.per.aa.ptm.stoichiomtery.dt[aa == "K"]

    mm.per.aa.ptm.stoichiomtery.with.annotation.dt <- merge(
        mm.per.aa.ptm.stoichiomtery.dt,
        protein.feature.per.pos.dt,
        by = c("protein_accession", "aa_pos")
    )

    fwrite(mm.per.aa.ptm.stoichiomtery.with.annotation.dt, paste0(output.prefix, "stoichiometry_per_position_with_annotation.csv"))

    return(mm.per.aa.ptm.stoichiomtery.with.annotation.dt)
}
