#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# SUSIE_COLOCALIZATION worker script. Lives in bin/ (not inline in the .nf file)
# because the JVM caps a Groovy string constant at 65,535 bytes. All Nextflow
# inputs arrive as NF_* environment variables exported by the process block.
# ----------------------------------------------------------------------------
nf_env <- function(k) { v <- Sys.getenv(paste0("NF_", k), unset = NA); if (is.na(v)) stop(sprintf("environment variable NF_%s not set (export it in the process block)", k)); v }
NF_susie_lambda_warn <- nf_env("susie_lambda_warn")
NF_susie_abf_fallback <- nf_env("susie_abf_fallback")
NF_susie_hybrid <- nf_env("susie_hybrid")
NF_exposure_sdY <- nf_env("exposure_sdY")
NF_outcome_sdY <- nf_env("outcome_sdY")
NF_exposure_ld_weights <- nf_env("exposure_ld_weights")
NF_ref_index_file <- nf_env("ref_index_file")
NF_ref_data_file <- nf_env("ref_data_file")
NF_ref_pop_desc_file <- nf_env("ref_pop_desc_file")
NF_af1_cutoff <- nf_env("af1_cutoff")
NF_gene_map <- nf_env("gene_map")
NF_outcome_raw <- nf_env("outcome_raw")
NF_do_liftover <- nf_env("do_liftover")
NF_chain_file <- nf_env("chain_file")
NF_gene_conversion_table <- nf_env("gene_conversion_table")
NF_proxy_genes <- nf_env("proxy_genes")
NF_sample_size <- nf_env("sample_size")
NF_window_size <- suppressWarnings(as.numeric(nf_env("window_size")))
NF_coloc_p12 <- suppressWarnings(as.numeric(nf_env("coloc_p12")))
NF_susie_h4_threshold <- suppressWarnings(as.numeric(nf_env("susie_h4_threshold")))
# v2.2 inputs (optional: a process block that does not export them gets the defaults)
nf_env_opt <- function(k, default) { v <- Sys.getenv(paste0("NF_", k), unset = NA); if (is.na(v) || !nzchar(v)) default else v }
NF_bin_dir <- nf_env_opt("bin_dir", "")
NF_coloc_rule <- nf_env_opt("coloc_rule", "H4 > 0.8")
NF_coloc_cond_threshold <- suppressWarnings(as.numeric(nf_env_opt("coloc_cond_threshold", "0.7")))
NF_coloc_power_min <- suppressWarnings(as.numeric(nf_env_opt("coloc_power_min", "0.5")))
NF_coloc_h3_strong <- suppressWarnings(as.numeric(nf_env_opt("coloc_h3_strong", "0.8")))
NF_coloc_p12_grid <- nf_env_opt("coloc_p12_grid", "1e-5,5e-6,1e-6")
NF_ld_check_r2 <- suppressWarnings(as.numeric(nf_env_opt("ld_check_r2", "0.8")))
NF_coloc_window_center <- tolower(trimws(nf_env_opt("coloc_window_center", "lead_instrument")))
NF_outcome_type <- tolower(trimws(nf_env_opt("outcome_type", "quant")))
NF_outcome_case_prop <- suppressWarnings(as.numeric(nf_env_opt("outcome_case_prop", "")))
    log_con <- file("susie_coloc_log.txt", open = "w")
    lg <- function(...) { msg <- sprintf(...); cat(msg, "\n"); writeLines(msg, log_con) }

    suppressMessages({
        library(coloc); library(susieR); library(TwoSampleMR); library(arrow)
        library(data.table); library(dplyr); library(tidyr); library(stringr); library(reshape2)
    })

    # Package provenance for this run. Written FIRST so the file exists even when
    # one of the hard setup checks below quits early.
    writeLines(capture.output(sessionInfo()), "sessionInfo_susie.txt")

    # ---- shared decision metrics / rule (bin/coloc_decision.R) -------------
    # Same functions as the classic COLOCALIZATION process, so both tables carry
    # identical derived columns (H3_plus_H4, H4_cond, ...) and the same rule.
    helper <- file.path(NF_bin_dir, "coloc_decision.R")
    if (!nzchar(NF_bin_dir) || !file.exists(helper)) helper <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "coloc_decision.R")
    if (!file.exists(helper)) helper <- Sys.which("coloc_decision.R")
    if (!nzchar(helper) || !file.exists(helper)) stop("bin/coloc_decision.R not found next to susie_colocalization.R")
    source(helper)
    p12_grid <- coloc_parse_p12_grid(NF_coloc_p12_grid)

    # Prior sensitivity for a coloc summary table: re-run the SAME coloc call at every
    # p12 of the grid (rerun_fn(p) returns its $summary) and attach H4_p12_<tag> /
    # H4cond_p12_<tag> per row, keyed by (idx1, idx2) when present, else by row order.
    # p12_min_pass / p12_max_pass come from coloc::sensitivity()'s p12 scan of `obj`.
    attach_p12_cols <- function(sm, rerun_fn, obj = NULL, rule = NF_coloc_rule) {
        # keying mode is decided ONCE from `sm` and applied to every re-run
        use_idx <- all(c("idx1","idx2") %in% names(sm)) && !all(is.na(sm$idx1))
        key_of <- function(d) if (use_idx && all(c("idx1","idx2") %in% names(d))) paste(d$idx1, d$idx2) else as.character(seq_len(nrow(d)))
        k0 <- key_of(sm)
        for (p in p12_grid) {
            tag <- coloc_p12_tag(p)
            s2 <- tryCatch(rerun_fn(p), error = function(e) { lg("  p12=%g re-run failed: %s", p, e$message); NULL })
            h4 <- rep(NA_real_, nrow(sm)); hc <- rep(NA_real_, nrow(sm))
            if (!is.null(s2) && nrow(s2) > 0 && all(c("PP.H3.abf","PP.H4.abf") %in% names(s2))) {
                s2 <- as.data.frame(s2)
                m <- match(k0, key_of(s2))
                h3v <- suppressWarnings(as.numeric(s2$PP.H3.abf[m])); h4v <- suppressWarnings(as.numeric(s2$PP.H4.abf[m]))
                h4 <- h4v
                hc <- ifelse(is.finite(h3v + h4v) & (h3v + h4v) > 1e-12, h4v / (h3v + h4v), NA_real_)
            }
            sm[[paste0("H4_p12_", tag)]] <- h4
            sm[[paste0("H4cond_p12_", tag)]] <- hc
        }
        sm$p12_min_pass <- rep(NA_real_, nrow(sm)); sm$p12_max_pass <- rep(NA_real_, nrow(sm))
        if (!is.null(obj) && "sensitivity" %in% getNamespaceExports("coloc")) {
            for (r in seq_len(nrow(sm))) {
                sens <- tryCatch(suppressMessages(coloc::sensitivity(obj, rule = rule, row = r, doplot = FALSE)),
                                 error = function(e) NULL)
                if (!is.null(sens) && all(c("p12","pass") %in% names(sens)) && any(sens$pass %in% TRUE)) {
                    sm$p12_min_pass[r] <- min(sens$p12[sens$pass %in% TRUE])
                    sm$p12_max_pass[r] <- max(sens$p12[sens$pass %in% TRUE])
                }
            }
        }
        sm
    }
    if (!(NF_coloc_window_center %in% c("lead_instrument", "first_instrument"))) {
        lg("WARNING: unknown coloc_window_center '%s'; using 'lead_instrument'", NF_coloc_window_center)
        NF_coloc_window_center <- "lead_instrument"
    }
    if (!(NF_outcome_type %in% c("quant", "cc"))) { lg("WARNING: unknown outcome_type '%s'; using 'quant'", NF_outcome_type); NF_outcome_type <- "quant" }
    if (NF_outcome_type == "cc" && (!is.finite(NF_outcome_case_prop) || NF_outcome_case_prop <= 0 || NF_outcome_case_prop >= 1))
        stop("--outcome_type cc requires --outcome_case_prop (fraction of cases, 0 < s < 1)")
    lg("Decision rule: %s | strong tier PP.H4 > %s | conditional tier H4/(H3+H4) >= %s with H3+H4 >= %s | distinct H3 > %s | LD check r2 >= %s | window centre: %s | outcome type: %s",
       NF_coloc_rule, NF_susie_h4_threshold, NF_coloc_cond_threshold, NF_coloc_power_min, NF_coloc_h3_strong,
       NF_ld_check_r2, NF_coloc_window_center, NF_outcome_type)

    # Runtime knobs (rendered by Nextflow)
    lam_warn <- suppressWarnings(as.numeric(NF_susie_lambda_warn))
    if (!is.finite(lam_warn)) lam_warn <- 0.2
    abf_fb <- suppressWarnings(as.logical(NF_susie_abf_fallback))
    if (is.na(abf_fb)) abf_fb <- TRUE
    # Hybrid coloc.susie_bf (SuSiE credible sets of the resolvable trait x single-signal
    # ABF of the other) is tried before the plain coloc.abf fallback.
    hyb_on <- suppressWarnings(as.logical(NF_susie_hybrid))
    if (is.na(hyb_on)) hyb_on <- TRUE

    # ---- sdY handling ------------------------------------------------
    #   exposure : GTEx expression is inverse-normal transformed, so the effect
    #              sizes are already on a unit-variance scale -> sdY = exposure_sdY.
    #   outcome  : a finite --outcome_sdY is used as given ("fixed:<v>"); otherwise
    #              coloc estimates sdY from the GWAS MAF + N ("MAF+N"). The chosen
    #              mode is logged and recorded per gene in susie_gene_status.tsv.
    # use_maf_n / sdy_mode are finalised once the outcome file is loaded (the
    # allele-frequency column has to exist for the MAF+N route).
    sdY_exp <- suppressWarnings(as.numeric(NF_exposure_sdY))
    if (!is.finite(sdY_exp) || sdY_exp <= 0) {
        lg("WARNING: exposure_sdY='%s' is not a positive number; using sdY = 1", NF_exposure_sdY)
        sdY_exp <- 1
    }
    sdY_out   <- suppressWarnings(as.numeric(NF_outcome_sdY))
    use_maf_n <- !is.finite(sdY_out) || sdY_out <= 0
    sdy_mode  <- if (use_maf_n) "MAF+N" else sprintf("fixed:%g", sdY_out)
    if (NF_outcome_type == "cc") {
        # a case/control dataset carries `s` (case fraction) instead of sdY
        use_maf_n <- FALSE; sdy_mode <- sprintf("cc:s=%g", NF_outcome_case_prop)
    }

    # ---- LD weighting mode -------------------------------------------
    # "outcome"  : one ancestry-weighted LD matrix (weights from the GWAS allele
    #              frequencies) used for BOTH traits - the previous behaviour.
    # "exposure" : additionally infer ancestry weights from the GTEx allele
    #              frequencies and use that LD matrix for the eQTL trait, while the
    #              GWAS keeps the outcome-weighted one.
    ld_wgt_mode <- tolower(trimws(NF_exposure_ld_weights))
    if (!(ld_wgt_mode %in% c("outcome","exposure"))) {
        lg("WARNING: unknown exposure_ld_weights '%s'; falling back to 'outcome'", ld_wgt_mode)
        ld_wgt_mode <- "outcome"
    }
    pop_wgt_exp <- NULL   # exposure-ancestry weights (only built in "exposure" mode)

    # Column schema shared by coloc.susie rows and ABF-fallback rows.
    std_cols <- c("nsnps","hit1","hit2","PP.H0.abf","PP.H1.abf","PP.H2.abf",
                  "PP.H3.abf","PP.H4.abf","idx1","idx2")

    # ---- per-gene status record -------------------------------------
    # One row per INPUT gene, written to susie_gene_status.tsv, so no gene is ever
    # silently dropped: it records why a gene produced (or failed to produce) rows,
    # plus the SuSiE diagnostics needed to tell a purity/coverage problem from an
    # LD-mismatch problem.
    new_status <- function(gene_id = NA_character_, gene_name = NA_character_,
                           reason = NA_character_) {
        data.frame(
            Gene_ID = as.character(gene_id), Gene_Name = as.character(gene_name),
            reason = as.character(reason), n_snps = NA_integer_,
            minp_exp = NA_real_, minp_out = NA_real_,
            lambda_exp = NA_real_, lambda_out = NA_real_, ld_mismatch_flag = NA, flip_applied = NA,
            n_cs_exp = NA_integer_, n_cs_out = NA_integer_,
            n_cs_unpruned_exp = NA_integer_, n_cs_unpruned_out = NA_integer_,
            n_cs_cov90_exp = NA_integer_, n_cs_cov90_out = NA_integer_,
            min_purity_exp = NA_real_, min_purity_out = NA_real_,
            max_pip_exp = NA_real_, max_pip_out = NA_real_,
            converged_exp = NA, converged_out = NA,
            n_susie_rows = 0L, hybrid_susie_bf = FALSE, abf_fallback = FALSE,
            ld_check_max_r2 = NA_real_, ld_check_pass = NA, ld_check_n_instr = 0L,
            window_center_snp = NA_character_, window_center_pos = NA_real_,
            outcome_sdY_mode = as.character(sdy_mode),
            stringsAsFactors = FALSE)
    }
    empty_status <- new_status()[0, , drop = FALSE]

    # Gene-level summary schema (one row per gene with rows; written by
    # write_gene_summary() below and, empty, by empty_out()).
    gene_summary_cols <- c("Gene_ID","Gene_Name","LD_type","nsnps","hit1","hit2",
                           "PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf",
                           "H3_plus_H4","H4_cond","H4_H3_odds","log2_H4_H3","dominant_hyp",
                           "pass_rule","coloc_class","n_pairs","n_pairs_pass","any_pass",
                           "any_pass_susie","max_H4_susie","max_H4_any",
                           "reason","lambda_exp","lambda_out","ld_mismatch_flag",
                           "n_cs_exp","n_cs_out","ld_check_max_r2","ld_check_pass")

    empty_out <- function(status = empty_status) {
        d <- data.frame(Gene_Name=character(), Gene_ID=character(), LD_type=character(),
                        nsnps=numeric(), hit1=character(), hit2=character(),
                        PP.H0.abf=numeric(), PP.H1.abf=numeric(), PP.H2.abf=numeric(),
                        PP.H3.abf=numeric(), PP.H4.abf=numeric(),
                        idx1=numeric(), idx2=numeric(), stringsAsFactors=FALSE)
        d <- coloc_derive_metrics(d)
        d$pass_rule <- logical(0); d$coloc_class <- character(0)
        saveRDS(d, "susie_coloc_results_summary.rds")
        data.table::fwrite(d, "susie_coloc_results_summary.tsv", sep = "\t")
        saveRDS(d[0,], "susie_passing_genes.rds")
        gs <- as.data.frame(setNames(replicate(length(gene_summary_cols), character(0), simplify = FALSE), gene_summary_cols))
        data.table::fwrite(gs, "susie_gene_summary.tsv", sep = "\t")
        # The status table is a declared output, so it must always exist.
        data.table::fwrite(status, "susie_gene_status.tsv", sep = "\t")
    }

    # ---- GAUSS multi-ethnic LD setup ----
    if (!requireNamespace("gauss", quietly = TRUE)) {
        lg("ERROR: package 'gauss' not installed (github: statsleelab/gauss). Cannot compute multi-ethnic LD.")
        empty_out(); close(log_con); quit(save="no", status=1)
    }
    ref_ok <- all(file.exists(NF_ref_index_file), file.exists(NF_ref_data_file), file.exists(NF_ref_pop_desc_file))
    if (!ref_ok) {
        lg("ERROR: GAUSS reference panel files not found:")
        lg("  index: %s", NF_ref_index_file); lg("  data:  %s", NF_ref_data_file); lg("  desc:  %s", NF_ref_pop_desc_file)
        empty_out(); close(log_con); quit(save="no", status=1)
    }

    # Fixed ancestry proportions (from the GWAS cohort Ns), e.g.
    # "EUR:0.86,AMR:0.08,AFR:0.05,EAS:0.002,SAS:0.0". Labels must match pop_desc.
    parse_props <- function(s) {
        parts <- strsplit(trimws(s), ",")[[1]]
        pops <- c(); wts <- c()
        for (p in parts) {
            kv <- strsplit(trimws(p), ":")[[1]]
            if (length(kv) == 2) { pops <- c(pops, trimws(kv[1])); wts <- c(wts, as.numeric(kv[2])) }
        }
        df <- data.frame(pop = pops, wgt = wts, stringsAsFactors = FALSE)
        df <- df[!is.na(df$wgt) & df$wgt > 0, ]
        df$wgt <- df$wgt / sum(df$wgt)
        df
    }
    # pop_wgt_df (ancestry proportions) is computed automatically via gauss::afmix
    # AFTER the outcome GWAS and rsID map are loaded (afmix needs rsid/chr/bp/a1/a2/af1).
    # If --ancestry_props is provided as a non-empty manual override, that is used instead.
    pop_wgt_df <- NULL

    af1c <- suppressWarnings(as.numeric(NF_af1_cutoff))
    if (is.na(af1c) || af1c < 0) af1c <- NULL

    # -----------------------------------------------------------------
    # rsID lookup built FROM THE GAUSS INDEX ITSELF.
    # Neither the GTEx association files nor the GWAS carry rsIDs (both use
    # chr_pos_ref_alt), but GAUSS matches its reference on rsID. Building the
    # map from the panel's own index guarantees the rsIDs AND the genome build
    # match the panel by construction.
    # Index columns (no header): rsid chr pos a1 a2 af ...
    # -----------------------------------------------------------------
    lg("Loading rsID map from GAUSS index (this may take a minute)...")
    idx <- tryCatch({
        # Avoid fread's gz path (needs R.utils). Prefer piping via zcat, which is
        # fast and dependency-free; fall back to base gzfile() if zcat is absent.
        idx_path <- NF_ref_index_file
        if (grepl("\\.gz$", idx_path) && nzchar(Sys.which("zcat"))) {
            as.data.frame(data.table::fread(cmd = paste("zcat", shQuote(idx_path)),
                           header = FALSE, select = 1:5,
                           col.names = c("rsid","chr","pos","a1","a2"),
                           showProgress = FALSE))
        } else if (grepl("\\.gz$", idx_path)) {
            # Base-R fallback: read all columns, keep the first five.
            d <- read.table(gzfile(idx_path), header = FALSE, stringsAsFactors = FALSE,
                            comment.char = "", quote = "")
            d <- d[, 1:5]
            names(d) <- c("rsid","chr","pos","a1","a2"); d
        } else {
            as.data.frame(data.table::fread(idx_path, header = FALSE, select = 1:5,
                           col.names = c("rsid","chr","pos","a1","a2"),
                           showProgress = FALSE))
        }
    }, error = function(e) { lg("ERROR reading index: %s", e$message); NULL })
    if (is.null(idx)) { empty_out(); close(log_con); quit(save="no", status=1) }
    idx <- idx[idx$rsid != "." & !is.na(idx$rsid), ]
    idx$chr <- sub("^chr", "", as.character(idx$chr))
    idx$pos <- suppressWarnings(as.numeric(idx$pos))
    idx <- idx[!is.na(idx$pos), ]
    idx$a1 <- toupper(idx$a1); idx$a2 <- toupper(idx$a2)
    # Two keys: alleles as-listed and swapped (our ref/alt orientation may differ)
    idx$posi <- format(as.integer(idx$pos), scientific = FALSE, trim = TRUE)
    idx$key1 <- paste(idx$chr, idx$posi, idx$a1, idx$a2, sep = ":")
    idx$key2 <- paste(idx$chr, idx$posi, idx$a2, idx$a1, sep = ":")
    rs_map <- c(setNames(idx$rsid, idx$key1), setNames(idx$rsid, idx$key2))
    lg("rsID map entries: %d (from %d index rows)", length(rs_map), nrow(idx))
    rm(idx); invisible(gc())

    lookup_rsid <- function(chr_i, pos_v, a1_v, a2_v) {
        pos_i <- format(as.integer(pos_v), scientific = FALSE, trim = TRUE)
        k <- paste(sub("^chr", "", as.character(chr_i)), pos_i, toupper(a1_v), toupper(a2_v), sep = ":")
        unname(rs_map[k])
    }

    # ---- Load inputs ----
    # Exposure = GTEx all-pairs association files (PARQUET), read directly as a lazy
    # arrow Dataset so each gene pulls only its own rows (no full-file load).
    assoc_paths <- list.files(".", pattern = "\\.parquet$", full.names = TRUE)
    if (length(assoc_paths) == 0) {
        lg("ERROR: no .parquet association files staged"); empty_out(); close(log_con); quit(save="no", status=1)
    }
    lg("Association parquet files: %d", length(assoc_paths))
    expo_ds <- arrow::open_dataset(assoc_paths, format = "parquet")
    lg("Association columns: %s", paste(names(expo_ds), collapse=","))
    gcol <- if ("gene_id" %in% names(expo_ds)) "gene_id" else "phenotype_id"

    # -----------------------------------------------------------------
    # Gene-ID bridge between the exposure (eGene) GENCODE release and the
    # all-pairs release, built by GENE_RELEASE_MAP. Keyed by version-stripped
    # exposure ENSG; the value is the EXACT versioned id used in the parquet,
    # so the read is an exact match rather than a grepl on a stripped ENSG.
    # -----------------------------------------------------------------
    gmap <- tryCatch(readRDS(NF_gene_map), error = function(e) NULL)
    if (is.null(gmap)) lg("WARNING: gene_release_map could not be read; every gene will look unmapped")
    assoc_id_of   <- if (!is.null(gmap)) setNames(as.character(gmap$assoc_gene_id),   as.character(gmap$exp_core)) else character(0)
    map_method_of <- if (!is.null(gmap)) setNames(as.character(gmap$map_method),      as.character(gmap$exp_core)) else character(0)
    assoc_name_of <- if (!is.null(gmap)) setNames(as.character(gmap$assoc_gene_name), as.character(gmap$exp_core)) else character(0)

    method_for <- function(ens) {
        if (ens %in% names(map_method_of)) map_method_of[[ens]] else NA_character_
    }

    read_gene_assoc <- function(ens) {
        aid <- if (ens %in% names(assoc_id_of)) assoc_id_of[[ens]] else NA_character_
        if (is.null(aid) || length(aid) == 0 || is.na(aid) || !nzchar(aid)) return(data.frame())
        mm <- method_for(ens)
        if (!is.na(mm) && mm != "id_match") {
            lg("  gene id bridged v11->v10 via %s: %s -> %s (%s)", mm, ens, aid,
               if (ens %in% names(assoc_name_of) && !is.na(assoc_name_of[[ens]])) assoc_name_of[[ens]] else "NA")
        }
        as.data.frame(expo_ds %>% filter(.data[[gcol]] == aid) %>% collect())
    }
    outcome  <- readRDS(NF_outcome_raw)
    if (!"se"   %in% colnames(outcome) && "standard_error" %in% colnames(outcome)) outcome$se   <- outcome$standard_error
    if (!"pval" %in% colnames(outcome) && "p_value"        %in% colnames(outcome)) outcome$pval <- outcome$p_value
    if (!"sample_size" %in% colnames(outcome)) {
        for (alt in c("TotalSampleSize","n","N","SampleSize","n_complete_samples")) {
            if (alt %in% colnames(outcome)) { outcome$sample_size <- outcome[[alt]]; break }
        }
    }
    if (!"sample_size" %in% colnames(outcome)) outcome$sample_size <- NA_real_

    # ---- Outcome allele frequency ------------------------------------
    # Used for the ancestry (afmix) fit and, when --outcome_sdY is not a number,
    # for the MAF coloc needs to estimate sdY. Detected once and carried into every
    # cis-window as eaf.outcome (harmonise_data keeps that column).
    out_af_col <- intersect(c("effect_allele_frequency","eaf","af","EAF","AF","freq"),
                            colnames(outcome))[1]
    outcome$eaf_out_raw <- if (!is.na(out_af_col)) suppressWarnings(as.numeric(outcome[[out_af_col]])) else NA_real_
    if (use_maf_n && is.na(out_af_col)) {
        lg("WARNING: --outcome_sdY not given and the outcome GWAS has no allele-frequency column")
        lg("         (looked for effect_allele_frequency/eaf/af/EAF/AF/freq); falling back to sdY = 1.")
        use_maf_n <- FALSE; sdY_out <- 1; sdy_mode <- "fixed:1(no_af)"
    }
    lg("Outcome sdY mode: %s%s", sdy_mode,
       if (use_maf_n) sprintf(" (allele frequency column '%s')", out_af_col) else "")
    lg("Exposure sdY: %g", sdY_exp)

    # -----------------------------------------------------------------
    # b38 -> b37 liftover (for LD lookup ONLY).
    # The GAUSS 33k panel is GRCh37 while GTEx/GWAS here are GRCh38, so panel
    # lookups need b37 coordinates. Summary stats, z-scores and all reported
    # results stay in b38 - we translate coordinates purely to find the right
    # variants in the reference panel.
    # -----------------------------------------------------------------
    do_lift <- as.logical(NF_do_liftover)
    if (is.na(do_lift)) do_lift <- FALSE
    chain <- NULL
    if (do_lift) {
        if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            lg("ERROR: rtracklayer not installed but liftover requested.")
            lg("       BiocManager::install('rtracklayer')")
            empty_out(); close(log_con); quit(save="no", status=1)
        }
        chain_path <- NF_chain_file
        if (!file.exists(chain_path)) {
            lg("ERROR: chain file not found: %s", chain_path)
            lg("       Get hg38ToHg19.over.chain.gz from UCSC and set --chain_file")
            empty_out(); close(log_con); quit(save="no", status=1)
        }
        # rtracklayer::import.chain needs an uncompressed chain
        if (grepl("\\.gz$", chain_path)) {
            unz <- file.path(getwd(), "chain_unzipped.chain")
            system2("zcat", args = shQuote(chain_path), stdout = unz)
            chain_path <- unz
        }
        chain <- tryCatch(rtracklayer::import.chain(chain_path),
                          error = function(e) { lg("ERROR importing chain: %s", e$message); NULL })
        if (is.null(chain)) { empty_out(); close(log_con); quit(save="no", status=1) }
        lg("Liftover enabled: GRCh38 -> GRCh37 for LD panel lookup")
    } else {
        lg("Liftover disabled: assuming data and LD panel share a genome build")
    }

    # Lift a vector of b38 positions on one chromosome -> b37.
    # Returns a vector the same length, NA where the position does not lift.
    lift_positions <- function(chr_i, pos_v) {
        if (!do_lift) return(pos_v)
        gr <- GenomicRanges::GRanges(
            seqnames = paste0("chr", sub("^chr", "", as.character(chr_i))),
            ranges = IRanges::IRanges(start = as.integer(pos_v), width = 1)
        )
        gr$idx <- seq_along(pos_v)
        lifted <- tryCatch(rtracklayer::liftOver(gr, chain), error = function(e) NULL)
        if (is.null(lifted)) return(rep(NA_real_, length(pos_v)))
        out <- rep(NA_real_, length(pos_v))
        # liftOver returns a GRangesList; keep only 1:1 unique mappings
        lens <- S4Vectors::elementNROWS(lifted)
        ok <- which(lens == 1)
        if (length(ok) > 0) {
            flat <- unlist(lifted[ok])
            out[gr$idx[ok]] <- GenomicRanges::start(flat)
        }
        out
    }

    # -----------------------------------------------------------------
    # Ancestry proportions via gauss::afmix (always inferred from the data;
    # the manual --ancestry_props override was removed by request).
    # afmix needs an input file with columns: rsid chr bp a1 a2 af1, matched to the
    # reference panel. A genome-wide thinned subset is enough and much faster.
    # dd must carry: chr (no "chr" prefix), bp (b38), a1 = REF, a2 = ALT,
    # af1 = frequency of a1(REF). Returns data.frame(pop, wgt) or NULL.
    # -----------------------------------------------------------------
    afmix_weights <- function(dd, label, infile, max_snps = 50000L) {
        dd <- dd[!is.na(dd$bp) & !is.na(dd$af1) & dd$af1 > 0 & dd$af1 < 1, , drop = FALSE]
        # thin to a manageable genome-wide sample for the fit
        set.seed(1)
        if (nrow(dd) > max_snps) dd <- dd[sample(nrow(dd), max_snps), , drop = FALSE]
        # Lift b38 -> b37 FIRST (the panel is GRCh37), then look up rsIDs on b37
        # positions. The afmix input must carry b37 positions AND their matching
        # rsIDs - mixing b38 bp with b37 rsIDs makes afmix segfault.
        dd$bp37 <- NA_real_
        for (cc in unique(dd$chr)) {
            sel <- dd$chr == cc
            dd$bp37[sel] <- lift_positions(cc, dd$bp[sel])
        }
        dd <- dd[!is.na(dd$bp37), , drop = FALSE]
        dd$rsid <- lookup_rsid(dd$chr, dd$bp37, dd$a1, dd$a2)
        dd <- dd[!is.na(dd$rsid) & grepl("^rs[0-9]+$", dd$rsid), , drop = FALSE]
        dd <- dd[!duplicated(dd$rsid), , drop = FALSE]
        lg("  afmix[%s] input SNPs (rsID-matched, b37): %d", label, nrow(dd))
        if (nrow(dd) < 100) {
            lg("  afmix[%s]: too few matched SNPs (%d)", label, nrow(dd)); return(NULL)
        }
        # afmix input columns: rsid chr bp a1 a2 af1  (bp = b37 to match the panel)
        afmix_in <- data.frame(rsid = dd$rsid, chr = dd$chr,
                               bp = as.integer(dd$bp37),
                               a1 = dd$a1, a2 = dd$a2, af1 = dd$af1,
                               stringsAsFactors = FALSE)
        afmix_file <- file.path(getwd(), infile)
        write.table(afmix_in, afmix_file, sep = " ", row.names = FALSE, quote = FALSE)
        pw <- tryCatch(
            gauss::afmix(input_file = afmix_file,
                         reference_index_file = NF_ref_index_file,
                         reference_data_file  = NF_ref_data_file,
                         reference_pop_desc_file = NF_ref_pop_desc_file,
                         interval = 1000L),
            error = function(e) { lg("  afmix[%s] error: %s", label, e$message); NULL })
        if (is.null(pw)) return(NULL)
        # afmix returns a data frame of population IDs + weights, but the column
        # names/types vary (pop may be a factor, wgt may be character). computeLD
        # needs pop=character, wgt=numeric, or it fails with a type error. Coerce
        # unconditionally rather than only when columns are missing.
        cn <- colnames(pw)
        pcol <- if ("pop" %in% cn) "pop" else cn[which(sapply(pw, function(x) is.character(x) || is.factor(x)))[1]]
        wcol <- if ("wgt" %in% cn) "wgt" else cn[which(sapply(pw, is.numeric))[1]]
        pw <- data.frame(
            pop = as.character(pw[[pcol]]),
            wgt = suppressWarnings(as.numeric(as.character(pw[[wcol]]))),
            stringsAsFactors = FALSE)
        pw <- pw[!is.na(pw$wgt) & pw$wgt > 0, , drop = FALSE]
        if (nrow(pw) == 0) { lg("  afmix[%s]: no positive weights returned", label); return(NULL) }
        pw$wgt <- pw$wgt / sum(pw$wgt)
        pw
    }

    log_weights <- function(pw, label) {
        lg("Ancestry proportions [%s] (%d populations):", label, nrow(pw))
        for (r in seq_len(nrow(pw))) lg("  %s: %.4f", pw$pop[r], pw$wgt[r])
    }

    # ---- Outcome (GWAS) ancestry weights: always needed ----
    {
        lg("Inferring ancestry proportions via gauss::afmix (outcome GWAS) ...")
        # locate outcome AF + position columns
        af_col <- out_af_col
        ea_col <- intersect(c("effect_allele","EA","ALT","alt"), colnames(outcome))[1]
        oa_col <- intersect(c("other_allele","OA","REF","ref"), colnames(outcome))[1]
        ch_col <- intersect(c("chromosome","chr","CHR"), colnames(outcome))[1]
        bp_col <- intersect(c("base_pair_location","bp","BP","pos","position"), colnames(outcome))[1]
        if (any(is.na(c(af_col,ea_col,oa_col,ch_col,bp_col)))) {
            lg("  afmix: outcome missing needed columns (af/alleles/chr/bp); cannot infer ancestry.")
            empty_out(); close(log_con); quit(save="no", status=1)
        }
        od <- data.frame(
            chr = sub("^chr","", as.character(outcome[[ch_col]])),
            bp  = suppressWarnings(as.numeric(outcome[[bp_col]])),
            a1  = toupper(as.character(outcome[[oa_col]])),   # a1 = REF (panel convention)
            a2  = toupper(as.character(outcome[[ea_col]])),   # a2 = ALT
            af1 = 1 - suppressWarnings(as.numeric(outcome[[af_col]])),  # af of a1(REF) = 1 - EAF(ALT)
            stringsAsFactors = FALSE)
        pop_wgt_df <- afmix_weights(od, "outcome", "afmix_input.txt")
        if (is.null(pop_wgt_df)) {
            lg("  afmix failed for the outcome GWAS; cannot compute ancestry-weighted LD.")
            empty_out(); close(log_con); quit(save="no", status=1)
        }
    }
    saveRDS(pop_wgt_df, "ancestry_proportions_used.rds")
    log_weights(pop_wgt_df, "outcome")

    # ---- Exposure (GTEx) ancestry weights: only for --susie_exposure_ld_weights exposure ----
    # GTEx donors are ~85% EUR while the GWAS may be far more admixed, so the LD that
    # matches the eQTL statistics is not the LD that matches the GWAS. In "exposure"
    # mode a second afmix fit is run on the GTEx allele frequencies and the resulting
    # weights give the eQTL trait its own LD matrix.
    if (ld_wgt_mode == "exposure") {
        lg("Inferring ancestry proportions via gauss::afmix (exposure / GTEx allele frequencies) ...")
        if (!all(c("variant_id","af") %in% names(expo_ds))) {
            lg("  association files lack variant_id/af; cannot infer exposure ancestry - using the outcome weights for both traits")
            ld_wgt_mode <- "outcome"
        } else {
            ea <- tryCatch(
                as.data.frame(expo_ds %>% select(variant_id, af) %>% head(2e6) %>% collect()),
                error = function(e) { lg("  exposure AF read failed: %s", e$message); NULL })
            if (is.null(ea) || nrow(ea) == 0) {
                lg("  no exposure allele frequencies read; using the outcome weights for both traits")
                ld_wgt_mode <- "outcome"
            } else {
                ea <- ea[!is.na(ea$af), , drop = FALSE]
                ea <- ea[!duplicated(ea$variant_id), , drop = FALSE]
                set.seed(1)
                if (nrow(ea) > 50000) ea <- ea[sample(nrow(ea), 50000), , drop = FALSE]
                # variant_id is "chr1_100_REF_ALT_b38": a1 = REF, a2 = ALT, af = ALT AF
                vs <- reshape2::colsplit(ea$variant_id, "_",
                        c("v_chr","v_pos","v_ref","v_alt","v_build"))
                ed <- data.frame(
                    chr = sub("^chr","", as.character(vs$v_chr)),
                    bp  = suppressWarnings(as.numeric(vs$v_pos)),
                    a1  = toupper(as.character(vs$v_ref)),
                    a2  = toupper(as.character(vs$v_alt)),
                    af1 = 1 - suppressWarnings(as.numeric(ea$af)),
                    stringsAsFactors = FALSE)
                pop_wgt_exp <- afmix_weights(ed, "exposure", "afmix_input_exposure.txt")
                if (is.null(pop_wgt_exp)) {
                    lg("  exposure afmix failed; using the outcome ancestry weights for both traits")
                    ld_wgt_mode <- "outcome"
                } else {
                    saveRDS(pop_wgt_exp, "ancestry_proportions_exposure.rds")
                    log_weights(pop_wgt_exp, "exposure")
                }
            }
        }
    }
    lg("LD weighting mode: %s (exposure trait uses %s-derived ancestry weights)",
       ld_wgt_mode, ld_wgt_mode)


    gene_conv <- tryCatch(readRDS(NF_gene_conversion_table), error = function(e) NULL)
    map_gene_name <- function(ids) {
        if (is.null(gene_conv) || !all(c("id","name") %in% colnames(gene_conv))) return(as.character(ids))
        ens <- sub(".*?(ENSG[0-9]+).*", "\\1", as.character(ids))
        nm <- gene_conv$name[match(ens, gene_conv$id)]
        nm[is.na(nm)] <- as.character(ids)[is.na(nm)]
        nm
    }

    # Combine harmonized data (used only to find each gene's chr + lead position)
    har_files <- list.files(".", pattern = "harmonized_.*\\.rds$")
    har_total <- data.frame()
    for (f in har_files) {
        df <- readRDS(f)
        if (!is.null(df) && nrow(df) > 0) {
            chr_num <- as.numeric(gsub("harmonized_chr_(.*)\\.rds", "\\1", f))
            df$chromosome <- sprintf("%s", chr_num)
            har_total <- rbind(har_total, df)
        }
    }

    # SuSiE now runs on the SAME FDR-significant gene set as classic coloc (not only
    # coloc failures). The input may carry gene ids in either 'Gene_ID' (old proxy_genes
    # format) or 'id.exposure' (the significant-results format) - accept both.
    gene_list <- tryCatch(readRDS(NF_proxy_genes), error = function(e) data.frame())
    if (is.null(gene_list) || nrow(gene_list) == 0) {
        lg("No input genes for SuSiE; nothing to do."); empty_out(); close(log_con); quit(save="no", status=0)
    }
    id_col <- if ("Gene_ID" %in% colnames(gene_list)) "Gene_ID" else
              if ("id.exposure" %in% colnames(gene_list)) "id.exposure" else NA
    if (is.na(id_col)) {
        lg("Input gene list has neither Gene_ID nor id.exposure column."); empty_out(); close(log_con); quit(save="no", status=1)
    }
    fail_ids <- unique(as.character(gene_list[[id_col]]))
    merged_df <- har_total[har_total$id.exposure %in% fail_ids, ]
    lg("SuSiE input genes (FDR-significant set): %d", length(fail_ids))

    # -----------------------------------------------------------------
    # LD via GAUSS computeLD: ancestry-WEIGHTED multi-ethnic correlation
    # matrix (signed R) for one locus. Keyed by chr:pos to match harm_data$SNP.
    # gauss_in: data.frame(rsid, chr, bp, a1, a2, af1, z) for the window.
    # -----------------------------------------------------------------
    # pop_wgt selects the ancestry weighting: the outcome-derived pop_wgt_df by
    # default, or pop_wgt_exp for the exposure-weighted matrix. tag only labels the
    # written input file so both runs stay inspectable.
    get_ld_matrix_gauss <- function(gauss_in, chr_i, start_bp, end_bp,
                                    pop_wgt = pop_wgt_df, tag = "out") {
        if (nrow(gauss_in) < 3) return(NULL)
        # Write into the task work dir (not tempfile) so the input is inspectable
        # after a failure, and matches GAUSS's example format exactly:
        #   space-delimited, WITH header, columns: rsid chr bp a1 a2 af1 z
        gfile <- file.path(getwd(), sprintf("gauss_input_chr%s_%s_%s.txt", chr_i, as.integer(start_bp), tag))
        write.table(gauss_in, gfile, sep = " ", row.names = FALSE, quote = FALSE)
        lg("  wrote GAUSS input: %s (%d rows)", basename(gfile), nrow(gauss_in))
        # GAUSS's shipped examples (PGC2_Chr22_ilmn1M_Z.txt / _AF1.txt) are
        # space-delimited WITH a header. The docs say input_file needs all of
        # rsid,chr,bp,a1,a2,af1,z - but the examples split z and af1 across two
        # 6-column files. Try the documented single-file layout first, then
        # fall back to a headerless variant if the loader rejects it.
        attempt <- function(path) {
            tryCatch({
                gauss::computeLD(
                    chr = as.integer(chr_i),
                    start_bp = as.integer(start_bp), end_bp = as.integer(end_bp),
                    input_file = path,
                    reference_index_file = NF_ref_index_file,
                    reference_data_file = NF_ref_data_file,
                    reference_pop_desc_file = NF_ref_pop_desc_file,
                    af1_cutoff = af1c,
                    pop_wgt_df = pop_wgt
                )
            }, error = function(e) { lg("  computeLD error: %s", e$message); NULL })
        }

        LD <- attempt(gfile)

        if (is.null(LD)) {
            # Fallback: same columns, no header row
            gfile2 <- sub("\\.txt$", "_nohdr.txt", gfile)
            write.table(gauss_in, gfile2, sep = " ", row.names = FALSE,
                        col.names = FALSE, quote = FALSE)
            lg("  retrying computeLD without header row")
            LD <- attempt(gfile2)
        }

        if (is.null(LD)) return(NULL)
        # Per ?computeLD the return is a LIST: a data frame (rsid, chr, bp, a1, a2,
        # af1mix) AND a correlation matrix. Pull out the matrix element.
        M <- NULL
        if (is.matrix(LD)) {
            M <- LD
        } else if (is.list(LD)) {
            # prefer an explicitly named element, else the first matrix-like element
            for (nm in c("cor","LD","ld","R","ld_r","corr","cor_matrix")) {
                if (!is.null(LD[[nm]]) && is.matrix(LD[[nm]])) { M <- LD[[nm]]; break }
            }
            if (is.null(M)) {
                mats <- Filter(function(x) is.matrix(x) && nrow(x) == ncol(x), LD)
                if (length(mats) > 0) M <- mats[[1]]
            }
            # rsid labels usually live on the accompanying data frame
            gdf <- NULL
            if (!is.null(M)) {
                dfs <- Filter(function(x) is.data.frame(x) && "rsid" %in% names(x), LD)
                if (length(dfs) > 0 && nrow(dfs[[1]]) == nrow(M)) {
                    gdf <- dfs[[1]]
                    if (is.null(rownames(M))) { rownames(M) <- gdf$rsid; colnames(M) <- gdf$rsid }
                }
            }
        }
        if (is.null(M)) {
            lg("  computeLD returned no usable matrix (names: %s)",
               paste(names(LD), collapse=","))
            return(NULL)
        }
        storage.mode(M) <- "numeric"
        # Attach GAUSS's per-SNP alleles (a1/a2 it actually used) so the caller can
        # verify/correct sign orientation against the harmonised betas.
        if (exists("gdf") && !is.null(gdf)) {
            aa <- gdf[, intersect(c("rsid","a1","a2"), names(gdf)), drop = FALSE]
            attr(M, "gauss_alleles") <- aa
        }
        M
    }


    tab_res <- data.frame()
    status_all <- empty_status
    genes <- unique(merged_df$id.exposure)

    # Input genes with no harmonized rows at all never reach the loop; record them
    # so the status table really does cover every input gene.
    missing_genes <- setdiff(fail_ids, genes)
    for (mg in missing_genes) {
        status_all <- rbind(status_all,
            new_status(mg, map_gene_name(mg), "gene not in harmonized data"))
    }
    if (length(missing_genes) > 0)
        lg("%d input gene(s) absent from the harmonized data; recorded in status table", length(missing_genes))

    for (gi in seq_along(genes)) {
        gene_id <- genes[gi]
        prot <- map_gene_name(gene_id)
        lg("[%d/%d] SuSiE gene %s (%s)", gi, length(genes), gene_id, prot)

        res_gene <- tryCatch({
            entry_all <- merged_df[merged_df$id.exposure == gene_id, , drop = FALSE]
            # Window centre: the LEAD instrument (smallest exposure p) by default;
            # "first_instrument" restores the v2.1 behaviour (first harmonized row).
            # Same rule as the classic COLOCALIZATION process.
            if (NF_coloc_window_center == "lead_instrument" && "pval.exposure" %in% names(entry_all)) {
                pe <- suppressWarnings(as.numeric(entry_all$pval.exposure))
                entry <- if (any(is.finite(pe))) entry_all[which.min(pe), , drop = FALSE] else entry_all[1, , drop = FALSE]
            } else {
                entry <- entry_all[1, , drop = FALSE]
            }
            chrI <- as.character(entry$chromosome)
            position <- if ("position" %in% names(entry)) as.numeric(entry$position) else as.numeric(entry$base_pair_location)
            if (is.na(position)) stop("no lead position")
            center_snp <- if ("SNP" %in% names(entry)) as.character(entry$SNP) else NA_character_
            # MR instruments of this gene (chr:pos ids, same coding as the window's SNP column)
            instr_snps <- if ("SNP" %in% names(entry_all)) unique(as.character(entry_all$SNP)) else character(0)
            pos_min <- position - NF_window_size; pos_max <- position + NF_window_size
            chr_label <- paste0("chr", chrI)
            lg("  window centre: %s (%s) at chr%s:%.0f; %d instrument(s)",
               NF_coloc_window_center, center_snp, chrI, position, length(instr_snps))

            # ---- Exposure cis-window from the FULL raw exposure file (John's way) ----
            # v11 all-pairs columns: gene_id, variant_id="chr1_100_REF_ALT_b38",
            # tss_distance, af, ma_samples, ma_count, pval_nominal, slope, slope_se.
            # (No maf/ref_factor: allele freq is `af` = ALT AF; effect allele = ALT.)
            ens <- sub(".*?(ENSG[0-9]+).*", "\\1", gene_id)
            mm_g <- method_for(ens)
            if (is.na(mm_g) || mm_g == "unmapped")
                stop("unmapped in all-pairs release (see gene_release_map.tsv)")
            expo_gene <- read_gene_assoc(ens)
            if (nrow(expo_gene) == 0) stop("gene not found in association files")

            vsplit <- reshape2::colsplit(expo_gene$variant_id, "_",
                        c("v_chr","v_pos","v_ref","v_alt","v_build"))
            expo_gene$v_chr <- as.character(vsplit$v_chr)
            expo_gene$v_pos <- suppressWarnings(as.numeric(vsplit$v_pos))
            expo_gene$v_ref <- toupper(as.character(vsplit$v_ref))
            expo_gene$v_alt <- toupper(as.character(vsplit$v_alt))

            # Effect-allele frequency: prefer 'af' (ALT AF, v11 all-pairs); fall back to 'maf'.
            if ("af" %in% colnames(expo_gene)) {
                expo_gene$eaf_exp <- suppressWarnings(as.numeric(expo_gene$af))
            } else if ("maf" %in% colnames(expo_gene)) {
                expo_gene$eaf_exp <- suppressWarnings(as.numeric(expo_gene$maf))
                # legacy ref_factor MAF correction only applies to older independent files
                if ("ref_factor" %in% colnames(expo_gene)) {
                    flip <- which(expo_gene$ref_factor == -1)
                    if (length(flip)) expo_gene$eaf_exp[flip] <- 1 - expo_gene$eaf_exp[flip]
                }
            } else {
                expo_gene$eaf_exp <- 0.3
            }

            expo_subset <- expo_gene %>%
                filter(v_chr == chr_label, v_pos >= pos_min, v_pos <= pos_max) %>%
                transmute(
                    SNP = paste0(chrI, ":", v_pos),
                    chromosome = chrI, base_pair_location = v_pos,
                    effect_allele.exposure = v_alt, other_allele.exposure = v_ref,
                    beta.exposure = slope, se.exposure = slope_se, pval.exposure = pval_nominal,
                    eaf.exposure = eaf_exp,
                    samplesize.exposure = if ("ma_samples" %in% names(.)) as.numeric(ma_samples) else NA_real_,
                    id.exposure = gene_id, exposure = gene_id
                )

            # ---- Outcome cis-window from the FULL raw outcome file ----
            # Raw outcome chromosome is like "chr17"; match to chr_label.
            outcome_subset <- outcome %>%
                filter(chromosome == chr_label,
                       base_pair_location >= pos_min, base_pair_location <= pos_max) %>%
                transmute(
                    SNP = paste0(chrI, ":", base_pair_location),
                    base_pair_location = base_pair_location,
                    effect_allele.outcome = toupper(effect_allele), other_allele.outcome = toupper(other_allele),
                    beta.outcome = beta, se.outcome = se, pval.outcome = pval,
                    eaf.outcome = eaf_out_raw, samplesize.outcome = sample_size,
                    id.outcome = "GWAS", outcome = "GWAS"
                )

            lg("  window SNPs: exposure=%d outcome=%d", nrow(expo_subset), nrow(outcome_subset))
            if (nrow(expo_subset) == 0 || nrow(outcome_subset) == 0) stop("empty cis-window")

            is_pal <- function(a1,a2) (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")
            expo_subset <- expo_subset[!is_pal(expo_subset$effect_allele.exposure, expo_subset$other_allele.exposure), ]

            # Keep the REAL outcome frequency keyed by SNP. harmonise_data may flip
            # eaf.outcome onto the exposure allele coding, but MAF = pmin(f, 1-f) is
            # flip-invariant, so the raw value is all that is needed afterwards.
            af_lookup <- setNames(suppressWarnings(as.numeric(outcome_subset$eaf.outcome)),
                                  as.character(outcome_subset$SNP))
            # harmonise_data(action = 2) only consults eaf for palindromic SNPs, and
            # those were just dropped on the exposure side - so a placeholder where the
            # GWAS has no frequency is safe. The real value is restored from af_lookup.
            na_eaf <- !is.finite(outcome_subset$eaf.outcome)
            if (any(na_eaf)) outcome_subset$eaf.outcome[na_eaf] <- 0.3

            harm_data <- harmonise_data(expo_subset, outcome_subset, action = 2)
            lg("  harmonise returned %d rows; cols: %s", nrow(harm_data),
               paste(head(names(harm_data), 25), collapse=","))
            harm_data <- harm_data[!duplicated(harm_data$SNP), ]
            if ("mr_keep" %in% names(harm_data)) harm_data <- harm_data[harm_data$mr_keep == TRUE, ]
            if (nrow(harm_data) < 3) stop(sprintf("fewer than 3 harmonized SNPs (%d)", nrow(harm_data)))

            # ---- Outcome MAF (only when coloc has to estimate sdY) ----
            # Done here, before any LD work, so harm_data / LD / flip vectors stay
            # aligned: coloc needs a usable MAF for every SNP in the dataset.
            if (use_maf_n) {
                maf_chk <- unname(af_lookup[as.character(harm_data$SNP)])
                maf_chk <- pmin(maf_chk, 1 - maf_chk)
                keep_maf <- is.finite(maf_chk) & maf_chk > 0 & maf_chk < 1
                if (any(!keep_maf)) lg("  dropped %d SNP(s) with NA/0/1 outcome MAF", sum(!keep_maf))
                harm_data <- harm_data[keep_maf, ]
                if (nrow(harm_data) < 3)
                    stop(sprintf("fewer than 3 SNPs with usable outcome MAF (%d)", nrow(harm_data)))
            }

            # Derive position from the SNP id ("chr:pos") because harmonise_data does
            # not preserve non-standard columns like base_pair_location.
            bp_vec <- suppressWarnings(as.numeric(sub("^.*:", "", harm_data$SNP)))
            if (all(is.na(bp_vec)) && "base_pair_location" %in% names(harm_data)) {
                bp_vec <- suppressWarnings(as.numeric(harm_data$base_pair_location))
            }
            if (all(is.na(bp_vec))) stop("cannot derive SNP positions after harmonise")

            # Order by position
            ord <- order(bp_vec)
            harm_data <- harm_data[ord, ]; bp_vec <- bp_vec[ord]

            # eaf.exposure may be dropped/NA by harmonise_data - guard it
            eaf_vec <- if ("eaf.exposure" %in% names(harm_data)) suppressWarnings(as.numeric(harm_data$eaf.exposure)) else rep(NA_real_, nrow(harm_data))
            eaf_vec[is.na(eaf_vec)] <- 0.3

            z_vec <- harm_data$beta.exposure / harm_data$se.exposure

            # ---- Lift b38 -> b37 for panel lookup (stats stay in b38) ----
            bp37 <- lift_positions(chrI, bp_vec)
            n_lift <- sum(!is.na(bp37))
            lg("  lifted b38->b37: %d / %d positions", n_lift, length(bp_vec))
            if (n_lift < 3) stop(sprintf("only %d positions lifted to b37", n_lift))

            # ---- Attach real rsIDs from the GAUSS panel index (b37 coords) ----
            # Panel index is keyed a1=REF, a2=ALT. Our lookup builds both orders,
            # but pass REF first to match the panel's convention.
            rs_vec <- lookup_rsid(chrI, bp37,
                                  harm_data$other_allele.exposure,
                                  harm_data$effect_allele.exposure)
            n_matched <- sum(!is.na(rs_vec))
            lg("  rsID matched: %d / %d SNPs", n_matched, length(rs_vec))
            if (n_matched < 3) stop(sprintf("only %d SNPs matched the LD panel by rsID", n_matched))

            keep_rs <- !is.na(rs_vec)
            harm_data <- harm_data[keep_rs, ]
            bp_vec <- bp_vec[keep_rs]; eaf_vec <- eaf_vec[keep_rs]
            z_vec <- z_vec[keep_rs]; rs_vec <- rs_vec[keep_rs]
            bp37 <- bp37[keep_rs]

            # ---- Ancestry-weighted LD via GAUSS computeLD ----
            # GAUSS convention (per ?computeLD and the panel index):
            #   a1 = REFERENCE allele, a2 = alternate, af1 = frequency of a1.
            # GTEx gives effect=ALT, other=REF and `af` = ALT allele frequency,
            # so we must SWAP: a1<-other(REF), a2<-effect(ALT), af1 <- 1 - af(ALT).
            # NOTE: bp passed to GAUSS must be b37 to match the panel.
            af_alt <- eaf_vec                       # ALT allele frequency (GTEx `af`)
            af1_ref <- 1 - af_alt                   # frequency of the REFERENCE allele
            gauss_in <- data.frame(
                rsid = as.character(rs_vec),
                chr  = rep(as.integer(chrI), nrow(harm_data)),
                bp   = as.integer(bp37),
                a1   = toupper(as.character(harm_data$other_allele.exposure)),   # REF
                a2   = toupper(as.character(harm_data$effect_allele.exposure)),  # ALT
                af1  = af1_ref,
                z    = as.numeric(z_vec),
                stringsAsFactors = FALSE
            )
            # Drop duplicate rsIDs and unusable values before LD
            gauss_in <- gauss_in[!duplicated(gauss_in$rsid), ]
            gauss_in <- gauss_in[is.finite(gauss_in$z) & !is.na(gauss_in$bp), ]
            # af1 must be a valid frequency in (0,1) for GAUSS's a1 filter
            gauss_in <- gauss_in[is.finite(gauss_in$af1) & gauss_in$af1 > 0 & gauss_in$af1 < 1, ]
            # Force clean numeric types and drop any residual non-finite values in ANY
            # numeric column, so computeLD's reader never hits a character where it
            # expects a double ("Not compatible with requested type").
            gauss_in$chr <- suppressWarnings(as.integer(gauss_in$chr))
            gauss_in$bp  <- suppressWarnings(as.integer(gauss_in$bp))
            gauss_in$af1 <- suppressWarnings(as.numeric(gauss_in$af1))
            gauss_in$z   <- suppressWarnings(as.numeric(gauss_in$z))
            gauss_in$rsid <- as.character(gauss_in$rsid)
            gauss_in$a1  <- toupper(as.character(gauss_in$a1))
            gauss_in$a2  <- toupper(as.character(gauss_in$a2))
            good <- is.finite(gauss_in$chr) & is.finite(gauss_in$bp) &
                    is.finite(gauss_in$af1) & is.finite(gauss_in$z) &
                    !is.na(gauss_in$rsid) & nzchar(gauss_in$rsid) &
                    gauss_in$a1 %in% c("A","C","G","T") & gauss_in$a2 %in% c("A","C","G","T")
            gauss_in <- gauss_in[good, , drop = FALSE]
            if (nrow(gauss_in) < 3) stop(sprintf("fewer than 3 clean SNPs for LD (%d)", nrow(gauss_in)))
            lg("  gauss_in: n=%d, af1 range %.3f-%.3f, bp37 range %d-%d",
               nrow(gauss_in), min(gauss_in$af1), max(gauss_in$af1),
               min(gauss_in$bp), max(gauss_in$bp))
            harm_data <- harm_data[match(gauss_in$rsid, rs_vec), ]
            # keep an rsid column so we can align to the LD matrix afterwards
            harm_data$rsid <- gauss_in$rsid
            lg("  harmonized=%d, LD input=%d", nrow(harm_data), nrow(gauss_in))

            ld_r <- get_ld_matrix_gauss(gauss_in, chrI,
                                        min(gauss_in$bp, na.rm = TRUE),
                                        max(gauss_in$bp, na.rm = TRUE))
            if (is.null(ld_r)) stop("no LD matrix returned from GAUSS")

            # Matrix subsetting drops custom attributes, so capture GAUSS's per-SNP
            # alleles BEFORE subsetting (this is why the orientation check used to
            # report "GAUSS alleles unavailable").
            gall <- attr(ld_r, "gauss_alleles")

            # ---- Optional second LD matrix, weighted by the EXPOSURE ancestry ----
            # Same window and same SNP input, different population weights. Every
            # subsetting / allele-flip decision taken for ld_r below is applied to it
            # too, so both matrices always describe the same SNPs in the same order.
            ld_r_exp <- NULL
            if (ld_wgt_mode == "exposure" && !is.null(pop_wgt_exp)) {
                ld_r_exp <- tryCatch(
                    get_ld_matrix_gauss(gauss_in, chrI,
                                        min(gauss_in$bp, na.rm = TRUE),
                                        max(gauss_in$bp, na.rm = TRUE),
                                        pop_wgt = pop_wgt_exp, tag = "exp"),
                    error = function(e) { lg("  exposure computeLD error: %s", e$message); NULL })
                if (is.null(ld_r_exp))
                    lg("  exposure-weighted LD unavailable; using the outcome-weighted LD for both traits")
            }

            # GAUSS keys the matrix by the rsid we supplied
            common <- intersect(harm_data$rsid, rownames(ld_r))
            if (length(common) < 3) stop(sprintf("fewer than 3 SNPs shared with LD matrix (%d)", length(common)))
            harm_data <- harm_data[match(common, harm_data$rsid), ]
            ld_r <- ld_r[common, common, drop = FALSE]
            if (!is.null(ld_r_exp)) {
                if (is.null(rownames(ld_r_exp)) || !all(common %in% rownames(ld_r_exp))) {
                    lg("  exposure-weighted LD covers a different SNP set; using the outcome-weighted LD for both traits")
                    ld_r_exp <- NULL
                } else {
                    ld_r_exp <- ld_r_exp[common, common, drop = FALSE]
                }
            }

            # ---- ALLELE-ORIENTATION CHECK (signed R) ----
            # susie_rss needs sign(z) and sign(R) on the same allele coding. GAUSS
            # reports the a1/a2 it used per SNP; where that is the reverse of our
            # REF/ALT the SNP *may* need its LD row/column negated. Whether GAUSS's R
            # is already oriented to the supplied alleles is not documented, so the
            # flip is NOT applied blindly: the candidate flip vector is kept and, after
            # LD conditioning, estimate_s_rss lambda is computed with and without it
            # and the orientation with the lower lambda (better z/R consistency) wins.
            # In practice the un-flipped matrix has given lambda ~0.01, i.e. consistent.
            flip_vec <- rep(FALSE, nrow(harm_data))
            if (!is.null(gall) && all(c("rsid","a1","a2") %in% names(gall))) {
                gal_map <- gall[match(common, gall$rsid), ]
                our_ref <- toupper(harm_data$other_allele.exposure)   # REF in our coding
                our_alt <- toupper(harm_data$effect_allele.exposure)  # ALT (z sign allele)
                g_a1 <- toupper(gal_map$a1); g_a2 <- toupper(gal_map$a2)
                same <- (g_a1 == our_ref & g_a2 == our_alt)
                flip <- (g_a1 == our_alt & g_a2 == our_ref)
                flip[is.na(flip)] <- FALSE
                # SNPs whose alleles match neither orientation (multiallelic / strand)
                # are unreliable for signed R -> drop them.
                bad <- !(same | flip)
                bad[is.na(bad)] <- TRUE
                if (any(bad)) {
                    keepb <- !bad
                    if (sum(keepb) < 3) stop(sprintf("fewer than 3 allele-consistent SNPs (%d)", sum(keepb)))
                    ld_r <- ld_r[keepb, keepb, drop = FALSE]
                    if (!is.null(ld_r_exp)) ld_r_exp <- ld_r_exp[keepb, keepb, drop = FALSE]
                    harm_data <- harm_data[keepb, ]
                    flip <- flip[keepb]
                    lg("  dropped %d allele-inconsistent SNPs", sum(bad))
                }
                flip_vec <- flip
                lg("  GAUSS allele labels reversed vs ours for %d / %d SNPs (flip decided by lambda below)",
                   sum(flip_vec), length(flip_vec))
            } else {
                lg("  WARNING: GAUSS alleles unavailable; cannot test LD sign orientation")
            }

            # ---- LD matrix conditioning ----
            # Cells GAUSS cannot compute come back NA. They used to be replaced by the
            # column mean, which invents correlations that susie_rss then treats as
            # data; the affected SNPs are DROPPED instead. Columns carrying more than
            # 10 NAs go first (a SNP the panel barely covers), then any SNP that still
            # has a single NA. Removing those rows/columns cannot create new NAs, so
            # one pass is enough. (cov2cor was a no-op here - computeLD already
            # returns a correlation matrix - and has been removed.)
            ld_keep <- function(M, max_na = 10L) {
                storage.mode(M) <- "numeric"
                k1 <- apply(M, 2, function(x) sum(is.na(x)) <= max_na)
                if (sum(k1) < 3) return(k1 & FALSE)
                M2 <- M[k1, k1, drop = FALSE]
                diag(M2) <- 1                 # the diagonal says nothing about coverage
                k2 <- apply(M2, 2, function(x) !any(is.na(x)))
                k1[k1] <- k2
                k1
            }
            condition_ld <- function(M, keep_idx) {
                M <- M[keep_idx, keep_idx, drop = FALSE]
                storage.mode(M) <- "numeric"
                diag(M) <- 1
                # susie_rss requires a symmetric R; GAUSS output can carry tiny FP
                # asymmetries. Symmetrize defensively.
                if (!isSymmetric(unname(M))) M <- (M + t(M)) / 2
                if (anyNA(M)) stop(sprintf("LD matrix still holds %d NA cell(s) after dropping", sum(is.na(M))))
                M
            }

            keep_ld <- ld_keep(ld_r)
            # Both traits must end up on the same SNP set, so a SNP unusable in either
            # matrix is dropped from both.
            if (!is.null(ld_r_exp)) keep_ld <- keep_ld & ld_keep(ld_r_exp)
            if (sum(keep_ld) < 3)
                stop(sprintf("fewer than 3 SNPs with a complete LD row/column (%d)", sum(keep_ld)))
            if (any(!keep_ld))
                lg("  dropped %d SNP(s) with NA cells in the LD matrix (dropped, not imputed)", sum(!keep_ld))
            harm_data <- harm_data[keep_ld, ]
            flip_vec <- flip_vec[keep_ld]
            ld_r <- condition_ld(ld_r, keep_ld)
            if (!is.null(ld_r_exp)) {
                ld_r_exp <- tryCatch(condition_ld(ld_r_exp, keep_ld),
                    error = function(e) {
                        lg("  exposure LD conditioning failed (%s); using the outcome-weighted LD for both traits", e$message)
                        NULL
                    })
            }

            z_expo <- harm_data$beta.exposure / harm_data$se.exposure
            z_outc <- harm_data$beta.outcome  / harm_data$se.outcome

            # ---- Decide the allele-orientation flip empirically ----
            # Only when some SNPs are candidates: compare z/R consistency (lambda from
            # estimate_s_rss, exposure trait - the stronger signal) with and without
            # negating those SNPs' rows/columns, and keep whichever is more consistent.
            flip_applied <- FALSE
            if (any(flip_vec) && !all(flip_vec)) {   # flipping ALL SNPs is a no-op
                n_tmp <- suppressWarnings(as.numeric(NF_sample_size)); if (!is.finite(n_tmp) || n_tmp <= 0) n_tmp <- 500
                sgn <- ifelse(flip_vec, -1, 1)
                apply_flip <- function(M) sweep(sweep(M, 1, sgn, `*`), 2, sgn, `*`)
                ld_flip <- apply_flip(ld_r)
                lam_nf <- tryCatch(susieR:::estimate_s_rss(z_expo, ld_r,    n = n_tmp), error = function(e) NA_real_)
                lam_fl <- tryCatch(susieR:::estimate_s_rss(z_expo, ld_flip, n = n_tmp), error = function(e) NA_real_)
                lg("  orientation test (exposure lambda): as-is=%.4g  flipped(%d SNPs)=%.4g",
                   lam_nf, sum(flip_vec), lam_fl)
                if (is.finite(lam_fl) && is.finite(lam_nf) && lam_fl < 0.8 * lam_nf) {
                    ld_r <- ld_flip; flip_applied <- TRUE
                    # the same decision must hold for the exposure-weighted matrix
                    if (!is.null(ld_r_exp)) ld_r_exp <- apply_flip(ld_r_exp)
                    lg("  -> applied sign flip to %d SNPs", sum(flip_vec))
                } else {
                    lg("  -> kept GAUSS R as-is (no flip)")
                }
            } else if (all(flip_vec) && length(flip_vec) > 0) {
                lg("  GAUSS labels reversed for ALL SNPs: uniform relabel, no sign change needed")
            }
            # Exposure N: GTEx tissue sample size OF THE ALL-PAIRS RELEASE (the window
            # is built from those statistics), i.e. --assoc_sample_size when set, else
            # --sample_size. ma_samples counts only minor-allele carriers, so it
            # UNDERSTATES N; prefer the value passed in.
            n_expo <- suppressWarnings(as.numeric(NF_sample_size))
            if (!is.finite(n_expo) || n_expo <= 0) {
                n_expo <- suppressWarnings(max(harm_data$samplesize.exposure, na.rm = TRUE))
            }
            if (!is.finite(n_expo) || n_expo <= 0) stop("cannot determine exposure sample size")
            n_outc <- suppressWarnings(mean(harm_data$samplesize.outcome, na.rm=TRUE))
            if (!is.finite(n_outc) || n_outc <= 0) stop("cannot determine outcome sample size")
            lg("  N: exposure=%.0f outcome=%.0f", n_expo, n_outc)

            # Save the per-gene harmonised window (both traits, full cis-window with
            # p-values/positions) for downstream diagnostic locus/scatter plots. This
            # is the exact SNP set fed to coloc.susie, so the plots reflect what the
            # method actually saw. One file per gene: window_<ENSG>.rds
            tryCatch({
                win <- harm_data
                # ensure p-values present (derive from beta/se if needed)
                if (!"pval.exposure" %in% names(win) || all(is.na(win$pval.exposure)))
                    win$pval.exposure <- 2 * pnorm(-abs(win$beta.exposure / win$se.exposure))
                if (!"pval.outcome" %in% names(win) || all(is.na(win$pval.outcome)))
                    win$pval.outcome <- 2 * pnorm(-abs(win$beta.outcome / win$se.outcome))
                win$Gene_ID <- gene_id; win$Gene_Name <- prot; win$chromosome <- chrI
                saveRDS(win, sprintf("window_%s.rds", ens))
                # Also save the SIGNED-R LD matrix (rsid-named), aligned to this window,
                # so downstream methods (e.g. SharePro) reuse the exact same LD SuSiE used.
                ld_save <- ld_r
                if (is.null(rownames(ld_save))) { rownames(ld_save) <- harm_data$rsid; colnames(ld_save) <- harm_data$rsid }
                saveRDS(ld_save, sprintf("ld_%s.rds", ens))
                # ld_<ENSG>.rds always stays the OUTCOME-weighted matrix (SharePro and
                # the diagnostics read it by that name); the exposure-weighted one, when
                # it exists, is saved alongside it.
                if (!is.null(ld_r_exp)) {
                    ld_save_e <- ld_r_exp
                    if (is.null(rownames(ld_save_e))) { rownames(ld_save_e) <- harm_data$rsid; colnames(ld_save_e) <- harm_data$rsid }
                    saveRDS(ld_save_e, sprintf("ld_exp_%s.rds", ens))
                }
            }, error = function(e) lg("  (window/LD save skipped: %s)", e$message))

            # ---- SuSiE-RSS on the SIGNED multi-ethnic LD matrix ----
            # Only signed R is valid here: susie_rss models signed correlations, so a
            # squared-R run (previously done for comparison) violates the model - it
            # inflated estimate_s_rss lambda and produced degenerate H4==1 rows. It has
            # been removed. The LD_type column now says where a row came from:
            #   "susie"        - a real coloc.susie credible-set pair
            #   "abf_fallback" - coloc.abf on the same window, because coloc.susie
            #                    could not form a credible-set pair
            # R_out is the outcome-weighted LD; R_exp is the matrix the EXPOSURE trait
            # is fine-mapped on (the exposure-weighted one when
            # --susie_exposure_ld_weights = "exposure" and it could be built, else the
            # same matrix as the outcome).
            run_coloc <- function(R_out, R_exp = NULL) {
                if (is.null(R_exp)) R_exp <- R_out
                st <- new_status(gene_id, prot, "ok")
                st$n_snps <- nrow(harm_data)
                st$minp_exp <- suppressWarnings(min(2 * pnorm(-abs(z_expo)), na.rm = TRUE))
                st$minp_out <- suppressWarnings(min(2 * pnorm(-abs(z_outc)), na.rm = TRUE))
                st$flip_applied <- flip_applied
                st$window_center_snp <- center_snp; st$window_center_pos <- position
                lg("  min p in window: exposure=%.2e outcome=%.2e", st$minp_exp, st$minp_out)

                # ---- LD check (Zheng 2020 Nat Genet; Zuber 2022 review): max r2 between
                # any MR instrument and the 30 strongest GWAS SNPs of the window, on the
                # same signed-R matrix SuSiE uses (R_out is aligned to harm_data).
                ldc <- tryCatch(coloc_ld_check(R_out, snp_ids = as.character(harm_data$SNP),
                                               instrument_ids = instr_snps,
                                               outcome_p = 2 * pnorm(-abs(z_outc)),
                                               top_n = 30, r2_thr = NF_ld_check_r2),
                                error = function(e) { lg("  LD check failed: %s", e$message); NULL })
                if (!is.null(ldc)) {
                    st$ld_check_max_r2 <- ldc$ld_check_max_r2; st$ld_check_pass <- ldc$ld_check_pass
                    st$ld_check_n_instr <- ldc$ld_check_n_instr
                    lg("  LD check: max r2(instrument, top-30 GWAS SNPs) = %s (%d instrument(s) in window; pass at r2 >= %s: %s)",
                       format(ldc$ld_check_max_r2, digits = 3), ldc$ld_check_n_instr, NF_ld_check_r2, ldc$ld_check_pass)
                }

                # each trait's z/R consistency is measured on the matrix it will use
                lam_e <- tryCatch(susieR:::estimate_s_rss(z_expo, R_exp, n = n_expo), error = function(e) NA_real_)
                lam_o <- tryCatch(susieR:::estimate_s_rss(z_outc, R_out, n = n_outc), error = function(e) NA_real_)
                st$lambda_exp <- lam_e; st$lambda_out <- lam_o
                st$ld_mismatch_flag <- isTRUE(lam_e > lam_warn) || isTRUE(lam_o > lam_warn)
                lg("  estimate_s_rss lambda: exposure=%.4g outcome=%.4g", lam_e, lam_o)
                if (isTRUE(st$ld_mismatch_flag))
                    lg("  WARNING: lambda above %.3g - the LD panel may not match these summary statistics", lam_warn)

                # Build coloc datasets and fine-map with coloc::runsusie (the coloc-
                # recommended wrapper around susie_rss). Both traits are quantitative
                # (eQTL expression and a quantitative GWAS).
                #   exposure: sdY = --exposure_sdY (GTEx expression is inverse-normal
                #             transformed, so 1 is the right scale).
                #   outcome : sdY = --outcome_sdY when that is a number; otherwise coloc
                #             estimates it from the GWAS MAF + N (the MAF comes from the
                #             outcome allele-frequency column, flip-invariant).
                # The *_noLD datasets are the SAME window without the LD element, and
                # are what the coloc.abf fallback below runs on.
                rs <- harm_data$rsid
                D_exp_noLD <- list(beta = harm_data$beta.exposure,
                                   varbeta = harm_data$se.exposure^2,
                                   snp = rs, type = "quant", N = n_expo, sdY = sdY_exp)
                D_out_noLD <- list(beta = harm_data$beta.outcome,
                                   varbeta = harm_data$se.outcome^2,
                                   snp = rs, type = NF_outcome_type, N = n_outc)
                if (NF_outcome_type == "cc") {
                    # case/control GWAS: coloc needs the case fraction; sdY is not used
                    D_out_noLD$s <- NF_outcome_case_prop
                } else if (use_maf_n) {
                    maf_o <- unname(af_lookup[as.character(harm_data$SNP)])
                    D_out_noLD$MAF <- pmin(maf_o, 1 - maf_o)
                } else {
                    D_out_noLD$sdY <- sdY_out
                }
                D_exp <- c(D_exp_noLD, list(LD = R_exp))
                D_out <- c(D_out_noLD, list(LD = R_out))

                # coloc.abf on the same window, run whenever coloc.susie cannot produce
                # a credible-set pair, so a gene is never silently dropped. The row is
                # tagged LD_type="abf_fallback" and st$abf_fallback=TRUE so downstream
                # code knows it is NOT a SuSiE result.
                do_fallback <- function(st, why) {
                    if (!abf_fb) {
                        lg("  no rows for this gene (%s); ABF fallback disabled", why)
                        return(list(sm = NULL, st = st))
                    }
                    r <- tryCatch(coloc::coloc.abf(D_exp_noLD, D_out_noLD, p12 = NF_coloc_p12),
                                  error = function(e) { lg("  ABF fallback failed: %s", e$message); NULL })
                    if (is.null(r) || is.null(r$summary)) return(list(sm = NULL, st = st))
                    sy <- as.list(r$summary)
                    h1 <- NA_character_; h2 <- NA_character_
                    rr <- tryCatch(as.data.frame(r$results), error = function(e) NULL)
                    if (!is.null(rr) && "snp" %in% names(rr) && nrow(rr) > 0) {
                        if ("lABF.df1" %in% names(rr)) h1 <- as.character(rr$snp[which.max(rr$lABF.df1)])
                        if ("lABF.df2" %in% names(rr)) h2 <- as.character(rr$snp[which.max(rr$lABF.df2)])
                    }
                    fb <- data.frame(
                        nsnps = suppressWarnings(as.numeric(sy$nsnps)),
                        hit1 = h1, hit2 = h2,
                        PP.H0.abf = suppressWarnings(as.numeric(sy$PP.H0.abf)),
                        PP.H1.abf = suppressWarnings(as.numeric(sy$PP.H1.abf)),
                        PP.H2.abf = suppressWarnings(as.numeric(sy$PP.H2.abf)),
                        PP.H3.abf = suppressWarnings(as.numeric(sy$PP.H3.abf)),
                        PP.H4.abf = suppressWarnings(as.numeric(sy$PP.H4.abf)),
                        idx1 = NA_real_, idx2 = NA_real_, stringsAsFactors = FALSE)
                    fb <- fb[, std_cols, drop = FALSE]
                    fb$LD_type   <- "abf_fallback"
                    fb$Gene_Name <- as.character(prot); fb$Gene_ID <- as.character(gene_id)
                    # prior sensitivity: coloc.abf on the same window under the p12 grid
                    fb <- attach_p12_cols(fb, function(p) {
                        z <- NULL
                        utils::capture.output(z <- suppressMessages(coloc::coloc.abf(D_exp_noLD, D_out_noLD, p12 = p)))
                        as.data.frame(as.list(z$summary))
                    }, obj = r)
                    st$abf_fallback <- TRUE
                    lg("  ABF fallback (%s): nsnps=%.0f H4=%.3f hit1=%s hit2=%s",
                       why, fb$nsnps[1], fb$PP.H4.abf[1], fb$hit1[1], fb$hit2[1])
                    list(sm = fb, st = st)
                }

                # ---- Hybrid: SuSiE credible sets x single-signal ABF (coloc.susie_bf) ----
                # When exactly ONE trait has a credible set (almost always the eQTL), plain
                # coloc.abf discards the fine-mapping we do have. coloc::coloc.susie_bf
                # colocalises each SuSiE credible set of the resolvable trait against the
                # other trait's per-SNP log Bayes factors (one causal variant), so with
                # several eQTL signals it says WHICH one the (weak) GWAS signal sits on.
                # Rows are tagged LD_type = "susie_bf"; never presented as coloc.susie.
                do_hybrid <- function(S_have, D_other_noLD, have_is_exposure) {
                    if (!exists("coloc.susie_bf", where = asNamespace("coloc"), inherits = FALSE)) {
                        lg("  coloc.susie_bf not available in this coloc version; using ABF fallback")
                        return(NULL)
                    }
                    fm <- tryCatch(coloc::finemap.abf(D_other_noLD),
                                   error = function(e) { lg("  finemap.abf failed: %s", e$message); NULL })
                    if (is.null(fm)) return(NULL)
                    fm <- as.data.frame(fm)
                    lcol <- grep("^lABF", names(fm), value = TRUE)[1]
                    if (is.na(lcol) || !"snp" %in% names(fm)) { lg("  finemap.abf output lacks lABF/snp columns"); return(NULL) }
                    fm <- fm[!is.na(fm$snp) & fm$snp != "null", , drop = FALSE]
                    bf2 <- setNames(as.numeric(fm[[lcol]]), as.character(fm$snp))
                    bf2 <- bf2[is.finite(bf2)]
                    hy <- tryCatch(coloc::coloc.susie_bf(S_have, bf2, p12 = NF_coloc_p12),
                                   error = function(e) { lg("  coloc.susie_bf failed: %s", e$message); NULL })
                    if (is.null(hy) || is.null(hy$summary) || nrow(hy$summary) == 0) {
                        lg("  coloc.susie_bf returned no summary"); return(NULL)
                    }
                    hs <- as.data.frame(hy$summary)
                    for (cc in std_cols) if (!(cc %in% colnames(hs))) hs[[cc]] <- NA
                    hs <- hs[, std_cols, drop = FALSE]
                    if (!have_is_exposure) {
                        # trait order was (outcome, exposure): swap back to (exposure, outcome)
                        hs <- hs[, c("nsnps","hit2","hit1","PP.H0.abf","PP.H2.abf","PP.H1.abf",
                                     "PP.H3.abf","PP.H4.abf","idx2","idx1"), drop = FALSE]
                        names(hs) <- std_cols
                    }
                    for (cc in c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","idx1","idx2"))
                        hs[[cc]] <- suppressWarnings(as.numeric(hs[[cc]]))
                    hs$hit1 <- as.character(hs$hit1); hs$hit2 <- as.character(hs$hit2)
                    hs$LD_type   <- "susie_bf"
                    hs$Gene_Name <- as.character(prot); hs$Gene_ID <- as.character(gene_id)
                    # prior sensitivity: same credible sets x BFs under the p12 grid. Rows of a
                    # re-run come back in the same credible-set order, so key by row order
                    # (idx columns were swapped above when the trait order was reversed).
                    hs_key <- hs; hs_key$idx1 <- NA_real_; hs_key$idx2 <- NA_real_
                    hs2 <- attach_p12_cols(hs_key, function(p) {
                        z <- coloc::coloc.susie_bf(S_have, bf2, p12 = p)
                        if (is.null(z) || is.null(z$summary)) NULL else as.data.frame(z$summary)
                    }, obj = hy)
                    hs <- cbind(hs, hs2[, setdiff(names(hs2), names(hs)), drop = FALSE])
                    lg("  susie_bf (%s): %d credible set(s) x single-signal ABF, max H4=%.3f",
                       if (have_is_exposure) "SuSiE exposure x ABF outcome" else "ABF exposure x SuSiE outcome",
                       nrow(hs), max(hs$PP.H4.abf, na.rm = TRUE))
                    hs
                }

                S_exp <- tryCatch(coloc::runsusie(D_exp, suffix = "exp"),
                                  error = function(e) { lg("  runsusie(exposure) failed: %s", e$message); NULL })
                S_out <- tryCatch(coloc::runsusie(D_out, suffix = "out"),
                                  error = function(e) { lg("  runsusie(outcome) failed: %s", e$message); NULL })

                # Per-trait SuSiE diagnostics. Comparing the three credible-set counts
                # tells the failure modes apart:
                #   unpruned > n_cs   -> sets exist but are killed by the purity filter
                #                        (min_abs_corr = 0.5) under this LD matrix
                #   cov90    > n_cs   -> coverage sits just under 0.95
                #   all zero          -> no signal SuSiE can resolve at all
                fill_diag <- function(st, S, sfx, Rx = R_out) {
                    if (is.null(S)) return(st)
                    st[[paste0("n_cs_", sfx)]] <- tryCatch(
                        if (is.null(S$sets) || is.null(S$sets$cs)) 0L else length(S$sets$cs),
                        error = function(e) NA_integer_)
                    # NO Xcorr => susie_get_cs applies no purity filter
                    st[[paste0("n_cs_unpruned_", sfx)]] <- tryCatch({
                        u <- susieR::susie_get_cs(S, coverage = 0.95)
                        if (is.null(u$cs)) 0L else length(u$cs)
                    }, error = function(e) NA_integer_)
                    # looser coverage, purity filter still applied
                    st[[paste0("n_cs_cov90_", sfx)]] <- tryCatch({
                        u <- susieR::susie_get_cs(S, Xcorr = Rx, coverage = 0.90, min_abs_corr = 0.5)
                        if (is.null(u$cs)) 0L else length(u$cs)
                    }, error = function(e) NA_integer_)
                    st[[paste0("min_purity_", sfx)]] <- tryCatch({
                        pu <- S$sets$purity
                        if (!is.null(pu) && "min.abs.corr" %in% colnames(pu) && nrow(pu) > 0)
                            min(as.numeric(pu$min.abs.corr), na.rm = TRUE) else NA_real_
                    }, error = function(e) NA_real_)
                    st[[paste0("max_pip_", sfx)]] <- tryCatch(max(as.numeric(S$pip), na.rm = TRUE),
                                                              error = function(e) NA_real_)
                    st[[paste0("converged_", sfx)]] <- isTRUE(S$converged)
                    st
                }
                st <- fill_diag(st, S_exp, "exp", R_exp)
                st <- fill_diag(st, S_out, "out", R_out)
                lg("  credible sets exp/out: %s/%s | unpruned %s/%s | cov90 %s/%s | min purity %s/%s | max PIP %s/%s | converged %s/%s",
                   st$n_cs_exp, st$n_cs_out,
                   st$n_cs_unpruned_exp, st$n_cs_unpruned_out,
                   st$n_cs_cov90_exp, st$n_cs_cov90_out,
                   format(st$min_purity_exp, digits = 3), format(st$min_purity_out, digits = 3),
                   format(st$max_pip_exp, digits = 3), format(st$max_pip_out, digits = 3),
                   st$converged_exp, st$converged_out)

                if (is.null(S_exp) || is.null(S_out)) {
                    st$reason <- if (is.null(S_exp)) "runsusie_failed_exposure" else "runsusie_failed_outcome"
                    return(do_fallback(st, st$reason))
                }

                ncs_e <- if (is.na(st$n_cs_exp)) 0L else st$n_cs_exp
                ncs_o <- if (is.na(st$n_cs_out)) 0L else st$n_cs_out
                if (ncs_e == 0 || ncs_o == 0) {
                    st$reason <- if (ncs_e == 0 && ncs_o == 0) "no_cs_both" else
                                 if (ncs_e == 0) "no_cs_exposure" else "no_cs_outcome"
                    lg("  no credible-set pair possible (%s)", st$reason)
                    if (hyb_on && xor(ncs_e == 0, ncs_o == 0)) {
                        hs <- if (ncs_o == 0) do_hybrid(S_exp, D_out_noLD, TRUE) else do_hybrid(S_out, D_exp_noLD, FALSE)
                        if (!is.null(hs)) {
                            st$hybrid_susie_bf <- TRUE; st$n_susie_rows <- nrow(hs)
                            return(list(sm = hs, st = st))
                        }
                    }
                    return(do_fallback(st, st$reason))
                }

                cs <- tryCatch(coloc::coloc.susie(S_exp, S_out, p12 = NF_coloc_p12),
                               error = function(e) { lg("  coloc.susie failed: %s", e$message); NULL })
                # A NULL result, or one carrying no $summary, simply means no credible-set
                # pair came back - that is a reportable outcome, not an error.
                if (is.null(cs) || is.null(cs$summary) || nrow(cs$summary) == 0) {
                    lg("  no coloc.susie summary (no credible-set pair returned)")
                    st$reason <- "coloc_susie_failed"
                    return(do_fallback(st, st$reason))
                }
                sm <- as.data.frame(cs$summary)
                need <- c("hit1","hit2","PP.H4.abf")
                if (!all(need %in% colnames(sm))) {
                    lg("  coloc.susie summary missing columns: %s",
                       paste(setdiff(need, colnames(sm)), collapse = ","))
                    st$reason <- "coloc_susie_failed"
                    return(do_fallback(st, st$reason))
                }

                for (cc in std_cols) if (!(cc %in% colnames(sm))) sm[[cc]] <- NA
                sm <- sm[, std_cols, drop = FALSE]
                for (cc in c("nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf","idx1","idx2"))
                    sm[[cc]] <- suppressWarnings(as.numeric(sm[[cc]]))
                sm$hit1 <- as.character(sm$hit1); sm$hit2 <- as.character(sm$hit2)
                sm$LD_type   <- "susie"
                sm$Gene_Name <- as.character(prot); sm$Gene_ID <- as.character(gene_id)
                # prior sensitivity: coloc.susie re-run on the SAME SuSiE fits under the
                # p12 grid (cheap), rows keyed by credible-set pair (idx1, idx2)
                sm <- attach_p12_cols(sm, function(p) {
                    z <- coloc::coloc.susie(S_exp, S_out, p12 = p)
                    if (is.null(z) || is.null(z$summary)) NULL else as.data.frame(z$summary)
                }, obj = cs)
                saveRDS(cs, sprintf("coloc_susie_%s_%s.rds", chrI, prot))
                st$reason <- "ok"; st$n_susie_rows <- nrow(sm)
                lg("  OK: %d credible-set pair(s), max H4=%.3f", nrow(sm), max(sm$PP.H4.abf, na.rm = TRUE))
                list(sm = sm, st = st)
            }

            tryCatch(run_coloc(ld_r, ld_r_exp), error = function(e) {
                lg("  coloc.susie stage error: %s", e$message)
                list(sm = NULL, st = new_status(gene_id, prot, "coloc_susie_failed"))
            })
        }, error = function(e) {
            lg("  skip %s: %s", gene_id, e$message)
            if (grepl("unmapped in all-pairs release", e$message, fixed = TRUE))
                lg("  hint: provide --gtf_exposure/--gtf_assoc to enable name/coordinate bridging")
            # Pre-LD failures reuse the stop() text as the status reason.
            list(sm = NULL, st = new_status(gene_id, prot, e$message))
        })

        status_all <- rbind(status_all, res_gene$st)
        if (!is.null(res_gene$sm) && nrow(res_gene$sm) > 0) {
            # Row types share the same schema; rbindlist(fill = TRUE) guards the
            # per-run prior-sensitivity columns all the same.
            tab_res <- if (nrow(tab_res) == 0) res_gene$sm else
                       as.data.frame(data.table::rbindlist(list(tab_res, res_gene$sm), fill = TRUE))
        }
    }

    # ---- per-gene status table (ALWAYS written: one row per input gene) ----
    if (is.null(status_all) || nrow(status_all) == 0) status_all <- empty_status
    data.table::fwrite(status_all, "susie_gene_status.tsv", sep = "\t")
    if (nrow(status_all) > 0) {
        rt <- table(status_all$reason, useNA = "ifany")
        lg("reasons: %s", paste(sprintf("%s=%d", names(rt), as.integer(rt)), collapse = " "))
    }
    n_fb <- if (!is.null(tab_res) && nrow(tab_res) > 0)
                sum(tab_res$LD_type == "abf_fallback", na.rm = TRUE) else 0L
    lg("abf_fallback rows: %d", n_fb)

    if (is.null(tab_res) || nrow(tab_res) == 0) {
        lg("SuSiE produced no results."); empty_out(status_all); close(log_con); quit(save="no", status=0)
    }

    lead <- c("Gene_Name","Gene_ID","LD_type")
    tab_res <- tab_res[, c(lead, setdiff(colnames(tab_res), lead))]

    # ---- decision metrics + rule + tiered class, per credible-set PAIR ----------
    # (bin/coloc_decision.R; identical to the classic COLOCALIZATION process)
    #   H3_plus_H4, H4_cond = H4/(H3+H4), H4_H3_odds, log2_H4_H3, dominant_hyp,
    #   pass_rule (NF_coloc_rule), coloc_class (strong / conditional / inconclusive /
    #   underpowered / distinct_signals). The strong tier uses NF_susie_h4_threshold.
    tab_res <- coloc_derive_metrics(tab_res)
    tab_res$pass_rule   <- coloc_apply_rule(tab_res, NF_coloc_rule)
    tab_res$coloc_class <- coloc_classify(tab_res, h4_strong = NF_susie_h4_threshold,
                                          cond_thr = NF_coloc_cond_threshold,
                                          power_min = NF_coloc_power_min,
                                          h3_strong = NF_coloc_h3_strong)
    tab_res$coloc_rule  <- rep(NF_coloc_rule, nrow(tab_res))
    saveRDS(tab_res, "susie_coloc_results_summary.rds")
    data.table::fwrite(tab_res, "susie_coloc_results_summary.tsv", sep = "\t")

    # Passing rows = rows (credible-set pairs) that satisfy the rule, any LD_type
    # (LD_type says which method produced the row). A gene passes if ANY pair passes.
    passing <- tab_res[tab_res$pass_rule %in% TRUE, , drop = FALSE]
    saveRDS(passing, "susie_passing_genes.rds")

    # ---- gene-level summary: one row per gene = its best pair -----------------
    # Priority: a real coloc.susie pair ("susie") > hybrid ("susie_bf") > coloc.abf
    # fallback; within a type the pair with the highest PP.H4. The ratio of the BEST
    # pair is reported (never the max of H4_cond across pairs, which a weak pair can
    # inflate). any_pass / any_pass_susie tell whether any pair / any real
    # coloc.susie pair passed the rule.
    prio <- c(susie = 1L, susie_bf = 2L, abf_fallback = 3L)
    pr <- unname(prio[as.character(tab_res$LD_type)]); pr[is.na(pr)] <- 9L
    o <- order(tab_res$Gene_ID, pr, -suppressWarnings(as.numeric(tab_res$PP.H4.abf)))
    best <- tab_res[o, , drop = FALSE]
    best <- best[!duplicated(best$Gene_ID), , drop = FALSE]
    agg <- function(x, g, f) { r <- tapply(x, g, f); unname(r[as.character(best$Gene_ID)]) }
    gid <- as.character(tab_res$Gene_ID)
    is_susie <- tab_res$LD_type == "susie"
    best$n_pairs        <- as.integer(agg(rep(1L, nrow(tab_res)), gid, sum))
    best$n_pairs_pass   <- as.integer(agg(as.integer(tab_res$pass_rule %in% TRUE), gid, sum))
    best$any_pass       <- best$n_pairs_pass > 0
    best$any_pass_susie <- as.logical(agg(tab_res$pass_rule %in% TRUE & is_susie, gid, any))
    best$max_H4_susie   <- agg(ifelse(is_susie, suppressWarnings(as.numeric(tab_res$PP.H4.abf)), NA_real_), gid,
                               function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    best$max_H4_any     <- agg(suppressWarnings(as.numeric(tab_res$PP.H4.abf)), gid,
                               function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    st_cols <- intersect(c("Gene_ID","reason","lambda_exp","lambda_out","ld_mismatch_flag",
                           "n_cs_exp","n_cs_out","ld_check_max_r2","ld_check_pass"), colnames(status_all))
    st_one <- status_all[!duplicated(status_all$Gene_ID), st_cols, drop = FALSE]
    gs <- merge(best, st_one, by = "Gene_ID", all.x = TRUE, sort = FALSE)
    for (cc in setdiff(gene_summary_cols, colnames(gs))) gs[[cc]] <- NA
    gs <- gs[, gene_summary_cols, drop = FALSE]
    data.table::fwrite(gs, "susie_gene_summary.tsv", sep = "\t")

    susie_tab <- tab_res[tab_res$LD_type == "susie", , drop = FALSE]
    hy_tab    <- tab_res[tab_res$LD_type == "susie_bf", , drop = FALSE]
    fb_tab    <- tab_res[tab_res$LD_type == "abf_fallback", , drop = FALSE]
    lg("SuSiE coloc complete. Rule: %s", NF_coloc_rule)
    lg("  susie       : %d credible-set pair(s), %d passing the rule, %d with PP.H4 > %s",
       nrow(susie_tab), sum(susie_tab$pass_rule), sum(susie_tab$PP.H4.abf > NF_susie_h4_threshold, na.rm = TRUE), NF_susie_h4_threshold)
    lg("  susie_bf    : %d credible set x ABF row(s), %d passing the rule, %d with PP.H4 > %s",
       nrow(hy_tab), sum(hy_tab$pass_rule), sum(hy_tab$PP.H4.abf > NF_susie_h4_threshold, na.rm = TRUE), NF_susie_h4_threshold)
    lg("  abf_fallback: %d row(s), %d passing the rule, %d with PP.H4 > %s",
       nrow(fb_tab), sum(fb_tab$pass_rule), sum(fb_tab$PP.H4.abf > NF_susie_h4_threshold, na.rm = TRUE), NF_susie_h4_threshold)
    cls_levels <- c("colocalized_strong","colocalized_conditional","inconclusive","underpowered","distinct_signals")
    cls <- table(factor(gs$coloc_class, levels = cls_levels))
    lg("  genes: %d with rows; %d pass the rule (any pair), %d via a real coloc.susie pair; best-pair classes: %s",
       nrow(gs), sum(gs$any_pass %in% TRUE), sum(gs$any_pass_susie %in% TRUE),
       paste(sprintf("%s=%d", names(cls), as.integer(cls)), collapse = " "))
    if ("ld_check_pass" %in% colnames(gs))
        lg("  LD check (r2 >= %s): %d pass, %d fail, %d NA", NF_ld_check_r2,
           sum(gs$ld_check_pass %in% TRUE), sum(gs$ld_check_pass %in% FALSE), sum(is.na(gs$ld_check_pass)))
    if (nrow(passing) > 0) {
        for (r in seq_len(nrow(passing))) {
            lg("  PASS [%s]: %s (%s) H4=%.3f H4/(H3+H4)=%.3f H3+H4=%.3f class=%s hit1=%s hit2=%s",
               passing$LD_type[r], passing$Gene_Name[r], passing$Gene_ID[r],
               passing$PP.H4.abf[r], passing$H4_cond[r], passing$H3_plus_H4[r], passing$coloc_class[r],
               passing$hit1[r], passing$hit2[r])
        }
    }
    # Rewrite sessionInfo now that every namespace (gauss, rtracklayer, ...) is loaded.
    writeLines(capture.output(sessionInfo()), "sessionInfo_susie.txt")
    close(log_con)
    