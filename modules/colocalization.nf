process COLOCALIZATION {
    container params.r_container
    publishDir "${params.outdir}/${params.tissue}/colocalization", mode: 'copy'

    maxRetries 1
    // Was hard-coded 'ignore', which silently hid every failure of this process
    // (an empty coloc table then looked like "no colocalising genes"). Use the
    // shared knob instead; the main workflow defaults it to 'terminate'.
    errorStrategy { params.coloc_error_strategy ?: 'terminate' }

    input:
    path harmonized_files
    path significant_results
    path assoc_files              // GTEx all-pairs association file(s), parquet
    path outcome_raw
    path gene_conversion_table
    path gene_map                 // gene_release_map.rds from GENE_RELEASE_MAP
    val tissue
    val window_size
    val h4_threshold
    val coloc_p12
    val outcome_sdY               // numeric -> fixed sdY for the GWAS; null/NA -> estimate from MAF+N
    val exposure_sdY              // GTEx expression is inverse-normal transformed -> 1
    val outdir

    output:
    path "coloc_results_summary.rds", emit: coloc_table
    path "coloc_results_summary.tsv", emit: coloc_tsv
    path "coloc_h4_fdr_table.tsv", emit: fdr_table
    path "passing_coloc_genes.rds", emit: passing_genes
    path "coloc_genes_for_proxy.rds", emit: proxy_genes
    path "coloc_skipped_genes.tsv", optional: true, emit: skipped
    path "sessionInfo_coloc.txt", emit: sessioninfo

    script:
    """
    #!/usr/bin/env Rscript

    # Only the packages this process actually uses. (LDlinkR / reshape / gridExtra /
    # ieugwasr / httr / purrr / jsonlite / susieR and the whole-tidyverse attach were
    # loaded but never called; they only slowed startup and widened the container
    # requirements.)
    suppressMessages({
        library(coloc)
        library(arrow)
        library(data.table)
        library(dplyr)
        library(tidyr)
        library(stringr)
        library(reshape2)
        library(TwoSampleMR)
    })

    # Package provenance for this run. Written FIRST so the file exists even if the
    # analysis below stops early.
    writeLines(capture.output(sessionInfo()), "sessionInfo_coloc.txt")

    # ---------------------------------------------------------------
    # Shared decision metrics (bin/coloc_decision.R): H3_plus_H4, H4_cond =
    # H4/(H3+H4), H4_H3_odds, log2_H4_H3, dominant_hyp, the configurable rule
    # (params.coloc_rule, coloc::sensitivity() syntax) and the tiered class.
    # bin/ is on PATH inside the task (Nextflow mounts it), so fall back to
    # Sys.which() if the projectDir path is not visible in the container.
    # ---------------------------------------------------------------
    helper <- "${projectDir}/bin/coloc_decision.R"
    if (!file.exists(helper)) helper <- Sys.which("coloc_decision.R")
    if (!nzchar(helper) || !file.exists(helper)) stop("bin/coloc_decision.R not found (projectDir/bin must sit next to nextflow_mr_pipeline.nf)")
    source(helper)
    rule_str  <- "${params.coloc_rule}"
    p12_grid  <- coloc_parse_p12_grid("${params.coloc_p12_grid}")
    min_snps  <- suppressWarnings(as.integer("${params.coloc_min_snps}")); if (!is.finite(min_snps) || min_snps < 3) min_snps <- 3L
    cond_thr  <- ${params.coloc_cond_threshold}
    power_min <- ${params.coloc_power_min}
    h3_strong <- ${params.coloc_h3_strong}
    fdr_target <- ${params.coloc_fdr_target}
    win_center <- tolower(trimws("${params.coloc_window_center}"))
    if (!(win_center %in% c("lead_instrument", "first_instrument"))) {
        cat(sprintf("WARNING: unknown coloc_window_center '%s'; using 'lead_instrument'\\n", win_center)); win_center <- "lead_instrument"
    }
    cat(sprintf("Decision rule: %s | strong tier H4 > %s | conditional tier H4/(H3+H4) >= %s with H3+H4 >= %s | distinct H3 > %s\\n",
                rule_str, ${h4_threshold}, cond_thr, power_min, h3_strong))
    cat(sprintf("p12 grid for prior sensitivity: %s | min SNPs per window: %d | window centre: %s\\n",
                paste(p12_grid, collapse = ","), min_snps, win_center))

    # ---------------------------------------------------------------
    # Outcome trait type: "quant" (default) or "cc" (case/control GWAS with logOR
    # betas -> coloc needs the case fraction `s`; sdY is not used for cc traits).
    # ---------------------------------------------------------------
    out_type <- tolower(trimws("${params.outcome_type}"))
    if (!(out_type %in% c("quant", "cc"))) { cat(sprintf("WARNING: unknown outcome_type '%s'; using 'quant'\\n", out_type)); out_type <- "quant" }
    s_cc <- suppressWarnings(as.numeric("${params.outcome_case_prop ?: ''}"))
    if (out_type == "cc" && (!is.finite(s_cc) || s_cc <= 0 || s_cc >= 1))
        stop("--outcome_type cc requires --outcome_case_prop (fraction of cases, 0 < s < 1)")

    # ---------------------------------------------------------------
    # sdY handling
    #   exposure : GTEx expression is inverse-normal transformed, so the effect
    #              sizes are already on a unit-variance scale -> sdY = exposure_sdY.
    #   outcome  : a finite --outcome_sdY is used as given ("fixed:<v>"); otherwise
    #              coloc estimates sdY itself, which needs MAF + N ("MAF+N").
    # The mode actually used is logged and written to coloc_results_summary.rds
    # as the column outcome_sdY_mode.
    # ---------------------------------------------------------------
    sdY_exp <- suppressWarnings(as.numeric("${exposure_sdY}"))
    if (!is.finite(sdY_exp) || sdY_exp <= 0) {
        cat(sprintf("WARNING: exposure_sdY='%s' is not a positive number; using sdY = 1\\n", "${exposure_sdY}"))
        sdY_exp <- 1
    }
    sdY_out  <- suppressWarnings(as.numeric("${outcome_sdY}"))
    use_maf_n <- !is.finite(sdY_out) || sdY_out <= 0
    sdy_mode <- if (use_maf_n) "MAF+N" else sprintf("fixed:%g", sdY_out)
    if (out_type == "cc") {
        # a cc dataset carries `s` instead of sdY; the MAF+N sdY route is not used
        use_maf_n <- FALSE; sdy_mode <- sprintf("cc:s=%g", s_cc)
    }

    # ---------------------------------------------------------------
    # Exposure = GTEx all-pairs association files (PARQUET).
    # Opened as a lazy arrow Dataset so each gene reads only its own
    # rows (avoids loading the full multi-GB file into memory).
    # Columns (v11 all-pairs): gene_id, variant_id, tss_distance, af,
    #   ma_samples, ma_count, pval_nominal, slope, slope_se
    # ---------------------------------------------------------------
    assoc_paths <- list.files(".", pattern = "\\\\.parquet\$", full.names = TRUE)
    if (length(assoc_paths) == 0) stop("no .parquet association files staged")
    cat(sprintf("Association parquet files: %d\\n", length(assoc_paths)))
    expo_ds <- arrow::open_dataset(assoc_paths, format = "parquet")
    cat("Association columns:", paste(names(expo_ds), collapse=","), "\\n")

    gcol <- if ("gene_id" %in% names(expo_ds)) "gene_id" else "phenotype_id"

    # ---------------------------------------------------------------
    # Gene-ID bridge between the exposure (eGene) GENCODE release and the
    # all-pairs release, built by GENE_RELEASE_MAP. Keyed by version-stripped
    # exposure ENSG; the value is the EXACT versioned id used in the parquet.
    # ---------------------------------------------------------------
    gmap <- tryCatch(readRDS("${gene_map}"), error = function(e) NULL)
    if (is.null(gmap)) cat("WARNING: gene_release_map could not be read; every gene will look unmapped\\n")
    assoc_id_of <- if (!is.null(gmap)) setNames(as.character(gmap\$assoc_gene_id), as.character(gmap\$exp_core)) else character(0)
    map_method_of <- if (!is.null(gmap)) setNames(as.character(gmap\$map_method), as.character(gmap\$exp_core)) else character(0)
    assoc_name_of <- if (!is.null(gmap)) setNames(as.character(gmap\$assoc_gene_name), as.character(gmap\$exp_core)) else character(0)

    method_for <- function(ens) {
        if (ens %in% names(map_method_of)) map_method_of[[ens]] else NA_character_
    }

    # Genes skipped for any reason are recorded here and written to
    # coloc_skipped_genes.tsv so no gene disappears silently.
    skipped_genes <- data.frame(Gene_ID = character(), reason = character(), stringsAsFactors = FALSE)
    add_skip <- function(gid, why) {
        skipped_genes <<- rbind(skipped_genes,
            data.frame(Gene_ID = as.character(gid), reason = as.character(why),
                       stringsAsFactors = FALSE))
    }

    # Result schema, used for the empty case and as the rbind target. Rows are
    # appended with rbindlist(fill = TRUE), so the prior-sensitivity columns
    # (H4_p12_<tag>, added per run from params.coloc_p12_grid) need not be listed.
    empty_res <- function() {
        data.frame(Gene_Name = character(), Gene_ID = character(), nSNP = numeric(),
                   H0 = numeric(), H1 = numeric(), H2 = numeric(),
                   H3 = numeric(), H4 = numeric(),
                   hit_exp = character(), hit_out = character(),
                   top_snp_H4 = character(), top_snp_PP_H4 = numeric(), cs95_size_H4 = integer(),
                   minp_exp = numeric(), minp_out = numeric(),
                   window_center_snp = character(), window_center_pos = numeric(),
                   p12_min_pass = numeric(), p12_max_pass = numeric(),
                   outcome_sdY_mode = character(), stringsAsFactors = FALSE)
    }

    # Pull one gene's cis rows from the parquet dataset (exact match on the
    # versioned id the map resolved - no grepl on a version-stripped ENSG).
    read_gene_assoc <- function(ens) {
        aid <- if (ens %in% names(assoc_id_of)) assoc_id_of[[ens]] else NA_character_
        if (is.null(aid) || length(aid) == 0 || is.na(aid) || !nzchar(aid)) return(data.frame())
        mm <- method_for(ens)
        if (!is.na(mm) && mm != "id_match") {
            cat(sprintf("  gene id bridged across releases via %s: %s -> %s (%s)\\n",
                        mm, ens, aid,
                        if (ens %in% names(assoc_name_of) && !is.na(assoc_name_of[[ens]])) assoc_name_of[[ens]] else "NA"))
        }
        d <- expo_ds %>%
            filter(.data[[gcol]] == aid) %>%
            collect()
        as.data.frame(d)
    }

    outcome <- readRDS("${outcome_raw}")
    # Robust outcome column standardization (handles n / TotalSampleSize / etc.)
    if (!"se" %in% colnames(outcome) && "standard_error" %in% colnames(outcome)) outcome\$se <- outcome\$standard_error
    if (!"pval" %in% colnames(outcome) && "p_value" %in% colnames(outcome)) outcome\$pval <- outcome\$p_value
    if (!"sample_size" %in% colnames(outcome)) {
        for (alt in c("TotalSampleSize","n","N","SampleSize","n_complete_samples")) {
            if (alt %in% colnames(outcome)) { outcome\$sample_size <- outcome[[alt]]; break }
        }
    }
    if (!"sample_size" %in% colnames(outcome)) outcome\$sample_size <- NA_real_

    # ---- Outcome allele frequency (needed only when sdY is estimated) ----
    # Detected once; carried into every cis-window as eaf.outcome so harmonise_data
    # keeps it, and MAF = pmin(eaf, 1-eaf) is recovered per gene afterwards.
    out_af_col <- intersect(c("effect_allele_frequency","eaf","af","EAF","AF","freq"),
                            colnames(outcome))[1]
    outcome\$eaf_out_raw <- if (!is.na(out_af_col)) suppressWarnings(as.numeric(outcome[[out_af_col]])) else NA_real_
    if (use_maf_n && is.na(out_af_col)) {
        cat("WARNING: --outcome_sdY not given and the outcome GWAS has no allele-frequency column\\n")
        cat("         (looked for effect_allele_frequency/eaf/af/EAF/AF/freq); falling back to sdY = 1.\\n")
        use_maf_n <- FALSE; sdY_out <- 1; sdy_mode <- "fixed:1(no_af)"
    }
    cat(sprintf("Outcome sdY mode: %s%s\\n", sdy_mode,
                if (use_maf_n) sprintf(" (allele frequency column '%s')", out_af_col) else ""))
    cat(sprintf("Exposure sdY: %g\\n", sdY_exp))

    # Load gene-name conversion table (columns: id, name) for labeling coloc output
    gene_conv <- tryCatch(readRDS("${gene_conversion_table}"), error = function(e) NULL)
    map_gene_name <- function(ids) {
        if (is.null(gene_conv) || !all(c("id","name") %in% colnames(gene_conv))) {
            return(as.character(ids))
        }
        ens <- sub(".*?(ENSG[0-9]+).*", "\\\\1", as.character(ids))
        nm <- gene_conv\$name[match(ens, gene_conv\$id)]
        nm[is.na(nm)] <- as.character(ids)[is.na(nm)]
        nm
    }

    # Load and combine harmonized data from all chromosomes
    har_files <- list.files(".", pattern = "harmonized_.*\\\\.rds\$")
    har_total <- data.frame()

    for (i in seq_along(har_files)) {
        df <- readRDS(har_files[i])
        if (nrow(df) > 0) {
            chr_num <- as.numeric(gsub("harmonized_chr_(.*)\\\\.rds", "\\\\1", har_files[i]))
            df\$chromosome <- sprintf("%s", chr_num)
            har_total <- rbind(har_total, df)
        }
    }

    # Load significant results (genome-wide, MHC already excluded, FDR-significant)
    signif_total <- tryCatch({
        s <- readRDS("${significant_results}")
        if (is.null(s)) data.frame() else s
    }, error = function(e) data.frame())

    cat(sprintf("Harmonized SNPs: %d\\n", nrow(har_total)))
    cat(sprintf("Significant genes: %d\\n", length(unique(signif_total\$id.exposure))))

    tab_res <- empty_res()

    if (nrow(signif_total) == 0) {
        cat("No significant results found for colocalization\\n")
    } else {
        merged_df <- inner_join(har_total, signif_total,
                               by = c("id.exposure" = "id.exposure"))

        # -----------------------------------------------------------------
        # One gene's coloc.abf run.
        # Returns list(row = <one-row data.frame>) on success, or
        # list(skip = "<reason>") for any expected reason to drop the gene.
        # This used to be the body of the loop, where every skip was a `next`
        # INSIDE a tryCatch: R cannot `next` out of a tryCatch frame (it raises
        # "no loop for break/next"), so the handler caught those and recorded
        # them as analysis errors. With a function the skips are plain returns.
        # -----------------------------------------------------------------
        run_gene <- function(gene_id) {

            entry_all <- merged_df[merged_df\$id.exposure == gene_id, , drop = FALSE]
            if (nrow(entry_all) == 0) return(list(skip = "no harmonized entry"))
            # Window centre. v2.1 took the FIRST harmonized row, which is an arbitrary
            # instrument when a gene has several clumped instruments; the default is
            # now the LEAD instrument (smallest exposure p-value), so classic coloc
            # and SuSiE see the same, best-motivated window.
            if (win_center == "lead_instrument" && "pval.exposure" %in% names(entry_all)) {
                pe <- suppressWarnings(as.numeric(entry_all\$pval.exposure))
                entry <- if (any(is.finite(pe))) entry_all[which.min(pe), , drop = FALSE] else entry_all[1, , drop = FALSE]
            } else {
                entry <- entry_all[1, , drop = FALSE]
            }

            chr <- as.character(entry\$chromosome)
            if (is.na(chr) || !nzchar(chr)) return(list(skip = "no chromosome on the harmonized record"))

            position <- NA
            if ("position" %in% names(entry)) {
                position <- entry\$position
            } else if ("base_pair_location" %in% names(entry)) {
                position <- entry\$base_pair_location
            }
            position <- suppressWarnings(as.numeric(position))
            if (!is.finite(position)) return(list(skip = "no lead position"))
            center_snp <- if ("SNP" %in% names(entry)) as.character(entry\$SNP) else NA_character_

            pos_min <- position - ${window_size}
            pos_max <- position + ${window_size}
            chr_label <- paste0("chr", chr)

            # ---- Exposure (GTEx all-pairs, parquet): read this gene's cis rows ----
            ens_g <- sub(".*?(ENSG[0-9]+).*", "\\\\1", gene_id)
            mm_g <- method_for(ens_g)
            if (is.na(mm_g) || mm_g == "unmapped")
                return(list(skip = "unmapped in all-pairs release (see gene_release_map.tsv)"))

            expo_gene <- read_gene_assoc(ens_g)
            if (nrow(expo_gene) == 0) return(list(skip = "gene not found in association files"))

            vsplit <- reshape2::colsplit(expo_gene\$variant_id, "_",
                        c("v_chr","v_pos","v_ref","v_alt","v_build"))
            expo_gene\$v_chr <- as.character(vsplit\$v_chr)
            expo_gene\$v_pos <- suppressWarnings(as.numeric(vsplit\$v_pos))
            expo_gene\$v_ref <- toupper(as.character(vsplit\$v_ref))
            expo_gene\$v_alt <- toupper(as.character(vsplit\$v_alt))

            expo_subset <- expo_gene %>%
                filter(v_chr == chr_label, v_pos >= pos_min, v_pos <= pos_max) %>%
                transmute(
                    SNP = paste0(chr, ":", v_pos),
                    effect_allele.exposure = v_alt,
                    other_allele.exposure  = v_ref,
                    beta.exposure = slope,
                    se.exposure   = slope_se,
                    pval.exposure = pval_nominal,
                    id.exposure = gene_id,
                    exposure = gene_id
                )

            # ---- Outcome (GWAS): filter to cis-window + format for TwoSampleMR ----
            outcome_subset <- outcome %>%
                filter(chromosome == chr_label,
                       base_pair_location >= pos_min, base_pair_location <= pos_max) %>%
                transmute(
                    SNP = paste0(chr, ":", base_pair_location),
                    base_pair_location = base_pair_location,
                    effect_allele.outcome = toupper(effect_allele),
                    other_allele.outcome  = toupper(other_allele),
                    beta.outcome = beta,
                    se.outcome   = se,
                    pval.outcome = pval,
                    samplesize.outcome = sample_size,
                    eaf.outcome = eaf_out_raw,
                    id.outcome = "GWAS",
                    outcome = "GWAS"
                )

            if (nrow(expo_subset) == 0 || nrow(outcome_subset) == 0)
                return(list(skip = "empty cis-window"))

            # Drop palindromic SNPs + set dummy exposure eaf (same trick that fixed the
            # MR step; prevents the harmonise NA-eaf crash). harmonise_data with
            # action=2 needs eaf only for palindromic SNPs, and those are gone.
            is_pal <- function(a1, a2) (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")
            expo_subset <- expo_subset[!is_pal(expo_subset\$effect_allele.exposure, expo_subset\$other_allele.exposure), ]
            if (nrow(expo_subset) == 0) return(list(skip = "no non-palindromic exposure SNPs"))
            expo_subset\$eaf.exposure <- 0.3

            # Keep the REAL outcome frequency keyed by SNP: harmonise_data may flip
            # eaf.outcome to the exposure's allele coding, but MAF = pmin(f, 1-f) is
            # flip-invariant, so the raw value is all we need afterwards.
            af_lookup <- setNames(suppressWarnings(as.numeric(outcome_subset\$eaf.outcome)),
                                  as.character(outcome_subset\$SNP))
            na_eaf <- !is.finite(outcome_subset\$eaf.outcome)
            if (any(na_eaf)) outcome_subset\$eaf.outcome[na_eaf] <- 0.3

            harm_err <- NULL
            harm_data <- tryCatch({
                harmonise_data(expo_subset, outcome_subset, action = 2)
            }, error = function(e) { harm_err <<- e\$message; NULL })

            if (is.null(harm_data)) {
                cat(sprintf("  harmonise err: %s\\n", harm_err))
                return(list(skip = sprintf("harmonise error: %s", harm_err)))
            }
            if (nrow(harm_data) == 0) return(list(skip = "no harmonized SNPs"))

            # harmonise_data flags SNPs it could not align (mr_keep = FALSE); those
            # rows carry unusable betas and must not reach coloc.
            if ("mr_keep" %in% names(harm_data)) harm_data <- harm_data[harm_data\$mr_keep == TRUE, ]
            harm_data <- harm_data[!duplicated(harm_data\$SNP), ]

            if (nrow(harm_data) < min_snps)
                return(list(skip = sprintf("fewer than %d harmonized SNPs (%d)", min_snps, nrow(harm_data))))

            # ---- Outcome MAF + N (only when sdY is estimated) ----
            maf_out <- NULL; n_out <- NULL
            if (use_maf_n) {
                maf_out <- unname(af_lookup[as.character(harm_data\$SNP)])
                maf_out <- pmin(maf_out, 1 - maf_out)
                keep_maf <- is.finite(maf_out) & maf_out > 0 & maf_out < 1
                n_bad <- sum(!keep_maf)
                if (n_bad > 0) cat(sprintf("  dropped %d SNP(s) with NA/0/1 outcome MAF\\n", n_bad))
                # Both datasets must describe the same SNPs, so drop them everywhere.
                harm_data <- harm_data[keep_maf, ]
                maf_out <- maf_out[keep_maf]
                if (nrow(harm_data) < min_snps)
                    return(list(skip = sprintf("fewer than %d SNPs with usable outcome MAF (%d)", min_snps, nrow(harm_data))))
                n_vec <- suppressWarnings(as.numeric(harm_data\$samplesize.outcome))
                if (all(is.finite(n_vec)) && all(n_vec > 0)) {
                    n_out <- n_vec
                } else {
                    n_ok <- n_vec[is.finite(n_vec) & n_vec > 0]
                    n_out <- if (length(n_ok) > 0) mean(n_ok) else NA_real_
                }
                if (!all(is.finite(n_out)) || any(n_out <= 0))
                    return(list(skip = "no usable outcome sample size for sdY estimation"))
            }

            # Positions are optional for coloc.abf; only attach them when every SNP
            # has one (harmonise_data does not always carry base_pair_location).
            bp_vec <- if ("base_pair_location" %in% names(harm_data))
                          suppressWarnings(as.numeric(harm_data\$base_pair_location)) else rep(NA_real_, nrow(harm_data))
            if (all(is.na(bp_vec)))
                bp_vec <- suppressWarnings(as.numeric(sub("^.*:", "", as.character(harm_data\$SNP))))

            D1_expo <- list(
                type = "quant",
                snp = as.character(harm_data\$SNP),
                beta = as.numeric(harm_data\$beta.exposure),
                varbeta = as.numeric(harm_data\$se.exposure)^2,
                sdY = sdY_exp
            )

            D2_outc <- list(
                type = out_type,
                snp = as.character(harm_data\$SNP),
                beta = as.numeric(harm_data\$beta.outcome),
                varbeta = as.numeric(harm_data\$se.outcome)^2
            )
            n_gwas <- suppressWarnings(mean(as.numeric(harm_data\$samplesize.outcome), na.rm = TRUE))
            if (out_type == "cc") {
                # case/control: coloc needs the case fraction; sdY is not used
                D2_outc\$s <- s_cc
                if (is.finite(n_gwas) && n_gwas > 0) D2_outc\$N <- n_gwas
            } else if (use_maf_n) {
                # coloc estimates sdY itself from MAF + N.
                D2_outc\$MAF <- as.numeric(maf_out)
                D2_outc\$N   <- as.numeric(n_out)
            } else {
                D2_outc\$sdY <- sdY_out
            }
            if (all(is.finite(bp_vec))) {
                D1_expo\$position <- bp_vec
                D2_outc\$position <- bp_vec
            }

            res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc, p12 = ${coloc_p12})

            # Lead SNP per trait, most likely shared SNP and 95% credible set under H4;
            # min p per trait (an H1 window = strong eQTL, no GWAS signal);
            # H4 / H4_cond under the other priors + the p12 range over which the rule holds.
            hits <- coloc_abf_hits(res)
            prof <- coloc_p12_profile(D1_expo, D2_outc, rule = rule_str, p12_grid = p12_grid, res = res)
            minp_exp <- suppressWarnings(min(2 * pnorm(-abs(D1_expo\$beta / sqrt(D1_expo\$varbeta))), na.rm = TRUE))
            minp_out <- suppressWarnings(min(2 * pnorm(-abs(D2_outc\$beta / sqrt(D2_outc\$varbeta))), na.rm = TRUE))

            # Prefer an explicit geneName on the record; otherwise map ENSG -> symbol
            gene_name <- if ("geneName" %in% names(entry) && !is.na(entry\$geneName)) {
                as.character(entry\$geneName)
            } else {
                as.character(map_gene_name(gene_id))
            }

            sy <- res\$summary
            nmv <- names(sy)
            pick <- function(want, idx) {
                if (!is.null(nmv) && want %in% nmv) {
                    return(suppressWarnings(as.numeric(sy[[want]])))
                }
                if (length(sy) >= idx) {
                    return(suppressWarnings(as.numeric(sy[[idx]])))
                }
                NA_real_
            }

            list(row = cbind(
                data.frame(
                    Gene_Name = gene_name,
                    Gene_ID   = as.character(gene_id),
                    nSNP = pick("nsnps", 1),
                    H0 = pick("PP.H0.abf", 2), H1 = pick("PP.H1.abf", 3),
                    H2 = pick("PP.H2.abf", 4), H3 = pick("PP.H3.abf", 5),
                    H4 = pick("PP.H4.abf", 6),
                    stringsAsFactors = FALSE),
                hits,
                data.frame(minp_exp = minp_exp, minp_out = minp_out,
                           window_center_snp = center_snp, window_center_pos = position,
                           stringsAsFactors = FALSE),
                prof,
                data.frame(outcome_sdY_mode = sdy_mode, stringsAsFactors = FALSE)))
        }

        # Test EVERY significant gene (an earlier 100-gene cap silently truncated the set).
        all_genes <- unique(merged_df\$id.exposure)
        cat(sprintf("Genes to test with coloc.abf: %d\\n", length(all_genes)))

        for (i in seq_along(all_genes)) {
            gene_id <- all_genes[i]
            cat(sprintf("Processing gene %d/%d: %s\\n", i, length(all_genes), gene_id))

            r <- tryCatch(run_gene(gene_id), error = function(e) {
                cat(sprintf("Error processing gene %d (%s): %s\\n", i, gene_id, e\$message))
                list(skip = sprintf("error: %s", e\$message))
            })

            if (!is.null(r\$skip)) {
                cat(sprintf("  skip %s: %s\\n", gene_id, r\$skip))
                add_skip(gene_id, r\$skip)
                next
            }
            # rbindlist(fill = TRUE): rows carry the per-run H4_p12_<tag> columns
            if (!is.null(r\$row)) tab_res <- as.data.frame(data.table::rbindlist(list(tab_res, r\$row), fill = TRUE))
        }
    }

    # ---------------------------------------------------------------
    # Decision metrics + rule + tiered class (bin/coloc_decision.R).
    #   pass_rule   : params.coloc_rule evaluated on the row -> passing_coloc_genes.rds
    #   coloc_class : colocalized_strong / colocalized_conditional / inconclusive /
    #                 underpowered / distinct_signals
    # ---------------------------------------------------------------
    tab_res <- coloc_derive_metrics(tab_res)
    tab_res\$pass_rule   <- coloc_apply_rule(tab_res, rule_str)
    tab_res\$coloc_class <- coloc_classify(tab_res, h4_strong = ${h4_threshold}, cond_thr = cond_thr,
                                          power_min = power_min, h3_strong = h3_strong)
    tab_res\$coloc_rule  <- rep(rule_str, nrow(tab_res))
    saveRDS(tab_res, "coloc_results_summary.rds")
    data.table::fwrite(tab_res, "coloc_results_summary.tsv", sep = "\\t")

    passing <- tab_res[tab_res\$pass_rule %in% TRUE, , drop = FALSE]
    failing <- tab_res[!(tab_res\$pass_rule %in% TRUE), , drop = FALSE]

    saveRDS(passing, "passing_coloc_genes.rds")
    saveRDS(failing, "coloc_genes_for_proxy.rds")

    # Reales et al. 2026 (PLoS Genet): if PP.H4 is calibrated, the expected FDR of the
    # calls with H4 > alpha is mean(1 - H4) over those calls. One row per alpha, plus
    # the smallest alpha whose expected FDR is below params.coloc_fdr_target.
    fdr_tab <- coloc_fdr_table(tab_res\$H4)
    fdr_tab\$fdr_target <- fdr_target
    fdr_tab\$alpha_for_target_FDR <- coloc_alpha_for_fdr(tab_res\$H4, fdr_target)
    data.table::fwrite(fdr_tab, "coloc_h4_fdr_table.tsv", sep = "\\t")

    # Always write the skip table (possibly empty) so dropped genes are visible.
    data.table::fwrite(skipped_genes[, c("Gene_ID","reason"), drop = FALSE], "coloc_skipped_genes.tsv", sep = "\\t")
    if (nrow(skipped_genes) > 0) {
        cat(sprintf("Genes skipped: %d\\n", nrow(skipped_genes)))
        rt <- table(skipped_genes\$reason)
        for (nmr in names(rt)) cat(sprintf("  %s: %d\\n", nmr, as.integer(rt[[nmr]])))
    }

    cls_levels <- c("colocalized_strong","colocalized_conditional","inconclusive","underpowered","distinct_signals")
    cls <- table(factor(tab_res\$coloc_class, levels = cls_levels))
    cat(sprintf("Colocalization complete: %d genes tested; rule '%s' passed by %d; strong (H4 > %s): %d; classes: %s; outcome sdY mode %s\\n",
                nrow(tab_res), rule_str, nrow(passing), ${h4_threshold},
                sum(tab_res\$H4 > ${h4_threshold}, na.rm = TRUE),
                paste(sprintf("%s=%d", names(cls), as.integer(cls)), collapse = " "), sdy_mode))
    if (is.finite(fdr_tab\$alpha_for_target_FDR[1]))
        cat(sprintf("Expected-FDR calibration (Reales 2026): H4 > %.2f keeps the expected FDR below %g\\n",
                    fdr_tab\$alpha_for_target_FDR[1], fdr_target))
    """
}
