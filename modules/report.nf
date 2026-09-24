process GENERATE_REPORT {
    container params.r_container
    publishDir "${params.outdir}/${params.tissue}", mode: 'copy'

    input:
    path 'exposure_clumped_*.rds'
    path missing_counts
    path proxy_merge_counts
    path mr_threshold_info
    path gene_map_tsv
    path coloc_table              // coloc_results_summary.rds (classic, with v2.2 metrics)
    path coloc_fdr_table          // coloc_h4_fdr_table.tsv (Reales 2026 expected-FDR calibration)
    path susie_gene_summary       // susie_gene_summary.tsv, or the NO_FILE_SUSIE_GS placeholder
    path pwcoco_table             // pwcoco_results_summary.rds, or the NO_FILE_PWCOCO placeholder
    val tissue
    val fstat_threshold
    val outdir

    output:
    path "analysis_report_${tissue}.txt", emit: report

    script:
    """
    #!/usr/bin/env Rscript

    safe_read <- function(f) tryCatch(readRDS(f), error = function(e) NULL)

    miss  <- safe_read("${missing_counts}")
    prox  <- safe_read("${proxy_merge_counts}")
    mrinf <- safe_read("${mr_threshold_info}")

    getv <- function(x, field, default = NA) {
        if (is.null(x)) return(default)
        if (is.list(x) && field %in% names(x)) return(x[[field]])
        if (is.numeric(x) && length(x) == 1 && is.na(field)) return(x)
        default
    }

    # ---- Exposure instruments: total rows across the clumped per-chromosome files ----
    # These are the instruments actually carried into MR.
    clumped_files <- list.files(".", pattern = "^exposure_clumped_.*\\\\.rds\$")
    n_exposure <- NA
    if (length(clumped_files) > 0) {
        n_exposure <- 0L
        for (cf in clumped_files) {
            d <- safe_read(cf)
            if (!is.null(d) && is.data.frame(d)) n_exposure <- n_exposure + nrow(d)
        }
    }
    cat(sprintf("Exposure instruments counted from %d clumped file(s)\\n", length(clumped_files)))

    n_outcome  <- getv(miss, "n_outcome_snps")
    n_missing  <- getv(miss, "n_missing_snps")
    n_proxy    <- getv(prox, "n_proxies_merged")

    fdr_alpha  <- getv(mrinf, "fdr_alpha")
    n_mr       <- getv(mrinf, "n_mr_total")     # MHC already excluded upstream
    n_sig_fdr  <- getv(mrinf, "n_sig_fdr")
    n_primary  <- getv(mrinf, "n_primary")
    prim_multi <- getv(mrinf, "mr_primary_multi", "Inverse variance weighted")
    n_st_fail  <- getv(mrinf, "n_steiger_failed")
    n_st_na    <- getv(mrinf, "n_steiger_na")
    st_filter  <- getv(mrinf, "steiger_filter", FALSE)
    mhc_build  <- getv(mrinf, "mhc_build", "GRCh38")
    mhc_chr    <- getv(mrinf, "mhc_chr", "6")
    mhc_start  <- getv(mrinf, "mhc_start", 28510120)
    mhc_end    <- getv(mrinf, "mhc_end", 33480577)

    # Gene-ID release bridging summary (from GENE_RELEASE_MAP); tolerate absence.
    gmap <- tryCatch(read.delim("${gene_map_tsv}", stringsAsFactors = FALSE),
                     error = function(e) NULL)
    gm_count <- function(m) {
        if (is.null(gmap) || !("map_method" %in% colnames(gmap))) return(NA)
        sum(gmap\$map_method == m, na.rm = TRUE)
    }
    n_id_match  <- gm_count("id_match")
    n_name_map  <- gm_count("name_match")
    n_coord_map <- gm_count("coord_overlap")
    n_unmapped  <- gm_count("unmapped")
    n_bridged   <- if (is.na(n_name_map) || is.na(n_coord_map)) NA else n_name_map + n_coord_map
    n_type_mm   <- NA
    if (!is.null(gmap) && ("type_mismatch" %in% colnames(gmap)))
        n_type_mm <- sum(as.logical(gmap\$type_mismatch), na.rm = TRUE)

    fmt <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "NA" else format(x, scientific = FALSE, big.mark = "")
    fmt_g <- function(x) if (is.null(x) || length(x) == 0 || is.na(x)) "NA" else formatC(x, format = "g", digits = 3)

    steiger_txt <- if (isTRUE(as.logical(st_filter))) "removed before FDR (steiger_filter = true)" else "kept and reported (steiger_filter = false)"

    # Package versions available in this container (r_container)
    pv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not available")

    # ---- Colocalization (v2.2 decision metrics) ----
    cls_levels <- c("colocalized_strong","colocalized_conditional","inconclusive","underpowered","distinct_signals")
    class_counts <- function(cls) { tb <- table(factor(cls, levels = cls_levels)); paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " ") }
    coloc <- safe_read("${coloc_table}")
    if (is.null(coloc) || !is.data.frame(coloc)) coloc <- data.frame()
    n_coloc <- nrow(coloc)
    n_rule  <- if (n_coloc && "pass_rule" %in% colnames(coloc)) sum(coloc\$pass_rule %in% TRUE) else NA
    n_strong <- if (n_coloc && "H4" %in% colnames(coloc)) sum(coloc\$H4 > ${params.coloc_h4_threshold}, na.rm = TRUE) else NA
    n_cond  <- if (n_coloc && all(c("H4_cond","H3_plus_H4") %in% colnames(coloc)))
                   sum(coloc\$H4_cond >= ${params.coloc_cond_threshold} & coloc\$H3_plus_H4 >= ${params.coloc_power_min}, na.rm = TRUE) else NA
    cls_txt <- if (n_coloc && "coloc_class" %in% colnames(coloc)) class_counts(coloc\$coloc_class) else "NA"
    fdr_tab <- tryCatch(read.delim("${coloc_fdr_table}", stringsAsFactors = FALSE), error = function(e) NULL)
    alpha_fdr <- if (!is.null(fdr_tab) && "alpha_for_target_FDR" %in% colnames(fdr_tab) && nrow(fdr_tab) > 0) fdr_tab\$alpha_for_target_FDR[1] else NA
    fdr_at_thr <- if (!is.null(fdr_tab) && all(c("alpha","expected_FDR") %in% colnames(fdr_tab))) {
        i <- which(abs(fdr_tab\$alpha - ${params.coloc_h4_threshold}) < 1e-9); if (length(i)) fdr_tab\$expected_FDR[i[1]] else NA } else NA
    susie_gs_path <- "${susie_gene_summary}"
    susie_lines <- "SuSiE colocalization: not run (run_susie = false)"
    if (!grepl("^NO_FILE", basename(susie_gs_path))) {
        gs <- tryCatch(read.delim(susie_gs_path, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(gs) && nrow(gs) > 0) {
            susie_lines <- c(
                sprintf("SuSiE colocalization: %s gene(s) with rows; %s pass the rule (any pair); %s via a real coloc.susie pair; best-pair classes: %s",
                        fmt(nrow(gs)), fmt(sum(gs\$any_pass %in% TRUE)), fmt(sum(gs\$any_pass_susie %in% TRUE)), class_counts(gs\$coloc_class)),
                sprintf("LD check (max r2 between an MR instrument and the top-30 GWAS SNPs >= %s): pass=%s fail=%s NA=%s",
                        "${params.ld_check_r2}", fmt(sum(gs\$ld_check_pass %in% TRUE)), fmt(sum(gs\$ld_check_pass %in% FALSE)), fmt(sum(is.na(gs\$ld_check_pass)))))
        } else {
            susie_lines <- "SuSiE colocalization: ran, but produced no rows (see susie_gene_status.tsv)"
        }
    }
    # PWCoCo (pair-wise conditional + coloc.abf): best signal pair per gene, same rule/classes
    pw_path <- "${pwcoco_table}"
    pwcoco_line <- "PWCoCo: not run (no --pwcoco_bfile, or run_susie = false)"
    if (!grepl("^NO_FILE", basename(pw_path))) {
        pw <- safe_read(pw_path)
        if (!is.null(pw) && is.data.frame(pw) && nrow(pw) > 0 && "H4" %in% colnames(pw)) {
            pw\$H4 <- suppressWarnings(as.numeric(pw\$H4))
            ok <- pw[is.finite(pw\$H4), , drop = FALSE]
            if (nrow(ok) > 0) {
                o <- order(ok\$Gene_ID, -ok\$H4); best <- ok[o, , drop = FALSE]; best <- best[!duplicated(best\$Gene_ID), , drop = FALSE]
                n_any <- if ("pass_rule" %in% colnames(ok)) length(unique(ok\$Gene_ID[ok\$pass_rule %in% TRUE])) else NA
                pwcoco_line <- sprintf("PWCoCo: %s gene(s) scored ; %s with best-pair H4 > %s ; %s pass the rule (any signal pair) ; best-pair classes: %s ; status: %s",
                                       fmt(nrow(best)), fmt(sum(best\$H4 > ${params.coloc_h4_threshold}, na.rm = TRUE)), "${params.coloc_h4_threshold}",
                                       fmt(n_any), if ("coloc_class" %in% colnames(best)) class_counts(best\$coloc_class) else "NA",
                                       if ("pwcoco_status" %in% colnames(pw)) { tb <- table(pw\$pwcoco_status[!duplicated(pw\$Gene_ID)]); paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " ") } else "NA")
            } else {
                pwcoco_line <- sprintf("PWCoCo: ran but scored no gene (statuses: %s)",
                                       if ("pwcoco_status" %in% colnames(pw)) paste(unique(pw\$pwcoco_status), collapse = ",") else "NA")
            }
        }
    }

    report_lines <- c(
        sprintf("Significance method: FDR (Benjamini-Hochberg), q < %s", fmt_g(fdr_alpha)),
        sprintf("F-statistics threshold: %s", "${fstat_threshold}"),
        sprintf("MHC region (%s, chr%s:%s-%s): excluded from all calculations",
                as.character(mhc_build), as.character(mhc_chr), fmt(mhc_start), fmt(mhc_end)),
        sprintf("Number of variants present in the exposure: %s", fmt(n_exposure)),
        sprintf("Number of variants present in the outcome: %s", fmt(n_outcome)),
        sprintf("Number of variants present in the exposure but missing from the outcome: %s", fmt(n_missing)),
        sprintf("Number of variants in the proxy search that meet criteria to merge with the outcome: %s", fmt(n_proxy)),
        sprintf("Number of MR results (MHC excluded): %s", fmt(n_mr)),
        sprintf("Primary MR tests carrying the multiple-testing burden (Wald ratio if 1 SNP, else '%s'): %s",
                as.character(prim_multi), fmt(n_primary)),
        sprintf("Genes failing Steiger directionality (exposure -> outcome): %s - %s", fmt(n_st_fail), steiger_txt),
        sprintf("Genes with undetermined Steiger direction (NA): %s", fmt(n_st_na)),
        sprintf("Number of MR results after FDR significant threshold (q < %s): %s", fmt_g(fdr_alpha), fmt(n_sig_fdr)),
        sprintf("Gene-ID release bridging: id_match=%s name_match=%s coord_overlap=%s unmapped=%s",
                fmt(n_id_match), fmt(n_name_map), fmt(n_coord_map), fmt(n_unmapped)),
        sprintf("Genes bridged across GENCODE releases (name or coordinate): %s", fmt(n_bridged)),
        sprintf("Bridged genes with a gene_type mismatch: %s", fmt(n_type_mm)),
        sprintf("Colocalization decision rule (params.coloc_rule): %s ; coloc p12 = %s ; window centre: %s ; outcome type: %s",
                "${params.coloc_rule}", "${params.coloc_p12}", "${params.coloc_window_center}", "${params.outcome_type}"),
        sprintf("Classic coloc.abf: %s gene(s) tested ; %s pass the rule ; %s with H4 > %s (strong tier) ; %s with H4/(H3+H4) >= %s and H3+H4 >= %s (conditional tier)",
                fmt(n_coloc), fmt(n_rule), fmt(n_strong), "${params.coloc_h4_threshold}", fmt(n_cond),
                "${params.coloc_cond_threshold}", "${params.coloc_power_min}"),
        sprintf("Classic coloc.abf classes (strong / conditional / inconclusive / underpowered [H0-H2 dominate] / distinct_signals [H3 > %s]): %s",
                "${params.coloc_h3_strong}", cls_txt),
        sprintf("Expected FDR of the H4 > %s call set (Reales 2026, mean of 1 - H4): %s ; smallest H4 cut with expected FDR < %s: %s",
                "${params.coloc_h4_threshold}", fmt_g(fdr_at_thr), "${params.coloc_fdr_target}", fmt_g(alpha_fdr)),
        susie_lines,
        pwcoco_line,
        sprintf("Software: R %s ; coloc %s ; TwoSampleMR %s ; susieR %s",
                as.character(getRversion()), pv("coloc"), pv("TwoSampleMR"), pv("susieR"))
    )

    writeLines(report_lines, "analysis_report_${tissue}.txt")
    cat(paste(report_lines, collapse = "\\n"), "\\n")
    """
}
