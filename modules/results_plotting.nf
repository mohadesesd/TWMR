process RESULTS_PLOTTING {
    container params.r_container
    publishDir "${params.outdir}/${params.tissue}/plots", mode: 'copy'

    input:
    path 'mr_results_*.rds'
    path coloc_results
    path susie_table                 // susie_coloc_results_summary.rds, or a NO_FILE* placeholder
    path gene_conversion_table
    path threshold_info
    val tissue
    val outdir

    output:
    path "*.pdf", optional: true, emit: plots
    path "*.png", optional: true, emit: png_plots
    path "*_summary_stats.txt", emit: summary_stats
    path "significant_results_with_names.rds", optional: true, emit: named_results
    path "significant_results_with_names.csv", optional: true, emit: named_results_csv

    script:
    """
    #!/usr/bin/env Rscript

    library(ggplot2)
    library(dplyr)
    library(gridExtra)
    library(data.table)
    library(tidyverse)

    # Load all MR results
    mr_files <- list.files(".", pattern = "mr_results_.*\\\\.rds\$")
    mr_list <- list()
    for (file in mr_files) {
        df <- tryCatch(readRDS(file), error = function(e) NULL)
        if (!is.null(df) && nrow(df) > 0) mr_list[[length(mr_list) + 1]] <- df
    }
    mr_combined <- if (length(mr_list) > 0) as.data.frame(data.table::rbindlist(mr_list, fill = TRUE)) else data.frame()

    # ------------------------------------------------------------------
    # Only the PRIMARY test of each gene is plotted and counted: one test per
    # gene (Wald ratio when nsnp == 1, otherwise params.mr_primary_multi). The
    # other methods are sensitivity analyses and carry pval_fdr = NA, because the
    # BH correction in COMBINE_MR_RESULTS is applied to the primary rows only.
    # ------------------------------------------------------------------
    if (nrow(mr_combined) > 0 && "primary" %in% colnames(mr_combined)) {
        n_all_rows <- nrow(mr_combined)
        mr_combined <- mr_combined[which(as.logical(mr_combined\$primary)), , drop = FALSE]
        cat(sprintf("Restricted to primary MR rows: %d of %d\\n", nrow(mr_combined), n_all_rows))
    } else if (nrow(mr_combined) > 0) {
        cat("No 'primary' column in the MR results; plotting all rows\\n")
    }

    # Shared decision metrics (bin/coloc_decision.R): derive the v2.2 columns when a
    # table lacks them, and the Reales 2026 expected-FDR calibration of the H4 cut.
    helper <- "${projectDir}/bin/coloc_decision.R"
    if (!file.exists(helper)) helper <- Sys.which("coloc_decision.R")
    if (nzchar(helper) && file.exists(helper)) source(helper) else cat("WARNING: bin/coloc_decision.R not found; class summaries skipped\\n")
    have_helper <- exists("coloc_derive_metrics")
    rule_str  <- "${params.coloc_rule}"
    h4_thr    <- ${params.coloc_h4_threshold}
    susie_thr <- ${params.susie_h4_threshold}
    cond_thr  <- ${params.coloc_cond_threshold}
    power_min <- ${params.coloc_power_min}
    h3_strong <- ${params.coloc_h3_strong}
    fdr_target <- ${params.coloc_fdr_target}
    cls_levels <- c("colocalized_strong","colocalized_conditional","inconclusive","underpowered","distinct_signals")
    cls_colors <- c(colocalized_strong = "#1B7837", colocalized_conditional = "#7FBF7B", inconclusive = "#BDBDBD",
                    underpowered = "#92C5DE", distinct_signals = "#D6604D")
    complete_metrics <- function(d, strong_thr) {
        if (!have_helper || nrow(d) == 0) return(d)
        if (!"H4_cond" %in% names(d)) d <- coloc_derive_metrics(d)
        if (!"pass_rule" %in% names(d)) d\$pass_rule <- coloc_apply_rule(d, rule_str)
        if (!"coloc_class" %in% names(d)) d\$coloc_class <- coloc_classify(d, h4_strong = strong_thr, cond_thr = cond_thr,
                                                                           power_min = power_min, h3_strong = h3_strong)
        d
    }
    class_counts <- function(cls) {
        tb <- table(factor(cls, levels = cls_levels))
        paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " | ")
    }

    # Load colocalization results (classic coloc.abf)
    coloc <- tryCatch(readRDS("${coloc_results}"), error = function(e) {
        cat(sprintf("Could not read coloc table (%s); treating as empty\\n", conditionMessage(e)))
        NULL
    })
    if (is.null(coloc) || !is.data.frame(coloc)) coloc <- data.frame()
    coloc <- complete_metrics(coloc, h4_thr)

    # Load the SuSiE coloc summary when it exists (a NO_FILE* placeholder, an
    # unreadable file or an empty table all mean "SuSiE did not run").
    susie_path <- "${susie_table}"
    susie <- NULL
    if (!grepl("^NO_FILE", basename(susie_path))) {
        susie <- tryCatch(readRDS(susie_path), error = function(e) {
            cat(sprintf("Could not read SuSiE table (%s); treating as absent\\n", conditionMessage(e)))
            NULL
        })
        if (!is.null(susie) && !is.data.frame(susie)) susie <- NULL
        if (!is.null(susie)) susie <- complete_metrics(susie, susie_thr)
    } else {
        cat("SuSiE table not provided (placeholder file)\\n")
    }

    # Load the SAME significance rule used in the report: FDR (BH).
    # all_mr_results_combined.rds is already MHC-free and carries pval_fdr.
    thr_info <- tryCatch(readRDS("${threshold_info}"), error = function(e) NULL)
    fdr_alpha <- if (!is.null(thr_info) && !is.null(thr_info\$fdr_alpha)) thr_info\$fdr_alpha else 0.05
    # Safety net: (re)compute BH if the column is missing for any reason.
    if (nrow(mr_combined) > 0 && !("pval_fdr" %in% colnames(mr_combined))) {
        mr_combined\$pval_fdr <- p.adjust(mr_combined\$pval, method = "BH")
    }
    cat(sprintf("Plotting significance: FDR (BH) q < %g [MHC excluded upstream, primary tests only]\\n", fdr_alpha))

    # Load gene name conversion table (columns: id, name)
    gtf_reduced <- readRDS("${gene_conversion_table}")

    # Helper: map id.exposure (e.g. ENSG00000172260.15) -> gene symbol
    add_gene_names <- function(df) {
        ens <- sub(".*?(ENSG[0-9]+).*", "\\\\1", df\$id.exposure)
        nm <- gtf_reduced\$name[match(ens, gtf_reduced\$id)]
        nm[is.na(nm)] <- df\$id.exposure[is.na(nm)]
        df\$geneName <- nm
        df
    }

    # ---- Summary statistics ----
    summary_lines <- c(
        "",
        "===== MR Analysis Summary for ${params.tissue} =====",
        "",
        sprintf("Total Number of Genes Tested (primary tests): %d", length(unique(mr_combined\$id.exposure))),
        sprintf("Number of Primary MR Tests: %d", nrow(mr_combined)),
        sprintf("Number of Significant Associations (FDR q < %g): %d",
                fdr_alpha, sum(mr_combined\$pval_fdr < fdr_alpha, na.rm = TRUE)),
        sprintf("MR Methods Used: %s", paste(unique(mr_combined\$method), collapse = ", ")),
        "",
        "===== Colocalization Summary (classic coloc.abf) =====",
        "",
        sprintf("Decision rule (params.coloc_rule): %s", rule_str),
        sprintf("Genes Tested: %d", nrow(coloc)),
        sprintf("Genes passing the rule: %d", if (nrow(coloc) && "pass_rule" %in% colnames(coloc)) sum(coloc\$pass_rule %in% TRUE) else NA_integer_),
        sprintf("Genes with H4 > %s (strong tier): %d", h4_thr, sum(coloc\$H4 > h4_thr, na.rm = TRUE)),
        sprintf("Genes with H4/(H3+H4) >= %s and H3+H4 >= %s (conditional tier, Wu 2024 / Guo 2015): %d", cond_thr, power_min,
                if (nrow(coloc) && "H4_cond" %in% colnames(coloc)) sum(coloc\$H4_cond >= cond_thr & coloc\$H3_plus_H4 >= power_min, na.rm = TRUE) else NA_integer_),
        sprintf("Classes: %s", if (nrow(coloc) && "coloc_class" %in% colnames(coloc)) class_counts(coloc\$coloc_class) else "NA"),
        sprintf("Expected-FDR calibration (Reales 2026): smallest H4 cut with expected FDR < %g: %s", fdr_target,
                if (have_helper && nrow(coloc)) format(coloc_alpha_for_fdr(coloc\$H4, fdr_target), digits = 2) else "NA"),
        ""
    )

    # ---- SuSiE lines: rows and passing genes per LD_type ----
    susie_lines <- c("===== SuSiE Colocalization Summary =====", "")
    if (!is.null(susie) && nrow(susie) > 0 && "PP.H4.abf" %in% colnames(susie)) {
        lt  <- if ("LD_type" %in% colnames(susie)) as.character(susie\$LD_type) else rep("susie", nrow(susie))
        h4  <- suppressWarnings(as.numeric(susie\$PP.H4.abf))
        pr  <- if ("pass_rule" %in% colnames(susie)) susie\$pass_rule %in% TRUE else rep(NA, nrow(susie))
        gid <- if ("Gene_ID" %in% colnames(susie)) as.character(susie\$Gene_ID) else as.character(seq_len(nrow(susie)))
        susie_lines <- c(susie_lines, sprintf("Decision rule (params.coloc_rule): %s (per credible-set pair)", rule_str))
        for (ty in c("susie", "susie_bf", "abf_fallback", setdiff(unique(lt), c("susie", "susie_bf", "abf_fallback")))) {
            sel <- which(lt == ty)
            n_strong <- length(unique(gid[sel][!is.na(h4[sel]) & h4[sel] > susie_thr]))
            n_rule   <- length(unique(gid[sel][pr[sel] %in% TRUE]))
            susie_lines <- c(susie_lines,
                sprintf("LD_type %-13s : %d row(s), %d gene(s) passing the rule, %d gene(s) with PP.H4.abf > %s",
                        ty, length(sel), n_rule, n_strong, susie_thr))
        }
        susie_lines <- c(susie_lines, "",
            sprintf("Total SuSiE table rows: %d across %d gene(s); genes passing the rule (any LD_type): %d; via a real coloc.susie pair: %d",
                    nrow(susie), length(unique(gid)), length(unique(gid[pr %in% TRUE])), length(unique(gid[pr %in% TRUE & lt == "susie"]))))
        if ("coloc_class" %in% colnames(susie)) {
            # class of each gene's best real coloc.susie pair (highest PP.H4)
            sr <- susie[lt == "susie", , drop = FALSE]
            if (nrow(sr) > 0) {
                sr <- sr[order(sr\$Gene_ID, -suppressWarnings(as.numeric(sr\$PP.H4.abf))), , drop = FALSE]
                sr <- sr[!duplicated(sr\$Gene_ID), , drop = FALSE]
                susie_lines <- c(susie_lines, sprintf("Classes (best real coloc.susie pair per gene): %s", class_counts(sr\$coloc_class)))
            }
        }
    } else {
        susie_lines <- c(susie_lines, "SuSiE results not available (table absent or empty)")
    }

    writeLines(c(summary_lines, susie_lines, ""), "${params.tissue}_summary_stats.txt")

    # ---- Plot 1: MR volcano (effect size vs FDR q-value) ----
    if (nrow(mr_combined) > 0) {
        mr_combined\$is_signif <- mr_combined\$pval_fdr < fdr_alpha
        p1 <- ggplot(mr_combined, aes(x = b, y = -log10(pval_fdr), color = is_signif)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                               name = sprintf("FDR q < %g", fdr_alpha)) +
            geom_hline(yintercept = -log10(fdr_alpha), linetype = "dashed", color = "red") +
            labs(title = sprintf("MR Volcano Plot (primary tests) - %s", "${params.tissue}"),
                 x = "Effect Size (Beta)", y = "-log10(FDR q-value)") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("01_mr_volcano_plot.pdf", p1, width = 10, height = 6)
        ggsave("01_mr_volcano_plot.png", p1, width = 10, height = 6, dpi = 300)
    }

    # ---- Plot 2: coloc H4 distribution ----
    if (nrow(coloc) > 0 && "H4" %in% colnames(coloc)) {
        p2 <- ggplot(coloc, aes(x = H4)) +
            geom_histogram(binwidth = 0.05, fill = "steelblue", alpha = 0.7, boundary = 0) +
            geom_vline(xintercept = h4_thr, linetype = "dashed", color = "red") +
            labs(title = sprintf("Colocalization H4 Distribution - %s", "${params.tissue}"),
                 subtitle = sprintf("dashed: strong tier H4 > %s; rule in use: %s", h4_thr, rule_str),
                 x = "H4 Posterior Probability", y = "Count") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("02_coloc_h4_distribution.pdf", p2, width = 10, height = 6)
        ggsave("02_coloc_h4_distribution.png", p2, width = 10, height = 6, dpi = 300)
    }

    # ---- Plot 6: tiered colocalization classes (classic; + best real coloc.susie pair) ----
    if (nrow(coloc) > 0 && "coloc_class" %in% colnames(coloc)) {
        cc <- data.frame(method = "coloc.abf (classic)", class = factor(coloc\$coloc_class, levels = cls_levels))
        if (!is.null(susie) && nrow(susie) > 0 && all(c("coloc_class", "LD_type", "PP.H4.abf") %in% colnames(susie))) {
            sr <- susie[susie\$LD_type == "susie", , drop = FALSE]
            if (nrow(sr) > 0) {
                sr <- sr[order(sr\$Gene_ID, -suppressWarnings(as.numeric(sr\$PP.H4.abf))), , drop = FALSE]
                sr <- sr[!duplicated(sr\$Gene_ID), , drop = FALSE]
                cc <- rbind(cc, data.frame(method = "coloc.susie (best real pair)", class = factor(sr\$coloc_class, levels = cls_levels)))
            }
        }
        cc <- cc[!is.na(cc\$class), , drop = FALSE]
        cnt <- as.data.frame(table(method = cc\$method, class = cc\$class))
        p6 <- ggplot(cnt, aes(x = class, y = Freq, fill = class)) +
            geom_col(width = 0.7) + geom_text(aes(label = Freq), vjust = -0.3, size = 3.2) +
            facet_wrap(~ method) +
            scale_fill_manual(values = cls_colors, drop = FALSE, guide = "none") +
            scale_x_discrete(labels = function(x) gsub("_", "\\n", x)) +
            labs(title = sprintf("Colocalization classes - %s", "${params.tissue}"),
                 subtitle = sprintf("strong: H4 > %s | conditional: H4/(H3+H4) >= %s & H3+H4 >= %s | distinct: H3 > %s | underpowered: H3+H4 < %s",
                                    h4_thr, cond_thr, power_min, h3_strong, power_min),
                 x = NULL, y = "genes") +
            theme_minimal() + theme(panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("06_coloc_class_summary.pdf", p6, width = 10, height = 5)
        ggsave("06_coloc_class_summary.png", p6, width = 10, height = 5, dpi = 300)
    }

    # ---- Plot 7: classic H4 vs conditional H4/(H3+H4), coloured by H3+H4 (power) ----
    if (nrow(coloc) > 0 && all(c("H4", "H4_cond", "H3_plus_H4", "coloc_class") %in% colnames(coloc))) {
        cl <- coloc[is.finite(coloc\$H4) & is.finite(coloc\$H4_cond), , drop = FALSE]
        if (nrow(cl) > 0) {
            cl\$coloc_class <- factor(cl\$coloc_class, levels = cls_levels)
            p7 <- ggplot(cl, aes(x = H4, y = H4_cond, colour = H3_plus_H4, shape = coloc_class)) +
                geom_vline(xintercept = h4_thr, linetype = "dashed", colour = "grey50") +
                geom_hline(yintercept = cond_thr, linetype = "dashed", colour = "grey50") +
                geom_point(size = 2.6, alpha = 0.85) +
                scale_colour_viridis_c(name = "H3 + H4\\n(power)", limits = c(0, 1)) +
                scale_shape_manual(values = c(colocalized_strong = 16, colocalized_conditional = 17, inconclusive = 1,
                                              underpowered = 4, distinct_signals = 15), drop = FALSE, name = "class") +
                scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
                scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
                labs(title = sprintf("coloc.abf: H4 vs conditional H4/(H3+H4) - %s", "${params.tissue}"),
                     subtitle = sprintf("top-left quadrant = called only by the ratio; dark (low H3+H4) points there are H1-dominated = underpowered\\ndashed: H4 > %s (strong tier), H4/(H3+H4) >= %s (conditional tier)", h4_thr, cond_thr),
                     x = "PP.H4", y = "PP.H4 / (PP.H3 + PP.H4)") +
                theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            ggsave("07_coloc_h4_vs_conditional.pdf", p7, width = 9, height = 6.5)
            ggsave("07_coloc_h4_vs_conditional.png", p7, width = 9, height = 6.5, dpi = 300)
        }
    }

    # ---- Plot 3: method comparison (primary tests: Wald ratio vs multi-SNP method) ----
    if (nrow(mr_combined) > 0 && length(unique(mr_combined\$method)) > 1) {
        p3 <- ggplot(mr_combined, aes(x = method, y = -log10(pval_fdr))) +
            geom_boxplot(fill = "lightblue", alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.4) +
            geom_hline(yintercept = -log10(fdr_alpha), linetype = "dashed", color = "red") +
            labs(title = sprintf("MR Method Comparison (primary tests) - %s", "${params.tissue}"),
                 x = "Method", y = "-log10(FDR q-value)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("03_mr_method_comparison.pdf", p3, width = 10, height = 6)
        ggsave("03_mr_method_comparison.png", p3, width = 10, height = 6, dpi = 300)
    }

    # ---- Plot 4: coloc hypothesis summary ----
    h_cols <- c("H0", "H1", "H2", "H3", "H4")
    if (nrow(coloc) > 0 && all(h_cols %in% colnames(coloc))) {
        h_data <- coloc %>%
            pivot_longer(all_of(h_cols), names_to = "Hypothesis", values_to = "Probability") %>%
            mutate(Hypothesis = factor(Hypothesis, levels = h_cols))
        p4 <- ggplot(h_data, aes(x = Hypothesis, y = Probability, fill = Hypothesis)) +
            geom_boxplot(alpha = 0.7) +
            scale_fill_brewer(palette = "Set2") +
            labs(title = sprintf("Colocalization Hypotheses - %s", "${params.tissue}"),
                 x = "Hypothesis", y = "Posterior Probability") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
        ggsave("04_coloc_hypothesis_summary.pdf", p4, width = 10, height = 6)
        ggsave("04_coloc_hypothesis_summary.png", p4, width = 10, height = 6, dpi = 300)
    }

    # ---- Plot 5: FOREST PLOT of significant genes (with gene names) ----
    if (nrow(mr_combined) > 0) {
        # Same rule as the report: FDR (BH) on the MHC-free PRIMARY results.
        signif <- mr_combined %>% filter(!is.na(pval_fdr) & pval_fdr < fdr_alpha)

        cat(sprintf("Forest plot: %d primary MR results at FDR q < %g\\n", nrow(signif), fdr_alpha))

        if (nrow(signif) > 0) {
            # OR and 95% CI (same formulas as your original script)
            signif\$OR <- exp(signif\$b)
            signif\$CI_lower <- signif\$b - signif\$se * qnorm(0.975)
            signif\$CI_upper <- signif\$b + signif\$se * qnorm(0.975)

            # Map to gene names
            signif <- add_gene_names(signif)

            # Order by effect size; make labels unique so duplicate symbols don't collapse
            signif <- signif[order(signif\$b), ]
            signif\$geneLabel <- make.unique(as.character(signif\$geneName))
            signif\$geneLabel <- factor(signif\$geneLabel, levels = signif\$geneLabel)

            p_forest <- ggplot(signif, aes(x = b, y = geneLabel, xmin = CI_lower, xmax = CI_upper)) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
                geom_errorbarh(height = 0.25, color = "steelblue") +
                geom_point(size = 2, color = "darkblue") +
                labs(title = sprintf("Significant MR Genes (FDR q < %g) - %s", fdr_alpha, "${params.tissue}"),
                     x = "Effect (Beta)", y = "Gene") +
                theme_minimal() +
                theme(plot.title = element_text(hjust = 0.5, face = "bold"))

            h <- max(4, nrow(signif) * 0.35)
            ggsave("05_forest_plot_significant.pdf", p_forest, width = 10, height = h, limitsize = FALSE)
            ggsave("05_forest_plot_significant.png", p_forest, width = 10, height = h, dpi = 300, limitsize = FALSE)

            saveRDS(signif, "significant_results_with_names.rds")
            write.csv(signif, "significant_results_with_names.csv", row.names = FALSE)
        } else {
            cat("No FDR-significant genes for forest plot\\n")
        }
    }

    cat("Plotting complete\\n")
    """
}
