#!/usr/bin/env Rscript
# ----------------------------------------------------------------------------
# COLOC_DIAGNOSTICS worker (v2.2). Called by modules/coloc_diagnostics.nf; all
# inputs arrive as NF_* environment variables. Can also be run by hand inside a
# directory holding window_*.rds / ld_*.rds and the two coloc summary tables.
#
#   PART B  classic (coloc.abf) vs coloc.susie comparison: one table with H4,
#           H4/(H3+H4), H3+H4, tiered class and pass_rule for both methods, the
#           H4 scatter, the H4-vs-conditional plot and the class bar chart.
#   PART A  LocusZoom-style locus plots (exposure over outcome) for the genes
#           selected by NF_locus_plot_set, coloured by r2 with the index SNP
#           (LocusZoom bins), plus the cross-trait scatter with the same colours.
# ----------------------------------------------------------------------------
nf_env <- function(k, default = NULL) {
    v <- Sys.getenv(paste0("NF_", k), unset = NA)
    if (is.na(v) || !nzchar(v)) { if (is.null(default)) stop(sprintf("environment variable NF_%s not set", k)) else default } else v
}

log_con <- file("coloc_diagnostics_log.txt", open = "w")
lg <- function(...) { m <- sprintf(...); cat(m, "\n"); writeLines(m, log_con) }

suppressMessages({
    library(data.table); library(ggplot2)
    have_patch <- requireNamespace("patchwork", quietly = TRUE)
    if (have_patch) library(patchwork)
    have_repel <- requireNamespace("ggrepel", quietly = TRUE)
})

# ---- shared decision metrics (bin/coloc_decision.R) -------------------------
bin_dir <- nf_env("bin_dir", "")
helper <- file.path(bin_dir, "coloc_decision.R")
if (!nzchar(bin_dir) || !file.exists(helper))
    helper <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "coloc_decision.R")
if (!file.exists(helper)) helper <- Sys.which("coloc_decision.R")
if (!nzchar(helper) || !file.exists(helper)) stop("bin/coloc_decision.R not found next to coloc_diagnostics.R")
source(helper)

coloc_summary <- nf_env("coloc_summary")
susie_summary <- nf_env("susie_summary")
pwcoco_summary <- nf_env("pwcoco_summary", "")       # optional: NO_FILE_PWCOCO placeholder when PWCoCo did not run
gene_conv_file <- nf_env("gene_conversion_table", "")
tissue    <- nf_env("tissue", "tissue")
rule_str  <- nf_env("coloc_rule", "H4 > 0.8")
h4_thr    <- as.numeric(nf_env("coloc_h4_threshold", "0.8"))
susie_thr <- as.numeric(nf_env("susie_h4_threshold", "0.8"))
cond_thr  <- as.numeric(nf_env("coloc_cond_threshold", "0.7"))
power_min <- as.numeric(nf_env("coloc_power_min", "0.5"))
h3_strong <- as.numeric(nf_env("coloc_h3_strong", "0.8"))
plot_set  <- tolower(trimws(nf_env("locus_plot_set", "both")))
if (!(plot_set %in% c("both", "classic", "susie", "pwcoco", "all_methods", "any", "all"))) {
    lg("WARNING: unknown locus_plot_set '%s'; using 'both'", plot_set); plot_set <- "both"
}
cls_levels <- c("colocalized_strong", "colocalized_conditional", "inconclusive", "underpowered", "distinct_signals")
cls_colors <- c(colocalized_strong = "#1B7837", colocalized_conditional = "#7FBF7B", inconclusive = "#BDBDBD",
                underpowered = "#92C5DE", distinct_signals = "#D6604D")

lg("=== Coloc diagnostics (v2.2): method comparison + LocusZoom-style locus plots ===")
lg("rule: %s | strong tier H4 > %s (classic) / %s (SuSiE) | conditional >= %s with H3+H4 >= %s | distinct H3 > %s | locus plots: %s",
   rule_str, h4_thr, susie_thr, cond_thr, power_min, h3_strong, plot_set)

# ---- gene-name map --------------------------------------------------------
gene_conv <- tryCatch(readRDS(gene_conv_file), error = function(e) NULL)
gene_name_of <- function(ens) {
    if (is.null(gene_conv) || !all(c("id", "name") %in% colnames(gene_conv))) return(ens)
    e <- sub(".*?(ENSG[0-9]+).*", "\\1", ens)
    nm <- gene_conv$name[match(e, gene_conv$id)]
    ifelse(is.na(nm), ens, nm)
}
ens_of <- function(x) sub(".*?(ENSG[0-9]+).*", "\\1", as.character(x))

# ---- per-gene windows (SuSiE; full cis-window, both traits) ----------------
wf <- list.files(".", pattern = "^window_.*\\.rds$", full.names = TRUE)
if (length(wf) == 0) lg("No window_*.rds from SuSiE; locus/scatter plots cannot be made. (Did SuSiE run?)")
har <- data.table()
for (f in wf) {
    d <- tryCatch(as.data.table(readRDS(f)), error = function(e) NULL)
    if (!is.null(d) && nrow(d) > 0) har <- rbind(har, d, fill = TRUE)
}
if (nrow(har) > 0) {
    if (!"base_pair_location" %in% names(har) || all(is.na(har$base_pair_location)))
        har$base_pair_location <- suppressWarnings(as.numeric(sub("^.*:", "", har$SNP)))
    har$ens_core <- ens_of(har$id.exposure)
    if (!"pval.exposure" %in% names(har) || all(is.na(suppressWarnings(as.numeric(har$pval.exposure)))))
        har$pval.exposure <- 2 * pnorm(-abs(as.numeric(har$beta.exposure) / as.numeric(har$se.exposure)))
    if (!"pval.outcome" %in% names(har) || all(is.na(suppressWarnings(as.numeric(har$pval.outcome)))))
        har$pval.outcome <- 2 * pnorm(-abs(as.numeric(har$beta.outcome) / as.numeric(har$se.outcome)))
    lg("Loaded %d SNPs across %d gene windows", nrow(har), length(unique(har$ens_core)))
}
# LD matrices are read one gene at a time (they can be large); index them by ENSG.
ldf <- list.files(".", pattern = "^ld_ENSG[0-9]+\\.rds$")
ld_path_of <- setNames(ldf, sub("^ld_(ENSG[0-9]+)\\.rds$", "\\1", ldf))
lg("LD matrices available for %d gene(s)%s", length(ld_path_of),
   if (length(ld_path_of) == 0) " - locus plots will be grey (no r2 colouring)" else "")

# ---- coloc tables ------------------------------------------------------------
read_tab <- function(f) { x <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(x) || !is.data.frame(x)) data.frame() else as.data.frame(x) }
complete_metrics <- function(d, strong_thr) {
    # tables from v2.1 lack the metrics; derive them so the comparison never breaks
    if (nrow(d) == 0) return(d)
    if (!"H4_cond" %in% names(d)) d <- coloc_derive_metrics(d)
    if (!"pass_rule" %in% names(d)) d$pass_rule <- coloc_apply_rule(d, rule_str)
    if (!"coloc_class" %in% names(d)) d$coloc_class <- coloc_classify(d, h4_strong = strong_thr, cond_thr = cond_thr,
                                                                       power_min = power_min, h3_strong = h3_strong)
    d
}
classic <- complete_metrics(read_tab(coloc_summary), h4_thr)
susie   <- complete_metrics(read_tab(susie_summary), susie_thr)
if (nrow(classic) > 0) { classic$ens_core <- ens_of(classic$Gene_ID); classic$H4 <- as.numeric(if ("H4" %in% names(classic)) classic$H4 else classic$PP.H4.abf) }
if (nrow(susie) > 0) {
    if (!"LD_type" %in% names(susie)) susie$LD_type <- "susie"
    susie$ens_core <- ens_of(susie$Gene_ID); susie$PP.H4.abf <- as.numeric(susie$PP.H4.abf)
}
lg("classic table: %d gene(s); SuSiE table: %d row(s) (%d real coloc.susie pairs)",
   nrow(classic), nrow(susie), if (nrow(susie)) sum(susie$LD_type == "susie") else 0L)

# best REAL coloc.susie pair per gene (highest PP.H4); the ratio of THAT pair is
# reported, never the max of H4_cond across pairs (a weak pair can inflate it)
su_real <- if (nrow(susie)) susie[susie$LD_type == "susie", , drop = FALSE] else susie
su_best <- data.frame(); su_any <- data.frame()
if (nrow(su_real) > 0) {
    o <- order(su_real$ens_core, -su_real$PP.H4.abf)
    su_best <- su_real[o, , drop = FALSE]; su_best <- su_best[!duplicated(su_best$ens_core), , drop = FALSE]
    keys <- unique(su_real$ens_core)
    su_any <- data.frame(ens_core = keys,
                         susie_any_pass = as.logical(tapply(su_real$pass_rule %in% TRUE, su_real$ens_core, any)[keys]),
                         susie_n_pairs  = as.integer(tapply(rep(1L, nrow(su_real)), su_real$ens_core, sum)[keys]),
                         stringsAsFactors = FALSE)
}
max_by <- function(d, col, type) {
    d <- d[d$LD_type == type & is.finite(d$PP.H4.abf), , drop = FALSE]
    if (nrow(d) == 0) return(data.frame(ens_core = character(0), x = numeric(0)) |> setNames(c("ens_core", col)))
    a <- aggregate(PP.H4.abf ~ ens_core, data = d, FUN = max); names(a)[2] <- col; a
}

# ---- PWCoCo (optional): one row per signal pair, H0..H4 like the classic table ----
pwcoco <- data.frame(); pw_best <- data.frame(); pwcoco_ran <- FALSE
if (nzchar(pwcoco_summary) && !grepl("^NO_FILE", basename(pwcoco_summary)) && file.exists(pwcoco_summary) && file.size(pwcoco_summary) > 0) {
    pwcoco <- read_tab(pwcoco_summary)
    if (nrow(pwcoco) > 0 && all(c("Gene_ID", "H4") %in% names(pwcoco))) {
        pwcoco_ran <- TRUE
        pwcoco$H4 <- suppressWarnings(as.numeric(pwcoco$H4))
        pwcoco <- complete_metrics(pwcoco, h4_thr)
        pwcoco$ens_core <- ens_of(pwcoco$Gene_ID)
        if (!"pwcoco_status" %in% names(pwcoco)) pwcoco$pwcoco_status <- NA_character_
        pwv <- pwcoco[is.finite(pwcoco$H4), , drop = FALSE]
        if (nrow(pwv) > 0) {
            o <- order(pwv$ens_core, -pwv$H4)
            pw_best <- pwv[o, , drop = FALSE]; pw_best <- pw_best[!duplicated(pw_best$ens_core), , drop = FALSE]
            keys <- unique(pwv$ens_core)
            pw_best$pwcoco_any_pass <- as.logical(tapply(pwv$pass_rule %in% TRUE, pwv$ens_core, any)[pw_best$ens_core])
            pw_best$pwcoco_n_pairs  <- as.integer(tapply(rep(1L, nrow(pwv)), pwv$ens_core, sum)[pw_best$ens_core])
        }
        lg("PWCoCo table: %d row(s) over %d gene(s); %d gene(s) with a passing pair",
           nrow(pwcoco), length(unique(pwcoco$ens_core)), if (nrow(pw_best)) sum(pw_best$pwcoco_any_pass %in% TRUE) else 0L)
    } else {
        lg("PWCoCo table present but empty / without H4 column; ignored")
    }
} else {
    lg("PWCoCo table not supplied (placeholder): no pwcoco_* columns")
}

# =================================================================
# PART B - classic-vs-SuSiE method comparison
# =================================================================
comparison <- data.frame()
if (nrow(classic) > 0) {
    getc <- function(cn, default = NA) if (cn %in% names(classic)) classic[[cn]] else rep(default, nrow(classic))
    comparison <- data.frame(
        Gene_ID = classic$Gene_ID, ens_core = classic$ens_core, Gene_Name = gene_name_of(classic$Gene_ID),
        classic_H4 = classic$H4, classic_H3_plus_H4 = getc("H3_plus_H4"), classic_H4_cond = getc("H4_cond"),
        classic_class = getc("coloc_class"), classic_pass = getc("pass_rule", FALSE) %in% TRUE,
        classic_top_snp_H4 = getc("top_snp_H4"), classic_minp_out = getc("minp_out"),
        stringsAsFactors = FALSE)
    if (nrow(su_best) > 0) {
        sb <- data.frame(ens_core = su_best$ens_core, susie_H4 = su_best$PP.H4.abf,
                         susie_H3_plus_H4 = su_best$H3_plus_H4, susie_H4_cond = su_best$H4_cond,
                         susie_class = su_best$coloc_class, susie_pass = su_best$pass_rule %in% TRUE,
                         susie_hit1 = as.character(su_best$hit1), susie_hit2 = as.character(su_best$hit2),
                         stringsAsFactors = FALSE)
        comparison <- merge(comparison, sb, by = "ens_core", all.x = TRUE, sort = FALSE)
        comparison <- merge(comparison, su_any, by = "ens_core", all.x = TRUE, sort = FALSE)
    } else {
        comparison$susie_H4 <- NA_real_; comparison$susie_H3_plus_H4 <- NA_real_; comparison$susie_H4_cond <- NA_real_
        comparison$susie_class <- NA_character_; comparison$susie_pass <- FALSE
        comparison$susie_hit1 <- NA_character_; comparison$susie_hit2 <- NA_character_
        comparison$susie_any_pass <- FALSE; comparison$susie_n_pairs <- 0L
    }
    comparison$susie_any_pass[is.na(comparison$susie_any_pass)] <- FALSE
    comparison$susie_pass[is.na(comparison$susie_pass)] <- FALSE
    comparison$susie_n_pairs[is.na(comparison$susie_n_pairs)] <- 0L
    hy <- max_by(susie, "susie_bf_H4", "susie_bf"); fb <- max_by(susie, "fallback_H4", "abf_fallback")
    comparison <- merge(comparison, hy, by = "ens_core", all.x = TRUE, sort = FALSE)
    comparison <- merge(comparison, fb, by = "ens_core", all.x = TRUE, sort = FALSE)
    comparison$susie_status <- ifelse(!is.na(comparison$susie_H4), "ok",
                               ifelse(!is.na(comparison$susie_bf_H4), "susie_bf_only",
                               ifelse(!is.na(comparison$fallback_H4), "abf_fallback_only", "no_susie_result")))
    comparison$both_pass <- comparison$classic_pass & comparison$susie_any_pass
    # PWCoCo: best signal pair (highest H4) of each gene + any-pair flag; NA when it did not run
    if (nrow(pw_best) > 0) {
        pb <- data.frame(ens_core = pw_best$ens_core, pwcoco_H4 = pw_best$H4, pwcoco_H3_plus_H4 = pw_best$H3_plus_H4,
                         pwcoco_H4_cond = pw_best$H4_cond, pwcoco_class = pw_best$coloc_class,
                         pwcoco_pass = pw_best$pass_rule %in% TRUE, pwcoco_any_pass = pw_best$pwcoco_any_pass %in% TRUE,
                         pwcoco_n_pairs = pw_best$pwcoco_n_pairs,
                         pwcoco_best_pair = paste(if ("SNP1" %in% names(pw_best)) pw_best$SNP1 else NA, if ("SNP2" %in% names(pw_best)) pw_best$SNP2 else NA, sep = "|"),
                         pwcoco_status = as.character(pw_best$pwcoco_status), stringsAsFactors = FALSE)
        comparison <- merge(comparison, pb, by = "ens_core", all.x = TRUE, sort = FALSE)
        comparison$pwcoco_pass[is.na(comparison$pwcoco_pass)] <- FALSE
        comparison$pwcoco_any_pass[is.na(comparison$pwcoco_any_pass)] <- FALSE
        # genes PWCoCo saw but could not score (no_signal etc.) keep their status
        if (nrow(pwcoco) > 0) {
            st_of <- tapply(as.character(pwcoco$pwcoco_status), pwcoco$ens_core, function(v) v[1])
            miss <- is.na(comparison$pwcoco_status)
            comparison$pwcoco_status[miss] <- unname(st_of[comparison$ens_core[miss]])
        }
        comparison$pwcoco_status[is.na(comparison$pwcoco_status)] <- "no_pwcoco_result"
    } else {
        comparison$pwcoco_H4 <- NA_real_; comparison$pwcoco_H3_plus_H4 <- NA_real_; comparison$pwcoco_H4_cond <- NA_real_
        comparison$pwcoco_class <- NA_character_; comparison$pwcoco_pass <- NA; comparison$pwcoco_any_pass <- NA
        comparison$pwcoco_n_pairs <- NA_integer_; comparison$pwcoco_best_pair <- NA_character_
        comparison$pwcoco_status <- if (pwcoco_ran) "no_pwcoco_result" else "not_run"
    }
    comparison$all_methods_pass <- if (pwcoco_ran) comparison$both_pass & (comparison$pwcoco_any_pass %in% TRUE) else NA
    comparison$coloc_rule <- rule_str
    comparison <- comparison[order(-comparison$both_pass, -comparison$classic_H4), , drop = FALSE]
    out_cols <- c("Gene_ID", "Gene_Name", "classic_H4", "classic_H3_plus_H4", "classic_H4_cond", "classic_class", "classic_pass",
                  "classic_top_snp_H4", "classic_minp_out", "susie_H4", "susie_H3_plus_H4", "susie_H4_cond", "susie_class",
                  "susie_pass", "susie_any_pass", "susie_n_pairs", "susie_hit1", "susie_hit2", "susie_bf_H4", "fallback_H4",
                  "susie_status", "both_pass",
                  "pwcoco_H4", "pwcoco_H3_plus_H4", "pwcoco_H4_cond", "pwcoco_class", "pwcoco_pass", "pwcoco_any_pass",
                  "pwcoco_n_pairs", "pwcoco_best_pair", "pwcoco_status", "all_methods_pass", "coloc_rule")
    fwrite(comparison[, out_cols], "coloc_method_comparison_table.tsv", sep = "\t")
    lg("Comparison table: %d classic genes | pass classic rule: %d | real coloc.susie pair passes: %d | both: %d | PWCoCo pair passes: %s | all three: %s | status: %s",
       nrow(comparison), sum(comparison$classic_pass), sum(comparison$susie_any_pass), sum(comparison$both_pass),
       if (pwcoco_ran) sum(comparison$pwcoco_any_pass %in% TRUE) else "not run",
       if (pwcoco_ran) sum(comparison$all_methods_pass %in% TRUE) else "not run",
       paste(sprintf("%s=%d", names(table(comparison$susie_status)), as.integer(table(comparison$susie_status))), collapse = " "))

    # ---- classic H4 vs coloc.susie H4 -------------------------------------
    sub_ok <- comparison[comparison$susie_status == "ok", , drop = FALSE]
    n_hyonly <- sum(comparison$susie_status == "susie_bf_only"); n_fbonly <- sum(comparison$susie_status == "abf_fallback_only")
    n_nosusie <- sum(comparison$susie_status == "no_susie_result")
    c_pos <- sub_ok$classic_H4 > h4_thr; s_pos <- sub_ok$susie_H4 > susie_thr
    sub_txt <- sprintf("both H4 > thr: %d | both below: %d | disagree: %d | susie_bf only: %d | ABF-fallback only: %d | no result: %d",
                       sum(c_pos & s_pos, na.rm = TRUE), sum(!c_pos & !s_pos, na.rm = TRUE), sum(c_pos != s_pos, na.rm = TRUE),
                       n_hyonly, n_fbonly, n_nosusie)
    lg("  susie: %s", sub_txt)
    if (nrow(sub_ok) > 0) {
        sub_ok$agreement <- ifelse(sub_ok$both_pass, "pass both (rule)", ifelse(sub_ok$classic_pass, "classic only", ifelse(sub_ok$susie_any_pass, "coloc.susie only", "neither")))
        p_cmp <- ggplot(sub_ok, aes(classic_H4, susie_H4, colour = agreement)) +
            geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40") +
            geom_hline(yintercept = susie_thr, linetype = "dashed", colour = "grey60") +
            geom_vline(xintercept = h4_thr, linetype = "dashed", colour = "grey60") +
            geom_point(alpha = .8, size = 2.4) +
            scale_colour_manual(values = c("pass both (rule)" = "#1B7837", "classic only" = "#E6550D",
                                           "coloc.susie only" = "#3182BD", "neither" = "grey55"), name = rule_str) +
            scale_x_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
            scale_y_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
            labs(title = "Classic coloc.abf vs coloc.susie (best credible-set pair)", subtitle = sub_txt,
                 x = "coloc.abf  PP.H4", y = "coloc.susie  PP.H4 (best pair)",
                 caption = sprintf("%d genes scored by both methods; dashed = H4 thresholds (%s / %s); not shown: %d susie_bf-only, %d ABF-fallback-only",
                                   nrow(sub_ok), h4_thr, susie_thr, n_hyonly, n_fbonly)) +
            theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
        if (have_repel) {
            p_cmp <- p_cmp + ggrepel::geom_text_repel(aes(label = Gene_Name), size = 2.5, max.overlaps = 100, min.segment.length = 0,
                                                     segment.size = 0.2, segment.colour = "grey70", box.padding = 0.4, point.padding = 0.2, show.legend = FALSE)
        } else {
            p_cmp <- p_cmp + geom_text(aes(label = Gene_Name), size = 2.5, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE)
        }
        ggsave("method_comparison_susie.png", p_cmp, width = 7.5, height = 7, dpi = 150)
    }

    # ---- classic H4 vs H4/(H3+H4), coloured by H3+H4 (power) ----------------
    cl <- comparison[is.finite(comparison$classic_H4) & is.finite(comparison$classic_H4_cond), , drop = FALSE]
    if (nrow(cl) > 0) {
        p_r <- ggplot(cl, aes(classic_H4, classic_H4_cond, colour = classic_H3_plus_H4, shape = classic_class)) +
            geom_vline(xintercept = h4_thr, linetype = "dashed", colour = "grey50") +
            geom_hline(yintercept = cond_thr, linetype = "dashed", colour = "grey50") +
            geom_point(size = 2.6, alpha = .85) +
            scale_colour_viridis_c(name = "H3 + H4\n(power)", limits = c(0, 1)) +
            scale_shape_manual(values = c(colocalized_strong = 16, colocalized_conditional = 17, inconclusive = 1,
                                          underpowered = 4, distinct_signals = 15), drop = FALSE, name = "class") +
            scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
            labs(title = "coloc.abf: H4 vs conditional H4/(H3+H4)",
                 subtitle = sprintf("top-left quadrant = called only by the ratio; dark (low H3+H4) points there are H1-dominated = underpowered\ndashed: H4 > %s (strong tier), H4/(H3+H4) >= %s (conditional tier)", h4_thr, cond_thr),
                 x = "PP.H4", y = "PP.H4 / (PP.H3 + PP.H4)") +
            theme_minimal(base_size = 12) + theme(legend.position = "right")
        ggsave("coloc_h4_vs_conditional.png", p_r, width = 8.5, height = 6, dpi = 150)
    }

    # ---- tiered class counts, classic vs coloc.susie ---------------------------
    cc <- rbind(data.frame(method = "coloc.abf (classic)", class = factor(comparison$classic_class, levels = cls_levels)),
                if (nrow(su_best)) data.frame(method = "coloc.susie (best real pair)", class = factor(su_best$coloc_class, levels = cls_levels)) else NULL,
                if (nrow(pw_best)) data.frame(method = "PWCoCo (best signal pair)", class = factor(pw_best$coloc_class, levels = cls_levels)) else NULL)
    cc <- cc[!is.na(cc$class), , drop = FALSE]
    if (nrow(cc) > 0) {
        cc$method <- factor(cc$method, levels = c("coloc.abf (classic)", "coloc.susie (best real pair)", "PWCoCo (best signal pair)"))
        cnt <- as.data.frame(table(method = cc$method, class = cc$class))
        cnt <- cnt[cnt$method %in% unique(cc$method), , drop = FALSE]
        p_c <- ggplot(cnt, aes(class, Freq, fill = class)) +
            geom_col(width = 0.7) + geom_text(aes(label = Freq), vjust = -0.3, size = 3.2) +
            facet_wrap(~ method) + scale_fill_manual(values = cls_colors, drop = FALSE, guide = "none") +
            scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
            labs(title = sprintf("Colocalization classes - %s", tissue),
                 subtitle = sprintf("strong: H4 > %s | conditional: H4/(H3+H4) >= %s & H3+H4 >= %s | distinct: H3 > %s | underpowered: H3+H4 < %s",
                                    h4_thr, cond_thr, power_min, h3_strong, power_min),
                 x = NULL, y = "genes") +
            theme_minimal(base_size = 11) + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 8))
        n_meth <- length(unique(cc$method))
        ggsave("coloc_class_summary.png", p_c, width = 4.8 * n_meth + 0.6, height = 4.8, dpi = 150)
    }
} else {
    lg("Classic coloc table is empty; method comparison skipped")
}

# =================================================================
# PART A - LocusZoom-style locus plots for the selected genes
# =================================================================
sel_classic <- if (nrow(classic)) unique(classic$ens_core[classic$pass_rule %in% TRUE]) else character(0)
sel_susie   <- if (nrow(su_real)) unique(su_real$ens_core[su_real$pass_rule %in% TRUE]) else character(0)
sel_pwcoco  <- if (nrow(pw_best)) unique(pw_best$ens_core[pw_best$pwcoco_any_pass %in% TRUE]) else character(0)
with_window <- unique(har$ens_core)
genes <- switch(plot_set,
    both        = intersect(sel_classic, sel_susie),
    classic     = sel_classic,
    susie       = sel_susie,
    pwcoco      = sel_pwcoco,
    all_methods = if (pwcoco_ran) Reduce(intersect, list(sel_classic, sel_susie, sel_pwcoco)) else intersect(sel_classic, sel_susie),
    any         = Reduce(union, list(sel_classic, sel_susie, sel_pwcoco)),
    all         = with_window)
genes <- intersect(genes, with_window)
lg("Locus plot selection '%s': classic rule passed by %d gene(s), real coloc.susie pair by %d, PWCoCo pair by %s, classic+susie by %d; %d selected gene(s) have a window",
   plot_set, length(sel_classic), length(sel_susie), if (pwcoco_ran) as.character(length(sel_pwcoco)) else "n/a (not run)",
   length(intersect(sel_classic, sel_susie)), length(genes))

pick_index <- function(g, ens) {
    # returns list(i = row index in g, source = text)
    crow <- if (nrow(classic)) classic[classic$ens_core == ens, , drop = FALSE] else data.frame()
    if (nrow(crow) > 0 && "top_snp_H4" %in% names(crow)) {
        ts <- as.character(crow$top_snp_H4[1])
        if (!is.na(ts) && ts %in% g$SNP)
            return(list(i = which(g$SNP == ts)[1],
                        source = sprintf("coloc.abf top shared SNP (SNP.PP.H4 = %s)", format(crow$top_snp_PP_H4[1], digits = 2))))
    }
    sb <- if (nrow(su_best)) su_best[su_best$ens_core == ens, , drop = FALSE] else data.frame()
    if (nrow(sb) > 0 && "rsid" %in% names(g)) {
        for (h in c("hit2", "hit1")) {
            hv <- as.character(sb[[h]][1])
            if (!is.na(hv) && hv %in% g$rsid)
                return(list(i = which(g$rsid == hv)[1],
                            source = sprintf("coloc.susie %s (%s credible-set lead)", h, if (h == "hit2") "outcome" else "exposure")))
        }
    }
    list(i = which.max(g$logP_out), source = "outcome lead SNP (no coloc index available)")
}

plotted <- data.frame()
make_locus_plots <- function(ens) {
    g <- as.data.frame(har[har$ens_core == ens, ])
    if (nrow(g) < 3) { lg("  skip %s: <3 harmonised SNPs (%d)", ens, nrow(g)); return(invisible(NULL)) }
    prot <- gene_name_of(ens)
    g$bp <- as.numeric(g$base_pair_location)
    g$logP_exp <- -log10(pmax(as.numeric(g$pval.exposure), 1e-300))
    g$logP_out <- -log10(pmax(as.numeric(g$pval.outcome), 1e-300))
    g <- g[is.finite(g$bp) & is.finite(g$logP_exp) & is.finite(g$logP_out), , drop = FALSE]
    if (nrow(g) < 3) { lg("  skip %s: <3 SNPs with finite bp/p", ens); return(invisible(NULL)) }
    chr <- unique(as.character(g$chromosome))[1]

    # ---- index SNP + r2 with it (LocusZoom colouring) ----
    idx <- pick_index(g, ens)
    g$is_index <- seq_len(nrow(g)) == idx$i
    ld <- if (ens %in% names(ld_path_of)) tryCatch(readRDS(ld_path_of[[ens]]), error = function(e) { lg("  %s: cannot read LD (%s)", ens, e$message); NULL }) else NULL
    g$r2 <- NA_real_; ld_ids <- NULL
    if (!is.null(ld) && !is.null(rownames(ld))) {
        # the LD matrix is rsid-named (SuSiE saves it that way); fall back to chr:pos ids
        if ("rsid" %in% names(g) && mean(g$rsid %in% rownames(ld), na.rm = TRUE) > 0.5) ld_ids <- as.character(g$rsid)
        else if (mean(g$SNP %in% rownames(ld)) > 0.5) ld_ids <- as.character(g$SNP)
        if (!is.null(ld_ids)) g$r2 <- coloc_r2_with_index(ld, ld_ids, ld_ids[idx$i])
    }
    rm(ld)
    g$ld_bin <- coloc_r2_bin(g$r2, g$is_index)
    # lead MR instrument (window centre recorded by the coloc modules)
    crow <- if (nrow(classic)) classic[classic$ens_core == ens, , drop = FALSE] else data.frame()
    center_snp <- if (nrow(crow) && "window_center_snp" %in% names(crow)) as.character(crow$window_center_snp[1]) else NA_character_
    g$is_instr <- !is.na(center_snp) & g$SNP == center_snp
    index_label <- if ("rsid" %in% names(g) && !is.na(g$rsid[idx$i]) && nzchar(g$rsid[idx$i])) g$rsid[idx$i] else g$SNP[idx$i]
    n_r2 <- sum(is.finite(g$r2))

    cH4 <- if (nrow(crow)) crow$H4[1] else NA; cCond <- if (nrow(crow) && "H4_cond" %in% names(crow)) crow$H4_cond[1] else NA
    sb <- if (nrow(su_best)) su_best[su_best$ens_core == ens, , drop = FALSE] else data.frame()
    sH4 <- if (nrow(sb)) sb$PP.H4.abf[1] else NA
    pwr <- if (nrow(pw_best)) pw_best[pw_best$ens_core == ens, , drop = FALSE] else data.frame()
    pH4 <- if (nrow(pwr)) pwr$H4[1] else NA
    sub_txt <- sprintf("Exposure (eQTL) association | index %s (chr%s:%s): %s | r2 for %d/%d SNPs\ncoloc.abf H4 = %s, H4/(H3+H4) = %s | coloc.susie best pair H4 = %s%s",
                       index_label, chr, format(g$bp[idx$i], big.mark = ","), idx$source, n_r2, nrow(g),
                       format(cH4, digits = 2), format(cCond, digits = 2), format(sH4, digits = 2),
                       if (pwcoco_ran) sprintf(" | PWCoCo best pair H4 = %s", format(pH4, digits = 2)) else "")

    add_points <- function(p, ycol) {
        p +
            geom_point(data = g[!g$is_index & !g$is_instr, ], aes(x = bp / 1e6, y = .data[[ycol]], fill = ld_bin),
                       shape = 21, colour = "grey25", size = 1.9, stroke = 0.25, alpha = 0.9) +
            geom_point(data = g[g$is_instr & !g$is_index, ], aes(x = bp / 1e6, y = .data[[ycol]], fill = ld_bin),
                       shape = 24, colour = "black", size = 3.2, stroke = 0.7) +
            geom_point(data = g[g$is_index, ], aes(x = bp / 1e6, y = .data[[ycol]]),
                       shape = 23, fill = coloc_r2_colors[["index SNP"]], colour = "black", size = 4, stroke = 0.7) +
            scale_fill_manual(values = coloc_r2_colors, drop = FALSE, name = expression(r^2 ~ "with index SNP"),
                              guide = guide_legend(override.aes = list(shape = 21, size = 3.2))) +
            theme_minimal(base_size = 11) +
            theme(panel.grid.minor = element_blank(), legend.position = "right")
    }
    label_index <- function(p, ycol) {
        d <- g[g$is_index, ]; d$lab <- index_label
        if (have_repel) p + ggrepel::geom_text_repel(data = d, aes(x = bp / 1e6, y = .data[[ycol]], label = lab), size = 3,
                                                     min.segment.length = 0, box.padding = 0.6, point.padding = 0.3, show.legend = FALSE)
        else p + geom_text(data = d, aes(x = bp / 1e6, y = .data[[ycol]], label = lab), size = 3, vjust = -1, show.legend = FALSE)
    }
    p_exp <- add_points(ggplot(), "logP_exp") +
        labs(title = sprintf("%s (%s) - chr%s", prot, ens, chr), subtitle = sub_txt, x = NULL,
             y = expression(-log[10](P)[exposure]), caption = NULL)
    p_exp <- label_index(p_exp, "logP_exp") + theme(plot.subtitle = element_text(size = 8.5))
    p_out <- add_points(ggplot(), "logP_out") +
        geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "grey50") +
        labs(subtitle = "Outcome (GWAS) association", x = sprintf("chr%s position (Mb)", chr), y = expression(-log[10](P)[outcome]),
             caption = "diamond = index SNP; triangle = lead MR instrument (window centre); colours = LocusZoom r2 bins from the GAUSS LD used by SuSiE")
    p_out <- label_index(p_out, "logP_out")

    f1 <- sprintf("locus_%s_%s.png", ens, prot)
    if (have_patch) {
        ggsave(f1, (p_exp / p_out) + plot_layout(guides = "collect") & theme(legend.position = "right"), width = 9, height = 8, dpi = 150)
    } else {
        ggsave(sub("\\.png$", "_exposure.png", f1), p_exp, width = 9, height = 4.2, dpi = 150)
        ggsave(sub("\\.png$", "_outcome.png", f1), p_out, width = 9, height = 4.2, dpi = 150)
    }

    # ---- cross-trait scatter with the same colouring ----
    r <- suppressWarnings(cor(g$logP_exp, g$logP_out, use = "complete.obs"))
    p_sc <- ggplot() +
        geom_point(data = g[!g$is_index & !g$is_instr, ], aes(logP_exp, logP_out, fill = ld_bin), shape = 21, colour = "grey25", size = 2, stroke = 0.25, alpha = .9) +
        geom_point(data = g[g$is_instr & !g$is_index, ], aes(logP_exp, logP_out, fill = ld_bin), shape = 24, colour = "black", size = 3.2, stroke = 0.7) +
        geom_point(data = g[g$is_index, ], aes(logP_exp, logP_out), shape = 23, fill = coloc_r2_colors[["index SNP"]], colour = "black", size = 4, stroke = 0.7) +
        scale_fill_manual(values = coloc_r2_colors, drop = FALSE, name = expression(r^2 ~ "with index SNP"),
                          guide = guide_legend(override.aes = list(shape = 21, size = 3.2))) +
        labs(title = sprintf("%s (%s): cross-trait association", prot, ens),
             subtitle = sprintf("Pearson r = %.2f  (n = %d SNPs) | index %s", r, nrow(g), index_label),
             x = expression(-log[10](P)[exposure]), y = expression(-log[10](P)[outcome])) +
        theme_minimal(base_size = 11)
    ggsave(sprintf("scatter_%s_%s.png", ens, prot), p_sc, width = 7, height = 5.6, dpi = 150)

    plotted <<- rbind(plotted, data.frame(Gene_ID = ens, Gene_Name = prot, chromosome = chr, n_snps = nrow(g), n_with_r2 = n_r2,
                                          index_snp = index_label, index_pos = g$bp[idx$i], index_source = idx$source,
                                          lead_instrument = center_snp, classic_H4 = cH4, classic_H4_cond = cCond, susie_best_H4 = sH4,
                                          pwcoco_best_H4 = pH4, selection = plot_set, stringsAsFactors = FALSE))
    lg("  locus plots: %s (%s) n=%d r2=%d/%d index=%s [%s]", prot, ens, nrow(g), n_r2, nrow(g), index_label, idx$source)
    invisible(TRUE)
}

if (nrow(har) == 0) {
    lg("No per-gene window data available; skipping locus/scatter plots.")
} else if (length(genes) == 0) {
    lg("No gene satisfies locus_plot_set = '%s'; no locus plots made (set --locus_plot_set any|all to widen).", plot_set)
} else {
    lg("Making LocusZoom-style plots for %d gene(s)", length(genes))
    for (e in genes) tryCatch(make_locus_plots(e), error = function(err) lg("  %s: locus plot failed (%s)", e, err$message))
}
if (nrow(plotted) > 0) fwrite(plotted, "locus_plot_genes.tsv", sep = "\t")

lg("Done.")
close(log_con)
