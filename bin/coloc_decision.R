#!/usr/bin/env Rscript
# =============================================================================
#  bin/coloc_decision.R  -  shared colocalization decision metrics
# -----------------------------------------------------------------------------
#  Sourced by COLOCALIZATION (classic coloc.abf), SUSIE_COLOCALIZATION
#  (coloc.susie / coloc.susie_bf / abf_fallback rows) and the diagnostics /
#  plotting processes, so every table carries the SAME derived columns and the
#  SAME decision rule.
#
#  Works on either column naming scheme:
#      classic table : H0 H1 H2 H3 H4
#      SuSiE table   : PP.H0.abf ... PP.H4.abf
#
#  Metrics added by coloc_derive_metrics():
#      H3_plus_H4   P(both traits have a causal variant in the window) = H3 + H4.
#                   This is the "power" mass: when it is small, H0/H1/H2 dominate
#                   and NEITHER H4 nor any ratio says much (Guo et al. 2015 HMG
#                   only looked at pairs with PP3+PP4 >= 0.8).
#      H4_cond      H4 / (H3 + H4): the conditional (on both traits being
#                   associated) probability that the causal variant is shared.
#                   The "colocalization probability > 70%" of Wu et al. 2024
#                   Nat Commun and Liu et al. 2026 Commun Biol (TWMR papers).
#      H4_H3_odds   H4 / H3: posterior odds of H4 vs H3 (Guo et al. 2015 used
#                   PP4/PP3 > 5 "convincing", > 3 "weaker").
#                   H4_cond > 0.7 <=> odds > 2.33 ; 0.75 <=> 3 ; 0.833 <=> 5 ; 0.9 <=> 9
#      log2_H4_H3   log2(H4/H3), the form Open Targets Genetics reports.
#      dominant_hyp which of H0..H4 carries the most posterior mass.
#
#  Decision rule (coloc_apply_rule) uses the SAME syntax as coloc::sensitivity():
#      "H4 > 0.8"
#      "H4/(H3+H4) > 0.7 & H3+H4 > 0.5"
#      "H4 > 0.9 & H4/H3 > 3"           (example from the coloc documentation)
#  so the rule string in params.coloc_rule can be handed unchanged to
#  coloc::sensitivity() for the p12 sensitivity plot.
# =============================================================================

# ---- internal: map either naming scheme onto H0..H4 ------------------------
.coloc_std_h <- function(df) {
    pick <- function(a, b) if (a %in% names(df)) df[[a]] else if (b %in% names(df)) df[[b]] else NA_real_
    out <- data.frame(
        H0 = suppressWarnings(as.numeric(pick("H0", "PP.H0.abf"))),
        H1 = suppressWarnings(as.numeric(pick("H1", "PP.H1.abf"))),
        H2 = suppressWarnings(as.numeric(pick("H2", "PP.H2.abf"))),
        H3 = suppressWarnings(as.numeric(pick("H3", "PP.H3.abf"))),
        H4 = suppressWarnings(as.numeric(pick("H4", "PP.H4.abf"))),
        stringsAsFactors = FALSE)
    if (nrow(df) == 0) out <- out[0, , drop = FALSE]
    out
}

# ---- derived metrics --------------------------------------------------------
coloc_derive_metrics <- function(df, eps = 1e-12) {
    if (is.null(df) || nrow(df) == 0) {
        for (cc in c("H3_plus_H4", "H4_cond", "H4_H3_odds", "log2_H4_H3")) df[[cc]] <- numeric(0)
        df$dominant_hyp <- character(0)
        return(df)
    }
    h <- .coloc_std_h(df)
    both <- h$H3 + h$H4
    df$H3_plus_H4 <- both
    df$H4_cond    <- ifelse(is.finite(both) & both > eps, h$H4 / both, NA_real_)
    df$H4_H3_odds <- ifelse(is.finite(h$H3) & h$H3 > eps, h$H4 / h$H3, ifelse(is.finite(h$H4) & h$H4 > eps, Inf, NA_real_))
    df$log2_H4_H3 <- ifelse(is.finite(df$H4_H3_odds) & df$H4_H3_odds > 0, log2(df$H4_H3_odds), NA_real_)
    hm <- as.matrix(h)
    df$dominant_hyp <- apply(hm, 1, function(r) if (all(is.na(r))) NA_character_ else paste0("H", which.max(r) - 1L))
    df
}

# ---- decision rule ----------------------------------------------------------
# rule: string in coloc::sensitivity() syntax using H0..H4 (H3_plus_H4, H4_cond,
# H4_H3_odds and log2_H4_H3 are also visible once coloc_derive_metrics() ran).
# NA posteriors evaluate to FALSE (never "pass" on a missing value).
coloc_apply_rule <- function(df, rule) {
    if (is.null(df) || nrow(df) == 0) return(logical(0))
    rule <- trimws(as.character(rule))
    if (!nzchar(rule)) stop("coloc_apply_rule: empty rule")
    env <- c(.coloc_std_h(df), df[setdiff(names(df), c("H0","H1","H2","H3","H4"))])
    env <- as.data.frame(env, stringsAsFactors = FALSE, check.names = FALSE)
    res <- tryCatch(with(env, eval(parse(text = rule))), error = function(e)
        stop(sprintf("coloc_apply_rule: cannot evaluate rule '%s': %s", rule, e$message)))
    res <- as.logical(res)
    if (length(res) == 1L && nrow(df) > 1L) res <- rep(res, nrow(df))
    if (length(res) != nrow(df)) stop("coloc_apply_rule: rule did not return one value per row")
    res[is.na(res)] <- FALSE
    res
}

# ---- tiered classification ----------------------------------------------------
# One label per row, in decreasing order of evidence FOR a shared variant:
#   colocalized_strong       H4 > h4_strong                       (classic call)
#   colocalized_conditional  H4 <= h4_strong, H3+H4 >= power_min and H4/(H3+H4) >= cond_thr
#                            (both traits associated; shared variant is the better explanation)
#   distinct_signals         H3 > h3_strong: both associated but different causal variants
#                            -> the MR estimate is at risk of LD confounding
#   underpowered             H3+H4 < power_min: H0/H1/H2 dominate. Inconclusive,
#                            NOT evidence against colocalization (dominant_hyp says which)
#   inconclusive             everything else (e.g. H3 ~ H4 with both moderate)
coloc_classify <- function(df, h4_strong = 0.8, cond_thr = 0.7, power_min = 0.5, h3_strong = 0.8) {
    if (is.null(df) || nrow(df) == 0) return(character(0))
    if (!"H4_cond" %in% names(df)) df <- coloc_derive_metrics(df)
    h <- .coloc_std_h(df)
    cls <- rep("inconclusive", nrow(df))
    cls[is.na(h$H4)] <- NA_character_
    under <- is.finite(df$H3_plus_H4) & df$H3_plus_H4 < power_min
    cls[under] <- "underpowered"
    dist <- is.finite(h$H3) & h$H3 > h3_strong
    cls[dist] <- "distinct_signals"
    cond <- is.finite(df$H4_cond) & df$H4_cond >= cond_thr & is.finite(df$H3_plus_H4) & df$H3_plus_H4 >= power_min
    cls[cond] <- "colocalized_conditional"
    strong <- is.finite(h$H4) & h$H4 > h4_strong
    cls[strong] <- "colocalized_strong"
    cls
}

# ---- expected FDR of a call set (Reales et al. 2026 PLoS Genet) --------------
# If PP.H4 is calibrated, the expected FDR among calls with H4 > alpha is
# mean(1 - H4) over those calls. Same idea for the conditional rule, using
# 1 - H4_cond over the calls (FDR conditional on both traits being associated).
coloc_fdr_table <- function(h4, alphas = c(0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)) {
    h4 <- suppressWarnings(as.numeric(h4)); h4 <- h4[is.finite(h4)]
    out <- lapply(alphas, function(a) {
        called <- h4[h4 > a]
        data.frame(alpha = a, n_called = length(called),
                   expected_FDR = if (length(called)) mean(1 - called) else NA_real_,
                   expected_false = if (length(called)) sum(1 - called) else NA_real_)
    })
    do.call(rbind, out)
}
coloc_alpha_for_fdr <- function(h4, target = 0.05, grid = seq(0.50, 0.99, by = 0.01)) {
    tab <- coloc_fdr_table(h4, grid)
    ok <- tab[is.finite(tab$expected_FDR) & tab$expected_FDR < target, , drop = FALSE]
    if (nrow(ok) == 0) return(NA_real_)
    min(ok$alpha)
}

# ---- column tag for a p12 value: 1e-05 -> "1e05", 5e-06 -> "5e06" -------------
coloc_p12_tag <- function(p) gsub("[^0-9a-zA-Z]", "", formatC(as.numeric(p), format = "e", digits = 0))

# parse "1e-5,5e-6,1e-6" (params.coloc_p12_grid) into a numeric vector
coloc_parse_p12_grid <- function(s, default = c(1e-5, 5e-6, 1e-6)) {
    v <- suppressWarnings(as.numeric(trimws(strsplit(as.character(s), ",")[[1]])))
    v <- v[is.finite(v) & v > 0]
    if (length(v) == 0) default else v
}

# ---- prior (p12) sensitivity for ONE coloc result ------------------------------
# Re-runs coloc.abf on the same two datasets for each p12 in p12_grid (cheap) and
# reports H4 and H4_cond at each, plus - via coloc::sensitivity() - the smallest
# p12 at which `rule` still passes. Returns a one-row data.frame.
coloc_p12_profile <- function(D1, D2, rule, p12_grid = c(1e-5, 5e-6, 1e-6), p1 = 1e-4, p2 = 1e-4, res = NULL) {
    out <- list()
    for (p in p12_grid) {
        # capture.output(): coloc.abf print()s its summary on every call
        r <- tryCatch({ z <- NULL; utils::capture.output(z <- suppressMessages(coloc::coloc.abf(D1, D2, p1 = p1, p2 = p2, p12 = p))); z },
                      error = function(e) NULL)
        tag <- coloc_p12_tag(p)   # 1e-05 -> 1e05
        if (is.null(r)) { out[[paste0("H4_p12_", tag)]] <- NA_real_; out[[paste0("H4cond_p12_", tag)]] <- NA_real_; next }
        s <- as.list(r$summary)
        h3 <- as.numeric(s$PP.H3.abf); h4 <- as.numeric(s$PP.H4.abf)
        out[[paste0("H4_p12_", tag)]]     <- h4
        out[[paste0("H4cond_p12_", tag)]] <- if ((h3 + h4) > 1e-12) h4 / (h3 + h4) else NA_real_
    }
    # smallest p12 for which the rule passes (sensitivity() scans p1*p2 .. min(p1,p2))
    p12_min <- NA_real_; p12_max <- NA_real_
    if (!is.null(res) && "sensitivity" %in% getNamespaceExports("coloc")) {
        sens <- tryCatch(suppressMessages(coloc::sensitivity(res, rule = rule, doplot = FALSE)), error = function(e) NULL)
        if (!is.null(sens) && all(c("p12", "pass") %in% names(sens)) && any(sens$pass, na.rm = TRUE)) {
            p12_min <- min(sens$p12[sens$pass %in% TRUE]); p12_max <- max(sens$p12[sens$pass %in% TRUE])
        }
    }
    out$p12_min_pass <- p12_min
    out$p12_max_pass <- p12_max
    as.data.frame(out, stringsAsFactors = FALSE)
}

# ---- per-trait lead SNPs / credible set from a coloc.abf result ---------------
coloc_abf_hits <- function(res, cs_coverage = 0.95) {
    rr <- tryCatch(as.data.frame(res$results), error = function(e) NULL)
    out <- data.frame(hit_exp = NA_character_, hit_out = NA_character_,
                      top_snp_H4 = NA_character_, top_snp_PP_H4 = NA_real_, cs95_size_H4 = NA_integer_,
                      stringsAsFactors = FALSE)
    if (is.null(rr) || !"snp" %in% names(rr) || nrow(rr) == 0) return(out)
    if ("lABF.df1" %in% names(rr)) out$hit_exp <- as.character(rr$snp[which.max(rr$lABF.df1)])
    if ("lABF.df2" %in% names(rr)) out$hit_out <- as.character(rr$snp[which.max(rr$lABF.df2)])
    if ("SNP.PP.H4" %in% names(rr)) {
        o <- order(rr$SNP.PP.H4, decreasing = TRUE)
        out$top_snp_H4    <- as.character(rr$snp[o[1]])
        out$top_snp_PP_H4 <- as.numeric(rr$SNP.PP.H4[o[1]])
        out$cs95_size_H4  <- as.integer(which(cumsum(rr$SNP.PP.H4[o]) >= cs_coverage)[1])
    }
    out
}

# ---- LD check (Zheng et al. 2020 Nat Genet; reviewed in Zuber et al. 2022) -----
# max r2 between any MR instrument and the top_n outcome SNPs (by p) in the window.
# ld_r: signed correlation matrix with row/col names = snp ids in `snp_ids`.
coloc_ld_check <- function(ld_r, snp_ids, instrument_ids, outcome_p, top_n = 30, r2_thr = 0.8) {
    out <- data.frame(ld_check_max_r2 = NA_real_, ld_check_pass = NA, ld_check_n_instr = 0L, stringsAsFactors = FALSE)
    if (is.null(ld_r) || length(snp_ids) == 0) return(out)
    ins <- which(snp_ids %in% instrument_ids)
    if (length(ins) == 0) return(out)
    op <- suppressWarnings(as.numeric(outcome_p)); op[!is.finite(op)] <- 1
    top <- order(op)[seq_len(min(top_n, length(op)))]
    sub <- ld_r[ins, top, drop = FALSE]
    m <- suppressWarnings(max(sub^2, na.rm = TRUE)); if (!is.finite(m)) m <- NA_real_
    out$ld_check_max_r2 <- m; out$ld_check_pass <- if (is.na(m)) NA else m >= r2_thr; out$ld_check_n_instr <- length(ins)
    out
}

# ---- LocusZoom-style r2 bins (for the locus plots in COLOC_DIAGNOSTICS) --------
# r2 with the index SNP, binned exactly as LocusZoom does; the index SNP itself is
# drawn as a purple diamond. Bins are ordered so legends read top-down from high LD.
coloc_r2_levels <- c("index SNP", "0.8 - 1.0", "0.6 - 0.8", "0.4 - 0.6", "0.2 - 0.4", "< 0.2", "no LD info")
coloc_r2_colors <- c("index SNP" = "#7B2CBF", "0.8 - 1.0" = "#D7191C", "0.6 - 0.8" = "#F58220",
                     "0.4 - 0.6" = "#2E9E44", "0.2 - 0.4" = "#74B9E6", "< 0.2" = "#1B2F75",
                     "no LD info" = "#BDBDBD")
coloc_r2_bin <- function(r2, is_index = rep(FALSE, length(r2))) {
    r2 <- suppressWarnings(as.numeric(r2))
    b <- ifelse(is.na(r2), "no LD info",
         ifelse(r2 >= 0.8, "0.8 - 1.0",
         ifelse(r2 >= 0.6, "0.6 - 0.8",
         ifelse(r2 >= 0.4, "0.4 - 0.6",
         ifelse(r2 >= 0.2, "0.2 - 0.4", "< 0.2")))))
    b[is_index %in% TRUE] <- "index SNP"
    factor(b, levels = coloc_r2_levels)
}
# r2 of every SNP in `snp_ids` with `index_id`, read off a signed-R matrix whose
# dimnames are SNP ids (ld_<ENSG>.rds from SUSIE_COLOCALIZATION). NA where absent.
coloc_r2_with_index <- function(ld_r, snp_ids, index_id) {
    out <- rep(NA_real_, length(snp_ids))
    if (is.null(ld_r) || is.null(rownames(ld_r)) || !(index_id %in% rownames(ld_r))) return(out)
    hit <- match(snp_ids, rownames(ld_r))
    ok <- !is.na(hit)
    out[ok] <- as.numeric(ld_r[hit[ok], index_id])^2
    out
}
