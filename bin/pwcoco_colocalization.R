#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# PWCoCo helper for the TWMR pipeline (called by modules/pwcoco_colocalization.nf)
#
#   Rscript pwcoco_colocalization.R prep   <exposure_N> <match_mode>
#   Rscript pwcoco_colocalization.R align  <match_mode> <min_snps>
#   Rscript pwcoco_colocalization.R parse  <gene_conversion_table.rds> <gene_map.rds>
#
# prep   : window_<ENSG>.rds (SuSiE windows)  ->  pwcoco_inputs/<ENSG>.window.tsv
#          + pwcoco_regions.tsv (ens chr start end) + pwcoco_prep_status.tsv
# align  : per-gene PLINK extract (pwcoco_panel/<ENSG>.bim/.frq)  ->  rewritten .bim
#          with unique chr:bp:A1:A2 ids + the two COJO-format .ma files PWCoCo reads
#          + pwcoco_manifest.tsv (ens chr) + pwcoco_align_status.tsv
# parse  : pwcoco_per_gene/<ENSG>.coloc + logs  ->  pwcoco_results_summary.rds,
#          pwcoco_gene_status.tsv
#
# Everything PWCoCo needs is derived from the same harmonised window SuSiE and
# SharePro saw, so the three methods are compared on one SNP set.
# ---------------------------------------------------------------------------
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
step <- if (length(args) >= 1) args[1] else stop("usage: pwcoco_colocalization.R <prep|align|parse> ...")

log_con <- file("pwcoco_run_log.txt", open = "a")
lg <- function(...) { m <- sprintf(...); cat(m, "\n"); writeLines(m, log_con) }
comp <- function(a) chartr("ACGT", "TGCA", a)
num  <- function(x) suppressWarnings(as.numeric(x))
ens_core <- function(x) sub(".*?(ENSG[0-9]+).*", "\\1", as.character(x))

dir.create("pwcoco_inputs",   showWarnings = FALSE)
dir.create("pwcoco_panel",    showWarnings = FALSE)
dir.create("pwcoco_per_gene", showWarnings = FALSE)

# ===========================================================================
if (step == "prep") {
    n_exp_param <- num(args[2])
    mode <- if (length(args) >= 3) args[3] else "position"
    wf <- list.files(".", pattern = "^window_.*\\.rds$")
    lg("prep: %d SuSiE window files, match mode = %s", length(wf), mode)
    regions <- list(); status <- list()
    for (f in wf) {
        ens <- sub("^window_(.*)\\.rds$", "\\1", f)
        win <- tryCatch(as.data.frame(readRDS(f)), error = function(e) NULL)
        if (is.null(win) || nrow(win) < 3 ||
            !all(c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome",
                   "effect_allele.exposure", "other_allele.exposure") %in% names(win))) {
            status[[ens]] <- data.table(Gene_ID = ens, prep_status = "bad_window", n_window = 0L)
            lg("  %s: unreadable window or missing columns", ens); next
        }
        chr <- num(sub("^chr", "", as.character(win$chromosome[1])))
        if (!is.finite(chr)) chr <- num(sub("^chr", "", sub(":.*$", "", win$SNP[1])))
        pos <- num(sub("^.*:", "", win$SNP))
        Nexp <- if (is.finite(n_exp_param) && n_exp_param > 0) n_exp_param else
                    suppressWarnings(max(num(win$samplesize.exposure), na.rm = TRUE))
        Nout <- suppressWarnings(mean(num(win$samplesize.outcome), na.rm = TRUE))
        if (!is.finite(Nexp) || Nexp <= 0 || !is.finite(Nout) || Nout <= 0) {
            status[[ens]] <- data.table(Gene_ID = ens, prep_status = "bad_window", n_window = 0L)
            lg("  %s: cannot determine sample sizes", ens); next
        }
        d <- data.table(
            key_snp  = as.character(win$SNP),
            rsid     = if ("rsid" %in% names(win)) as.character(win$rsid) else NA_character_,
            chr      = as.character(as.integer(chr)),
            pos      = as.integer(pos),
            A1       = toupper(as.character(win$effect_allele.exposure)),
            A2       = toupper(as.character(win$other_allele.exposure)),
            freq_exp = if ("eaf.exposure" %in% names(win)) num(win$eaf.exposure) else NA_real_,
            beta_exp = num(win$beta.exposure), se_exp = num(win$se.exposure),
            p_exp    = if ("pval.exposure" %in% names(win)) num(win$pval.exposure) else NA_real_,
            N_exp    = round(Nexp),
            beta_out = num(win$beta.outcome), se_out = num(win$se.outcome),
            p_out    = if ("pval.outcome" %in% names(win)) num(win$pval.outcome) else NA_real_,
            N_out    = round(Nout),
            freq_out = if ("eaf.outcome" %in% names(win)) num(win$eaf.outcome) else NA_real_
        )
        d[is.na(p_exp), p_exp := 2 * pnorm(-abs(beta_exp / se_exp))]
        d[is.na(p_out), p_out := 2 * pnorm(-abs(beta_out / se_out))]
        ok <- is.finite(d$pos) & is.finite(d$beta_exp) & is.finite(d$se_exp) & d$se_exp > 0 &
              is.finite(d$beta_out) & is.finite(d$se_out) & d$se_out > 0 &
              d$A1 %in% c("A", "C", "G", "T") & d$A2 %in% c("A", "C", "G", "T") & d$A1 != d$A2
        d <- d[ok]
        d <- d[!duplicated(key_snp)]
        if (nrow(d) < 3) {
            status[[ens]] <- data.table(Gene_ID = ens, prep_status = "bad_window", n_window = nrow(d))
            lg("  %s: fewer than 3 usable SNPs", ens); next
        }
        if (mode == "rsid") {
            rs <- unique(d$rsid[!is.na(d$rsid) & nzchar(d$rsid)])
            if (length(rs) < 3) {
                status[[ens]] <- data.table(Gene_ID = ens, prep_status = "no_rsid", n_window = nrow(d))
                lg("  %s: rsid matching requested but the window carries no rsIDs", ens); next
            }
            writeLines(rs, sprintf("pwcoco_inputs/%s.rsids.txt", ens))
        }
        fwrite(d, sprintf("pwcoco_inputs/%s.window.tsv", ens), sep = "\t")
        regions[[ens]] <- data.table(ens = ens, chr = d$chr[1], start = min(d$pos), end = max(d$pos))
        status[[ens]] <- data.table(Gene_ID = ens, prep_status = "ok", n_window = nrow(d))
        lg("  %s: chr%s:%d-%d, %d SNPs, N_exp=%d N_out=%d", ens, d$chr[1], min(d$pos), max(d$pos),
           nrow(d), as.integer(round(Nexp)), as.integer(round(Nout)))
    }
    reg <- if (length(regions)) rbindlist(regions) else data.table(ens = character(), chr = character(), start = integer(), end = integer())
    fwrite(reg, "pwcoco_regions.tsv", sep = "\t", col.names = FALSE)
    st <- if (length(status)) rbindlist(status) else data.table(Gene_ID = character(), prep_status = character(), n_window = integer())
    fwrite(st, "pwcoco_prep_status.tsv", sep = "\t")
    lg("prep: %d of %d windows ready for panel extraction", nrow(reg), length(wf))

# ===========================================================================
} else if (step == "align") {
    mode <- if (length(args) >= 2) args[2] else "position"
    min_snps <- if (length(args) >= 3) as.integer(args[3]) else 10L
    ext <- if (file.exists("pwcoco_extract_status.tsv") && file.size("pwcoco_extract_status.tsv") > 0)
               fread("pwcoco_extract_status.tsv", header = FALSE, col.names = c("ens", "extract_status", "chr")) else
               data.table(ens = character(), extract_status = character(), chr = character())
    lg("align: %d genes extracted from the panel", sum(ext$extract_status == "extracted"))
    manifest <- list(); status <- list()
    for (i in seq_len(nrow(ext))) {
        ens <- ext$ens[i]; chr <- as.character(ext$chr[i])
        rec <- function(st, n_pan = NA_integer_, n_mat = NA_integer_, n_cmp = NA_integer_, n_mis = NA_integer_)
            data.table(Gene_ID = ens, align_status = st, n_panel = n_pan, n_matched = n_mat,
                       n_complemented = n_cmp, n_allele_mismatch = n_mis)
        if (ext$extract_status[i] != "extracted") { status[[ens]] <- rec(ext$extract_status[i]); next }
        bimf <- sprintf("pwcoco_panel/%s.bim", ens)
        bim <- tryCatch(fread(bimf, header = FALSE, colClasses = "character"), error = function(e) NULL)
        win <- tryCatch(fread(sprintf("pwcoco_inputs/%s.window.tsv", ens), colClasses = list(character = c("chr", "rsid", "key_snp"))),
                        error = function(e) NULL)
        if (is.null(bim) || ncol(bim) < 6 || nrow(bim) == 0 || is.null(win)) {
            status[[ens]] <- rec("bad_panel_extract"); lg("  %s: unreadable .bim/window", ens); next
        }
        setnames(bim, 1:6, c("chr", "id", "cm", "bp", "a1", "a2"))
        bim[, chr := sub("^chr", "", chr)]
        bim[, a1 := toupper(a1)]; bim[, a2 := toupper(a2)]
        bim[, bp := as.integer(bp)]
        bim[, row := .I]
        # ---- panel allele frequency of a1 (plink 1.9 .frq or plink 2 .afreq) ----
        bim[, fa1 := NA_real_]
        frqf <- sprintf("pwcoco_panel/%s.frq", ens); afrf <- sprintf("pwcoco_panel/%s.afreq", ens)
        if (file.exists(frqf)) {
            fq <- tryCatch(fread(frqf, colClasses = list(character = c("SNP", "A1", "A2"))), error = function(e) NULL)
            if (!is.null(fq) && nrow(fq) == nrow(bim)) {
                bim[, fa1 := ifelse(toupper(fq$A1) == a1, num(fq$MAF),
                             ifelse(toupper(fq$A2) == a1, 1 - num(fq$MAF), NA_real_))]
            }
        } else if (file.exists(afrf)) {
            fq <- tryCatch(fread(afrf, colClasses = list(character = c("ID", "REF", "ALT"))), error = function(e) NULL)
            if (!is.null(fq) && nrow(fq) == nrow(bim) && all(c("REF", "ALT", "ALT_FREQS") %in% names(fq))) {
                bim[, fa1 := ifelse(toupper(fq$ALT) == a1, num(fq$ALT_FREQS),
                             ifelse(toupper(fq$REF) == a1, 1 - num(fq$ALT_FREQS), NA_real_))]
            }
        }
        if (all(is.na(bim$fa1))) lg("  %s: WARNING no usable panel frequencies (.frq/.afreq); GTEx eaf used for both traits", ens)
        # A constant eaf (the pipeline's 0.3 placeholder when a file has no frequency
        # column) is not a frequency: treat it as missing so the panel value is used.
        # The windows are TwoSampleMR-harmonised, so both eaf columns refer to A1.
        if (length(unique(round(win$freq_exp[is.finite(win$freq_exp)], 4))) <= 1) {
            lg("  %s: exposure eaf is a placeholder/constant; panel frequencies used for the exposure", ens)
            win[, freq_exp := NA_real_]
        }
        if (!"freq_out" %in% names(win)) win[, freq_out := NA_real_]
        if (length(unique(round(win$freq_out[is.finite(win$freq_out)], 4))) <= 1) {
            lg("  %s: outcome eaf is a placeholder/constant; panel frequencies used for the outcome", ens)
            win[, freq_out := NA_real_]
        } else {
            # SNPs the GWAS had no frequency for carry the 0.3 placeholder: fall back per SNP
            win[abs(freq_out - 0.3) < 1e-9, freq_out := NA_real_]
        }
        # unique, build-free ids: chr:bp:a1:a2
        bim[, id_new := sprintf("%s:%d:%s:%s", chr, bp, a1, a2)]
        # ---- match window SNPs to panel rows ----
        if (mode == "rsid") {
            cand <- merge(win[!is.na(rsid) & nzchar(rsid)], bim, by.x = "rsid", by.y = "id", allow.cartesian = TRUE)
        } else {
            cand <- merge(win, bim, by.x = c("chr", "pos"), by.y = c("chr", "bp"), allow.cartesian = TRUE)
        }
        n_pan <- nrow(bim)
        if (nrow(cand) == 0) { status[[ens]] <- rec("no_panel_overlap", n_pan, 0L, 0L, 0L); lg("  %s: no window SNP found in the panel extract", ens); next }
        cand[, mtype := fifelse((A1 == a1 & A2 == a2) | (A1 == a2 & A2 == a1), "same",
                        fifelse((comp(A1) == a1 & comp(A2) == a2) | (comp(A1) == a2 & comp(A2) == a1), "complement", "mismatch"))]
        n_mis <- sum(cand$mtype == "mismatch")
        cand <- cand[mtype != "mismatch"]
        cand[, rank := fifelse(mtype == "same", 1L, 2L)]
        setorder(cand, key_snp, rank)
        cand <- cand[!duplicated(key_snp)]      # one panel row per window SNP (prefer same-strand)
        cand <- cand[!duplicated(row)]          # one window SNP per panel row
        cand <- cand[!duplicated(id_new)]       # exact duplicate panel rows: PWCoCo renames them
        n_cmp <- sum(cand$mtype == "complement")
        cand[mtype == "complement", `:=`(A1 = comp(A1), A2 = comp(A2))]
        # frequency of OUR effect allele (A1) in the panel
        cand[, fpan := fifelse(A1 == a1, fa1, 1 - fa1)]
        # Each trait keeps its own allele frequency when it is a real value (GTEx af
        # for the eQTL, the GWAS file's EAF for the outcome); otherwise the panel
        # frequency is used, and as a last resort the other trait's value.
        cand[, freq_exp := fifelse(is.finite(freq_exp) & freq_exp > 0 & freq_exp < 1, freq_exp, fpan)]
        cand[, freq_out := fifelse(is.finite(freq_out) & freq_out > 0 & freq_out < 1, freq_out, fpan)]
        cand[!(is.finite(freq_exp) & freq_exp > 0 & freq_exp < 1), freq_exp := freq_out]
        cand[!(is.finite(freq_out) & freq_out > 0 & freq_out < 1), freq_out := freq_exp]
        cand <- cand[is.finite(freq_exp) & freq_exp > 0 & freq_exp < 1 & is.finite(freq_out) & freq_out > 0 & freq_out < 1]
        n_mat <- nrow(cand)
        if (n_mat < min_snps) {
            status[[ens]] <- rec("no_panel_overlap", n_pan, n_mat, n_cmp, n_mis)
            lg("  %s: only %d window SNPs matched the panel (min %d)", ens, n_mat, min_snps); next
        }
        setorder(cand, pos)
        fwrite(cand[, .(SNP = id_new, A1, A2, freq = signif(freq_exp, 6), b = beta_exp, se = se_exp, p = p_exp, N = N_exp)],
               sprintf("pwcoco_inputs/%s_exposure.ma", ens), sep = "\t")
        fwrite(cand[, .(SNP = id_new, A1, A2, freq = signif(freq_out, 6), b = beta_out, se = se_out, p = p_out, N = N_out)],
               sprintf("pwcoco_inputs/%s_outcome.ma", ens), sep = "\t")
        # rewrite the extracted .bim with the new ids (numeric chr, as PWCoCo requires)
        file.copy(bimf, paste0(bimf, ".orig"), overwrite = TRUE)
        fwrite(bim[, .(chr, id_new, cm, bp, a1, a2)], bimf, sep = "\t", col.names = FALSE)
        manifest[[ens]] <- data.table(ens = ens, chr = chr)
        status[[ens]] <- rec("ok", n_pan, n_mat, n_cmp, n_mis)
        lg("  %s: panel %d SNPs, matched %d of %d window SNPs (%d strand-complemented, %d allele mismatches)",
           ens, n_pan, n_mat, nrow(win), n_cmp, n_mis)
    }
    man <- if (length(manifest)) rbindlist(manifest) else data.table(ens = character(), chr = character())
    fwrite(man, "pwcoco_manifest.tsv", sep = "\t", col.names = FALSE)
    st <- if (length(status)) rbindlist(status) else
              data.table(Gene_ID = character(), align_status = character(), n_panel = integer(),
                         n_matched = integer(), n_complemented = integer(), n_allele_mismatch = integer())
    fwrite(st, "pwcoco_align_status.tsv", sep = "\t")
    lg("align: %d genes ready for PWCoCo", nrow(man))

# ===========================================================================
} else if (step == "parse") {
    gene_conv <- if (length(args) >= 2) tryCatch(readRDS(args[2]), error = function(e) NULL) else NULL
    gmap      <- if (length(args) >= 3) tryCatch(readRDS(args[3]), error = function(e) NULL) else NULL
    name_of <- function(ids) {
        if (is.null(gene_conv) || !all(c("id", "name") %in% colnames(gene_conv))) return(as.character(ids))
        nm <- as.character(gene_conv$name)[match(ens_core(ids), ens_core(gene_conv$id))]
        nm[is.na(nm)] <- as.character(ids)[is.na(nm)]
        nm
    }
    rd <- function(p) if (file.exists(p) && file.size(p) > 0) fread(p) else NULL
    prep  <- rd("pwcoco_prep_status.tsv")
    align <- rd("pwcoco_align_status.tsv")
    run   <- if (file.exists("pwcoco_run_status.tsv") && file.size("pwcoco_run_status.tsv") > 0)
                 fread("pwcoco_run_status.tsv", header = FALSE, col.names = c("Gene_ID", "run_status", "seconds")) else NULL
    ids <- unique(c(if (!is.null(prep)) prep$Gene_ID, if (!is.null(align)) align$Gene_ID, if (!is.null(run)) run$Gene_ID))
    rows <- list(); gst <- list()
    for (ens in ids) {
        p_st <- if (!is.null(prep))  prep[Gene_ID == ens]  else NULL
        a_st <- if (!is.null(align)) align[Gene_ID == ens] else NULL
        r_st <- if (!is.null(run))   run[Gene_ID == ens]   else NULL
        n_window <- if (!is.null(p_st) && nrow(p_st)) p_st$n_window[1] else NA_integer_
        n_panel <- n_matched <- n_cmp <- NA_integer_
        if (!is.null(a_st) && nrow(a_st)) { n_panel <- a_st$n_panel[1]; n_matched <- a_st$n_matched[1]; n_cmp <- a_st$n_complemented[1] }
        # ---- status before PWCoCo even ran ----
        status <- NA_character_
        if (!is.null(p_st) && nrow(p_st) && p_st$prep_status[1] != "ok") status <- p_st$prep_status[1]
        else if (!is.null(a_st) && nrow(a_st) && a_st$align_status[1] != "ok") status <- a_st$align_status[1]
        else if (is.null(r_st) || nrow(r_st) == 0) status <- "not_run"
        else if (r_st$run_status[1] != "ok") status <- r_st$run_status[1]
        cf <- sprintf("pwcoco_per_gene/%s.coloc", ens)
        lf <- sprintf("pwcoco_per_gene/%s.pwcoco_log.txt", ens)
        logtxt <- if (file.exists(lf)) readLines(lf, warn = FALSE) else character(0)
        res <- if (file.exists(cf) && file.size(cf) > 0) tryCatch(fread(cf), error = function(e) NULL) else NULL
        n_sig_exp <- n_sig_out <- NA_integer_
        m <- regmatches(logtxt, regexec("There are ([0-9]+) selected SNPs in the exposure dataset and ([0-9]+) in the outcome dataset", logtxt))
        m <- Filter(function(x) length(x) == 3, m)
        if (length(m)) { n_sig_exp <- as.integer(m[[1]][2]); n_sig_out <- as.integer(m[[1]][3]) }
        n_inc_exp <- n_inc_out <- NA_integer_
        cnt <- function(pat) { f <- list.files("pwcoco_per_gene", pattern = pat, full.names = TRUE); if (length(f)) max(0L, length(readLines(f[1], warn = FALSE)) - 1L) else NA_integer_ }
        if (is.na(status)) {
            n_inc_exp <- cnt(sprintf("^%s\\..*_exposure\\.ma\\.included$", ens))
            n_inc_out <- cnt(sprintf("^%s\\..*_outcome\\.ma\\.included$", ens))
            if (is.null(res) || nrow(res) == 0 || !"H4" %in% names(res)) {
                status <- if (any(grepl("no SNPs included in initial", logtxt))) "no_overlap" else
                          if (any(grepl("Included list of SNPs is empty", logtxt))) "no_panel_match" else "pwcoco_failed"
            } else if (any(grepl("already at or above threshold", logtxt))) {
                status <- "ok_initial_h4"     # marginal coloc H4 >= init_h4: PWCoCo stops before conditioning
            } else if (any(grepl("Both conditional analyses failed", logtxt))) {
                status <- "no_signal"         # no SNP passed p_cutoff in either trait: only the marginal row exists
                n_sig_exp <- 0L; n_sig_out <- 0L
            } else {
                status <- "ok"
                if (is.na(n_sig_exp)) n_sig_exp <- length(setdiff(unique(res$SNP1), "unconditioned"))
                if (is.na(n_sig_out)) n_sig_out <- length(setdiff(unique(res$SNP2), "unconditioned"))
            }
        }
        if (!is.null(res) && nrow(res) > 0 && "H4" %in% names(res)) {
            res <- as.data.table(res)
            res[, Gene_ID := ens]
            res[, n_signals_exposure := n_sig_exp]; res[, n_signals_outcome := n_sig_out]
            res[, pwcoco_status := status]
            rows[[ens]] <- res
            j <- which.max(num(res$H4))
            if (length(j) == 1) {
                bH4 <- num(res$H4)[j]; bH3 <- num(res$H3)[j]; bS1 <- as.character(res$SNP1)[j]; bS2 <- as.character(res$SNP2)[j]
            } else { bH4 <- bH3 <- NA_real_; bS1 <- bS2 <- NA_character_ }
            n_pairs <- nrow(res)
        } else {
            rows[[ens]] <- data.table(Dataset1 = NA_character_, Dataset2 = NA_character_, SNP1 = NA_character_, SNP2 = NA_character_,
                                      nsnps = NA_integer_, H0 = NA_real_, H1 = NA_real_, H2 = NA_real_, H3 = NA_real_, H4 = NA_real_,
                                      log_abf_all = NA_real_, Gene_ID = ens, n_signals_exposure = n_sig_exp,
                                      n_signals_outcome = n_sig_out, pwcoco_status = status)
            bH4 <- bH3 <- NA_real_; bS1 <- bS2 <- NA_character_; n_pairs <- 0L
        }
        gst[[ens]] <- data.table(Gene_ID = ens, pwcoco_status = status, n_window_snps = n_window, n_panel_snps = n_panel,
                                 n_matched = n_matched, n_strand_complemented = n_cmp,
                                 n_included_exposure = n_inc_exp, n_included_outcome = n_inc_out,
                                 n_signals_exposure = n_sig_exp, n_signals_outcome = n_sig_out, n_coloc_rows = n_pairs,
                                 best_H4 = bH4, best_H3 = bH3, best_SNP1 = bS1, best_SNP2 = bS2,
                                 seconds = if (!is.null(r_st) && nrow(r_st)) r_st$seconds[1] else NA_real_)
    }
    out <- if (length(rows)) rbindlist(rows, fill = TRUE) else
               data.table(Dataset1 = character(), Dataset2 = character(), SNP1 = character(), SNP2 = character(), nsnps = integer(),
                          H0 = numeric(), H1 = numeric(), H2 = numeric(), H3 = numeric(), H4 = numeric(), log_abf_all = numeric(),
                          Gene_ID = character(), n_signals_exposure = integer(), n_signals_outcome = integer(), pwcoco_status = character())
    out[, Gene_Name := name_of(Gene_ID)]
    if (!is.null(gmap) && all(c("exp_core", "map_method") %in% colnames(gmap)))
        out[, gene_map_method := as.character(gmap$map_method)[match(ens_core(Gene_ID), as.character(gmap$exp_core))]]
    setcolorder(out, c("Gene_ID", "Gene_Name"))
    g <- if (length(gst)) rbindlist(gst) else data.table(Gene_ID = character(), pwcoco_status = character())
    g[, Gene_Name := name_of(Gene_ID)]
    setcolorder(g, c("Gene_ID", "Gene_Name"))

    # ---- v2.2 decision metrics, per (exposure signal x outcome signal) pair ------
    # Same functions as the other coloc processes (bin/coloc_decision.R, base R only):
    # H3_plus_H4, H4_cond = H4/(H3+H4), H4_H3_odds, log2_H4_H3, dominant_hyp,
    # pass_rule (NF_coloc_rule) and coloc_class. Applied to every row, the marginal
    # "unconditioned" pair included; a gene passes when ANY pair passes (the PWCoCo
    # authors call a locus colocalised when any pair reaches H4 >= 0.8).
    nf_opt <- function(k, default) { v <- Sys.getenv(paste0("NF_", k), unset = NA); if (is.na(v) || !nzchar(v)) default else v }
    helper <- file.path(nf_opt("bin_dir", ""), "coloc_decision.R")
    if (!file.exists(helper)) helper <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "coloc_decision.R")
    if (!file.exists(helper)) helper <- Sys.which("coloc_decision.R")
    out_df <- as.data.frame(out)
    if (nzchar(helper) && file.exists(helper)) {
        source(helper)
        rule_str  <- nf_opt("coloc_rule", "H4 > 0.8")
        h4_strong <- num(nf_opt("coloc_h4_threshold", "0.8"))
        cond_thr  <- num(nf_opt("coloc_cond_threshold", "0.7"))
        power_min <- num(nf_opt("coloc_power_min", "0.5"))
        h3_strong <- num(nf_opt("coloc_h3_strong", "0.8"))
        out_df <- coloc_derive_metrics(out_df)
        out_df$pass_rule   <- coloc_apply_rule(out_df, rule_str)
        out_df$coloc_class <- coloc_classify(out_df, h4_strong = h4_strong, cond_thr = cond_thr,
                                             power_min = power_min, h3_strong = h3_strong)
        out_df$coloc_rule  <- rep(rule_str, nrow(out_df))
        # per-gene: metrics of the best pair (highest H4) + any-pair flags
        if (nrow(out_df) > 0) {
            o <- order(out_df$Gene_ID, -ifelse(is.finite(out_df$H4), out_df$H4, -Inf))
            b <- out_df[o, , drop = FALSE]; b <- b[!duplicated(b$Gene_ID), , drop = FALSE]
            gid <- as.character(out_df$Gene_ID)
            agg <- function(x, f) { r <- tapply(x, gid, f); unname(r[as.character(b$Gene_ID)]) }
            b$n_pairs_pass <- as.integer(agg(as.integer(out_df$pass_rule %in% TRUE), sum))
            b$any_pass     <- b$n_pairs_pass > 0
            bm <- data.frame(Gene_ID = b$Gene_ID, best_H4_cond = b$H4_cond, best_H3_plus_H4 = b$H3_plus_H4,
                             best_class = b$coloc_class, best_pass = b$pass_rule %in% TRUE,
                             n_pairs_pass = b$n_pairs_pass, any_pass = b$any_pass, stringsAsFactors = FALSE)
            g <- merge(g, as.data.table(bm), by = "Gene_ID", all.x = TRUE, sort = FALSE)
        }
        lg("parse: rule '%s' -> %d of %d rows pass; genes with any passing pair: %d; classes (best pair): %s",
           rule_str, sum(out_df$pass_rule %in% TRUE), nrow(out_df), sum(g$any_pass %in% TRUE),
           paste(sprintf("%s=%d", names(table(g$best_class)), as.integer(table(g$best_class))), collapse = " "))
    } else {
        lg("parse: bin/coloc_decision.R not found - no v2.2 metric columns (set NF_bin_dir)")
    }
    saveRDS(out_df, "pwcoco_results_summary.rds")
    fwrite(out_df, "pwcoco_results_summary.tsv", sep = "\t")
    fwrite(g, "pwcoco_gene_status.tsv", sep = "\t")
    lg("parse: %d genes, %d coloc rows; %d genes with best H4 > 0.8", nrow(g), sum(!is.na(out_df$H4)), sum(g$best_H4 > 0.8, na.rm = TRUE))
    if (nrow(g)) { tt <- table(g$pwcoco_status); lg("  status: %s", paste(sprintf("%s=%d", names(tt), as.integer(tt)), collapse = " ")) }
} else {
    stop("unknown step: ", step)
}
close(log_con)
