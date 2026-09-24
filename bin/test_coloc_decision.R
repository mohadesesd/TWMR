suppressMessages(library(coloc))
helper <- if (file.exists("coloc_decision.R")) "coloc_decision.R" else "bin/coloc_decision.R"
source(helper)   # run from bin/ or from the pipeline root
data(coloc_test_data)
attach(coloc_test_data)

rule_cond <- "H4/(H3+H4) > 0.7 & H3+H4 > 0.5"

# --- 1. real coloc.abf results on the package test data ----------------------
pairs <- list(D1_D2 = list(D1, D2), D3_D4 = list(D3, D4), D1_D3 = list(D1, D3), D2_D4 = list(D2, D4))
rows <- list()
for (nm in names(pairs)) {
    res <- suppressMessages(coloc.abf(pairs[[nm]][[1]], pairs[[nm]][[2]], p12 = 1e-5))
    sy <- as.list(res$summary)
    row <- data.frame(pair = nm, nSNP = sy$nsnps, H0 = sy$PP.H0.abf, H1 = sy$PP.H1.abf, H2 = sy$PP.H2.abf,
                      H3 = sy$PP.H3.abf, H4 = sy$PP.H4.abf, stringsAsFactors = FALSE)
    row <- cbind(row, coloc_abf_hits(res),
                 coloc_p12_profile(pairs[[nm]][[1]], pairs[[nm]][[2]], rule = rule_cond, res = res))
    rows[[nm]] <- row
}
tab <- do.call(rbind, rows)
tab <- coloc_derive_metrics(tab)
tab$pass_h4   <- coloc_apply_rule(tab, "H4 > 0.8")
tab$pass_cond <- coloc_apply_rule(tab, rule_cond)
tab$pass_odds <- coloc_apply_rule(tab, "H4 > 0.9 & H4/H3 > 3")
tab$coloc_class <- coloc_classify(tab)
print(tab[, c("pair","nSNP","H0","H1","H2","H3","H4","H3_plus_H4","H4_cond","H4_H3_odds","log2_H4_H3","dominant_hyp",
              "pass_h4","pass_cond","pass_odds","coloc_class")], digits = 3)
print(tab[, grep("p12|hit_|top_snp|cs95", names(tab))], digits = 3)

# --- 2. synthetic rows covering every class -----------------------------------
syn <- data.frame(
    Gene = c("strong","cond_ok","underpowered_H1","distinct_H3","cond_but_underpowered","borderline"),
    H0 = c(0.00, 0.02, 0.05, 0.00, 0.10, 0.02),
    H1 = c(0.01, 0.10, 0.85, 0.01, 0.75, 0.10),
    H2 = c(0.01, 0.05, 0.02, 0.01, 0.05, 0.05),
    H3 = c(0.08, 0.20, 0.02, 0.90, 0.02, 0.45),
    H4 = c(0.90, 0.63, 0.06, 0.08, 0.08, 0.38))
syn <- coloc_derive_metrics(syn)
syn$pass_h4   <- coloc_apply_rule(syn, "H4 > 0.8")
syn$pass_cond <- coloc_apply_rule(syn, rule_cond)
syn$class     <- coloc_classify(syn)
print(syn, digits = 3)
stopifnot(identical(syn$class, c("colocalized_strong","colocalized_conditional","underpowered",
                                 "distinct_signals","underpowered","inconclusive")))
stopifnot(identical(syn$pass_cond, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)))
# the ratio alone (no power guard) would 'pass' the underpowered H1 rows: show it
syn$pass_ratio_only <- coloc_apply_rule(syn, "H4/(H3+H4) > 0.7")
cat("\nratio-only passes:", paste(syn$Gene[syn$pass_ratio_only], collapse = ", "), "\n")
stopifnot(all(c("underpowered_H1","cond_but_underpowered") %in% syn$Gene[syn$pass_ratio_only]))

# --- 3. SuSiE-style column names ------------------------------------------------
ss <- data.frame(Gene_ID = "g", LD_type = "susie", PP.H0.abf = 0.01, PP.H1.abf = 0.02, PP.H2.abf = 0.02,
                 PP.H3.abf = 0.25, PP.H4.abf = 0.70)
ss <- coloc_derive_metrics(ss)
stopifnot(abs(ss$H4_cond - 0.70/0.95) < 1e-9)
stopifnot(coloc_apply_rule(ss, rule_cond), !coloc_apply_rule(ss, "H4 > 0.8"))
cat("SuSiE-style columns OK; class =", coloc_classify(ss), "\n")

# --- 4. FDR table (Reales et al. 2026) --------------------------------------------
set.seed(1); h4 <- c(runif(40, 0, 0.5), runif(30, 0.5, 0.85), runif(30, 0.85, 1))
print(coloc_fdr_table(h4), digits = 3)
cat("alpha for FDR<0.05:", coloc_alpha_for_fdr(h4, 0.05), "\n")

# --- 5. LD check -------------------------------------------------------------------
set.seed(2); n <- 50; R <- diag(n); R[1, 2] <- R[2, 1] <- 0.95; R[3, 4] <- R[4, 3] <- 0.5
ids <- paste0("rs", 1:n); pv <- runif(n); pv[2] <- 1e-9; pv[4] <- 1e-8
print(coloc_ld_check(R, ids, instrument_ids = c("rs1"), outcome_p = pv, top_n = 5, r2_thr = 0.8))
print(coloc_ld_check(R, ids, instrument_ids = c("rs3"), outcome_p = pv, top_n = 5, r2_thr = 0.8))

# --- 6. empty inputs ----------------------------------------------------------------
e <- coloc_derive_metrics(data.frame(H0=numeric(),H1=numeric(),H2=numeric(),H3=numeric(),H4=numeric()))
stopifnot(nrow(e) == 0, length(coloc_apply_rule(e, "H4 > 0.8")) == 0, length(coloc_classify(e)) == 0)
cat("\nALL CHECKS PASSED\n")
