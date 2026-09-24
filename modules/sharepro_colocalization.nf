process SHAREPRO_COLOCALIZATION {
    // NO container directive on purpose. This process needs BOTH R (to turn the
    // SuSiE .rds windows into SharePro's text inputs and to parse the results) and
    // Python 3 with numpy/scipy/pandas (to run SharePro itself). The old
    // 'python:3.9-slim' image has no R at all, so every Rscript call inside it
    // failed. No image in this pipeline ships both, so the process runs in the
    // host / module environment: the script loads a python module and builds an
    // isolated venv, and expects Rscript on PATH.
    // To run it in an image that has both, set one in nextflow.config, e.g.
    //   process { withName: SHAREPRO_COLOCALIZATION { container = '<image>' } }
    publishDir "${params.outdir}/${params.tissue}/colocalization/sharepro", mode: 'copy'

    maxRetries 1
    errorStrategy { params.coloc_error_strategy ?: 'ignore' }
    memory { params.sharepro_memory ?: '16 GB' }
    time   { params.sharepro_time   ?: '6h' }

    input:
    path window_files                 // window_<ENSG>.rds from SuSiE (beta/se/N both traits)
    path ld_files                     // ld_<ENSG>.rds signed-R matrices from SuSiE
    path coloc_summary                // classic coloc_results_summary.rds (for comparison)
    path susie_summary                // susie_coloc_results_summary.rds (for comparison)
    path pwcoco_summary               // pwcoco_results_summary.rds, or the empty NO_FILE_PWCOCO placeholder
    path gene_conversion_table        // GeneCode_ConversionTable.rds (id -> name)
    path gene_map                     // gene_release_map.rds (annotation only)
    val tissue
    val sample_size                   // exposure N of the ALL-PAIRS release
    val sharepro_k                    // max causal signals (--K)
    val sharepro_dir                  // existing SharePro_coloc checkout; "" / null -> clone
    val sharepro_commit               // git ref to check out when cloning (e.g. "main" or a SHA)
    val outdir

    output:
    path "sharepro_results_summary.rds", emit: sharepro_table
    path "sharepro_per_gene/*.sharepro.txt", optional: true, emit: per_gene
    path "sharepro_vs_others.tsv", optional: true, emit: comparison
    path "sharepro_log.txt", emit: log
    path "sharepro_provenance.txt", emit: provenance

    script:
    """
    set -e
    log=sharepro_log.txt
    prov=sharepro_provenance.txt
    echo "=== SharePro colocalization ===" | tee \$log
    # Provenance is a declared output, so start it immediately: even the early
    # exits below then leave a file behind saying how far the run got.
    echo "=== SharePro provenance ===" > \$prov
    echo "date_utc=\$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> \$prov
    echo "requested_sharepro_dir=${sharepro_dir}" >> \$prov
    echo "requested_sharepro_commit=${sharepro_commit}" >> \$prov
    echo "R=\$(Rscript -e 'cat(R.version.string)' 2>/dev/null || echo NA)" >> \$prov

    # ---- 0. Get a Python WITH scientific stack ----
    # The cluster profile only loads R, so load Python explicitly and build an
    # isolated venv with the needed packages. Fail loudly if the stack is missing
    # rather than silently producing 0 results.
    module load StdEnv/2023 python/3.11 2>>\$log || module load python 2>>\$log || true
    if ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: no python3 after 'module load python'. Check available modules with 'module spider python'." | tee -a \$log
        echo "status=no_python" >> \$prov
        # still write empty outputs so the process produces its declared files
        Rscript -e 'saveRDS(data.frame(Gene_ID=character(), Gene_Name=character(), sharepro_share=numeric(), sharepro_status=character()), "sharepro_results_summary.rds")'
        exit 1
    fi
    python3 -m venv sp_env 2>>\$log
    . sp_env/bin/activate
    export PIP_ROOT_USER_ACTION=ignore
    pip install --quiet --no-index numpy scipy pandas 2>>\$log \
        || pip install --quiet numpy scipy pandas 2>>\$log
    # verify the stack actually imports; if not, stop with a clear message
    if ! python3 -c "import numpy, scipy, pandas" 2>>\$log; then
        echo "ERROR: numpy/scipy/pandas not importable in venv. On Alliance clusters try:" | tee -a \$log
        echo "  module load scipy-stack   (bundles numpy/scipy/pandas)" | tee -a \$log
        echo "status=no_python_stack" >> \$prov
        Rscript -e 'saveRDS(data.frame(Gene_ID=character(), Gene_Name=character(), sharepro_share=numeric(), sharepro_status=character()), "sharepro_results_summary.rds")'
        exit 1
    fi
    PY=python3
    echo "Python stack OK: \$(python3 -c 'import pandas; print(pandas.__version__)')" | tee -a \$log
    {
      echo "python=\$(python3 -V 2>&1)"
      echo "numpy=\$(python3 -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo NA)"
      echo "scipy=\$(python3 -c 'import scipy; print(scipy.__version__)' 2>/dev/null || echo NA)"
      echo "pandas=\$(python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)"
    } >> \$prov

    # ---- 1. Get SharePro ----
    # A local checkout (--sharepro_dir) is used as-is, which is what a cluster with
    # no outbound network needs. Otherwise clone and pin the requested commit, and
    # record the SHA actually used so a run can be reproduced.
    SP_DIR="${sharepro_dir}"
    SP_COMMIT="${sharepro_commit}"
    SP_SOURCE=unknown
    SP_COMMIT_RESOLVED=unknown
    if [ -n "\$SP_DIR" ] && [ "\$SP_DIR" != "null" ] && [ -d "\$SP_DIR" ]; then
        SP="\$SP_DIR/src/SharePro/sharepro_coloc.py"
        SP_SOURCE=user_supplied
        SP_COMMIT_RESOLVED=\$(git -C "\$SP_DIR" rev-parse HEAD 2>/dev/null || echo unknown)
        echo "Using SharePro checkout supplied via --sharepro_dir: \$SP_DIR" | tee -a \$log
    else
        if [ -n "\$SP_DIR" ] && [ "\$SP_DIR" != "null" ]; then
            echo "WARNING: --sharepro_dir '\$SP_DIR' does not exist; falling back to a clone" | tee -a \$log
        fi
        if [ ! -d SharePro_coloc ]; then
            if command -v git >/dev/null 2>&1; then
                # full clone (not --depth 1): a shallow clone cannot check out an
                # arbitrary commit
                git clone https://github.com/zhwm/SharePro_coloc.git 2>>\$log && SP_SOURCE=git_clone \
                    || echo "ERROR: git clone failed" | tee -a \$log
            else
                \$PY -c "import urllib.request,zipfile,io; d=urllib.request.urlopen('https://github.com/zhwm/SharePro_coloc/archive/refs/heads/main.zip').read(); zipfile.ZipFile(io.BytesIO(d)).extractall('.'); import os; os.rename('SharePro_coloc-main','SharePro_coloc')" 2>>\$log \
                    && SP_SOURCE=zip_download_main \
                    || echo "ERROR: SharePro download failed" | tee -a \$log
            fi
        else
            SP_SOURCE=existing_checkout
        fi
        if [ "\$SP_SOURCE" = "git_clone" ] && [ -n "\$SP_COMMIT" ] && [ "\$SP_COMMIT" != "null" ]; then
            git -C SharePro_coloc checkout --quiet "\$SP_COMMIT" 2>>\$log \
                || echo "WARNING: could not check out '\$SP_COMMIT'; staying on the default branch" | tee -a \$log
        fi
        SP=SharePro_coloc/src/SharePro/sharepro_coloc.py
        SP_COMMIT_RESOLVED=\$(git -C SharePro_coloc rev-parse HEAD 2>/dev/null || echo unknown)
    fi
    {
      echo "sharepro_script=\$SP"
      echo "sharepro_source=\$SP_SOURCE"
      echo "sharepro_commit_resolved=\$SP_COMMIT_RESOLVED"
      echo "sharepro_K=${sharepro_k}"
    } >> \$prov
    echo "SharePro: \$SP (source=\$SP_SOURCE commit=\$SP_COMMIT_RESOLVED)" | tee -a \$log
    if [ ! -f "\$SP" ]; then echo "ERROR: SharePro script not found at \$SP" | tee -a \$log; fi

    # ---- 2. Prepare per-gene inputs from the SuSiE window + LD, run SharePro ----
    # An R step reads the .rds windows/LD and writes SharePro's SNP/BETA/SE/N files
    # and the .ld matrix; then Python runs SharePro per gene; then R parses results.
    Rscript - <<'RS' 2>>\$log
    suppressMessages({ library(data.table) })
    log_con <- file("sharepro_log.txt", open = "a")
    lg <- function(...) { m <- sprintf(...); cat(m, "\\n"); writeLines(m, log_con) }

    wf <- list.files(".", pattern = "^window_.*\\\\.rds\$")
    ens_of <- function(f) sub("^window_(.*)\\\\.rds\$", "\\\\1", f)
    dir.create("sharepro_inputs", showWarnings = FALSE)
    dir.create("sharepro_per_gene", showWarnings = FALSE)

    # Every gene that has a window is accounted for here, so a blank share in the
    # summary always has a stated reason: no_ld / bad_window / (later) failed.
    prep_status <- list()
    manifest <- list()
    n_ok <- 0
    for (f in wf) {
        ens <- ens_of(f)
        ldf <- sprintf("ld_%s.rds", ens)
        if (!file.exists(ldf)) {
            lg("  %s: no LD file, skip", ens)
            prep_status[[ens]] <- data.frame(Gene_ID=ens, sharepro_status="no_ld", stringsAsFactors=FALSE)
            next
        }
        win <- tryCatch(as.data.frame(readRDS(f)), error=function(e) NULL)
        ld  <- tryCatch(as.matrix(readRDS(ldf)),    error=function(e) NULL)
        if (is.null(win) || is.null(ld) || nrow(win) < 3) {
            lg("  %s: bad window/LD", ens)
            prep_status[[ens]] <- data.frame(Gene_ID=ens, sharepro_status="bad_window", stringsAsFactors=FALSE)
            next
        }
        # align window rows to LD order by rsid
        rs <- if ("rsid" %in% names(win)) win\$rsid else rownames(ld)
        keep <- rs %in% rownames(ld)
        win <- win[keep, ]; rs <- rs[keep]
        if (length(rs) < 3) {
            lg("  %s: fewer than 3 SNPs shared between window and LD", ens)
            prep_status[[ens]] <- data.frame(Gene_ID=ens, sharepro_status="bad_window", stringsAsFactors=FALSE)
            next
        }
        ld  <- ld[rs, rs, drop=FALSE]
        n_exp <- suppressWarnings(as.numeric("${sample_size}"))
        # exposure and outcome N
        Nexp <- if (is.finite(n_exp) && n_exp>0) n_exp else suppressWarnings(max(win\$samplesize.exposure, na.rm=TRUE))
        Nout <- suppressWarnings(mean(win\$samplesize.outcome, na.rm=TRUE))
        if (!is.finite(Nexp) || Nexp<=0) Nexp <- 500
        if (!is.finite(Nout) || Nout<=0) Nout <- 100000
        # SharePro GWAS files: SNP BETA SE N  (one per trait)
        exp_tab <- data.frame(SNP=rs, BETA=win\$beta.exposure, SE=win\$se.exposure, N=round(Nexp))
        out_tab <- data.frame(SNP=rs, BETA=win\$beta.outcome,  SE=win\$se.outcome,  N=round(Nout))
        fe <- sprintf("sharepro_inputs/%s_exposure.txt", ens)
        fo <- sprintf("sharepro_inputs/%s_outcome.txt", ens)
        fl <- sprintf("sharepro_inputs/%s.ld", ens)
        write.table(exp_tab, fe, row.names=FALSE, quote=FALSE, sep="\\t")
        write.table(out_tab, fo, row.names=FALSE, quote=FALSE, sep="\\t")
        # LD as plain correlation matrix, whitespace-delimited, no header/rownames
        write.table(format(ld, digits=6), fl, row.names=FALSE, col.names=FALSE, quote=FALSE)
        manifest[[ens]] <- data.frame(ens=ens, fe=fe, fo=fo, fl=fl, stringsAsFactors=FALSE)
        prep_status[[ens]] <- data.frame(Gene_ID=ens, sharepro_status="ok", stringsAsFactors=FALSE)
        n_ok <- n_ok + 1
    }
    if (length(manifest) > 0)
        write.table(do.call(rbind, manifest), "sharepro_manifest.txt",
                    sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    ps <- if (length(prep_status) > 0) do.call(rbind, prep_status) else
              data.frame(Gene_ID=character(), sharepro_status=character(), stringsAsFactors=FALSE)
    write.table(ps, "sharepro_prep_status.tsv", sep="\\t", row.names=FALSE, quote=FALSE)
    lg("Prepared SharePro inputs for %d of %d genes", n_ok, length(wf))
    close(log_con)
RS

    # ---- 3. Run SharePro per gene ----
    mkdir -p sharepro_per_gene
    : > sharepro_run_status.tsv
    if [ -f sharepro_manifest.txt ]; then
        while IFS=\$'\\t' read -r ens fe fo fl; do
            echo "  running SharePro: \$ens" >> \$log
            if \$PY "\$SP" --z "\$fe" "\$fo" --ld "\$fl" \
                   --save "sharepro_per_gene/\$ens" --K ${sharepro_k} >> \$log 2>&1; then
                printf '%s\\t%s\\n' "\$ens" "ok" >> sharepro_run_status.tsv
            else
                echo "  SharePro failed for \$ens" >> \$log
                printf '%s\\t%s\\n' "\$ens" "failed" >> sharepro_run_status.tsv
            fi
        done < sharepro_manifest.txt
    fi

    # ---- 4. Parse results + compare to classic/SuSiE ----
    Rscript - <<'RS' 2>>\$log
    suppressMessages({ library(data.table) })
    log_con <- file("sharepro_log.txt", open = "a")
    lg <- function(...) { m <- sprintf(...); cat(m, "\\n"); writeLines(m, log_con) }

    # ---- shares from whatever SharePro actually wrote ----
    res_files <- list.files("sharepro_per_gene", pattern="\\\\.sharepro\\\\.txt\$", full.names=TRUE)
    share_map <- setNames(rep(NA_real_, length(res_files)),
                          sub("\\\\.sharepro\\\\.txt\$", "", basename(res_files)))
    for (rf in res_files) {
        ens <- sub("\\\\.sharepro\\\\.txt\$", "", basename(rf))
        d <- tryCatch(fread(rf), error=function(e) NULL)
        # 'share' column holds colocalization probability per effect group; take the max
        if (!is.null(d) && "share" %in% names(d) && nrow(d) > 0) {
            v <- suppressWarnings(max(as.numeric(d\$share), na.rm=TRUE))
            share_map[[ens]] <- if (is.finite(v)) v else NA_real_
        }
    }
    # SharePro has no H3: its per-effect-group 'share' is P(shared causal variant).
    # The closest analogue to coloc's H3 (distinct causal variants) is 1 - share of
    # the best effect group, reported as sharepro_H3_equiv.
    h3eq_map <- 1 - share_map

    # ---- one row per gene, with WHY a share is blank ----
    #   ok          - SharePro ran and produced a share
    #   failed      - inputs were prepared but SharePro errored or wrote no share
    #   no_ld       - SuSiE emitted no ld_<ENSG>.rds for that window
    #   bad_window  - window/LD unreadable, <3 SNPs, or no overlap between them
    rd <- function(p, hdr) {
        if (!file.exists(p) || file.size(p) == 0) return(NULL)
        tryCatch(as.data.frame(fread(p, header=hdr,
                 col.names=c("Gene_ID","sharepro_status"))), error=function(e) NULL)
    }
    runst <- rd("sharepro_run_status.tsv", FALSE)
    prep  <- rd("sharepro_prep_status.tsv", TRUE)
    st <- rbind(runst, prep)                       # run status wins over prep status
    if (!is.null(st)) st <- st[!duplicated(as.character(st\$Gene_ID)), , drop=FALSE]

    ids <- unique(c(if (!is.null(st)) as.character(st\$Gene_ID) else character(0),
                    names(share_map)))
    if (length(ids) == 0) {
        sp <- data.frame(Gene_ID=character(), Gene_Name=character(),
                         sharepro_share=numeric(), sharepro_H3_equiv=numeric(),
                         sharepro_status=character(), stringsAsFactors=FALSE)
    } else {
        status <- setNames(rep(NA_character_, length(ids)), ids)
        if (!is.null(st)) status[as.character(st\$Gene_ID)] <- as.character(st\$sharepro_status)
        share <- unname(share_map[ids])
        status[is.na(status)] <- "ok"
        status[is.na(share) & status == "ok"] <- "failed"
        sp <- data.frame(Gene_ID=ids, sharepro_share=share,
                         sharepro_H3_equiv=unname(h3eq_map[ids]),
                         sharepro_status=unname(status), stringsAsFactors=FALSE)
    }

    # ---- map ENSG -> gene name ----
    gene_conv <- tryCatch(readRDS("${gene_conversion_table}"), error=function(e) NULL)
    name_of <- function(ids) {
        if (is.null(gene_conv) || !all(c("id","name") %in% colnames(gene_conv))) return(as.character(ids))
        e  <- sub(".*?(ENSG[0-9]+).*", "\\\\1", as.character(ids))
        ci <- sub("(ENSG[0-9]+).*", "\\\\1", as.character(gene_conv\$id))  # strip version on map side too
        nm <- gene_conv\$name[match(e, ci)]
        nm[is.na(nm)] <- as.character(ids)[is.na(nm)]
        nm
    }
    sp\$Gene_Name <- name_of(sp\$Gene_ID)
    sp <- sp[, c("Gene_ID","Gene_Name","sharepro_share","sharepro_H3_equiv","sharepro_status")]
    # Evidence tiers of the SharePro paper (Zhang et al. 2024 Bioinformatics) and of the
    # Zhang et al. 2026 GPB benchmark: > 0.8 strong, > 0.5 suggestive, < 0.2 against.
    sp\$sharepro_tier <- ifelse(is.na(sp\$sharepro_share), NA_character_,
                        ifelse(sp\$sharepro_share > 0.8, "strong",
                        ifelse(sp\$sharepro_share > 0.5, "suggestive",
                        ifelse(sp\$sharepro_share < 0.2, "against", "inconclusive"))))
    saveRDS(sp, "sharepro_results_summary.rds")
    lg("SharePro results: %d genes; share > 0.8 (strong): %d; > 0.5 (suggestive): %d; < 0.2 (against): %d", nrow(sp),
       sum(sp\$sharepro_share > 0.8, na.rm=TRUE), sum(sp\$sharepro_share > 0.5, na.rm=TRUE), sum(sp\$sharepro_share < 0.2, na.rm=TRUE))
    if (nrow(sp) > 0) {
        stt <- table(sp\$sharepro_status, useNA="ifany")
        lg("  status: %s", paste(sprintf("%s=%d", names(stt), as.integer(stt)), collapse=" "))
    }

    # comparison table vs classic + SuSiE
    classic <- tryCatch(as.data.table(readRDS("${coloc_summary}")), error=function(e) NULL)
    susie   <- tryCatch(as.data.table(readRDS("${susie_summary}")),  error=function(e) NULL)
    cmp <- as.data.table(sp)
    # normalize gene id (strip version) for joins
    cmp\$ens_core <- sub(".*?(ENSG[0-9]+).*", "\\\\1", cmp\$Gene_ID)
    # annotate how each gene's id was bridged onto the all-pairs release
    gmap <- tryCatch(readRDS("${gene_map}"), error=function(e) NULL)
    gm_method <- function(cores) {
        if (is.null(gmap) || !all(c("exp_core","map_method") %in% colnames(gmap)))
            return(rep(NA_character_, length(cores)))
        as.character(gmap\$map_method)[match(cores, as.character(gmap\$exp_core))]
    }
    cmp\$gene_map_method <- gm_method(cmp\$ens_core)
    if (!is.null(classic) && nrow(classic) > 0) {
        ch <- if ("H4" %in% names(classic)) "H4" else if ("PP.H4.abf" %in% names(classic)) "PP.H4.abf" else NA
        if (!is.na(ch)) {
            classic[, ens_core := sub(".*?(ENSG[0-9]+).*", "\\\\1", Gene_ID)]
            c3 <- if ("H3" %in% names(classic)) "H3" else if ("PP.H3.abf" %in% names(classic)) "PP.H3.abf" else NA
            cc <- classic[, {
                h4 <- suppressWarnings(as.numeric(get(ch)))
                h3 <- if (is.na(c3)) rep(NA_real_, length(h4)) else suppressWarnings(as.numeric(get(c3)))
                j  <- if (all(is.na(h4))) NA_integer_ else which.max(h4)
                .(classic_H4 = if (is.na(j)) NA_real_ else h4[j],
                  classic_H3 = if (is.na(j)) NA_real_ else h3[j])
            }, by=ens_core]
            # v2.2 decision metrics of the classic table, when present
            if (all(c("H4_cond","H3_plus_H4","coloc_class","pass_rule") %in% names(classic))) {
                cc2 <- classic[, .(classic_H4_cond = as.numeric(H4_cond)[1], classic_H3_plus_H4 = as.numeric(H3_plus_H4)[1],
                                   classic_class = as.character(coloc_class)[1], classic_pass = as.logical(pass_rule)[1]), by=ens_core]
                cc <- merge(cc, cc2, by="ens_core", all.x=TRUE)
            }
            cmp <- merge(cmp, cc, by="ens_core", all.x=TRUE)
        }
    }
    if (!is.null(susie) && nrow(susie) > 0 && "PP.H4.abf" %in% names(susie)) {
        susie[, ens_core := sub(".*?(ENSG[0-9]+).*", "\\\\1", Gene_ID)]
        lt <- if ("LD_type" %in% names(susie)) susie\$LD_type else rep("susie", nrow(susie))
        susie[, LDt := lt]
        # class / rule of each gene's best REAL coloc.susie pair (v2.2 columns), when present
        if (all(c("coloc_class","pass_rule","H4_cond") %in% names(susie)) && any(lt == "susie")) {
            sr <- susie[LDt == "susie"]
            setorder(sr, ens_core, -PP.H4.abf)
            sb <- sr[!duplicated(ens_core), .(ens_core, susie_best_class = as.character(coloc_class), susie_best_pass = as.logical(pass_rule),
                                              susie_best_H4_cond = as.numeric(H4_cond))]
            sa <- sr[, .(susie_any_pass = any(pass_rule %in% TRUE)), by=ens_core]
            cmp <- merge(cmp, merge(sb, sa, by="ens_core"), by="ens_core", all.x=TRUE)
        }
        if (!"PP.H3.abf" %in% names(susie)) susie[, PP.H3.abf := NA_real_]
        # H4 and H3 are taken from the SAME row (the credible-set pair with the
        # highest H4), so the two numbers describe one hypothesis test.
        ss <- susie[, {
            h4 <- suppressWarnings(as.numeric(PP.H4.abf))
            h3 <- suppressWarnings(as.numeric(PP.H3.abf))
            j  <- if (all(is.na(h4))) NA_integer_ else which.max(h4)
            .(H4 = if (is.na(j)) NA_real_ else h4[j],
              H3 = if (is.na(j)) NA_real_ else h3[j])
        }, by=.(ens_core, LDt)]
        # LD_type values are "susie", "susie_bf" and "abf_fallback"; dcast yields
        # H4_<type> / H3_<type>, renamed to <type>_H4 / <type>_H3 so the meaning is
        # explicit ("susie" = coloc.susie, "susie_bf" = credible sets x single-signal
        # ABF, "abf_fallback" = coloc.abf on the same window).
        ss <- dcast(ss, ens_core ~ LDt, value.var=c("H4","H3"))
        vc <- setdiff(names(ss), "ens_core")
        if (length(vc) > 0) setnames(ss, vc, sub("^(H[34])_(.*)\$", "\\\\2_\\\\1", vc))
        cmp <- merge(cmp, ss, by="ens_core", all.x=TRUE)
    }
    # ---- PWCoCo (pair-wise conditional + coloc.abf) ----
    # pwcoco_results_summary.rds holds one row per (exposure signal x outcome
    # signal) pair, the marginal "unconditioned" pair included. The gene's
    # PWCoCo answer is the pair with the highest H4 (the authors call a locus
    # colocalised when any pair reaches H4 >= 0.8); H3 comes from that same row.
    # pwcoco_signals = "<n exposure signals>x<n outcome signals>" from the
    # stepwise selection; pwcoco_status says why a value is blank (tool_missing_*,
    # no_panel_overlap, no_signal, ok_initial_h4 = stopped after the marginal
    # coloc because H4 was already >= init_h4, ...).
    pwf <- "${pwcoco_summary}"
    pw <- if (!grepl("^NO_FILE", basename(pwf)) && file.exists(pwf) && file.size(pwf) > 0)
              tryCatch(as.data.table(readRDS(pwf)), error=function(e) NULL) else NULL
    if (!is.null(pw) && nrow(pw) > 0 && all(c("Gene_ID","H4","H3","pwcoco_status") %in% names(pw))) {
        pw[, ens_core := sub(".*?(ENSG[0-9]+).*", "\\\\1", Gene_ID)]
        for (cc in c("SNP1","SNP2")) if (!cc %in% names(pw)) pw[, (cc) := NA_character_]
        for (cc in c("n_signals_exposure","n_signals_outcome")) if (!cc %in% names(pw)) pw[, (cc) := NA_integer_]
        has_v22 <- all(c("H4_cond","H3_plus_H4","coloc_class","pass_rule") %in% names(pw))
        pwb <- pw[, {
            h4 <- suppressWarnings(as.numeric(H4)); h3 <- suppressWarnings(as.numeric(H3))
            j  <- if (all(is.na(h4))) NA_integer_ else which.max(h4)
            ne <- suppressWarnings(as.integer(n_signals_exposure[1])); no <- suppressWarnings(as.integer(n_signals_outcome[1]))
            # v2.2 metrics of that same best pair, when the PWCoCo table carries them
            .(pwcoco_H4 = if (is.na(j)) NA_real_ else h4[j],
              pwcoco_H3 = if (is.na(j)) NA_real_ else h3[j],
              pwcoco_H4_cond = if (is.na(j) || !has_v22) NA_real_ else suppressWarnings(as.numeric(H4_cond))[j],
              pwcoco_H3_plus_H4 = if (is.na(j) || !has_v22) NA_real_ else suppressWarnings(as.numeric(H3_plus_H4))[j],
              pwcoco_class = if (is.na(j) || !has_v22) NA_character_ else as.character(coloc_class)[j],
              pwcoco_best_pass = if (is.na(j) || !has_v22) NA else as.logical(pass_rule)[j],
              pwcoco_any_pass = if (!has_v22) NA else any(as.logical(pass_rule) %in% TRUE),
              pwcoco_signals = if (is.na(ne) && is.na(no)) NA_character_ else sprintf("%sx%s", ne, no),
              pwcoco_best_pair = if (is.na(j)) NA_character_ else paste(SNP1[j], SNP2[j], sep="|"),
              pwcoco_status = as.character(pwcoco_status[1]))
        }, by=ens_core]
        cmp <- merge(cmp, pwb, by="ens_core", all.x=TRUE)
        lg("PWCoCo columns added: %d genes with a value, %d with H4 > 0.8, %d passing the rule (any pair)",
           sum(!is.na(cmp\$pwcoco_H4)), sum(cmp\$pwcoco_H4 > 0.8, na.rm=TRUE), sum(cmp\$pwcoco_any_pass %in% TRUE))
    } else {
        lg("PWCoCo table not supplied (or empty): no pwcoco_* columns")
    }
    # order: each method's H4 immediately followed by its H3, PWCoCo details last
    ord <- unlist(lapply(c("susie","susie_bf","abf_fallback","pwcoco"),
                         function(m) intersect(c(paste0(m,"_H4"), paste0(m,"_H3")), names(cmp))))
    tail_cols <- intersect(c("pwcoco_H4_cond","pwcoco_H3_plus_H4","pwcoco_class","pwcoco_best_pass","pwcoco_any_pass",
                             "pwcoco_signals","pwcoco_best_pair","pwcoco_status"), names(cmp))
    lead <- setdiff(names(cmp), c(ord, tail_cols))
    setcolorder(cmp, c(lead, ord, tail_cols))
    fwrite(cmp, "sharepro_vs_others.tsv", sep="\\t")
    lg("Wrote comparison table: %d rows", nrow(cmp))
    close(log_con)
RS
    echo "SharePro done." >> \$log
    """
}
