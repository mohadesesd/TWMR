process PWCOCO_COLOCALIZATION {
    // PWCoCo = pair-wise conditional analysis (GCTA-COJO style) + coloc.abf on
    // every (exposure signal x outcome signal) pair. Robinson et al. 2022,
    // https://github.com/jwr-git/pwcoco
    //
    // NO container directive on purpose (same reasoning as SharePro): it needs
    // host Rscript + data.table, a PLINK binary (1.9 or 2) and the compiled
    // pwcoco executable. Build pwcoco once on the cluster:
    //   module load StdEnv/2023 gcc cmake
    //   git clone https://github.com/jwr-git/pwcoco.git && cd pwcoco
    //   mkdir build && cd build && cmake .. && make          # -> build/pwcoco
    // then point --pwcoco_bin at that binary.
    //
    // Per gene it (1) cuts the cis-window out of the PLINK panel, (2) rewrites
    // the panel ids to chr:bp:A1:A2 and writes the two COJO-format .ma files from
    // the SAME harmonised window SuSiE/SharePro used, (3) runs pwcoco, and
    // (4) collects every H0-H4 row plus a per-gene status.
    //
    // Inputs the panel must satisfy (checked in the log, not assumed):
    //   * PLINK 1 bed/bim/fam; "{chr}" in the prefix is replaced by the chromosome
    //   * same population as the GWAS; >= 4,000 people recommended by COJO
    //   * pwcoco_match = "position": panel on GRCh38 (like the GTEx/GWAS data)
    //     pwcoco_match = "rsid":     any build, but the .bim must carry rsIDs
    publishDir "${params.outdir}/${params.tissue}/colocalization/pwcoco", mode: 'copy'

    maxRetries 1
    errorStrategy { params.pwcoco_error_strategy ?: 'ignore' }
    memory { params.pwcoco_memory ?: '16 GB' }
    cpus   { params.pwcoco_cpus   ?: 4 }
    time   { params.pwcoco_time   ?: '24h' }

    input:
    path window_files                 // window_<ENSG>.rds from SuSiE (beta/se/N/alleles, both traits)
    path gene_conversion_table        // GeneCode_ConversionTable.rds (id -> name)
    path gene_map                     // gene_release_map.rds (annotation only)
    val tissue
    val sample_size                   // exposure N of the ALL-PAIRS release
    val coloc_p12                     // shared prior, same value the other methods use
    val bfile                         // PLINK prefix (no extension); may contain {chr}
    val outdir

    output:
    path "pwcoco_results_summary.rds", emit: pwcoco_table
    path "pwcoco_results_summary.tsv", optional: true, emit: pwcoco_tsv   // v2.2: + H4_cond, H3_plus_H4, pass_rule, coloc_class
    path "pwcoco_gene_status.tsv", emit: gene_status
    path "pwcoco_per_gene/*", optional: true, emit: per_gene
    path "pwcoco_inputs/*.ma", optional: true, emit: ma_files
    path "pwcoco_run_log.txt", emit: log
    path "pwcoco_provenance.txt", emit: provenance

    script:
    def p_cut1   = params.pwcoco_p_cutoff1      ?: 5e-8
    def p_cut2   = params.pwcoco_p_cutoff2      ?: 5e-8
    def maf      = params.pwcoco_maf            != null ? params.pwcoco_maf : 0.01
    def freq_thr = params.pwcoco_freq_threshold ?: 0.2
    def collin   = params.pwcoco_collinear      ?: 0.9
    def init_h4  = params.pwcoco_init_h4        != null ? params.pwcoco_init_h4 : 80
    def match    = params.pwcoco_match          ?: 'position'
    def min_snps = params.pwcoco_min_snps       ?: 10
    def verbose  = params.pwcoco_verbose ? '--verbose' : ''
    def extra    = params.pwcoco_extra_args     ?: ''
    def mem_mb   = task.memory ? (task.memory.toMega() as int) : 8000
    def plink_mb = Math.max(2000, mem_mb - 2048)
    """
    log=pwcoco_run_log.txt
    prov=pwcoco_provenance.txt
    echo "=== PWCoCo colocalization ===" | tee \$log
    echo "=== PWCoCo provenance ===" > \$prov
    echo "date_utc=\$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> \$prov
    echo "bfile=${bfile}" >> \$prov
    echo "match_mode=${match}" >> \$prov
    echo "p_cutoff1=${p_cut1} p_cutoff2=${p_cut2} maf=${maf} freq_threshold=${freq_thr} collinear=${collin} init_h4=${init_h4}" >> \$prov
    echo "coloc_pp=1e-4 1e-4 ${coloc_p12}" >> \$prov
    echo "R=\$(Rscript -e 'cat(R.version.string)' 2>/dev/null || echo NA)" >> \$prov
    mkdir -p pwcoco_inputs pwcoco_panel pwcoco_per_gene

    # Every exit path below leaves the declared outputs behind with a stated
    # reason, so a missing tool shows up as a status in the comparison table
    # instead of silently removing PWCoCo from it.
    bail() {
        echo "ERROR: \$1" | tee -a \$log
        echo "status=\$2" >> \$prov
        printf 'Gene_ID\\tpwcoco_status\\n' > pwcoco_gene_status.tsv
        for f in window_*.rds; do [ -e "\$f" ] || continue; e=\${f#window_}; printf '%s\\t%s\\n' "\${e%.rds}" "\$2" >> pwcoco_gene_status.tsv; done
        # one NA row per gene, carrying the reason, so the comparison table shows it
        Rscript -e 'g <- read.delim("pwcoco_gene_status.tsv", colClasses="character"); n <- nrow(g);
                    saveRDS(data.frame(Gene_ID=g\$Gene_ID, Gene_Name=g\$Gene_ID, SNP1=rep(NA_character_,n), SNP2=rep(NA_character_,n),
                            H3=rep(NA_real_,n), H4=rep(NA_real_,n), n_signals_exposure=rep(NA_integer_,n), n_signals_outcome=rep(NA_integer_,n),
                            pwcoco_status=g\$pwcoco_status, stringsAsFactors=FALSE), "pwcoco_results_summary.rds")'
        exit 0
    }

    # ---- 0. Tools ----
    PLINK="${params.plink_bin ?: 'plink'}"
    if ! command -v "\$PLINK" >/dev/null 2>&1; then
        module load plink/1.9b_6.21-x86_64 2>>\$log || module load plink 2>>\$log || true
        for cand in plink plink1.9 plink2; do
            if command -v "\$cand" >/dev/null 2>&1; then PLINK="\$cand"; break; fi
        done
    fi
    command -v "\$PLINK" >/dev/null 2>&1 || bail "PLINK not found (--plink_bin '${params.plink_bin ?: 'plink'}'; tried 'module load plink')" tool_missing_plink
    PLINK_VER=\$("\$PLINK" --version 2>/dev/null | head -1)
    echo "plink=\$PLINK (\$PLINK_VER)" | tee -a \$prov \$log
    if echo "\$PLINK_VER" | grep -q "v2"; then
        PLINK_OPTS="--silent --threads ${task.cpus} --memory ${plink_mb}"
    else
        PLINK_OPTS="--allow-no-sex --silent --threads ${task.cpus} --memory ${plink_mb}"
    fi

    PWCOCO="${params.pwcoco_bin ?: 'pwcoco'}"
    if [ ! -x "\$PWCOCO" ]; then
        command -v "\$PWCOCO" >/dev/null 2>&1 && PWCOCO=\$(command -v "\$PWCOCO") || bail "pwcoco executable not found at '\$PWCOCO' (set --pwcoco_bin to the compiled binary)" tool_missing_pwcoco
    fi
    ( mkdir -p ver_tmp && cd ver_tmp && "\$PWCOCO" --version 2>/dev/null | grep -o 'Version: [0-9.]*' | head -1 ) > pwcoco_version.txt || true
    echo "pwcoco=\$PWCOCO (\$(cat pwcoco_version.txt 2>/dev/null || echo unknown))" | tee -a \$prov \$log
    rm -rf ver_tmp
    Rscript -e 'suppressMessages(library(data.table))' 2>>\$log || bail "R package data.table not available to Rscript" tool_missing_r

    # ---- 1. Windows -> regions + per-gene tables ----
    Rscript "${projectDir}/bin/pwcoco_colocalization.R" prep "${sample_size}" "${match}" 2>>\$log \
        || bail "prep step failed (see log)" prep_failed
    [ -s pwcoco_regions.tsv ] || bail "no usable SuSiE window (window_*.rds) to run PWCoCo on" no_windows

    # ---- 2. Cut each cis-window out of the PLINK panel ----
    BFILE_PAT='${bfile}'
    : > pwcoco_extract_status.tsv
    while IFS=\$'\\t' read -r ens chr start end; do
        prefix=\$(printf '%s' "\$BFILE_PAT" | sed "s/{chr}/\$chr/g")
        if [ ! -f "\$prefix.bed" ] || [ ! -f "\$prefix.bim" ] || [ ! -f "\$prefix.fam" ]; then
            echo "  \$ens: panel files \$prefix.{bed,bim,fam} not found" >> \$log
            printf '%s\\tpanel_missing\\t%s\\n' "\$ens" "\$chr" >> pwcoco_extract_status.tsv
            continue
        fi
        if [ "${match}" = "rsid" ]; then
            sel="--extract pwcoco_inputs/\$ens.rsids.txt"
        else
            sel="--from-bp \$start --to-bp \$end"
        fi
        if "\$PLINK" --bfile "\$prefix" --chr "\$chr" \$sel --make-bed --out "pwcoco_panel/\$ens" \$PLINK_OPTS \
                > "pwcoco_panel/\$ens.extract.stdout" 2>&1 && [ -s "pwcoco_panel/\$ens.bim" ]; then
            "\$PLINK" --bfile "pwcoco_panel/\$ens" --freq --out "pwcoco_panel/\$ens" \$PLINK_OPTS \
                > "pwcoco_panel/\$ens.freq.stdout" 2>&1 || echo "  \$ens: plink --freq failed (frequencies fall back to GTEx eaf)" >> \$log
            printf '%s\\textracted\\t%s\\n' "\$ens" "\$chr" >> pwcoco_extract_status.tsv
        else
            # plink 1.9 says "All variants excluded", plink 2 "No variants remaining"
            if grep -qi "variants excluded\\|No variants remain\\|no variants remain" "pwcoco_panel/\$ens.extract.stdout" "pwcoco_panel/\$ens.log" 2>/dev/null; then
                echo "  \$ens: the panel has no variant inside chr\$chr:\$start-\$end" >> \$log
                printf '%s\\tno_panel_overlap\\t%s\\n' "\$ens" "\$chr" >> pwcoco_extract_status.tsv
            else
                echo "  \$ens: plink extraction failed:" >> \$log
                tail -n 5 "pwcoco_panel/\$ens.extract.stdout" >> \$log 2>/dev/null || true
                printf '%s\\tplink_failed\\t%s\\n' "\$ens" "\$chr" >> pwcoco_extract_status.tsv
            fi
        fi
    done < pwcoco_regions.tsv
    echo "panel extraction: \$(grep -c 'extracted' pwcoco_extract_status.tsv) of \$(wc -l < pwcoco_regions.tsv) windows" | tee -a \$log

    # ---- 3. Align alleles, write COJO-format .ma files, rewrite panel ids ----
    Rscript "${projectDir}/bin/pwcoco_colocalization.R" align "${match}" "${min_snps}" 2>>\$log \
        || bail "align step failed (see log)" align_failed

    # ---- 4. Run PWCoCo per gene ----
    : > pwcoco_run_status.tsv
    if [ -s pwcoco_manifest.tsv ]; then
        while IFS=\$'\\t' read -r ens chr; do
            t0=\$(date +%s)
            rm -f "pwcoco_per_gene/\$ens.coloc"           # pwcoco APPENDS to an existing .coloc
            if "\$PWCOCO" --bfile "pwcoco_panel/\$ens" \
                    --sum_stats1 "pwcoco_inputs/\${ens}_exposure.ma" \
                    --sum_stats2 "pwcoco_inputs/\${ens}_outcome.ma" \
                    --out "pwcoco_per_gene/\$ens" --log "pwcoco_per_gene/\$ens.pwcoco_log" \
                    --chr "\$chr" --p_cutoff1 ${p_cut1} --p_cutoff2 ${p_cut2} \
                    --maf ${maf} --freq_threshold ${freq_thr} --collinear ${collin} \
                    --init_h4 ${init_h4} --coloc_pp 1e-4 1e-4 ${coloc_p12} \
                    --threads ${task.cpus} ${verbose} ${extra} > /dev/null 2>> \$log; then
                st=ok
            else
                st=pwcoco_failed
                echo "  \$ens: pwcoco exited with an error" >> \$log
            fi
            printf '%s\\t%s\\t%s\\n' "\$ens" "\$st" "\$(( \$(date +%s) - t0 ))" >> pwcoco_run_status.tsv
        done < pwcoco_manifest.tsv
    fi
    echo "pwcoco runs: \$(grep -c 'ok' pwcoco_run_status.tsv || true) ok" | tee -a \$log

    # ---- 5. Collect results (+ v2.2 decision metrics via bin/coloc_decision.R) ----
    export NF_bin_dir='${projectDir}/bin'
    export NF_coloc_rule='${params.coloc_rule ?: "H4 > 0.8"}'
    export NF_coloc_h4_threshold='${params.coloc_h4_threshold ?: 0.8}'
    export NF_coloc_cond_threshold='${params.coloc_cond_threshold ?: 0.7}'
    export NF_coloc_power_min='${params.coloc_power_min ?: 0.5}'
    export NF_coloc_h3_strong='${params.coloc_h3_strong ?: 0.8}'
    Rscript "${projectDir}/bin/pwcoco_colocalization.R" parse "${gene_conversion_table}" "${gene_map}" 2>>\$log \
        || bail "parse step failed (see log)" parse_failed
    echo "status=done" >> \$prov
    echo "PWCoCo done." >> \$log
    """
}
