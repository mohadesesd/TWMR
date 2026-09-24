process SUSIE_COLOCALIZATION {
    container params.r_container   // must also contain gauss + susieR + coloc
    publishDir "${params.outdir}/${params.tissue}/colocalization/susie", mode: 'copy'

    maxRetries 1
    errorStrategy { params.coloc_error_strategy ?: 'terminate' }
    // afmix + computeLD load the 44GB GAUSS reference panel; give generous memory
    // to avoid segfaults from OOM. Adjust via the params if your cluster differs.
    memory { params.susie_memory ?: '96 GB' }
    cpus   { params.susie_cpus   ?: 4 }
    time   { params.susie_time   ?: '24h' }

    input:
    path harmonized_files
    path proxy_genes              // the FDR-significant gene set (same input as classic coloc)
    path assoc_files              // GTEx all-pairs association file(s), parquet
    path outcome_raw
    path gene_conversion_table
    path gene_map                 // gene_release_map.rds from GENE_RELEASE_MAP
    val tissue
    val sample_size
    val window_size
    val susie_h4_threshold
    val coloc_p12                 // shared coloc prior (also used by classic coloc.abf)
    val susie_abf_fallback        // run coloc.abf when coloc.susie yields no credible-set pair
    val susie_lambda_warn         // estimate_s_rss lambda above this flags an LD mismatch
    val outcome_sdY               // numeric -> fixed sdY for the GWAS; null/NA -> estimate from MAF+N
    val exposure_sdY              // GTEx expression is inverse-normal transformed -> 1
    val exposure_ld_weights       // "outcome" (one LD matrix) | "exposure" (own afmix weights for the eQTL)
    val ref_index_file
    val ref_data_file
    val ref_pop_desc_file
    val af1_cutoff
    val do_liftover
    val chain_file
    val ancestry_props
    val outdir

    output:
    path "susie_coloc_results_summary.rds", emit: susie_table
    path "susie_coloc_results_summary.tsv", emit: susie_tsv
    path "susie_passing_genes.rds", emit: susie_passing
    path "susie_gene_summary.tsv", emit: gene_summary     // one row per gene: best pair, metrics, class, LD check
    path "susie_gene_status.tsv", emit: gene_status
    path "ancestry_proportions_used.rds", emit: ancestry_props, optional: true
    path "ancestry_proportions_exposure.rds", emit: ancestry_props_exposure, optional: true
    path "susie_coloc_log.txt", emit: log
    path "sessionInfo_susie.txt", emit: sessioninfo
    path "window_*.rds", optional: true, emit: windows
    path "ld_*.rds", optional: true, emit: ld_matrices
    // Exposure-ancestry LD, only written when --susie_exposure_ld_weights = "exposure".
    // NOTE: these files also match the "ld_*.rds" glob above, so they appear in
    // ld_matrices as well; SharePro looks up ld_<ENSG>.rds by name and ignores them.
    path "ld_exp_*.rds", optional: true, emit: ld_exp
    path "gauss_input_*.txt", emit: gauss_inputs, optional: true

    script:
    """
    # The R worker lives in bin/ because a Groovy string constant is capped at
    # 65,535 bytes and this script is larger. Inputs are passed as NF_* env vars.
    export NF_susie_lambda_warn='${susie_lambda_warn}'
    export NF_susie_abf_fallback='${susie_abf_fallback}'
    export NF_exposure_sdY='${exposure_sdY}'
    export NF_outcome_sdY='${outcome_sdY}'
    export NF_exposure_ld_weights='${exposure_ld_weights}'
    export NF_ref_index_file='${ref_index_file}'
    export NF_ref_data_file='${ref_data_file}'
    export NF_ref_pop_desc_file='${ref_pop_desc_file}'
    export NF_af1_cutoff='${af1_cutoff}'
    export NF_gene_map='${gene_map}'
    export NF_outcome_raw='${outcome_raw}'
    export NF_do_liftover='${do_liftover}'
    export NF_chain_file='${chain_file}'
    export NF_gene_conversion_table='${gene_conversion_table}'
    export NF_proxy_genes='${proxy_genes}'
    export NF_sample_size='${sample_size}'
    export NF_susie_hybrid='${params.susie_hybrid == null ? true : params.susie_hybrid}'
    export NF_window_size='${window_size}'
    export NF_coloc_p12='${coloc_p12}'
    export NF_susie_h4_threshold='${susie_h4_threshold}'
    # ---- v2.2: shared decision metrics / rule (bin/coloc_decision.R), LD check,
    #      outcome trait type, window centre; read with nf_env_opt() in the R worker
    export NF_bin_dir='${projectDir}/bin'
    export NF_coloc_rule='${params.coloc_rule}'
    export NF_coloc_cond_threshold='${params.coloc_cond_threshold}'
    export NF_coloc_power_min='${params.coloc_power_min}'
    export NF_coloc_h3_strong='${params.coloc_h3_strong}'
    export NF_coloc_p12_grid='${params.coloc_p12_grid}'
    export NF_ld_check_r2='${params.ld_check_r2}'
    export NF_coloc_window_center='${params.coloc_window_center}'
    export NF_outcome_type='${params.outcome_type}'
    export NF_outcome_case_prop='${params.outcome_case_prop ?: ''}'
    Rscript "${projectDir}/bin/susie_colocalization.R"
    """
}
