process COLOC_DIAGNOSTICS {
    container params.r_container
    publishDir "${params.outdir}/${params.tissue}/coloc_diagnostics", mode: 'copy'

    // v2.2: the R code lives in bin/coloc_diagnostics.R (same mechanism as the SuSiE
    // worker) so it can be run and tested outside Nextflow. Inputs arrive as NF_* vars.
    //
    // What it produces
    //   coloc_method_comparison_table.tsv  one row per classic-coloc gene: H4, H4/(H3+H4),
    //                                      H3+H4, class, pass_rule for classic, for the best
    //                                      real coloc.susie pair (+ susie_bf / fallback H4) and,
    //                                      when PWCoCo ran, for its best signal pair
    //   method_comparison_susie.png        classic H4 vs coloc.susie H4
    //   coloc_h4_vs_conditional.png        classic H4 vs H4/(H3+H4), coloured by H3+H4 (power)
    //   coloc_class_summary.png            tiered class counts, classic vs coloc.susie
    //   locus_<ENSG>_<gene>.png            LocusZoom-style stacked locus plots (exposure over
    //                                      outcome), points coloured by r2 with the index SNP
    //                                      (LocusZoom bins 0.8/0.6/0.4/0.2), index SNP = purple
    //                                      diamond, lead MR instrument = black triangle. Only
    //                                      for the genes selected by params.locus_plot_set
    //                                      ("both" = passed the classic rule AND a real
    //                                      coloc.susie pair passed it).
    //   scatter_<ENSG>_<gene>.png          cross-trait -log10(p) scatter, same colouring
    //   locus_plot_genes.tsv               which genes were plotted, why, and their index SNP

    input:
    path window_files                     // window_<ENSG>.rds from SuSiE (full cis-window, both traits)
    path ld_files                         // ld_<ENSG>.rds from SuSiE (signed R, rsid-named) or the NO_FILE_LD placeholder
    path coloc_summary                    // coloc_results_summary.rds (classic; H0..H4 + v2.2 metrics)
    path susie_summary                    // susie_coloc_results_summary.rds (PP.H4.abf, LD_type + v2.2 metrics)
    path pwcoco_summary                   // pwcoco_results_summary.rds (H0..H4 per signal pair + v2.2 metrics) or NO_FILE_PWCOCO
    path gene_conversion_table
    val tissue
    val outdir

    output:
    path "locus_*.png", optional: true, emit: locus_plots
    path "scatter_*.png", optional: true, emit: scatter_plots
    path "method_comparison_*.png", optional: true, emit: comparison_plots
    path "coloc_h4_vs_conditional.png", optional: true, emit: conditional_plot
    path "coloc_class_summary.png", optional: true, emit: class_plot
    path "coloc_method_comparison_table.tsv", optional: true, emit: comparison_table
    path "locus_plot_genes.tsv", optional: true, emit: locus_plot_genes
    path "coloc_diagnostics_log.txt", emit: log

    script:
    """
    export NF_bin_dir='${projectDir}/bin'
    export NF_coloc_summary='${coloc_summary}'
    export NF_susie_summary='${susie_summary}'
    export NF_pwcoco_summary='${pwcoco_summary}'
    export NF_gene_conversion_table='${gene_conversion_table}'
    export NF_tissue='${tissue}'
    export NF_coloc_rule='${params.coloc_rule}'
    export NF_coloc_h4_threshold='${params.coloc_h4_threshold}'
    export NF_susie_h4_threshold='${params.susie_h4_threshold}'
    export NF_coloc_cond_threshold='${params.coloc_cond_threshold}'
    export NF_coloc_power_min='${params.coloc_power_min}'
    export NF_coloc_h3_strong='${params.coloc_h3_strong}'
    export NF_locus_plot_set='${params.locus_plot_set}'
    Rscript "${projectDir}/bin/coloc_diagnostics.R"
    """
}
