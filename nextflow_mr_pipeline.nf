#!/usr/bin/env nextflow

// =============================================================================
//  Cis-MR + colocalization pipeline (GTEx eQTL exposure -> GWAS outcome)
// -----------------------------------------------------------------------------
//  REQUIRED ENVIRONMENT VARIABLES  (secrets are never hardcoded in this file)
//
//    export OPENGWAS_JWT=...   OpenGWAS token. Read INSIDE the R process that
//                              clumps (CLUMP_EXPOSURE -> TwoSampleMR::clump_data
//                              -> ieugwasr::ld_clump, which calls the OpenGWAS
//                              API and authenticates with this JWT).
//    export LDLINK_TOKEN=...   LDlink token for the proxy search. Read on the
//                              HOST here (params.lddlink_token below) and passed
//                              into PROXY_BATCH_LDLINK as a value, so it does
//                              not need to be present inside the container.
//
//  CONTAINERS DO NOT INHERIT THE SHELL ENVIRONMENT. Because OPENGWAS_JWT is read
//  inside the container, it must be forwarded explicitly, e.g. in nextflow.config
//      docker.envWhitelist      = 'OPENGWAS_JWT'
//      singularity.envWhitelist = 'OPENGWAS_JWT'
//  or per process:  containerOptions '-e OPENGWAS_JWT'
//  See nextflow.config.example next to this file.
//
//  DEPENDENCY NOTE (changed): the proxy branch no longer starts from the raw
//  instruments. The chain is now
//      EXPOSURE_PREPROCESSING -> MAP_RSID_GAUSS -> CLUMP_EXPOSURE -> FIND_MISSING_SNPS
//  so the GAUSS reference index (--gauss_reference_index_file) and, when
//  --do_liftover is true, the chain file (--chain_file) are REQUIRED for EVERY
//  run, including runs with --run_susie false. Clumping is no longer optional.
// =============================================================================

nextflow.enable.dsl = 2

params.tissue = "Whole_Blood"
params.exposure_file = null
params.outcome_file = null
params.outdir = "results"
params.sample_size = 670
params.fstat_threshold = 10
params.p_threshold = 0.05
// Secrets come from the environment - never hardcode them here.
// Before running: export LDLINK_TOKEN=...   and   export OPENGWAS_JWT=...
// (OPENGWAS_JWT is read inside the clumping container - see the header note.)
params.lddlink_token = System.getenv('LDLINK_TOKEN') ?: ''
params.ld_threshold = 0.8
params.genome_build = "grch38"   // LDlink genome build: grch37 | grch38 | grch38_high_coverage
params.population = "ALL"        // LDlink reference population (e.g. EUR, ALL)
params.coloc_h4_threshold = 0.8
params.window_size = 500000
// Plain Groovy list, not a range: OUTCOME_SPLIT receives it as one val and renders
// it into the script as "[1, 2, ...]" (it also parses the "1..22" range form).
params.chromosomes = (1..22).toList()
params.gene_conversion_table = "../GeneCode_ConversionTable.rds"

// Image used by EVERY R process. Must contain TwoSampleMR, coloc, susieR, arrow,
// data.table, LDlinkR, rtracklayer and gauss (see environment.yml).
params.r_container = 'rocker/tidyverse:4.4.1'

params.run_susie = true
params.susie_h4_threshold = 0.8

// ---- Colocalization decision rule + metrics (v2.2; shared by classic coloc and SuSiE) ----
// Every coloc table now carries H3_plus_H4 (power mass), H4_cond = H4/(H3+H4) (the
// "colocalization probability" of Wu 2024 Nat Commun / Liu 2026 Commun Biol),
// H4_H3_odds, log2_H4_H3, dominant_hyp, a pass_rule flag and a tiered coloc_class
// (see bin/coloc_decision.R and COLOC_DECISION_METRICS.md).
// coloc_rule uses coloc::sensitivity() syntax on H0..H4:
//   "H4 > 0.8"                          classic (default; comparable with the literature)
//   "H4/(H3+H4) > 0.7 & H3+H4 > 0.5"    conditional probability + power guard (Guo 2015 / Wu 2024)
//   "H4 > 0.9 & H4/H3 > 3"              example from the coloc documentation
// The classic H4 thresholds above still define the 'colocalized_strong' tier.
params.coloc_rule           = "H4 > 0.8"
params.coloc_cond_threshold = 0.7     // H4/(H3+H4) needed for the 'colocalized_conditional' tier
params.coloc_power_min      = 0.5     // minimum H3+H4 before the ratio is trusted (Guo 2015 used 0.8)
params.coloc_h3_strong      = 0.8     // H3 above this -> 'distinct_signals' (LD-confounded MR)
params.coloc_p12_grid       = "1e-5,5e-6,1e-6"   // prior-sensitivity columns H4_p12_<tag>
params.coloc_fdr_target     = 0.05    // Reales 2026: smallest H4 alpha with expected FDR below this
params.coloc_min_snps       = 50      // coloc.abf needs a dense map; was 3
params.ld_check_r2          = 0.8     // Zheng 2020 LD check (SuSiE step): max r2(instrument, top-30 GWAS SNPs)
// Window centre for BOTH coloc modules: "lead_instrument" (min pval.exposure among the
// gene's clumped instruments) or "first_instrument" (v2.1 behaviour: first harmonized row).
params.coloc_window_center  = "lead_instrument"
// Outcome trait type for coloc: "quant" (sdY from --outcome_sdY or MAF+N) or "cc"
// (case/control GWAS: logOR betas; needs the case fraction in --outcome_case_prop).
params.outcome_type         = "quant"
params.outcome_case_prop    = null
// Which genes get LocusZoom-style locus plots in COLOC_DIAGNOSTICS:
//   "both"        passed the classic rule AND a real coloc.susie credible-set pair (default)
//   "classic"     passed the classic rule; "susie" passed a real coloc.susie pair;
//   "pwcoco"      a PWCoCo signal pair passed the rule (needs --pwcoco_bfile)
//   "all_methods" classic AND coloc.susie AND PWCoCo (PWCoCo ignored when it did not run)
//   "any"         any of the three; "all" every gene with a window
params.locus_plot_set       = "both"
// Shared coloc prior for a variant affecting BOTH traits. Passed explicitly to both
// coloc.abf (classic) and coloc.susie so the two methods use the same prior.
params.coloc_p12 = 1e-5
// When coloc.susie finds no credible-set pair (or runsusie fails), fall back to
// coloc.abf on the same window so genes are never silently dropped.
params.susie_abf_fallback = true
// Before the ABF fallback, try coloc.susie_bf: SuSiE credible sets of the trait that
// resolves x single-signal ABF of the other (rows tagged LD_type = "susie_bf").
params.susie_hybrid = true
// estimate_s_rss lambda above this flags a likely LD/summary-stat mismatch.
params.susie_lambda_warn = 0.2
// GTEx all-pairs association files (PARQUET) used by BOTH coloc.abf and SuSiE.
// Accepts a glob, e.g. "/data/GTEx/Whole_Blood.v11.cis_qtl_pairs.chr*.parquet"
params.assoc_files = null
// ---- Gene-ID bridging between GENCODE releases ----
// The eGene/exposure file and the all-pairs parquet files often come from
// different GTEx releases (hence different GENCODE annotations), so a gene may
// exist in one and not the other under the same ENSG. Supplying the "genes.gtf"
// GENCODE file that matches EACH release lets GENE_RELEASE_MAP bridge the ids by
// gene name and by coordinate overlap instead of dropping those genes.
params.gtf_exposure = null          // GTF matching the eGene/exposure release (optional)
params.gtf_assoc    = null          // GTF matching the all-pairs release (optional)
params.gene_map_min_overlap = 0.5   // minimum overlap fraction for coordinate bridging
// Exposure N inside the SuSiE/SharePro windows, which are filled with ALL-PAIRS
// statistics. params.sample_size is the eGene-release tissue N used for the MR
// F-statistics; set this when the all-pairs release has a different tissue N.
params.assoc_sample_size = null
// Coloc/SuSiE failures terminate by default so they surface. Set to 'ignore' to skip.
params.coloc_error_strategy = 'terminate'

// ---- LD clumping of the exposure instruments (CLUMP_EXPOSURE) ----
// Clumping ALWAYS runs now: with an all-pairs exposure file exposure_by_chr holds
// many correlated SNPs per gene, so MR needs the clumped, independent instruments.
params.clump_kb = 10000        // clumping window, kb (per gene cis-window)
params.clump_r2 = 0.001        // r2 threshold for independence
params.clump_p  = 1            // p-value threshold (1 = keep every instrument)
params.clump_pop = "EUR"       // single population for clump_data (LD reference)
params.clump_retries = 3       // attempts against the OpenGWAS API before giving up

// ---- Resource knobs for the two memory-hungry single-task processes ----
params.rsid_map_memory = '64 GB'   // MAP_RSID_GAUSS holds the whole GAUSS index
params.outcome_memory  = '32 GB'   // OUTCOME_SPLIT reads the GWAS once, in full

// ---- MR result combination (COMBINE_MR_RESULTS) ----
// Steiger is always computed; genes failing it are only DROPPED when this is true.
params.steiger_filter = false
// Primary estimate for genes with more than one instrument (nsnp == 1 -> Wald ratio).
params.mr_primary_multi = "Inverse variance weighted"

// ---- coloc trait variance (sdY) ----
// outcome_sdY: null -> coloc estimates sdY from the GWAS MAF + N; a number -> fixed.
// exposure_sdY: GTEx expression is inverse-normal transformed, so 1.
params.outcome_sdY  = null
params.exposure_sdY = 1

// ---- SuSiE LD via GAUSS (ancestry-weighted multi-ethnic) ----
params.gauss_reference_index_file    = "/path/to/33kg_index.gz"
params.gauss_reference_data_file     = "/path/to/33kg_geno.gz"
params.gauss_reference_pop_desc_file = "/path/to/33kg_pop_desc.txt"
params.gauss_af1_cutoff              = 0.001
// GAUSS 33k panel is GRCh37 while GTEx/GWAS here are GRCh38: lift coords for LD lookup.
params.do_liftover = true
params.chain_file  = "/path/to/hg38ToHg19.over.chain.gz"
// Fixed ancestry proportions from the GWAS cohort Ns (labels must match pop_desc):
params.susie_ancestry_props = ""  // ignored; ancestry is always inferred via afmix
// Which trait's ancestry weights build the LD used by SuSiE: "outcome" (one shared
// matrix) or "exposure" (a second, eQTL-ancestry matrix for the exposure).
params.susie_exposure_ld_weights = "outcome"

// ---- SharePro (method comparison) ----
params.sharepro_k = 10          // SharePro max causal signals (--K)
params.sharepro_dir = null      // existing SharePro_coloc checkout; null/"" -> clone
params.sharepro_commit = "main" // git ref checked out when cloning

// ---- PWCoCo (pair-wise conditional analysis + coloc.abf per signal pair) ----
// Needs an INDIVIDUAL-LEVEL reference panel in PLINK 1 format (bed/bim/fam),
// same population as the GWAS, ideally >= 4,000 people (GCTA-COJO guidance).
// "{chr}" in the prefix is replaced by the chromosome number, e.g.
//   --pwcoco_bfile /project/1000G/plink/EUR/Chr{chr}.EUR
// Leave null to skip PWCoCo (the comparison tables then have no pwcoco_* columns).
// v2.2: PWCoCo rows carry the same decision metrics / rule / class as the other
// methods (bin/coloc_decision.R) and feed COLOC_DIAGNOSTICS, SharePro and the report.
params.pwcoco_bfile = null
params.pwcoco_bin = "pwcoco"           // compiled binary (github.com/jwr-git/pwcoco, cmake + make)
params.plink_bin = "plink"             // PLINK 1.9 or 2; falls back to 'module load plink'
params.pwcoco_match = "position"       // "position": panel on GRCh38 (matches GTEx/GWAS here)
                                       // "rsid": any build, .bim ids must be rsIDs
params.pwcoco_p_cutoff1 = 5e-8         // stepwise-selection threshold, exposure (eQTL)
params.pwcoco_p_cutoff2 = 5e-8         // stepwise-selection threshold, outcome (GWAS)
params.pwcoco_maf = 0.01               // panel MAF filter (PWCoCo's own default is 0.1; use 0.05 with a ~500-person panel)
params.pwcoco_freq_threshold = 0.2     // drop SNPs whose freq differs from the panel by more
params.pwcoco_collinear = 0.9          // COJO collinearity cut-off
params.pwcoco_init_h4 = 80             // stop after the marginal coloc if H4 >= this (%); 100 = always condition
params.pwcoco_min_snps = 10            // minimum window SNPs found in the panel to run a gene
params.pwcoco_verbose = false          // PWCoCo --verbose (writes .badfreq / conditioned .cojo files)
params.pwcoco_extra_args = ""          // anything else to pass to pwcoco verbatim
params.pwcoco_error_strategy = "ignore"

include { EXPOSURE_PREPROCESSING } from './modules/exposure_preprocessing.nf'
include { MAP_RSID_GAUSS } from './modules/rsid_map_gauss.nf'
include { CLUMP_EXPOSURE } from './modules/clump_exposure.nf'
include { GENE_RELEASE_MAP } from './modules/gene_release_map.nf'
include { FIND_MISSING_SNPS; PROXY_BATCH_LDLINK; PROCESS_PROXY_RESULTS; MERGE_PROXIED_DATA } from './modules/proxy_batch.nf'
include { OUTCOME_SPLIT } from './modules/outcome_preprocessing.nf'
include { MR_ANALYSIS; COMBINE_MR_RESULTS } from './modules/mr_analysis.nf'
include { COLOCALIZATION } from './modules/colocalization.nf'
include { SUSIE_COLOCALIZATION } from './modules/susie_colocalization.nf'
include { COLOC_DIAGNOSTICS } from './modules/coloc_diagnostics.nf'
include { SHAREPRO_COLOCALIZATION } from './modules/sharepro_colocalization.nf'
include { PWCOCO_COLOCALIZATION } from './modules/pwcoco_colocalization.nf'
include { RESULTS_PLOTTING } from './modules/results_plotting.nf'
include { GENERATE_REPORT } from './modules/report.nf'

workflow {
    main:
        if (!params.exposure_file || !params.outcome_file) {
            error "ERROR: --exposure_file and --outcome_file are required"
        }

        log.info "Starting MR Pipeline for ${params.tissue}..."

        // A `val` process input bound to a Groovy null becomes an UNBOUND value
        // channel (Channel.value(null) is empty), and a process whose input never
        // gets a value simply never runs. Nullable params are therefore passed as
        // "" and each module maps "" to its documented null behaviour:
        //   outcome_sdY  -> as.numeric("") is NA -> coloc estimates sdY from MAF+N
        //   sharepro_dir -> empty -> clone SharePro_coloc at --sharepro_commit
        outcome_sdY_val  = params.outcome_sdY  != null ? params.outcome_sdY  : ''
        sharepro_dir_val = params.sharepro_dir != null ? params.sharepro_dir : ''

        // ---- Exposure instruments -------------------------------------------
        // EXPOSURE_PREPROCESSING is FILTERING ONLY now: it stops at the
        // per-chromosome, indel/palindrome-free instrument set.
        exposure_preprocessed = EXPOSURE_PREPROCESSING(
            params.exposure_file, params.tissue, params.sample_size,
            params.fstat_threshold, params.outdir
        )

        // rsIDs for every instrument, built once from the GAUSS index (the
        // OpenGWAS clumping API matches on rsID; GTEx carries none).
        log.info "Building the exposure rsID map from the GAUSS reference index..."
        rsid_map = MAP_RSID_GAUSS(
            exposure_preprocessed.exposure_by_chr.collect(),
            params.gauss_reference_index_file,
            params.do_liftover,
            params.chain_file
        )

        // Per-chromosome LD clumping. The rsID map is turned into a VALUE channel
        // (a single-task output is already a value channel, so every chromosome task receives it).
        exposure_by_chr_tuples = exposure_preprocessed.exposure_by_chr
            .flatten()
            .map { f -> def m = (f.getName() =~ /exposure_by_chr_(\d+)\.rds/); tuple(m[0][1] as Integer, f) }

        log.info "LD-clumping the exposure instruments, one chromosome per task..."
        clumped = CLUMP_EXPOSURE(
            exposure_by_chr_tuples,
            rsid_map.rsid_map,
            params.clump_kb,
            params.clump_r2,
            params.clump_p,
            params.clump_pop,
            params.clump_retries
        )

        // Parquet association files, shared by coloc.abf, SuSiE and the gene-ID map
        assoc_ch = Channel.fromPath(params.assoc_files, checkIfExists: true).collect()

        // Optional GENCODE GTFs. Nextflow can only stage files that exist, so when a
        // GTF is not supplied we stage an empty placeholder the process recognises by
        // name (two distinct names so the staged files never collide).
        def no_gtf_exp = file("NO_FILE")
        def no_gtf_ass = file("NO_FILE2")
        if (!params.gtf_exposure && !no_gtf_exp.exists()) no_gtf_exp.text = ''
        if (!params.gtf_assoc    && !no_gtf_ass.exists()) no_gtf_ass.text = ''
        gtf_exp_ch   = params.gtf_exposure ? file(params.gtf_exposure, checkIfExists: true) : no_gtf_exp
        gtf_assoc_ch = params.gtf_assoc    ? file(params.gtf_assoc,    checkIfExists: true) : no_gtf_ass

        log.info "Mapping exposure gene IDs onto the all-pairs release..."
        gene_map = GENE_RELEASE_MAP(
            exposure_preprocessed.exposure_raw,
            assoc_ch,
            gtf_exp_ch,
            gtf_assoc_ch,
            file(params.gene_conversion_table),
            params.gene_map_min_overlap
        )

        // Exposure N for windows built from the ALL-PAIRS release (SuSiE/SharePro).
        assoc_n = params.assoc_sample_size ?: params.sample_size
        if (!params.assoc_sample_size) {
            log.warn "assoc_sample_size not set; using sample_size (${params.sample_size}) as the all-pairs exposure N — set --assoc_sample_size to the tissue N of the all-pairs release"
        }

        // ---- Outcome SNPs missing for the (clumped) instruments -> LD proxies ----
        missing_snps = FIND_MISSING_SNPS(
            clumped.clumped.collect(),
            params.outcome_file, params.tissue
        )

        proxy_results = PROXY_BATCH_LDLINK(
            missing_snps.missing_snps, params.lddlink_token, params.genome_build, params.population
        )

        proxies_processed = PROCESS_PROXY_RESULTS(
            proxy_results.proxy_results, params.outcome_file, params.ld_threshold
        )

        proxied_outcome_data = MERGE_PROXIED_DATA(
            proxies_processed.proxies_processed,
            missing_snps.missing_snps,
            params.outcome_file,
            params.ld_threshold
        )

        // ---- Outcome: read the GWAS once, split it per chromosome ----
        outcome_split = OUTCOME_SPLIT(
            params.outcome_file,
            proxied_outcome_data.outcome_with_proxies,
            params.tissue,
            params.chromosomes,
            params.outdir
        )

        log.info "Running MR Analysis..."
        exposure_files = clumped.clumped
            .flatten()
            .map { f -> def m = (f.getName() =~ /exposure_clumped_(\d+)\.rds/); tuple(m[0][1] as Integer, f) }
        outcome_files = outcome_split.outcome_by_chr
            .flatten()
            .map { f -> def m = (f.getName() =~ /outcome_chr_(\d+)\.rds/); tuple(m[0][1] as Integer, f) }
        mr_input = exposure_files.join(outcome_files)

        mr_results = MR_ANALYSIS(
            mr_input, params.tissue, params.sample_size, params.p_threshold, params.outdir
        )

        log.info "Combining MR results, computing threshold from harmonized data, removing MHC..."
        combined_mr = COMBINE_MR_RESULTS(
            mr_results.all_results.collect(),
            mr_results.harmonized_data.collect(),
            mr_results.gene_count.collect(),
            params.p_threshold,
            params.steiger_filter,
            params.mr_primary_multi,
            params.tissue
        )

        log.info "Running colocalization (decision rule: ${params.coloc_rule}; strong tier H4 > ${params.coloc_h4_threshold})..."
        coloc_results = COLOCALIZATION(
            mr_results.harmonized_data.collect(),
            combined_mr.significant,
            assoc_ch,
            outcome_split.outcome_raw,
            file(params.gene_conversion_table),
            gene_map.map_rds,
            params.tissue,
            params.window_size,
            params.coloc_h4_threshold,
            params.coloc_p12,
            outcome_sdY_val,
            params.exposure_sdY,
            params.outdir
        )

        // Placeholders for the optional SuSiE outputs. They are declared BEFORE the
        // if-block so they stay in scope afterwards, and they carry three DISTINCT
        // names so two empty inputs of the same process never collide when staged.
        // The modules recognise them by the "NO_FILE" prefix / by globbing only
        // window_*.rds and ld_*.rds.
        def no_win   = file("NO_FILE_WINDOWS")
        def no_ld    = file("NO_FILE_LD")
        def no_susie = file("NO_FILE_SUSIE")
        def no_susie_gs = file("NO_FILE_SUSIE_GS")   // susie_gene_summary.tsv placeholder (report)
        def no_pwcoco = file("NO_FILE_PWCOCO")       // pwcoco_results_summary.rds placeholder
        if (!no_win.exists())   no_win.text   = ''
        if (!no_ld.exists())    no_ld.text    = ''
        if (!no_susie.exists()) no_susie.text = ''
        if (!no_susie_gs.exists()) no_susie_gs.text = ''
        if (!no_pwcoco.exists()) no_pwcoco.text = ''
        // The PWCoCo table reaches COLOC_DIAGNOSTICS, SharePro and the report; it is
        // the placeholder unless SuSiE ran (PWCoCo uses its windows) AND a panel was given.
        pwcoco_table = no_pwcoco

        if (params.run_susie) {
            log.info "Running SuSiE colocalization (GAUSS ancestry-weighted LD) on the FDR-significant gene set..."
            susie_results = SUSIE_COLOCALIZATION(
                mr_results.harmonized_data.collect(),
                combined_mr.significant,
                assoc_ch,
                outcome_split.outcome_raw,
                file(params.gene_conversion_table),
                gene_map.map_rds,
                params.tissue,
                assoc_n,
                params.window_size,
                params.susie_h4_threshold,
                params.coloc_p12,
                params.susie_abf_fallback,
                params.susie_lambda_warn,
                outcome_sdY_val,
                params.exposure_sdY,
                params.susie_exposure_ld_weights,
                params.gauss_reference_index_file,
                params.gauss_reference_data_file,
                params.gauss_reference_pop_desc_file,
                params.gauss_af1_cutoff,
                params.do_liftover,
                params.chain_file,
                params.susie_ancestry_props,
                params.outdir
            )

            // window_*.rds / ld_*.rds are OPTIONAL outputs: when SuSiE colocalises
            // nothing those channels stay empty and the two consumers below would
            // wait for input for ever. ifEmpty() feeds them the placeholder instead
            // so the processes still run (and say so in their logs).
            windows_ch = susie_results.windows.collect().ifEmpty(no_win)
            ld_ch      = susie_results.ld_matrices.collect().ifEmpty(no_ld)

            // ---- PWCoCo on the same SuSiE windows (only with a PLINK panel) ----
            // Without --pwcoco_bfile the empty placeholder goes to the comparison
            // steps, which then simply omit the pwcoco_* columns.
            if (params.pwcoco_bfile) {
                log.info "Running PWCoCo (conditional analysis + coloc.abf per signal pair)..."
                pwcoco_results = PWCOCO_COLOCALIZATION(
                    windows_ch,
                    file(params.gene_conversion_table),
                    gene_map.map_rds,
                    params.tissue,
                    assoc_n,
                    params.coloc_p12,
                    params.pwcoco_bfile,
                    params.outdir
                )
                pwcoco_table = pwcoco_results.pwcoco_table
            } else {
                log.warn "pwcoco_bfile not set: PWCoCo skipped (no pwcoco_* columns in the comparison tables)"
            }

            log.info "Generating coloc diagnostic + method-comparison plots (LocusZoom-style locus plots for genes passing '${params.locus_plot_set}')..."
            COLOC_DIAGNOSTICS(
                windows_ch,
                ld_ch,                       // ld_<ENSG>.rds: r2 colouring of the locus plots
                coloc_results.coloc_table,
                susie_results.susie_table,
                pwcoco_table,                // pwcoco_results_summary.rds or the NO_FILE_PWCOCO placeholder
                file(params.gene_conversion_table),
                params.tissue,
                params.outdir
            )

            log.info "Running SharePro colocalization for comparison..."
            SHAREPRO_COLOCALIZATION(
                windows_ch,
                ld_ch,
                coloc_results.coloc_table,
                susie_results.susie_table,
                pwcoco_table,
                file(params.gene_conversion_table),
                gene_map.map_rds,
                params.tissue,
                assoc_n,
                params.sharepro_k,
                sharepro_dir_val,
                params.sharepro_commit,
                params.outdir
            )
        }

        log.info "Generating plots..."
        // RESULTS_PLOTTING always runs; without SuSiE it gets the NO_FILE_SUSIE
        // placeholder and simply omits the SuSiE panels.
        susie_tbl_ch = params.run_susie ? susie_results.susie_table : no_susie
        final_plots = RESULTS_PLOTTING(
            combined_mr.combined.collect(),
            coloc_results.coloc_table,
            susie_tbl_ch,
            file(params.gene_conversion_table),
            combined_mr.threshold_info,
            params.tissue,
            params.outdir
        )

        log.info "Generating summary report..."
        // Colocalization lines in the report: classic table + FDR table always exist;
        // the SuSiE gene summary is replaced by a placeholder when SuSiE did not run.
        susie_gs_ch = params.run_susie ? susie_results.gene_summary : no_susie_gs
        GENERATE_REPORT(
            clumped.clumped.collect(),
            missing_snps.counts,
            proxied_outcome_data.counts,
            combined_mr.threshold_info,
            gene_map.map_tsv,
            coloc_results.coloc_table,
            coloc_results.fdr_table,
            susie_gs_ch,
            pwcoco_table,
            params.tissue,
            params.fstat_threshold,
            params.outdir
        )

        log.info "Pipeline completed!"
}
