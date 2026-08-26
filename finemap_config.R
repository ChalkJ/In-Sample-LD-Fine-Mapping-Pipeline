#!/usr/bin/env Rscript
# =============================================================================
# Fine-mapping config, sourced (not executed) by run_finemapping.R -- see
# --finemap-config=<path> on 05_run_finemapping.sh / submit_full_pipeline.sh
# to point at a different copy of this file instead of editing this one.
#
# Edit the values below; the field names/structure must stay as-is (
# run_finemapping.R reads CONFIG$method, CONFIG$susie$lambda, etc.).
# =============================================================================

CONFIG <- list(

  # "susie"   -- run SuSiE only
  # "finemap" -- run FINEMAP only
  # "both"    -- run both and cross-check top-SNP agreement (default,
  #              matches every existing published result on this project)
  # "none"    -- skip fine-mapping entirely (e.g. you only want this
  #              pipeline's QC/clump/LD stages); run_finemapping.R exits
  #              immediately without error
  method = "both",

  # SuSiE's L (max number of credible sets) and FINEMAP's --n-causal-snps
  # (k) at the initial pass. Loci that saturate this (SuSiE non-converged
  # with L credible sets, or FINEMAP's Pr(k = max_causal_variants) exceeds
  # adaptive_rerun$k_posterior_threshold) get automatically re-run at
  # adaptive_rerun$max_causal_variants_rerun -- see below.
  max_causal_variants = 5,

  # Skip a locus entirely if its output already exists (SuSiE .rds / FINEMAP
  # .snp file, matching the method(s) selected above). Lets you resubmit
  # this stage after a partial failure without redoing finished loci.
  skip_existing = TRUE,

  susie = list(
    estimate_residual_variance = TRUE,
    lambda                     = 0.1
  ),

  finemap = list(
    # REQUIRED (non-blank, executable) if method is "finemap" or "both".
    # Snellius has no FINEMAP module -- place the binary somewhere in your
    # own scratch/home space first and set the full path here, e.g.:
    #   binary = "/gpfs/home3/<you>/finemapping/bin/finemap"
    # run_finemapping.R checks this up front and fails loudly before doing
    # any SuSiE work if it's needed but missing/not executable.
    binary = ""
  ),

  # Loci that saturate the initial max_causal_variants pass are
  # automatically re-run once at a higher L/k, then (if still not
  # converged) again with more SuSiE iterations. Only applies to whichever
  # method(s) `method` above actually runs.
  adaptive_rerun = list(
    enabled                   = TRUE,
    k_posterior_threshold     = 0.5,  # FINEMAP: Pr(k = max_causal_variants) above this = saturated
    max_causal_variants_rerun = 10,   # L/k used for the rerun pass
    rerun_max_iter            = 500   # SuSiE max_iter if still unconverged after the rerun pass
  )
)
