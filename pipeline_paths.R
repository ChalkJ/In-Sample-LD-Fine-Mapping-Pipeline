#!/usr/bin/env Rscript
# =============================================================================
# Shared PHENO -> path mapping, sourced by run_finemapping.R and
# generate_reports.R. Single source of truth for this pipeline's on-disk
# layout -- avoids the path-drift bugs already hit twice on this project
# (ld/bipvscon vs ld/bpvscon, hardcoded constants out of sync across
# copy-pasted scripts).
#
# Requires the FINEMAP_ROOT environment variable (e.g.
# `export FINEMAP_ROOT=/path/to/your/finemapping`) -- see README.md
# "Configuration".
# =============================================================================

get_pipeline_paths <- function(pheno) {
  finemap_root <- Sys.getenv("FINEMAP_ROOT")
  if (finemap_root == "") {
    stop("FINEMAP_ROOT environment variable must be set -- e.g. export FINEMAP_ROOT=/path/to/your/finemapping")
  }

  # Combined QC'd sumstats filename is a fixed, phenotype-agnostic
  # convention -- qc_filter_sumstats.R writes this same name, so no
  # per-phenotype mapping is needed (or kept in sync) here.
  sumstats_name <- paste0("daner_", pheno, "_qc.gz")

  datadir <- file.path(finemap_root, pheno)

  list(
    pheno               = pheno,
    datadir             = datadir,
    lddir               = file.path(datadir, "output", "ld"),
    susiedir            = file.path(datadir, "susie"),
    fmdir               = file.path(datadir, "finemap"),
    sumstats            = file.path(datadir, sumstats_name),
    finemapping_summary = file.path(datadir, paste0("finemapping_", pheno, "_summary.tsv")),
    resultsdir          = file.path(datadir, "results"),
    plotsdir            = file.path(datadir, "results", "plots")
  )
}
