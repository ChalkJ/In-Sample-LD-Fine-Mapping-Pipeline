#!/usr/bin/env Rscript
# =============================================================================
# Sumstats QC: optional duplicate-position + HetPVa filter, then combine 22
# per-chromosome sumstats files into the single combined daner-format .gz
# that run_ld_pipeline.sh hard-requires (SUMSTATS_NAME below must stay in
# sync with that script's case statement -- confirmed from its actual code,
# not assumed).
#
# Both filters are OPTIONAL and tunable -- not everyone wants them, and the
# threshold is a judgment call, not a fixed constant. Defaults match the
# rule confirmed for this project (2026-08-05):
#   - Duplicate CHR:BP within a chromosome -> keep the row with the lowest
#     P-value, drop the rest. Disable with --dedup=off.
#   - HetPVa < 0.05 -> drop the row. Change the cutoff with
#     --hetpva-cutoff=<value>, or disable entirely with --hetpva-cutoff=off.
#   - HetPVa missing/non-numeric (heterogeneity test not computed for that
#     SNP, e.g. only tested in one cohort) is always KEPT, never dropped --
#     a SNP can't fail a test that was never run.
#
# Usage:
#   Rscript qc_filter_sumstats.R <sczvscon|bipvscon>
#   Rscript qc_filter_sumstats.R <sczvscon|bipvscon> --dedup=off
#   Rscript qc_filter_sumstats.R <sczvscon|bipvscon> --hetpva-cutoff=0.01
#   Rscript qc_filter_sumstats.R <sczvscon|bipvscon> --dedup=off --hetpva-cutoff=off
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args  <- commandArgs(trailingOnly = TRUE)
PHENO <- args[1]
if (is.na(PHENO) || !(PHENO %in% c("sczvscon", "bipvscon"))) {
  stop("usage: Rscript qc_filter_sumstats.R <sczvscon|bipvscon> [--dedup=on|off] [--hetpva-cutoff=<value>|off]")
}

opt_args <- args[-1]
get_opt  <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), opt_args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[length(hit)])
}

DEDUP_ON     <- tolower(get_opt("dedup", "on")) != "off"
HETPVA_RAW   <- get_opt("hetpva-cutoff", "0.05")
HETPVA_ON    <- tolower(HETPVA_RAW) != "off"
HETPVA_CUTOFF <- if (HETPVA_ON) as.numeric(HETPVA_RAW) else NA_real_
if (HETPVA_ON && is.na(HETPVA_CUTOFF)) {
  stop("--hetpva-cutoff must be a number or 'off', got: ", HETPVA_RAW)
}

# Per-chromosome file prefix and combined-output filename are DIFFERENT
# strings on this project (confirmed in run_ld_pipeline.sh: PREFIX=
# sc_vs_allcontrols but SUMSTATS_NAME uses "scz_vs_allcontrols") -- kept as
# two explicit case statements rather than derived from one, to avoid
# silently drifting from what run_ld_pipeline.sh actually expects.
chr_prefix   <- switch(PHENO, sczvscon = "sc_vs_allcontrols",  bipvscon = "bp_vs_allcontrols")
out_name     <- switch(PHENO, sczvscon = "daner_scz_vs_allcontrols_1025_noduppos_hetpva.gz",
                               bipvscon = "daner_bip_vs_allcontrols_1025_noduppos_hetpva.gz")

FINEMAP_ROOT <- Sys.getenv("FINEMAP_ROOT")
if (FINEMAP_ROOT == "") {
  stop("FINEMAP_ROOT environment variable must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping")
}

BASE    <- file.path(FINEMAP_ROOT, PHENO)
OUT_GZ  <- file.path(BASE, out_name)
REPORT  <- file.path(BASE, "qc_sumstats_report.txt")

# --- Skip-if-done, matching the skip-if-done idiom used throughout this
# pipeline (01_extract_merge_clump.sh, run_ld_pipeline.sh) so resubmitting
# after a partial failure is always safe and cheap. Delete OUT_GZ to force
# a rerun (e.g. after changing --dedup/--hetpva-cutoff).
if (file.exists(OUT_GZ) && file.size(OUT_GZ) > 0) {
  cat(OUT_GZ, "already exists, skipping. Delete it to force a rerun with different settings.\n")
  quit(save = "no", status = 0)
}

report_lines <- character(0)
add_report   <- function(...) {
  line <- paste0(...)
  cat(line, "\n")
  report_lines <<- c(report_lines, line)
}

add_report("QC filter report for ", PHENO, " (", chr_prefix, ")")
add_report("Settings: dedup=", if (DEDUP_ON) "on (keep lowest P per CHR:BP)" else "off",
           "; hetpva_filter=", if (HETPVA_ON) paste0("on (drop HetPVa < ", HETPVA_CUTOFF, ")") else "off")
add_report("NA/non-numeric HetPVa is always kept, regardless of hetpva_filter setting.")
add_report("")

chr_tables   <- vector("list", 22)
header_cols  <- NULL
totals <- list(rows_in = 0L, dup_removed = 0L, hetpva_removed = 0L, hetpva_na_kept = 0L, rows_out = 0L)

for (CHR in 1:22) {
  f <- file.path(BASE, paste0(chr_prefix, "_chr", CHR, ".txt"))
  if (!file.exists(f)) {
    add_report("MISSING: ", f, " -- skipping chr", CHR)
    next
  }

  dt <- fread(f, sep = "\t", header = TRUE, colClasses = list(character = "HetPVa"))
  if (is.null(header_cols)) header_cols <- names(dt)

  n_in <- nrow(dt)

  # fread silently drops malformed rows (e.g. a truncated line with fewer
  # fields than the header, treated as a "footer") with only an R warning --
  # easy to miss in a SLURM .err log. Cross-check against a raw line count
  # so any such loss is explicit in the QC report, not just a suppressed
  # warning. (Hit once for real: sczvscon chr5 had one line truncated after
  # its Direction field, silently dropped by fread -- a raw-data artifact
  # upstream of this script, not something to guess is fine without checking.)
  raw_line_count <- length(readLines(f)) - 1L
  if (raw_line_count != n_in) {
    add_report(sprintf("WARNING: chr%d -- fread read %d rows but the file has %d data lines; %d row(s) silently dropped by fread (likely malformed/truncated lines) -- inspect %s directly",
                        CHR, n_in, raw_line_count, raw_line_count - n_in, f))
  }

  hetpva_num <- suppressWarnings(as.numeric(dt$HetPVa))
  na_mask    <- is.na(hetpva_num)
  n_na_kept  <- sum(na_mask)

  if (HETPVA_ON) {
    het_fail      <- !na_mask & hetpva_num < HETPVA_CUTOFF
    n_het_removed <- sum(het_fail)
    dt            <- dt[!het_fail]
  } else {
    n_het_removed <- 0L
  }

  if (DEDUP_ON) {
    # Dedup CHR:BP keeping lowest P. Ties on P keep the first-encountered row
    # (data.table's order() is stable), which is an arbitrary but deterministic
    # tiebreak -- not specified by the user, flagged here rather than hidden.
    dt[, P_num := suppressWarnings(as.numeric(P))]
    setorder(dt, BP, P_num)
    n_before_dedup <- nrow(dt)
    dt <- dt[!duplicated(BP)]
    n_dup_removed <- n_before_dedup - nrow(dt)
    dt[, P_num := NULL]
  } else {
    n_dup_removed <- 0L
  }

  n_out <- nrow(dt)
  add_report(sprintf("chr%-2d  in=%7d  dup_removed=%5d  hetpva_removed=%5d  hetpva_na_kept=%6d  out=%7d",
                      CHR, n_in, n_dup_removed, n_het_removed, n_na_kept, n_out))

  totals$rows_in        <- totals$rows_in + n_in
  totals$dup_removed     <- totals$dup_removed + n_dup_removed
  totals$hetpva_removed  <- totals$hetpva_removed + n_het_removed
  totals$hetpva_na_kept  <- totals$hetpva_na_kept + n_na_kept
  totals$rows_out        <- totals$rows_out + n_out

  chr_tables[[CHR]] <- dt
}

combined <- rbindlist(chr_tables[!sapply(chr_tables, is.null)])
setcolorder(combined, header_cols)

add_report("")
add_report(sprintf("TOTAL  in=%d  dup_removed=%d  hetpva_removed=%d  hetpva_na_kept=%d  out=%d",
                    totals$rows_in, totals$dup_removed, totals$hetpva_removed,
                    totals$hetpva_na_kept, totals$rows_out))

tmp_txt <- tempfile(fileext = ".txt")
fwrite(combined, tmp_txt, sep = "\t")

ok <- system2("gzip", c("-c", shQuote(tmp_txt)), stdout = OUT_GZ)
if (ok != 0 || !file.exists(OUT_GZ) || file.size(OUT_GZ) == 0) {
  stop("gzip failed writing ", OUT_GZ)
}
file.remove(tmp_txt)

writeLines(report_lines, REPORT)

add_report("")
add_report("Wrote ", OUT_GZ)
add_report("Wrote ", REPORT)
