#!/usr/bin/env Rscript
# =============================================================================
# Sumstats QC: optional duplicate-position + HetPVa filter, then write the
# combined daner-format .gz that run_ld_pipeline.sh and pipeline_paths.R
# expect (out_name below must stay in sync with those -- it's the same
# fixed, phenotype-agnostic convention in all three places).
#
# Takes a single whole-genome sumstats file (see README.md "Prerequisites")
# -- there is no per-chromosome file requirement.
#
# Both filters are OPTIONAL and tunable -- not everyone wants them, and the
# threshold is a judgment call, not a fixed constant. Defaults match the
# rule confirmed for this project (2026-08-05):
#   - Duplicate CHR:BP -> keep the row with the lowest P-value, drop the
#     rest. Disable with --dedup=off.
#   - HetPVa < 0.05 -> drop the row. Change the cutoff with
#     --hetpva-cutoff=<value>, or disable entirely with --hetpva-cutoff=off.
#   - HetPVa missing/non-numeric (heterogeneity test not computed for that
#     SNP, e.g. only tested in one cohort) is always KEPT, never dropped --
#     a SNP can't fail a test that was never run.
#
# Usage:
#   Rscript qc_filter_sumstats.R <phenotype> --sumstats-file=<filename>
#   Rscript qc_filter_sumstats.R <phenotype> --sumstats-file=<filename> --dedup=off
#   Rscript qc_filter_sumstats.R <phenotype> --sumstats-file=<filename> --hetpva-cutoff=0.01
#   Rscript qc_filter_sumstats.R <phenotype> --sumstats-file=<filename> --dedup=off --hetpva-cutoff=off
# =============================================================================

suppressPackageStartupMessages(library(data.table))

args  <- commandArgs(trailingOnly = TRUE)
PHENO <- args[1]
if (is.na(PHENO) || PHENO == "") {
  stop("usage: Rscript qc_filter_sumstats.R <phenotype> --sumstats-file=<filename> [--dedup=on|off] [--hetpva-cutoff=<value>|off]")
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

# Whole-genome sumstats filename is user-supplied (however your own upstream
# GWAS pipeline named it) -- not something this pipeline can derive on its
# own.
sumstats_file <- get_opt("sumstats-file", "")
if (sumstats_file == "") {
  stop("--sumstats-file=<filename> is required -- your whole-genome sumstats file")
}

# Combined-output filename is a fixed, phenotype-agnostic convention (this
# is pipeline-internal -- qc_filter_sumstats.R generates this file itself,
# unlike sumstats_file above which describes pre-existing user data) -- must
# stay in sync with pipeline_paths.R and run_ld_pipeline.sh, which both
# read this same name.
out_name <- paste0("daner_", PHENO, "_qc.gz")

FINEMAP_ROOT <- Sys.getenv("FINEMAP_ROOT")
if (FINEMAP_ROOT == "") {
  stop("FINEMAP_ROOT environment variable must be set -- e.g. export FINEMAP_ROOT=/path/to/your/finemapping")
}

BASE    <- file.path(FINEMAP_ROOT, PHENO)
IN_FILE <- file.path(BASE, sumstats_file)
OUT_GZ  <- file.path(BASE, out_name)
REPORT  <- file.path(BASE, "qc_sumstats_report.txt")

if (!file.exists(IN_FILE)) {
  stop("Sumstats file not found: ", IN_FILE)
}

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

add_report("QC filter report for ", PHENO, " (", sumstats_file, ")")
add_report("Settings: dedup=", if (DEDUP_ON) "on (keep lowest P per CHR:BP)" else "off",
           "; hetpva_filter=", if (HETPVA_ON) paste0("on (drop HetPVa < ", HETPVA_CUTOFF, ")") else "off")
add_report("NA/non-numeric HetPVa is always kept, regardless of hetpva_filter setting.")
add_report("")

cat("Loading", IN_FILE, "...\n")
# fill=TRUE is required here, not optional: with a single multi-million-row
# whole-genome file, fread's default behaviour on hitting a short/ragged
# line *mid-file* is not to drop just that one line -- it stops parsing
# entirely at that point and silently truncates everything after it.
# Confirmed for real during validation: one truncated line partway through
# a 7.18M-row real sumstats file caused fread to return only 2.66M rows
# (63% of the data silently lost) without fill=TRUE. This is a materially
# different (and far worse) failure mode than the old per-chromosome-file
# version of this script ever saw, where the same kind of truncated line
# only ever appeared at/near the end of a much smaller single-chromosome
# file and was merely dropped as a "footer". fill=TRUE pads missing
# trailing fields with NA instead of aborting, which is what we want.
dt <- fread(IN_FILE, sep = "\t", header = TRUE, colClasses = list(character = "HetPVa"), fill = TRUE)
header_cols <- names(dt)
n_in_total  <- nrow(dt)

# With fill=TRUE, a short/ragged line is no longer dropped, so raw line
# count vs. fread row count will usually match even when a row was
# malformed -- this check now mainly catches a genuine parse failure
# (fread erroring/dying) rather than the ragged-row case.
raw_line_count <- length(readLines(IN_FILE)) - 1L
if (raw_line_count != n_in_total) {
  add_report(sprintf("WARNING: fread read %d rows but the file has %d data lines; %d row(s) may have been dropped -- inspect %s directly",
                      n_in_total, raw_line_count, raw_line_count - n_in_total, IN_FILE))
}

# Catch ragged/short rows explicitly instead: fill=TRUE pads a short row's
# missing trailing columns with NA, so a real GWAS meta-analysis row that
# always populates its last column (here, Neff_half) but comes back NA is
# a strong signal of a truncated line, not a real missing value. Flagged
# as a heads-up, not a hard error -- some other sumstats format could
# legitimately leave a trailing column blank for a subset of rows.
last_col <- header_cols[length(header_cols)]
n_ragged <- sum(is.na(dt[[last_col]]))
if (n_ragged > 0) {
  add_report(sprintf("WARNING: %d row(s) have NA in the last column (%s) after fill=TRUE padding -- likely truncated/ragged lines in the raw file, not necessarily real missing data. Affected SNPs: %s",
                      n_ragged, last_col, paste(head(dt[is.na(get(last_col)), SNP], 20), collapse = ", ")))
}

# Per-chromosome row counts, captured before any filtering, purely for the
# report below.
n_in_by_chr <- dt[, .N, by = CHR]

hetpva_num <- suppressWarnings(as.numeric(dt$HetPVa))
na_mask    <- is.na(hetpva_num)
n_na_by_chr <- dt[na_mask, .N, by = CHR]

if (HETPVA_ON) {
  het_fail <- !na_mask & hetpva_num < HETPVA_CUTOFF
} else {
  het_fail <- rep(FALSE, nrow(dt))
}
n_het_removed_by_chr <- dt[het_fail, .N, by = CHR]
dt <- dt[!het_fail]

if (DEDUP_ON) {
  # Dedup CHR:BP keeping lowest P. Ties on P keep the first-encountered row
  # (data.table's order() is stable), which is an arbitrary but deterministic
  # tiebreak -- not specified by the user, flagged here rather than hidden.
  # Must key on CHR *and* BP together, not BP alone -- unlike the old
  # per-chromosome-file version of this script (where each file was
  # necessarily one chromosome, so BP alone was an unambiguous key), a
  # single whole-genome table can have the same BP recur across different
  # chromosomes, and those are not duplicates of each other.
  dt[, P_num := suppressWarnings(as.numeric(P))]
  setorder(dt, CHR, BP, P_num)
  n_before_dedup_by_chr <- dt[, .N, by = CHR]
  dt <- dt[!duplicated(dt, by = c("CHR", "BP"))]
  n_after_dedup_by_chr <- dt[, .N, by = CHR]
  dt[, P_num := NULL]
  dup_removed_by_chr <- merge(n_before_dedup_by_chr, n_after_dedup_by_chr, by = "CHR", all.x = TRUE)
  dup_removed_by_chr[is.na(N.y), N.y := 0L]
  dup_removed_by_chr[, N := N.x - N.y]
  dup_removed_by_chr[, c("N.x", "N.y") := NULL]
} else {
  dup_removed_by_chr <- data.table(CHR = integer(0), N = integer(0))
}

n_out_by_chr <- dt[, .N, by = CHR]

get_n <- function(tbl, chr) { v <- tbl[CHR == chr, N]; if (length(v) == 0) 0L else v }

totals <- list(rows_in = 0L, dup_removed = 0L, hetpva_removed = 0L, hetpva_na_kept = 0L, rows_out = 0L)
for (CHR_i in sort(unique(n_in_by_chr$CHR))) {
  n_in          <- get_n(n_in_by_chr, CHR_i)
  n_dup_removed <- get_n(dup_removed_by_chr, CHR_i)
  n_het_removed <- get_n(n_het_removed_by_chr, CHR_i)
  n_na_kept     <- get_n(n_na_by_chr, CHR_i)
  n_out         <- get_n(n_out_by_chr, CHR_i)

  add_report(sprintf("chr%-2d  in=%7d  dup_removed=%5d  hetpva_removed=%5d  hetpva_na_kept=%6d  out=%7d",
                      CHR_i, n_in, n_dup_removed, n_het_removed, n_na_kept, n_out))

  totals$rows_in        <- totals$rows_in + n_in
  totals$dup_removed     <- totals$dup_removed + n_dup_removed
  totals$hetpva_removed  <- totals$hetpva_removed + n_het_removed
  totals$hetpva_na_kept  <- totals$hetpva_na_kept + n_na_kept
  totals$rows_out        <- totals$rows_out + n_out
}

combined <- dt
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
