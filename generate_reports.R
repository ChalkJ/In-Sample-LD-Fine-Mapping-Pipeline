#!/usr/bin/env Rscript
# =============================================================================
# Reports stage: per-SNP results table, genome-wide summary-stats rollup, and
# a PIP-colored locus-zoom scatter + LD heatmap per locus. No VEP/gene
# annotation -- deliberately out of scope (see README.md "Reports stage").
#
# Existence-driven, not config-driven: reads whatever SuSiE/FINEMAP output is
# actually on disk under this phenotype's susie/finemap dirs, regardless of
# what finemap_config.R's method= currently says -- so re-running this after
# the config changed (or against a partial run) still reflects reality
# rather than a stale assumption. No --finemap-config flag needed as a
# result.
#
# Usage: Rscript generate_reports.R <sczvscon|bipvscon>
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) return(".")
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}
SCRIPT_DIR <- get_script_dir()
source(file.path(SCRIPT_DIR, "pipeline_paths.R"))

args  <- commandArgs(trailingOnly = TRUE)
PHENO <- args[1]
if (is.na(PHENO) || !(PHENO %in% c("sczvscon", "bipvscon"))) {
  stop("usage: Rscript generate_reports.R <sczvscon|bipvscon>")
}

P <- get_pipeline_paths(PHENO)
dir.create(P$resultsdir, showWarnings = FALSE, recursive = TRUE)
dir.create(P$plotsdir,   showWarnings = FALSE, recursive = TRUE)

if (!file.exists(P$sumstats)) {
  stop("Combined sumstats not found: ", P$sumstats, " -- run qc_filter_sumstats.sh first.")
}

loci_input_file <- file.path(P$datadir, "loci_input.txt")
if (!file.exists(loci_input_file)) {
  stop("loci_input.txt not found: ", loci_input_file, " -- run 03_make_chunk_lookup.sh first.")
}
loci_input <- fread(loci_input_file, header = FALSE,
                     col.names = c("loci_id", "chr", "start", "end"))

snp_log_files <- list.files(P$lddir, pattern = "^[0-9]{3}\\.snp\\.log$", full.names = FALSE)
loci_all      <- sort(sub("\\.snp\\.log$", "", snp_log_files))
cat("Loci found in", P$lddir, ":", length(loci_all), "\n")
if (length(loci_all) == 0) stop("No loci found in ", P$lddir, " -- has the LD stage completed?")

n_susie   <- length(list.files(P$susiedir, pattern = "\\.susie_l[0-9]+\\.rds$"))
n_finemap <- length(list.files(P$fmdir, recursive = TRUE, pattern = "\\.snp$"))
cat("SuSiE output files found:", n_susie, " | FINEMAP output files found:", n_finemap, "\n\n")

cat("Loading sumstats...\n")
sumstats <- fread(P$sumstats)
sumstats[, zscore := log(OR) / SE]

# =============================================================================
# Per-locus "best output" resolution
# =============================================================================

# SuSiE reruns write a NEW L-suffixed file rather than overwriting (e.g. both
# NNN.susie_l5.rds and NNN.susie_l10.rds can exist for a locus that saturated
# and got re-run) -- the rerun only ever moves L up, so the highest-L file
# present is always the final fit.
find_best_susie_file <- function(susiedir, locus) {
  files <- list.files(susiedir, pattern = paste0("^", locus, "\\.susie_l([0-9]+)\\.rds$"), full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  Ls <- as.integer(sub(paste0("^", locus, "\\.susie_l([0-9]+)\\.rds$"), "\\1", basename(files)))
  files[which.max(Ls)]
}

# FINEMAP's .credK files are NOT "the Nth credible set" -- each is a
# separate hypothesis assuming exactly K causal SNPs, and a locus can have
# .cred1, .cred2, ... simultaneously on disk. The correct one is whichever
# has the highest posterior probability in its own header line, not the
# highest K (a real bug hit and fixed once already on this project, in
# summarize_polyfun_tables.R -- ported here from master_results.Rmd's
# make_fm_plot(), not re-derived from scratch).
resolve_finemap_cred <- function(fm_locus_dir, locus) {
  cred_files <- sort(list.files(fm_locus_dir, pattern = paste0("^", locus, "\\.cred[0-9]+$"), full.names = TRUE))
  best_cred <- NULL; best_pr <- -1
  for (cf in cred_files) {
    hdr <- readLines(cf, n = 2)
    pr  <- suppressWarnings(as.numeric(regmatches(hdr[1], regexpr("[0-9.eE+-]+$", hdr[1]))))
    if (length(pr) > 0 && !is.na(pr) && pr > best_pr) { best_pr <- pr; best_cred <- cf }
  }
  cs_map <- setNames(integer(0), character(0))
  if (!is.null(best_cred)) {
    clines <- readLines(best_cred)
    skip   <- grep("^index", clines) - 1
    if (length(skip) > 0) {
      cdat <- tryCatch(fread(best_cred, skip = skip[1]), error = function(e) NULL)
      if (!is.null(cdat)) {
        cs_cols <- grep("^cred[0-9]+$", names(cdat), value = TRUE)
        for (i in seq_along(cs_cols)) {
          snps <- cdat[[cs_cols[i]]]
          snps <- snps[!is.na(snps) & snps != ""]
          if (length(snps) > 0) cs_map[snps] <- i
        }
      }
    }
  }
  list(cred_file = best_cred, cs_map = cs_map)
}

# =============================================================================
# Per-locus SNP table (sumstats + PIPs + CS membership) -- also the direct
# input to the plots below, so PIP/CS lookup logic isn't duplicated between
# the table and the plotting code.
# =============================================================================
build_snp_table <- function(locus) {
  snplog_file <- file.path(P$lddir, paste0(locus, ".snp.log"))
  if (!file.exists(snplog_file)) return(NULL)

  ld_snps     <- fread(snplog_file, header = FALSE)[[1]]
  ss_locus    <- sumstats[SNP %in% ld_snps]
  common_snps <- intersect(ld_snps, ss_locus$SNP)
  if (length(common_snps) < 2) return(NULL)
  ss_locus <- ss_locus[match(common_snps, SNP)]

  frq_a_col <- grep("^FRQ_A_", colnames(ss_locus), value = TRUE)[1]
  frq_u_col <- grep("^FRQ_U_", colnames(ss_locus), value = TRUE)[1]

  dt <- data.table(
    locus  = locus,
    CHR    = ss_locus$CHR, SNP = ss_locus$SNP, BP = ss_locus$BP,
    A1     = ss_locus$A1,  A2  = ss_locus$A2,
    OR     = ss_locus$OR,  SE  = ss_locus$SE, P = ss_locus$P,
    FRQ_A  = ss_locus[[frq_a_col]], FRQ_U = ss_locus[[frq_u_col]],
    zscore = ss_locus$zscore
  )

  susie_file <- find_best_susie_file(P$susiedir, locus)
  if (!is.na(susie_file)) {
    fit     <- readRDS(susie_file)
    pip_map <- setNames(fit$pip, fit$snp_ids)
    cs_map  <- setNames(rep(NA_integer_, length(fit$snp_ids)), fit$snp_ids)
    if (length(fit$sets$cs) > 0) {
      for (i in seq_along(fit$sets$cs)) cs_map[fit$snp_ids[fit$sets$cs[[i]]]] <- i
    }
    dt[, susie_pip := unname(pip_map[SNP])]
    dt[, susie_cs  := unname(cs_map[SNP])]
  } else {
    dt[, susie_pip := NA_real_]
    dt[, susie_cs  := NA_integer_]
  }

  fm_locus_dir <- file.path(P$fmdir, locus)
  fm_snp_file  <- file.path(fm_locus_dir, paste0(locus, ".snp"))
  if (file.exists(fm_snp_file)) {
    fm_snps <- fread(fm_snp_file)
    pip_map <- setNames(fm_snps$prob, fm_snps$rsid)
    cred    <- resolve_finemap_cred(fm_locus_dir, locus)
    dt[, finemap_pip := unname(pip_map[SNP])]
    dt[, finemap_cs  := unname(cred$cs_map[SNP])]
  } else {
    dt[, finemap_pip := NA_real_]
    dt[, finemap_cs  := NA_integer_]
  }

  dt
}

# =============================================================================
# Plots -- PIP-colored scatter (no gene track) + LD heatmap. Both skip if
# their PNG already exists (resumable, matches master_results.Rmd's
# embed_gg() idiom: a killed/re-run job never redraws finished plots).
# =============================================================================
plot_pip_scatter <- function(snp_dt, pip_col, cs_col, locus, chr_val, method_label, out_png) {
  if (file.exists(out_png)) return(invisible(NULL))

  d <- data.table(BP = snp_dt$BP, P = snp_dt$P,
                   PIP = snp_dt[[pip_col]], cs = snp_dt[[cs_col]])
  d[, log10p := -log10(P)]
  d[, in_cs  := !is.na(cs)]
  d[is.na(PIP), PIP := 0]

  p <- ggplot(d, aes(x = BP, y = log10p)) +
    geom_point(data = d[!(in_cs)], aes(size = PIP), color = "#BBBBBB", alpha = 0.5) +
    geom_point(data = d[(in_cs)],  aes(color = PIP, size = PIP, alpha = PIP)) +
    scale_color_gradient(low = "#5D3A9B", high = "#E66100", limits = c(0, 1), name = "PIP") +
    scale_size_continuous(range  = c(0.5, 5), limits = c(0, 1)) +
    scale_alpha_continuous(range = c(0.4, 1), limits = c(0, 1)) +
    scale_x_continuous(labels = scales::comma) +
    labs(x = "Position", y = expression(-log[10](italic(p))),
         title = sprintf("%s — Locus %s  Chr%s", method_label, locus, chr_val)) +
    guides(size = "none", alpha = "none") +
    theme_bw()

  ggsave(out_png, p, width = 8, height = 4.5, dpi = 150)
}

plot_ld_heatmap <- function(locus, out_png) {
  if (file.exists(out_png)) return(invisible(NULL))
  ld_file <- file.path(P$lddir, paste0(locus, ".ld.gz"))
  if (!file.exists(ld_file)) return(invisible(NULL))

  LD <- as.matrix(fread(ld_file, header = FALSE))
  n  <- nrow(LD)
  ld_pal <- colorRampPalette(c("#313695", "#74add1", "#e0f3f8", "#ffffbf",
                                "#fee090", "#f46d43", "#a50026"))(256)

  png(out_png, width = 900, height = 800, res = 120)
  old <- par(no.readonly = TRUE)
  on.exit({ par(old); dev.off() })
  layout(matrix(c(1L, 2L), 1, 2), widths = c(8, 1))
  par(mar = c(2, 2, 3, 0.5), pty = "s")
  image(seq_len(n), seq_len(n), LD, col = ld_pal, zlim = c(-1, 1),
        axes = FALSE, useRaster = TRUE,
        main = sprintf("LD r — locus %s  (%d x %d SNPs)", locus, n, n))
  par(mar = c(2, 0.5, 3, 2.5), pty = "m")
  image(1L, seq(-1, 1, length.out = 256),
        matrix(seq(-1, 1, length.out = 256), nrow = 1L),
        col = ld_pal, zlim = c(-1, 1), xaxt = "n", yaxt = "n")
  axis(4, at = c(-1, 0, 1), labels = c("-1", "0", "1"), las = 1, cex.axis = 0.85)
  mtext("r", side = 4, line = 2, cex = 0.85)
}

# =============================================================================
# Main loop: build the per-SNP table and draw plots for every locus
# =============================================================================
cat(strrep("=", 60), "\n")
cat("BUILDING RESULTS TABLE + PLOTS\n")
cat(strrep("=", 60), "\n\n")

snp_rows <- list()

for (locus in loci_all) {
  snp_dt <- build_snp_table(locus)
  if (is.null(snp_dt)) {
    cat("[", locus, "] WARNING: could not build SNP table -- skipping\n")
    next
  }
  snp_rows[[locus]] <- snp_dt

  li      <- loci_input[loci_id == as.integer(locus)]
  chr_val <- if (nrow(li) > 0) li$chr[1] else NA

  has_susie   <- any(!is.na(snp_dt$susie_pip))
  has_finemap <- any(!is.na(snp_dt$finemap_pip))

  cat("[", locus, "]", nrow(snp_dt), "SNPs | susie:", has_susie, "| finemap:", has_finemap, "\n")

  if (has_susie) {
    plot_pip_scatter(snp_dt, "susie_pip", "susie_cs", locus, chr_val, "SuSiE",
                      file.path(P$plotsdir, paste0(locus, "_susie_scatter.png")))
  }
  if (has_finemap) {
    plot_pip_scatter(snp_dt, "finemap_pip", "finemap_cs", locus, chr_val, "FINEMAP",
                      file.path(P$plotsdir, paste0(locus, "_finemap_scatter.png")))
  }
  plot_ld_heatmap(locus, file.path(P$plotsdir, paste0(locus, "_ld_heatmap.png")))
}

all_snps <- rbindlist(snp_rows)
all_snps_file <- file.path(P$resultsdir, paste0("all_snps_", PHENO, ".tsv"))
fwrite(all_snps, all_snps_file, sep = "\t")
cat("\nWrote", all_snps_file, "(", nrow(all_snps), "rows)\n")

# =============================================================================
# Genome-wide summary-stats rollup
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("SUMMARY STATS\n")
cat(strrep("=", 60), "\n\n")

summary_rows <- list(n_loci = length(loci_all))

if (file.exists(P$finemapping_summary)) {
  locus_summary <- fread(P$finemapping_summary)

  susie_present   <- !is.na(locus_summary$s_max_pip)
  finemap_present <- !is.na(locus_summary$f_max_pip)

  summary_rows$n_loci_with_susie   <- sum(susie_present)
  summary_rows$n_loci_with_finemap <- sum(finemap_present)
  summary_rows$mean_max_pip_susie    <- if (any(susie_present))   round(mean(locus_summary$s_max_pip[susie_present]), 4) else NA_real_
  summary_rows$median_max_pip_susie  <- if (any(susie_present))   round(median(locus_summary$s_max_pip[susie_present]), 4) else NA_real_
  summary_rows$mean_max_pip_finemap  <- if (any(finemap_present)) round(mean(locus_summary$f_max_pip[finemap_present]), 4) else NA_real_
  summary_rows$median_max_pip_finemap <- if (any(finemap_present)) round(median(locus_summary$f_max_pip[finemap_present]), 4) else NA_real_
} else {
  cat("NOTE:", P$finemapping_summary, "not found -- max-PIP rollups skipped (run 05 first).\n")
}

# CS-size rollups come from the freshly-built all_snps table (actual SNP
# membership per credible set), not from locus_summary's comma-joined
# s_cs_sizes string -- avoids re-parsing that string and gives FINEMAP a
# real CS-size rollup too (locus_summary only has FINEMAP's credible-set
# *count*, f_n_cred, not their sizes).
add_cs_rollup <- function(prefix, cs_col) {
  if (!(cs_col %in% names(all_snps))) return(invisible(NULL))
  has_cs <- all_snps[!is.na(get(cs_col))]
  if (nrow(has_cs) == 0) {
    summary_rows[[paste0("mean_cs_size_", prefix)]] <<- NA_real_
    summary_rows[[paste0("mean_n_cs_", prefix)]]     <<- 0
    summary_rows[[paste0("n_loci_zero_cs_", prefix)]] <<- length(loci_all)
    return(invisible(NULL))
  }
  cs_sizes <- has_cs[, .N, by = .(locus, cs = get(cs_col))]
  n_cs_per_locus <- has_cs[, uniqueN(get(cs_col)), by = locus]

  summary_rows[[paste0("mean_cs_size_", prefix)]]  <<- round(mean(cs_sizes$N), 2)
  summary_rows[[paste0("mean_n_cs_", prefix)]]      <<- round(mean(n_cs_per_locus$V1), 2)
  summary_rows[[paste0("n_loci_zero_cs_", prefix)]] <<- length(loci_all) - uniqueN(has_cs$locus)
}
add_cs_rollup("susie", "susie_cs")
add_cs_rollup("finemap", "finemap_cs")

summary_dt <- as.data.table(summary_rows)
summary_file <- file.path(P$resultsdir, paste0("summary_stats_", PHENO, ".tsv"))
fwrite(summary_dt, summary_file, sep = "\t")
print(summary_dt)
cat("\nWrote", summary_file, "\n")

cat("\nDone.\n")
