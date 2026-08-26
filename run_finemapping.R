#!/usr/bin/env Rscript
# =============================================================================
# Config-driven fine-mapping (SuSiE / FINEMAP / both / none), Snellius-native.
#
# Snellius-side port of scripts/run_finemapping_sczvscon.r (a Windows-local
# script that invoked FINEMAP through WSL). On Snellius, FINEMAP runs
# natively -- no WSL, no path translation -- and this reads the QC'd
# combined sumstats this pipeline's own qc_filter_sumstats.R stage produces,
# rather than re-reading the 22 raw per-chromosome files.
#
# Usage: Rscript run_finemapping.R <sczvscon|bipvscon> [--finemap-config=<path>]
#   Defaults --finemap-config to "finemap_config.R" alongside this script if
#   omitted -- 05_run_finemapping.sh always passes an explicit resolved
#   path, so this default only matters for a manual/standalone invocation.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

# Directory this script itself lives in -- used to locate pipeline_paths.R
# regardless of the caller's working directory (Rscript's own commandArgs
# includes a "--file=<path>" entry, the standard way to self-locate).
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) return(".")
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}
SCRIPT_DIR <- get_script_dir()
source(file.path(SCRIPT_DIR, "pipeline_paths.R"))

args <- commandArgs(trailingOnly = TRUE)

CONFIG_PATH <- ""
POSITIONAL  <- character(0)
for (a in args) {
  if (startsWith(a, "--finemap-config=")) {
    CONFIG_PATH <- sub("^--finemap-config=", "", a)
  } else {
    POSITIONAL <- c(POSITIONAL, a)
  }
}

PHENO <- POSITIONAL[1]
if (is.na(PHENO) || !(PHENO %in% c("sczvscon", "bipvscon"))) {
  stop("usage: Rscript run_finemapping.R <sczvscon|bipvscon> [--finemap-config=<path>]")
}

if (CONFIG_PATH == "") CONFIG_PATH <- file.path(SCRIPT_DIR, "finemap_config.R")
if (!file.exists(CONFIG_PATH)) {
  stop("Fine-mapping config not found: ", CONFIG_PATH,
       " -- pass --finemap-config=<path> or run from the pipeline directory.")
}
source(CONFIG_PATH)  # defines CONFIG
cat("Fine-mapping config loaded from:", CONFIG_PATH, "\n")
cat("method =", CONFIG$method, "\n\n")

if (identical(CONFIG$method, "none")) {
  cat("method = \"none\" -- nothing to do, exiting.\n")
  quit(save = "no", status = 0)
}
if (!(CONFIG$method %in% c("susie", "finemap", "both"))) {
  stop("CONFIG$method must be \"susie\", \"finemap\", \"both\", or \"none\" -- got: ", CONFIG$method)
}

RUN_SUSIE   <- CONFIG$method %in% c("susie", "both")
RUN_FINEMAP <- CONFIG$method %in% c("finemap", "both")

# Fail loudly and early -- before any SuSiE work is wasted -- if FINEMAP is
# needed but Snellius has no FINEMAP module and the user hasn't pointed the
# config at a binary they've placed themselves.
if (RUN_FINEMAP) {
  fbin <- CONFIG$finemap$binary
  if (is.null(fbin) || fbin == "" || !file.exists(fbin) || file.access(fbin, mode = 1) != 0) {
    stop("method=\"", CONFIG$method, "\" requires FINEMAP, but finemap$binary in ", CONFIG_PATH,
         " is blank or not executable (\"", fbin, "\"). ",
         "Snellius has no FINEMAP module -- place the binary in your own scratch/home space ",
         "and set finemap$binary to its full path before running.")
  }
}

if (RUN_SUSIE) suppressPackageStartupMessages(library(susieR))

# ── Paths (shared with generate_reports.R via pipeline_paths.R -- one
# source of truth, not duplicated per script) ───────────────────────────────
P        <- get_pipeline_paths(PHENO)
LDDIR    <- P$lddir
SUSIEDIR <- P$susiedir
FMDIR    <- P$fmdir
SUMSTATS <- P$sumstats
OUT_FILE <- P$finemapping_summary

# ── Configuration (from CONFIG, not hardcoded) ──────────────────────────────
N_CAUSAL              <- CONFIG$max_causal_variants
SKIP_EXISTING         <- isTRUE(CONFIG$skip_existing)
SUSIE_EST_RES_VAR     <- isTRUE(CONFIG$susie$estimate_residual_variance)
SUSIE_LAMBDA          <- CONFIG$susie$lambda
ADAPTIVE_ON           <- isTRUE(CONFIG$adaptive_rerun$enabled)
K_POSTERIOR_THRESHOLD <- CONFIG$adaptive_rerun$k_posterior_threshold
N_CAUSAL_RERUN        <- CONFIG$adaptive_rerun$max_causal_variants_rerun
RERUN_MAX_ITER        <- CONFIG$adaptive_rerun$rerun_max_iter

if (RUN_SUSIE)   dir.create(SUSIEDIR, showWarnings = FALSE, recursive = TRUE)
if (RUN_FINEMAP) dir.create(FMDIR,    showWarnings = FALSE, recursive = TRUE)

if (!file.exists(SUMSTATS)) {
  stop("Combined sumstats not found: ", SUMSTATS, " -- run qc_filter_sumstats.sh first.")
}

# Discover loci from output/ld (whatever run_ld_pipeline.sh's array has
# actually completed so far -- this stage can run against a partial LD
# array, though 04_finalize_loci_and_launch_ld.sh's --dependency=afterok
# normally means the whole array finished first).
snp_log_files <- list.files(LDDIR, pattern = "^[0-9]{3}\\.snp\\.log$", full.names = FALSE)
loci_all      <- sort(sub("\\.snp\\.log$", "", snp_log_files))
cat("Loci found in", LDDIR, ":", length(loci_all), "\n\n")
if (length(loci_all) == 0) {
  stop("No loci found in ", LDDIR, " -- has run_ld_pipeline.sh's array completed?")
}

# =============================================================================
# 1. Load sumstats
# =============================================================================
cat("Loading sumstats...\n")
sumstats <- fread(SUMSTATS)
sumstats[, zscore := log(OR) / SE]

N <- round(mean(4 / (1/sumstats$Nca + 1/sumstats$Nco), na.rm = TRUE))
cat("Effective sample size N =", N, "\n\n")

# =============================================================================
# 2. Per-locus: SuSiE + FINEMAP input prep
# =============================================================================
cat(strrep("=", 60), "\n")
cat("FINE-MAPPING (method =", CONFIG$method, ")\n")
cat(strrep("=", 60), "\n\n")

fm_master_rows <- list()
loci_processed <- character()

for (locus in loci_all) {

  susie_out    <- file.path(SUSIEDIR, paste0(locus, ".susie_l", N_CAUSAL, ".rds"))
  fm_locus_dir <- file.path(FMDIR, locus)
  fm_z_file    <- file.path(fm_locus_dir, paste0(locus, ".z"))
  fm_ld_file   <- file.path(fm_locus_dir, paste0(locus, ".ld"))
  fm_snp_file  <- file.path(fm_locus_dir, paste0(locus, ".snp"))

  susie_done   <- !RUN_SUSIE   || (SKIP_EXISTING && file.exists(susie_out))
  finemap_done <- !RUN_FINEMAP || (SKIP_EXISTING && file.exists(fm_snp_file))

  if (susie_done && finemap_done) {
    cat("[", locus, "] complete, skipping\n")
    if (RUN_FINEMAP) {
      fm_master_rows[[locus]] <- data.table(
        z = fm_z_file, ld = fm_ld_file, snp = fm_snp_file,
        config = file.path(fm_locus_dir, paste0(locus, ".config")),
        cred   = file.path(fm_locus_dir, paste0(locus, ".cred")),
        log    = file.path(fm_locus_dir, paste0(locus, ".log")),
        n_samples = N
      )
    }
    loci_processed <- c(loci_processed, locus)
    next
  }

  snplog_file <- file.path(LDDIR, paste0(locus, ".snp.log"))
  ldfile      <- file.path(LDDIR, paste0(locus, ".ld.gz"))

  if (!file.exists(snplog_file) || !file.exists(ldfile)) {
    cat("[", locus, "] WARNING: LD files not found -- skipping\n")
    next
  }

  ld_snps <- fread(snplog_file, header = FALSE)[[1]]
  LD      <- as.matrix(fread(ldfile, header = FALSE))
  rownames(LD) <- colnames(LD) <- ld_snps
  diag(LD) <- 1.0

  ss_locus    <- sumstats[SNP %in% ld_snps]
  common_snps <- intersect(ld_snps, ss_locus$SNP)

  if (length(common_snps) < 2) {
    cat("[", locus, "] WARNING: <2 overlapping SNPs -- skipping\n")
    next
  }

  ss_locus <- ss_locus[match(common_snps, SNP)]
  LD       <- LD[common_snps, common_snps]

  cat("[", locus, "]", length(common_snps), "SNPs |")

  # ── SuSiE ──
  if (RUN_SUSIE && !susie_done) {
    z   <- ss_locus$zscore
    fit <- tryCatch(
      susie_rss(
        z = z, R = LD, n = N, L = N_CAUSAL,
        lambda = SUSIE_LAMBDA,
        estimate_residual_variance = SUSIE_EST_RES_VAR
      ),
      error = function(e) {
        cat(" SuSiE ERROR:", conditionMessage(e), " |")
        return(NULL)
      }
    )
    if (!is.null(fit)) {
      fit$snp      <- common_snps
      fit$snp_ids  <- common_snps
      fit$zscore   <- z
      fit$position <- ss_locus$BP
      fit$chr      <- ss_locus$CHR[1]
      saveRDS(fit, susie_out)
      cat(sprintf(" SuSiE: %d CS, max_pip=%.3f |", length(fit$sets$cs), max(fit$pip)))
    }
  } else if (RUN_SUSIE) {
    cat(" SuSiE: cached |")
  }

  # ── FINEMAP input prep ──
  if (RUN_FINEMAP && !finemap_done) {
    dir.create(fm_locus_dir, showWarnings = FALSE)

    frq_col <- grep("^FRQ_A_", colnames(ss_locus), value = TRUE)[1]

    z_dt <- data.table(
      rsid       = ss_locus$SNP,
      chromosome = ss_locus$CHR,
      position   = ss_locus$BP,
      allele1    = ss_locus$A1,
      allele2    = ss_locus$A2,
      maf        = pmin(ss_locus[[frq_col]], 1 - ss_locus[[frq_col]]),
      beta       = log(ss_locus$OR),
      se         = ss_locus$SE
    )
    fwrite(z_dt, fm_z_file, sep = " ", quote = FALSE, eol = "\n")

    LD_fm <- LD; diag(LD_fm) <- 1.0
    fwrite(as.data.table(LD_fm), fm_ld_file, sep = " ",
           col.names = FALSE, quote = FALSE, eol = "\n")

    cat(" FINEMAP inputs written |")
  } else if (RUN_FINEMAP) {
    cat(" FINEMAP: cached |")
  }

  if (RUN_FINEMAP) {
    fm_master_rows[[locus]] <- data.table(
      z = fm_z_file, ld = fm_ld_file, snp = fm_snp_file,
      config = file.path(fm_locus_dir, paste0(locus, ".config")),
      cred   = file.path(fm_locus_dir, paste0(locus, ".cred")),
      log    = file.path(fm_locus_dir, paste0(locus, ".log")),
      n_samples = N
    )
  }
  loci_processed <- c(loci_processed, locus)
  cat("\n")
}

# =============================================================================
# 3. Run FINEMAP (master multi-locus batch -- one invocation covers every
#    locus whose output isn't already complete, not one call per locus)
# =============================================================================
master_run_dt <- NULL

if (RUN_FINEMAP) {
  master_dt <- rbindlist(fm_master_rows)

  # Only submit loci whose FINEMAP output is missing or a truncated/empty
  # stub (e.g. left behind by a killed prior run) -- a non-empty .snp file
  # means this locus's --sss output already completed and shouldn't be
  # recomputed.
  snp_incomplete <- vapply(master_dt$snp, function(f)
    !file.exists(f) || length(readLines(f)) == 0, logical(1))
  master_run_dt <- master_dt[snp_incomplete]
  cat(sprintf("\n%d / %d loci already have complete FINEMAP output; running FINEMAP on the remaining %d\n",
              sum(!snp_incomplete), nrow(master_dt), nrow(master_run_dt)))

  master_path <- file.path(FMDIR, "master")
  if (nrow(master_run_dt) > 0) {
    fwrite(master_run_dt, master_path, sep = ";", quote = FALSE, eol = "\n")
    cat("FINEMAP master written:", master_path, "(", nrow(master_run_dt), "loci)\n")
  }

  cat("\n", strrep("=", 60), "\n")
  cat("RUNNING FINEMAP --sss (k =", N_CAUSAL, ")\n")
  cat(strrep("=", 60), "\n\n")

  if (nrow(master_run_dt) == 0) {
    cat("Nothing to run -- all loci already have complete FINEMAP output.\n\n")
  } else {
    finemap_cmd <- sprintf("%s --sss --in-files %s --n-causal-snps %d --log",
                            CONFIG$finemap$binary, master_path, N_CAUSAL)
    cat("Command:", finemap_cmd, "\n\n")
    exit_code <- system(finemap_cmd)
    if (exit_code != 0) {
      cat("WARNING: FINEMAP exited with code", exit_code, "\n\n")
    } else {
      cat("FINEMAP completed.\n\n")
    }
  }
}

# =============================================================================
# 3b. Adaptive re-run (only for whichever method(s) are enabled)
# =============================================================================
load_locus_data <- function(loc) {
  snplog_file <- file.path(LDDIR, paste0(loc, ".snp.log"))
  ldfile      <- file.path(LDDIR, paste0(loc, ".ld.gz"))
  if (!file.exists(snplog_file) || !file.exists(ldfile)) return(NULL)

  ld_snps <- fread(snplog_file, header = FALSE)[[1]]
  LD      <- as.matrix(fread(ldfile, header = FALSE))
  rownames(LD) <- colnames(LD) <- ld_snps
  diag(LD) <- 1.0

  ss_locus    <- sumstats[SNP %in% ld_snps]
  common_snps <- intersect(ld_snps, ss_locus$SNP)
  if (length(common_snps) < 2) return(NULL)

  ss_locus <- ss_locus[match(common_snps, SNP)]
  LD       <- LD[common_snps, common_snps]
  list(ss_locus = ss_locus, LD = LD, common_snps = common_snps)
}

run_susie_locus <- function(loc, dat, L, max_iter = 100) {
  z   <- dat$ss_locus$zscore
  fit <- tryCatch(
    susie_rss(
      z = z, R = dat$LD, n = N, L = L,
      lambda = SUSIE_LAMBDA,
      estimate_residual_variance = SUSIE_EST_RES_VAR,
      max_iter = max_iter
    ),
    error = function(e) { cat(" SuSiE ERROR:", conditionMessage(e)); NULL }
  )
  if (!is.null(fit)) {
    fit$snp      <- dat$common_snps
    fit$snp_ids  <- dat$common_snps
    fit$zscore   <- z
    fit$position <- dat$ss_locus$BP
    fit$chr      <- dat$ss_locus$CHR[1]
    saveRDS(fit, file.path(SUSIEDIR, paste0(loc, ".susie_l", L, ".rds")))
  }
  fit
}

run_finemap_loci <- function(loci_to_run, k) {
  fm_rows <- lapply(loci_to_run, function(loc) {
    fm_locus_dir <- file.path(FMDIR, loc)
    data.table(
      z      = file.path(fm_locus_dir, paste0(loc, ".z")),
      ld     = file.path(fm_locus_dir, paste0(loc, ".ld")),
      snp    = file.path(fm_locus_dir, paste0(loc, ".snp")),
      config = file.path(fm_locus_dir, paste0(loc, ".config")),
      cred   = file.path(fm_locus_dir, paste0(loc, ".cred")),
      log    = file.path(fm_locus_dir, paste0(loc, ".log")),
      n_samples = N
    )
  })
  rerun_dt   <- rbindlist(fm_rows)
  rerun_path <- file.path(FMDIR, paste0("master_rerun_k", k))
  fwrite(rerun_dt, rerun_path, sep = ";", quote = FALSE, eol = "\n")
  cmd <- sprintf("%s --sss --in-files %s --n-causal-snps %d --log",
                 CONFIG$finemap$binary, rerun_path, k)
  cat("  FINEMAP rerun:", cmd, "\n")
  system(cmd)
}

pass1_loci <- character()

if (ADAPTIVE_ON && (RUN_SUSIE || RUN_FINEMAP)) {
  cat(strrep("=", 60), "\n")
  cat("ADAPTIVE RE-RUN: PASS 1\n")
  cat(strrep("=", 60), "\n\n")

  for (loc in loci_processed) {
    needs_rerun <- FALSE; reasons <- character()

    if (RUN_SUSIE) {
      susie_file <- file.path(SUSIEDIR, paste0(loc, ".susie_l", N_CAUSAL, ".rds"))
      if (file.exists(susie_file)) {
        fit <- readRDS(susie_file)
        if (!fit$converged && length(fit$sets$cs) >= N_CAUSAL) {
          needs_rerun <- TRUE
          reasons <- c(reasons, sprintf("SuSiE: %d CS, not converged", length(fit$sets$cs)))
        }
      }
    }

    if (RUN_FINEMAP) {
      log_file <- file.path(FMDIR, loc, paste0(loc, ".log_sss"))
      if (file.exists(log_file)) {
        log_lines  <- readLines(log_file)
        post_start <- grep("Post-Pr\\(# of causal", log_lines)
        post_lines <- if (length(post_start) > 0) log_lines[post_start[1]:length(log_lines)] else character(0)
        k_pattern  <- sprintf("^\\s+%d -> ", N_CAUSAL)
        k_line     <- grep(k_pattern, post_lines, value = TRUE)
        if (length(k_line) > 0) {
          k_prob <- as.numeric(sub(".*-> ", "", trimws(k_line[1])))
          if (!is.na(k_prob) && k_prob > K_POSTERIOR_THRESHOLD) {
            needs_rerun <- TRUE
            reasons <- c(reasons, sprintf("FINEMAP: Pr(k=%d) = %.3g", N_CAUSAL, k_prob))
          }
        }
      }
    }

    if (needs_rerun) {
      cat(sprintf("[%s] SATURATED -- %s\n", loc, paste(reasons, collapse = "; ")))
      pass1_loci <- c(pass1_loci, loc)
    }
  }

  if (length(pass1_loci) > 0) {
    cat(sprintf("\nRe-running %d loci with L/k = %d\n", length(pass1_loci), N_CAUSAL_RERUN))
    if (RUN_SUSIE) {
      for (loc in pass1_loci) {
        dat <- load_locus_data(loc)
        if (is.null(dat)) next
        cat(sprintf("[%s] %d SNPs | L: %d -> %d |", loc, length(dat$common_snps), N_CAUSAL, N_CAUSAL_RERUN))
        fit <- run_susie_locus(loc, dat, L = N_CAUSAL_RERUN)
        if (!is.null(fit)) cat(sprintf(" SuSiE: %d CS, conv=%s", length(fit$sets$cs), fit$converged))
        cat("\n")
      }
    }
    if (RUN_FINEMAP) run_finemap_loci(pass1_loci, k = N_CAUSAL_RERUN)
  } else {
    cat("No saturated loci.\n\n")
  }
}

# =============================================================================
# 3c. Adaptive re-run pass 2 -- SuSiE-iteration-specific, so only applies
#     when SuSiE is one of the enabled methods.
# =============================================================================
if (ADAPTIVE_ON && RUN_SUSIE && length(pass1_loci) > 0) {
  cat(strrep("=", 60), "\n")
  cat("ADAPTIVE RE-RUN: PASS 2\n")
  cat(strrep("=", 60), "\n\n")

  needs_more_iter <- character()
  poorly_resolved <- character()

  for (loc in pass1_loci) {
    susie_file <- file.path(SUSIEDIR, paste0(loc, ".susie_l", N_CAUSAL_RERUN, ".rds"))
    if (!file.exists(susie_file)) next
    fit <- readRDS(susie_file)
    if (fit$converged) {
      cat(sprintf("[%s] Converged at L=%d with %d CS\n", loc, N_CAUSAL_RERUN, length(fit$sets$cs)))
    } else if (length(fit$sets$cs) >= N_CAUSAL_RERUN) {
      cat(sprintf("[%s] HIT CEILING: %d CS at L=%d -- poorly resolved\n",
                  loc, length(fit$sets$cs), N_CAUSAL_RERUN))
      poorly_resolved <- c(poorly_resolved, loc)
    } else {
      cat(sprintf("[%s] NEEDS ITERATIONS: re-running with max_iter=%d\n", loc, RERUN_MAX_ITER))
      needs_more_iter <- c(needs_more_iter, loc)
    }
  }

  if (length(needs_more_iter) > 0) {
    for (loc in needs_more_iter) {
      dat <- load_locus_data(loc)
      if (is.null(dat)) next
      fit <- run_susie_locus(loc, dat, L = N_CAUSAL_RERUN, max_iter = RERUN_MAX_ITER)
      if (!is.null(fit)) cat(sprintf("[%s] %d CS, conv=%s\n", loc, length(fit$sets$cs), fit$converged))
    }
    if (RUN_FINEMAP) run_finemap_loci(needs_more_iter, k = N_CAUSAL_RERUN)
  }

  if (length(poorly_resolved) > 0) {
    cat("\nPOORLY RESOLVED:", paste(poorly_resolved, collapse = ", "), "\n")
    cat("Complex LD -- interpret results with caution.\n\n")
  }

  cat("Adaptive re-run complete.\n\n")
}

# =============================================================================
# 4. Diagnostics
# =============================================================================
cat(strrep("=", 60), "\n")
cat("DIAGNOSTICS\n")
cat(strrep("=", 60), "\n\n")

best_L <- setNames(rep(N_CAUSAL, length(loci_all)), loci_all)
for (loc in pass1_loci) {
  if (file.exists(file.path(SUSIEDIR, paste0(loc, ".susie_l", N_CAUSAL_RERUN, ".rds"))))
    best_L[loc] <- N_CAUSAL_RERUN
}

build_summary <- function() {
  rows <- lapply(loci_all, function(locus) {

    L_used     <- best_L[locus]
    susie_file <- file.path(SUSIEDIR, paste0(locus, ".susie_l", L_used, ".rds"))

    if (file.exists(susie_file)) {
      fit <- readRDS(susie_file)
      cs_lbf <- if (length(fit$sets$cs) > 0)
        sapply(seq_along(fit$sets$cs), function(i)
          round(max(fit$lbf_variable[i, fit$sets$cs[[i]]], na.rm = TRUE), 2))
        else numeric(0)
      cs_size   <- sapply(fit$sets$cs, length)
      cs_purity <- if (!is.null(fit$sets$purity) && nrow(fit$sets$purity) > 0)
        round(fit$sets$purity[, "min.abs.corr"], 3)
        else rep(NA_real_, length(fit$sets$cs))

      susie_dt <- data.table(
        s_L         = L_used,
        s_n_cs      = length(fit$sets$cs),
        s_max_pip   = round(max(fit$pip), 4),
        s_top_snp   = fit$snp[which.max(fit$pip)],
        s_converged = fit$converged,
        s_sigma2    = round(fit$sigma2, 4),
        s_cs_sizes  = paste(cs_size, collapse = ","),
        s_cs_lbf    = paste(cs_lbf, collapse = ","),
        s_cs_purity = paste(cs_purity, collapse = ",")
      )
    } else {
      susie_dt <- data.table(
        s_L = NA_integer_, s_n_cs = NA_integer_, s_max_pip = NA_real_,
        s_top_snp = NA_character_, s_converged = NA,
        s_sigma2 = NA_real_, s_cs_sizes = "", s_cs_lbf = "", s_cs_purity = ""
      )
    }

    fm_locus_dir  <- file.path(FMDIR, locus)
    fm_snp_file   <- file.path(fm_locus_dir, paste0(locus, ".snp"))
    fm_cred_files <- list.files(fm_locus_dir, pattern = paste0("^", locus, "\\.cred\\d+$"),
                                full.names = TRUE)

    if (file.exists(fm_snp_file)) {
      fm_snps <- fread(fm_snp_file)

      log_file   <- file.path(fm_locus_dir, paste0(locus, ".log_sss"))
      post_exp_k <- NA_real_; post_pr_k5 <- NA_real_
      if (file.exists(log_file)) {
        ll         <- readLines(log_file)
        post_start <- grep("Post-Pr\\(# of causal", ll)
        post_lines <- if (length(post_start) > 0) ll[post_start[1]:length(ll)] else character(0)
        exp_line   <- grep("Post-expected", ll, value = TRUE)
        if (length(exp_line) > 0)
          post_exp_k <- as.numeric(regmatches(exp_line[1], regexpr("[0-9\\.]+$", exp_line[1])))
        k5_line <- grep(sprintf("^\\s+%d -> ", N_CAUSAL), post_lines, value = TRUE)
        if (length(k5_line) > 0)
          post_pr_k5 <- as.numeric(sub(".*-> ", "", trimws(k5_line[1])))
      }

      n_cred <- length(fm_cred_files)

      fm_config_file <- file.path(fm_locus_dir, paste0(locus, ".config"))
      config_dt      <- if (file.exists(fm_config_file)) fread(fm_config_file) else NULL
      best_prob      <- if (!is.null(config_dt) && nrow(config_dt) > 0) config_dt[1, prob] else NA_real_

      fm_dt <- data.table(
        f_n_cred    = n_cred,
        f_max_pip   = round(max(fm_snps$prob, na.rm = TRUE), 4),
        f_top_snp   = fm_snps$rsid[which.max(fm_snps$prob)],
        f_exp_k     = round(post_exp_k, 2),
        f_pr_k5     = round(post_pr_k5, 4),
        f_best_prob = round(best_prob, 4)
      )
    } else {
      fm_dt <- data.table(
        f_n_cred = NA_integer_, f_max_pip = NA_real_, f_top_snp = NA_character_,
        f_exp_k = NA_real_, f_pr_k5 = NA_real_, f_best_prob = NA_real_
      )
    }

    cbind(data.table(locus = locus), susie_dt, fm_dt)
  })
  rbindlist(rows, fill = TRUE)
}

summary_dt <- build_summary()

if (RUN_SUSIE) {
  cat("-- SuSiE (L=", N_CAUSAL, ") summary --\n")
  print(summary_dt[, .(locus, s_L, s_n_cs, s_max_pip, s_top_snp, s_converged,
                        s_sigma2, s_cs_lbf, s_cs_purity)])
}

if (RUN_FINEMAP) {
  cat("\n-- FINEMAP (k=", N_CAUSAL, ") summary --\n")
  print(summary_dt[, .(locus, f_n_cred, f_max_pip, f_top_snp, f_exp_k, f_pr_k5, f_best_prob)])
}

if (RUN_SUSIE && RUN_FINEMAP) {
  cat("\n-- Top-SNP agreement --\n")
  both <- summary_dt[!is.na(s_max_pip) & !is.na(f_max_pip)]
  if (nrow(both) > 0) {
    both[, agree := s_top_snp == f_top_snp]
    cat(sprintf("Agreement: %d / %d\n", sum(both$agree, na.rm = TRUE), nrow(both)))
    if (any(!both$agree, na.rm = TRUE))
      print(both[agree == FALSE, .(locus, s_top_snp, s_max_pip, f_top_snp, f_max_pip)])
  }
}

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)
fwrite(summary_dt, OUT_FILE, sep = "\t")
cat("\nSummary written to:", OUT_FILE, "\nDone.\n")
