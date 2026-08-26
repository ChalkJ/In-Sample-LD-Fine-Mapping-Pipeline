#!/bin/bash

#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --job-name=pgc4-qc-sumstats
#SBATCH --output=logs/qc_sumstats_%A.out
#SBATCH --error=logs/qc_sumstats_%A.err

# =============================================================================
# Thin SBATCH wrapper around qc_filter_sumstats.R -- dedups per-chromosome
# sumstats (keep lowest P per CHR:BP) and filters HetPVa, then combines all
# 22 chromosomes into the single daner-format .gz run_ld_pipeline.sh expects.
# --sumstats-prefix=<prefix> is required; both filters are optional/tunable.
# Any extra args after PHENO are forwarded straight through to the R script,
# e.g.:
#
#   sbatch --job-name=my_pheno-qc \
#     --output=$FINEMAP_ROOT/my_pheno/logs/qc_%A.out \
#     --error=$FINEMAP_ROOT/my_pheno/logs/qc_%A.err \
#     qc_filter_sumstats.sh my_pheno --sumstats-prefix=my_pheno_gwas
#
#   sbatch qc_filter_sumstats.sh my_pheno --sumstats-prefix=my_pheno_gwas --dedup=off
#   sbatch qc_filter_sumstats.sh my_pheno --sumstats-prefix=my_pheno_gwas --hetpva-cutoff=0.01
#   sbatch qc_filter_sumstats.sh my_pheno --sumstats-prefix=my_pheno_gwas --dedup=off --hetpva-cutoff=off
# =============================================================================

module load 2025
module load R/4.5.1-gfbf-2025a

set -euo pipefail
PHENO=${1:?"usage: sbatch qc_filter_sumstats.sh <phenotype> --sumstats-prefix=<prefix> [--dedup=on|off] [--hetpva-cutoff=<value>|off]"}
shift

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Rscript --vanilla "$SCRIPT_DIR/qc_filter_sumstats.R" "$PHENO" "$@"
