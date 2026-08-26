#!/bin/bash

#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --job-name=pgc4-reports
#SBATCH --output=logs/reports_%A.out
#SBATCH --error=logs/reports_%A.err

# =============================================================================
# Reports stage: per-SNP results table, genome-wide summary-stats rollup,
# PIP-colored locus-zoom scatter + LD heatmap per locus -- see
# generate_reports.R. No VEP/gene annotation.
#
# Existence-driven, not config-driven -- reads whatever SuSiE/FINEMAP output
# is actually on disk, so it doesn't need --finemap-config at all (unlike
# stages 01/03/05, which need to know the *dataset*/fine-mapping config
# before they run something; this stage only ever reads what's already
# there).
#
# Usage: sbatch 06_generate_reports.sh <sczvscon|bipvscon>
# =============================================================================

module load 2025
module load R/4.5.1-gfbf-2025a

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PHENO=${1:?"usage: sbatch 06_generate_reports.sh <sczvscon|bipvscon>"}

case "$PHENO" in
  sczvscon|bipvscon) ;;
  *) echo "unknown PHENO: $PHENO" >&2; exit 1 ;;
esac

Rscript --vanilla "$SCRIPT_DIR/generate_reports.R" "$PHENO"
