#!/bin/bash

#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --job-name=pgc4-finemap
#SBATCH --output=logs/finemap_%A.out
#SBATCH --error=logs/finemap_%A.err

# =============================================================================
# Fine-mapping stage: config-driven SuSiE / FINEMAP / both / none, single job
# (not an array) -- see run_finemapping.R and finemap_config.R. Single job
# because FINEMAP's --in-files batches every locus into one invocation and
# the adaptive-rerun passes reason about saturated loci across the whole
# set; a per-locus array would need those restructured.
#
# Generous default time limit: this loops sequentially over every locus and
# may include FINEMAP re-runs -- same risk profile as the MHC-region loci
# that hit run_ld_pipeline.sh's 8hr array limit. Override at submission
# time same as the other stages, e.g.:
#
#   sbatch --job-name=my_pheno-finemap \
#     --output=$FINEMAP_ROOT/my_pheno/logs/finemap_%A.out \
#     --error=$FINEMAP_ROOT/my_pheno/logs/finemap_%A.err \
#     --time=48:00:00 \
#     05_run_finemapping.sh my_pheno
#
# Add --finemap-config=<path> to use a config other than the default
# finemap_config.R alongside this script (method toggle, max_causal_variants,
# FINEMAP binary path -- see that file's comments).
#
# Once run_finemapping.R succeeds, self-submits 06_generate_reports.sh (no
# --dependency needed -- that submission only happens after the Rscript call
# below has already completed, thanks to set -euo pipefail below; unlike
# 04->05, there's nothing async left to wait on here).
# =============================================================================

module load 2025
module load R/4.5.1-gfbf-2025a

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FINEMAP_CONFIG=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --finemap-config=*) FINEMAP_CONFIG="${arg#--finemap-config=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"
FINEMAP_CONFIG="${FINEMAP_CONFIG:-$SCRIPT_DIR/finemap_config.R}"

PHENO=${1:?"usage: sbatch 05_run_finemapping.sh <phenotype> [--finemap-config=<path>]"}

Rscript --vanilla "$SCRIPT_DIR/run_finemapping.R" "$PHENO" "--finemap-config=${FINEMAP_CONFIG}"

# generate_reports.R is existence-driven (reads whatever SuSiE/FINEMAP
# output is actually on disk), so it doesn't need --finemap-config at all.
REPORTS_JOBID=$(sbatch --parsable \
  --job-name="${PHENO}-reports" \
  "$SCRIPT_DIR/06_generate_reports.sh" "$PHENO")

echo "Submitted 06_generate_reports.sh: $REPORTS_JOBID"
