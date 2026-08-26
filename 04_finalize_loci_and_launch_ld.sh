#!/bin/bash

#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --job-name=pgc4-finalize-loci
#SBATCH --output=logs/finalize_loci_%A.out
#SBATCH --error=logs/finalize_loci_%A.err

# =============================================================================
# Runs after 01_extract_merge_clump.sh's array completes (submitted with
# --dependency=afterok by submit_full_pipeline.sh). Turns the per-chromosome
# clump reports into loci_input.txt + chunk_lookup.txt, then self-submits
# run_ld_pipeline.sh as a correctly-sized SLURM array -- this is the step
# that used to require watching squeue, counting loci by hand, and
# submitting stage 5 yourself. Once that's submitted, also self-submits
# 05_run_finemapping.sh (--dependency=afterok on the LD array) -- it always
# gets queued; run_finemapping.R itself no-ops immediately if the fine-
# mapping config's method="none", so there's no need to branch on that here.
#
# Usage: sbatch 04_finalize_loci_and_launch_ld.sh <sczvscon|bipvscon> [ref_cohort] [--dataset-config=<path>] [--finemap-config=<path>]
#   ref_cohort defaults to grp10neu3 (see 03_make_chunk_lookup.sh) -- only
#   needs to be a cohort with a complete dosage-chunk directory for whatever
#   dataset --dataset-config points at.
#
# NOTE: known post-processing not automated here (still manual, same as
# before) -- see README.md "Known gotchas": MHC-region loci (chr6:30-33Mb)
# and loci landing exactly on a chunk boundary have historically needed
# manual review/adjustment of loci_input.txt/chunk_lookup.txt before this
# step. If you want those loci flagged automatically, do that before
# resubmitting stage 5.
# =============================================================================

set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DATASET_CONFIG=""
FINEMAP_CONFIG=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --dataset-config=*) DATASET_CONFIG="${arg#--dataset-config=}" ;;
    --finemap-config=*) FINEMAP_CONFIG="${arg#--finemap-config=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"
DATASET_CONFIG="${DATASET_CONFIG:-$SCRIPT_DIR/datasets/ricopili_cross_bcs.sh}"
FINEMAP_CONFIG="${FINEMAP_CONFIG:-$SCRIPT_DIR/finemap_config.R}"

PHENO=${1:?"usage: sbatch 04_finalize_loci_and_launch_ld.sh <sczvscon|bipvscon> [ref_cohort] [--dataset-config=<path>] [--finemap-config=<path>]"}
REF_COHORT=${2:-grp10neu3}

case "$PHENO" in
  sczvscon|bipvscon) ;;
  *) echo "unknown PHENO: $PHENO" >&2; exit 1 ;;
esac

BASE="$FINEMAP_ROOT/$PHENO"
LOCI_FILE="$BASE/loci_input.txt"

bash "$SCRIPT_DIR/02_make_loci_report.sh" "$PHENO"
bash "$SCRIPT_DIR/03_make_chunk_lookup.sh" "$PHENO" "$REF_COHORT" "--dataset-config=${DATASET_CONFIG}"

if [ ! -s "$LOCI_FILE" ]; then
  echo "ERROR: $LOCI_FILE missing/empty after 03_make_chunk_lookup.sh -- nothing to submit for stage 5" >&2
  exit 1
fi

N_LOCI=$(wc -l < "$LOCI_FILE")
if [ "$N_LOCI" -eq 0 ]; then
  echo "ERROR: $LOCI_FILE has 0 loci -- nothing to submit for stage 5" >&2
  exit 1
fi

echo "$N_LOCI loci in $LOCI_FILE -- submitting run_ld_pipeline.sh as an array of that size"

# SLURM does not create --output/--error directories itself -- must exist
# before the array starts writing to them.
mkdir -p "$BASE/output"

LD_JOBID=$(sbatch --parsable \
  --job-name="${PHENO}-ld" \
  --output="$BASE/output/LD_dosage_%A_%a.out" \
  --error="$BASE/output/LD_dosage_%A_%a.err" \
  --array="1-${N_LOCI}%8" \
  "$SCRIPT_DIR/run_ld_pipeline.sh" "$PHENO" "--dataset-config=${DATASET_CONFIG}")

echo "Submitted run_ld_pipeline.sh array: $LD_JOBID (1-${N_LOCI})"

mkdir -p "$BASE/logs"

FINEMAP_JOBID=$(sbatch --parsable \
  --dependency=afterok:$LD_JOBID \
  --job-name="${PHENO}-finemap" \
  --output="$BASE/logs/finemap_%A.out" \
  --error="$BASE/logs/finemap_%A.err" \
  "$SCRIPT_DIR/05_run_finemapping.sh" "$PHENO" "--finemap-config=${FINEMAP_CONFIG}")

echo "Submitted 05_run_finemapping.sh: $FINEMAP_JOBID (depends on LD array $LD_JOBID)"
echo "  (no-ops immediately if the fine-mapping config's method=\"none\")"
echo "Monitor with: squeue -u \$USER"
