#!/bin/bash
# =============================================================================
# One-command orchestrator for the full clump -> loci/chunk -> LD pipeline.
# Chains the stages as SLURM jobs with --dependency=afterok, so they queue
# and run in order automatically without watching squeue and submitting each
# stage by hand -- same idiom as submit_snpvar_pipeline.sh in
# r_files/snellius_pipeline/, not a new pattern:
#
#   qc_filter_sumstats.sh        (single job)  -- sumstats QC + combine
#   01_extract_merge_clump.sh    (22-array)    -- per-chromosome clump
#         |                             |
#         +----------- both must succeed ------+
#                        v
#   04_finalize_loci_and_launch_ld.sh (single job)
#     -> runs 02_make_loci_report.sh + 03_make_chunk_lookup.sh
#     -> self-submits run_ld_pipeline.sh as a correctly-sized array
#     -> self-submits 05_run_finemapping.sh (--dependency=afterok on that
#        array) -- config-driven SuSiE/FINEMAP/both/none, see
#        finemap_config.R. Always queued; no-ops immediately if method="none"
#
# Requires FINEMAP_ROOT to be set in the environment, e.g.:
#   export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping
# sbatch inherits the submitting shell's environment by default, so setting
# it once before running this script propagates automatically through the
# whole self-chaining job tree this submits. See README.md "Configuration".
#
# Usage:
#   bash submit_full_pipeline.sh <sczvscon|bipvscon>
#   bash submit_full_pipeline.sh sczvscon --dedup=off
#   bash submit_full_pipeline.sh sczvscon --hetpva-cutoff=0.01
#   bash submit_full_pipeline.sh sczvscon --hetpva-cutoff=off --dedup=off
#   bash submit_full_pipeline.sh sczvscon --ref-cohort=grp10neu3
#   bash submit_full_pipeline.sh sczvscon --dataset-config=datasets/my_other_dataset.sh
#   bash submit_full_pipeline.sh sczvscon --finemap-config=finemap_config_susie_only.R
#
# Check progress with: squeue -u $USER
# If the QC job or any clump-array task fails, the finalize/LD-launch job is
# cancelled automatically (afterok) rather than running on incomplete input.
# =============================================================================
set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping}"

DATASET_CONFIG=""
FINEMAP_CONFIG=""
DEDUP=""
HETPVA_CUTOFF=""
REF_COHORT=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --dataset-config=*) DATASET_CONFIG="${arg#--dataset-config=}" ;;
    --finemap-config=*) FINEMAP_CONFIG="${arg#--finemap-config=}" ;;
    --dedup=*)          DEDUP="${arg#--dedup=}" ;;
    --hetpva-cutoff=*)  HETPVA_CUTOFF="${arg#--hetpva-cutoff=}" ;;
    --ref-cohort=*)     REF_COHORT="${arg#--ref-cohort=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"

PHENO=${1:?"usage: bash submit_full_pipeline.sh <sczvscon|bipvscon> [--dataset-config=<path>] [--finemap-config=<path>] [--dedup=on|off] [--hetpva-cutoff=<value>|off] [--ref-cohort=<cohort>]"}

case "$PHENO" in
  sczvscon|bipvscon) ;;
  *) echo "unknown PHENO: $PHENO" >&2; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATASET_CONFIG="${DATASET_CONFIG:-$SCRIPT_DIR/datasets/ricopili_cross_bcs.sh}"
FINEMAP_CONFIG="${FINEMAP_CONFIG:-$SCRIPT_DIR/finemap_config.R}"

# SLURM does not create --output/--error directories itself -- must exist
# before any job starts writing to them.
mkdir -p "$FINEMAP_ROOT/logs" "$FINEMAP_ROOT/$PHENO/logs"

QC_ARGS=()
[ -n "$DEDUP" ]         && QC_ARGS+=("--dedup=${DEDUP}")
[ -n "$HETPVA_CUTOFF" ] && QC_ARGS+=("--hetpva-cutoff=${HETPVA_CUTOFF}")

QC_JOBID=$(sbatch --parsable \
  --job-name="${PHENO}-qc" \
  "$SCRIPT_DIR/qc_filter_sumstats.sh" "$PHENO" "${QC_ARGS[@]}")
echo "Submitted sumstats QC:        $QC_JOBID"

CLUMP_JOBID=$(sbatch --parsable \
  --job-name="${PHENO}-clump" \
  --output="$FINEMAP_ROOT/$PHENO/logs/clump_%A_%a.out" \
  --error="$FINEMAP_ROOT/$PHENO/logs/clump_%A_%a.err" \
  "$SCRIPT_DIR/01_extract_merge_clump.sh" "$PHENO" "--dataset-config=${DATASET_CONFIG}")
echo "Submitted clump array:        $CLUMP_JOBID"

FINALIZE_ARGS=("$PHENO")
[ -n "$REF_COHORT" ] && FINALIZE_ARGS+=("$REF_COHORT")
FINALIZE_ARGS+=("--dataset-config=${DATASET_CONFIG}")
FINALIZE_ARGS+=("--finemap-config=${FINEMAP_CONFIG}")

FINALIZE_JOBID=$(sbatch --parsable \
  --dependency=afterok:$CLUMP_JOBID:$QC_JOBID \
  --job-name="${PHENO}-finalize" \
  "$SCRIPT_DIR/04_finalize_loci_and_launch_ld.sh" "${FINALIZE_ARGS[@]}")
echo "Submitted finalize+LD-launch: $FINALIZE_JOBID (depends on clump $CLUMP_JOBID and QC $QC_JOBID)"

echo ""
echo "All stages queued for $PHENO. Monitor with: squeue -u \$USER"
echo "04_finalize_loci_and_launch_ld.sh submits run_ld_pipeline.sh's array, then 05_run_finemapping.sh,"
echo "itself once it runs -- check squeue again after $FINALIZE_JOBID completes to see those job IDs."
