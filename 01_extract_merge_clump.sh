#!/bin/bash

#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --job-name=pgc4-clump
#SBATCH --output=logs/pgc4_clump_%A_%a.out
#SBATCH --error=logs/pgc4_clump_%A_%a.err
#SBATCH --array=1-22%8

# Requires FINEMAP_ROOT to be set in the environment (e.g.
# `export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping`) -- sbatch inherits
# the submitting shell's environment by default, so exporting it once
# before submitting is enough; see README.md "Configuration".
#
# Override --job-name / --output / --error at submission time, and pass the
# phenotype as a plain script argument (more reliable than --export across
# Slurm configs), e.g.:
#
#   sbatch --job-name=my_pheno-clump \
#     --output=$FINEMAP_ROOT/my_pheno/logs/clump_%A_%a.out \
#     --error=$FINEMAP_ROOT/my_pheno/logs/clump_%A_%a.err \
#     01_extract_merge_clump.sh my_pheno --sumstats-prefix=my_pheno_gwas
#
# Add --dataset-config=<path> to point at a different individual-level
# dataset's config (see datasets/ricopili_cross_bcs.sh for the required
# interface); defaults to that file if omitted, e.g.:
#
#     01_extract_merge_clump.sh my_pheno --sumstats-prefix=my_pheno_gwas --dataset-config=datasets/my_other_dataset.sh

set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CHR=$SLURM_ARRAY_TASK_ID

DATASET_CONFIG=""
SUMSTATS_PREFIX=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --dataset-config=*)   DATASET_CONFIG="${arg#--dataset-config=}" ;;
    --sumstats-prefix=*)  SUMSTATS_PREFIX="${arg#--sumstats-prefix=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"
DATASET_CONFIG="${DATASET_CONFIG:-$SCRIPT_DIR/datasets/ricopili_cross_bcs.sh}"
source "$DATASET_CONFIG"

PHENO=${1:?"usage: sbatch 01_extract_merge_clump.sh <phenotype> --sumstats-prefix=<prefix> [--dataset-config=<path>]"}
SUMSTATS_PREFIX="${SUMSTATS_PREFIX:?--sumstats-prefix=<prefix> is required -- the per-chromosome sumstats filename prefix}"

BASE=$FINEMAP_ROOT/$PHENO
COHORT_LIST=$BASE/EUR_cohorts.txt
PLINK2=$FINEMAP_ROOT/plink/plink2/plink2
PLINK19=$FINEMAP_ROOT/plink/plink1.9/plink

SUMSTATS=$BASE/${SUMSTATS_PREFIX}_chr${CHR}.txt
SNPLIST=$BASE/snplists/${SUMSTATS_PREFIX}_chr${CHR}.snplist

# All bulky intermediates (per-cohort bed, merged bed, merged pgen) live in
# $TMPDIR (node-local scratch, wiped when the job ends). Only the small
# clump report gets written back to the home-quota'd $BASE. This is an
# intentional choice, not just a size limit: the merged whole-chromosome
# pgen is only ever consumed by the --clump call below, and it is not
# needed by the later per-locus fine-mapping steps (those re-extract a
# small, targeted region from the raw cohort bg files once locus
# boundaries are known from the clump report). If that assumption changes,
# add a `cp -r "$WORKDIR/pgen_merged" "$BASE/"` line before the job's
# $TMPDIR is cleaned up.
WORKDIR=$TMPDIR/pgc4_${PHENO}_chr${CHR}
mkdir -p "$WORKDIR/pgen_cohort" "$WORKDIR/pgen_merged" "$BASE/clump" "$BASE/logs"

# Skip-if-done: lets you resubmit the whole 1-22 array after a partial
# failure without redoing chromosomes that already finished successfully.
# Because this script has `set -euo pipefail`, any failure at any earlier
# step (extraction, merge) aborts immediately -- the --clump call is the
# last thing that runs, so chrN.clumps only ever exists if the ENTIRE
# pipeline succeeded for that chromosome. Delete clump/chrN.clumps to
# force that chromosome to rerun even if this check would otherwise skip it.
CLUMP_DONE=$BASE/clump/chr${CHR}.clumps
if [ -s "$CLUMP_DONE" ]; then
  echo "chr$CHR already completed ($CLUMP_DONE exists), skipping. Delete that file to force a rerun."
  exit 0
fi

if [ ! -s "$SNPLIST" ]; then
  echo "SNP list missing: $SNPLIST (run 00_make_snplists.sh first)" >&2
  exit 1
fi

MERGELIST=$WORKDIR/pgen_merged/chr${CHR}_mergelist.txt
: > "$MERGELIST"

while read -r COHORT; do
  [ -z "$COHORT" ] && continue

  mapfile -t BG_MATCHES < <(dataset_find_cohort_bed "$COHORT")

  if [ "${#BG_MATCHES[@]}" -eq 0 ]; then
    echo "WARNING: no .bg.bed found for cohort $COHORT (dataset_find_cohort_bed, config: $DATASET_CONFIG)" >&2
    continue
  fi
  if [ "${#BG_MATCHES[@]}" -gt 1 ]; then
    echo "WARNING: multiple .bg.bed matches for cohort $COHORT, using first: ${BG_MATCHES[0]}" >&2
  fi

  BGPREFIX=${BG_MATCHES[0]%.bed}
  OUT=$WORKDIR/pgen_cohort/${COHORT}_chr${CHR}

  # NOTE: --make-bed, not --make-pgen. plink2's --pmerge-list only supports
  # "concatenating" merges (same samples, different variants). Our case is
  # the opposite (same variants, different samples per cohort), which plink2
  # explicitly does not support yet ("Non-concatenating --pmerge[-list] is
  # under development"). So the sample-merge is done with plink1.9's
  # --merge-list below, which requires bed/bim/fam input.
  "$PLINK2" --bfile "$BGPREFIX" \
    --chr "$CHR" \
    --extract "$SNPLIST" \
    --rm-dup exclude-all \
    --make-bed \
    --memory 18000 \
    --out "$OUT"

  if [ -f "${OUT}.bed" ]; then
    echo "$OUT" >> "$MERGELIST"
  else
    echo "WARNING: plink2 extract produced no .bed for cohort $COHORT chr $CHR" >&2
  fi
done < "$COHORT_LIST"

if [ ! -s "$MERGELIST" ]; then
  echo "ERROR: no cohort .bed files produced for chr $CHR, aborting merge/clump" >&2
  exit 1
fi

MERGED_BED=$WORKDIR/pgen_merged/chr${CHR}_merged_bed
MERGED=$WORKDIR/pgen_merged/chr${CHR}_merged

BASE_COHORT=$(head -n1 "$MERGELIST")
REST_LIST=$WORKDIR/pgen_merged/chr${CHR}_mergelist_rest.txt
tail -n +2 "$MERGELIST" > "$REST_LIST"

if [ -s "$REST_LIST" ]; then
  "$PLINK19" --bfile "$BASE_COHORT" \
    --merge-list "$REST_LIST" \
    --allow-no-sex \
    --make-bed \
    --memory 18000 \
    --out "$MERGED_BED"
else
  # only one cohort had data for this chromosome, nothing to merge
  "$PLINK19" --bfile "$BASE_COHORT" \
    --allow-no-sex \
    --make-bed \
    --memory 18000 \
    --out "$MERGED_BED"
fi

# Free the per-cohort bed files now -- the merge is done, they're no longer
# needed, and on a big chromosome (chr1/chr2, ~40 cohorts) they're the
# single largest thing sitting in $TMPDIR at once.
rm -rf "$WORKDIR/pgen_cohort"

"$PLINK2" --bfile "$MERGED_BED" \
  --make-pgen \
  --memory 18000 \
  --out "$MERGED"

# Merged bed is superseded by the merged pgen now; drop it before clumping.
rm -f "${MERGED_BED}".bed "${MERGED_BED}".bim "${MERGED_BED}".fam "${MERGED_BED}".log

CLUMPOUT=$BASE/clump/chr${CHR}

"$PLINK2" --pfile "$MERGED" \
  --clump "$SUMSTATS" \
  --clump-p1 5e-8 \
  --clump-p2 1e-4 \
  --clump-r2 0.1 \
  --clump-kb 3000 \
  --clump-id-field SNP \
  --clump-field P \
  --memory 18000 \
  --out "$CLUMPOUT"
