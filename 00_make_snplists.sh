#!/bin/bash

#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=pgc4-snplists
#SBATCH --output=logs/snplists_%A.out
#SBATCH --error=logs/snplists_%A.err

set -euo pipefail

# Extracts the SNP column from each of a phenotype's 22 per-chromosome
# sumstats files, for use as a plink --extract list in stage 01. Lightweight
# awk pass over 22 files -- the small default resource request above is
# intentional, not a placeholder to size up.
#
# Usage: 00_make_snplists.sh <phenotype> --sumstats-prefix=<prefix>
#   <prefix> is whatever your per-chromosome sumstats files are named:
#   expects "$FINEMAP_ROOT/<phenotype>/<prefix>_chr<N>.txt" for N = 1..22.

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/path/to/your/finemapping}"

SUMSTATS_PREFIX=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --sumstats-prefix=*) SUMSTATS_PREFIX="${arg#--sumstats-prefix=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"

PHENO=${1:?"usage: 00_make_snplists.sh <phenotype> --sumstats-prefix=<prefix>"}
SUMSTATS_PREFIX="${SUMSTATS_PREFIX:?--sumstats-prefix=<prefix> is required -- the per-chromosome sumstats filename prefix}"

BASE=$FINEMAP_ROOT/$PHENO
mkdir -p "$BASE/snplists"

for CHR in $(seq 1 22); do
  SUMSTATS="$BASE/${SUMSTATS_PREFIX}_chr${CHR}.txt"
  SNPLIST="$BASE/snplists/${SUMSTATS_PREFIX}_chr${CHR}.snplist"

  if [ ! -f "$SUMSTATS" ]; then
    echo "MISSING: $SUMSTATS" >&2
    continue
  fi

  awk -F'\t' 'NR>1{print $2}' "$SUMSTATS" > "$SNPLIST"
  echo "$PHENO chr$CHR: $(wc -l < "$SNPLIST") SNPs -> $SNPLIST"
done
