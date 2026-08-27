#!/bin/bash

#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --job-name=pgc4-snplists
#SBATCH --output=logs/snplists_%A.out
#SBATCH --error=logs/snplists_%A.err

set -euo pipefail

# Extracts the SNP column from the phenotype's whole-genome sumstats file,
# for use as a plink --extract list in stage 01. A single, genome-wide list
# is sufficient -- 01_extract_merge_clump.sh already restricts genotype data
# with --chr per array task before applying --extract, so this list never
# needed to be split by chromosome.
#
# Usage: 00_make_snplists.sh <phenotype> --sumstats-file=<filename>
#   <filename> is your whole-genome sumstats file, expected directly under
#   "$FINEMAP_ROOT/<phenotype>/".

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/path/to/your/finemapping}"

SUMSTATS_FILE=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --sumstats-file=*) SUMSTATS_FILE="${arg#--sumstats-file=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"

PHENO=${1:?"usage: 00_make_snplists.sh <phenotype> --sumstats-file=<filename>"}
SUMSTATS_FILE="${SUMSTATS_FILE:?--sumstats-file=<filename> is required -- your whole-genome sumstats file}"

BASE=$FINEMAP_ROOT/$PHENO
SUMSTATS="$BASE/$SUMSTATS_FILE"
SNPLIST="$BASE/snplist.txt"

if [ ! -f "$SUMSTATS" ]; then
  echo "ERROR: $SUMSTATS not found" >&2
  exit 1
fi

awk -F'\t' 'NR>1{print $2}' "$SUMSTATS" > "$SNPLIST"
echo "$PHENO: $(wc -l < "$SNPLIST") SNPs -> $SNPLIST"
