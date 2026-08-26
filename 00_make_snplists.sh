#!/bin/bash
set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping}"

declare -A PREFIXES=( ["sczvscon"]="sc_vs_allcontrols" ["bipvscon"]="bp_vs_allcontrols" )

for PHENO in sczvscon bipvscon; do
  PREFIX=${PREFIXES[$PHENO]}
  BASE=$FINEMAP_ROOT/$PHENO
  mkdir -p "$BASE/snplists"

  for CHR in $(seq 1 22); do
    SUMSTATS="$BASE/${PREFIX}_chr${CHR}.txt"
    SNPLIST="$BASE/snplists/${PREFIX}_chr${CHR}.snplist"

    if [ ! -f "$SUMSTATS" ]; then
      echo "MISSING: $SUMSTATS" >&2
      continue
    fi

    awk -F'\t' 'NR>1{print $2}' "$SUMSTATS" > "$SNPLIST"
    echo "$PHENO chr$CHR: $(wc -l < "$SNPLIST") SNPs -> $SNPLIST"
  done
done
