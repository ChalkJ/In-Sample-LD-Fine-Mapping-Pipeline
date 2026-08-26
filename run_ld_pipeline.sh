#!/bin/bash

#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --job-name=LD-dosage
#SBATCH --output=logs/LD_dosage_%A_%a.out
#SBATCH --error=logs/LD_dosage_%A_%a.err
#SBATCH --array=1-10%8

# Requires FINEMAP_ROOT to be set in the environment -- see README.md
# "Configuration". PHENO is passed as a script argument, not --export (see
# prior pipeline scripts for why). Override --job-name/--output/--error/
# --array at submission time, e.g.:
#
#   sbatch --job-name=my_pheno-ld \
#     --output=$FINEMAP_ROOT/my_pheno/output/LD_dosage_%A_%a.out \
#     --error=$FINEMAP_ROOT/my_pheno/output/LD_dosage_%A_%a.err \
#     --array=1-273%8 \
#     run_ld_pipeline.sh my_pheno
#
# (273 = number of loci in my_pheno's loci_input.txt, i.e.
#  `wc -l < loci_input.txt` after running 03_make_chunk_lookup.sh, since the
#  array bound can't be computed dynamically inside the script itself --
#  04_finalize_loci_and_launch_ld.sh computes and passes this automatically
#  when it self-submits this script.)

module load 2025
module load R/4.5.1-gfbf-2025a

set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/path/to/your/finemapping}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Add --dataset-config=<path> to point at a different individual-level
# dataset's config (see datasets/ricopili_cross_bcs.sh); defaults to that
# file if omitted.
DATASET_CONFIG=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --dataset-config=*) DATASET_CONFIG="${arg#--dataset-config=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"
DATASET_CONFIG="${DATASET_CONFIG:-$SCRIPT_DIR/datasets/ricopili_cross_bcs.sh}"
# Resolve to an absolute path now -- this script cd's into $TMPDIR below, and
# make_dataset_ld.sh (run from there) needs to source this same file.
case "$DATASET_CONFIG" in
  /*) ;;
  *) DATASET_CONFIG="$SCRIPT_DIR/$DATASET_CONFIG" ;;
esac

PHENO=${1:?"usage: sbatch run_ld_pipeline.sh <phenotype> [--dataset-config=<path>]"}

# --- Paths ---
# Combined sumstats filename is a fixed convention -- see qc_filter_sumstats.R
# and pipeline_paths.R, which both write/read this same name; not phenotype-
# specific beyond substituting PHENO itself.
DATADIR="$FINEMAP_ROOT/$PHENO"
SUMSTATS="$DATADIR/daner_${PHENO}_qc.gz"
PLINK2="$FINEMAP_ROOT/plink/plink2/plink2"
PLINK1="$FINEMAP_ROOT/plink/plink1.9/plink"
OUTDIR="$DATADIR/output/ld"
LOCIFILE="loci_input.txt"
COHORTLIST="EUR_cohorts.txt"
CHUNKLOOKUP="chunk_lookup.txt"

mkdir -p "$OUTDIR"

# --- Copy files to scratch ---
cd "$TMPDIR"
cp "$DATADIR/make_dataset_ld.sh" .
cp "$DATADIR/LDmerge_v2.R" .
cp "$DATADIR/$LOCIFILE" .          # loci_input_new.txt
cp "$DATADIR/$COHORTLIST" .
cp "$DATADIR/$CHUNKLOOKUP" .
chmod +x make_dataset_ld.sh

# --- Parse locus (task IDs 1-76 map directly to rows 1-76, no header) ---
# loci_input_new.txt columns: locus_no | CHR | left_bp | right_bp
locusinfo=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR == n' "$LOCIFILE" | tr -d '\r')
if [ -z "$locusinfo" ]; then
    echo "ERROR: no locus found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi
locusarray=($locusinfo)
locus=$(printf "%03d" "${locusarray[0]}")
chr_num="${locusarray[1]}"
start="${locusarray[2]}"
end="${locusarray[3]}"
chr="chr${chr_num}"
echo "Locus ${locus}: ${chr}:${start}-${end}"

# --- Skip if already done ---
if [ -f "$OUTDIR/${locus}.ld.gz" ]; then
    echo "Locus ${locus} already exists in $OUTDIR, skipping."
    exit 0
fi

# --- Extend narrow loci (<10 kb) to reduce edge effects ---
locuslength=$(( end - start ))
if [ "$locuslength" -le 10000 ]; then
    start=$(( start - 5000 ))
    end=$(( end + 5000 ))
    echo "  Narrow locus extended to ${chr}:${start}-${end}"
fi

# --- Look up which dosage chunk(s) cover this locus ---
chunkinfo=$(awk -v loc="$locus" '$1 == loc' "$CHUNKLOOKUP")
if [ -z "$chunkinfo" ]; then
    echo "ERROR: locus ${locus} not found in ${CHUNKLOOKUP}."
    echo "  Check that chunk_lookup.txt has an entry for every locus."
    exit 1
fi
chunk1=$(echo "$chunkinfo" | awk '{print $2}')
chunk2=$(echo "$chunkinfo" | awk '{print $3}')
echo "  Chunk(s): ${chunk1}  ${chunk2}"

# --- Build SNP reference file from sumstats ---
# Daner-style format: col 1 = CHR (numeric), col 2 = SNP, col 3 = BP,
#   col 4 = A1 (effect allele), col 5 = A2 (non-effect/reference allele)
# Output: SNP_ID and A2 sorted by position — used by PLINK2 --ref-allele
# to orient all cohorts relative to the GWAS effect allele before LD is computed.
if [ ! -f "$SUMSTATS" ]; then
    echo "ERROR: sumstats file not found: ${SUMSTATS}"
    exit 1
fi
zcat "$SUMSTATS" \
    | awk -v CHR="$chr_num" -v S="$start" -v E="$end" \
        'NR > 1 && $1 == CHR && $3 >= S && $3 <= E { print $2, $5, $3 }' \
    | sort -nk3 \
    | awk '{ print $1, $2 }' \
    > "${locus}.ref"

if [ ! -s "${locus}.ref" ]; then
    echo "ERROR: no SNPs extracted from sumstats for locus ${locus} (${chr}:${start}-${end})"
    exit 1
fi

# A1 file for PLINK1.9 --a1-allele: col 4 = A1 (effect allele), ensures correct LD sign
zcat "$SUMSTATS" \
    | awk -v CHR="$chr_num" -v S="$start" -v E="$end" \
        'NR > 1 && $1 == CHR && $3 >= S && $3 <= E { print $2, $4, $3 }' \
    | sort -nk3 \
    | awk '{ print $1, $2 }' \
    > "${locus}.a1"

# One-column SNP list for PLINK2 --extract
awk '{ print $1 }' "${locus}.ref" > "${locus}.snplist"
echo "  SNPs in reference: $(wc -l < "${locus}.ref")"

# --- Per-cohort LD matrices ---
while IFS= read -r cohort; do
    echo "--- cohort: ${cohort} ---"
    bash make_dataset_ld.sh \
        "$cohort" "$locus" "$chr" "$chr_num" "$start" "$end" \
        "$PLINK2" "$PLINK1" "$chunk1" "$chunk2" "$DATASET_CONFIG"
    rc=$?
    [ "$rc" -ne 0 ] && echo "  WARNING: non-zero exit ($rc) for ${cohort}"
done < "$COHORTLIST"

# --- Save per-cohort LD files for inspection ---
for f in *_"${locus}".ld.gz; do
    [ -f "$f" ] || continue
    cohort_tag="${f%_${locus}.ld.gz}"
    cp "$f" "$OUTDIR/locus_${locus}_${cohort_tag}.ld.gz"
    echo "  Saved per-cohort LD: $OUTDIR/locus_${locus}_${cohort_tag}.ld.gz"
done

# --- Require at least 2 matrices (LDmerge_v2.R enforces this too) ---
nfiles=$(ls *_"${locus}".ld.gz 2>/dev/null | wc -l)
echo "${nfiles} per-cohort LD matrices found for locus ${locus}"
if [ "$nfiles" -lt 2 ]; then
    echo "ERROR: need ≥2 LD matrices to merge; found ${nfiles}. Exiting."
    exit 1
fi

# --- Merge with N-effective weighting ---
Rscript --vanilla LDmerge_v2.R "$locus" PLINK
if [ $? -ne 0 ]; then
    echo "ERROR: LDmerge_v2.R failed for locus ${locus}"
    exit 1
fi

# --- Compress and archive ---
gzip -9 "${locus}.ld"
cp "${locus}.ld.gz"        "$OUTDIR/"
cp "${locus}.snp.log"      "$OUTDIR/"
cp "${locus}.samples.log"  "$OUTDIR/"
echo "Locus ${locus} complete."
