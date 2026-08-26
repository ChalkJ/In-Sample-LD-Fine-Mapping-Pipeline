#!/bin/bash
# Per-cohort helper: dosage -> pgen -> BED -> pairwise LD (.ld.gz)
# Called by run_ld_pipeline.sh for each cohort at each locus.

cohort="$1"
locus="$2"
chr="$3"       # e.g. chr4
chr_num="$4"   # e.g. 4
start="$5"
end="$6"
PLINK2="$7"
PLINK1="$8"
chunk1="${9}"
chunk2="${10}"
DATASET_CONFIG="${11}"

# Sourced for dataset_find_dosage_gz() -- see datasets/ricopili_cross_bcs.sh
# for the required interface. Passed down (as an absolute path) by
# run_ld_pipeline.sh rather than hardcoded here, so a different
# individual-level dataset's naming convention doesn't need this file edited.
source "$DATASET_CONFIG"

# ---------------------------------------------------------------------------
# find_dosage_gz COHORT CHUNK
#   Resolves the wildcard QC-version string in the directory path and returns
#   the full path to the .dosage.gz file, or empty string if not found.
# ---------------------------------------------------------------------------
find_dosage_gz() {
    local c="$1" chunk="$2"
    dataset_find_dosage_gz "$c" "$chunk"
}

# ---------------------------------------------------------------------------
# import_chunk COHORT CHUNK OUT_PREFIX
#   Imports one dosage chunk as a PLINK2 pgen dataset.
#   Returns 0 on success, 1 on failure.
# ---------------------------------------------------------------------------
import_chunk() {
    local c="$1" chunk="$2" out="$3"

    local gz
    gz=$(find_dosage_gz "$c" "$chunk")
    if [ -z "$gz" ]; then
        echo "  WARNING: dosage file not found for cohort=${c} chunk=${chunk}"
        return 1
    fi

    # The .fam and .map share the same base path as the .gz
    local base="${gz%.dosage.gz}.dosage"
    if [ ! -f "${base}.fam" ] || [ ! -f "${base}.map" ]; then
        echo "  WARNING: .fam or .map missing alongside ${gz}"
        return 1
    fi

    $PLINK2 \
        --import-dosage "$gz" \
        --fam "${base}.fam" \
        --map "${base}.map" \
        --memory 4000 --threads 1 \
        --make-pgen \
        --out "$out"

    if [ ! -f "${out}.pgen" ]; then
        echo "  WARNING: PLINK2 import failed for ${c} chunk ${chunk}"
        return 1
    fi
    return 0
}

# ---------------------------------------------------------------------------
# 1. Import chunk 1
# ---------------------------------------------------------------------------
import_chunk "$cohort" "$chunk1" "${cohort}_${locus}_c1" || exit 0

# ---------------------------------------------------------------------------
# 2. If locus spans two chunks, import chunk 2 and merge
# ---------------------------------------------------------------------------
if [ "$chunk2" != "NA" ]; then
    import_chunk "$cohort" "$chunk2" "${cohort}_${locus}_c2"
    if [ -f "${cohort}_${locus}_c2.pgen" ]; then
        printf '%s\n' "${cohort}_${locus}_c1" "${cohort}_${locus}_c2" \
            > "${cohort}_${locus}_mergelist.txt"

        $PLINK2 \
            --pmerge-list "${cohort}_${locus}_mergelist.txt" pfile \
            --memory 4000 --threads 1 \
            --make-pgen \
            --out "${cohort}_${locus}_merged"

        if [ -f "${cohort}_${locus}_merged.pgen" ]; then
            mv "${cohort}_${locus}_merged.pgen" "${cohort}_${locus}_c1.pgen"
            mv "${cohort}_${locus}_merged.pvar" "${cohort}_${locus}_c1.pvar"
            mv "${cohort}_${locus}_merged.psam" "${cohort}_${locus}_c1.psam"
            echo "  Two chunks merged for ${cohort}"
        else
            echo "  WARNING: chunk merge failed for ${cohort}; using chunk1 only"
        fi

        rm -f "${cohort}_${locus}_c2".* \
               "${cohort}_${locus}_mergelist.txt" \
               "${cohort}_${locus}_merged".*
    fi
fi

# ---------------------------------------------------------------------------
# 3. Extract locus region, restrict to GWAS SNPs, remove duplicate IDs
#    --rm-dup force-first: for any variant seen more than once (imputation
#    window overlap between adjacent chunks), keep the first occurrence.
#
#    --ref-allele: sets the A2 (reference/non-effect) allele from the GWAS
#    .ref file so that all cohorts have consistent allele orientation before
#    PLINK1.9 computes R. This ensures the sign of R is relative to the GWAS
#    effect allele, which SuSiE requires (R must be consistent with z-scores).
# ---------------------------------------------------------------------------
$PLINK2 \
    --pfile "${cohort}_${locus}_c1" \
    --extract "${locus}.snplist" \
    --ref-allele "${locus}.ref" \
    --rm-dup force-first \
    --memory 4000 --threads 1 \
    --make-bed \
    --out "${cohort}_${locus}"

rm -f "${cohort}_${locus}_c1".*

if [ ! -f "${cohort}_${locus}.bed" ]; then
    echo "  WARNING: no variants remain for ${cohort} at locus ${locus} after filtering"
    exit 0
fi

nvars=$(wc -l < "${cohort}_${locus}.bim")
if [ "$nvars" -lt 2 ]; then
    echo "  WARNING: only ${nvars} variant(s) for ${cohort} at locus ${locus} — skipping"
    rm -f "${cohort}_${locus}".{bed,bim,fam}
    exit 0
fi
echo "  ${cohort}: ${nvars} variants in region"

# ---------------------------------------------------------------------------
# 4. Compute pairwise signed LD (R) with PLINK 1.9
#    --ld-window / --ld-window-kb: large values so all pairs are considered.
#    --ld-window-r2 0: output pairs with any R² (not just ≥ 0.2 default).
#    Output: ${cohort}_${locus}.ld  (columns: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R)
#    This matches the "PLINK" format expected by LDmerge_v2.R.
# ---------------------------------------------------------------------------
$PLINK1 \
    --bfile "${cohort}_${locus}" \
    --a1-allele "${locus}.a1" 2 1 \
    --r \
    --ld-window 999999 \
    --ld-window-kb 99999 \
    --ld-window-r2 0 \
    --memory 4000 --threads 1 \
    --out "${cohort}_${locus}"

if [ ! -f "${cohort}_${locus}.ld" ]; then
    echo "  WARNING: PLINK1.9 --r failed for ${cohort} at locus ${locus}"
    rm -f "${cohort}_${locus}".{bed,bim,fam}
    exit 0
fi

# ---------------------------------------------------------------------------
# 5. Compress LD file; keep .fam for LDmerge_v2.R N-effective calculation
# ---------------------------------------------------------------------------
gzip -9 "${cohort}_${locus}.ld"
rm -f "${cohort}_${locus}".{bed,bim,log}

npairs=$(zcat "${cohort}_${locus}.ld.gz" | wc -l)
echo "  ${cohort}: LD done (${npairs} pairs)"
