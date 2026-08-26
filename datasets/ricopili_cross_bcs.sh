#!/bin/bash
# =============================================================================
# Dataset config: RICOPILI cross_bcs individual-level data (the dataset this
# pipeline was originally built against). Sourced (not executed) by
# 01_extract_merge_clump.sh, 03_make_chunk_lookup.sh, run_ld_pipeline.sh, and
# make_dataset_ld.sh, via --dataset-config=<path>.
#
# This file is also the template for a DIFFERENT individual-level dataset:
# copy it, rename it, and rewrite the variable/functions below to match that
# dataset's actual path and file-naming convention. Every consuming script
# only calls these functions/variables -- none of them hardcode a path or
# filename pattern of their own anymore.
#
# Required interface (must be preserved in any new dataset config file):
#   DATASET_BG_ROOT           -- root directory containing one subdirectory
#                                 per cohort.
#   dataset_find_cohort_bed <cohort>
#                             -- echo the full path (with .bed extension) to
#                                 cohort <cohort>'s best-guess plink1
#                                 bed/bim/fam file. Echo nothing (no error) if
#                                 not found -- callers check for an empty
#                                 result themselves. If more than one file
#                                 could match, echo all of them, one per
#                                 line (some callers warn and pick the first).
#   dataset_find_qc1_dir <cohort>
#                             -- echo the directory holding cohort <cohort>'s
#                                 RICOPILI-style dosage-chunk files (used to
#                                 discover chunk Mb boundaries from disk).
#                                 Echo nothing if not found.
#   DATASET_CHUNK_NAME_REGEX  -- a `sed -E` pattern that extracts the
#                                 "chrN_STARTMB_ENDMB" token (e.g.
#                                 "chr4_172_196") from one dosage filename in
#                                 the directory above, as capture group \1.
#   dataset_find_dosage_gz <cohort> <chunk>
#                             -- echo the full path to cohort <cohort>'s
#                                 dosage .gz file for chunk <chunk> (chunk is
#                                 e.g. "chr4_172_196", matching
#                                 DATASET_CHUNK_NAME_REGEX's token format).
#                                 Echo nothing if not found.
# =============================================================================

DATASET_BG_ROOT="/gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/cross_bcs"

dataset_find_cohort_bed() {
  local cohort="$1"
  ls "${DATASET_BG_ROOT}/${cohort}/cobg_dir_genome_wide"/*.bg.bed 2>/dev/null
}

dataset_find_qc1_dir() {
  local cohort="$1"
  ls -d "${DATASET_BG_ROOT}/${cohort}/dasuqc1_"*/qc1 2>/dev/null
}

DATASET_CHUNK_NAME_REGEX='.*\.(chr[0-9]+_[0-9]{3}_[0-9]{3})\.out\.dosage\.gz$'

dataset_find_dosage_gz() {
  local cohort="$1" chunk="$2"
  ls "${DATASET_BG_ROOT}/${cohort}/dasuqc1_scb_${cohort}_eur_sa-qc"*.hg38.ch.fl/qc1/"dos_scb_${cohort}_eur_sa-qc"*.hg38.ch.fl."${chunk}".out.dosage.gz \
      2>/dev/null | head -1
}
