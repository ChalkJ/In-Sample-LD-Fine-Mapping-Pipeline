#!/bin/bash
set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --dataset-config=<path> is pulled out of the args wherever it appears;
# the remaining positional args are PHENO and (optionally) REF_COHORT, e.g.:
#   03_make_chunk_lookup.sh sczvscon
#   03_make_chunk_lookup.sh sczvscon grp10neu3 --dataset-config=datasets/my_other_dataset.sh
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
source "$DATASET_CONFIG"

PHENO=${1:?"usage: 03_make_chunk_lookup.sh <sczvscon|bipvscon> [ref_cohort] [--dataset-config=<path>]"}
REF_COHORT=${2:-grp10neu3}

case "$PHENO" in
  sczvscon) PREFIX=sc_vs_allcontrols ;;
  bipvscon) PREFIX=bp_vs_allcontrols ;;
  *) echo "unknown PHENO: $PHENO" >&2; exit 1 ;;
esac

BASE=$FINEMAP_ROOT/$PHENO
LOCI_REPORT=$BASE/clump/${PREFIX}_loci_report.txt
CHUNK_OUT=$BASE/chunk_lookup.txt
LOCI_OUT=$BASE/loci_input.txt

if [ ! -s "$LOCI_REPORT" ]; then
  echo "ERROR: $LOCI_REPORT missing/empty -- run 02_make_loci_report.sh first" >&2
  exit 1
fi

# Skip-if-done: loci_input.txt/chunk_lookup.txt are hand-edited after this
# script runs (MHC-region loci dropped, chunk-boundary loci forced into
# two-chunk entries -- see README "Known gotchas"). Regenerating them
# unconditionally on every rerun would silently overwrite those edits and
# renumber every locus after the first dropped one, desyncing existing
# output/ld/NNN.* and fine-mapping files from their real locus definitions.
# Delete BOTH files to force a genuine rerun (e.g. after the clump report
# actually changed).
if [ -s "$LOCI_OUT" ] && [ -s "$CHUNK_OUT" ]; then
  echo "$LOCI_OUT and $CHUNK_OUT already exist, skipping (preserves any hand edits, e.g. MHC exclusion)."
  echo "Delete both files to force a rerun."
  exit 0
fi

# Chunk boundaries are the same dosage chunking for every cohort on this
# build, so we only need to read one cohort's directory to build the
# chr -> chunk(start_mb, end_mb) table, rather than hardcoding the file list.
REF_QC1_DIR=$(dataset_find_qc1_dir "$REF_COHORT")
if [ -z "$REF_QC1_DIR" ]; then
  echo "ERROR: dataset_find_qc1_dir found nothing for cohort $REF_COHORT (config: $DATASET_CONFIG)" >&2
  exit 1
fi

CHUNK_TABLE=$(mktemp)
ls "$REF_QC1_DIR"/*.out.dosage.gz \
  | sed -E "s#${DATASET_CHUNK_NAME_REGEX}#\1#" \
  | sort -u \
  | awk -F'_' -v OFS='\t' '{ chrnum=substr($1,4); print chrnum, $0, $2+0, $3+0 }' \
  > "$CHUNK_TABLE"

# CHR, START, STOP from the loci report, sorted and numbered.
# Two locus-number forms are carried through:
#   plain  (1, 2, ... 10, ...)   -> loci_input.txt (run_ld_pipeline.sh does its
#                                    own printf "%03d" padding on this; if this
#                                    were already zero-padded, e.g. "008", bash's
#                                    printf %d would misparse it as octal and error)
#   padded (001, 002, ...)       -> chunk_lookup.txt, to match the padded string
#                                    run_ld_pipeline.sh looks up with
LOCI_NUMBERED=$(mktemp)
tail -n +2 "$LOCI_REPORT" \
  | sort -k1,1n -k4,4n \
  | awk -v OFS='\t' '{ printf "%d\t%03d\t%s\t%s\t%s\n", NR, NR, $1, $4, $5 }' \
  > "$LOCI_NUMBERED"

awk -v OFS='\t' '{ print $1, $3, $4, $5 }' "$LOCI_NUMBERED" > "$LOCI_OUT"

awk -v OFS=' ' '
  NR==FNR {
    n[$1]++
    i=n[$1]
    chr_s[$1,i]=$3
    chr_e[$1,i]=$4
    chr_lab[$1,i]=$2
    next
  }
  {
    locus_padded=$2; chr=$3; start=$4; stop=$5
    smb=int(start/1000000)
    emb=int(stop/1000000)
    sidx=find_chunk(chr, smb, locus_padded)
    eidx=find_chunk(chr, emb, locus_padded)
    if (chr_lab[chr,sidx] == chr_lab[chr,eidx]) {
      print locus_padded, chr_lab[chr,sidx], "NA"
    } else {
      print locus_padded, chr_lab[chr,sidx], chr_lab[chr,eidx]
    }
  }
  function find_chunk(chr, mb, locus,   i, s, e, d, best, bestd) {
    best=1; bestd=999999
    for (i=1; i<=n[chr]; i++) {
      s=chr_s[chr,i]; e=chr_e[chr,i]
      if (mb>=s && mb<=e) return i
      d = (mb<s) ? s-mb : mb-e
      if (d<bestd) { bestd=d; best=i }
    }
    print "WARNING: locus " locus " chr" chr " pos " mb "Mb falls between chunks, using nearest chunk " chr_lab[chr,best] > "/dev/stderr"
    return best
  }
' "$CHUNK_TABLE" "$LOCI_NUMBERED" > "$CHUNK_OUT"

rm -f "$CHUNK_TABLE" "$LOCI_NUMBERED"

echo "Wrote $LOCI_OUT"
wc -l "$LOCI_OUT"
echo "Wrote $CHUNK_OUT"
wc -l "$CHUNK_OUT"
