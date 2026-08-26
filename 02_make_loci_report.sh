#!/bin/bash
set -euo pipefail

FINEMAP_ROOT="${FINEMAP_ROOT:?FINEMAP_ROOT must be set -- e.g. export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping}"

SUMSTATS_PREFIX=""
POSITIONAL=()
for arg in "$@"; do
  case "$arg" in
    --sumstats-prefix=*) SUMSTATS_PREFIX="${arg#--sumstats-prefix=}" ;;
    *) POSITIONAL+=("$arg") ;;
  esac
done
set -- "${POSITIONAL[@]}"

PHENO=${1:?"usage: 02_make_loci_report.sh <phenotype> --sumstats-prefix=<prefix>"}
SUMSTATS_PREFIX="${SUMSTATS_PREFIX:?--sumstats-prefix=<prefix> is required -- the per-chromosome sumstats filename prefix}"

BASE=$FINEMAP_ROOT/$PHENO
REPORT=$BASE/clump/${SUMSTATS_PREFIX}_loci_report.txt

echo -e "CHR\tSNP\tP\tSTART\tSTOP" > "$REPORT"

for CHR in $(seq 1 22); do
  SUMSTATS=$BASE/${SUMSTATS_PREFIX}_chr${CHR}.txt
  CLUMP=$BASE/clump/chr${CHR}.clumps

  if [ ! -f "$CLUMP" ]; then
    # plink2 doesn't write a .clumps file at all when nothing passes
    # --clump-p1 (it prints "Warning: No significant --clump results.
    # Skipping." instead of writing a header-only file). That's 0 loci,
    # not a failure -- confirm via the .log, and only warn if there's no
    # such confirmation (which would mean the job never got this far).
    if grep -q "No significant --clump results" "$BASE/clump/chr${CHR}.log" 2>/dev/null; then
      echo "chr$CHR: 0 genome-wide-significant loci, nothing to add to report"
    else
      echo "WARNING: $CLUMP missing and chr${CHR}.log doesn't confirm a clean 'no significant results' run -- check that chromosome's job log" >&2
    fi
    continue
  fi

  awk -v OFS='\t' '
    NR==FNR {
      if (FNR>1) bp[$2]=$3
      next
    }
    FNR==1 { next }
    NF==0 { next }
    {
      chr=$1; pos=$2; id=$3; p=$4; sp2=$11
      start=pos; stop=pos
      if (sp2 != "NONE" && sp2 != "." && sp2 != "") {
        n=split(sp2, ids, ",")
        for (i=1; i<=n; i++) {
          other=ids[i]
          if (other == ".") continue
          if (other in bp) {
            if (bp[other] < start) start=bp[other]
            if (bp[other] > stop) stop=bp[other]
          } else {
            print "WARNING: BP not found for SP2 member " other " (chr " chr ")" > "/dev/stderr"
          }
        }
      }
      print chr, id, p, start, stop
    }
  ' "$SUMSTATS" "$CLUMP" >> "$REPORT"
done

echo "Wrote $REPORT"
wc -l "$REPORT"
