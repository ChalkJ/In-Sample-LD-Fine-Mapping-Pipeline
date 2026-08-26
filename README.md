# In-Sample LD & Fine-Mapping Pipeline

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)

*A self-contained pipeline, performing multi-method fine-mapping and
generating in-sample LD matrices from individual-level data.*

A self-contained SLURM/HPC pipeline: takes per-chromosome GWAS summary
statistics through sumstats QC, clumping, locus definition, in-sample 
linkage-disequilibrium (LD) generation, and config-driven fine-mapping 
(SuSiE/FINEMAP/both/none — all native to the cluster, no local machine or
WSL involved), finishing with a results / reporting stage (per-SNP table,
summary statistics, locus-zoom plots). See "Quick start" below.

VEP/gene annotation is **not** part of this pipeline due to mirror instability,
it's a manual seperate task (submitting to the VEP web tool) done at a later phase.
The locus-zoom plots this pipeline does produce are plain PIP-colored SNP scatters 
with no gene track. The gene track can be added once positional values have been acquired
through a VEP call.

Originally built for the PGC SCZvsBP case-case GWAS fine-mapping
project, but the pipeline itself is phenotype- and dataset-agnostic.
See "Configuration" and "Dataset config" below.

## Quick start

```bash
export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping   # see "Configuration" below
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas
```

That's the whole thing. It queues every stage below as SLURM jobs chained
with `--dependency=afterok`, so they run in order automatically.
Monitor with:

```bash
squeue -u $USER
```

`<phenotype>` (`my_pheno` above) can be anything you like, it's just the
subdirectory name under `$FINEMAP_ROOT`. `--sumstats-prefix` is required:
it's whatever your own per-chromosome sumstats files are named (see
"Prerequisites" below). There is no fixed list of accepted phenotype names —
any identifier works, nothing in the pipeline needs editing to add a new one.

Optional flags (all can be combined):

```bash
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas --dedup=off
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas --hetpva-cutoff=0.01
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas --hetpva-cutoff=off --dedup=off
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas --ref-cohort=grp10neu3
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas --dataset-config=datasets/my_other_dataset.sh
bash submit_full_pipeline.sh my_pheno --sumstats-prefix=my_pheno_gwas --finemap-config=finemap_config_susie_only.R
```

## Configuration

**`FINEMAP_ROOT`** (required environment variable) the working-directory
root every stage reads and writes under (per-phenotype subdirectories,
sumstats, LD output, fine-mapping output, results). Every script checks
this is set and fails immediately with a (hopefully) clear message if it isn't:

```bash
export FINEMAP_ROOT=/gpfs/home3/<you>/finemapping
```

Set once in your shell before running `submit_full_pipeline.sh` (or in
your shell profile) — `sbatch` inherits the submitting shell's environment
by default, this export propagates automatically through the whole
self-chaining job tree the pipeline submits. It also determines where the
plink2/plink1.9 binaries are expected, this assumes pre-placed binaries, not a HPC module
(`$FINEMAP_ROOT/plink/plink2/plink2`, `$FINEMAP_ROOT/plink/plink1.9/plink`).

Two more, separate config layers, covered in their own sections below:
**`--dataset-config=<path>`** (which individual-level genotype dataset to
use) and **`--finemap-config=<path>`** (fine-mapping method/parameters).

## Prerequisites

Before running, for the phenotype you're processing (`$PHENO`), under
`$FINEMAP_ROOT/$PHENO/`:

- 22 per-chromosome sumstats files, one per autosome, named
  `<prefix>_chr<N>.txt` where `<prefix>` is whatever you pass as
  `--sumstats-prefix=<prefix>` entirely up to you, however your own
  upstream GWAS pipeline happened to name its per-chromosome output.
- `EUR_cohorts.txt` — one cohort name per line, e.g.:
  ```
  grp10neu3
  grp11uk2
  grp12ger3
  ...
  ```
  These names must match directories under the individual-level dataset's
  root (see "Dataset config" below) — one cohort's genotype/dosage data per
  line.

Plink binaries, expected under `$FINEMAP_ROOT` (derived from `FINEMAP_ROOT`
in `01_extract_merge_clump.sh` and `run_ld_pipeline.sh` — not
dataset-config-driven, since these are software installs, not data):
`$FINEMAP_ROOT/plink/plink2/plink2` and `$FINEMAP_ROOT/plink/plink1.9/plink`
(note: `plink`, not `plink1.9`, is the actual binary name in that directory).

## Input format

Each per-chromosome sumstats file is expected to be in DANER format ie:

```
CHR  SNP  BP  A1  A2  FRQ_A_x  FRQ_U_x  INFO  OR  SE  P  ngt  Direction  HetISqt  HetDf  HetPVa  Nca  Nco  Neff_half
```

`CHR`, `SNP`, `BP`, `A1`, `A2` in columns 1–5 specifically must stay in that
order — `run_ld_pipeline.sh` and `make_dataset_ld.sh` read these by column
position when building the SNP reference/allele-orientation files.

## Stages

Each stage is idempotent/skip-if-done (safe and cheap to resubmit the whole
pipeline after a partial failure) — the exact skip check for each is listed
below so a rerun's behavior is predictable rather than something to
reverse-engineer from the script.

| # | Script | What it does | Skip-if-done check |
|---|---|---|---|
| QC | `qc_filter_sumstats.sh` → `qc_filter_sumstats.R` | Dedups per-chromosome sumstats (default: keep lowest P per CHR:BP) and filters HetPVa (default: drop < 0.05; NA/non-numeric always kept), then combines all 22 chromosomes into one gzipped daner-format file. Both filters are optional/tunable — `--dedup=off`, `--hetpva-cutoff=<value>\|off`. Writes a plain-text QC report alongside (row counts in/out per chromosome, per filter) so filtering is auditable, not a black box. | Skips entirely if the combined output `.gz` already exists and is non-empty. Delete it to force a rerun (e.g. after changing the filter settings). |
| 00 | `00_make_snplists.sh` | Extracts the SNP column from each per-chromosome sumstats file, for use as a plink `--extract` list in stage 01. | None — always regenerates (cheap, deterministic). |
| 01 | `01_extract_merge_clump.sh` | SLURM array (chromosomes 1–22). Per cohort: plink2 extract of that chromosome's GWAS SNPs from the cohort's best-guess genotypes; then plink1.9 `--merge-list` across cohorts (not plink2 `--pmerge-list` — plink2 doesn't support merging different cohorts' samples over the same SNP set); then plink2 `--clump` (p1=5e-8, p2=1e-4, r2=0.1, kb=3000). | Per-chromosome: skips if `clump/chr<N>.clumps` exists and is non-empty. Delete that file to force a rerun of that one chromosome. |
| 02 | `02_make_loci_report.sh` | Concatenates the 22 chromosomes' clump output into one `CHR SNP P START STOP` report, deriving each locus's boundaries from the min/max BP across the index SNP and all its secondary-clumped (`SP2`) SNPs. | None — always regenerates from the existing per-chromosome clump files (cheap). A missing `.clumps` file for a chromosome is treated as "0 genome-wide-significant loci there" (confirmed via that chromosome's `.log`), not an error. |
| 03 | `03_make_chunk_lookup.sh` | Derives the dosage-chunk boundary table from one reference cohort's actual files on disk (via the dataset config, not hardcoded), then writes `loci_input.txt` (plain locus numbers, unpadded — see "Known gotchas") and `chunk_lookup.txt` (zero-padded locus numbers, one or two chunk IDs). | Skips entirely if **both** `loci_input.txt` and `chunk_lookup.txt` already exist and are non-empty — deliberately, since these files are hand-edited afterward (MHC exclusion, chunk-boundary fixes) and regenerating them unconditionally would silently overwrite those edits and renumber every locus after the first dropped one. Delete both files to force a rerun. |
| 04 | `04_finalize_loci_and_launch_ld.sh` | Runs 02 then 03, counts the resulting loci, then self-submits `run_ld_pipeline.sh` as a `--array=1-N` job sized to that exact count — the step that used to mean watching `squeue`, counting loci by hand, and submitting stage 5 yourself. | None itself (cheap to rerun — 02/03 regenerate deterministically and 05's own per-locus skip check makes a duplicate array submission harmless, just wasteful). |
| LD | `run_ld_pipeline.sh` + `make_dataset_ld.sh` + `LDmerge_v2.R` | SLURM array, one task per locus. Imports the relevant dosage chunk(s), extracts that locus's GWAS SNPs with consistent allele orientation, computes per-cohort pairwise LD (plink1.9 `--r`), then merges across cohorts with N-effective weighting (R). | Per-locus: skips if `output/ld/<locus>.ld.gz` already exists. Delete that file to force a rerun of that one locus. |
| 05 | `05_run_finemapping.sh` → `run_finemapping.R` | Single SLURM job (not an array — see "Fine-mapping" below for why), self-submitted by 04 once the LD array completes. Config-driven SuSiE / FINEMAP / both / none (`finemap_config.R`), with an adaptive re-run pass for loci that saturate the initial `max_causal_variants`. | Per-locus, per enabled method: skips a locus if its SuSiE `.rds` / FINEMAP `.snp` output already exists (whichever method(s) are enabled). `method="none"` skips the whole stage immediately. |
| 06 | `06_generate_reports.sh` → `generate_reports.R` | Single SLURM job, self-submitted by 05 right after its `Rscript` call succeeds (no `--dependency` needed for that submission — see "Reports" below). Existence-driven: reads whatever SuSiE/FINEMAP output is actually on disk and writes a per-SNP results table, a genome-wide summary-stats rollup, and per-locus plots. No VEP/gene annotation. | Table/summary TSVs always regenerate (cheap, deterministic). Each plot PNG skips if it already exists — a partial prior run never redraws finished plots. |

## Output structure

Under `$FINEMAP_ROOT/$PHENO/`:

```
qc_sumstats_report.txt                       # QC stage: what was dropped and why
daner_<phenotype>_qc.gz                      # QC stage: combined, filtered sumstats (fixed naming convention)
snplists/<prefix>_chr<N>.snplist             # stage 00
clump/chr<N>.clumps, chr<N>.log              # stage 01
clump/<prefix>_loci_report.txt               # stage 02
loci_input.txt                               # stage 03: locus_id  chr  start  stop  (unpadded, no header)
chunk_lookup.txt                             # stage 03: locus_id(padded)  chunk1  chunk2-or-NA
output/ld/<locus>.ld.gz                      # LD stage: merged LD matrix (gzipped)
output/ld/<locus>.snp.log                    # LD stage: SNPs entering that locus's LD
output/ld/<locus>.samples.log                # LD stage: sample info for N-effective weighting
output/ld/locus_<locus>_<cohort>.ld.gz       # LD stage: per-cohort LD (kept for inspection)
output/LD_dosage_<jobid>_<taskid>.{out,err}  # LD stage: SLURM array logs
susie/<locus>.susie_l<L>.rds                 # stage 05: SuSiE fit object (L = max_causal_variants used)
finemap/<locus>/<locus>.{z,ld,snp,config,cred,log}  # stage 05: FINEMAP inputs/outputs
finemap/master, finemap/master_rerun_k<N>    # stage 05: FINEMAP multi-locus --in-files batch lists
finemapping_<phenotype>_summary.tsv          # stage 05: per-locus SuSiE/FINEMAP diagnostics summary
results/all_snps_<phenotype>.tsv             # stage 06: one row per SNP -- sumstats + PIPs + locus + CS
results/summary_stats_<phenotype>.tsv        # stage 06: genome-wide rollup (mean CS size, mean PIP, ...)
results/plots/<locus>_susie_scatter.png      # stage 06: PIP-colored SNP scatter (SuSiE)
results/plots/<locus>_finemap_scatter.png    # stage 06: same, FINEMAP
results/plots/<locus>_ld_heatmap.png         # stage 06: LD r heatmap
```

Example real content (from an existing analysis, format reference only):

```
# chunk_lookup.txt
001 chr1_172_196 NA
002 chr8_137_146 NA

# loci_input.txt
1	1	173461333	175049333
2	8	143902410	144039410
```

## Dataset config (portability to a different individual-level dataset)

The pipeline's individual-level genotype/dosage data (paths, cohort file
naming, dosage-chunk naming) is described by one sourced config file, not
hardcoded into the scripts — see `datasets/ricopili_cross_bcs.sh`, the
default (this project's actual RICOPILI `cross_bcs` layout) and template.
Every stage that touches individual-level data (`01_extract_merge_clump.sh`,
`03_make_chunk_lookup.sh`, `run_ld_pipeline.sh` → `make_dataset_ld.sh`)
accepts `--dataset-config=<path>`, defaulting to that file if omitted.

To point this pipeline at a different dataset: copy
`datasets/ricopili_cross_bcs.sh` to a new file and rewrite its variable and
three functions to match the new dataset's actual layout — the required
interface (what each function must be given and must echo) is documented in
that file's header comment. Nothing else needs to change; different chunk
*boundaries* on a different dataset are already handled automatically (stage
03 always re-derives them from whatever's actually on disk, never
hardcoded).

## Fine-mapping (stage 05)

Config-driven via `finemap_config.R` (sourced, not parsed — copy it and pass
`--finemap-config=<path>` to use a different one instead of editing the
shipped default). Key fields, all with defaults matching every existing
published result on this project:

- `method`: `"susie"` | `"finemap"` | `"both"` | `"none"`. `"both"` runs
  both and cross-checks their top SNP per locus (today's default).
  `"none"` skips fine-mapping entirely — the stage is still submitted by
  `04` every time (simpler than teaching bash to parse this R config file),
  it just exits immediately.
- `max_causal_variants`: SuSiE's `L` and FINEMAP's `--n-causal-snps`/`k` for
  the initial pass. Loci that saturate this (SuSiE doesn't converge with
  this many credible sets, or FINEMAP's posterior probability for
  `k = max_causal_variants` exceeds `adaptive_rerun$k_posterior_threshold`)
  are automatically re-run once at `adaptive_rerun$max_causal_variants_rerun`,
  then (if still unconverged) again with `adaptive_rerun$rerun_max_iter`
  SuSiE iterations.
- `finemap$binary`: **blank by default and must be set if `method` is
  `"finemap"` or `"both"`** — Snellius has no FINEMAP module, so you need to
  place a FINEMAP binary in your own scratch/home space first and point
  this at its full path. `run_finemapping.R` checks this up front and fails
  loudly, before doing any SuSiE work, if it's needed but blank or not
  executable — you shouldn't discover this partway through a long run.
- `susie$estimate_residual_variance`, `susie$lambda`: passed straight
  through to `susieR::susie_rss()`.

**Not yet verified**: whether `susieR` is installed in Snellius's
`R/4.5.1-gfbf-2025a` module. Check this (`Rscript -e 'library(susieR)'`
after `module load R/4.5.1-gfbf-2025a`) before a first real run with
`method` including `"susie"`.

Runs as **one single SLURM job, not a per-locus array** — it loops over
every locus sequentially, and FINEMAP's `--in-files` batches every
not-yet-complete locus into one invocation rather than one call per locus
(this is FINEMAP's own efficient multi-locus mode, not something this
pipeline could easily split across array tasks without restructuring both
that batching and the multi-locus adaptive-rerun logic). This means it can
run long — the default `--time=24:00:00` in `05_run_finemapping.sh` may
need to be raised for a phenotype with many loci; override the same way as
any other stage's SBATCH directives (see that script's header comment).

## Reports (stage 06)

`generate_reports.R` — no config, no flags beyond the phenotype. It's
**existence-driven**: it looks at whatever SuSiE `.rds` / FINEMAP `.snp`
files actually exist under `susie/` and `finemap/` and reports on those,
rather than trusting what `finemap_config.R`'s `method=` currently says —
so re-running this stage after the config changed (or against a partial
fine-mapping run) still reflects reality.

- **`results/all_snps_<phenotype>.tsv`**: one row per SNP that entered
  fine-mapping at each locus — `locus, CHR, SNP, BP, A1, A2, OR, SE, P,
  FRQ_A, FRQ_U, zscore, susie_pip, susie_cs, finemap_pip, finemap_cs`. The
  `*_pip`/`*_cs` columns are `NA` wherever that method has no output for
  that locus (never fabricated). `*_cs` is the credible-set index (`1`,
  `2`, …) that SNP belongs to, or `NA` if it's not in any CS.
- **`results/summary_stats_<phenotype>.tsv`**: genome-wide rollup —
  `n_loci`, `n_loci_with_susie`/`n_loci_with_finemap`, mean/median max-PIP
  per method (from stage 05's per-locus summary), and mean CS
  size/mean number of CS per locus/count of zero-CS loci per method
  (computed directly from `all_snps`'s actual CS membership, not from the
  per-locus summary's comma-joined size string — this also gives FINEMAP a
  real CS-size rollup, which the per-locus summary alone doesn't have).
  "Top SNP by PIP per locus" isn't repeated here — it already lives in
  `finemapping_<phenotype>_summary.tsv`'s `s_top_snp`/`f_top_snp` columns.
- **Plots**: `<locus>_susie_scatter.png`/`<locus>_finemap_scatter.png` —
  PIP-colored SNP scatter (color/size/alpha by PIP, credible-set members
  highlighted, everything else grey), one per locus per method that has
  output there. `<locus>_ld_heatmap.png` — LD r heatmap, reusing the LD
  matrix the LD stage already computed (no new computation). None of these
  include a gene-boundary track (that's a VEP/biomaRt-dependent addition in
  `master_results.Rmd`, deliberately dropped here).
- **A locus that went through fine-mapping's adaptive re-run has two SuSiE
  files on disk** (e.g. both `NNN.susie_l5.rds` and `NNN.susie_l10.rds` —
  the rerun writes a new L-suffixed file rather than overwriting). This
  stage always uses the highest-L file present, since the rerun only ever
  moves L up.
- **FINEMAP's `.credK` files are not "the Nth credible set"** — a locus can
  have `.cred1`, `.cred2`, … simultaneously on disk, each a separate
  hypothesis assuming exactly `K` causal SNPs. This stage picks whichever
  has the highest posterior probability in its own header line, not the
  highest `K` (a real bug already hit and fixed once elsewhere on this
  project, in `summarize_polyfun_tables.R` — ported here correctly from the
  start, not re-derived).

**Not yet verified**: whether `ggplot2`/`scales` are installed in Snellius's
`R/4.5.1-gfbf-2025a` module (same caveat as `susieR` above) — needed for the
scatter plots; the LD heatmap only needs base R.

## Known gotchas

- **`plink2`/`plink1.9` binary paths**: `finemap/plink/plink2/plink2` and
  `finemap/plink/plink1.9/plink` — note the second one is a directory named
  `plink1.9` containing a binary literally named `plink`, not `plink1.9`.
- **plink2 vs plink1.9 for the cohort merge** (stage 01): plink2's
  `--pmerge-list` only supports "concatenating" merges (same samples,
  different variants). This pipeline's merge is the opposite — same
  variants, different samples per cohort — which plink2 explicitly doesn't
  support yet, hence plink1.9's `--merge-list` for that one step.
  `--allow-no-sex` is required on that call (ambiguous-sex samples plus a
  loaded phenotype otherwise blocks combining `--merge-list` with
  `--make-bed`).
- **`loci_input.txt`'s locus numbers must stay unpadded** (plain `1`, `2`,
  … not `001`, `002`). `run_ld_pipeline.sh` does its own `printf "%03d"`
  padding on this column; a pre-padded value like `"008"` would be
  misparsed by bash's `printf %d` as octal and crash.
- **MHC region (chr6:30,000,000–33,000,000)**: historically caused
  multi-hour+ LD computation (hits SLURM time limits) due to extreme LD
  complexity. Not automatically excluded by this pipeline — if you hit this,
  the established practice on this project has been to drop MHC-overlapping
  loci from `loci_input.txt`/`chunk_lookup.txt` by hand after stage 03,
  before stage 05 runs. **Stage 03 now skips entirely if those two files
  already exist** (see the stages table above), specifically so a rerun of
  the pipeline never silently undoes this kind of hand edit.
- **Loci landing exactly on a chunk's rounded-Mb boundary**: dosage chunk
  filenames encode *rounded* Mb boundaries, so a locus whose start/stop
  lands exactly on a chunk's own boundary number is at risk of the wrong or
  incomplete chunk being assigned. Not automatically detected by this
  pipeline — established practice has been to manually force such loci into
  genuine two-chunk loci (add the adjacent chunk instead of `NA` in
  `chunk_lookup.txt`) after stage 03. Same skip-if-done protection as above.
- **`fread` silently drops malformed/truncated sumstats rows** — found for
  real during validation: `sczvscon`'s `sc_vs_allcontrols_chr5.txt` has one
  line (SNP `rs163017`) truncated after the `Direction` field, missing
  `HetISqt` onward. `data.table::fread` treats a short trailing line like
  this as a "footer" and drops it with only an R warning — easy to miss in
  a SLURM `.err` log. `qc_filter_sumstats.R` now cross-checks each
  chromosome's raw line count against what `fread` actually returned and
  writes an explicit `WARNING:` line in the QC report if they differ, so
  this shows up somewhere you'll actually see it.
- **A missing `.clumps` file after stage 01 is not necessarily a failure** —
  plink2 doesn't write one at all when zero variants pass `--clump-p1`
  (prints "No significant --clump results" instead of a header-only file).
  Stage 02 checks the chromosome's `.log` for that message before treating a
  missing `.clumps` as a real error.

## Status

Written and locally verified every shell script bash-syntax-checks, every
R script parse-checks, and the algorithmically trickiest pieces (sumstats
QC filter logic, dataset-config plumbing, fine-mapping's method-toggle and
fail-fast paths, the highest-L SuSiE and best-posterior FINEMAP `.credK`
resolution in the reports stage, plot generation) were verified against
real or realistic synthetic data. **Not yet run end-to-end on a live SLURM
cluster** — confirm `susieR`/`ggplot2`/`scales` are available in your R
environment before a first real run (see "Fine-mapping" and "Reports"
above), and set `finemap_config.R`'s `finemap$binary` or use `method="susie"`
if no FINEMAP binary is available.

## Citation

If you use this pipeline, please cite it :) see [`CITATION.cff`](CITATION.cff)
for the machine-readable citation metadata. If this repository has been archived on Zenodo,
prefer citing the archived version's DOI:

[![DOI](https://zenodo.org/badge/1347179228.svg)](https://doi.org/10.5281/zenodo.22112840)

## License

Licensed under the Apache License, Version 2.0 — see [`LICENSE`](LICENSE)
for the full text.
