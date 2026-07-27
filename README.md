# multiPPIseq_DMS_central

Central code + methods home for the **multiPPIseq DMS** pipeline — Python
scripts, SGE/qsub job scripts, shared reference libraries, and the analysis plan
/ lab notebook. Sequencing data lives outside this repo (per-batch
`multiPPIseq_DMS_dict_*` / `multiPPIseq_DMS_editRate_*` folders).

The pipeline has two parts:

1. **Dictionary building** (`*/dict_build/`) — from raw amplicon reads, build a
   barcode → variant dictionary per gene/block, plus build-time diagnostics and
   (optional) combined / concatenated dictionaries.
2. **cDNA edit count** (`*/cdna_edit_count/`) — apply a dictionary to cDNA reads
   and count base edits per variant, with cDNA-sample diagnostics.

## Layout

```
multiPPIseq_DMS_central/
  README.md                       # this file
  analysis_plan/                  # methods specs (start at analysis_plan/README.md)
    dict_build_analysis_plan.md   #   methods + onboarding — dictionary building
    cdna_edit_count_analysis_plan.md  # methods + onboarding — cDNA edit count
  script/
    dict_build/                   # dictionary-build Python (step1_v2/step2/step3 + diagnostics)
    cdna_edit_count/              # cDNA parse Step1/Step2 + run-stats
  qsub/
    dict_build/                   # step-wise SGE job scripts + config.sh
    cdna_edit_count/              # step-wise SGE job scripts + config.sh
  ref/dict_build/                 # shared refs: r2 fixed-seq, block-ref fasta + bowtie2 index, variant lib
```

A working checkout also carries local-only material that is **not tracked** (see
`.gitignore`): the `analysis_plan/*_decisions.md` logs, and a `metadata/`
directory holding the fastq sample inventory. Neither is needed to run the
pipeline — the methods specs above are self-sufficient — and neither is backed
up by this repo.

Each pipeline's `config.sh` holds the batch-specific knobs (paths, genes,
parameters); the step scripts source it so they stay generic. Override any value
by exporting it before `qsub` (e.g. `DATA_PATH=… qsub …`).

## Part 1 — dictionary building (`qsub/dict_build/`)

Edit `qsub/dict_build/config.sh` for the batch (`DATA_PATH`, `LIBTYPE`, `GENES`),
then submit the steps in order.

| Step | Script | What |
|------|--------|------|
| 00 | `00_fastqc.qsub` | FastQC on raw reads → inspect for adapter read-through |
| 01 | `01_cutadapt.qsub` | *(optional)* adapter trim; set `USE_NOADAPTER=1` after |
| 02 | `02_extract_barcodes.qsub` | extract bc1/bc2/oligo → FASTA (`step1_v2.py`) |
| 03 | `03_bowtie2_align.qsub` | align oligo to block refs → BAM |
| 04 | `04_variant_map.qsub` | map reads to variants (`step2.py`) |
| 05 | `05_barcode_to_variant.qsub` | collapse to per-block barcode→variant dict (`step3.py`), per-block cutoffs from `dictionary_parameters.csv` |
| 06 | `06_dict_stats.qsub` | dictionary-build diagnostic stats *(optional)* |
| 07 | `07_dict_sweep.qsub` | support/fidelity threshold sweep *(optional)* |
| 08 | `08_support_hist.qsub` | per-barcode read-support histogram *(optional)* |
| 09 | `09_support_accuracy_scatter.qsub` | per-barcode support × fidelity scatter *(optional)* |
| 10 | `10_combine_variant_mapped.qsub` | pool variant-mapped reads across a gene group *(optional, `BUILD_COMBINED=1`)* |
| 11 | `11_combined_dictionary.qsub` | build collapsed combined dictionary per group *(optional, `BUILD_COMBINED=1`)* |
| 12 | `12_concatenate_dictionary.qsub` | concatenate per-block dicts per group → `concatenated/` (default) |

Dictionary forms (under `barcode_to_variant/<LIBTYPE>/`): **`blocks/`** (per-block,
canonical), **`combined/`** (per-group, collapsed with one cutoff — optional), and
**`concatenated/`** (per-group, per-block dicts pooled, keeping each block's cutoff
and cross-block bc1 rows — what part 2 consumes).

Steps 06–09 are **dictionary-build** diagnostics (QC of the dictionary itself),
distinct from the cDNA-sample diagnostics in part 2.

## Part 2 — cDNA edit count (`qsub/cdna_edit_count/`)

Applies a dictionary (the per-block-derived **concatenated** dictionary) to cDNA
reads. Edit `qsub/cdna_edit_count/config.sh` for the sample/batches.

| Step | Script | What |
|------|--------|------|
| 01 | `01_edit_count.qsub` | Step1 (expand edited barcodes) → Step2 (anchor template, call edits per position, match bc1 → variant, credit a multi-mapping bc1 to every variant). Writes `<rep>.cDNA_parsed.csv` per batch×rep, then a per-sample edit-position bar plot under `plot/` |
| 02 | `02_cDNA_run_stats.qsub` | summarize the per-rep `Report_line`s into a `cDNA_diagnosis_stats_*.csv` per batch (chain with `qsub -hold_jid <01_jobid>`) |

Outputs go to a separate `cDNA_parsed_<…>/` folder per batch so the canonical
`cDNA_parsed/` is untouched.

**R2 contamination screen.** Parsing is R1-only by default. Set `R2_FILTER=1` for a
batch whose cDNA library is contaminated with dictionary reads (currently
`editRate_071526` only): step 01 then reads R2 alongside R1 and keeps a pair only if
R2 is R1's reverse mate — no dictionary block fixed sequence, and starting with the
reverse complement of `R1_POST_BC2` (derived in-script, with edit-variable positions
wildcarded so edited reads are not preferentially discarded). Contamination and
adapter-dimer rejects are counted separately and reported as extra stats columns.

**Primary edit position.** `PRIMARY_POS` (1-based in the edit template; pinned to 14, the
on-target site) is surfaced per variant as `Primary_pos`/`Primary_edit_ct` and highlighted in
the per-sample plot. Leave `PRIMARY_POS=` empty to auto-detect (most-edited position); it is
tied to `TEMPLATE_SEQ`, so re-derive it if the template changes.

## Conventions

- **bc1-only keying** (`--ignore_bc2`) is the committed Dual_NEBuilder design.
- Two conda envs cover everything: `access` (most steps) and `cutadapt` (step 01).
- **Onboarding:** load the relevant `analysis_plan/*_analysis_plan.md` "Status /
  orientation" block first; the `*_decisions.md` logs hold the *why/when*
  (local-only — gitignored, so not present in a fresh clone).
