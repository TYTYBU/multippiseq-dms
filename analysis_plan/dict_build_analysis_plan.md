# Analysis Plan — multiPPIseq DMS dictionary building

## Status / orientation   (read first)

*Load this file to bring a new session up to speed; this block is the current state, the rest is the detailed method.*

- **Canonical pipeline:** `multiPPIseq_DMS_central/qsub/dict_build/` — step-wise SGE
  job scripts (`00_fastqc` → `05_barcode_to_variant`) sharing `config.sh`, calling the
  Python lineage in `multiPPIseq_DMS_central/script/dict_build/` (`step1_v2.py` →
  bowtie2 → `step2.py` → `step3.py`). Decomposed from the legacy
  `qsub_script/multiPPIseq_DMS_build_dict/multiPPIseq_DMS_build_dict_master.qsub`.
- **Output:** under `<DATA_PATH>/<stage>/<LIBTYPE>/` — e.g. `fastqc/`, `noAdapter/`,
  `processed_fa/`, `bowtie2_out/`, `variant_mapped/`, `barcode_to_variant/`. Final
  dictionaries: `barcode_to_variant/<LIBTYPE>/*.barcode_to_variant.csv.gz` (per block under
  `blocks/`, per gene at top level). **Dictionary-build diagnostics** (steps 06–09, which
  assess the *dictionary itself*) go in `barcode_to_variant/<LIBTYPE>/stats/`. These are
  distinct from **cDNA-sample diagnostics** — the QC of *applying* a finished dictionary to
  a cDNA sample for variant counting — which is a later, separate pipeline (part 2) with its
  own outputs; do not conflate the two.
- **Conda envs:** two cover the pipeline (set in `config.sh`): `access`
  (`/net/apps/conda/miniconda3/envs/access`) provides fastqc, bowtie2, samtools, and
  python+pandas+pysam — used by steps 00/02/03/04/05; `cutadapt`
  (`/net/apps/conda/miniconda3/envs/cutadapt`, v4.6) — step 01 only. Each step script
  `source`s the right one. (No dedicated env exists and this host has no outbound internet
  to build one.)
- **Reference set:** dated extended copies live in
  `multiPPIseq_DMS_central/ref/dict_build/`, all extended with RAD51D blk7/8 + XRCC2
  blk6/7: `r2_fixed_seq_061626.csv` (step 02); `block_ref_seqs_061626.fa` + `.fai` +
  `bowtie2_index/block_ref_seqs_061626.*` (steps 03–04); `variant_sequence_library_lite_061626.csv`
  (step 04, 7917 rows). All four genes now covered end-to-end.
- **Locked parameters:** bc1 len 12–15; bc2 **range mode** 2–7 with `r1_post_bc2`
  anchor `ACCTTTCGACCT` (the r1_barcode2 library is *not* used); r1_linker
  `TATCATAC…TGTAGTACC`. Step3 keying is **bc1-only, committed** (`IGNORE_BC=--ignore_bc2`
  in config; bc1+bc2 and grid sweeps were evaluated and dropped). Thresholds: per-block
  support/fidelity come from `DICT_PARAMS_CSV` (`<DATA_PATH>/dictionary_parameters.csv`)
  when present, else the global `MIN_SUPPORT_READS=3`/`MIN_ACCURACY=0.9`. All in
  `qsub/dict_build/config.sh`; override per-batch via env or by editing config.
- **Library types:** reads come in two construct types — **regular** and
  **Dual_NEBuilder** — processed separately and namespaced by `LIBTYPE`. (Final dicts of
  the two types can be merged downstream via `merge_regular_NEbuilder_dicts.pl`.)
- **Filename convention:** raw fastq must parse as `<GENE>_DMS_<block>_Dictionary_R{1,2}.fastq.gz`
  (a person/type prefix is allowed). `step1_v2.py` parses GENE (uppercase) + block via
  `([^_]+)_DMS_(\d+)_Dictionary`.
- **Decision history & rationale:** `dict_build_decisions.md`.
- **Status:** Active. Current batch = **Ampseq pilot** (`multiPPIseq_DMS_dict_Ampseq_pilot_062626`,
  Dual_NEBuilder; genes RAD51C, RAD51D, XRCC2, XRCC3). Steps 00–04 done. Step 05 **block-level
  only** (`BLOCK_LEVEL_ONLY=1`), bc1-only, regenerated with **tentative per-block thresholds**
  from `dictionary_parameters.csv` (support 3 for all blocks; fidelity 0.05–0.125 per block).
  Yield now thousands of barcodes/block (e.g. XRCC3_DMS_3 19,220; RAD51C_DMS_2 4,736).
  Dictionary diagnostics (06/08/09) rerun against these dictionaries.
- **Last updated:** 2026-06-23

---

## Pipeline steps

Each step is one SGE job script under `qsub/dict_build/`, sourcing `config.sh`.
Run in order; step 01 is conditional. The preflight runs once per batch before step 00.

### Preflight — raw-data filename check & layout
- **Inputs:** delivered raw `*.fastq.gz` for the batch.
- **Method:** verify every file parses as `<GENE>_DMS_<block>_Dictionary_R{1,2}.fastq.gz`
  with `step1_v2.py`'s regex `([^_]+)_DMS_(\d+)_Dictionary` — GENE uppercase (the bc2
  library match, when used, is case-sensitive), `block` digits-only, `Dictionary`
  capitalized. A person/library-type prefix before the gene is allowed (e.g.
  `NR-TV__XRCC3_DMS_5_Dictionary_R1.fastq.gz`). Rename any non-conforming files. Place
  files under `raw_data/<LIBTYPE>/` (`Dual_NEBuilder` or `regular`).
- **Output:** conforming, correctly-located raw reads; nothing else parses downstream.
- **Rationale:** `step1_v2.py` derives gene/block from the R1 filename; a bad name is
  silently mis-parsed or skipped. (Candidate to become a dedicated skill/script later.)

### Step 00 — FastQC
- **Inputs:** `raw_data/<LIBTYPE>/*.fastq.gz`.
- **Method:** `fastqc` per file (`00_fastqc.qsub`).
- **Output:** `fastqc/<LIBTYPE>/*_fastqc.html|zip`.
- **Rationale:** QC, and specifically to read the **Adapter Content** module to decide
  whether step 01 is needed. The cutadapt decision is based **solely on Adapter Content** —
  the Overrepresented-sequences module is ignored here (for these amplicon libraries it
  flags the construct's own constant regions, not sequencing adapters).

### Step 01 — Adapter trim (optional)
- **Inputs:** raw paired reads.
- **Method:** `cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 20` (`01_cutadapt.qsub`).
  Run **only** if FastQC Adapter Content shows read-through. **This step sets the
  downstream input directory:** after running it you must set `USE_NOADAPTER=1` so step 02
  reads trimmed reads from `noAdapter/<LIBTYPE>/`; if skipped, leave `USE_NOADAPTER=0`
  (default) and step 02 reads raw reads from `raw_data/<LIBTYPE>/`.
- **Output:** `noAdapter/<LIBTYPE>/*_R{1,2}.noAdapter.fastq.gz`.
- **Rationale:** regular-type libraries historically need no trimming; Dual_NEBuilder is
  evaluated per batch from FastQC.

### Step 02 — Barcode / oligo extraction → FASTA
- **Inputs:** paired reads (raw, or trimmed if `USE_NOADAPTER=1`); `r2_fixed_seq.csv`.
- **Method:** `step1_v2.py` (`02_extract_barcodes.qsub`). Extracts bc1 (12–15 nt) and bc2
  via **range mode** (2–7 nt delimited by `r1_post_bc2` anchor `ACCTTTCGACCT`), splits R1
  on `r1_linker`, validates R2 against the per-gene/block fixed-sequence library.
- **Output:** `processed_fa/<LIBTYPE>/<prefix>.fasta.gz` (read header carries barcodes;
  sequence is the oligo).
- **Rationale:** range mode avoids depending on the r1_barcode2 library, which does not
  cover all genes/blocks in newer batches.

### Step 03 — Align to block references
- **Inputs:** `processed_fa/<LIBTYPE>/*.fasta.gz`; bowtie2 index of `block_ref_seqs.fa`.
- **Method:** `bowtie2 --local --very-sensitive-local --no-unal` → `samtools sort`
  (`03_bowtie2_align.qsub`).
- **Output:** `bowtie2_out/<LIBTYPE>/<prefix>.block_ref_seqs.bam`.

### Step 04 — Map reads to variants
- **Inputs:** per-block BAM; `block_ref_seqs.fa`; `variant_sequence_library_lite.csv`.
- **Method:** `step2.py` (`04_variant_map.qsub`; `access` env / pysam).
- **Output:** `variant_mapped/<LIBTYPE>/blocks/<prefix>.block_ref_seqs.csv.gz`
  (read-level barcode↔variant).

### Step 05 — Collapse to barcode→variant dictionary
- **Inputs:** read-level variant maps (per block, and per gene after concatenating a
  gene's blocks).
- **Method:** `step3.py $IGNORE_BC --min_support_reads <S> --min_accuracy <A>`
  (`05_barcode_to_variant.qsub`). Keying is **bc1-only** (`IGNORE_BC=--ignore_bc2`,
  committed). A barcode is assigned its **max-support variant** if that variant has ≥ S
  reads (support) and accounts for ≥ A of the barcode's reads (fidelity = best_ct/total_ct).
  **Per-block thresholds:** when `DICT_PARAMS_CSV` exists, each block's S/A are read from it
  (cols `Gene,Block_number,proposed_fidelity_cutoff,proposed_support_read_cutoff`; awk strips
  CR/BOM); otherwise the global `MIN_SUPPORT_READS`/`MIN_ACCURACY` apply to all blocks. For
  the current Dual_NEBuilder batch `BLOCK_LEVEL_ONLY=1`: only per-block dictionaries are
  built; the per-gene merge is skipped while thresholds are still being explored. Per-block
  values are tentative and may change.
- **Output:** `barcode_to_variant/<LIBTYPE>/blocks/*.barcode_to_variant.csv.gz` (per
  block) and `barcode_to_variant/<LIBTYPE>/*.barcode_to_variant.csv.gz` (per gene, when not
  block-level-only), plus `*.raw_dict.csv.gz`.

## Dictionary-build diagnostics (steps 06–09)

These steps assess the **dictionary itself** (read-survival funnel, threshold sensitivity,
per-barcode support depth, support×fidelity structure). They are QC for *building* the
dictionary and are kept separate from **cDNA-sample diagnostics** (part 2: QC of applying a
finished dictionary to a cDNA sample to count variant edits), which has its own scripts and
output location. All step 06–09 outputs live under `barcode_to_variant/<LIBTYPE>/stats/`.

### Step 06 — Dictionary-build diagnostic stats (optional)
- **Inputs:** raw R1 fastqs; step-02 `processed_fa`; step-04 `variant_mapped/blocks`;
  step-05 `barcode_to_variant/blocks` (`*.raw_dict.csv.gz` + `*.barcode_to_variant.csv.gz`).
- **Method:** `multiPPIseq_DMS_dict_stats.py` (`06_dict_stats.qsub`; Python port of the
  legacy `DMS_dict_stats.pl`). Discovers blocks from the raw R1 names and counts, per block:
  `total_reads`, `reads_with_barcode`, `reads_with_barcode_and_variant`,
  `reads_in_faithful_barcode_pairs`, `faithful_barcode_pairs`.
- **Output:** one CSV `barcode_to_variant/<LIBTYPE>/stats/dict_stats.csv` (one row per block).
- **Rationale:** read-survival funnel per stage + final dictionary yield, for QC and
  parameter tuning. `reads_in_faithful_barcode_pairs` (was the Perl's ambiguous
  `treads_passed_QC`) = the raw reads belonging to the barcode pairs that passed the
  step-05 support/accuracy filters, i.e. the reads captured by the final dictionary.

### Step 07 — Threshold sweep (parameter tuning; optional)
- **Inputs:** step-04 `variant_mapped/<LIBTYPE>/blocks` per-read CSVs only.
- **Method:** `multiPPIseq_DMS_dict_sweep.py` (`07_dict_sweep.qsub`). Reuses the step3
  functions to load + redistribute counts **once per block**, then re-counts the full
  `--support_values` × `--accuracy_values` grid in memory. Uses `IGNORE_BC` from config
  (bc1-only).
- **Output:** one CSV `barcode_to_variant/<LIBTYPE>/stats/dict_param_sweep_bc1.csv` with columns
  `block, min_support_reads, min_accuracy, variant_mapped_reads,
  reads_in_faithful_barcode_pairs, frac_reads_retained, observed_barcode_pairs,
  faithful_barcode_pairs, frac_pairs_retained`.
- **Rationale:** choose Dual_NEBuilder support/fidelity thresholds. `variant_mapped_reads` =
  reads assigned a variant (read-fraction denom); `observed_barcode_pairs` = distinct bc
  pairs seen (pair-fraction denom). The sweep reproduces the step-05 dictionary exactly at
  matching thresholds (validated). (bc1-only is committed; the earlier bc1+bc2/grid sweeps
  were dropped — see decisions log.)

### Step 08 — Per-barcode read-support distribution (QC; optional)
- **Inputs:** step-04 `variant_mapped/<LIBTYPE>/blocks` per-read CSVs only.
- **Method:** `multiPPIseq_DMS_dict_support_hist.py` (`08_support_hist.qsub`). Counts reads
  per unique barcode (keying from `IGNORE_BC`, bc1-only) and histograms the per-barcode
  support depth.
- **Output:** under `barcode_to_variant/<LIBTYPE>/stats/support_hist/`: `support_hist.hist.csv`
  (long: block, support_reads, n_barcodes; plus an `ALL` block), `support_hist.summary.csv`
  (n_unique_barcodes, n_reads, mean/median/max, barcodes/reads ≥ 1..5,10 thresholds), and
  per-block + `ALL` PNG histograms (log-y, bins 1..xmax then `≥xmax+1`).
- **Rationale:** quantify under-sequencing. The pilot is deliberately under-sequenced; this
  step measures the read-support depth per bc1. (Pilot ALL: ~214k unique bc1, bimodal — ~90k
  singletons vs ~70k with ≥31 reads.)

### Step 09 — Per-barcode support × fidelity scatter (QC / parameter context; optional)
- **Inputs:** step-04 `variant_mapped/<LIBTYPE>/blocks` per-read CSVs only.
- **Method:** `multiPPIseq_DMS_dict_support_accuracy_scatter.py` (`09_support_accuracy_scatter.qsub`).
  Reuses the step3 load+redistribute (so per-barcode numbers match the dictionary), then for
  every barcode plots one dot on log2–log2 axes: x = support reads of the max-support variant
  (`best_ct`), y = fidelity (`best_ct/total_ct`). Dots are colored by local barcode density
  (an 80×80 log2 `histogram2d`; raw overplotting would be unreadable). Tie barcodes are kept
  for the cloud but flagged via `n_tie` in the CSV (`--drop_ties` to exclude). Uses `IGNORE_BC`
  (bc1-only). **The plotted cloud is parameter-independent** (the input is the step-04
  `variant_mapped` reads, not the dictionary; per-barcode `(best_ct, accuracy)` is recomputed
  here). Only the dashed reference lines depend on parameters: when `DICT_PARAMS_CSV` is passed
  (`--params_csv`), each per-block panel draws **that block's chosen support/fidelity cutoff**
  so the decision can be seen against the cloud (the `ALL` panel omits the box since cutoffs
  differ per block); otherwise a generic `--min_support`/`--min_accuracy` guide (default 3 × 0.9)
  is drawn. Because the lines are a pure overlay on a parameter-independent cloud, annotating the
  chosen cutoffs is informative, not circular.
- **Output:** under `barcode_to_variant/<LIBTYPE>/stats/support_accuracy_scatter/`:
  `support_accuracy.perbarcode.csv.gz` (`block, bc_pair, best_ct, total_ct, accuracy, n_tie`)
  and per-block + `ALL` PNG scatters.
- **Rationale:** read off, per block, where the barcode mass sits so support/fidelity cutoffs
  can be picked. Pilot bc1-only: the cloud sits at high support but low fidelity (reads smeared
  across many variants per bc1); some blocks (e.g. RAD51D_DMS_8) are bimodal in fidelity,
  reflecting two bc1 collision-multiplicity populations (~6 vs ~tens of variants per bc1).

## Combined dictionaries (steps 10–11; OPTIONAL — `BUILD_COMBINED=1`)
The cDNA pipeline keys on the per-block dictionaries (step 05); the combined dictionaries are
an optional convenience and the whole branch is gated behind `BUILD_COMBINED` (default 0).

### Step 10 — Combine variant-mapped reads across blocks/genes (optional)
- **Inputs:** step-04 `variant_mapped/<LIBTYPE>/blocks` per-read CSVs.
- **Method:** `10_combine_variant_mapped.qsub` concatenates (gzip multi-member; valid) all
  blocks of the genes in each `GENE_GROUPS` entry (COMMD* together; RAD51D+XRCC2; RAD51C+XRCC3)
  into one combined variant-mapped file per group, keeping the per-read format and the
  `_Dictionary.block_ref_seqs.csv.gz` suffix so steps 09/11 consume it unchanged. The variant
  column keeps its gene prefix (`GENE---variant`) so genes stay distinguishable after merge.
- **Output:** `variant_mapped/<LIBTYPE>/combined/<group>_DMS_Dictionary.block_ref_seqs.csv.gz`.
- **QC:** run step 09 on the combined dir (`SCATTER_IN_DIR=…/combined`,
  `SCATTER_OUT_DIR=…/support_accuracy_scatter_combined`, `SCATTER_PARAMS_CSV=` to suppress the
  box) to choose per-group cutoffs before building.

### Step 11 — Build combined dictionaries (optional)
- **Inputs:** step-10 combined variant-mapped files; per-group cutoffs from
  `DICT_PARAMS_COMBINED_CSV` (`Group,proposed_fidelity_cutoff,proposed_support_read_cutoff`).
- **Method:** `11_combined_dictionary.qsub` runs step3.py per group, looking up its cutoff by
  `Group` label (falls back to global `MIN_SUPPORT_READS`/`MIN_ACCURACY` if absent), with
  `IGNORE_BC` (bc1-only).
- **Output:** `barcode_to_variant/<LIBTYPE>/combined/<group>.barcode_to_variant.csv.gz`
  (+ `.raw_dict.csv.gz`).
- **Pilot result (support≥3):** `RAD51C_XRCC3` fidelity≥0.125 → 42,475 barcodes;
  `RAD51D_XRCC2` fidelity≥0.0625 → 41,289 barcodes.

## Notes
- `--ignore_bc2` (bc1-only) is the **committed Dual_NEBuilder design**, not a defect; the
  bc1+bc2 and grid sweeps were evaluated and dropped. Low pilot yield is expected from
  under-sequencing. See decisions log.
- **Dictionary diagnostics vs cDNA-sample diagnostics:** steps 06–09 QC the dictionary build
  itself. A later, separate pipeline (part 2) QCs *applying* a dictionary to a cDNA sample
  for variant counting — different scripts, different outputs. Keep the two sets distinct.
- Legacy Perl lineage (`step1.pl`, `step3a/b/c.pl`, `partition_barcodes.pl`, the
  array-job `step2_arrayJob.qsub`, `DMS_dict_stats.pl`) is superseded by the Python lineage
  above. See the decisions log.
