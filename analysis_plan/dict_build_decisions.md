# Decisions log — multiPPIseq DMS dictionary building

Dated, append-only. Newest entries on top. One entry per real decision.

---

## 2026-06-23 — Combined dictionaries built but kept optional (cDNA uses per-block)

**Context.** By convention reads are pooled across blocks and across paired genes into
combined groups (COMMD* together; RAD51D+XRCC2; RAD51C+XRCC3). We built combined
variant-mapped files (step 10) and combined dictionaries (step 11) from them. But the cDNA
counting pipeline downstream will key on the **per-block** dictionaries, not the combined ones.

**Decision.** Everything touching the combined dictionaries — step 10 (combine
variant-mapped), the combined-files scatter, and the new step 11 (combined dictionary build)
— is **gated behind `BUILD_COMBINED` (default 0)**. Per-block dictionaries (step 05) run
unconditionally. Combined per-group cutoffs come from `dictionary_parameters_combined.csv`
(cols `Group,proposed_fidelity_cutoff,proposed_support_read_cutoff`), falling back to the
global 3 × 0.9 if a group is absent.

**Evidence / pointer.** Built with Tian's chosen combined cutoffs (support≥3):
`RAD51C_XRCC3` fidelity≥0.125 → 42,475 barcodes; `RAD51D_XRCC2` fidelity≥0.0625 → 41,289
barcodes, at
`barcode_to_variant/Dual_NEBuilder/combined/`. Gene prefix retained in the variant column
(`RAD51C---K194N`) so genes stay distinguishable post-merge. New script
`qsub/dict_build/11_combined_dictionary.qsub`.

---

## 2026-06-23 — Step-09 scatter: show the chosen per-block cutoffs (supersedes prior entry)

**Context.** The prior entry removed the per-block threshold box to avoid "circularity."
On reflection that reasoning was wrong: step 09's input is the step-04 `variant_mapped`
reads, **not** the dictionary, and the per-barcode `(best_ct, accuracy)` cloud is recomputed
independently of any threshold. So the dots/density/`perbarcode.csv` are identical regardless
of parameters — the only parameter-dependent element is the dashed reference line. Drawing the
chosen cutoffs is therefore a pure annotation on a parameter-independent cloud, which is
informative (you see how the decision sits on the data), not circular.

**Decision.** Restored the `--params_csv` feature (recovered from commit 77ef2e9). Step 09
now draws **each block's chosen support/fidelity cutoff** (from `DICT_PARAMS_CSV`) on its
panel; the `ALL` panel omits the box since cutoffs differ per block. Falls back to the generic
`--min_support`/`--min_accuracy` guide (3 × 0.9) when no params CSV is given.

**Status / pointer.** Applied; scatter regenerated with per-block boxes. Supersedes
"2026-06-23 — Step-09 scatter stays parameter-agnostic" below.

---

## 2026-06-23 — Step-09 scatter stays parameter-agnostic (it *selects* the parameters) [SUPERSEDED]

**Context.** While wiring per-block thresholds I briefly made step 09 draw each block's
chosen cutoff (from `dictionary_parameters.csv`) as its threshold box. That inverts the
workflow: the support×fidelity scatter is the **input** Tian uses to *choose* the per-block
cutoffs, so drawing the decided cutoffs back onto it is circular and misleading.

**Decision.** Step 09 is produced **independently of the decided parameters**. Reverted the
`--params_csv` per-block-box feature (removed from the script and `09_support_accuracy_scatter.qsub`).
The dashed reference lines are a generic reading guide only (script defaults 3 × 0.9); the
qsub no longer passes any chosen thresholds. The scatter is run first; per-block cutoffs are
read off it and recorded in `dictionary_parameters.csv`.

**Status / pointer.** Applied. Scatter regenerated parameter-agnostic. Methods file Step 09
re-labeled "parameter selection" and notes the independence explicitly.

---

## 2026-06-23 — Per-block dictionary thresholds from CSV; pilot block dicts regenerated

**Context.** With bc1-only committed, the pilot needed real per-block support/fidelity
cutoffs rather than one global 0.9 fidelity (which yielded 0–8 barcodes/block — far too
strict for the smeared bc1 clouds seen in steps 08/09). Tian supplied tentative per-block
values in `<DATA_PATH>/dictionary_parameters.csv` (cols `Gene, Block_number,
proposed_fidelity_cutoff, proposed_support_read_cutoff`): support 3 for all 9 blocks;
fidelity 0.0625 (RAD51C/RAD51D, XRCC2_6), 0.05 (XRCC2_7), 0.125 (XRCC3 all). Values are
tentative and may be updated.

**Decision.** Step 05 now reads per-block thresholds from `DICT_PARAMS_CSV`
(config default `<DATA_PATH>/dictionary_parameters.csv`) when it exists, falling back to the
global `MIN_SUPPORT_READS`/`MIN_ACCURACY` otherwise. The per-block dictionaries in
`barcode_to_variant/Dual_NEBuilder/blocks/` were regenerated in place with these values.
Yield rose to thousands of barcodes/block (XRCC3_DMS_3 19,220; RAD51D_DMS_7 14,944;
RAD51C_DMS_2 4,736). These cutoffs are deliberately permissive (top variant need only hold
~5–12.5% of a barcode's reads), so many called bc1 are multi-molecule collisions —
acceptable for a tentative pilot dictionary.

**Tooling / fix.** `config.sh` `DICT_PARAMS_CSV`; `05_barcode_to_variant.qsub` per-block
loop with awk lookup that strips CR/BOM (the CSV arrived UTF-8-BOM + CRLF). The CSV was
subsequently cleaned to plain ASCII/LF. Diagnostics steps 06/08/09 rerun against the new
dictionaries.

**Naming note.** Steps 06–09 are **dictionary-build diagnostics** (assessing the dictionary
itself), explicitly distinct from future **cDNA-sample diagnostics** (part 2: QC of applying
a dictionary to a cDNA sample for variant counting). Documented in the methods file.

**Status / pointer.** Applied. Per-block thresholds tentative; revisit as the pilot is
analyzed. Per-gene merge still skipped (`BLOCK_LEVEL_ONLY=1`).

---

## 2026-06-17 — Committed to bc1-only keying (`--ignore_bc2`) going forward

**Context.** Evaluated three keyings for the Dual_NEBuilder dictionary — bc1+bc2, the
standard grid sweep, and bc1-only — plus a bc1-only low-fidelity sweep (support 3, accuracy
0.1–0.5). Even at accuracy 0.1 only ~36% of bc1 are retained and retention collapses to
<0.3% by 0.5; no usable plateau. bc1-only is the intended Dual_NEBuilder design and the path
forward; low yield is expected from undersequencing.

**Decision.** bc1-only committed. Centralized the flag as `IGNORE_BC=--ignore_bc2` in
`config.sh`; steps 05/07/08/09 reference `$IGNORE_BC`. The bc1+bc2 and grid sweep result
CSVs were deleted; the surviving sweep is `stats/dict_param_sweep_bc1.csv`. Set `IGNORE_BC=""`
only to revisit bc1+bc2.

**Status / pointer.** Applied in config.sh + 05/07/08/09 qsubs. Support/fidelity thresholds
themselves still being explored (now per-block via CSV; see 2026-06-23 entry).

---

## 2026-06-17 — Threshold sweep (step 07); bc1+bc2 vs ignore_bc2 for Dual_NEBuilder

**Context.** Dual_NEBuilder step-05 yield was near-zero at S=3/A=0.9. Built a sweep to see
whether any (min_support_reads × min_accuracy) choice rescues it, and to compare bc1-only
vs bc1+bc2 keying.

**Tooling.** `script/dict_build/multiPPIseq_DMS_dict_sweep.py` + `qsub/dict_build/07_dict_sweep.qsub`.
Reuses step3's load/neighbor/redistribute (imported, registered in `sys.modules["step3"]`
so multiprocessing workers resolve) — redistribute once per block, re-count the grid.
New column names: `variant_mapped_reads` (was `reads_with_barcode_and_variant`; read-frac
denom) and `observed_barcode_pairs` (distinct bc pairs; pair-frac denom). At S=3/A=0.9 it
reproduces the step-05 dictionary exactly (XRCC3_DMS_3: valid=8, reads_in_valid=25).

**Findings (9 blocks, summed).**
- `--ignore_bc2` (the intended design): low yield across the whole grid — even the most
  permissive S=1/A=0.6 keeps only 6.1% of pairs and 0.11% of reads; S≥2 → ~0. **This is
  expected**: the pilot is deliberately under-sequenced, so per-barcode read support is low.
- bc1+bc2 (exploratory comparison only, *not* the design): higher at S=1 (40.7% pairs, 5.4%
  reads) but still collapses at S≥3 (S=3/A=0.6: 0.36% pairs; S=3/A=0.9: 0.04%).
- Accuracy is a *minor* lever (A=0.6→1.0 barely moves yield at fixed S); **support is the
  dominant factor** — most barcodes are singletons/low-coverage, consistent with
  undersequencing rather than a threshold or keying problem.

**Decision.** `--ignore_bc2` stays — it is what the experiment design intends; do NOT treat
the low yield as a defect to fix by switching to bc1+bc2. Keep step 05 outputs as-is
(block-level, ignore_bc2). Next: characterize the per-barcode read-support distribution
(histogram, step 08) to quantify the undersequencing. Sweep CSVs: `dict_param_sweep.csv`
(ignore_bc2, the design), `dict_param_sweep_bc1bc2.csv` (bc1+bc2, exploratory only).

---

## 2026-06-17 — Port DMS_dict_stats.pl to Python; rename ambiguous column; fix Perl bug

**Context.** Diagnostic stats for the per-block dictionaries were produced by the legacy
`script/multiPPIseq_DMS_dict_build/DMS_dict_stats.pl` (wrapped by
`qsub_script/.../multiPPIseq_DMS_dict_stats.qsub`). Wanted a Python version in central,
CSV output, and a clearer name for the ambiguous `treads_passed_QC` column.

**Decision.**
- New `script/dict_build/multiPPIseq_DMS_dict_stats.py` + `qsub/dict_build/06_dict_stats.qsub`
  (step 06; `access` env). Discovers blocks from raw R1 names, computes the per-stage read
  funnel, parallelizes over blocks, writes a single CSV
  `barcode_to_variant/<LIBTYPE>/dict_stats.csv` (one row per block). Columns: `block`,
  `total_reads`, `reads_with_barcode`, `reads_with_barcode_and_variant`,
  `reads_in_faithful_barcode_pairs`, `faithful_barcode_pairs`.
- **Renamed** `treads_passed_QC` → **`reads_in_faithful_barcode_pairs`**. It was never a
  per-read QC pass count; it is the sum of raw read support over the barcode pairs that
  passed the step-05 support/accuracy filters (reads captured by the final dictionary).

**Bug found + fixed (both scripts).** The raw dict has one row per `(bc_pair, variant)`
(verified on a real file: 604,758 lines vs 13,685 unique bc_pairs; one pair on up to 351
lines; col4 = read count, sum of col4 = total reads). The Perl built
`$barcode_pair_hash{$bc_pair} = $raw_read_ct` with plain `=`, so it kept only the *last*
variant's count per pair instead of the pair total — undercounting `reads_passed_QC` for
any multi-variant pair. Fixed the Perl to `+=`; the Python port sums correctly by design.
Effect is invisible on the current Dual_NEBuilder batch (QC-passing pairs are nearly all
single-variant, accuracy ~1.0) but matters once bc2 is used / accuracy is loosened. The two
implementations now agree.

---

## 2026-06-17 — Step 04/05 run; step2.py basename fix; block-level-only step 05

**Context.** Variant library and all 061626 references complete; steps 02/03 done for all 9
Dual_NEBuilder samples. Ready to run step 04 (variant map) + step 05 (barcode→variant dict).

**step2.py path-parsing bug (fixed).** `step2.py` parsed gene/block via
`re.compile('([^_]+)_DMS_(\\d+)_Dictionary').search(bam_path)` on the *full path*. Because
`[^_]+` spans `/` and `.`, the real BAM dir `.../bowtie2_out/Dual_NEBuilder/` would make it
parse gene as `NEBuilder/XRCC2` and crash on `fasta.fetch`. Same class of bug as the
step1_v2 fix. Changed to parse `os.path.basename(bam_path)`. Ported to the original
`script/multiPPIseq_DMS_dict_build/multiPPIseq_DMS_dict_build_step2.py` (kept byte-identical
to the central copy). Smoke-tested on a 20k-read XRCC2_blk_6 subset placed under a
`Dual_NEBuilder/` dir: parsed gene XRCC2 block 6, loaded 102 variant oligos, mapped reads to
variants (e.g. `XRCC2---A63R`, `P68*`) — confirms the padded variant_seq also works.

**Decision: block-level-only step 05 for now.** Dual_NEBuilder step3 parameters are still
being explored, so step 05 builds **only per-block dictionaries** and skips the per-gene
merge, gated by `BLOCK_LEVEL_ONLY=1` (new config default; set 0 to also build per-gene).
Kept `--min_support_reads 3 --min_accuracy 0.9` but flagged as provisional/subject to change.
Submitted step 04 (job 234334) → step 05 (job 234335, `-hold_jid 234334`).

---

## 2026-06-16 — Extend variant library for RAD51D/XRCC2 from validation-library xlsx

**Context.** Step 04 (`step2.py`) needs `variant_sequence_library_lite.csv`
(`gene_aa_str, block_number, variant_seq`); the 121725 copy lacks RAD51D/XRCC2. Source:
`0326_RAD51D_XRCC2_validation_library.xlsx`, sheets "XRCC2 Validation Library" /
"RAD51D Validation Library", columns `gene_aa_str`, `Block Number`, `Oligo Block for mutation`.

**Decision.** Built `multiPPIseq_DMS_central/ref/dict_build/variant_sequence_library_lite_061626.csv`
= 121725 rows + 445 new. For each new variant, `variant_seq = "Oligo Block for mutation" +
block_ref 3' tail` padded to block-ref length. Reason: step2.py compares
`variant_seq[:len(ref)]` position-by-position from index 0, so variant_seq must align at
pos 0 and be ≥ ref length. The oligos already start at the block ref's position 0 but are
20–24 nt shorter than the (061626) block refs; the ref region beyond the oligo is constant,
so appending it is correct and adds no spurious mismatches. `config.sh` `VARIANT_LIB`
updated.

**Evidence.** All 445 new variants: 2–3 mismatches vs block ref (the variant codon),
mean identity ~98.7%, none shorter than ref. Counts: XRCC2 blk6=102, blk7=110; RAD51D
blk7=100, blk8=133. Matches the existing library's convention (variant aligned at pos 0,
1–4 mismatches).

**Status / pointer.** Applied. Reference set now covers all four pilot genes end-to-end
(steps 02–05).

---

## 2026-06-16 — Fix step1_v2.py gene parsing to use basename (Dual_NEBuilder dir bug)

**Context.** Smoke-testing steps 02→03 on peek data, step 02 failed with "R2 fixed
sequence not found for gene first1000/RAD51C". `step1_v2.py` ran its gene/block regex
`([^_]+)_DMS_(\d+)_Dictionary` on the **full input path**; `[^_]+` stops at underscores, so
an underscore in a parent directory name leaks into the gene. The real input dir
`raw_data/Dual_NEBuilder/` would have parsed gene as `NEBuilder/RAD51C` — i.e. the full
production run would have failed identically, not just the peek.

**Decision.** Patched the central copy
(`script/dict_build/multiPPIseq_DMS_dict_build_step1_v2.py`) to parse
`os.path.basename(infile)` instead of `infile` (added `import os`). Same fix ported back to
the original `script/multiPPIseq_DMS_dict_build/multiPPIseq_DMS_dict_build_step1_v2.py`
(now byte-identical to the central copy).

**Evidence.** `re.search(rx, full_path)` → `NEBuilder/RAD51C`; `re.search(rx, basename)` →
`RAD51C`. After the fix, smoke test parsed RAD51C/XRCC2/RAD51D correctly, extracted
~80–82% of pairs, and aligned at 99.6–100% each to its own block.

**Status / pointer.** Applied; steps 02+03 then submitted on full data (jobs 234332/234333).

---

## 2026-06-16 — Extend block_ref_seqs (+ faidx + bowtie2 index) for RAD51D/XRCC2

**Context.** Steps 03–04 align/variant-map against `block_ref_seqs.fa`; the 121725 copy
lacks RAD51D blk7/8 and XRCC2 blk6/7. Tian provided the four block reference sequences, but
each still carried the R2 fixed sequence and a 5' NNN prefix — whereas existing block refs
begin **immediately 3' of** the R2 fixed seq (e.g. RAD51C_blk_1 starts `ACCTTCTGTTCAGCACTAGATGATA`,
right after fixed `GCATACCCAGGGCTTCATAATC`).

**Decision.** Stripped the R2 fixed seq + everything 5' of it from each provided sequence,
then copied `block_ref_seqs.fa` → `multiPPIseq_DMS_central/ref/dict_build/block_ref_seqs_061626.fa`
and appended the four trimmed entries (`>RAD51D_blk_7/8`, `>XRCC2_blk_6/7`). Built
`samtools faidx` and `bowtie2-build` (env `access`). `config.sh` `BLOCK_REF_FA` and
`BOWTIE2_INDEX` now point to the dated copies.

**Evidence.** For all four, the fixed seq occurred exactly once, preceded only by NNN(N).
Trimmed lengths: RAD51D_blk_7=174, RAD51D_blk_8=172, XRCC2_blk_6=167, XRCC2_blk_7=171
(comparable to existing ~125–200 nt blocks). Fasta now has 22 entries; `.fai` and 6 `.bt2`
index files generated.

**Status / pointer.** Applied. Remaining gap: the variant library still lacks RAD51D/XRCC2
(needed for step 04).

---

## 2026-06-16 — Extend r2 fixed-seq lib for RAD51D/XRCC2 → dated central copy

**Context.** Step 02 (`step1_v2.py`) needs an R2 fixed sequence per gene/block; reads with
no match are skipped. The canonical `data/multiPPIseq_DMS_dict_121725/r2_fixed_seq.csv`
covers COMMD*, RAD51C, XRCC3 but not RAD51D (blk7/8) or XRCC2 (blk6/7), which are in the
Ampseq pilot. Tian provided candidate fixed seqs for the four missing gene/blocks.

**Decision.** Verified each provided sequence against its gene×block R2 reads, then copied
the canonical CSV into the workspace as
`multiPPIseq_DMS_central/ref/dict_build/r2_fixed_seq_061626.csv` and appended the four
rows. `config.sh` `R2_FIXED_SEQ_LIB` now points to this dated copy.

**Evidence.** Fraction of R2 reads (first 100k) containing the fixed seq:
| gene/block | fixed seq | hit rate |
|---|---|---|
| RAD51D/7 | AGTTCTCGGCTTTCCCCGTGAAT | 98.3% |
| RAD51D/8 | GCAGGTGACTGGTTCTTCAGGAACT | 94.6% |
| XRCC2/6 | CCTGTGCATGGTGATATTCTTGAATTTCAT | 97.0% |
| XRCC2/7 | AGAAATAATCAAATACTGCCTGGGA | 97.8% |

**Status / pointer.** Applied. Note `block_ref_seqs.fa`, `variant_sequence_library_lite.csv`,
and the bowtie2 index still lack RAD51D/XRCC2 coverage — must be extended before steps
03–05 run for those genes.

---

## 2026-06-16 — Tool envs: `access` + `cutadapt` (no new env built)

**Context.** Wanted one conda env with fastqc + cutadapt + bowtie2 + samtools + pysam. No
single existing env had all; tried to `conda create` a new `multiPPIseq` env but this host
(muhee) has no outbound internet (conda/pypi/google all time out, no proxy), so the build
failed. fastqc and bowtie2 also aren't conda packages here — only system binaries.

**Decision.** Use two existing apps-conda envs instead of building a new one:
`access` (fastqc, bowtie2, samtools, python 3.12 + pandas 2.2.2 + pysam 0.22.1) for steps
00/02/03/04/05, and `cutadapt` (v4.6) for step 01. Each step script activates the right
env via `source "$CONDA_BASE/bin/activate" "$ENV_*"`; `CONDA_BASE/ENV_MAIN/ENV_CUTADAPT`
live in `config.sh`. Dropped the earlier reliance on the home `bio` env and hardcoded
binary paths.

**Evidence.** With `access` active: fastqc 0.11.5, bowtie2 (PATH), samtools 1.21,
`import pandas, pysam` OK. With `cutadapt` active: cutadapt 4.6. Smoke-tested activation in
a clean login shell.

**Status / pointer.** Applied in `config.sh` + all `qsub/dict_build/*.qsub`. If a unified
env is ever wanted, build it where there is internet (see decisions on host networking).

---

## 2026-06-16 — Ampseq pilot (Dual_NEBuilder): skip cutadapt — no adapter read-through

**Context.** New-batch policy is to run FastQC (step 00) first and only run cutadapt
(step 01) if adapter read-through is detected. Ran step 00 on the 18 Dual_NEBuilder pilot
files (job 234331).

**Decision.** **Skip step 01** for this batch; run step 02 on raw reads (`USE_NOADAPTER=0`).
Policy: the cutadapt decision is based **solely on the Adapter Content module**;
Overrepresented sequences are ignored (per Tian).

**Evidence.** FastQC "Adapter Content" = PASS on all 18 files; max standard-Illumina
adapter fraction ≤ 0.09% across every file — far below any trimming threshold. (The
Overrepresented-sequence flags are the construct's own constant regions — R2 fixed
sequence, R1 bc1+r1_linker — not adapters, and are not part of this decision.) Consistent
with the legacy master.qsub "NO NEED TO TRIM ADAPTERS" note.

**Downstream impact.** Because step 01 is skipped, `USE_NOADAPTER` stays 0 and step 02
reads `raw_data/Dual_NEBuilder/`. If cutadapt were run, step 02's input would instead be
`noAdapter/Dual_NEBuilder/` (`USE_NOADAPTER=1`).

**Status / pointer.** Applied. `fastqc/Dual_NEBuilder/` holds the reports. If a future
batch shows real adapter read-through, run `01_cutadapt.qsub` then set `USE_NOADAPTER=1`.

---

## 2026-06-16 — Decompose master.qsub into step-wise scripts; create central workspace

**Context.** Dictionary building lived in one monolithic
`qsub_script/multiPPIseq_DMS_build_dict/multiPPIseq_DMS_build_dict_master.qsub` mixing
FastQC, cutadapt, extraction, alignment, variant mapping, and dictionary collapse — hard
to run a single stage (e.g. FastQC only) for a new batch. Starting a new Ampseq pilot
batch and wanting a reusable home for the whole multiPPIseq DMS pipeline (dictionary
building + cDNA edit count).

**Decision.** Created `multiPPIseq_DMS_central/` as the pipeline workspace
(`script/`, `qsub/`, `analysis_plan/`). Split the master qsub into six step scripts
`00_fastqc` … `05_barcode_to_variant` under `qsub/dict_build/`, all sourcing a shared
`config.sh`. **Copied** (not moved) the canonical Python scripts `step1_v2.py`,
`step2.py`, `step3.py` into `script/dict_build/` — originals in
`script/multiPPIseq_DMS_dict_build/` left in place so existing callers don't break.

**Evidence.** master.qsub's active (uncommented) path is the Python lineage with bc2
range mode and step3 `--ignore_bc2 --min_support_reads 3 --min_accuracy 0.9`; that is what
was carried into the step scripts and `config.sh`.

**Status / pointer.** Applied. See `dict_build_analysis_plan.md` for the methods.

---

## 2026-06-16 — Canonical lineage = Python (step1_v2/step2/step3), not legacy Perl

**Context.** Two lineages exist: an older Perl lineage (`step1.pl`, `step3a/b/c.pl` with
`partition_barcodes.pl` + an SGE array job for Hamming-distance barcode pairing) and a
newer Python lineage. The legacy `step1/2/3.qsub` scripts are almost entirely commented
out; the master qsub runs the Python lineage.

**Decision.** Treat the **Python lineage** as canonical for dictionary building and build
the decomposed step scripts around it. Leave the Perl scripts in the old script dir for
reference only.

**Status / pointer.** Applied. Perl lineage noted as superseded in the methods file.

---

## 2026-06-16 — bc2 via range mode, not the r1_barcode2 library

**Context.** `step1_v2.py` can delimit bc2 either from a fixed-length `r1_barcode2`
library (matched by gene/block, case-sensitive uppercase gene) or in "range mode"
(`--bc2_min_len/--bc2_max_len` + `--r1_post_bc2` anchor). The available
`r1_barcode2.csv` (121725) only covers COMMD genes, RAD51C, XRCC3 — not RAD51D or XRCC2,
which appear in current batches.

**Decision.** Use **range mode** (bc2 2–7 nt, anchor `ACCTTTCGACCT`), matching the active
master.qsub. This sidesteps the incomplete barcode2 library and works for all genes.

**Evidence.** `r1_barcode2.csv` gene column = {COMMD1, COMMD5, COMMD6, COMMD10, RAD51C,
XRCC3}; current Ampseq pilot also has RAD51D and XRCC2.

**Status / pointer.** Applied in `config.sh` (BC2_MIN_LEN=2, BC2_MAX_LEN=7,
R1_POST_BC2=ACCTTTCGACCT); no `--r1_barcode2_lib` passed.

---

## 2026-06-16 — Ampseq pilot raw filenames renamed to convention

**Context.** Pilot raw files arrived as `rad51c_block1_dictionary_R1.fastq.gz` etc.,
which `step1_v2.py`'s parser `([^_]+)_DMS_(\d+)_Dictionary` cannot read (lowercase gene,
`block1` instead of `DMS_1`, lowercase `dictionary`), and the library match is
case-sensitive on uppercase gene.

**Decision.** Renamed all 18 files to `<GENE>_DMS_<block>_Dictionary_R{1,2}.fastq.gz`
(uppercase gene, `DMS_<n>`, `Dictionary`) and moved them under
`raw_data/Dual_NEBuilder/` (this batch is the Dual_NEBuilder construct type). A 1000-read
uncompressed peek for RAD51C_DMS_1 is in `raw_data/Dual_NEBuilder/peek_first1000/`.

**Status / pointer.** Applied. Files verified to parse to correct gene/block.

---

## 2026-06-17 — Committed to bc1-only keying (--ignore_bc2) going forward

**Context.** Evaluated three keyings for the Dual_NEBuilder dictionary: bc1+bc2,
the standard grid sweep, and bc1-only. Also ran a bc1-only, support=3 sweep over low
accuracy thresholds (0.1–0.5); even at accuracy 0.1 only ~36% of bc1 barcodes are
retained (each bc1 is smeared across many variants), and retention collapses to <0.3%
by 0.5 — no usable plateau. bc1-only is the intended Dual_NEBuilder design and remains
the path forward; low yield is expected from undersequencing in this pilot.

**Decision.** bc1-only is committed. Centralized the flag as `IGNORE_BC=--ignore_bc2`
in `config.sh`; steps 05/07/08 now reference `$IGNORE_BC` instead of per-script defaults.
The bc1+bc2 and grid sweep result CSVs were deleted; the surviving sweep is
`stats/dict_param_sweep_bc1.csv` (step 07 now writes that name). Set `IGNORE_BC=""`
only to revisit bc1+bc2.

**Status / pointer.** Applied in config.sh + 05/07/08 qsubs. Support/accuracy thresholds
themselves still being explored (BLOCK_LEVEL_ONLY=1 stays).
