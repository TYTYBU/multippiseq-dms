# Analysis Plan — cDNA edit count (part 2)

## Status / orientation   (read first)

*Load this file to bring a new session up to speed; this block is the current state, the rest is the detailed method.*

- **Canonical pipeline:** `multiPPIseq_DMS_central/qsub/cdna_edit_count/` (`config.sh`,
  `01_edit_count.qsub`, `02_cDNA_run_stats.qsub`) driving
  `multiPPIseq_DMS_central/script/cdna_edit_count/` (`cDNA_parse_step1.py`,
  `cDNA_parse_step2_allPos.py`, `extract_cDNA_run_stats.py`). Applies a barcode→variant
  dictionary (part 1 output) to cDNA reads and counts base edits per variant.
- **Output:** per batch, a **separate** pilot folder `data/<batch>/<OUT_SUBDIR>/`
  (default `cDNA_parsed_dictPilot062626/`) with one `<rep>.cDNA_parsed.csv` per rep
  (variant × edit-position count table) plus a `cDNA_diagnosis_stats_*.csv`.
  The canonical `cDNA_parsed/` and any existing stats files are left untouched.
- **Locked parameters:** bc1-only (`--ignore_bc2`); `EDIT_MODE=T2C`;
  `TEMPLATE_SEQ=TATCATAC…TGTAGTACC`; `R1_POST_BC2=ACCTTTCGACCTGCA`;
  `REGEX_MATCH_LEN=8`; `BC1` length 12–15. Override via env before `qsub` (see `config.sh`).
- **R2 contamination screen:** `R2_FILTER=0` by default (R1-only, the historical behaviour).
  Set `R2_FILTER=1` only for batches contaminated with dictionary reads — currently
  **`multiPPIseq_DMS_editRate_071526` only**. See step 01 below.
- **Primary edit position:** `PRIMARY_POS=14` (1-based, in the edit template) — the on-target
  base-editing site, surfaced per variant as `Primary_pos`/`Primary_edit_ct` and highlighted in
  the per-sample plot. Verified: across all 32 cDNA parses the auto-detected primary was 14 in
  **100%** of samples (44–62% of edits vs ≤13% runner-up). **`PRIMARY_POS` is tied to the current
  `TEMPLATE_SEQ`** — if the template changes, re-derive it (run once with `PRIMARY_POS=` empty and
  read the logged auto-detected value) before re-pinning.
- **Dictionary:** the **concatenated** dictionary built by dict_build step 12
  (`barcode_to_variant/<LIBTYPE>/concatenated/<LABEL>.barcode_to_variant.csv.gz`) — per-block
  dicts pooled, NOT the collapsed *combined* dict. Step 12 emits both a **per-gene** and a
  **per-group** label (`COMMD_all`, `RAD51C_XRCC3`, `RAD51D_XRCC2`); pick the label matching
  the cDNA sample's pooled content. A bc1 reused across blocks for different variants is
  counted toward **every** such variant (15 such bc1 for RAD51C_XRCC3 in the pilot).
  Sanity check when choosing a file: a concatenated dict has **no** `*.raw_dict.csv.gz`
  sibling; if one is present you are looking at a combined dict.
- **Other outputs:** Step1 bc-alt temp → `…dict_Ampseq_pilot_062626/temp/`; per-rep
  `Report_line` files → `data/<batch>/<OUT_SUBDIR>/reports/`.
- **Decision history & rationale:** `cdna_edit_count_decisions.md` (local-only, not in the repo).
- **Input integrity:** cDNA fastqs are **not** assumed intact — two files in `editRate_041726`
  were silently corrupt in transit and one killed a rep for months without failing the job
  (fixed; see decisions log 2026-07-22). `gzip -t` a batch before parsing it, and note that
  step 01 now exits non-zero when any rep dies.
- **Status:** all 24 biological cDNA samples reprocessed 2026-07-27 with the updated pipeline
  (primary position + plots). Canonical (sample × dictionary) parses, per the committed mapping:
  - **regular `RAD51C_XRCC3`** — 041726 (HCT+U2OS) + 051526 (HCT; no U2OS in that batch) →
    `cDNA_parsed/regular_primaryPos/` (12 samples).
  - **NEBuilder `RAD51C_XRCC3`** — 041726 + 051526, **HCT only** →
    `cDNA_parsed/Dual_NEBuilder_primaryPos/` (8 samples).
  - **NEBuilder `RAD51D_XRCC2`** — 071526, all X2RD samples, `R2_FILTER=1` →
    `cDNA_parsed/Dual_NEBuilder_r2filtered/` (12 samples). 071526 contamination is **per-sample**:
    PPIseq_regular rep1–3 and anti-myc_regular rep1 are 45–56% dictionary reads; the other 8 are
    ≤0.5%. Use the `_r2filtered/` tables; the unfiltered `Dual_NEBuilder/` is the baseline only.
  - RAD51C_XRCC3 cDNA is **regular** libtype (regular dict matches better: 24% vs 18% barcode
    match); it is parsed against both dicts as a comparison, U2OS against regular only.
- **Last updated:** 2026-07-27

---

## Pipeline steps

### Step 01 — Count edits in cDNA samples (`01_edit_count.qsub`)
- **Inputs:** the concatenated dictionary
  `barcode_to_variant/<LIBTYPE>/concatenated/<SAMPLE>.barcode_to_variant.csv.gz` (dict_build
  step 12); cDNA R1 fastqs `data/<batch>/raw_data/*<SAMPLE_GLOB>*_R1.fastq.gz`. With
  `R2_FILTER=1`, the matching `_R2.fastq.gz` and the block fixed-seq table
  `ref/dict_build/r2_fixed_seq_061626.csv` as well.
- **Method:**
  1. **Step1** (`cDNA_parse_step1.py`) — expand every dictionary barcode into all edited
     (T→C) alt-barcode forms → `…dict_Ampseq_pilot_062626/temp/<SAMPLE>.bc1_alts.txt`
     (+ empty bc2 alts).
  2. **Step2** (`cDNA_parse_step2_allPos.py`) — for each batch×rep: anchor the template
     (start/end regex of `REGEX_MATCH_LEN`), require extracted cDNA == template length, read
     bc1 upstream of the template and call edits at each ref-base (T) position; match the read's
     bc1 (via the alt table) to the dictionary and credit the read to **all** variants that bc1
     maps to. Writes `<rep>.cDNA_parsed.csv`; prints one `Report_line` to stdout, captured under
     `data/<batch>/<OUT_SUBDIR>/reports/<rep>.report.txt`.
  3. **R2 screen** (`R2_FILTER=1` only) — read R2 in lockstep with R1 and drop the pair unless
     R2 is R1's reverse mate. A pair is kept when R2 (a) contains **no** dictionary block fixed
     sequence and (b) starts (offset 0–2) with the reverse complement of `R1_POST_BC2`. That
     anchor is **derived in-script** from `R1_POST_BC2` + `EDIT_MODE`, wildcarding the positions
     an edit can change on the reverse strand (T2C → `rc(ACCTTTCGACCTGCA)` = `TGCAGGTCGAAAGGT`
     → `TGC.GGTCG...GGT`); a literal RC would match only *unedited* reads and bias the edit rate
     down. Rejects are counted separately as `dict` (contamination) vs `no_anchor` (adapter
     dimer / short insert), so the ~5–11% dimer background is never reported as contamination.
     `Total reads` still counts every pair, so the screen is visible in the stats rather than
     hidden. Mismatched/truncated R1–R2 pairs abort the rep.
  4. **Per-sample plot** (`plot_edit_positions.py`) — after every rep succeeds, one QC bar plot
     per parsed CSV: total edit events per editable template position (summed over variants),
     with `Primary_pos` highlighted → `data/<batch>/<OUT_SUBDIR>/plot/<rep>.edit_positions.png`.
     Runs as a post-pass over the CSVs, so a plot folder never reflects a partially-failed batch.
- **Failure semantics:** reps run concurrently; each PID is waited on individually and the job
  exits non-zero if any rep fails. A bare `wait` returns 0 even when a child died, which once let
  a corrupt-fastq rep fail silently while step 02 summarized only the survivors.
- **Output:** `data/<batch>/<OUT_SUBDIR>/<rep>.cDNA_parsed.csv` — columns
  `Outcome  Total_reads  Is_edited  Any_edit_ct  Primary_pos  Primary_edit_ct  Pos<N>…` (one row
  per variant, one Pos column per editable template position). `Primary_edit_ct` is the read count
  edited at `Primary_pos` — i.e. it equals that row's `Pos<Primary_pos>` cell, surfaced up front.
  `Primary_pos` is constant within a file (from `PRIMARY_POS`, or auto-detected when unset).
- **Rationale:** edits are scored at every template T position so per-position editing can be
  profiled; bc1-only matching mirrors the committed dictionary keying.

### Step 02 — cDNA-sample diagnosis stats (`02_cDNA_run_stats.qsub`)
- **Inputs:** the `Report_line` files from step 01 (`data/<batch>/<OUT_SUBDIR>/reports/`).
- **Method:** `extract_cDNA_run_stats.py` summarizes per rep; written **per batch** (rep
  basenames repeat across batches, so one CSV each avoids index collisions). Chain after
  step 01 with `qsub -hold_jid <01_jobid>`.
- **Output:** `data/<batch>/<OUT_SUBDIR>/cDNA_diagnosis_stats_dictPilot062626.csv`:
  Dictionary type, Dictionary name, Total reads, Reads with barcode matched, Reads with template
  matched, Reads with barcode and template matched, Unique barcodes in cDNA samples, Unique
  barcodes matched in dictionary. Two further columns — Read pairs screened on R2, Read pairs
  rejected as dictionary contamination — appear **only** when the batch was run with
  `R2_FILTER=1`, so unscreened batches' CSVs keep their original shape.
  Dictionary type/name are read off `DICT_FILE`'s path (not `LIBTYPE`/`SAMPLE`) so the recorded
  columns cannot drift from the dictionary step 01 actually used.
- **Rationale:** QC of *applying* the dictionary to cDNA — distinct from the dictionary-**build**
  diagnostics (dict_build steps 06–09).

## Notes
- These are cDNA-**sample** diagnostics; do not conflate with dictionary-**build** diagnostics.
- The pilot deliberately writes to a separate `OUT_SUBDIR` so the canonical `cDNA_parsed/`
  (built from a different dictionary) is never overwritten.
