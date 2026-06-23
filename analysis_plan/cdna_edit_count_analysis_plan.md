# Analysis Plan — cDNA edit count (part 2)

## Status / orientation   (read first)

*Load this file to bring a new session up to speed; this block is the current state, the rest is the detailed method.*

- **Canonical pipeline:** `multiPPIseq_DMS_central/qsub/cdna_edit_count/` (`config.sh`,
  `01_edit_count.qsub`, `02_cDNA_run_stats.qsub`) driving
  `multiPPIseq_DMS_central/script/cdna_edit_count/` (`cDNA_parse_step1.py`,
  `cDNA_parse_step2_allPos.py`, `extract_cDNA_run_stats.py`). Applies a barcode→variant
  dictionary (part 1 output) to cDNA reads and counts base edits per variant.
- **Output:** per batch, a **separate** pilot folder `data/<batch>/<OUT_SUBDIR>/`
  (default `cDNA_parsed_dictPilot062626/`) with one `<rep>.cDNA_parsed.txt` per rep
  (variant × edit-position count table) plus a `cDNA_diagnosis_stats_*.csv`/`.xlsx`.
  The canonical `cDNA_parsed/` and any existing stats files are left untouched.
- **Locked parameters:** bc1-only (`--ignore_bc2`); `EDIT_MODE=T2C`;
  `TEMPLATE_SEQ=TATCATAC…TGTAGTACC`; `R1_POST_BC2=ACCTTTCGACCTGCA`;
  `REGEX_MATCH_LEN=8`; `BC1` length 12–15. Override via env before `qsub` (see `config.sh`).
- **Pilot dictionary:** the **per-block** dictionaries (NOT the combined dict). The pooled
  RAD51C_XRCC3 cDNA sample's dictionary is the concatenation of the RAD51C + XRCC3 per-block
  dicts; a bc1 reused across blocks for different variants is counted toward **every** such
  variant (15 such bc1 in the pilot).
- **Decision history & rationale:** `cdna_edit_count_decisions.md`.
- **Status:** pilot run complete for HCT_RAD51C_XRCC3 (editRate_041726 + editRate_051526).
- **Last updated:** 2026-06-23

---

## Pipeline steps

### Step 01 — Count edits in cDNA samples (`01_edit_count.qsub`)
- **Inputs:** per-block dicts under `…dict_Ampseq_pilot_062626/barcode_to_variant/<LIBTYPE>/blocks/`
  (globs in `SAMPLE_BLOCK_GLOBS`); cDNA R1 fastqs `data/<batch>/raw_data/*<SAMPLE_GLOB>*_R1.fastq.gz`.
- **Method:**
  1. **Pool dictionary** — concatenate the SAMPLE's per-block dicts (gzip multi-member concat)
     into `WORK_PATH/dict/<SAMPLE>.perblock_concat.barcode_to_variant.csv.gz`. Duplicate bc1
     rows (same bc1 across blocks) are kept on purpose.
  2. **Step1** (`cDNA_parse_step1.py`) — expand every dictionary barcode into all edited
     (T→C) alt-barcode forms → `WORK_PATH/temp/<SAMPLE>.bc1_alts.txt` (+ empty bc2 alts).
  3. **Step2** (`cDNA_parse_step2_allPos.py`) — for each batch×rep: anchor the template
     (start/end regex of `REGEX_MATCH_LEN`), require extracted cDNA == template length, read
     bc1 upstream of the template and call edits at each ref-base (T) position; match the read's
     bc1 (via the alt table) to the dictionary and credit the read to **all** variants that bc1
     maps to. Writes `<rep>.cDNA_parsed.txt`; prints one `Report_line` to stdout, captured under
     `WORK_PATH/reports/<batch>/<rep>.report.txt`.
- **Output:** `data/<batch>/<OUT_SUBDIR>/<rep>.cDNA_parsed.txt` — columns
  `Outcome  Total_reads  Is_edited  Any_edit_ct  Pos<N>…` (one row per variant, one Pos column
  per editable template position).
- **Rationale:** edits are scored at every template T position so per-position editing can be
  profiled; bc1-only matching mirrors the committed dictionary keying.

### Step 02 — cDNA-sample diagnosis stats (`02_cDNA_run_stats.qsub`)
- **Inputs:** the `Report_line` files from step 01 (`WORK_PATH/reports/<batch>/`).
- **Method:** `extract_cDNA_run_stats.py` summarizes per rep; written **per batch** (rep
  basenames repeat across batches, so one CSV each avoids index collisions). Chain after
  step 01 with `qsub -hold_jid <01_jobid>`.
- **Output:** `data/<batch>/<OUT_SUBDIR>/cDNA_diagnosis_stats_dictPilot062626.csv` (+ `.xlsx`):
  Total reads, Reads with barcode matched, Reads with template matched, Reads with barcode and
  template matched, Unique barcodes in cDNA samples, Unique barcodes matched in dictionary.
- **Rationale:** QC of *applying* the dictionary to cDNA — distinct from the dictionary-**build**
  diagnostics (dict_build steps 06–09).

## Notes
- These are cDNA-**sample** diagnostics; do not conflate with dictionary-**build** diagnostics.
- The pilot deliberately writes to a separate `OUT_SUBDIR` so the canonical `cDNA_parsed/`
  (built from a different dictionary) is never overwritten.
