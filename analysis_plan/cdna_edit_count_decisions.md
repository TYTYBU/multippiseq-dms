# Decisions log — multiPPIseq DMS cDNA edit count (part 2)

Dated, append-only. Newest entries on top. One entry per real decision.

---

## 2026-06-23 — Concatenated dictionary promoted to a dict_build artifact; temp/reports relocated

**Context.** The pooled per-block dictionary was originally built inline by cDNA step 01 into a
`cdna_edit_count/` work folder under the pilot dict dir. Tian wanted it to live alongside the
`combined/` dictionaries as a first-class artifact, and the scratch folder cleaned up.

**Decisions.**
1. **Concatenated dictionaries are now built by dict_build step 12**
   (`12_concatenate_dictionary.qsub`) into `barcode_to_variant/<LIBTYPE>/concatenated/`, one per
   `GENE_GROUPS` group, by default (cheap gzip concat; not gated by `BUILD_COMBINED`). Generated
   for `RAD51C_XRCC3` (47,961 rows) and `RAD51D_XRCC2` (44,168 rows); COMMD_all skipped (no
   blocks in batch).
2. **cDNA step 01 now consumes that prebuilt dict** (`DICT_FILE` → `concatenated/<SAMPLE>…`)
   instead of building it inline; errors out if the dict is missing.
3. **Relocations:** Step1 bc-alt temp → `…dict_Ampseq_pilot_062626/temp/`; per-rep Report_line
   files → `data/<batch>/<OUT_SUBDIR>/reports/` (reports are needed by step 02, so kept, but
   moved under the respective cDNA project folder). The `cdna_edit_count/` work folder was
   removed.

**Evidence / pointer.** `qsub/dict_build/12_concatenate_dictionary.qsub`; updated
`qsub/cdna_edit_count/{config.sh,01_edit_count.qsub,02_cDNA_run_stats.qsub}`. Re-ran the pilot;
parsed CSVs + stats unchanged.

---

## 2026-06-23 — Pilot cDNA edit count: per-block dict, multi-variant, separate output

**Context.** First run of part 2 on the pilot dictionary, for the HCT_RAD51C_XRCC3 cDNA
samples in `multiPPIseq_DMS_editRate_041726` and `_051526`. The legacy reference is
`qsub_script/multiPPIseq_DMS_editRate.qsub` (single per-gene-group dict, `--ignore_bc2`,
T2C). We had just built both per-block and (optional) combined pilot dictionaries.

**Decisions.**
1. **Dictionary = concatenated per-block dicts**, not the combined dict. The pooled
   RAD51C_XRCC3 sample needs one dict spanning RAD51C blk1–2 + XRCC3 blk1–3, so we `cat` the
   five per-block `barcode_to_variant.csv.gz` files. This honors the per-block support/fidelity
   cutoffs chosen in part 1.
2. **A bc1 mapping to multiple variants counts the read toward ALL of them.** Concatenating
   per-block dicts can put the same bc1 on multiple rows (15 such bc1 in the pilot,
   47,946 unique of 47,961 rows). The legacy Step2 built `barcode_pair_dict[(bc1,bc2)] = outcome`
   (last-wins under `--ignore_bc2`, where bc2=""), silently dropping the others. Fixed in the
   central copy `cDNA_parse_step2_allPos.py`: the pair maps to a **set** of outcomes and a
   matching read is credited to every variant in it.
3. **Separate pilot output folder.** Outputs go to `data/<batch>/cDNA_parsed_dictPilot062626/`
   (parsed txt + `cDNA_diagnosis_stats_dictPilot062626.csv/.xlsx`); the canonical `cDNA_parsed/`
   and pre-existing stats files are untouched.
4. **Stats per batch.** Rep basenames repeat across the two batches, so one diagnosis CSV is
   written per batch to avoid sample-index collisions.

**Evidence / pointer.** New module `multiPPIseq_DMS_central/{qsub,script}/cdna_edit_count/`.
Pilot result (barcode-and-template-matched reads): 041726 ≈ 1.79–2.31 M/rep; 051526
≈ 0.54–0.67 M/rep; ~27–30 k dictionary barcodes matched per rep. RAD51C_XRCC3 parsed tables
cover 1,192 variants.
