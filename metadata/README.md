# metadata — sample inventory

## `DMS_fastq_summary.csv`

One row per fastq **sample** (R1/R2 pair) across every multiPPIseq DMS project folder —
the provenance index for "which reads went into which analysis." Moved here from
`data/DMS_fastq_summary.csv` on 2026-07-21 so it sits with the pipeline rather than among
the reads and outputs.

> **The CSV itself is deliberately untracked** (see `.gitignore`) — it is a data record, not
> code. This README is tracked and describes it; the file lives on disk only. Back it up
> with the data, not with the repo, and do not expect `git` to restore it.

### Columns

| Column | Meaning |
|---|---|
| `project_folder` | folder under `PPIseq_proj/data/` |
| `project_type` | `DMS_dict` (dictionary building) or `DMS_editRate` (cDNA edit counting) |
| `sample_name` | fastq basename minus `_R{1,2}.fastq.gz` |
| `library_type` | `regular` or `Dual_NEBuilder` — the construct/prep type |
| `fastq_location` | path **relative to `project_folder`** |
| `is_merged` | `yes` if these reads are a concatenation of several batches |
| `merged_from` | **per-block** source batches (not the folder label — a folder named for four batches can still hold blocks that exist in only one) |
| `NGS batch` | delivery batch(es), matching `merged_from` for merged rows |
| `Original sample name` | human-readable name as delivered, e.g. `RAD51C Block 3 Dictionary` |
| `notes` | provenance caveats — added 2026-07-21, empty for pre-existing rows |

### Row order and what is kept

- **Sorted newest batch first.** Batch date is parsed from the `MMDDYY` suffix of
  `project_folder` (`071526` = 2026-07-15).
- **Current-year batches on top; everything older sits at the bottom** and is marked
  `RECORD ONLY` in `notes`. Those rows are kept for provenance only — they are not part of
  the current dictionary-building scheme and should not be fed to the pipeline.
- **Only the most recent merged set of a lineage is kept.** Each new merge supersedes the
  last because it is strictly more complete, so the older ones are removed rather than left
  to be picked up by mistake. Dropped on 2026-07-21:
  `dict_041726/032426_041726_merged/regular` (13 rows) and
  `dict_051526/032426_041726_051526_merged/regular` (4 rows), both superseded by
  `dict_071526/merged_data/regular_032426_041726_051526_071526`, verified a strict superset
  of their block sets beforehand.
  This applies to **successive merges of the same lineage**, not to unrelated experiments
  that happen to be merged — e.g. the `DMS_editRate` merged samples are distinct experiments
  and are all retained.

### Reading it correctly

- **`merged_from` is per block, deliberately.** In
  `merged_data/regular_032426_041726_051526_071526`, `COMMD5_DMS_1` came from 032426
  alone while `RAD51C_DMS_3` came from three batches. Do not infer a block's depth from
  the folder name.
- **Check `notes` before trusting a merged folder.** The legacy
  `032426_041726_051526_merged/regular` under `dict_051526` was built by an
  intersection-only merge that silently dropped all 9 COMMD blocks; it is superseded.
- **`library_type` for the 071526 batch was corrected on 2026-07-21** — the delivery had
  regular and Dual_NEBuilder swapped. Verify libtype from read content (R2 fixed-sequence
  match rate), never from a delivered filename.

### Maintenance

Append rows when a batch is delivered or a merged read set is created, then re-sort so the
new batch is on top. Batch dates are `MMDDYY`.

When a new merge supersedes an existing one, **verify the new set is a superset of the old
block set before dropping the old rows** — a merge script that intersects rather than unions
will silently produce a smaller set (this is exactly what happened to the COMMD blocks in
`032426_041726_051526_merged`).

At year rollover the `RECORD ONLY` boundary moves: the batches currently on top become
record-only and sink to the bottom. That is a deliberate manual step, not something the file
does on its own — re-tag and re-sort when it happens.
