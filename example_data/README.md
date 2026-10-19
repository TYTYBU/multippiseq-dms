# Example data

Small, self-contained datasets for running the multiPPIseq DMS pipeline
end-to-end as a smoke test. Both halves of the pipeline are exercised:

1. **Dictionary build** — barcode extraction, alignment, variant assignment, and
   collapse to a barcode-to-variant dictionary.
2. **cDNA edit counting** — expansion of dictionary barcodes to their edited
   forms, then per-variant edit counting from a cDNA sample.

## The two library types

Reads come in two construct types, processed separately and never pooled across
types. One example is bundled for each. The pipeline **scripts are the same** for
both, and the `step1` barcode parameters are identical; what differs is the
`step3` filtering (and, in production, adapter trimming).

| | **Dual_NEBuilder** | **regular** |
|---|---|---|
| example block | RAD51C block 1 (residues 109-113) | RAD51C block 3 (residues 232-236) |
| step3 cutoffs | per-block: support >= 5, fidelity >= 0.10 | global: support >= 3, fidelity >= 0.90 |
| construct | low-fidelity ("messy"); many variants per barcode | clean; ~1 variant per barcode (fidelity ~1.0) |
| adapter trimming (production) | trimmed (cutadapt) | untrimmed |
| barcode keying | variant barcode only (`--ignore_bc2`) | variant barcode only (`--ignore_bc2`) |

Both examples use the **same cDNA sample** (`NR-TV_HCT_RAD51C_XRCC3_Rep1`), which
spans the whole RAD51C ORF, so it drives either block.

Each example is **sized for overlap**: the cDNA reads are those whose variant
barcode is present in the bundled dictionary, so edit counting finds real matches
instead of coming up empty on a random subsample. All FASTQs are ~7-8k reads.

## Layout

```
example_data/
  Dual_NEBuilder/
    dict_build/        RAD51C_DMS_1_Dictionary_R{1,2}.fastq.gz
    cdna_edit_count/   NR-TV_HCT_RAD51C_XRCC3_Rep1__R1.fastq.gz
    expected_output/   built dictionary + per-variant edit counts
  regular/
    dict_build/        RAD51C_DMS_3_Dictionary_R{1,2}.fastq.gz
    cdna_edit_count/   NR-TV_HCT_RAD51C_XRCC3_Rep1__R1.fastq.gz
    expected_output/   built dictionary + per-variant edit counts
```

The end-to-end walkthrough lives at the repository root:
`../full_pipeline.ipynb` runs both halves of the pipeline on this data,
calling the pipeline scripts directly (the same commands the `qsub/` jobs issue,
run inline rather than submitted to the cluster).

The cDNA half uses only R1: the `NR-TV_HCT_RAD51C_XRCC3` sample does not require
the R2 dictionary-contamination screen, so no cDNA R2 is bundled.

## Provenance

- **Dictionary reads** come from the `multiPPIseq_DMS_dict_071526` build
  (merged `RAD51C_DMS_1` Dual_NEBuilder reads, merged `RAD51C_DMS_3` regular reads).
- **cDNA reads** come from `multiPPIseq_DMS_editRate_041726`, sample
  `NR-TV_HCT_RAD51C_XRCC3_Rep1`.

These are the datasets paired in the production run that the methods describe.
The FASTQs were subset from the full production data by keeping, for a few
consecutive residues with strong editing signal, the dictionary read pairs for a
set of well-covered variant barcodes and the cDNA reads whose (possibly edited)
variant barcode matches one of them &mdash; so edit counting finds real matches
rather than coming up empty on a random subsample. Reproducing the subset needs
the lab-internal production data mounted and is not needed to run the pipeline.

## How to run

Open `../full_pipeline.ipynb` (at the repository root) with a kernel that
has `pysam`, `pandas`, `numpy` and `matplotlib` (e.g. the `access` conda env),
with `bowtie2` and `samtools` available. Pick the library type in the config cell
(`LIBTYPE = "Dual_NEBuilder"` or `"regular"`) and run all cells. The notebook
builds the bowtie2 index from the committed reference FASTA on first run (it is
not checked in), runs both halves of the pipeline by invoking the scripts
directly, plots the per-variant edit rates, and checks the result against that
library type's `expected_output/`.

## Parameters used (production settings)

Shared, both halves (`qsub/dict_build/config.sh`, `qsub/cdna_edit_count/config.sh`):

- R1 linker / cDNA edit template: `TATCATACggccctgcagacatgggtgatcctcatgttggccTGTAGTACC`
- variant barcode (bc1) length 12-15; `--ignore_bc2` (keyed on the variant barcode only)
- bowtie2: `--local --very-sensitive-local --no-unal`
- edit mode T2C; primary (on-target) edit position: template position 14
- `--r1_post_bc2 ACCTTTCGACCTGCA` (cDNA parse)

Differs by library type (step3): Dual_NEBuilder uses per-block cutoffs
(support >= 5, fidelity >= 0.10 for RAD51C block 1); regular uses the global
support >= 3, fidelity >= 0.90.

## Expected result

- **Dual_NEBuilder**: dictionary build yields **62 barcodes** across residues
  109-113; parsing the cDNA matches ~5,970 of 7,156 reads (~83%).
- **regular**: dictionary build yields **75 barcodes** across residues 232-236,
  all at fidelity ~1.0 (clean construct); parsing the cDNA matches essentially
  all ~7,280 reads.

In both, the primary edit position (template position 14) is the most edited in
nearly every variant, and the per-variant edit rate is
`Primary_edit_ct / Total_reads`.
