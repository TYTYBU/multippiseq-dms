# multiPPIseq_DMS_central

Central workspace for the **multiPPIseq DMS** pipeline — scripts, qsub job
scripts, and the analysis plan / lab notebook. The pipeline has two parts:

1. **Dictionary building** — from raw amplicon reads, build a
   barcode → variant dictionary for each gene/block. *(in progress)*
2. **cDNA edit count** — use the dictionaries to count edits/variants in cDNA
   readout. *(to be filled in)*

## Layout

```
multiPPIseq_DMS_central/
  README.md                       # this file
  analysis_plan/                  # lab notebook (see analysis_plan/README.md)
    README.md
    dict_build_analysis_plan.md   # methods + onboarding for dictionary building
    dict_build_decisions.md       # dated decision log
  script/
    dict_build/                   # canonical pipeline scripts (copied here)
    cdna_edit_count/              # (to be filled in)
  qsub/
    dict_build/                   # step-wise SGE job scripts + config.sh
    cdna_edit_count/              # (to be filled in)
```

## Running the dictionary build

Edit `qsub/dict_build/config.sh` for the batch (DATA_PATH, LIBTYPE, GENES),
then submit the steps in order:

| Step | Script | What |
|------|--------|------|
| 00 | `00_fastqc.qsub` | FastQC on raw reads → inspect for adapter read-through |
| 01 | `01_cutadapt.qsub` | *(optional)* adapter trim; set `USE_NOADAPTER=1` after |
| 02 | `02_extract_barcodes.qsub` | extract bc1/bc2/oligo → FASTA (`step1_v2.py`) |
| 03 | `03_bowtie2_align.qsub` | align oligo to block refs → BAM |
| 04 | `04_variant_map.qsub` | map reads to variants (`step2.py`) |
| 05 | `05_barcode_to_variant.qsub` | collapse to barcode→variant dict (`step3.py`) |

```bash
qsub qsub/dict_build/00_fastqc.qsub      # then inspect fastqc/, decide on 01
qsub qsub/dict_build/02_extract_barcodes.qsub
...
```

Start with `analysis_plan/dict_build_analysis_plan.md` (its "Status / orientation"
block) to come up to speed.
