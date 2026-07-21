# =====================================================================
# Shared config for the multiPPIseq DMS dictionary-building pipeline.
# Sourced by every NN_*.qsub step script. Edit batch-specific values
# here (or export them before `qsub`) so the step scripts stay generic.
# =====================================================================

PROJ_PATH=/net/bgm/sherwood/PPIseq_proj
CENTRAL=$PROJ_PATH/multiPPIseq_DMS_central
SCRIPT_PATH=$CENTRAL/script/dict_build

# ---- batch-specific (override via env before qsub if needed) ----------
DATA_PATH=${DATA_PATH:-$PROJ_PATH/data/multiPPIseq_DMS_dict_Ampseq_pilot_062626}
LIBTYPE=${LIBTYPE:-Dual_NEBuilder}          # Dual_NEBuilder | regular
GENES=${GENES:-"RAD51C RAD51D XRCC2 XRCC3"} # genes present in this batch

# ---- shared reference libraries (121725 = canonical reference set) ----
REF_PATH=$PROJ_PATH/data/multiPPIseq_DMS_dict_121725
# r2 fixed-seq lib: dated copy in central, extended with RAD51D blk7/8 + XRCC2 blk6/7
# (verified present in 94.6-98.3% of the corresponding R2 reads, 2026-06-16).
R2_FIXED_SEQ_LIB=$CENTRAL/ref/dict_build/r2_fixed_seq_061626.csv
# block ref fasta + bowtie2 index: dated central copies, extended with RAD51D blk7/8 +
# XRCC2 blk6/7 (R2 fixed-seq + 5' prefix stripped; samtools faidx + bowtie2-build done).
BLOCK_REF_FA=$CENTRAL/ref/dict_build/block_ref_seqs_061626.fa
BOWTIE2_INDEX=$CENTRAL/ref/dict_build/bowtie2_index/block_ref_seqs_061626
# variant lib: dated central copy, extended with RAD51D blk7/8 + XRCC2 blk6/7 (445 variants
# from the validation library xlsx; each variant_seq = oligo + constant 3' ref tail, padded
# to block-ref length so step2's var[:len(ref)] compare works; ~98.7% identity to ref).
VARIANT_LIB=$CENTRAL/ref/dict_build/variant_sequence_library_lite_061626.csv

# ---- per-step I/O paths (all namespaced by library type) --------------
# RAW_PATH is overridable so a batch can be run off a merged read set
# (e.g. merged_data/<label>/) instead of its own raw_data/<libtype>.
RAW_PATH=${RAW_PATH:-$DATA_PATH/raw_data/$LIBTYPE}
FASTQC_PATH=$DATA_PATH/fastqc/$LIBTYPE
NOADAPTER_PATH=$DATA_PATH/noAdapter/$LIBTYPE
PROCESSED_FA_PATH=$DATA_PATH/processed_fa/$LIBTYPE
BAM_PATH=$DATA_PATH/bowtie2_out/$LIBTYPE
VARIANT_MAPPED_PATH=$DATA_PATH/variant_mapped/$LIBTYPE
BC2VAR_PATH=$DATA_PATH/barcode_to_variant/$LIBTYPE
STATS_PATH=$BC2VAR_PATH/stats                 # diagnostics: steps 06-08

# ---- adapter handling -------------------------------------------------
# Set USE_NOADAPTER=1 once 01_cutadapt has been run; step 02 will then read
# trimmed reads from NOADAPTER_PATH instead of raw reads from RAW_PATH.
# Whether cutadapt is needed at all is decided from the 00_fastqc report.
USE_NOADAPTER=${USE_NOADAPTER:-0}
CUTADAPT_ADAPTER=${CUTADAPT_ADAPTER:-AGATCGGAAGAG}   # 3' Illumina TruSeq

# ---- tool environments ------------------------------------------------
# Two conda envs cover the whole pipeline (no dedicated env exists / this host
# has no outbound internet to build one). Each step script activates the right
# one via:  source "$CONDA_BASE/bin/activate" "$ENV_<X>"
#   ENV_MAIN     (access)   -> fastqc, bowtie2, samtools, python+pandas+pysam
#                             (steps 00, 02, 03, 04, 05)
#   ENV_CUTADAPT (cutadapt) -> cutadapt 4.6 (step 01 only)
CONDA_BASE=/net/apps/conda/miniconda3
ENV_MAIN=access
ENV_CUTADAPT=cutadapt
NUM_THREADS=${NUM_THREADS:-32}

# ---- step 02 (barcode/oligo extraction; step1_v2.py) params -----------
R1_LINKER=TATCATACggccctgcagacatgggtgatcctcatgttggccTGTAGTACC
BC1_MIN_LEN=12
BC1_MAX_LEN=15
BC2_MIN_LEN=2
BC2_MAX_LEN=7
R1_POST_BC2=ACCTTTCGACCT

# ---- step 05 (barcode->variant dictionary; step3.py) params -----------
MIN_SUPPORT_READS=3
MIN_ACCURACY=0.9
# Optional per-block thresholds. If this CSV exists, step 05 reads each block's
# support/fidelity cutoff from it (columns: Gene,Block_number,
# proposed_fidelity_cutoff,proposed_support_read_cutoff) instead of the global
# MIN_SUPPORT_READS/MIN_ACCURACY above. Tentative pilot values; may be updated.
DICT_PARAMS_CSV=${DICT_PARAMS_CSV:-$DATA_PATH/dictionary_parameters.csv}
# Barcode keying: committed to bc1 only for Dual_NEBuilder (--ignore_bc2). The
# bc1+bc2 and grid sweeps were evaluated and dropped; bc1-only is the design
# going forward. Used by steps 05, 07, 08. Set to "" only to revisit bc1+bc2.
IGNORE_BC=${IGNORE_BC:---ignore_bc2}
# Dual_NEBuilder: step3 parameters are still being explored, so for now build
# only the per-block dictionaries and skip the per-gene merge (set to 0 to also
# build per-gene dictionaries once parameters are finalized).
BLOCK_LEVEL_ONLY=${BLOCK_LEVEL_ONLY:-1}

# ---- steps 10-11 (combined dictionaries) params -----------------------
# The combined dictionaries (pool reads across blocks AND across the genes in
# each group) are OPTIONAL: the cDNA pipeline uses the per-block dictionaries.
# Steps 10 (combine variant-mapped), the combined scatter, and 11 (combined
# dictionary build) only run when BUILD_COMBINED=1.
BUILD_COMBINED=${BUILD_COMBINED:-0}
# Group definitions, by convention: all COMMD* together; RAD51D+XRCC2 together;
# RAD51C+XRCC3 together. Format: "<label>:<gene,gene,...>" entries,
# space-separated. Groups with no matching blocks in the batch are skipped.
GENE_GROUPS=${GENE_GROUPS:-"RAD51C_XRCC3:RAD51C,XRCC3 RAD51D_XRCC2:RAD51D,XRCC2 COMMD_all:COMMD1,COMMD5,COMMD6,COMMD10"}
# Per-group support/fidelity cutoffs for the combined dictionaries
# (columns: Group,proposed_fidelity_cutoff,proposed_support_read_cutoff).
# Falls back to global MIN_SUPPORT_READS/MIN_ACCURACY if a group is absent.
DICT_PARAMS_COMBINED_CSV=${DICT_PARAMS_COMBINED_CSV:-$DATA_PATH/dictionary_parameters_combined.csv}
