# =====================================================================
# Shared config for the multiPPIseq DMS cDNA edit-count pipeline (part 2:
# apply a barcode->variant dictionary to cDNA samples and count edits).
# Sourced by every NN_*.qsub step script.
# =====================================================================

PROJ_PATH=/net/bgm/sherwood/PPIseq_proj
CENTRAL=$PROJ_PATH/multiPPIseq_DMS_central
SCRIPT_PATH=$CENTRAL/script/cdna_edit_count

# ---- tool environment (access: python+pandas) -------------------------
CONDA_BASE=/net/apps/conda/miniconda3
ENV_MAIN=access

# ---- dictionary source (PILOT per-block dictionaries) -----------------
# This pilot uses the per-block dictionaries (NOT the combined dictionary):
# the pooled RAD51C_XRCC3 cDNA sample's dictionary is the concatenation of the
# RAD51C and XRCC3 per-block dicts. A bc1 reused across blocks for different
# variants is kept (cDNA_parse_step2 counts a read toward every such variant).
# The pooled dictionary is the CONCATENATED dictionary built by dict_build
# step 12 (barcode_to_variant/<LIBTYPE>/concatenated/<SAMPLE>...).
LIBTYPE=${LIBTYPE:-Dual_NEBuilder}
DICT_PILOT=${DICT_PILOT:-$PROJ_PATH/data/multiPPIseq_DMS_dict_Ampseq_pilot_062626}
SAMPLE=${SAMPLE:-RAD51C_XRCC3}
DICT_FILE=${DICT_FILE:-$DICT_PILOT/barcode_to_variant/$LIBTYPE/concatenated/${SAMPLE}.barcode_to_variant.csv.gz}

# ---- bc-alt temp (Step1 output; under the pilot dict folder) ----------
TEMP_DIR=${TEMP_DIR:-$DICT_PILOT/temp}

# ---- cDNA inputs / outputs -------------------------------------------
# One or more editRate batches; for each, parse the reps whose filename
# contains SAMPLE_GLOB. Output goes to a SEPARATE pilot folder per batch
# (OUT_SUBDIR) so the canonical cDNA_parsed/ is left untouched. The per-rep
# Report_line files (consumed by step 02) live under <batch>/<OUT_SUBDIR>/reports.
EDITRATE_BATCHES=${EDITRATE_BATCHES:-"multiPPIseq_DMS_editRate_041726 multiPPIseq_DMS_editRate_051526"}
SAMPLE_GLOB=${SAMPLE_GLOB:-HCT_RAD51C_XRCC3}
OUT_SUBDIR=${OUT_SUBDIR:-cDNA_parsed_dictPilot062626}
REPORT_SUBDIR=${REPORT_SUBDIR:-reports}
STATS_CSV_NAME=${STATS_CSV_NAME:-cDNA_diagnosis_stats_dictPilot062626.csv}

# ---- cDNA parse params (from the editRate legacy script) --------------
EDIT_MODE=${EDIT_MODE:-T2C}
TEMPLATE_SEQ=${TEMPLATE_SEQ:-TATCATACggccctgcagacatgggtgatcctcatgttggccTGTAGTACC}
REGEX_MATCH_LEN=${REGEX_MATCH_LEN:-8}
R1_POST_BC2=${R1_POST_BC2:-ACCTTTCGACCTGCA}
BC1_MIN_LEN=${BC1_MIN_LEN:-12}
BC1_MAX_LEN=${BC1_MAX_LEN:-15}
# Barcode keying: bc1-only (committed Dual_NEBuilder design); matches the dict.
IGNORE_BC=${IGNORE_BC:---ignore_bc2}
# ---- R2 dictionary-contamination screen (step 01) ---------------------
# Off by default. Set R2_FILTER=1 for batches contaminated with dictionary reads
# (known: multiPPIseq_DMS_editRate_071526). Step 01 then also reads R2 and keeps
# only pairs whose R2 is R1's reverse mate: it must start with the reverse
# complement of R1_POST_BC2 (derived in-script from R1_POST_BC2 + EDIT_MODE, with
# the editable positions wildcarded) and must not contain any dictionary block
# fixed sequence. Measured on gold-standard data 2026-07-22: 0.00% of true
# dictionary reads pass the screen, 0.00% of true cDNA reads are rejected by the
# fixed-seq clause.
R2_FILTER=${R2_FILTER:-0}
R2_FIXED_SEQ_LIB=${R2_FIXED_SEQ_LIB:-$CENTRAL/ref/dict_build/r2_fixed_seq_061626.csv}

# ---- primary (on-target) edit position (step 01) ----------------------
# 1-based position IN THE EDIT TEMPLATE, surfaced per variant as Primary_pos /
# Primary_edit_ct and highlighted in the per-sample edit-position plot. Pinned to 14:
# across all 32 cDNA parses (2 cell lines, both dict libtypes, 3 assays, 3 batches) the
# auto-detected primary was 14 in 100% of samples, carrying 44-62% of edits vs <=13% for
# the runner-up. Set PRIMARY_POS= (empty) to fall back to per-sample auto-detection
# (argmax of edits over all variants) -- convenient for a new template, but NOT stable
# across samples, so pin it for any multi-sample comparison.
# CAVEAT: 14 is tied to the CURRENT TEMPLATE_SEQ. If TEMPLATE_SEQ changes, both the
# coordinate system and which position is on-target change -- re-derive PRIMARY_POS
# (run once with PRIMARY_POS= empty and read the auto-detected value) before re-pinning.
PRIMARY_POS=${PRIMARY_POS:-14}
