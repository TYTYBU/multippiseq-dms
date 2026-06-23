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
LIBTYPE=${LIBTYPE:-Dual_NEBuilder}
DICT_PILOT=${DICT_PILOT:-$PROJ_PATH/data/multiPPIseq_DMS_dict_Ampseq_pilot_062626}
BLOCKS_DICT_DIR=$DICT_PILOT/barcode_to_variant/$LIBTYPE/blocks
# per-block dict file globs that compose this sample's pooled dictionary
SAMPLE=${SAMPLE:-RAD51C_XRCC3}
SAMPLE_BLOCK_GLOBS=${SAMPLE_BLOCK_GLOBS:-"RAD51C_DMS_* XRCC3_DMS_*"}

# ---- work area (concat dict + bc-alt temp; under the pilot dict folder) -
WORK_PATH=${WORK_PATH:-$DICT_PILOT/cdna_edit_count}
DICT_OUT_DIR=$WORK_PATH/dict
TEMP_DIR=$WORK_PATH/temp
REPORT_DIR=$WORK_PATH/reports
DICT_FILE=$DICT_OUT_DIR/${SAMPLE}.perblock_concat.barcode_to_variant.csv.gz

# ---- cDNA inputs / outputs -------------------------------------------
# One or more editRate batches; for each, parse the reps whose filename
# contains SAMPLE_GLOB. Output goes to a SEPARATE pilot folder per batch
# (OUT_SUBDIR) so the canonical cDNA_parsed/ is left untouched.
EDITRATE_BATCHES=${EDITRATE_BATCHES:-"multiPPIseq_DMS_editRate_041726 multiPPIseq_DMS_editRate_051526"}
SAMPLE_GLOB=${SAMPLE_GLOB:-HCT_RAD51C_XRCC3}
OUT_SUBDIR=${OUT_SUBDIR:-cDNA_parsed_dictPilot062626}
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
