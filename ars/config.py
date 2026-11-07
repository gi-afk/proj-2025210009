"""Project-wide configuration: paths, constants and download endpoints."""

from pathlib import Path

# Repository root = two levels up from ars/config.py
ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
SAMPLE_DIR = DATA_DIR / "sample"
RESULTS_DIR = ROOT / "results"
FIGURE_DIR = RESULTS_DIR / "figures"

# ---------------------------------------------------------------------------
# Reference genome: Saccharomyces cerevisiae S288C, assembly R64 (sacCer3)
# NCBI RefSeq GCF_000146045.2. Coordinates match SGD / OriDB (R64).
# ---------------------------------------------------------------------------
GENOME_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/"
    "GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz"
)
GENOME_FILE = RAW_DIR / "GCF_000146045.2_R64_genomic.fna.gz"

ASSEMBLY_REPORT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/"
    "GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_report.txt"
)
ASSEMBLY_REPORT_FILE = RAW_DIR / "GCF_000146045.2_R64_assembly_report.txt"

# ---------------------------------------------------------------------------
# Experimental ARS coordinates: OriDB "all replication origins" table.
# Only rows with status == "Confirmed" are used as positives.
# ---------------------------------------------------------------------------
ORIDB_BASE = "http://cerevisiae.oridb.org/data_output.php"
ORIDB_PARAMS = {"main": "sc_ori", "table": "sc_ori", "ext_format": "", "format": "tab"}
ORIDB_FILE = RAW_DIR / "oridb_sc_ori.tsv"

# ---------------------------------------------------------------------------
# Landmark annotations (centromeres) from SGD. Telomeres are represented by
# the chromosome ends (position 0 and chromosome length).
# ---------------------------------------------------------------------------
SGD_FEATURES_URL = (
    "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
)
SGD_FEATURES_FILE = RAW_DIR / "SGD_features.tab"

# ---------------------------------------------------------------------------
# Dataset construction
# ---------------------------------------------------------------------------
WINDOW = 500              # bp, centred on the annotated ARS midpoint
MIN_POS_DISTANCE = 500    # merge positive midpoints closer than this
POS_EXCLUSION = 1000      # no negative window within this distance of a positive
SEED = 2024
N_BOOTSTRAP = 1000
N_FOLDS_GROUPS = 16       # leave-one-chromosome-out

NEG_RANDOM = "random"
NEG_MATCHED = "distance_matched"
NEG_TYPES = (NEG_RANDOM, NEG_MATCHED)

FEATURES = [
    "gc_content",
    "at_content",
    "acs_pwm_max",
    "kmer3_entropy",
    "dist_cen_tel_log",
]

FEATURE_LABELS = {
    "gc_content": "GC content",
    "at_content": "AT content",
    "acs_pwm_max": "ACS PWM score (bits)",
    "kmer3_entropy": "3-mer entropy (norm.)",
    "dist_cen_tel_log": "log1p distance to nearest\ncentromere / telomere",
}

# Fixed model settings (no hyper-parameter search by design)
LR_C = 1.0
RF_N_ESTIMATORS = 300
RF_MIN_SAMPLES_LEAF = 5
