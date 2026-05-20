"""REFOLD global constants. All magic numbers live here. Never inline them elsewhere."""

from pathlib import Path
from typing import Final, FrozenSet

# ─────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────

ROOT_DIR: Final[Path] = Path(__file__).parent.parent.parent
DATA_DIR: Final[Path] = ROOT_DIR / "data"
RAW_DIR: Final[Path] = DATA_DIR / "raw"
PROCESSED_DIR: Final[Path] = DATA_DIR / "processed"
RESULTS_DIR: Final[Path] = DATA_DIR / "results"
CHECKPOINT_DIR: Final[Path] = ROOT_DIR / "checkpoints"
LOG_DIR: Final[Path] = ROOT_DIR / "logs"

# ─────────────────────────────────────────────────────────────
# AMINO ACIDS
# ─────────────────────────────────────────────────────────────

STANDARD_AAS: Final[str] = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX: Final[dict[str, int]] = {aa: i for i, aa in enumerate(STANDARD_AAS)}
IDX_TO_AA: Final[dict[int, str]] = {i: aa for aa, i in AA_TO_IDX.items()}
NUM_AAS: Final[int] = 20

AA_ONE_TO_THREE: Final[dict[str, str]] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
AA_THREE_TO_ONE: Final[dict[str, str]] = {v: k for k, v in AA_ONE_TO_THREE.items()}

AA_PROPERTIES: Final[dict] = {
    "A": {"hydrophobicity": 1.8,  "charge": 0,  "mw": 89.1,  "volume": 88.6},
    "C": {"hydrophobicity": 2.5,  "charge": 0,  "mw": 121.2, "volume": 108.5},
    "D": {"hydrophobicity": -3.5, "charge": -1, "mw": 133.1, "volume": 111.1},
    "E": {"hydrophobicity": -3.5, "charge": -1, "mw": 147.1, "volume": 138.4},
    "F": {"hydrophobicity": 2.8,  "charge": 0,  "mw": 165.2, "volume": 189.9},
    "G": {"hydrophobicity": -0.4, "charge": 0,  "mw": 75.1,  "volume": 60.1},
    "H": {"hydrophobicity": -3.2, "charge": 0,  "mw": 155.2, "volume": 153.2},
    "I": {"hydrophobicity": 4.5,  "charge": 0,  "mw": 131.2, "volume": 166.7},
    "K": {"hydrophobicity": -3.9, "charge": 1,  "mw": 146.2, "volume": 168.6},
    "L": {"hydrophobicity": 3.8,  "charge": 0,  "mw": 131.2, "volume": 166.7},
    "M": {"hydrophobicity": 1.9,  "charge": 0,  "mw": 149.2, "volume": 162.9},
    "N": {"hydrophobicity": -3.5, "charge": 0,  "mw": 132.1, "volume": 114.1},
    "P": {"hydrophobicity": -1.6, "charge": 0,  "mw": 115.1, "volume": 112.7},
    "Q": {"hydrophobicity": -3.5, "charge": 0,  "mw": 146.1, "volume": 143.8},
    "R": {"hydrophobicity": -4.5, "charge": 1,  "mw": 174.2, "volume": 173.4},
    "S": {"hydrophobicity": -0.8, "charge": 0,  "mw": 105.1, "volume": 89.0},
    "T": {"hydrophobicity": -0.7, "charge": 0,  "mw": 119.1, "volume": 116.1},
    "V": {"hydrophobicity": 4.2,  "charge": 0,  "mw": 117.1, "volume": 140.0},
    "W": {"hydrophobicity": -0.9, "charge": 0,  "mw": 204.2, "volume": 227.8},
    "Y": {"hydrophobicity": -1.3, "charge": 0,  "mw": 181.2, "volume": 193.6},
}

# ─────────────────────────────────────────────────────────────
# STRUCTURE
# ─────────────────────────────────────────────────────────────

ATOM_NAMES: Final[list[str]] = [
    "N", "CA", "C", "CB", "O", "CG", "CG1", "CG2", "OG", "OG1",
    "SG", "CD", "CD1", "CD2", "ND1", "ND2", "OD1", "OD2", "SD",
    "CE", "CE1", "CE2", "CE3", "NE", "NE1", "NE2", "OE1", "OE2",
    "CH2", "NH1", "NH2", "OH", "CZ", "CZ2", "CZ3", "NZ", "OXT",
]
ATOM_TO_IDX: Final[dict[str, int]] = {a: i for i, a in enumerate(ATOM_NAMES)}
CA_IDX: Final[int] = 1       # Cα is ALWAYS index 1
N_IDX: Final[int] = 0        # backbone N
C_IDX: Final[int] = 2        # backbone C
O_IDX: Final[int] = 4        # backbone O
CB_IDX: Final[int] = 3       # Cβ
N_ATOM_TYPES: Final[int] = 37

INTERFACE_HEAVY_ATOM_CUTOFF: Final[float] = 5.0
CONTACT_CA_CUTOFF: Final[float] = 8.0
ENM_CUTOFF: Final[float] = 12.0
SURFACE_PROBE_RADIUS: Final[float] = 1.4
MIN_POCKET_VOLUME: Final[float] = 300.0
MAX_POCKET_VOLUME: Final[float] = 2000.0

DSSP_SS_MAP: Final[dict[str, int]] = {
    "H": 0,   # α-helix
    "B": 1,   # β-bridge
    "E": 2,   # β-strand
    "G": 3,   # 3₁₀-helix
    "I": 4,   # π-helix
    "T": 5,   # turn
    "S": 6,   # bend
    " ": 7,   # coil
}

MAX_ASA: Final[dict[str, float]] = {
    "A": 121.0, "R": 265.0, "N": 187.0, "D": 187.0, "C": 148.0,
    "Q": 214.0, "E": 214.0, "G": 97.0,  "H": 216.0, "I": 195.0,
    "L": 191.0, "K": 230.0, "M": 203.0, "F": 228.0, "P": 154.0,
    "S": 143.0, "T": 163.0, "W": 264.0, "Y": 255.0, "V": 165.0,
}

# ─────────────────────────────────────────────────────────────
# RESCUE CLASSIFIER THRESHOLDS
# ─────────────────────────────────────────────────────────────

DDG_DESTABILIZING_THRESHOLD: Final[float] = 0.5
DDG_SEVERELY_UNSTABLE: Final[float] = 5.0
DDG_ACTIVE_SITE_CUTOFF: Final[float] = 5.0

RESCUE_POSITIVE_THRESHOLD: Final[float] = 0.5
RESCUE_PROBABILITY_THRESHOLD: Final[float] = 0.5
HIGH_CONFIDENCE_THRESHOLD: Final[float] = 0.8

# ─────────────────────────────────────────────────────────────
# ENM / CONFORMATIONAL SAMPLING
# ─────────────────────────────────────────────────────────────

N_ENM_MODES: Final[int] = 20
N_CONFORMATIONS: Final[int] = 50
ENM_AMPLITUDE_SCALE: Final[float] = 1.0
ENM_N_SKIP_TRIVIAL: Final[int] = 6

# ─────────────────────────────────────────────────────────────
# POCKET DETECTION
# ─────────────────────────────────────────────────────────────

FPOCKET_MIN_ALPHA_SPHERE_RADIUS: Final[float] = 3.0
FPOCKET_MAX_ALPHA_SPHERE_RADIUS: Final[float] = 6.0
FPOCKET_MIN_DRUGGABILITY_SCORE: Final[float] = 0.5
POCKET_OVERLAP_THRESHOLD: Final[float] = 0.3
TRANSIENT_POCKET_FREQ_THRESHOLD: Final[float] = 0.20
DRUGGABILITY_THRESHOLD: Final[float] = 0.50

# ─────────────────────────────────────────────────────────────
# MOLECULE GENERATION
# ─────────────────────────────────────────────────────────────

MAX_MW: Final[float] = 500.0
MAX_LOGP: Final[float] = 5.0
MAX_HBD: Final[int] = 5
MAX_HBA: Final[int] = 10
MAX_ROTATABLE_BONDS: Final[int] = 10
MAX_TPSA: Final[float] = 140.0
MAX_SA_SCORE: Final[float] = 4.0
MIN_QED: Final[float] = 0.30

N_MOLECULES_PER_POCKET: Final[int] = 1000
N_DIFFUSION_STEPS_DEFAULT: Final[int] = 1000
N_DIFFUSION_STEPS_M1: Final[int] = 500
DIFFUSION_BETA_START: Final[float] = 1e-4
DIFFUSION_BETA_END: Final[float] = 0.02
N_ATOMS_MIN: Final[int] = 10
N_ATOMS_MAX: Final[int] = 35

MOL_ATOM_TYPES: Final[list[str]] = ["C", "N", "O", "S", "F", "P", "Cl", "Br", "I"]
N_MOL_ATOM_TYPES: Final[int] = len(MOL_ATOM_TYPES)

BOND_TYPES: Final[list[str]] = ["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC"]
N_BOND_TYPES: Final[int] = len(BOND_TYPES)

# ─────────────────────────────────────────────────────────────
# ADMET
# ─────────────────────────────────────────────────────────────

BBB_POSITIVE_THRESHOLD: Final[float] = 0.5
BIOAVAILABILITY_THRESHOLD: Final[float] = 0.5
HERG_INHIBITION_THRESHOLD: Final[float] = 0.5
AMES_POSITIVE_THRESHOLD: Final[float] = 0.5

# ─────────────────────────────────────────────────────────────
# DATA SOURCES — URLS
# ─────────────────────────────────────────────────────────────

CLINVAR_VCF_URL: Final[str] = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
    "clinvar.vcf.gz"
)
CLINVAR_XML_URL: Final[str] = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/"
    "ClinVarFullRelease_00-latest.xml.gz"
)
ALPHAFOLD_BASE_URL: Final[str] = "https://alphafold.ebi.ac.uk/files/"
UNIPROT_FASTA_URL: Final[str] = (
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
    "knowledgebase/reference_proteomes/Eukaryota/"
    "UP000005640/UP000005640_9606.fasta.gz"
)
THERMOMUTDB_URL: Final[str] = (
    "https://biosig.lab.uq.edu.au/thermomutdb/static/media/"
    "thermomutdb_json.zip"
)
PROTHERM_URL: Final[str] = (
    "https://web.iitm.ac.in/bioinfo2/prothermdb/"
    "rawdata.zip"
)
CHEMBL_URL: Final[str] = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    "chembl_34_sqlite.tar.gz"
)
BINDINGDB_URL: Final[str] = (
    "https://www.bindingdb.org/bind/downloads/"
    "BindingDB_All_2024m1.tsv.zip"
)

# ─────────────────────────────────────────────────────────────
# MODEL ARCHITECTURE DEFAULTS
# ─────────────────────────────────────────────────────────────

RESCUE_CLASSIFIER_HIDDEN_DIMS: Final[list[int]] = [512, 256, 128, 64]
RESCUE_CLASSIFIER_DROPOUT: Final[float] = 0.2

POCKET_GNN_HIDDEN_DIM: Final[int] = 256
POCKET_GNN_N_LAYERS: Final[int] = 6
POCKET_GNN_N_HEADS: Final[int] = 4

MOL_GEN_HIDDEN_DIM: Final[int] = 256
MOL_GEN_N_LAYERS: Final[int] = 6
MOL_GEN_N_HEADS: Final[int] = 8
MOL_GEN_POCKET_EMBED_DIM: Final[int] = 128

ADMET_HIDDEN_DIMS: Final[list[int]] = [512, 256, 128]
ADMET_N_TASKS: Final[int] = 8

# ─────────────────────────────────────────────────────────────
# FEATURE DIMENSIONS
# ─────────────────────────────────────────────────────────────

GNN_NODE_DIM: Final[int] = 30
GNN_EDGE_DIM: Final[int] = 20
ESM2_EMBED_DIM: Final[int] = 480
ESM2_LARGE_EMBED_DIM: Final[int] = 1280
ESM1V_MAX_SEQ_LEN: Final[int] = 1022
THERMO_FEAT_DIM: Final[int] = 16
EVO_FEAT_DIM: Final[int] = 32

POCKET_EMBED_DIM: Final[int] = 128
POCKET_NODE_FEAT_DIM: Final[int] = 25
POCKET_SCALAR_DIM: Final[int] = 8

# ─────────────────────────────────────────────────────────────
# TRAINING / OPTIMISATION
# ─────────────────────────────────────────────────────────────

EMA_DECAY: Final[float] = 0.9999
GRAD_CLIP_NORM: Final[float] = 1.0
FOCAL_LOSS_GAMMA: Final[float] = 2.0
FOCAL_LOSS_ALPHA: Final[float] = 0.25

ALPHAFOLD_RATE_LIMIT: Final[float] = 0.2

# ─────────────────────────────────────────────────────────────
# BENCHMARK MUTATIONS (known chaperone-responsive)
# ─────────────────────────────────────────────────────────────

BENCHMARK_MUTATIONS: FrozenSet[str] = frozenset({
    # Fabry disease (GLA, α-galactosidase A)
    "P06280_152_Y_C",
    "P06280_215_R_W",
    "P06280_231_R_H",
    # Pompe disease (GAA)
    "P10253_525_L_S",
    "P10253_600_P_L",
    # Gaucher disease (GBA)
    "P04062_444_N_S",
    "P04062_296_L_V",
    # CFTR (cystic fibrosis)
    "P13569_508_F_C",
    "P13569_551_R_H",
    # TTR (transthyretin amyloidosis)
    "P02766_30_V_M",
    "P02766_122_V_I",
    # p53 (Li-Fraumeni syndrome)
    "P04637_175_R_H",
    "P04637_248_R_W",
    "P04637_249_R_S",
})
