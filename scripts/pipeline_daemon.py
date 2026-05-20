#!/usr/bin/env python3
"""
REFOLD Proteome Pipeline Daemon
================================
Continuously processes pathogenic missense variants from ClinVar,
running the full REFOLD Stage 1→2→3 pipeline and injecting results
into PCD_global_atlas.json.

Auto-restarts on crash via macOS LaunchAgent.
Run manually:  python3 scripts/pipeline_daemon.py
Check status:  cat data/pipeline_state.json | python3 -m json.tool
"""

import json
import os
import sys
import time
import logging
import traceback
import shutil
import subprocess
from pathlib import Path
from datetime import datetime, timezone

# ── Paths ──────────────────────────────────────────────────────────────────
ROOT        = Path(__file__).parent.parent.resolve()
STATE_FILE  = ROOT / "data" / "pipeline_state.json"
ATLAS_FILE  = ROOT / "data" / "results" / "PCD_global_atlas.json"
PUBLIC_ATLAS = ROOT / "pcd-website" / "public" / "PCD_global_atlas.json"
PUBLIC_PDB   = ROOT / "pcd-website" / "public" / "structures"
CLINVAR_QUEUE = ROOT / "data" / "clinvar_pathogenic_missense_queue.json"
GITHUB_DATA_REPO = ROOT / "pcd-atlas-data"   # local clone of aesenthilvanan-coder/pcd-atlas-data
LOG_DIR     = ROOT / "logs"
LOG_FILE    = LOG_DIR / "pipeline_daemon.log"

LOG_DIR.mkdir(parents=True, exist_ok=True)
PUBLIC_PDB.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger("pcd_daemon")

# ── ClinVar variant queue ──────────────────────────────────────────────────
# Curated list of pathogenic missense variants across major misfolding diseases.
# Format: (gene, uniprot, mutation_mature, mutation_precursor, clinvar_id, disease, mechanism)
INITIAL_QUEUE = [
    # ── Lysosomal Storage ──────────────────────────────────────────────────
    ("GBA1",    "P04062", "N370S",  "N409S",  "4288",   "Gaucher Disease Type 1",           "ER retention / ERAD"),
    ("GBA1",    "P04062", "L444P",  "L483P",  "4288",   "Gaucher Disease Type 1",           "ER retention / ERAD"),
    ("HEXA",    "P06865", "G269S",  "G269S",  "3827",   "Tay-Sachs Disease",                "ER retention of β-hexosaminidase"),
    ("HEXB",    "P06865", "R211H",  "R211H",  "3828",   "Sandhoff Disease",                 "ER retention of β-hexosaminidase"),
    ("ARSA",    "P15289", "P426L",  "P426L",  "2740",   "Metachromatic Leukodystrophy",     "Misfolding of arylsulfatase A"),
    ("GLA",     "P06280", "A143T",  "A143T",  "9968",   "Fabry Disease",                    "ER retention of α-galactosidase"),
    ("GLA",     "P06280", "R301Q",  "R301Q",  "9969",   "Fabry Disease",                    "ER retention of α-galactosidase"),
    ("NPC1",    "O15118", "P1007A", "P1007A", "17828",  "Niemann-Pick Type C1",             "Cholesterol trafficking defect"),
    ("SMPD1",   "P17405", "L302P",  "L302P",  "9610",   "Niemann-Pick Type A",              "ER retention of sphingomyelinase"),
    ("GALC",    "P54803", "I562T",  "I562T",  "3793",   "Krabbe Disease",                   "Misfolding of galactocerebrosidase"),
    ("NAGLU",   "P20933", "R643C",  "R643C",  "3779",   "Sanfilippo Syndrome Type B",       "ER retention of α-N-acetylglucosaminidase"),
    # ── Amino Acid Metabolism ──────────────────────────────────────────────
    ("PAH",     "P00439", "R408W",  "R408W",  "699",    "Phenylketonuria",                  "Misfolding of phenylalanine hydroxylase"),
    ("PAH",     "P00439", "R261Q",  "R261Q",  "700",    "Phenylketonuria",                  "Misfolding of phenylalanine hydroxylase"),
    ("FAH",     "P16930", "W262R",  "W262R",  "2984",   "Tyrosinemia Type 1",               "ER retention of fumarylacetoacetase"),
    ("GALT",    "P07902", "Q188R",  "Q188R",  "3118",   "Galactosemia",                     "Misfolding of galactose-1-P uridyltransferase"),
    ("BCKDHA",  "P12694", "Y393N",  "Y393N",  "4981",   "Maple Syrup Urine Disease",        "Misfolding of BCKDH E1α"),
    # ── Respiratory ────────────────────────────────────────────────────────
    ("CFTR",    "P13569", "F508del","F508del", "7107",  "Cystic Fibrosis",                  "NBD1 misfolding / ΔF508"),
    ("CFTR",    "P13569", "G85E",   "G85E",    "7107",  "Cystic Fibrosis",                  "TM1 misfolding"),
    ("CFTR",    "P13569", "R117H",  "R117H",   "7107",  "Cystic Fibrosis",                  "Gating defect with misfolding component"),
    ("SERPINA1","P01009", "E342K",  "E366K",   "3764",  "Alpha-1 Antitrypsin Deficiency",   "Loop insertion polymerization"),
    # ── Neurological ───────────────────────────────────────────────────────
    ("TTR",     "P02766", "V30M",   "V30M",    "2101",  "Transthyretin Amyloidosis",        "Tetramer dissociation / amyloid"),
    ("TTR",     "P02766", "V122I",  "V122I",   "2102",  "Transthyretin Cardiomyopathy",     "Tetramer dissociation / amyloid"),
    ("SOD1",    "P00441", "A4V",    "A4V",     "4832",  "ALS (SOD1-related)",               "Misfolding of Cu/Zn superoxide dismutase"),
    ("SOD1",    "P00441", "G93A",   "G93A",    "4833",  "ALS (SOD1-related)",               "Misfolding of Cu/Zn superoxide dismutase"),
    ("LRRK2",  "Q5S007", "G2019S", "G2019S",  "9072",  "Parkinson Disease (LRRK2)",        "Kinase domain misfolding"),
    ("PINK1",  "Q9BXM7", "G309D",  "G309D",   "4327",  "Parkinson Disease (PINK1)",        "Mitochondrial kinase misfolding"),
    ("APP",    "P05067", "V717I",  "V717I",   "12345", "Alzheimer Disease (familial)",     "Altered APP cleavage / Aβ42 aggregation"),
    # ── Cardiac ────────────────────────────────────────────────────────────
    ("MYH7",   "P12883", "R403Q",  "R403Q",   "177",   "Hypertrophic Cardiomyopathy",      "Myosin motor domain misfolding"),
    ("KCNQ1",  "P51787", "R243H",  "R243H",   "3820",  "Long QT Syndrome Type 1",          "ER retention of Kv7.1 channel"),
    ("SCN5A",  "Q14524", "R1432G", "R1432G",  "3821",  "Brugada Syndrome",                 "ER retention of Nav1.5"),
    ("LDLR",   "P01130", "W23G",   "W23G",    "3766",  "Familial Hypercholesterolemia",    "ER retention of LDL receptor"),
    # ── Hematological ──────────────────────────────────────────────────────
    ("HBB",    "P68871", "E6V",    "E7V",     "15",    "Sickle Cell Disease",              "Deoxygenated Hb polymerization"),
    ("G6PD",   "P11413", "G202A",  "G202A",   "6649",  "G6PD Deficiency",                  "Misfolding reduces enzyme stability"),
    ("F8",     "P00451", "R2150H", "R2150H",  "4380",  "Hemophilia A",                     "ER retention of Factor VIII"),
    # ── Endocrine ──────────────────────────────────────────────────────────
    ("INS",    "P01308", "R46Q",   "R46Q",    "4082",  "Neonatal Diabetes Mellitus",       "Proinsulin misfolding / ER stress"),
    ("GCK",    "P35557", "G261R",  "G261R",   "3959",  "MODY Type 2",                      "Glucokinase misfolding"),
    ("PROP1",  "P51840", "R120C",  "R120C",   "9012",  "Combined Pituitary Hormone Def.",  "Homeodomain misfolding"),
    # ── Connective Tissue ──────────────────────────────────────────────────
    ("COL1A1", "P02452", "G277D",  "G277D",   "3872",  "Osteogenesis Imperfecta Type II",  "Triple helix destabilization"),
    ("FBN1",   "P35555", "C1117Y", "C1117Y",  "4018",  "Marfan Syndrome",                  "Fibrillin-1 domain misfolding"),
    # ── Immune / Other ─────────────────────────────────────────────────────
    ("BRCA1",  "P38398", "C61G",   "C61G",    "55466", "Hereditary Breast Cancer",         "RING domain misfolding"),
    ("TP53",   "P04637", "R175H",  "R175H",   "12375", "Li-Fraumeni Syndrome",             "DNA-binding domain misfolding"),
    ("TP53",   "P04637", "R248W",  "R248W",   "12376", "Li-Fraumeni Syndrome",             "DNA-binding domain misfolding"),
    ("ATP7B",  "P35670", "H1069Q", "H1069Q",  "3759",  "Wilson Disease",                   "ER retention of copper transporter"),
    ("PTEN",   "P60484", "H93R",   "H93R",    "9740",  "Cowden Syndrome",                  "PTEN phosphatase misfolding"),
    ("VHL",    "P40337", "R167Q",  "R167Q",   "8951",  "von Hippel-Lindau Syndrome",       "pVHL elongin complex misfolding"),
]


def load_state() -> dict:
    if STATE_FILE.exists():
        with open(STATE_FILE) as f:
            return json.load(f)

    # Build queue from real ClinVar data (178,597 pathogenic missense variants)
    queue = []
    if CLINVAR_QUEUE.exists():
        log.info(f"Loading ClinVar queue from {CLINVAR_QUEUE}...")
        with open(CLINVAR_QUEUE) as f:
            clinvar = json.load(f)
        for v in clinvar["queue"]:
            queue.append({
                "gene": v["gene"],
                "uniprot": "",          # fetched dynamically via UniProt API per gene
                "mutation_mature": v["mutation_mature"],
                "mutation_precursor": v["mutation_mature"],  # updated if signal peptide info available
                "clinvar_id": v.get("variation_id", v.get("clinvar_allele_id", "")),
                "disease": v.get("phenotype", "Unknown disease"),
                "mechanism": "Pathogenic missense — mechanism TBD",
                "wt_aa": v.get("wt_aa",""),
                "pos": v.get("pos", 0),
                "mut_aa": v.get("mut_aa",""),
            })
        log.info(f"Loaded {len(queue)} variants from ClinVar")
    else:
        # Fall back to curated list if ClinVar queue not yet generated
        log.warning("ClinVar queue file not found — using curated list")
        queue = [
            {"gene": g, "uniprot": u, "mutation_mature": mm, "mutation_precursor": mp,
             "clinvar_id": ci, "disease": dis, "mechanism": mech, "wt_aa":"","pos":0,"mut_aa":""}
            for g, u, mm, mp, ci, dis, mech in INITIAL_QUEUE
        ]

    return {
        "version": "2.0",
        "daemon_started": datetime.now(timezone.utc).isoformat(),
        "total_queued": len(queue),
        "total_processed": 0,
        "total_complete": 0,
        "total_failed": 0,
        "total_skipped": 0,
        "queue": queue,
        "completed": [],
        "failed": [],
        "skipped": [],
        "last_updated": datetime.now(timezone.utc).isoformat(),
    }


def save_state(state: dict):
    STATE_FILE.parent.mkdir(parents=True, exist_ok=True)
    state["last_updated"] = datetime.now(timezone.utc).isoformat()
    tmp = STATE_FILE.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(state, f, indent=2)
    tmp.replace(STATE_FILE)


def load_atlas() -> dict:
    if ATLAS_FILE.exists():
        with open(ATLAS_FILE) as f:
            return json.load(f)
    # Get actual total from ClinVar queue file
    total_target = 178597  # real ClinVar pathogenic missense count
    if CLINVAR_QUEUE.exists():
        try:
            with open(CLINVAR_QUEUE) as f:
                cq = json.load(f)
            total_target = cq.get("total", len(cq.get("queue", [])))
        except Exception:
            pass

    return {
        "database": "Pharmacological Chaperone Database (PCD)",
        "version": "1.0.0",
        "build_date": datetime.now(timezone.utc).date().isoformat(),
        "investigator": "Aaryan Senthilvanan",
        "institution": "S.Y.A.L.I.S Labs",
        "powered_by": "REFOLD",
        "total_entries": 0,
        "proteome_targets": {
            "total_clinvar_pathogenic_missense": total_target,
            "total_genes_affected": 5992,
            "total_processed": 0,
            "queue_remaining": total_target,
        },
        "entries": [],
    }


def save_atlas(atlas: dict):
    ATLAS_FILE.parent.mkdir(parents=True, exist_ok=True)
    tmp = ATLAS_FILE.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(atlas, f, indent=2)
    tmp.replace(ATLAS_FILE)
    # Sync to website public dir
    shutil.copy(ATLAS_FILE, PUBLIC_ATLAS)
    log.info(f"Atlas synced locally ({atlas['total_entries']} entries)")


ATLAS_RAW_BASE = "https://raw.githubusercontent.com/aesenthilvanan-coder/pcd-atlas-data/main"


def push_to_github(atlas: dict, pdb_src: Path | None = None, entry_id: str | None = None):
    """Commit and push atlas JSON (+ optional PDB) to pcd-atlas-data repo."""
    if not GITHUB_DATA_REPO.exists():
        log.warning("GitHub data repo not found — skipping push")
        return
    try:
        shutil.copy(ATLAS_FILE, GITHUB_DATA_REPO / "PCD_global_atlas.json")
        files_to_add = ["PCD_global_atlas.json"]

        if pdb_src and pdb_src.exists() and entry_id:
            struct_dir = GITHUB_DATA_REPO / "structures"
            struct_dir.mkdir(exist_ok=True)
            dest = struct_dir / f"{entry_id}.pdb"
            shutil.copy(pdb_src, dest)
            files_to_add.append(f"structures/{entry_id}.pdb")

        n = atlas["total_entries"]
        eid = entry_id or (atlas["entries"][-1]["entry_id"] if atlas["entries"] else "?")
        msg = f"Add {eid} — {n} total entries"
        subprocess.run(["git", "-C", str(GITHUB_DATA_REPO), "add"] + files_to_add,
                       check=True, capture_output=True)
        subprocess.run(["git", "-C", str(GITHUB_DATA_REPO), "commit", "-m", msg],
                       check=True, capture_output=True)
        subprocess.run(["git", "-C", str(GITHUB_DATA_REPO), "push"],
                       check=True, capture_output=True, timeout=30)
        log.info(f"✓ Pushed to GitHub: {msg}")
    except Exception as e:
        log.error(f"GitHub push failed: {e}")


def _update_proteome_count(atlas: dict, state: dict):
    """Update proteome_targets counters. total_processed = all variants attempted."""
    processed = state["total_processed"] + 1   # +1 because we haven't incremented yet
    total = atlas["proteome_targets"].get("total_clinvar_pathogenic_missense", 178597)
    atlas["proteome_targets"]["total_processed"] = processed
    atlas["proteome_targets"]["queue_remaining"] = max(0, total - processed)


def entry_exists(atlas: dict, entry_id: str) -> bool:
    return any(e["entry_id"] == entry_id for e in atlas["entries"])


def run_stage2(variant: dict, work_dir: Path) -> dict | None:
    """Run ANM conformational sampling + fpocket. Returns pocket summary or None."""
    gene = variant["gene"]
    mut  = variant["mutation_mature"]
    log.info(f"Stage 2: {gene} {mut} — ANM + fpocket")

    script = ROOT / "scripts" / f"analyze_{gene.lower()}.py"
    if not script.exists():
        # Use generic stage2 script
        script = ROOT / "scripts" / "analyze_generic.py"
    if not script.exists():
        log.warning(f"No Stage 2 script for {gene}. Skipping.")
        return None

    result = subprocess.run(
        [sys.executable, str(script),
         "--gene", gene,
         "--uniprot", variant["uniprot"],
         "--mutation", mut,
         "--outdir", str(work_dir)],
        capture_output=True, text=True, timeout=600
    )
    if result.returncode != 0:
        log.error(f"Stage 2 failed:\n{result.stderr[-1000:]}")
        return None

    summary_file = work_dir / "pocket_summary.json"
    if not summary_file.exists():
        log.warning("No pocket_summary.json produced.")
        return None

    with open(summary_file) as f:
        return json.load(f)


def run_stage3(variant: dict, pocket_summary: dict, work_dir: Path) -> dict | None:
    """Run chaperone generation. Returns chaperone dict or None."""
    gene = variant["gene"]
    mut  = variant["mutation_mature"]
    log.info(f"Stage 3: {gene} {mut} — chaperone generation")

    script = ROOT / "scripts" / f"generate_chaperone_{gene.lower()}.py"
    if not script.exists():
        script = ROOT / "scripts" / "generate_chaperone_generic.py"
    if not script.exists():
        log.warning(f"No Stage 3 script for {gene}. Skipping.")
        return None

    result = subprocess.run(
        [sys.executable, str(script),
         "--gene", gene,
         "--mutation", mut,
         "--pocket-summary", str(work_dir / "pocket_summary.json"),
         "--outdir", str(work_dir)],
        capture_output=True, text=True, timeout=1200
    )
    if result.returncode != 0:
        log.error(f"Stage 3 failed:\n{result.stderr[-1000:]}")
        return None

    cand_file = work_dir / "stage3_chaperone_candidates.json"
    if not cand_file.exists():
        return None

    with open(cand_file) as f:
        candidates = json.load(f)

    # Return the top-ranked candidate
    if not candidates:
        return None
    return candidates[0] if isinstance(candidates, list) else candidates


def build_atlas_entry(variant: dict, pocket_summary: dict, chaperone: dict,
                       entry_id: str, work_dir: Path) -> dict:
    """Assemble a complete PCDEntry dict matching the GBA1/CFTR schema."""
    now = datetime.now(timezone.utc).isoformat()
    pdb_src = next(work_dir.glob("conformations/mt_*.pdb"), None)
    pdb_public = f"{ATLAS_RAW_BASE}/structures/{entry_id}.pdb"
    if pdb_src:
        shutil.copy(pdb_src, PUBLIC_PDB / f"{entry_id}.pdb")

    eij_csv = work_dir / "E_ij_matrix.csv"
    eij_values = []
    if eij_csv.exists():
        import csv
        with open(eij_csv) as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                eij_values.append([round(float(x), 3) for x in row])

    # Pocket — merge raw pocket dict with computed fields
    raw_pocket = pocket_summary.get("pocket", {})

    # Sequence — use enriched sequence section from analyze_generic.py
    raw_seq = pocket_summary.get("sequence", {})
    # pocket_lining_residues at top level is a list of label strings e.g. ["A123","R124"]
    lining_labels = pocket_summary.get("pocket_lining_residues", [])
    if isinstance(lining_labels, list) and lining_labels and isinstance(lining_labels[0], dict):
        # Old format: list of dicts — convert
        lining_labels = [f"{r['aa']}{r['pos']}" for r in lining_labels]

    lining_str = raw_seq.get("pocket_lining_residues") or "-".join(lining_labels)
    lining_positions = raw_seq.get("pocket_lining_positions_precursor") or [
        int("".join(c for c in s if c.isdigit())) for s in lining_labels if any(c.isdigit() for c in s)
    ]
    fasta_seq = raw_seq.get("fasta_sequence") or raw_seq.get("full_sequence", "")
    full_seq   = raw_seq.get("full_sequence", fasta_seq)
    mut_pos    = raw_seq.get("mutation_pos_1indexed", 0)
    slc_start  = max(0, mut_pos - 11)
    slc_end    = min(len(full_seq), mut_pos + 9)
    seq_slice  = raw_seq.get("sequence_slice_around_pocket") or full_seq[slc_start:slc_end]

    pocket_section = {
        "target_conformation":               raw_pocket.get("target_conformation",
                                              f"ANM mutant conformation (of {raw_pocket.get('n_conformations_sampled', 20)})"),
        "n_conformations_sampled":           raw_pocket.get("n_conformations_sampled", 20),
        "fpocket_druggability":              raw_pocket.get("fpocket_druggability", 0.0),
        "wt_baseline_druggability":          raw_pocket.get("wt_baseline_druggability", 0.0),
        "volume_angstrom3":                  raw_pocket.get("volume_angstrom3", 0.0),
        "alpha_sphere_count":                raw_pocket.get("alpha_sphere_count", 0),
        "center_angstrom":                   raw_pocket.get("center_angstrom", [0.0, 0.0, 0.0]),
        "dist_mutation_to_pocket_angstrom":  raw_pocket.get("dist_mutation_to_pocket_angstrom", 0.0),
        "exceeds_threshold":                 raw_pocket.get("exceeds_threshold", True),
        "pocket_type":                       raw_pocket.get("pocket_type",
                                              "Cryptic" if raw_pocket.get("cryptic", True) else "Allosteric"),
        "detection_frequency":               raw_pocket.get("detection_frequency", 1.0),
        "cryptic":                           raw_pocket.get("cryptic", True),
        "pocket_residues":                   raw_pocket.get("pocket_residues", {}),
        "n_lining_residues":                 raw_pocket.get("n_lining_residues", len(lining_labels)),
        "shell_radius_angstrom":             raw_pocket.get("shell_radius_angstrom", 8.0),
    }

    sequence_section = {
        "fasta_header":                   raw_seq.get("fasta_header",
                                           f">sp|{variant['uniprot']}|{variant['gene']} variant {variant['mutation_mature']}"),
        "fasta_sequence":                 fasta_seq,
        "full_sequence":                  full_seq,
        "n_residues":                     raw_seq.get("n_residues", len(fasta_seq)),
        "mutation_pos_1indexed":          mut_pos,
        "wt_aa":                          raw_seq.get("wt_aa", ""),
        "mut_aa":                         raw_seq.get("mut_aa", ""),
        "pocket_lining_residues":         lining_str,
        "pocket_lining_positions_precursor": lining_positions,
        "pocket_lining_1letter":          raw_seq.get("pocket_lining_1letter", "".join(s[0] for s in lining_labels if s)),
        "sequence_slice_around_pocket":   seq_slice,
    }

    return {
        "entry_id":      entry_id,
        "pdb_structure": pdb_public,
        "investigator":  "Aaryan Senthilvanan",
        "institution":   "S.Y.A.L.I.S Labs",
        "metadata": {
            "gene":               variant["gene"],
            "uniprot":            variant["uniprot"],
            "clinvar_id":         variant["clinvar_id"],
            "disease":            variant["disease"],
            "mutation_precursor": variant["mutation_precursor"],
            "mutation_mature":    variant["mutation_mature"],
            "variant_class":      "pathogenic",
            "mechanism":          variant["mechanism"],
        },
        "pocket":   pocket_section,
        "sequence": sequence_section,
        "eij_matrix": {
            "file":           str(work_dir / "E_ij_matrix.npy"),
            "shape":          pocket_summary.get("eij_shape", [0, 0]),
            "residue_labels": pocket_summary.get("residue_labels", lining_labels),
            "values":         eij_values,
        },
        "chaperone": chaperone,
        "assets": {
            "eij_matrix_npy":           str(work_dir / "E_ij_matrix.npy"),
            "eij_matrix_csv":           str(eij_csv) if eij_csv.exists() else None,
            "transient_conformation_pdb": str(pdb_src) if pdb_src else "",
            "pocket_summary_json":      str(work_dir / "pocket_summary.json"),
            "stage3_candidates_json":   str(work_dir / "stage3_chaperone_candidates.json"),
            "pdb_structure":            pdb_public,
        },
        "status":    "COMPLETE",
        "created_at": now,
    }


def process_variant(variant: dict, state: dict, atlas: dict) -> str:
    """
    Returns 'complete', 'failed', or 'skipped'.
    """
    gene = variant["gene"]
    mut  = variant["mutation_mature"]
    entry_id = f"PCD-{gene}-{mut}"

    if entry_exists(atlas, entry_id):
        log.info(f"Already in atlas: {entry_id}. Skipping.")
        _update_proteome_count(atlas, state)
        # Push to GitHub every 50 skips so the progress bar still advances
        if state["total_skipped"] % 50 == 0:
            save_atlas(atlas)
            push_to_github(atlas)
        return "skipped"

    work_dir = ROOT / "data" / "results" / f"{gene}_{mut}"
    work_dir.mkdir(parents=True, exist_ok=True)

    try:
        pocket_summary = run_stage2(variant, work_dir)
        if pocket_summary is None:
            log.warning(f"{entry_id}: Stage 2 produced no valid pocket.")
            _update_proteome_count(atlas, state)
            # Push every 20 failures so progress bar advances
            if state["total_failed"] % 20 == 0:
                save_atlas(atlas)
                push_to_github(atlas)
            return "failed"

        drug_score = pocket_summary.get("pocket", {}).get("fpocket_druggability", 0)
        if drug_score < 0.70:
            log.info(f"{entry_id}: druggability={drug_score:.3f} < 0.70 threshold. Skipping.")
            return "skipped"

        chaperone = run_stage3(variant, pocket_summary, work_dir)
        if chaperone is None:
            log.warning(f"{entry_id}: Stage 3 failed.")
            return "failed"

        entry = build_atlas_entry(variant, pocket_summary, chaperone, entry_id, work_dir)
        atlas["entries"].append(entry)
        atlas["total_entries"] = len(atlas["entries"])
        _update_proteome_count(atlas, state)
        save_atlas(atlas)
        pdb_src = next(work_dir.glob("conformations/mt_*.pdb"), None)
        push_to_github(atlas, pdb_src=pdb_src, entry_id=entry_id)

        log.info(f"✓ {entry_id} — drug={drug_score:.3f} — composite={chaperone.get('composite_score',0):.3f}")
        return "complete"

    except Exception as e:
        log.error(f"{entry_id} error:\n{traceback.format_exc()}")
        return "failed"


def main():
    log.info("=" * 60)
    log.info("REFOLD Proteome Pipeline Daemon — starting")
    log.info(f"Root: {ROOT}")
    log.info("=" * 60)

    while True:
        try:
            state = load_state()
            atlas = load_atlas()

            if not state["queue"]:
                log.info("Queue empty. Daemon complete. Sleeping 1 hour.")
                time.sleep(3600)
                continue

            variant = state["queue"].pop(0)
            gene = variant["gene"]
            mut  = variant["mutation_mature"]
            log.info(f"\n── Processing {gene} {mut} ({variant['disease']}) ──")
            log.info(f"   Queue remaining: {len(state['queue'])}")

            result = process_variant(variant, state, atlas)

            state["total_processed"] += 1
            if result == "complete":
                state["total_complete"] += 1
                state["completed"].append({"gene": gene, "mutation": mut, "ts": datetime.now(timezone.utc).isoformat()})
            elif result == "failed":
                state["total_failed"] += 1
                state["failed"].append({"gene": gene, "mutation": mut, "ts": datetime.now(timezone.utc).isoformat()})
            elif result == "skipped":
                state["total_skipped"] += 1
                state["skipped"].append({"gene": gene, "mutation": mut})

            save_state(state)

            log.info(f"Progress: {state['total_complete']} complete, "
                     f"{state['total_skipped']} skipped, "
                     f"{state['total_failed']} failed, "
                     f"{len(state['queue'])} remaining")

            # Pause between variants (be kind to filesystem + external APIs)
            time.sleep(10)

        except KeyboardInterrupt:
            log.info("Daemon interrupted by user.")
            break
        except Exception as e:
            log.error(f"Daemon loop error:\n{traceback.format_exc()}")
            time.sleep(30)


if __name__ == "__main__":
    main()
