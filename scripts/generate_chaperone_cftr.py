"""
Stage 3: Pharmacological Chaperone Generation for CFTR G85E.

Target: ICL1-TM2 junction pocket (residues I70/L73/R74/F77/F78/F81/L195)
        fpocket druggability = 0.807
        Pocket center: (-5.46, 19.94, -22.07) Å
        VX-809 canonical site: drug=0.003 (destroyed)

Pharmacophore features:
  - Aromatic π-stacker for the F77-F78-F81 cage (three Phe in TM2)
  - H-bond acceptor for R74 Arg guanidinium
  - Hydrophobic contacts for I70, L73, L195
  - logP 2-4, MW 280-420 Da (TM-accessible corrector)
  - Basic N preferred (lysosomal/ER trapping aids CFTR ER-to-Golgi trafficking)

Seeds include known CFTR corrector scaffolds and TM-binding analogs.
"""

import sys, json, random, logging
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import (
        Descriptors, rdMolDescriptors, AllChem, RWMol, DataStructs
    )
    from rdkit.Chem.rdMolDescriptors import CalcTPSA
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
except ImportError:
    logger.error("RDKit not available"); sys.exit(1)

# ── Pocket pharmacophore parameters ──────────────────────────────────────────
POCKET = {
    "residues":      ["I70","L73","R74","F77","F78","F81","L195"],
    "center":        (-5.458, 19.942, -22.066),
    "druggability":  0.807,
    "volume":        484.0,
    "hydrophobicity": 2.29,        # highly hydrophobic (Kyte-Doolittle mean)
    "n_phe":         3,            # aromatic cage F77-F78-F81
    "has_arg":       True,         # R74 → wants H-bond acceptor
}

# ── Seed SMILES library ───────────────────────────────────────────────────────
# Curated for CFTR TMD1 correctors: aromatic cores, amide/sulfonamide HBA for R74,
# hydrophobic tails for Ile/Leu contacts, basic N for ER localisation
SEEDS = [
    # VX-809 / lumacaftor scaffold series
    "OC(=O)c1ccc(NC(=O)Nc2ccc(F)cc2)cc1",            # VX-809 core analog
    "OC(=O)c1cc(F)ccc1NC(=O)Nc1ccc(Cl)cc1",
    "Cc1ccc(NC(=O)c2ccncc2)cc1C(=O)N1CCCC1",         # pyridine amide + pip
    # Fluorobenzamide + piperidine (VX-661/tezacaftor class)
    "O=C(c1ccc(F)cc1)N1CCCCC1",
    "O=C(c1cc(F)ccc1F)N1CCN(Cc2ccccc2)CC1",
    "O=C(Nc1ccc(F)cc1)C1CCNCC1",
    "O=C(Nc1cc(F)cc(F)c1)C1CCN(C)CC1",
    # Indole / benzimidazole cores (planar → F77-F78-F81 stacking)
    "c1ccc2[nH]ccc2c1",                               # indole
    "O=C(Nc1nc2ccccc2[nH]1)C1CCNCC1",                # benzimidazole-pipC
    "Cc1nc2ccccc2n1CC(=O)N1CCCCC1",
    "O=C(c1ccc2[nH]ccc2c1)N1CCNCC1",                 # indole-3-carbonyl pip
    "c1cnc2ccccc2c1",                                  # quinoline
    "O=C(Nc1ccncc1)c1ccc2ccccc2c1",                   # naphthalene-pyridine amide
    # Phenethylamine / benzylamine (compact aromatic + basic N)
    "NCCc1ccccc1",
    "NCCc1ccc(F)cc1",
    "NC(Cc1ccccc1)C(=O)NC1CCCCC1",                   # Phe analog
    "NC(Cc1ccc(F)cc1)C(=O)N1CCCC1",
    # Sulfonamide series (strong H-bond acceptor for R74)
    "O=S(=O)(Nc1ccccc1)c1ccc(N)cc1",
    "O=S(=O)(N1CCCCC1)c1ccc(F)cc1",
    "CS(=O)(=O)Nc1ccc(C(=O)N2CCCCC2)cc1",
    # Diphenylmethane / diaryl (fill aromatic cage breadth)
    "O=C(NCc1ccccc1)c1ccccc1",
    "O=C(NCc1ccc(F)cc1)c1ccc(N)cc1",
    "Fc1ccc(CNCc2ccccn2)cc1",                         # fluorobenzyl-pyridyl
    # Heterocyclic correctors: pyrazole, triazole
    "Cc1cc(-c2ccc(F)cc2)n[nH]1",
    "Fc1ccc(-c2nnc(N3CCCC3)s2)cc1",                  # thiadiazole-F-phenyl
    "O=c1cc(-c2ccc(F)cc2)cc(=O)[nH]1",               # 3-Fphenyl-chromone
    # Amide + aromatic + basic N triad (CFTR corrector pharmacophore)
    "O=C(Nc1ccc(F)cc1)N1CCCCC1",                     # carbamate-F-phenyl
    "O=C(Nc1ccc(F)cc1)c1ccncc1",                     # nicotinamide-F-phenyl
    "CC(=O)Nc1ccc(Cc2ccc(N3CCCCC3)cc2)cc1",
    "O=C(c1ccc2c(c1)CCCC2)N1CCCCC1",                 # tetralin amide
]

# ── REFOLD scoring function ───────────────────────────────────────────────────
def lipinski_ok(mol) -> bool:
    return (Descriptors.MolWt(mol) <= 500
            and Descriptors.MolLogP(mol) <= 5
            and rdMolDescriptors.CalcNumHBD(mol) <= 5
            and rdMolDescriptors.CalcNumHBA(mol) <= 10)

def veber_ok(mol) -> bool:
    return (rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
            and CalcTPSA(mol) <= 140)

def sa_score(mol) -> float:
    try:
        from rdkit.Chem import RDConfig
        import os, importlib.util
        sas_path = os.path.join(RDConfig.RDContribDir, "SA_Score", "sascorer.py")
        spec = importlib.util.spec_from_file_location("sascorer", sas_path)
        sas = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sas)
        return sas.calculateScore(mol)
    except Exception:
        return 3.0   # default if not available

def pains_ok(mol) -> bool:
    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)
        return not catalog.HasMatch(mol)
    except Exception:
        return True

def pocket_affinity(mol) -> float:
    """
    Heuristic affinity score for CFTR ICL1-TM2 pocket pharmacophore:
      - Aromatic ring(s) for F77-F78-F81 cage (strong weight)
      - H-bond acceptor (amide/sulfonamide) for R74
      - logP 2–4 for TM accessibility
      - MW 280–420 (corrector sweet spot)
      - Basic N (pKa ~ 8–10, ER/lysosomal trapping)
    """
    score = 0.0

    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    naro = rdMolDescriptors.CalcNumAromaticRings(mol)
    tpsa = CalcTPSA(mol)

    # Aromatic rings (critical — F77/F78/F81 cage)
    if naro >= 1: score += 0.25
    if naro >= 2: score += 0.15   # two rings → sandwich F77-F81

    # logP 2–4 (TM corrector window)
    if 2.0 <= logp <= 4.0:
        score += 0.20
    elif 1.5 <= logp < 2.0 or 4.0 < logp <= 5.0:
        score += 0.08

    # MW 280–420
    if 280 <= mw <= 420:
        score += 0.15
    elif 240 <= mw < 280 or 420 < mw <= 450:
        score += 0.07

    # H-bond acceptors (for R74 guanidinium: amide, sulfonamide, pyridine N)
    if 2 <= hba <= 5: score += 0.12
    if hbd <= 2:      score += 0.08   # low HBD → better membrane permeability

    # Basic nitrogen (ER trapping aids CFTR ER-to-Golgi rescue)
    basic_n_patt = Chem.MolFromSmarts("[N;!$(N-C=O);!$(N=*);H0,H1;+0]")
    if basic_n_patt and mol.HasSubstructMatch(basic_n_patt):
        score += 0.05

    return min(score, 1.0)

def score_molecule(smi: str) -> tuple[float, dict] | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    smi_canon = Chem.MolToSmiles(mol)
    if len(smi_canon) < 5:
        return None
    if not lipinski_ok(mol):
        return None
    if not veber_ok(mol):
        return None
    sa = sa_score(mol)
    if sa > 5.0:
        return None
    qed = Descriptors.qed(mol)
    if qed < 0.2:
        return None
    if not pains_ok(mol):
        return None

    aff = pocket_affinity(mol)
    mw  = Descriptors.MolWt(mol)
    lp  = Descriptors.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)

    composite = 0.45 * aff + 0.30 * qed + 0.15 * max(0, (5 - sa) / 4) + 0.10 * min(lp / 4, 1)

    return composite, {
        "smiles": smi_canon, "score": composite, "mw": mw,
        "logp": lp, "qed": qed, "sa": sa, "hbd": hbd, "hba": hba,
        "affinity": aff,
    }

# ── Mutation operators ────────────────────────────────────────────────────────
_FRAGMENTS = [
    "c1ccccc1", "c1ccncc1", "c1ccc(F)cc1", "c1cc(F)cc(F)c1", "c1ccc(Cl)cc1",
    "C1CCNCC1", "C1CCCNC1", "C1CCN(C)CC1", "N1CCOCC1",
    "C(=O)N", "S(=O)(=O)N", "C(=O)O", "NC", "NC(=O)",
    "c1ccc2[nH]ccc2c1", "c1cnc2ccccc2c1",
]

def _add_substituent(smi: str, rng: random.Random) -> str | None:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    rwmol = RWMol(mol)
    atoms = [a.GetIdx() for a in rwmol.GetAtoms()
             if a.GetAtomicNum() in (6, 7, 8)
             and a.GetNumImplicitHs() > 0]
    if not atoms:
        return None
    sub = rng.choice([
        "[F]", "[Cl]", "[OH]", "C", "[NH2]", "OC", "NC",
        "C(=O)N", "S(=O)(=O)N", "CC",
    ])
    sub_mol = Chem.MolFromSmarts(sub)
    if sub_mol is None:
        return None
    idx = rng.choice(atoms)
    try:
        rwmol.GetAtomWithIdx(idx).SetNumExplicitHs(0)
        new_smi = smi + "." + sub  # simplified: try as fragment combine
        return None   # skip unsafe edit; use crossover instead
    except Exception:
        return None

def mutate(smi: str, rng: random.Random) -> str | None:
    """Apply a random mutation to a SMILES string."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    op = rng.randint(0, 5)
    try:
        if op == 0:
            # Add fluorine to aromatic ring
            patt = Chem.MolFromSmarts("c1ccccc1")
            if patt and mol.HasSubstructMatch(patt):
                new = smi.replace("c1ccccc1", "c1ccc(F)cc1", 1)
                if Chem.MolFromSmiles(new):
                    return new
        elif op == 1:
            # Swap benzene for pyridine
            if "c1ccccc1" in smi:
                new = smi.replace("c1ccccc1", "c1ccncc1", 1)
                if Chem.MolFromSmiles(new):
                    return new
        elif op == 2:
            # Replace OCC with NCC (add basic nitrogen)
            if "OCC" in smi:
                new = smi.replace("OCC", "NCC", 1)
                if Chem.MolFromSmiles(new):
                    return new
        elif op == 3:
            # Extend piperidine linker by one CH2
            if "C1CCNCC1" in smi:
                new = smi.replace("C1CCNCC1", "C1CCCNCC1", 1)
                if Chem.MolFromSmiles(new):
                    return new
        elif op == 4:
            # Add methyl to aromatic
            if "c1ccccc1" in smi or "c1ccncc1" in smi:
                new = smi.replace("cc1)", "c(C)c1)", 1)
                if Chem.MolFromSmiles(new):
                    return new
        elif op == 5:
            # Replace single-ring with naphthalene for deeper Phe cage penetration
            if "c1ccccc1" in smi:
                new = smi.replace("c1ccccc1", "c1ccc2ccccc2c1", 1)
                if Chem.MolFromSmiles(new):
                    return new
    except Exception:
        pass
    return None

def crossover(s1: str, s2: str, rng: random.Random) -> str | None:
    """Scaffold fragment swap: replace an aromatic ring from s2 into s1."""
    m1, m2 = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
    if m1 is None or m2 is None:
        return None
    # Simple approach: pick a random fragment from _FRAGMENTS and substitute
    # the first aromatic ring of s1 with it
    frag = rng.choice(_FRAGMENTS[:8])   # only validated ring fragments
    fmol = Chem.MolFromSmiles(frag)
    if fmol is None:
        return None
    # Replace first 'c1ccccc1' occurrence if present
    if "c1ccccc1" in s1:
        new = s1.replace("c1ccccc1", frag, 1)
        m = Chem.MolFromSmiles(new)
        if m is not None:
            return Chem.MolToSmiles(m)
    return None

# ── Evolutionary loop ─────────────────────────────────────────────────────────
def run_evolution(
    seeds: list[str],
    n_generations: int = 50,
    pop_size: int = 40,
    n_elite: int = 10,
    seed: int = 42,
) -> list[dict]:
    rng = random.Random(seed)
    np_rng = np.random.default_rng(seed)

    logger.info("Phase 1: Scoring seed SMILES library...")
    scored = []
    for smi in seeds:
        r = score_molecule(smi)
        if r:
            scored.append(r[1])
    logger.info(f"  {len(scored)} valid seeds")

    population = sorted(scored, key=lambda x: x["score"], reverse=True)[:pop_size]
    if not population:
        logger.error("No valid seed molecules")
        return []

    logger.info(f"Phase 2: Evolutionary optimization ({n_generations} generations, pop={pop_size})...")
    best_score = population[0]["score"]

    for gen in range(n_generations):
        elite = population[:n_elite]
        candidates = list(elite)

        attempts = 0
        while len(candidates) < pop_size * 2 and attempts < pop_size * 20:
            attempts += 1
            parent = rng.choice(elite)["smiles"]

            if rng.random() < 0.35 and len(elite) >= 2:
                p2 = rng.choice(elite)["smiles"]
                child = crossover(parent, p2, rng)
            else:
                child = mutate(parent, rng)

            if child is None:
                continue
            r = score_molecule(child)
            if r:
                candidates.append(r[1])

        population = sorted(candidates, key=lambda x: x["score"], reverse=True)[:pop_size]
        best_score = population[0]["score"]

        if (gen + 1) % 10 == 0:
            logger.info(f"  Gen {gen+1:3d}: pop={len(candidates)} | best_score={best_score:.4f}")

    return population

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    out_dir = Path("data/results/CFTR_G85E")
    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("REFOLD Stage 3: Pharmacological Chaperone Generation")
    logger.info("Target: CFTR G85E — ICL1-TM2 junction pocket")
    logger.info(f"  fpocket druggability: {POCKET['druggability']:.3f}")
    logger.info(f"  Volume: {POCKET['volume']:.0f} Å³")
    logger.info(f"  Pocket residues: {POCKET['residues']}")
    logger.info(f"  Aromatic cage: F77-F78-F81 (3× Phe in TM2)")
    logger.info(f"  H-bond anchor: R74 (Arg)")
    logger.info("=" * 70)

    final_pop = run_evolution(SEEDS, n_generations=50, pop_size=40, n_elite=10)

    if not final_pop:
        logger.error("No candidates generated.")
        return

    logger.info("\nPhase 3: Final ranking and selection...")
    top10 = final_pop[:10]

    logger.info("\n" + "=" * 70)
    logger.info("TOP 10 PHARMACOLOGICAL CHAPERONE CANDIDATES:")
    logger.info("=" * 70)
    logger.info(f"{'Rank':>4} {'Score':>7} {'MW':>7} {'logP':>6} {'QED':>6} {'SA':>5} "
                f"{'HBD':>4} {'HBA':>4}  SMILES")
    logger.info("-" * 70)
    for i, cand in enumerate(top10):
        logger.info(
            f"{i+1:4d}  {cand['score']:7.4f}  {cand['mw']:7.1f}  {cand['logp']:6.2f}"
            f"  {cand['qed']:6.3f}  {cand['sa']:5.2f}  {cand['hbd']:4d}  {cand['hba']:4d}"
            f"  {cand['smiles']}"
        )

    best = top10[0]
    logger.info("\n" + "=" * 70)
    logger.info("BEST CANDIDATE — FULL REPORT:")
    logger.info(f"  SMILES:              {best['smiles']}")
    logger.info(f"  Composite score:     {best['score']:.4f}")
    logger.info(f"  MW:                  {best['mw']:.2f} Da")
    logger.info(f"  logP:                {best['logp']:.3f}")
    logger.info(f"  HBD / HBA:           {best['hbd']} / {best['hba']}")
    logger.info(f"  SA score:            {best['sa']:.2f} (1=easy, 10=hard)")
    logger.info(f"  QED:                 {best['qed']:.4f}")
    logger.info(f"  Pocket affinity:     {best['affinity']:.4f}")
    logger.info(f"  Lipinski:            True")
    logger.info(f"  Veber:               True")
    logger.info("=" * 70)

    results = {
        "target": {
            "mutation": "CFTR G85E",
            "uniprot": "P13569",
            "pocket_residues": POCKET["residues"],
            "pocket_druggability": POCKET["druggability"],
            "pocket_volume_A3": POCKET["volume"],
            "pocket_center": list(POCKET["center"]),
            "vx809_site_drug": 0.003,
            "vx809_site_status": "DESTROYED",
            "key_pharmacophore": "3x Phe (F77-F78-F81) aromatic cage + R74 H-bond anchor",
        },
        "top10_candidates": top10,
        "best_smiles": best["smiles"],
        "best_score": best["score"],
    }

    with open(out_dir / "stage3_chaperone_candidates.json", "w") as f:
        json.dump(results, f, indent=2)
    logger.info(f"Saved: {out_dir}/stage3_chaperone_candidates.json")


if __name__ == "__main__":
    main()
