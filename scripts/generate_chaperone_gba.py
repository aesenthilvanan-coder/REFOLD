"""
Stage 3: De novo pharmacological chaperone generation for GBA1 L444P.

Uses pharmacophore-guided SMILES evolution:
  1. Read MT conformation 10 pocket geometry (39 alpha spheres, drug=0.926)
  2. Infer pharmacophore from pocket-lining residues
  3. Evolve SMILES population using RDKit mutations
  4. Score each candidate with REFOLD pipeline:
       - Binding affinity heuristic (pocket-pharmacophore complementarity)
       - Lipinski RO5 + Veber
       - SA score, QED, PAINS filter
       - Composite rescue_probability weighting
  5. Report top-ranked SMILES with all scores

Pocket residues: A1, R2, P3, C4 | D24, S25, F26 | R48, M49, E50, L51
Pharmacophore:
  - Hydrophobic region: F26, L51, M49 → needs aromatic or aliphatic hydrophobic
  - H-bond acceptors: D24(OD), E50(OE) → basic N or OH on molecule
  - H-bond donors: R48 guanidinium → anionic group on molecule
  - Cysteine C4: potential disulfide/SH interaction → thiol/halide optional
"""

import sys, random, math, logging
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

SEED_SMILES = [
    # Iminosugar/aminosugar scaffolds (known GBA chaperone class)
    "OC1CN[C@@H]2CCCC[C@@H]12",           # isofagomine-like
    "NC1CC(O)C(O)C(O)C1",                 # amino-cyclohexane polyol
    "OC[C@H]1CN[C@@H]2CC[C@@H]12",        # bicyclic amino alcohol
    # Aminopiperidine / piperazine variants
    "N1CC(O)(CC1)c1ccccc1",               # 4-phenyl-4-OH-piperidine
    "NC1CCC(F)CC1",                        # fluorinated aminocyclohexane
    "OC1CNCC(O)C1",                        # aminocyclohexanediol
    "N1CCN(CC1)c1ccc(O)cc1",              # piperazine-phenol
    "NC1CC(O)CCC1=O",                      # amino-hydroxycyclohexanone
    # Aromatic amines with H-bond groups
    "Nc1ccc(O)cc1",                        # p-aminophenol
    "NCc1ccc(O)cc1",                       # tyramine
    "NCC(O)c1ccccc1",                      # 2-amino-1-phenylethan-1-ol
    "NC(CC(=O)O)Cc1ccccc1",               # phenylalanine-like
    "NC1CCCCC1O",                          # amino-hydroxycyclohexane
    "OC1CN(C)CC1",                         # N-methyl aminocyclopentanol
    # Indole/benzimidazole series
    "NCc1c[nH]c2ccccc12",                  # tryptamine
    "Nc1nc2ccccc2[nH]1",                   # 2-aminobenzimidazole
    "OCC(N)c1ccc(O)cc1",                   # norepinephrine-like
    # Small fragments for growth
    "NC1CCCC1",                             # aminocyclopentane
    "OC1CCCCC1N",                           # trans-2-aminocyclohexanol
    "NC(=O)c1ccc(O)cc1",                   # 4-hydroxybenzamide
    # Fluorinated aromatic amines (good for LLE)
    "Nc1ccc(F)cc1",                         # 4-fluoroaniline
    "NCc1cccc(F)c1",                        # 3-fluorobenzylamine
    "NC1CC(=O)CC1",                         # 3-aminocyclopentanone
    "NC1CCCCC1=O",                          # 2-aminocyclohexanone
    # Sulfonamide variants
    "NS(=O)(=O)c1ccc(O)cc1",               # sulfanilamide
    "NS(=O)(=O)c1ccncc1",                  # pyridine-4-sulfonamide
    # Morpholine / oxetane fragments
    "N1CCOCC1",                             # morpholine
    "C1COC1N",                              # 3-aminooxetane
    "OCC1CCN(C1)C",                         # N-methyl prolinal
    "OC1CNCCO1",                            # morpholine-diol
]

# Functional group mutations for SMILES evolution
MUTATIONS = [
    # Add hydroxyl
    lambda s: s.replace("c1ccccc1", "c1ccc(O)cc1", 1),
    lambda s: s.replace("C1CCCCC1", "C1CC(O)CCC1", 1),
    # Add amine
    lambda s: s.replace("c1ccccc1", "c1ccc(N)cc1", 1),
    # Add fluorine (lipophilicity modifier)
    lambda s: s.replace("c1ccccc1", "c1ccc(F)cc1", 1),
    # Add methyl
    lambda s: s.replace("N", "NC", 1),
    lambda s: s.replace("O)", "OC)", 1),
    # Add carbonyl
    lambda s: s.replace("CC", "CC(=O)", 1),
    # Fragment combination: add piperidine
    lambda s: s + "N1CCCCC1" if len(s) < 40 else s,
    # Fragment combination: add pyridine
    lambda s: s.replace("c1ccccc1", "c1ccncc1", 1),
    # Add amide
    lambda s: s.replace("NC", "NC(=O)C", 1),
]


def smiles_to_mol(smi):
    from rdkit import Chem
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None: return None
        return mol
    except: return None


def compute_pocket_affinity(smi, pocket_residues, pocket_volume=462.6):
    """
    Heuristic binding affinity score based on pharmacophore complementarity.

    Pocket residues: A, R, P, C, D, S, F, R, M, E, L
    Pharmacophore requirements:
      - Hydrophobic heavy atoms → complement F26/L51/M49/P3
      - H-bond donors (NH/OH) → complement D24/E50 carboxylates
      - H-bond acceptors (N/O) → complement R48 guanidinium
      - Basic N → salt bridge with D24(pKa~3.7) or E50(pKa~4.2)
      - Aromatic ring → pi-stacking with F26
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Fragments

    mol = smiles_to_mol(smi)
    if mol is None: return 0.0

    score = 0.0

    # Hydrophobic atoms (C not adjacent to polar)
    n_atoms = mol.GetNumHeavyAtoms()
    logp = Descriptors.MolLogP(mol)
    # Optimal logP for this pocket: 1-3 (moderately hydrophobic)
    logp_score = 1.0 - abs(logp - 2.0) / 3.0
    score += 0.25 * max(0, logp_score)

    # H-bond donors (interact with D24, E50 carboxylates)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hbd_score = min(hbd / 3.0, 1.0)  # up to 3 donors optimal
    score += 0.20 * hbd_score

    # H-bond acceptors (interact with R48 NH)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    hba_score = min(hba / 4.0, 1.0)
    score += 0.15 * hba_score

    # Basic nitrogen (salt bridge with D24/E50) — key pharmacophore
    n_basic = Fragments.fr_NH2(mol) + Fragments.fr_NH1(mol) + Fragments.fr_NH0(mol)
    basic_score = min(n_basic / 2.0, 1.0)
    score += 0.20 * basic_score

    # Aromatic ring (pi-stacking with F26)
    n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    arom_score = min(n_arom / 1.0, 1.0)
    score += 0.10 * arom_score

    # Size matching to pocket volume (462 Å³ → MW ~300-400)
    mw = Descriptors.ExactMolWt(mol)
    # Optimal MW = pocket_volume * 0.7 ≈ 324 Da
    optimal_mw = pocket_volume * 0.7
    size_score = 1.0 - abs(mw - optimal_mw) / optimal_mw
    score += 0.10 * max(0, size_score)

    return float(np.clip(score, 0.0, 1.0))


def refold_score(smi, pocket_residues):
    """Full REFOLD scoring: affinity + drug-likeness."""
    from refold.scoring.filters import (
        compute_lipinski_properties,
        compute_sa_score,
        compute_qed,
        check_pains,
    )

    try:
        lip = compute_lipinski_properties(smi)
        sa = compute_sa_score(smi)
        qed = compute_qed(smi)
        pains = check_pains(smi)
        affinity = compute_pocket_affinity(smi, pocket_residues)

        # Hard filters
        if not lip["passes_lipinski"]: return 0.0
        if not lip["passes_veber"]: return 0.0
        if sa > 5.0: return 0.0
        if pains: return 0.0
        if qed < 0.2: return 0.0

        # Composite score (matches REFOLD rescue_probability weights)
        composite = (
            0.30 * affinity        # pocket binding
            + 0.25 * affinity      # affinity (same heuristic, both components)
            + 0.15 * min(1.0, lip["logp"] / 3.0 + 0.3)   # druggability proxy
            + 0.10 * (1.0 - sa / 10.0)   # synthesizability
            + 0.10 * float(lip["passes_lipinski"])
            + 0.05 * qed
            + 0.05 * (1.0 - float(pains))
        )
        return float(composite), {
            "mw": lip["mw"], "logp": lip["logp"], "hbd": lip["hbd"],
            "hba": lip["hba"], "sa_score": sa, "qed": qed,
            "passes_lipinski": lip["passes_lipinski"],
            "passes_veber": lip["passes_veber"],
            "pains": pains, "affinity_score": affinity,
            "composite_score": composite,
        }
    except Exception as e:
        return 0.0, {}


def canonicalize(smi):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smi)
    if mol is None: return None
    return Chem.MolToSmiles(mol)


def mutate(smi, rng):
    """Apply a random SMILES mutation."""
    from rdkit import Chem
    for _ in range(10):
        fn = rng.choice(MUTATIONS)
        try:
            new_smi = fn(smi)
            if new_smi == smi: continue
            mol = Chem.MolFromSmiles(new_smi)
            if mol is None: continue
            canon = Chem.MolToSmiles(mol)
            if canon and canon != smi:
                return canon
        except: pass
    return smi


def combine(smi_a, smi_b, rng):
    """Fragment combination: take prefix of a + suffix of b."""
    from rdkit import Chem
    # Simple fragment joining via '.' (mixture) then sanitize
    combined = f"{smi_a}.{smi_b}"
    # Try to form a real bond using common linkers
    linkers = ["", "C", "CC", "CCC", "c1ccc(cc1)", "NC", "OC"]
    linker = str(rng.choice(linkers))
    # Just return valid if combinable
    for l in linkers:
        try:
            test = f"{smi_a}{l}{smi_b}"
            mol = Chem.MolFromSmiles(test)
            if mol and 10 <= mol.GetNumHeavyAtoms() <= 40:
                return Chem.MolToSmiles(mol)
        except: pass
    return smi_a


def main():
    from rdkit import Chem

    POCKET_RESIDUES = {
        "A40": "A", "R41": "R", "P42": "P", "C43": "C",
        "D63": "D", "S64": "S", "F65": "F",
        "R87": "R", "M88": "M", "E89": "E", "L90": "L",
    }
    POCKET_DRUGGABILITY = 0.926
    POCKET_VOLUME = 462.6

    logger.info("=" * 70)
    logger.info("REFOLD Stage 3: Pharmacological Chaperone Generation")
    logger.info(f"Target pocket: GBA1 L444P MT conformation 10")
    logger.info(f"Pocket druggability: {POCKET_DRUGGABILITY}")
    logger.info(f"Pocket volume: {POCKET_VOLUME:.0f} Å³")
    logger.info(f"Pocket-lining residues: {list(POCKET_RESIDUES.values())}")
    logger.info("=" * 70)

    rng = np.random.default_rng(42)
    random.seed(42)

    # ── Phase 1: Score seed library ────────────────────────────────────────
    logger.info("\nPhase 1: Scoring seed SMILES library...")
    population = {}  # canon_smi → (score, props)

    for smi in SEED_SMILES:
        canon = canonicalize(smi)
        if canon is None: continue
        result = refold_score(canon, POCKET_RESIDUES)
        if isinstance(result, tuple):
            score, props = result
        else:
            score, props = result, {}
        if score > 0:
            population[canon] = (score, props)

    logger.info(f"  {len(population)} valid seeds")

    # ── Phase 2: Evolutionary SMILES optimization ─────────────────────────
    N_GENERATIONS = 50
    POP_SIZE = 40
    ELITE_K = 10

    logger.info(f"\nPhase 2: Evolutionary optimization ({N_GENERATIONS} generations, pop={POP_SIZE})...")

    # Expand population with mutations of seeds
    all_smiles = list(population.keys())
    for _ in range(200):
        parent = rng.choice(all_smiles)
        child = mutate(parent, rng)
        if child and child not in population:
            result = refold_score(child, POCKET_RESIDUES)
            if isinstance(result, tuple) and result[0] > 0:
                population[child] = result

    for gen in range(N_GENERATIONS):
        # Sort by score
        sorted_pop = sorted(population.items(), key=lambda x: x[1][0], reverse=True)

        # Keep elite
        elite = dict(sorted_pop[:ELITE_K])

        # Generate new candidates from elite
        elite_smiles = list(elite.keys())
        new_gen = {}
        attempts = 0
        while len(new_gen) < POP_SIZE - ELITE_K and attempts < 500:
            attempts += 1
            op = rng.random()
            if op < 0.6 or len(elite_smiles) < 2:
                # Mutation
                parent = rng.choice(elite_smiles)
                child = mutate(parent, rng)
            else:
                # Combination
                p1, p2 = rng.choice(elite_smiles, size=2, replace=False)
                child = combine(p1, p2, rng)

            if child is None: continue
            canon = canonicalize(child)
            if canon is None or canon in population: continue

            result = refold_score(canon, POCKET_RESIDUES)
            if isinstance(result, tuple) and result[0] > 0:
                new_gen[canon] = result

        population.update(elite)
        population.update(new_gen)

        if (gen + 1) % 10 == 0:
            best_score = max(v[0] for v in population.values())
            logger.info(f"  Gen {gen+1:3d}: pop={len(population)} | best_score={best_score:.4f}")

    # ── Phase 3: Final ranking ─────────────────────────────────────────────
    logger.info("\nPhase 3: Final ranking and selection...")
    final = sorted(population.items(), key=lambda x: x[1][0], reverse=True)

    logger.info(f"\n{'='*70}")
    logger.info("TOP 10 PHARMACOLOGICAL CHAPERONE CANDIDATES:")
    logger.info(f"{'='*70}")
    logger.info(f"{'Rank':4s} {'Score':7s} {'MW':7s} {'logP':6s} {'QED':6s} {'SA':5s} {'HBD':4s} {'HBA':4s}  SMILES")
    logger.info("-" * 70)

    top10 = []
    for i, (smi, (score, props)) in enumerate(final[:10]):
        logger.info(
            f"  {i+1:2d}  {score:.4f}  {props.get('mw',0):6.1f}  "
            f"{props.get('logp',0):5.2f}  {props.get('qed',0):.3f}  "
            f"{props.get('sa_score',0):.2f}  {props.get('hbd',0):3d}  "
            f"{props.get('hba',0):3d}  {smi}"
        )
        top10.append({
            "rank": i + 1,
            "smiles": smi,
            "composite_score": score,
            **props
        })

    # ── Best candidate full analysis ───────────────────────────────────────
    best_smi, (best_score, best_props) = final[0]
    logger.info(f"\n{'='*70}")
    logger.info("BEST CANDIDATE — FULL REPORT:")
    logger.info(f"  SMILES:            {best_smi}")
    logger.info(f"  Composite score:   {best_score:.4f}")
    logger.info(f"  MW:                {best_props.get('mw',0):.2f} Da")
    logger.info(f"  logP:              {best_props.get('logp',0):.3f}")
    logger.info(f"  HBD / HBA:         {best_props.get('hbd',0)} / {best_props.get('hba',0)}")
    logger.info(f"  Rot. bonds:        {best_props.get('rotatable_bonds', 'N/A')}")
    logger.info(f"  SA score:          {best_props.get('sa_score',0):.2f} (1=easy, 10=hard)")
    logger.info(f"  QED:               {best_props.get('qed',0):.4f}")
    logger.info(f"  Passes Lipinski:   {best_props.get('passes_lipinski', False)}")
    logger.info(f"  Passes Veber:      {best_props.get('passes_veber', False)}")
    logger.info(f"  PAINS:             {best_props.get('pains', False)}")
    logger.info(f"  Pocket affinity:   {best_props.get('affinity_score',0):.4f}")
    logger.info(f"{'='*70}")

    # Save results
    import json
    out_dir = Path("data/results/GBA_L444P")
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(out_dir / "stage3_chaperone_candidates.json", "w") as f:
        json.dump({
            "target": "GBA1 L444P",
            "pocket": "MT conformation 10, Pocket 1 (fpocket drug=0.926)",
            "pocket_volume_angstrom3": POCKET_VOLUME,
            "pocket_residues": POCKET_RESIDUES,
            "generation_method": "pharmacophore-guided SMILES evolution",
            "n_generations": N_GENERATIONS,
            "top_smiles": best_smi,
            "top_score": best_score,
            "top10_candidates": top10,
        }, f, indent=2)

    logger.info(f"Saved: data/results/GBA_L444P/stage3_chaperone_candidates.json")
    return best_smi, best_score, top10


if __name__ == "__main__":
    main()
