#!/usr/bin/env python3
"""
Stage 3 — Generic pharmacological chaperone generation.

Uses pharmacophore-guided SMILES evolution driven by pocket geometry
from pocket_summary.json produced by analyze_generic.py.

Usage:
    python3 scripts/generate_chaperone_generic.py \
        --gene GBA1 \
        --mutation L444P \
        --pocket-summary data/results/GBA1_L444P/pocket_summary.json \
        --outdir data/results/GBA1_L444P

Outputs:
    outdir/stage3_chaperone_candidates.json   (list, top candidate first)
"""

import sys, json, re, argparse, logging
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
log = logging.getLogger("chaperone_generic")

# ── Seed library: broadly active pharmacological chaperone scaffolds ──────────
# Covers iminosugar, aminoalcohol, benzimidazole, sulfonamide, pyrimidine
# classes with known track record in protein misfolding diseases.
SEED_SMILES = [
    # Iminosugar / aminosugar (GBA, GLA, HEXA, HEXB, GALC)
    "OC1CN[C@@H]2CCCC[C@@H]12",
    "NC1CC(O)C(O)C(O)C1",
    "OC[C@H]1CN[C@@H]2CC[C@@H]12",
    "OC1CNCC(O)C1",
    "N[C@@H]1CC[C@H](O)C1",
    # Aminoalcohol / amino-cyclohexane (CFTR, PAH, TTR)
    "NC1CCCCC1O",
    "OC1CNCCO1",
    "NC1CC(O)CCC1=O",
    "NC1CCCCC1=O",
    "OC1CC(N)CC1",
    # Benzimidazole / imidazole (TTR, SOD1, LRRK2)
    "Nc1nc2ccccc2[nH]1",
    "NCc1c[nH]c2ccccc12",
    "Cc1nc2ccc(N)cc2[nH]1",
    "Nc1ccnc(N)n1",
    "NC(=O)c1ccc(N)nc1",
    # Hydroxamic acid / amide (SERPINA1, MYH7)
    "NC(=O)c1ccc(O)cc1",
    "ONC(=O)c1cccc(F)c1",
    "NC(=O)c1ccncc1",
    # Aminopyrimidine / aminopyridine (KCNQ1, SCN5A, CFTR)
    "Nc1ccncc1",
    "Nc1ncc(F)cn1",
    "Nc1ccc(F)nc1",
    "CNc1ncnc2[nH]cnc12",
    # Sulfonamide (BRCA1, TP53, VHL)
    "NS(=O)(=O)c1ccc(N)cc1",
    "NS(=O)(=O)c1ccncc1",
    "NS(=O)(=O)c1ccc(O)cc1",
    # Small fragments for growth
    "NC1CCCC1",
    "N1CCOCC1",
    "NC(=O)c1ccccc1",
    "Nc1ccc(O)cc1",
    "NCc1ccc(O)cc1",
    "NC1CCC(F)CC1",
    # Hydroxypyridinone (iron-binding, applicable to metal-enzyme chaperones)
    "OC1=CC(=O)CC(O)=C1",
    "Oc1ccc(O)cc1C(=O)N",
    # Carbamate/urea fragments
    "NC(=O)Nc1ccccc1",
    "NC(=O)NC1CCCCC1",
]

# SMILES mutation operators (functional group modifications)
MUTATIONS = [
    lambda s: s.replace("c1ccccc1", "c1ccc(O)cc1", 1),
    lambda s: s.replace("c1ccccc1", "c1ccc(N)cc1", 1),
    lambda s: s.replace("c1ccccc1", "c1ccc(F)cc1", 1),
    lambda s: s.replace("c1ccccc1", "c1ccncc1", 1),
    lambda s: s.replace("c1ccccc1", "c1ccc(C(F)(F)F)cc1", 1),
    lambda s: s.replace("CC", "CC(=O)", 1),
    lambda s: s.replace("N", "NC", 1),
    lambda s: s.replace("O)", "OC)", 1),
    lambda s: s + "C(=O)N" if len(s) < 35 else s,
    lambda s: s + "N1CCCCC1" if len(s) < 30 else s,
    lambda s: s.replace("NC", "NCC", 1),
    lambda s: s.replace("OC", "OCC", 1),
    lambda s: s.replace("c1ccccc1", "c1ccc(Cl)cc1", 1),
    lambda s: s.replace("CC", "C(O)C", 1),
    lambda s: s.replace("N)", "NC(=O)", 1),
]


# ── RDKit scoring ──────────────────────────────────────────────────────────────

def canon(smi: str):
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except:
        return None


def score_molecule(smi: str, pocket: dict) -> tuple:
    """Returns (composite_score, props_dict). Returns (0.0, {}) on failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, QED, Fragments
        from rdkit.Chem.rdMolDescriptors import CalcNumHBD, CalcNumHBA

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return 0.0, {}

        mw    = Descriptors.ExactMolWt(mol)
        logp  = Descriptors.MolLogP(mol)
        hbd   = CalcNumHBD(mol)
        hba   = CalcNumHBA(mol)
        rotb  = rdMolDescriptors.CalcNumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        qed_v = QED.qed(mol)
        n_heavy = mol.GetNumHeavyAtoms()

        # SA score (synthetic accessibility, 1-10, lower=easier)
        try:
            from rdkit.Chem import RDConfig
            import os
            sa_path = os.path.join(RDConfig.RDContribDir, 'SA_Score', 'sascorer.py')
            import importlib.util
            spec = importlib.util.spec_from_file_location("sascorer", sa_path)
            sascorer = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(sascorer)
            sa = float(sascorer.calculateScore(mol))
        except:
            sa = 3.0  # default moderate SA score

        # ── Lipinski filters ─────────────────────────────────────────────────
        if mw > 500 or hbd > 5 or hba > 10 or logp > 5:
            return 0.0, {}
        if rotb > 10:
            return 0.0, {}
        if n_heavy < 5 or n_heavy > 50:
            return 0.0, {}
        if sa > 6.0:
            return 0.0, {}
        if qed_v < 0.15:
            return 0.0, {}

        # ── Pocket-specific pharmacophore score ──────────────────────────────
        pocket_vol   = pocket.get("volume_angstrom3", 400.0)
        pocket_drug  = pocket.get("fpocket_druggability", 0.8)
        pocket_res   = pocket.get("pocket_residues", {})

        # Count pharmacophore residues by type
        residue_seq = "".join(pocket_res.values()) if pocket_res else ""
        n_hydrophobic = sum(1 for r in residue_seq if r in "FYWLIMVA")
        n_pos_charged = sum(1 for r in residue_seq if r in "RKH")
        n_neg_charged = sum(1 for r in residue_seq if r in "DE")
        n_polar       = sum(1 for r in residue_seq if r in "STNQ")

        # Basic nitrogen score (needed for charged pockets with D/E)
        n_basic_N = (Fragments.fr_NH2(mol) + Fragments.fr_NH1(mol) + Fragments.fr_NH0(mol))
        basic_score = min(n_basic_N / max(n_neg_charged, 1), 1.0) if n_neg_charged > 0 else 0.3

        # Aromatic/hydrophobic matching
        n_arom = rdMolDescriptors.CalcNumAromaticRings(mol)
        hydrophobic_score = min((n_arom * 0.5 + max(logp, 0) * 0.1), 1.0) if n_hydrophobic > 2 else 0.3

        # H-bond donor score (for R/K-rich pockets)
        donor_score = min(hbd / max(n_pos_charged, 1), 1.0) if n_pos_charged > 0 else min(hbd / 3.0, 1.0)

        # Size matching: pocket volume ~0.7 × optimal MW
        optimal_mw = pocket_vol * 0.7
        size_score = max(0.0, 1.0 - abs(mw - optimal_mw) / max(optimal_mw, 1.0))

        # logP optimal range for druggability
        logp_score = max(0.0, 1.0 - abs(logp - 2.0) / 3.0)

        # ── Composite score ──────────────────────────────────────────────────
        composite = (
            0.30 * basic_score
            + 0.20 * hydrophobic_score
            + 0.15 * donor_score
            + 0.15 * size_score
            + 0.10 * logp_score
            + 0.05 * (1.0 - sa / 10.0)
            + 0.05 * qed_v
        )

        # Boost by pocket druggability
        composite *= (0.5 + 0.5 * pocket_drug)
        composite = float(np.clip(composite, 0.0, 1.0))

        props = {
            "smiles": smi,
            "composite_score": composite,
            "mw": round(mw, 2),
            "logp": round(logp, 3),
            "hbd": int(hbd),
            "hba": int(hba),
            "rotatable_bonds": int(rotb),
            "n_rings": int(rings),
            "qed": round(qed_v, 4),
            "sa_score": round(sa, 2),
            "n_heavy_atoms": int(n_heavy),
            "basic_n_score": round(basic_score, 3),
            "hydrophobic_score": round(hydrophobic_score, 3),
            "size_score": round(size_score, 3),
        }
        return composite, props
    except Exception as e:
        return 0.0, {}


# ── Evolutionary SMILES optimization ──────────────────────────────────────────

def mutate_smiles(smi: str, rng) -> str:
    from rdkit import Chem
    for _ in range(15):
        fn = rng.choice(MUTATIONS)
        try:
            new_smi = fn(smi)
            if new_smi == smi:
                continue
            mol = Chem.MolFromSmiles(new_smi)
            if mol is None:
                continue
            c = Chem.MolToSmiles(mol)
            if c and c != smi:
                return c
        except:
            pass
    return smi


def combine_smiles(smi_a: str, smi_b: str, rng) -> str:
    from rdkit import Chem
    linkers = ["C", "CC", "CCC", "NC", "OC", "c1ccc(cc1)"]
    for linker in rng.permutation(linkers):
        try:
            test = f"{smi_a}{linker}{smi_b}"
            mol = Chem.MolFromSmiles(test)
            if mol and 8 <= mol.GetNumHeavyAtoms() <= 45:
                return Chem.MolToSmiles(mol)
        except:
            pass
    return smi_a


def evolve(pocket: dict, n_gen: int = 40, pop_size: int = 40, elite_k: int = 10,
           seed: int = 42, rng=None) -> list:
    if rng is None:
        rng = np.random.default_rng(seed)

    log.info("Phase 1: Scoring seed library...")
    population = {}
    for smi in SEED_SMILES:
        c = canon(smi)
        if c is None:
            continue
        sc, props = score_molecule(c, pocket)
        if sc > 0:
            population[c] = (sc, props)
    log.info(f"  {len(population)} valid seeds")

    # Initial mutations of seeds
    seed_list = list(population.keys())
    for _ in range(150):
        parent = rng.choice(seed_list)
        child = mutate_smiles(parent, rng)
        if child and child not in population:
            sc, props = score_molecule(child, pocket)
            if sc > 0:
                population[child] = (sc, props)

    log.info(f"Phase 2: Evolutionary optimization ({n_gen} generations)...")
    for gen in range(n_gen):
        sorted_pop = sorted(population.items(), key=lambda x: x[1][0], reverse=True)
        elite = dict(sorted_pop[:elite_k])
        elite_list = list(elite.keys())

        new_gen = {}
        attempts = 0
        while len(new_gen) < pop_size - elite_k and attempts < 600:
            attempts += 1
            if rng.random() < 0.65 or len(elite_list) < 2:
                child = mutate_smiles(rng.choice(elite_list), rng)
            else:
                p1, p2 = rng.choice(elite_list, size=2, replace=False)
                child = combine_smiles(p1, p2, rng)

            c = canon(child) if child else None
            if c is None or c in population:
                continue
            sc, props = score_molecule(c, pocket)
            if sc > 0:
                new_gen[c] = (sc, props)

        population.update(elite)
        population.update(new_gen)

        if (gen + 1) % 10 == 0:
            best_score = max(v[0] for v in population.values()) if population else 0
            log.info(f"  Gen {gen+1:3d}: pop={len(population)} | best={best_score:.4f}")

    return sorted(population.items(), key=lambda x: x[1][0], reverse=True)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="REFOLD Stage 3 — generic chaperone generation")
    parser.add_argument("--gene",           required=True)
    parser.add_argument("--mutation",       required=True)
    parser.add_argument("--pocket-summary", required=True, type=Path)
    parser.add_argument("--outdir",         required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    # Load pocket summary
    with open(args.pocket_summary) as f:
        pocket_summary = json.load(f)

    pocket = pocket_summary.get("pocket", {})
    drug_score = pocket.get("fpocket_druggability", 0.8)
    volume     = pocket.get("volume_angstrom3", 400.0)

    log.info(f"Target: {args.gene} {args.mutation}")
    log.info(f"Pocket: druggability={drug_score:.3f}, volume={volume:.0f} Å³")
    log.info(f"Pocket residues: {list(pocket.get('pocket_residues', {}).values())}")

    # Deterministic seed from gene+mutation for reproducibility
    import hashlib
    seed_bytes = f"{args.gene}{args.mutation}".encode()
    seed_int   = int(hashlib.md5(seed_bytes).hexdigest()[:8], 16) % (2**31)
    rng        = np.random.default_rng(seed_int)

    # Run evolutionary optimization (20 gen × 40 pop keeps Stage 3 under ~3 min)
    final = evolve(pocket, n_gen=20, pop_size=40, elite_k=10, rng=rng)

    if not final:
        log.error("No valid candidates generated.")
        sys.exit(1)

    # Build top-10 output
    top10 = []
    for i, (smi, (score, props)) in enumerate(final[:10]):
        entry = {"rank": i + 1, **props}
        if "smiles" not in entry:
            entry["smiles"] = smi
        if "composite_score" not in entry:
            entry["composite_score"] = score
        top10.append(entry)

    best = top10[0]
    log.info(f"\nTop candidate: {best['smiles']}")
    log.info(f"  composite_score={best['composite_score']:.4f}")
    log.info(f"  MW={best['mw']:.1f}, logP={best['logp']:.2f}, QED={best['qed']:.3f}")

    # Save
    out = {
        "target":           f"{args.gene} {args.mutation}",
        "pocket_volume":    volume,
        "pocket_drug_score": drug_score,
        "generation_method": "pharmacophore-guided SMILES evolution (REFOLD Stage 3)",
        "n_generations":    40,
        "candidates":       top10,
    }
    out_path = args.outdir / "stage3_chaperone_candidates.json"
    with open(out_path, "w") as f:
        json.dump(top10, f, indent=2)  # daemon reads this as a list

    log.info(f"Saved {len(top10)} candidates → {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
