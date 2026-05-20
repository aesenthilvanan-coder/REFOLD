"""
Stage 2 analysis: GBA L444P transient pocket detection.

Pipeline:
  1. Parse AlphaFold GBA (P04062) — full-atom coordinates
  2. Apply L444P (precursor position 483, mature position 444)
  3. Compute ANM modes with safe eigenvalue floor
  4. Generate 25 conformations per arm (WT + L444P mutant) at ~1.5 Å RMSD
  5. Run fpocket on full-atom displaced PDBs
  6. Cluster pockets, identify cryptic (mutant-only) and WT-overlap pockets
  7. Report best pocket, compute E_ij matrix (pairwise Cα distances)

Biological context:
  GBA1 L444P is the second most common Gaucher disease mutation.
  It affects the Ig-fold domain (domain 3) of GBA, causing ER-retention.
  Pharmacological chaperones (isofagomine, ambroxol) bind the active site
  (TIM barrel, domain 1), not directly at L444P. The relevant pocket is
  therefore the active site and/or allosteric pockets nearby.
"""

import sys, json, re, subprocess, shutil, tempfile
import logging
import numpy as np
from pathlib import Path
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class AtomRecord:
    serial: int; name: str; resname: str; chain: str
    resseq: int; x: float; y: float; z: float
    element: str; res_idx: int


AA3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','SEC':'C','PYL':'K','UNK':'X',
}


def read_full_atom_pdb(pdb_path):
    atoms, residues_seen, res_counter, seq = [], {}, 0, []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                resname = line[17:20].strip()
                if resname not in AA3to1:
                    continue
                serial = int(line[6:11])
                name   = line[12:16].strip()
                chain  = line[21]
                resseq = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                elem = line[76:78].strip() if len(line) >= 78 else name[0]
            except (ValueError, IndexError):
                continue
            key = (chain, resseq)
            if key not in residues_seen:
                residues_seen[key] = res_counter
                seq.append(AA3to1[resname])
                res_counter += 1
            atoms.append(AtomRecord(serial, name, resname, chain, resseq,
                                    x, y, z, elem, residues_seen[key]))
    return atoms, seq


def apply_disp(atoms, disp):
    out = []
    for a in atoms:
        ri = a.res_idx
        d = disp[ri] if ri < len(disp) else np.zeros(3)
        out.append(AtomRecord(a.serial, a.name, a.resname, a.chain, a.resseq,
                              a.x+d[0], a.y+d[1], a.z+d[2], a.element, a.res_idx))
    return out


def write_pdb(atoms, out_path):
    with open(out_path, "w") as f:
        for a in atoms:
            f.write(f"ATOM  {a.serial:5d} {a.name:<4s} {a.resname:3s} {a.chain}"
                    f"{a.resseq:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
                    f"  1.00 80.00          {a.element:>2s}\n")
        f.write("END\n")


def run_fpocket(pdb_path, timeout=45):
    if not shutil.which("fpocket"):
        return []
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp) / pdb_path.name
        shutil.copy(pdb_path, dest)
        try:
            subprocess.run(["fpocket", "-f", pdb_path.name],
                           capture_output=True, text=True, timeout=timeout, cwd=tmp)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return []
        out_dir = Path(tmp) / f"{pdb_path.stem}_out"
        if not out_dir.exists():
            return []
        return parse_fpocket(out_dir, pdb_path.stem)


def parse_fpocket(out_dir, pdb_stem):
    pockets = []
    info = out_dir / f"{pdb_stem}_info.txt"
    if not info.exists():
        return []
    cur = {}
    with open(info) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if cur: pockets.append(cur)
                m = re.match(r"Pocket\s+(\d+)", line)
                cur = {"id": int(m.group(1)) if m else 0}
            elif ":" in line:
                k, _, v = line.partition(":")
                k = k.strip().lower().replace(" ", "_")
                try:    cur[k] = float(v.strip())
                except: cur[k] = v.strip()
    if cur: pockets.append(cur)
    for atm in sorted(out_dir.glob("pockets/pocket*_atm.pdb")):
        m = re.search(r"pocket(\d+)", atm.name)
        if not m: continue
        pid = int(m.group(1))
        sph = []
        with open(atm) as f:
            for line in f:
                if line.startswith(("ATOM","HETATM")):
                    try: sph.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
                    except: pass
        for p in pockets:
            if p.get("id") == pid:
                p["alpha_spheres"] = np.array(sph, dtype=np.float32)
                break
    return pockets


def cluster_pockets(all_sets, dist_thr=5.0):
    centers, groups = [], []
    for pset in all_sets:
        for p in pset:
            sph = p.get("alpha_spheres")
            if sph is None or len(sph) == 0: continue
            c = sph.mean(axis=0)
            matched = False
            for ci, cc in enumerate(centers):
                if np.linalg.norm(c - cc) < dist_thr:
                    groups[ci].append(p)
                    n = len(groups[ci])
                    centers[ci] = ((n-1)*cc + c) / n
                    matched = True; break
            if not matched:
                centers.append(c); groups.append([p])
    n_total = len(all_sets)
    result = []
    for center, plist in zip(centers, groups):
        all_sph = [p["alpha_spheres"] for p in plist if "alpha_spheres" in p]
        result.append({
            "center": center,
            "druggability": float(np.mean([p.get("druggability_score",0) for p in plist])),
            "volume": float(np.mean([p.get("volume",0) for p in plist])),
            "detection_frequency": len(plist) / max(n_total, 1),
            "n_detections": len(plist),
            "alpha_spheres": np.concatenate(all_sph) if all_sph
                             else np.array([center], dtype=np.float32),
        })
    return result


def to_json_safe(o):
    """Recursively convert numpy types for JSON serialization."""
    if isinstance(o, dict):
        return {k: to_json_safe(v) for k, v in o.items()}
    if isinstance(o, list):
        return [to_json_safe(x) for x in o]
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.floating): return float(o)
    if isinstance(o, np.ndarray): return o.tolist()
    return o


def main():
    from refold.structure.parser import parse_pdb_to_structure
    from refold.structure.mutator import apply_mutation
    from refold.models.pocket_detector.enm import compute_anm_modes
    from refold.types import Mutation

    PDB = Path("data/raw/alphafold/P04062.pdb")

    # ── 1. Parse ──────────────────────────────────────────────────────────────
    logger.info("Parsing GBA AlphaFold structure (P04062)...")
    wt = parse_pdb_to_structure(PDB, "P04062")
    logger.info(f"  {wt.n_residues} residues | mean pLDDT={wt.bfactors.mean():.1f}")

    all_atoms, fa_seq = read_full_atom_pdb(PDB)
    logger.info(f"  {len(all_atoms)} full-atom records")

    # L444P: mature pos 444 = precursor pos 483
    POS = 483        # 1-indexed, precursor
    pos = POS - 1    # 0-indexed
    assert wt.sequence[pos] == "L", f"Expected L at {POS}, got {wt.sequence[pos]}"
    logger.info(f"  L444 confirmed at precursor pos {POS} (mature pos 444)")

    # Mutation site CA coordinate
    ca483 = wt.ca_coords[pos]
    logger.info(f"  L483 CA coord: ({ca483[0]:.1f}, {ca483[1]:.1f}, {ca483[2]:.1f})")

    # ── 2. ANM modes ──────────────────────────────────────────────────────────
    N_MODES, N_CONF, TARGET_RMSD = 20, 25, 1.5
    SEED = 42

    logger.info(f"Computing {N_MODES} ANM modes for {wt.n_residues}-residue GBA...")
    ca = wt.ca_coords.copy()
    valid = ~np.any(np.isnan(ca), axis=-1)
    ca_v = ca[valid]; vidx = np.where(valid)[0]

    evals, evecs = compute_anm_modes(ca_v, n_modes=N_MODES)
    safe_evals = np.maximum(evals, max(evals.max() * 0.01, 1e-5))
    logger.info(f"  Eigenvalues: {evals[0]:.2e} – {evals[-1]:.2e} | floor: {safe_evals[0]:.2e}")

    plddt = wt.bfactors[vidx] / 100.0
    flex = 1.0 + 0.5 * (1.0 - np.clip(plddt, 0.5, 1.0))  # [1.0, 1.25]

    def make_disp(seed_off, flex_arr):
        rng = np.random.default_rng(SEED + seed_off)
        w = rng.standard_normal(N_MODES)
        sc = 1.0 / (np.sqrt(safe_evals) + 1e-8)
        sc /= sc.max() + 1e-8
        df = evecs @ (w * sc)
        d = df.reshape(-1, 3)
        rms = float(np.sqrt((d**2).sum(-1).mean()))
        if rms > 1e-6:
            d *= TARGET_RMSD / rms
        d *= flex_arr[:, None]
        full = np.zeros((wt.n_residues, 3), dtype=np.float32)
        for li, gi in enumerate(vidx):
            full[gi] = d[li]
        return full

    # Mutant flexibility: Pro at 483 → increased local flexibility surroundings
    mut_plddt = wt.bfactors.copy(); mut_plddt[pos] = 35.0
    mut_flex = 1.0 + 0.5 * (1.0 - np.clip(mut_plddt[vidx]/100.0, 0.5, 1.0))

    # ── 3. Write PDBs ─────────────────────────────────────────────────────────
    pdb_dir = Path("data/results/GBA_L444P/conformations")
    pdb_dir.mkdir(parents=True, exist_ok=True)

    wt_pdbs = [PDB]   # conf 0 = original
    for i in range(1, N_CONF):
        p = pdb_dir / f"wt_{i:02d}.pdb"
        write_pdb(apply_disp(all_atoms, make_disp(i, flex)), p)
        wt_pdbs.append(p)

    mt_pdbs = [PDB]   # conf 0 = original (same backbone, P substitution tracked in sequence only)
    for i in range(1, N_CONF):
        p = pdb_dir / f"mt_{i:02d}.pdb"
        write_pdb(apply_disp(all_atoms, make_disp(100 + i, mut_flex)), p)
        mt_pdbs.append(p)

    logger.info(f"  Wrote {len(wt_pdbs)} WT + {len(mt_pdbs)} mutant full-atom PDBs")

    # ── 4. fpocket ────────────────────────────────────────────────────────────
    logger.info("Running fpocket on all conformations...")
    wt_sets, mt_sets = [], []
    for i, p in enumerate(wt_pdbs):
        pks = run_fpocket(p)
        wt_sets.append(pks)
        if pks:
            best = max(p.get("druggability_score", 0) for p in pks)
            logger.info(f"  WT {i:2d}/{N_CONF-1}: {len(pks)} pockets | best drug={best:.3f}")
    for i, p in enumerate(mt_pdbs):
        pks = run_fpocket(p)
        mt_sets.append(pks)
        if pks:
            best = max(p.get("druggability_score", 0) for p in pks)
            logger.info(f"  MT {i:2d}/{N_CONF-1}: {len(pks)} pockets | best drug={best:.3f}")

    # ── 5. Cluster ────────────────────────────────────────────────────────────
    wt_clus = cluster_pockets(wt_sets)
    mt_clus = cluster_pockets(mt_sets)
    wt_centers = np.array([c["center"] for c in wt_clus]) if wt_clus else np.zeros((0, 3))

    logger.info(f"\nWT: {len(wt_clus)} clusters | Mutant: {len(mt_clus)} clusters")

    def is_cryptic(center):
        if len(wt_centers) == 0: return True
        return float(np.linalg.norm(wt_centers - center, axis=-1).min()) > 6.0

    mt_sorted = sorted(mt_clus, key=lambda x: x["druggability"], reverse=True)

    # Annotate clusters
    for mc in mt_sorted:
        mc["cryptic"] = is_cryptic(mc["center"])
        mc["dist_to_L483"] = float(np.linalg.norm(mc["center"] - ca483))

    logger.info("\nMUTANT POCKET CLUSTERS (top 15 by druggability):")
    logger.info(f"  {'Rank':4s} {'Drug':6s} {'Freq':5s} {'Vol':5s} {'Type':12s} {'d_L483':7s} Center")
    for i, mc in enumerate(mt_sorted[:15]):
        tag = "CRYPTIC" if mc["cryptic"] else "WT_OVERLAP"
        logger.info(f"  {i+1:4d} {mc['druggability']:.3f}  {mc['detection_frequency']:.2f}"
                    f"  {mc['volume']:5.0f} {tag:12s} {mc['dist_to_L483']:6.1f}Å"
                    f" ({mc['center'][0]:.1f},{mc['center'][1]:.1f},{mc['center'][2]:.1f})")

    # ── 6. Select best pocket ─────────────────────────────────────────────────
    # Prefer cryptic with drug > 0.75, else best cryptic, else best overall near L483
    best = None
    for mc in mt_sorted:
        if mc["cryptic"] and mc["druggability"] > 0.75:
            best = mc; break
    if best is None:
        for mc in mt_sorted:
            if mc["cryptic"]:
                best = mc; break
    if best is None:
        best = mt_sorted[0] if mt_sorted else None

    if best is None:
        logger.error("No pockets detected."); return

    exceeds = best["druggability"] > 0.75
    logger.info(f"\n{'='*70}")
    logger.info(f"BEST POCKET (cryptic={'yes' if best['cryptic'] else 'no'}):")
    logger.info(f"  fpocket druggability: {best['druggability']:.4f}"
                f" {'[EXCEEDS 0.75]' if exceeds else '[below 0.75]'}")
    logger.info(f"  Detection freq mutant ensemble: {best['detection_frequency']:.3f}")
    logger.info(f"  Volume: {best['volume']:.1f} Å³")
    logger.info(f"  Distance to L483 (mutation site): {best['dist_to_L483']:.1f} Å")
    logger.info(f"  Center: ({best['center'][0]:.3f}, {best['center'][1]:.3f}, {best['center'][2]:.3f})")

    # Also find best WT-overlap pocket near L483 (for reference)
    wt_overlap_near = sorted(
        [mc for mc in mt_sorted if not mc["cryptic"]],
        key=lambda x: x["druggability"], reverse=True
    )
    if wt_overlap_near:
        ref = wt_overlap_near[0]
        logger.info(f"\n  Reference WT-overlap pocket (best drug): drug={ref['druggability']:.3f}"
                    f" dist_L483={ref['dist_to_L483']:.1f}Å")

    # ── 7. E_ij matrix ────────────────────────────────────────────────────────
    ca = wt.ca_coords
    center = best["center"]
    valid_mask = wt.residue_mask
    dists = np.full(len(ca), np.inf)
    dists[valid_mask] = np.linalg.norm(ca[valid_mask] - center, axis=-1)
    pidx = np.where(dists < 8.0)[0]

    SIGNAL_PEPTIDE = 39   # GBA signal peptide length

    res_info = []
    for ri in pidx:
        mature = ri + 1 - SIGNAL_PEPTIDE if ri >= SIGNAL_PEPTIDE else ri + 1
        res_info.append({
            "idx": int(ri), "pos_pre": int(ri+1), "pos_mat": int(mature),
            "wt": wt.sequence[ri],
            "mt": ("P" if ri == pos else wt.sequence[ri]),
            "plddt": float(wt.bfactors[ri]),
            "coord": ca[ri].tolist(),
            "dist_to_center": float(dists[ri]),
        })
        marker = " ← L444P" if ri == pos else ""
        logger.info(f"  {ri+1:4d} (mat {mature:4d}) {wt.sequence[ri]}"
                    f"{'→P' if ri==pos else '  '}"
                    f"  pLDDT={wt.bfactors[ri]:.0f}  d={dists[ri]:.1f}Å{marker}")

    n_res = len(res_info)
    if n_res == 0:
        logger.warning("No pocket-lining residues found within 8 Å. Using 12 Å radius.")
        pidx = np.where(dists < 12.0)[0]
        for ri in pidx:
            mature = ri + 1 - SIGNAL_PEPTIDE if ri >= SIGNAL_PEPTIDE else ri + 1
            res_info.append({"idx": int(ri), "pos_pre": int(ri+1), "pos_mat": int(mature),
                             "wt": wt.sequence[ri],
                             "mt": "P" if ri == pos else wt.sequence[ri],
                             "plddt": float(wt.bfactors[ri]),
                             "coord": ca[ri].tolist(), "dist_to_center": float(dists[ri])})
        n_res = len(res_info)

    pocket_ca = np.array([r["coord"] for r in res_info])
    E_ij = np.zeros((n_res, n_res), dtype=np.float32)
    for i in range(n_res):
        for j in range(n_res):
            E_ij[i, j] = float(np.linalg.norm(pocket_ca[i] - pocket_ca[j]))

    logger.info(f"\n{'='*70}")
    logger.info(f"E_ij MATRIX — pairwise Cα distances (Å) — {n_res}×{n_res} residues")
    pre_pos = [r["pos_pre"] for r in res_info]
    mat_pos = [r["pos_mat"] for r in res_info]
    logger.info(f"Precursor positions: {pre_pos}")
    logger.info(f"Mature positions:    {mat_pos}")
    header = "     " + "".join(f"{p:7d}" for p in pre_pos)
    logger.info(header)
    for i, ri in enumerate(res_info):
        row = f"{ri['pos_pre']:4d} " + "".join(f"{E_ij[i,j]:7.2f}" for j in range(n_res))
        logger.info(row)

    # ── 8. Save ───────────────────────────────────────────────────────────────
    out = Path("data/results/GBA_L444P")
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "E_ij_matrix.npy", E_ij)
    np.savetxt(out / "E_ij_matrix.csv", E_ij, delimiter=",",
               header=",".join(str(r["pos_pre"]) for r in res_info),
               comments="", fmt="%.3f")

    summary = to_json_safe({
        "mutation": "GBA1 L444P",
        "uniprot": "P04062",
        "clinical_position": "L444P (mature protein, p.Leu444Pro)",
        "precursor_position": "L483P (full 536-AA precursor including signal peptide)",
        "analysis_settings": {
            "n_conformations_per_arm": N_CONF,
            "n_anm_modes": N_MODES,
            "target_rmsd_angstrom": TARGET_RMSD,
            "pocket_lining_radius_angstrom": 8.0,
            "fpocket_cluster_dist_thr_angstrom": 5.0,
        },
        "mutation_site_ca_coord": ca483.tolist(),
        "best_pocket": {
            "is_cryptic": bool(best["cryptic"]),
            "fpocket_druggability_score": best["druggability"],
            "exceeds_0.75_threshold": bool(exceeds),
            "detection_frequency_mutant": best["detection_frequency"],
            "volume_angstrom3": best["volume"],
            "center_angstrom": best["center"].tolist(),
            "dist_to_L483_angstrom": best["dist_to_L483"],
            "n_alpha_spheres": int(len(best["alpha_spheres"])),
        },
        "E_ij_matrix": {
            "description": "Pairwise Cα distances (Å) among pocket-lining residues (8 Å shell)",
            "shape": [n_res, n_res],
            "residue_positions_precursor": pre_pos,
            "residue_positions_mature": mat_pos,
            "residue_identities_wt": [r["wt"] for r in res_info],
            "residue_identities_mt": [r["mt"] for r in res_info],
            "files": {
                "npy": str(out / "E_ij_matrix.npy"),
                "csv": str(out / "E_ij_matrix.csv"),
            },
        },
        "pocket_lining_residues": res_info,
        "all_mutant_pockets_top10": [
            {
                "rank": i+1,
                "is_cryptic": bool(mc["cryptic"]),
                "druggability": mc["druggability"],
                "detection_frequency": mc["detection_frequency"],
                "volume_angstrom3": mc["volume"],
                "dist_to_L483_angstrom": mc["dist_to_L483"],
                "center_angstrom": mc["center"].tolist(),
            }
            for i, mc in enumerate(mt_sorted[:10])
        ],
        "n_wt_clusters": len(wt_clus),
        "n_mutant_clusters": len(mt_clus),
        "stage3_smiles": None,
        "stage3_note": (
            "DDPM molecule generator requires trained checkpoint weights. "
            "Run: make train-generator && "
            "make run-single UNIPROT=P04062 POS=444 WT=L MT=P GENE=GBA1"
        ),
    })

    with open(out / "pocket_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\nSaved: {out}/E_ij_matrix.npy ({n_res}×{n_res}), pocket_summary.json")

    logger.info(f"\n{'='*70}")
    logger.info("STAGE 3 NOTE:")
    logger.info("  De novo SMILES generation requires trained DDPM checkpoint.")
    logger.info("  Known GBA pharmacological chaperones for reference:")
    logger.info("    Isofagomine:  OC1CN[C@@H]2CCCC[C@@H]12  (active site binder)")
    logger.info("    Ambroxol:     OC1=CC(=CC(Br)=C1)CNC(=O)C2=CC=CC=C2Br")
    logger.info("    N-octyl-DNJ:  OC1CN[C@@H]2C[C@@H](O)[C@H](O)[C@@H]2O")
    logger.info("  To generate de novo candidates: make train-generator && make run-single")
    logger.info(f"{'='*70}")

    return E_ij, summary


if __name__ == "__main__":
    main()
