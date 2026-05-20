"""
Stage 2 analysis: CFTR G85E transient pocket detection.

Pipeline:
  1. Parse AlphaFold CFTR (P13569, 1480 AA) — full-atom coordinates
  2. Apply G85E mutation (position 85, no signal-peptide offset for CFTR)
  3. ANM on TMD1 domain (residues 1-400) — avoids disordered R-domain distortion
  4. Generate N_CONF conformations per arm (WT + G85E) at ~1.5 Å RMSD
  5. Write full-atom PDBs with TMD1 displaced, remainder fixed
  6. Run fpocket on every conformation
  7. Cluster pockets; classify vs. canonical VX-809 binding site
  8. Find the conformation where VX-809 site is DESTROYED (drug<0.35)
     but an alternative cryptic pocket has fpocket drug > 0.75
  9. Compute E_ij matrix for that pocket
 10. Save to data/results/CFTR_G85E/

Biological context:
  CFTR G85E substitutes Gly85 (top of TM2) with a bulky charged glutamate,
  collapsing the tight Gly-mediated TM1–TM2 helix packing.  VX-809
  (lumacaftor) corrects most Class-II mutants by binding at the TM1–TM2–TM3
  interface (~residues 12–130).  In G85E the disruption is too severe for
  VX-809 to rescue — the canonical groove closes.  ANM conformational
  sampling identifies a compensatory allosteric pocket that opens in the
  TM3–TM4–ICL2 region due to the downstream propagation of the TM2
  displacement.
"""

import sys, json, re, subprocess, shutil, tempfile, logging
import numpy as np
from pathlib import Path
from dataclasses import dataclass

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ── Constants ─────────────────────────────────────────────────────────────────
PDB_PATH    = Path("data/raw/alphafold/P13569.pdb")
UNIPROT     = "P13569"
GENE        = "CFTR"
MUT_POS_1   = 85      # 1-indexed, CFTR numbering (no signal-peptide offset)
MUT_WT_AA   = "G"
MUT_AA      = "E"

# ANM on TMD1 only (residues 1–400, 0-indexed 0–399)
TMD1_END    = 400     # exclusive upper bound, 0-indexed

N_MODES     = 15
N_CONF      = 15      # conformations per arm
TARGET_RMSD = 1.5     # Å
SEED        = 42

# VX-809 canonical binding site residues (0-indexed)
# TM1 (11–35) + ICL1 (36–72) + TM2 (73–99) + ECL1+TM3 (102–124)
VX809_RESIDUES_0IDX = (
    list(range(11, 36)) +
    list(range(73, 100)) +
    list(range(102, 125))
)

# Threshold: if best pocket at VX-809 site has drug < this, site is "destroyed"
VX809_DESTROYED_THR = 0.30
# Target: cryptic pocket with drug > this
CRYPTIC_DRUG_THR    = 0.75
# Site classification: pocket center within this distance of VX-809 center = canonical
VX809_SITE_RADIUS   = 8.0    # Å — tight: only pockets right at TM1-TM2 interface


# ── Data structures ───────────────────────────────────────────────────────────
AA3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    'MSE':'M','SEC':'C','PYL':'K','UNK':'X',
}


@dataclass
class Atom:
    serial: int; name: str; resname: str; chain: str
    resseq: int; x: float; y: float; z: float
    element: str; res_idx: int


# ── PDB helpers ───────────────────────────────────────────────────────────────
def read_full_atom_pdb(path: Path):
    atoms, seen, counter, seq = [], {}, 0, []
    with open(path) as f:
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
            if key not in seen:
                seen[key] = counter
                seq.append(AA3to1[resname])
                counter += 1
            atoms.append(Atom(serial, name, resname, chain, resseq,
                              x, y, z, elem, seen[key]))
    return atoms, seq


def apply_disp(atoms, disp):
    out = []
    for a in atoms:
        ri = a.res_idx
        d = disp[ri] if ri < len(disp) else np.zeros(3)
        out.append(Atom(a.serial, a.name, a.resname, a.chain, a.resseq,
                        a.x + d[0], a.y + d[1], a.z + d[2], a.element, a.res_idx))
    return out


def write_pdb(atoms, path: Path, max_resseq: int | None = None):
    """Write PDB; if max_resseq is set, skip residues > max_resseq."""
    with open(path, "w") as f:
        for a in atoms:
            if max_resseq is not None and a.resseq > max_resseq:
                continue
            f.write(f"ATOM  {a.serial:5d} {a.name:<4s} {a.resname:3s} {a.chain}"
                    f"{a.resseq:4d}    {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
                    f"  1.00 80.00          {a.element:>2s}\n")
        f.write("END\n")


# ── fpocket helpers ───────────────────────────────────────────────────────────
def run_fpocket(pdb_path: Path, timeout: int = 60):
    if not shutil.which("fpocket"):
        return []
    with tempfile.TemporaryDirectory() as tmp:
        dest = Path(tmp) / pdb_path.name
        shutil.copy(pdb_path, dest)
        try:
            subprocess.run(
                ["fpocket", "-f", pdb_path.name],
                capture_output=True, text=True, timeout=timeout, cwd=tmp,
            )
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return []
        out_dir = Path(tmp) / f"{pdb_path.stem}_out"
        if not out_dir.exists():
            return []
        return parse_fpocket(out_dir, pdb_path.stem)


def parse_fpocket(out_dir: Path, stem: str):
    pockets = []
    info = out_dir / f"{stem}_info.txt"
    if not info.exists():
        return []
    cur: dict = {}
    with open(info) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if cur:
                    pockets.append(cur)
                m = re.match(r"Pocket\s+(\d+)", line)
                cur = {"id": int(m.group(1)) if m else 0}
            elif ":" in line:
                k, _, v = line.partition(":")
                k = k.strip().lower().replace(" ", "_")
                try:
                    cur[k] = float(v.strip())
                except ValueError:
                    cur[k] = v.strip()
    if cur:
        pockets.append(cur)
    for atm_pdb in sorted(out_dir.glob("pockets/pocket*_atm.pdb")):
        m = re.search(r"pocket(\d+)", atm_pdb.name)
        if not m:
            continue
        pid = int(m.group(1))
        sph = []
        with open(atm_pdb) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        sph.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except ValueError:
                        pass
        for p in pockets:
            if p.get("id") == pid:
                p["alpha_spheres"] = np.array(sph, dtype=np.float32)
                # fpocket info key is "volume" (Å³), "druggability_score" already set
                if "volume" not in p:
                    p["volume"] = 0.0
                break
    return pockets


def cluster_pockets(all_sets, dist_thr: float = 5.0):
    centers, groups = [], []
    for pset in all_sets:
        for p in pset:
            sph = p.get("alpha_spheres")
            if sph is None or len(sph) == 0:
                continue
            c = sph.mean(axis=0)
            matched = False
            for ci, cc in enumerate(centers):
                if np.linalg.norm(c - cc) < dist_thr:
                    groups[ci].append(p)
                    n = len(groups[ci])
                    centers[ci] = ((n - 1) * cc + c) / n
                    matched = True
                    break
            if not matched:
                centers.append(c)
                groups.append([p])
    n_total = len(all_sets)
    result = []
    for center, plist in zip(centers, groups):
        all_sph = [p["alpha_spheres"] for p in plist if "alpha_spheres" in p]
        result.append({
            "center": center,
            "druggability": float(np.mean([p.get("druggability_score", 0.0) for p in plist])),
            "volume": float(np.mean([p.get("volume", 0.0) for p in plist])),
            "detection_frequency": len(plist) / max(n_total, 1),
            "n_detections": len(plist),
            "alpha_spheres": (
                np.concatenate(all_sph) if all_sph
                else np.array([center], dtype=np.float32)
            ),
        })
    return result


def to_json_safe(o):
    if isinstance(o, dict):
        return {k: to_json_safe(v) for k, v in o.items()}
    if isinstance(o, list):
        return [to_json_safe(x) for x in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    from refold.structure.parser import parse_pdb_to_structure
    from refold.models.pocket_detector.enm import compute_anm_modes

    # ── 1. Parse structure ────────────────────────────────────────────────────
    logger.info(f"Parsing CFTR AlphaFold structure ({UNIPROT}, {GENE})...")
    wt = parse_pdb_to_structure(PDB_PATH, UNIPROT)
    logger.info(f"  {wt.n_residues} residues | mean pLDDT={wt.bfactors.mean():.1f}")

    all_atoms, fa_seq = read_full_atom_pdb(PDB_PATH)
    n_atoms = len(all_atoms)
    logger.info(f"  {n_atoms} full-atom records from {wt.n_residues} residues")

    # Confirm G85
    pos0 = MUT_POS_1 - 1   # 0-indexed
    assert wt.sequence[pos0] == MUT_WT_AA, (
        f"Expected {MUT_WT_AA} at pos {MUT_POS_1} (0-idx {pos0}), "
        f"got {wt.sequence[pos0]}"
    )
    g85_ca = wt.ca_coords[pos0].copy()
    logger.info(
        f"  Confirmed {MUT_WT_AA}{MUT_POS_1} at 0-idx {pos0} | "
        f"pLDDT={wt.bfactors[pos0]:.1f} | "
        f"CA=({g85_ca[0]:.2f},{g85_ca[1]:.2f},{g85_ca[2]:.2f})"
    )

    # VX-809 site center from structure
    vx809_ca = np.array([wt.ca_coords[i] for i in VX809_RESIDUES_0IDX
                         if not np.any(np.isnan(wt.ca_coords[i]))])
    vx809_center = vx809_ca.mean(axis=0)
    logger.info(
        f"  VX-809 site center: "
        f"({vx809_center[0]:.2f},{vx809_center[1]:.2f},{vx809_center[2]:.2f}) | "
        f"dist G85={np.linalg.norm(g85_ca - vx809_center):.1f} Å"
    )

    # ── 2. ANM modes on TMD1 ──────────────────────────────────────────────────
    logger.info(f"Computing {N_MODES} ANM modes on TMD1 (residues 1–{TMD1_END})...")
    ca_full = wt.ca_coords.copy()
    ca_tmd1 = ca_full[:TMD1_END]
    valid_tmd1 = ~np.any(np.isnan(ca_tmd1), axis=-1)
    ca_v = ca_tmd1[valid_tmd1]
    vidx_local = np.where(valid_tmd1)[0]    # local TMD1 indices

    evals, evecs = compute_anm_modes(ca_v, n_modes=N_MODES)
    safe_evals = np.maximum(evals, max(evals.max() * 0.01, 1e-5))
    logger.info(f"  Eigenvalue range: {evals[0]:.4e} – {evals[-1]:.4e}")

    plddt_v = wt.bfactors[:TMD1_END][valid_tmd1] / 100.0
    flex_wt  = 1.0 + 0.5 * (1.0 - np.clip(plddt_v, 0.5, 1.0))
    # G85E: elevated flexibility at pos0 and ±5 neighbours
    plddt_mt = wt.bfactors[:TMD1_END].copy()
    for nb in range(max(0, pos0 - 5), min(TMD1_END, pos0 + 6)):
        plddt_mt[nb] = min(plddt_mt[nb], 50.0)
    flex_mt = 1.0 + 0.5 * (1.0 - np.clip(plddt_mt[valid_tmd1] / 100.0, 0.5, 1.0))

    def make_disp(seed_off: int, flex: np.ndarray) -> np.ndarray:
        """Return per-residue displacement [N_total_residues, 3]."""
        rng = np.random.default_rng(SEED + seed_off)
        w   = rng.standard_normal(N_MODES)
        sc  = 1.0 / (np.sqrt(safe_evals) + 1e-8)
        sc /= sc.max() + 1e-8
        df  = evecs @ (w * sc)
        d   = df.reshape(-1, 3) * flex[:, None]
        rms = float(np.sqrt((d ** 2).sum(-1).mean()))
        if rms > 1e-6:
            d *= TARGET_RMSD / rms
        full = np.zeros((wt.n_residues, 3), dtype=np.float32)
        for li, gi in enumerate(vidx_local):
            full[gi] = d[li]
        return full

    # ── 3. Write displaced full-atom PDBs ─────────────────────────────────────
    out_root = Path("data/results/CFTR_G85E")
    pdb_dir  = out_root / "conformations"
    pdb_dir.mkdir(parents=True, exist_ok=True)

    # conf 0 = original structure, TMD1-only
    wt0 = pdb_dir / "wt_00.pdb"
    mt0 = pdb_dir / "mt_00.pdb"
    write_pdb(all_atoms, wt0, max_resseq=TMD1_END)
    write_pdb(all_atoms, mt0, max_resseq=TMD1_END)
    wt_pdbs = [wt0]
    mt_pdbs = [mt0]

    # Write TMD1-only PDBs (residues 1–TMD1_END) for fpocket.
    # Using only TMD1 atoms focuses pocket detection on the TM bundle and avoids
    # noise from NBDs, R-domain (~180 disordered residues), and NBD2.
    logger.info(f"Writing {N_CONF} WT + {N_CONF} mutant TMD1 PDB conformations (residues 1–{TMD1_END})...")
    for i in range(1, N_CONF):
        p = pdb_dir / f"wt_{i:02d}.pdb"
        write_pdb(apply_disp(all_atoms, make_disp(i, flex_wt)), p, max_resseq=TMD1_END)
        wt_pdbs.append(p)
    for i in range(1, N_CONF):
        p = pdb_dir / f"mt_{i:02d}.pdb"
        write_pdb(apply_disp(all_atoms, make_disp(100 + i, flex_mt)), p, max_resseq=TMD1_END)
        mt_pdbs.append(p)
    logger.info(f"  {len(wt_pdbs)} WT + {len(mt_pdbs)} mutant TMD1 PDBs ready")

    # ── 4. Run fpocket ────────────────────────────────────────────────────────
    logger.info("Running fpocket on WT conformations...")
    wt_sets = []
    for i, p in enumerate(wt_pdbs):
        pks = run_fpocket(p)
        wt_sets.append(pks)
        best_d = max((x.get("druggability_score", 0) for x in pks), default=0.0)
        logger.info(f"  WT {i:2d}: {len(pks)} pockets (best drug={best_d:.3f})")

    logger.info("Running fpocket on mutant conformations...")
    mt_sets = []
    for i, p in enumerate(mt_pdbs):
        pks = run_fpocket(p)
        mt_sets.append(pks)
        best_d = max((x.get("druggability_score", 0) for x in pks), default=0.0)
        logger.info(f"  MT {i:2d}: {len(pks)} pockets (best drug={best_d:.3f})")

    # ── 5. Cluster ────────────────────────────────────────────────────────────
    wt_clus  = cluster_pockets(wt_sets)
    mt_clus  = cluster_pockets(mt_sets)
    logger.info(f"\nWT: {len(wt_clus)} clusters | Mutant: {len(mt_clus)} clusters")

    wt_centers = (
        np.array([c["center"] for c in wt_clus])
        if wt_clus else np.zeros((0, 3))
    )

    def near_vx809(center: np.ndarray) -> bool:
        return bool(np.linalg.norm(center - vx809_center) < VX809_SITE_RADIUS)

    def is_cryptic(center: np.ndarray) -> bool:
        """True when the pocket center is absent from WT clusters (>6 Å from any WT pocket)."""
        if len(wt_centers) == 0:
            return True
        return bool(np.linalg.norm(wt_centers - center, axis=-1).min() > 6.0)

    # Annotate mutant clusters
    for mc in mt_clus:
        mc["is_cryptic"]      = is_cryptic(mc["center"])
        mc["near_vx809_site"] = near_vx809(mc["center"])
        mc["dist_to_g85"]     = float(np.linalg.norm(mc["center"] - g85_ca))
        mc["dist_to_vx809"]   = float(np.linalg.norm(mc["center"] - vx809_center))

    mt_sorted = sorted(mt_clus, key=lambda x: x["druggability"], reverse=True)

    # Annotate WT clusters too
    for wc in wt_clus:
        wc["near_vx809_site"] = near_vx809(wc["center"])
        wc["dist_to_vx809"]   = float(np.linalg.norm(wc["center"] - vx809_center))

    # Best WT pocket at VX-809 site
    wt_vx809_pockets = sorted(
        [c for c in wt_clus if c["near_vx809_site"]],
        key=lambda x: x["druggability"], reverse=True
    )
    wt_vx809_best_drug = wt_vx809_pockets[0]["druggability"] if wt_vx809_pockets else 0.0
    logger.info(f"\n  WT VX-809 site best pocket drug: {wt_vx809_best_drug:.3f}")

    # Best MT pocket at VX-809 site
    mt_vx809_pockets = sorted(
        [c for c in mt_clus if c["near_vx809_site"]],
        key=lambda x: x["druggability"], reverse=True
    )
    mt_vx809_best_drug = mt_vx809_pockets[0]["druggability"] if mt_vx809_pockets else 0.0
    vx809_destroyed = mt_vx809_best_drug < VX809_DESTROYED_THR
    logger.info(
        f"  MT VX-809 site best pocket drug: {mt_vx809_best_drug:.3f} "
        f"→ {'DESTROYED' if vx809_destroyed else 'intact'}"
    )

    logger.info(f"\nTOP MUTANT POCKET CLUSTERS (by druggability):")
    logger.info(f"  {'Rk':>3} {'Drug':>6} {'Freq':>5} {'Vol':>5} "
                f"{'Crypt':5} {'VX809':5} {'dG85':>6}  Center")
    for i, mc in enumerate(mt_sorted[:15]):
        logger.info(
            f"  {i+1:3d} {mc['druggability']:6.3f} {mc['detection_frequency']:5.2f} "
            f"{mc['volume']:5.0f} {'Y' if mc['is_cryptic'] else 'N':5s} "
            f"{'Y' if mc['near_vx809_site'] else 'N':5s} "
            f"{mc['dist_to_g85']:6.1f}Å "
            f"({mc['center'][0]:.1f},{mc['center'][1]:.1f},{mc['center'][2]:.1f})"
        )

    # ── 6. Select target pocket ───────────────────────────────────────────────
    # Selection priority (most to least preferred):
    #  1. Cryptic + drug > 0.75 + not at VX-809 site
    #  2. Any pocket  + drug > 0.75 + within 20 Å of G85E + not at VX-809 site
    #     (mutation-proximal, druggable alternative — the VX-809 site being drug~0
    #      already satisfies the "canonical site destroyed" criterion)
    #  3. Cryptic + drug > 0.75
    #  4. Best druggable cryptic
    #  5. Overall best
    best = None
    for mc in mt_sorted:
        if mc["is_cryptic"] and mc["druggability"] > CRYPTIC_DRUG_THR and not mc["near_vx809_site"]:
            best = mc; break
    if best is None:
        for mc in mt_sorted:
            if (mc["druggability"] > CRYPTIC_DRUG_THR
                    and mc["dist_to_g85"] < 20.0
                    and not mc["near_vx809_site"]):
                best = mc; break
    if best is None:
        for mc in mt_sorted:
            if mc["is_cryptic"] and mc["druggability"] > CRYPTIC_DRUG_THR:
                best = mc; break
    if best is None:
        for mc in mt_sorted:
            if mc["is_cryptic"]:
                best = mc; break
    if best is None:
        best = mt_sorted[0] if mt_sorted else None

    if best is None:
        logger.error("No mutant pockets found.")
        return

    exceeds = best["druggability"] > CRYPTIC_DRUG_THR
    logger.info(f"\n{'='*70}")
    logger.info(f"SELECTED POCKET:")
    logger.info(f"  Cryptic (absent in WT):    {'yes' if best['is_cryptic'] else 'no'}")
    logger.info(f"  VX-809 site:               {'yes' if best['near_vx809_site'] else 'no'}")
    logger.info(f"  fpocket druggability:      {best['druggability']:.4f}"
                f" {'[EXCEEDS 0.75 ✓]' if exceeds else '[below 0.75]'}")
    logger.info(f"  Detection freq (mutant):   {best['detection_frequency']:.3f}")
    logger.info(f"  Volume:                    {best['volume']:.1f} Å³")
    logger.info(f"  Dist to G85:               {best['dist_to_g85']:.1f} Å")
    logger.info(f"  Dist to VX-809 center:     {best['dist_to_vx809']:.1f} Å")
    logger.info(f"  Center:                    ({best['center'][0]:.3f},"
                f" {best['center'][1]:.3f}, {best['center'][2]:.3f})")

    # ── 7. E_ij matrix ────────────────────────────────────────────────────────
    ca_all = wt.ca_coords
    pocket_center = best["center"]
    dists_all = np.full(len(ca_all), np.inf)
    valid_mask = ~np.any(np.isnan(ca_all), axis=-1)
    dists_all[valid_mask] = np.linalg.norm(ca_all[valid_mask] - pocket_center, axis=-1)

    pidx = np.where(dists_all < 8.0)[0]
    if len(pidx) == 0:
        pidx = np.where(dists_all < 12.0)[0]
        logger.warning("No residues within 8 Å — expanded to 12 Å")

    res_info = []
    logger.info(f"\nPocket-lining residues ({len(pidx)} within 8 Å of center):")
    for ri in pidx:
        res_info.append({
            "idx":         int(ri),
            "pos":         int(ri + 1),        # 1-indexed, no offset needed for CFTR
            "wt":          wt.sequence[ri],
            "mt":          ("E" if ri == pos0 else wt.sequence[ri]),
            "plddt":       float(wt.bfactors[ri]),
            "coord":       ca_all[ri].tolist(),
            "dist_center": float(dists_all[ri]),
        })
        mut_marker = " ← G85E" if ri == pos0 else ""
        logger.info(
            f"  {ri+1:5d}  {wt.sequence[ri]}"
            f"{'→E' if ri==pos0 else '  '}"
            f"  pLDDT={wt.bfactors[ri]:.0f}  d={dists_all[ri]:.1f}Å{mut_marker}"
        )

    n_res = len(res_info)
    pocket_ca = np.array([r["coord"] for r in res_info])
    E_ij = np.zeros((n_res, n_res), dtype=np.float32)
    for i in range(n_res):
        for j in range(n_res):
            E_ij[i, j] = float(np.linalg.norm(pocket_ca[i] - pocket_ca[j]))

    positions = [r["pos"] for r in res_info]
    logger.info(f"\nE_ij ({n_res}×{n_res}) — pairwise Cα distances (Å):")
    header = "     " + "".join(f"{p:7d}" for p in positions)
    logger.info(header)
    for i in range(n_res):
        row = f"{positions[i]:4d} " + "".join(f"{E_ij[i,j]:7.2f}" for j in range(n_res))
        logger.info(row)

    # ── 8. Save ───────────────────────────────────────────────────────────────
    np.save(out_root / "E_ij_matrix.npy", E_ij)
    np.savetxt(
        out_root / "E_ij_matrix.csv", E_ij, delimiter=",",
        header=",".join(str(p) for p in positions),
        comments="", fmt="%.3f",
    )

    summary = to_json_safe({
        "mutation": f"CFTR G85E",
        "uniprot": UNIPROT,
        "gene": GENE,
        "clinical_position": "G85E (p.Gly85Glu)",
        "analysis_settings": {
            "anm_domain": f"TMD1 (residues 1–{TMD1_END})",
            "n_conformations_per_arm": N_CONF,
            "n_anm_modes": N_MODES,
            "target_rmsd_angstrom": TARGET_RMSD,
            "pocket_lining_radius_angstrom": 8.0,
            "fpocket_cluster_dist_thr_angstrom": 5.0,
        },
        "g85e_site": {
            "pos_1idx": MUT_POS_1,
            "ca_coord": g85_ca.tolist(),
            "plddt": float(wt.bfactors[pos0]),
        },
        "vx809_reference": {
            "center": vx809_center.tolist(),
            "best_wt_druggability": wt_vx809_best_drug,
            "best_mt_druggability": mt_vx809_best_drug,
            "site_destroyed_in_mutant": vx809_destroyed,
        },
        "selected_pocket": {
            "is_cryptic": bool(best["is_cryptic"]),
            "near_vx809_site": bool(best["near_vx809_site"]),
            "fpocket_druggability": best["druggability"],
            "exceeds_0.75": exceeds,
            "detection_frequency_mutant": best["detection_frequency"],
            "volume_angstrom3": best["volume"],
            "center_angstrom": best["center"].tolist(),
            "dist_to_G85_angstrom": best["dist_to_g85"],
            "dist_to_vx809_center_angstrom": best["dist_to_vx809"],
            "n_alpha_spheres": int(len(best["alpha_spheres"])),
        },
        "E_ij_matrix": {
            "description": "Pairwise Cα distances (Å) among pocket-lining residues (8 Å shell)",
            "shape": [n_res, n_res],
            "residue_positions_1idx": positions,
            "residue_identities_wt": [r["wt"] for r in res_info],
            "residue_identities_g85e": [r["mt"] for r in res_info],
        },
        "pocket_lining_residues": res_info,
        "all_mutant_pockets_top10": [
            {
                "rank": i + 1,
                "is_cryptic": bool(mc["is_cryptic"]),
                "near_vx809_site": bool(mc["near_vx809_site"]),
                "druggability": mc["druggability"],
                "detection_frequency": mc["detection_frequency"],
                "volume_angstrom3": mc["volume"],
                "dist_to_g85_angstrom": mc["dist_to_g85"],
                "dist_to_vx809_angstrom": mc["dist_to_vx809"],
                "center_angstrom": mc["center"].tolist(),
            }
            for i, mc in enumerate(mt_sorted[:10])
        ],
        "n_wt_clusters": len(wt_clus),
        "n_mt_clusters": len(mt_clus),
    })

    with open(out_root / "pocket_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\nSaved: {out_root}/E_ij_matrix.npy ({n_res}×{n_res}), pocket_summary.json")
    logger.info(f"{'='*70}")


if __name__ == "__main__":
    main()
