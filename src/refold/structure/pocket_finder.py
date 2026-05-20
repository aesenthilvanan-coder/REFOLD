"""Pocket detection using fpocket on conformational ensembles.

Detects transient pockets that appear in the mutant ensemble
but are absent in the wildtype — classified as TRANSIENT_MISFOLDING.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from refold.constants import (
    CA_IDX, FPOCKET_MIN_DRUGGABILITY_SCORE, POCKET_OVERLAP_THRESHOLD,
    MIN_POCKET_VOLUME, MAX_POCKET_VOLUME, TRANSIENT_POCKET_FREQ_THRESHOLD,
)
from refold.types import Pocket, PocketType

if TYPE_CHECKING:
    from refold.types import ProteinStructure

logger = logging.getLogger(__name__)


def _write_pdb_from_ca(ca_coords: np.ndarray, pdb_path: Path) -> None:
    """Write minimal PDB file with only Cα atoms."""
    with open(pdb_path, "w") as f:
        for i, pos in enumerate(ca_coords):
            if np.any(np.isnan(pos)):
                continue
            f.write(
                f"ATOM  {i+1:5d}  CA  ALA A{i+1:4d}    "
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                f"  1.00  0.00           C\n"
            )
        f.write("END\n")


def _parse_fpocket_output(out_dir: Path) -> list[dict]:
    """Parse fpocket output directory into pocket dicts."""
    pockets = []
    pdb_stem = out_dir.stem.removesuffix("_out")
    info_file = out_dir / f"{pdb_stem}_info.txt"

    if not info_file.exists():
        return pockets

    current_pocket: dict = {}
    with open(info_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("Pocket"):
                if current_pocket:
                    pockets.append(current_pocket)
                m = re.match(r"Pocket\s+(\d+)", line)
                current_pocket = {"pocket_id": int(m.group(1)) if m else 0}
            elif ":" in line:
                key, _, val = line.partition(":")
                key = key.strip().lower().replace(" ", "_")
                try:
                    current_pocket[key] = float(val.strip())
                except ValueError:
                    current_pocket[key] = val.strip()

    if current_pocket:
        pockets.append(current_pocket)

    # Parse alpha sphere coordinates from pdb files
    for pocket_dir in sorted(out_dir.glob("pockets/pocket*_atm.pdb")):
        m = re.search(r"pocket(\d+)", pocket_dir.name)
        if not m:
            continue
        pid = int(m.group(1))
        sphere_coords = []
        with open(pocket_dir) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        sphere_coords.append([x, y, z])
                    except ValueError:
                        continue
        for p in pockets:
            if p.get("pocket_id") == pid:
                p["alpha_spheres"] = np.array(sphere_coords, dtype=np.float32)
                break

    return pockets


def run_fpocket_on_conformation(
    ca_coords: np.ndarray,
    full_pdb_path: Path | None = None,
    min_volume: float = MIN_POCKET_VOLUME,
    max_volume: float = MAX_POCKET_VOLUME,
) -> list[dict]:
    """Run fpocket on a single conformation. Returns list of pocket dicts."""
    if shutil.which("fpocket") is None:
        logger.debug("fpocket binary not found — skipping pocket detection")
        return []

    with tempfile.TemporaryDirectory() as tmp:
        pdb_path = Path(tmp) / "conformation.pdb"
        if full_pdb_path and full_pdb_path.exists():
            shutil.copy(full_pdb_path, pdb_path)
        else:
            _write_pdb_from_ca(ca_coords, pdb_path)

        try:
            result = subprocess.run(
                ["fpocket", "-f", str(pdb_path)],
                capture_output=True, text=True, timeout=120, cwd=tmp,
            )
            if result.returncode != 0:
                logger.debug(f"fpocket failed: {result.stderr[:200]}")
                return []
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return []

        out_dir = Path(tmp) / "conformation_out"
        if not out_dir.exists():
            return []

        raw_pockets = _parse_fpocket_output(out_dir)

    filtered = []
    for p in raw_pockets:
        vol = p.get("volume_score", p.get("real_volume", 0.0))
        drug = p.get("druggability_score", 0.0)
        if min_volume <= vol <= max_volume:
            p["volume"] = vol
            p["druggability"] = drug
            filtered.append(p)

    return filtered


def _cluster_pockets_across_conformations(
    all_pocket_sets: list[list[dict]],
    distance_threshold: float = 4.0,
) -> list[dict]:
    """Greedy centroid-distance clustering across conformations.

    Returns list of cluster dicts with detection_frequency.
    """
    cluster_centers: list[np.ndarray] = []
    cluster_data: list[list[dict]] = []

    for conformation_pockets in all_pocket_sets:
        for pocket in conformation_pockets:
            spheres = pocket.get("alpha_spheres")
            if spheres is None or len(spheres) == 0:
                continue
            center = spheres.mean(axis=0)

            matched = False
            for ci, cc in enumerate(cluster_centers):
                dist = np.linalg.norm(center - cc)
                if dist < distance_threshold:
                    cluster_data[ci].append(pocket)
                    # Update running center
                    n = len(cluster_data[ci])
                    cluster_centers[ci] = ((n - 1) * cc + center) / n
                    matched = True
                    break

            if not matched:
                cluster_centers.append(center)
                cluster_data.append([pocket])

    n_conformations = len(all_pocket_sets)
    clusters = []
    for ci, (center, pockets) in enumerate(zip(cluster_centers, cluster_data)):
        detection_freq = len(pockets) / max(n_conformations, 1)
        avg_vol = np.mean([p.get("volume", 0.0) for p in pockets])
        avg_drug = np.mean([p.get("druggability", 0.0) for p in pockets])

        all_spheres = []
        for p in pockets:
            if p.get("alpha_spheres") is not None:
                all_spheres.append(p["alpha_spheres"])

        alpha_spheres = (
            np.concatenate(all_spheres, axis=0)
            if all_spheres
            else np.array([center], dtype=np.float32)
        )

        clusters.append({
            "center": center,
            "volume": avg_vol,
            "druggability": avg_drug,
            "detection_frequency": detection_freq,
            "alpha_spheres": alpha_spheres,
            "n_detections": len(pockets),
        })

    return clusters


def detect_transient_pockets(
    structure: "ProteinStructure",
    wt_conformations: list[np.ndarray],
    mutant_conformations: list[np.ndarray],
    full_pdb_path: Path | None = None,
    freq_threshold: float = TRANSIENT_POCKET_FREQ_THRESHOLD,
    druggability_threshold: float = FPOCKET_MIN_DRUGGABILITY_SCORE,
) -> list[Pocket]:
    """Detect transient pockets by comparing WT and mutant conformational ensembles.

    Pockets detected in the mutant ensemble but absent in WT
    are classified as TRANSIENT_MISFOLDING.
    """
    logger.info(f"Running fpocket on {len(wt_conformations)} WT + {len(mutant_conformations)} mutant conformations")

    wt_pocket_sets = []
    for ca in wt_conformations:
        pockets = run_fpocket_on_conformation(ca, full_pdb_path)
        wt_pocket_sets.append(pockets)

    mutant_pocket_sets = []
    for ca in mutant_conformations:
        pockets = run_fpocket_on_conformation(ca)
        mutant_pocket_sets.append(pockets)

    wt_clusters = _cluster_pockets_across_conformations(wt_pocket_sets)
    mutant_clusters = _cluster_pockets_across_conformations(mutant_pocket_sets)

    wt_centers = np.array([c["center"] for c in wt_clusters]) if wt_clusters else np.zeros((0, 3))

    result_pockets: list[Pocket] = []
    pocket_counter = 0

    for mc in mutant_clusters:
        if mc["detection_frequency"] < freq_threshold:
            continue
        if mc["druggability"] < druggability_threshold:
            continue

        center = mc["center"]
        is_transient = True
        if len(wt_centers) > 0:
            dists_to_wt = np.linalg.norm(wt_centers - center, axis=-1)
            if dists_to_wt.min() < 6.0:
                is_transient = False

        pocket_type = (
            PocketType.TRANSIENT_MISFOLDING if is_transient
            else PocketType.WILDTYPE_ALLOSTERIC
        )

        ca = structure.ca_coords
        valid_ca = ca[structure.residue_mask]
        if len(valid_ca) > 0:
            dists_to_ca = np.linalg.norm(valid_ca - center, axis=-1)
            residue_indices = np.where(dists_to_ca < 8.0)[0].tolist()
        else:
            residue_indices = []

        hydrophobicity = _compute_hydrophobicity(
            structure, residue_indices
        )

        pocket_counter += 1
        pocket = Pocket(
            pocket_id=f"pocket_{pocket_counter:03d}",
            pocket_type=pocket_type,
            center=center.astype(np.float32),
            volume=mc["volume"],
            druggability_score=mc["druggability"],
            hydrophobicity=hydrophobicity,
            residue_indices=residue_indices,
            alpha_sphere_coords=mc["alpha_spheres"].astype(np.float32),
            detection_frequency=mc["detection_frequency"],
            parent_structure=structure.uniprot_id,
        )
        result_pockets.append(pocket)

    n_transient = sum(1 for p in result_pockets if p.is_transient)
    logger.info(
        f"Detected {len(result_pockets)} pockets "
        f"({n_transient} transient/misfolding-specific)"
    )
    return result_pockets


def _compute_hydrophobicity(
    structure: "ProteinStructure",
    residue_indices: list[int],
    hydrophobicity_scale: dict[str, float] | None = None,
) -> float:
    """Compute mean hydrophobicity of pocket-lining residues."""
    from refold.constants import AA_PROPERTIES, STANDARD_AAS

    if not residue_indices:
        return 0.5

    if hydrophobicity_scale is None:
        hydrophobicity_scale = {
            aa: AA_PROPERTIES.get(one, {}).get("hydrophobicity", 0.0)
            for aa, one in zip(STANDARD_AAS, "ACDEFGHIKLMNPQRSTVWY")
        }

    scores = []
    for ri in residue_indices:
        if ri >= len(structure.sequence):
            continue
        aa_one = structure.sequence[ri]
        h = hydrophobicity_scale.get(aa_one, 0.0)
        h_norm = (h + 4.5) / 9.0
        scores.append(np.clip(h_norm, 0.0, 1.0))

    return float(np.mean(scores)) if scores else 0.5
