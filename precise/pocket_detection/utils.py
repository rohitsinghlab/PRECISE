from typing import List, Dict
import warnings
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import numpy as np
import torch
from sklearn.cluster import MeanShift
from torch_geometric.nn import radius


def parse_pdb_residues(pdb_path: str) -> List[Dict]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBConstructionWarning)

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)

    residues = []

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                if residue.get_id()[0] != " ":
                    continue

                residue_id = residue.get_id()[1]
                residue_name = residue.get_resname()
                try:
                    ca_atom = residue["CA"]
                except KeyError:
                    continue
                coords = np.array(ca_atom.get_coord())

                residue_info = {
                    "residue_id": residue_id,
                    "residue_name": residue_name,
                    "chain_id": chain_id,
                    "coordinates": coords,
                }

                residues.append(residue_info)

    return residues


def find_closest_residues(
    patch_centers: torch.Tensor, receptor_residues: List[Dict]
) -> List[int]:
    if isinstance(patch_centers, torch.Tensor):
        patch_centers = patch_centers.detach().cpu().numpy()
    patch_coords = np.asarray(patch_centers)
    if patch_coords.size == 0:
        return []

    residue_coords = np.array([res["coordinates"] for res in receptor_residues])

    closest_residues = []

    for patch_coord in patch_coords:
        distances = np.linalg.norm(residue_coords - patch_coord, axis=1)

        closest_idx = np.argmin(distances)
        closest_residue = receptor_residues[closest_idx].copy()
        closest_residue["distance"] = distances[closest_idx]

        closest_residues.append(closest_residue)

    return sorted({int(res["residue_id"]) for res in closest_residues})


def find_hotspot_vertices(
    vertex_positions: torch.Tensor,
    vertex_scores: np.ndarray,
    threshold: float = 0.00,
    bandwidth: float = 12.0,
    min_samples: int = 15,
) -> List[Dict]:
    if isinstance(vertex_positions, torch.Tensor):
        vertex_positions = vertex_positions.detach().cpu().numpy()

    high_score_mask = vertex_scores >= threshold
    if not np.any(high_score_mask):
        return []

    high_score_vertices = vertex_positions[high_score_mask]
    high_score_indices = np.where(high_score_mask)[0]

    clustering = MeanShift(bandwidth=bandwidth).fit(high_score_vertices)
    cluster_labels = clustering.labels_

    if len(cluster_labels) == 0:
        return []

    clusters: List[Dict] = []
    for cluster_label in np.unique(cluster_labels):
        if cluster_label == -1:
            continue
        cluster_vertices = high_score_vertices[cluster_labels == cluster_label]
        cluster_indices = high_score_indices[cluster_labels == cluster_label]
        clusters.append(
            {
                "vertices": cluster_vertices,
                "indices": cluster_indices,
                "average_score": vertex_scores[cluster_indices].mean(),
            }
        )

    if not clusters:
        return []

    enriched_clusters = []
    for cluster in clusters:
        vertices = cluster["vertices"]
        centroid = vertices.mean(axis=0) if len(vertices) > 0 else None
        enriched_clusters.append(
            {
                "vertices": vertices,
                "indices": cluster["indices"],
                "average_score": float(cluster["average_score"]),
                "size": int(len(vertices)),
                "centroid": centroid if centroid is None else centroid.astype(float),
            }
        )

    enriched_clusters.sort(key=lambda c: c["average_score"], reverse=True)
    return enriched_clusters


def accumulate_vertex_scores_via_radius(
    rec_centers: torch.Tensor,
    rec_pos: torch.Tensor,
    patch_scores: torch.Tensor,
    rec_patch_indices: torch.Tensor,
    radius_r: float = 7.0,
    max_num_neighbors: int = 128,
    threshold: float = 0.0,
    weight_mode: str = "gaussian",
    sigma: float = 1.5,
    agg: str = "mean",
    softmax_tau: float = 0.15,
) -> np.ndarray:
    device = rec_pos.device
    V = rec_pos.shape[0]

    keep = patch_scores >= threshold
    if keep.sum() == 0:
        return np.zeros(V, dtype=np.float32)

    kept_centers = rec_centers[keep]
    kept_scores = patch_scores[keep].clamp(0, 1)  # [K]
    kept_center_vertices = rec_patch_indices[keep]

    rows, cols = radius(
        kept_centers,
        rec_pos,
        radius_r,
        torch.zeros(kept_centers.shape[0], device=device),
        torch.zeros(rec_pos.shape[0], device=device),
        max_num_neighbors=max_num_neighbors,
    )

    K = kept_centers.shape[0]
    Vmax = max(1, V - 1)
    Kmax = max(1, K - 1)

    def fits_A(r, c, V, K):
        return r.numel() == 0 or (
            int(r.max().item()) <= Vmax and int(c.max().item()) <= Kmax
        )

    def fits_B(r, c, V, K):
        return r.numel() == 0 or (
            int(r.max().item()) <= Kmax and int(c.max().item()) <= Vmax
        )

    if fits_A(rows, cols, V, K):
        vert_idx = rows
        patch_local_idx = cols
    elif fits_B(rows, cols, V, K):
        vert_idx = cols
        patch_local_idx = rows
    else:
        vert_idx = rows
        patch_local_idx = cols

    v_scores = np.zeros(V, dtype=np.float32)
    if vert_idx.numel() == 0:
        vi_cent = kept_center_vertices.detach().cpu().numpy()
        ps_cent = kept_scores.detach().cpu().numpy()
        np.maximum.at(v_scores, vi_cent, ps_cent)
        return v_scores

    vi = vert_idx.detach().cpu().numpy()  # [E] vertex indices
    pi = patch_local_idx.detach().cpu().numpy()
    ps = kept_scores.detach().cpu().numpy()  # [K] patch scores
    ps_edge = ps[pi]  # [E] patch score per edge

    vp = rec_pos[vert_idx]  # [E, 3] vertex positions
    pc = kept_centers[patch_local_idx]  # [E, 3] patch centers
    d = torch.linalg.norm(vp - pc, dim=1).detach().cpu().numpy()

    if weight_mode == "gaussian":
        sigma = max(1e-6, float(sigma))
        w = np.exp(-0.5 * (d / sigma) ** 2).astype(np.float32)  # [E]
    elif weight_mode == "uniform":
        w = np.ones_like(d, dtype=np.float32)  # [E]
    else:
        raise ValueError(f"Unknown weight_mode: {weight_mode}")

    sum_w = np.zeros(V, dtype=np.float32)
    np.add.at(sum_w, vi, w)
    sum_w = np.maximum(sum_w, 1e-8)

    if agg == "mean":
        contrib = (w * ps_edge) / sum_w[vi]  # [E]
        np.add.at(v_scores, vi, contrib)

    center_vi = kept_center_vertices.detach().cpu().numpy()
    center_ps = kept_scores.detach().cpu().numpy()
    np.maximum.at(v_scores, center_vi, center_ps)

    return np.clip(v_scores, 0.0, 1.0)
