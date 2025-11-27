from precise.dataset.precise_dataset import (
    smiles_to_morgan_fingerprint,
    build_single_example_batch,
)
from precise.pocket_detection.utils import parse_pdb_residues, find_closest_residues
from pathlib import Path
import torch
from precise.preprocess.msms_pipeline import process_protein
from precise.dataset.precise_dataset import get_data_from_ply
from precise.pocket_detection.utils import (
    find_hotspot_vertices,
    accumulate_vertex_scores_via_radius,
)


def inference(
    model,
    ligand_smi: str,
    receptor_pdb_path: str,
    quantize: bool = False,
    device: torch.device = torch.device("cpu"),
    bandwidth: float = 12.0,
    threshold: float = 0.35,
    num_patches: int = 2048,
) -> dict:
    payload = extract_patches(
        model, ligand_smi, receptor_pdb_path, device, bandwidth, threshold, num_patches
    )
    vertex_scores_d = payload["vertex_scores_d"]
    receptor_residues = parse_pdb_residues(receptor_pdb_path)

    cluster_key = "rec_clusters_qd" if quantize else "rec_clusters_d"
    clusters = payload.get(cluster_key, [])

    residue_clusters = []
    for cluster in clusters:
        residues = find_closest_residues(cluster["vertices"], receptor_residues)
        if not residues:
            continue
        residue_clusters.append(
            {
                "residue_ids": residues,
                "cluster_size": cluster["size"],
                "average_score": cluster["average_score"],
                "centroid": None
                if cluster["centroid"] is None
                else cluster["centroid"].astype(float).tolist(),
            }
        )

    primary_residues = residue_clusters[0]["residue_ids"] if residue_clusters else []

    return {
        "clusters": residue_clusters,
        "mode": "qd" if quantize else "d",
        "residue_ids": primary_residues,
        "vertex_scores_d": vertex_scores_d,
    }


def extract_patches(
    model,
    ligand_smi: str,
    receptor_pdb_path: str,
    device: torch.device,
    bandwidth: float = 12.0,
    threshold: float = 0.35,
    num_patches: int = 2048,
) -> dict:
    drug_fp = smiles_to_morgan_fingerprint(ligand_smi)
    surface_path = Path(receptor_pdb_path).parent / "surface.ply"
    process_protein(Path(receptor_pdb_path), "RECP", "A", surface_path)
    receptor_data = get_data_from_ply(surface_path)
    receptor_data = build_single_example_batch(receptor_data, drug_fp, device)
    payload = model.infer(receptor_data, quantize=False, num_patches=num_patches)

    rec_logits_qd = torch.sigmoid(payload["logits_qd"][0])  # [n_patches]
    rec_logits_d = torch.sigmoid(payload["logits_d"][0])  # [n_patches]
    rec_centers = payload["patch_centers"][0]  # [n_patches, 3]
    rec_patch_indices = payload["patch_indices"][0]  # [n_patches]

    vertex_scores_qd = accumulate_vertex_scores_via_radius(
        rec_centers, receptor_data.pos, rec_logits_qd, rec_patch_indices
    )
    vertex_scores_d = accumulate_vertex_scores_via_radius(
        rec_centers, receptor_data.pos, rec_logits_d, rec_patch_indices
    )

    hotspot_clusters_qd = find_hotspot_vertices(
        receptor_data.pos,
        vertex_scores_qd,
        bandwidth=bandwidth,
        threshold=threshold,
    )
    hotspot_clusters_d = find_hotspot_vertices(
        receptor_data.pos,
        vertex_scores_d,
        bandwidth=bandwidth,
        threshold=threshold,
    )

    return {
        "rec_clusters_qd": hotspot_clusters_qd,
        "rec_clusters_d": hotspot_clusters_d,
        "vertex_scores_d": vertex_scores_d,
        "rec": receptor_data,
    }
