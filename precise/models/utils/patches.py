from __future__ import annotations

from typing import Tuple

import torch

from .geometry import radius_pytorch


__all__ = ["create_patches_vectorized", "aggregate_labels_to_patches"]


def create_patches_vectorized(
    node_features: torch.Tensor,
    node_original_features: torch.Tensor,
    positions: torch.Tensor,
    batch: torch.Tensor,
    n_patches: int = 512,
    radius: float = 7.0,
) -> Tuple[
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
    torch.Tensor,
]:
    device = node_features.device
    batch_size = batch.max().item() + 1 if batch.numel() > 0 else 0
    feat_dim = node_features.shape[-1]
    original_feat_dim = node_original_features.shape[-1]
    patch_features = torch.zeros(batch_size, n_patches, feat_dim, device=device)
    patch_original_features = torch.zeros(
        batch_size, n_patches, original_feat_dim, device=device
    )
    patch_centers = torch.zeros(batch_size, n_patches, 3, device=device)
    patch_masks = torch.zeros(batch_size, n_patches, dtype=torch.bool, device=device)
    patch_center_indices = torch.zeros(
        batch_size, n_patches, dtype=torch.long, device=device
    )

    for b in range(batch_size):
        batch_mask = batch == b
        batch_indices = torch.where(batch_mask)[0]
        n_nodes = batch_indices.shape[0]

        if n_nodes == 0:
            # Fill with zeros for empty batches
            patch_center_indices[b] = torch.zeros(
                n_patches, dtype=torch.long, device=device
            )
            continue

        batch_positions = positions[batch_indices]
        batch_features = node_features[batch_indices]
        batch_original_features = node_original_features[batch_indices]
        centers, selected_indices = _sample_patch_centers(batch_positions, n_patches)
        patch_centers[b] = centers

        # Store the center indices (mapped back to original node indices)
        center_indices_in_batch = batch_indices[selected_indices]
        patch_center_indices[b] = center_indices_in_batch

        search_radius = radius
        edge_index = radius_pytorch(x=batch_positions, y=centers, r=search_radius)

        local_features, local_original_features, local_mask = _aggregate_patch_features(
            batch_features,
            batch_original_features,
            edge_index,
            selected_indices,
            n_nodes,
            n_patches,
        )

        patch_features[b] = local_features
        patch_original_features[b] = local_original_features
        patch_masks[b] = local_mask

    return (
        patch_features,
        patch_original_features,
        patch_centers,
        patch_masks,
        patch_center_indices,
    )


def aggregate_labels_to_patches(
    node_labels: torch.Tensor,
    batch: torch.Tensor,
    patch_center_indices: torch.Tensor,
) -> torch.Tensor:
    batch_size = patch_center_indices.shape[0]
    n_patches = patch_center_indices.shape[1]
    device = node_labels.device

    patch_labels = torch.zeros(batch_size, n_patches, device=device)

    for b in range(batch_size):
        center_indices = patch_center_indices[b]
        # Get labels for the center nodes
        patch_labels[b] = node_labels[center_indices].float()

    return patch_labels


def _sample_patch_centers(
    positions: torch.Tensor, n_patches: int
) -> Tuple[torch.Tensor, torch.Tensor]:
    from torch_cluster import fps

    n_nodes = positions.size(0)
    actual_patches = min(n_patches, n_nodes)
    if n_nodes > 1:
        fps_ratio = min(1.0, actual_patches / n_nodes)
        indices = fps(positions, batch=None, ratio=fps_ratio, random_start=True)
    else:
        indices = torch.tensor([0], device=positions.device, dtype=torch.long)

    indices = _normalize_patch_count(indices, n_patches)
    centers = positions[indices]
    return centers, indices


def _normalize_patch_count(indices: torch.Tensor, n_patches: int) -> torch.Tensor:
    if len(indices) < n_patches:
        pad_size = n_patches - len(indices)
        if len(indices) > 0:
            pad_indices = indices[
                torch.randint(0, len(indices), (pad_size,), device=indices.device)
            ]
            indices = torch.cat([indices, pad_indices])
        else:
            indices = torch.zeros(n_patches, dtype=torch.long, device=indices.device)
    elif len(indices) > n_patches:
        indices = indices[:n_patches]
    return indices


def _aggregate_patch_features(
    batch_features: torch.Tensor,
    batch_original_features: torch.Tensor,
    edge_index: torch.Tensor,
    selected_indices: torch.Tensor,
    n_nodes: int,
    n_patches: int,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    feat_dim = batch_features.size(-1)
    original_feat_dim = batch_original_features.size(-1)
    device = batch_features.device

    patch_features = torch.zeros(n_patches, feat_dim, device=device)
    patch_original_features = torch.zeros(n_patches, original_feat_dim, device=device)
    patch_counts = torch.zeros(n_patches, device=device)

    if edge_index.numel() > 0:
        source_nodes = edge_index[0]
        target_patches = edge_index[1]

        patch_features.index_add_(0, target_patches, batch_features[source_nodes])
        patch_original_features.index_add_(
            0, target_patches, batch_original_features[source_nodes]
        )
        patch_counts.index_add_(
            0, target_patches, torch.ones_like(target_patches, dtype=torch.float)
        )

    mask = patch_counts > 0
    patch_mask = mask.clone()

    if mask.any():
        normalizer = patch_counts[mask].unsqueeze(-1)
        patch_features[mask] = patch_features[mask] / normalizer
        patch_original_features[mask] = patch_original_features[mask] / normalizer

    empty_indices = torch.where(~mask)[0]
    for idx in empty_indices:
        local_idx = selected_indices[idx] % n_nodes
        patch_features[idx] = batch_features[local_idx]
        patch_original_features[idx] = batch_original_features[local_idx]
        patch_counts[idx] = 1
        patch_mask[idx] = True

    return patch_features, patch_original_features, patch_mask
