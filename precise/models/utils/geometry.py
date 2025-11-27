import torch


__all__ = ["radius_pytorch"]


def radius_pytorch(
    x: torch.Tensor,
    y: torch.Tensor,
    r: float,
    max_num_neighbors: int | None = None,
) -> torch.Tensor:
    dist_matrix = torch.cdist(x, y)
    mask = dist_matrix <= r
    edge_index = mask.nonzero().t()

    if max_num_neighbors is not None and edge_index.numel() > 0:
        edge_index, _ = _limit_neighbors(edge_index, max_num_neighbors)

    return edge_index


def _limit_neighbors(
    edge_index: torch.Tensor, max_num_neighbors: int
) -> tuple[torch.Tensor, torch.Tensor]:
    target_nodes = edge_index[1]
    sorted_indices = torch.argsort(target_nodes)
    edge_index_sorted = edge_index[:, sorted_indices]

    kept_mask = torch.ones(
        edge_index_sorted.size(1), dtype=torch.bool, device=edge_index.device
    )
    counts = torch.zeros(
        target_nodes.numel() if target_nodes.numel() > 0 else 0,
        dtype=torch.long,
        device=edge_index.device,
    )

    current_target = None
    current_count = 0
    for idx in range(edge_index_sorted.size(1)):
        target = edge_index_sorted[1, idx].item()
        if target != current_target:
            current_target = target
            current_count = 0
        if current_count >= max_num_neighbors:
            kept_mask[idx] = False
        else:
            current_count += 1
            counts[target] += 1

    filtered_edge_index = edge_index_sorted[:, kept_mask]
    return filtered_edge_index, counts
