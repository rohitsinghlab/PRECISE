import torch
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score


def compute_metrics(
    logits: torch.Tensor,
    labels: torch.Tensor,
    mask: torch.Tensor,
    max_neg_per_pos: int = 10,
    max_total: int = 5000,
    seed: int | None = None,
) -> dict[str, float]:
    metrics = {"auroc": 0.0, "auprc": 0.0, "auprc_gain": 0.0, "n_sampled": 0}

    if torch.isnan(logits).any():
        print("Warning: NaN detected in logits. Skipping metrics calculation.")
        return metrics

    sampled_logits, sampled_labels = sample_for_metrics(
        logits, labels, mask, max_neg_per_pos, max_total, seed
    )

    if sampled_logits is None or sampled_logits.numel() == 0:
        return metrics

    sampled_logits = torch.clamp(sampled_logits, min=-50, max=50)
    sampled_prob = torch.sigmoid(sampled_logits)
    sampled_prob_np = sampled_prob.to(torch.float32).cpu().numpy()

    if np.isnan(sampled_prob_np).any():
        print("Warning: NaN in probabilities after sigmoid. Skipping metrics.")
        return metrics

    sampled_labels_np = sampled_labels.cpu().numpy()
    if len(np.unique(sampled_labels_np)) < 2:
        return metrics

    try:
        auroc = roc_auc_score(sampled_labels_np, sampled_prob_np)
        auprc = average_precision_score(sampled_labels_np, sampled_prob_np)
        random_aupr = float(np.mean(sampled_labels_np))

        metrics["auroc"] = float(auroc)
        metrics["auprc"] = float(auprc)
        metrics["auprc_gain"] = float(auprc) / random_aupr if random_aupr > 0 else 0.0
        metrics["n_sampled"] = int(sampled_logits.numel())
    except Exception as exc:  # pragma: no cover - defensive
        print(f"Error calculating metrics: {exc}")

    return metrics


def sample_for_metrics(
    logits: torch.Tensor,
    true_labels: torch.Tensor,
    mask: torch.Tensor,
    max_neg_per_pos: int = 10,
    max_total: int = 5000,
    seed: int | None = None,
):
    valid = torch.where(mask)[0]
    if valid.numel() == 0:
        return None, None

    v_logits = logits[valid]
    v_labels = true_labels[valid]

    pos = torch.where(v_labels == 1)[0]
    neg = torch.where(v_labels == 0)[0]
    if pos.numel() == 0 or neg.numel() == 0:
        return None, None

    # keep all positives
    pos_idx = pos
    # cap negatives
    n_neg_target = min(neg.numel(), max_neg_per_pos * pos_idx.numel())
    if max_total:
        n_neg_target = min(n_neg_target, max_total - pos_idx.numel())

    g = torch.Generator(device=logits.device)
    if seed is not None:
        g.manual_seed(seed)
    neg_perm = torch.randperm(neg.numel(), generator=g, device=logits.device)[
        :n_neg_target
    ]
    neg_idx = neg[neg_perm]

    idx = torch.cat([pos_idx, neg_idx])
    return v_logits[idx], v_labels[idx]
