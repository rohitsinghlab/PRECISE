from __future__ import annotations

from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from einops import rearrange
from einops import einsum
from precise.models.gotennet.models.components.layers import CosineCutoff
from precise.models.gotennet.models.representation.gotennet import GotenNetV2Wrapper
from precise.models.utils import (
    BinaryFocalLoss,
    aggregate_labels_to_patches,
    create_patches_vectorized,
    compute_metrics,
)


def codes_to_embedding_index(codes: torch.Tensor) -> torch.Tensor:
    embedding_indices = codes[:, 0] * (32 * 32) + codes[:, 1] * 32 + codes[:, 2]
    return embedding_indices


class PreciseDTI(nn.Module):
    """Precise model for DTI"""

    def __init__(
        self,
        n_input_node_feats: int = 3,
        drug_dim: int = 2048,
        concise: Optional[nn.Module] = None,
        patch_radius: float = 7.0,
        n_atom_basis: int = 384,
        dropout: float = 0.1,
        num_heads: int = 8,
        n_patches: int = 512,
        mix: bool = False,
    ):
        super().__init__()

        cutoff_fn = CosineCutoff(4)

        self.encoder = GotenNetV2Wrapper(
            n_input_node_feats=n_input_node_feats,
            n_atom_basis=n_atom_basis,
            cutoff_fn=cutoff_fn,
            n_interactions=3,
            layernorm="yes",
            steerable_norm="yes",
        )

        # Node-level feature projectors
        self.rec_node_embed_projector = nn.Sequential(
            nn.Linear(n_atom_basis, n_atom_basis, bias=False),
            nn.GELU(),
            nn.Linear(n_atom_basis, n_atom_basis, bias=False),
            nn.Dropout(dropout),
        )

        self.drug_projector = nn.Sequential(
            nn.Linear(drug_dim, n_atom_basis),
            nn.SiLU(),
            nn.Linear(n_atom_basis, n_atom_basis),
        )

        self.concise = concise if concise is not None else nn.Identity()
        self.concise_dim = 256 * 3
        self.concise_projector = nn.Sequential(
            nn.Linear(self.concise_dim, n_atom_basis),
            nn.SiLU(),
            nn.Linear(n_atom_basis, n_atom_basis),
        )
        self.focal_loss_fn = BinaryFocalLoss(alpha=0.7, gamma=2.0)

        self.n_patches = n_patches
        self.patch_radius = patch_radius
        self.mix = mix

    def forward(
        self,
        receptor_data,
        return_loss=True,
        num_patches: Optional[int] = None,
    ):
        rec_batch = receptor_data.batch
        drug = receptor_data.drug_ecfp  # [B, D]
        drug = torch.stack(drug.chunk(receptor_data.batch_size, dim=0), dim=0)

        with torch.no_grad():
            code_emb = self.concise.emb(drug).squeeze(1)

        qd_emb = self.concise_projector(rearrange(code_emb, "b n d -> b (n d)"))
        d_emb = self.drug_projector(drug)
        qd_emb = F.normalize(qd_emb, p=2, dim=-1, eps=1e-8)
        d_emb = F.normalize(d_emb, p=2, dim=-1, eps=1e-8)

        rec_embed = self.encode_surface(receptor_data)

        if torch.isnan(rec_embed).any():
            return {"loss": None}

        (
            rec_patches,
            rec_original_patches,
            rec_centers,
            rec_patch_mask,
            rec_patch_indices,
        ) = create_patches_vectorized(
            rec_embed,
            receptor_data.z,
            receptor_data.pos,
            rec_batch,
            n_patches=num_patches if num_patches is not None else self.n_patches,
            radius=self.patch_radius,
        )  # [B, n_patches, D], [B, n_patches, 3], [B, n_patches], [B, n_patches]

        rec_logits_qd = einsum(rec_patches, qd_emb, "b n d, b d -> b n")
        rec_logits_d = einsum(rec_patches, d_emb, "b n d, b d -> b n")

        payload = {
            "logits_qd": rec_logits_qd,
            "logits_d": rec_logits_d,
            "mask": rec_patch_mask,
            "patch_centers": rec_centers,
            "patch_indices": rec_patch_indices,
        }

        if return_loss:
            loss_metrics = self._compute_loss_and_metrics(
                receptor_data,
                rec_batch,
                rec_logits_qd,
                rec_logits_d,
                rec_patch_indices,
                rec_patch_mask,
            )
            payload.update(loss_metrics)
            if loss_metrics["loss"] is None:
                return payload

        return payload

    @torch.no_grad()
    def find_best_codes(
        self,
        receptor_data: Data,
        center: tuple[float, float, float] | torch.Tensor,
        radius: float,
        topk: int = 20,
        batch_index: int = 0,
        chunk_size: int = 4096,
        return_scores: bool = False,
    ) -> dict[str, torch.Tensor]:
        device = next(self.parameters()).device

        pos = receptor_data.pos.to(device)
        batch_vec = (
            receptor_data.batch.to(device)
            if (hasattr(receptor_data, "batch") and receptor_data.batch is not None)
            else torch.zeros(pos.size(0), dtype=torch.long, device=device)
        )
        center = torch.as_tensor(center, dtype=pos.dtype, device=device)

        rec_embed = self.encode_surface(receptor_data).to(device)  # [N, D]

        mask_batch = batch_vec == int(batch_index)
        if not mask_batch.any():
            raise ValueError(f"No nodes found for batch_index={batch_index}.")

        pos_b = pos[mask_batch]  # [Nb, 3]
        emb_b = rec_embed[mask_batch]  # [Nb, D]
        d2 = torch.sum((pos_b - center) ** 2, dim=-1)
        in_sphere = d2 <= (radius * radius)

        if not in_sphere.any():
            raise ValueError(
                "No receptor nodes found within the specified (center, radius)."
            )

        region_embed = emb_b[in_sphere]  # [Nb, D]

        vals = torch.arange(32, device=device, dtype=torch.long)
        all_codes = torch.cartesian_prod(vals, vals, vals)  # [32768, 3]

        scores_list = []
        probs_list = []

        for start in range(0, all_codes.size(0), chunk_size):
            end = min(start + chunk_size, all_codes.size(0))
            drug_chunk = all_codes[start:end]  # [B, 3]

            code_emb = self.concise.d_encoder.embed(drug_chunk)
            code_emb = self.concise.d_project(code_emb)
            code_mixed, _ = self.concise.d_to_d_attention(
                rearrange(code_emb, "b n k -> n b k")
            )
            code_mixed = rearrange(code_mixed, "n b k -> b n k")
            code_emb = code_emb + code_mixed

            qd_emb = self.concise_projector(rearrange(code_emb, "b n d -> b (n d)"))
            qd_emb = F.normalize(qd_emb, p=2, dim=-1, eps=1e-8)

            # logits: [B]
            logits = (
                torch.einsum("bd,Nd->bN", qd_emb, region_embed).squeeze(-1).mean(dim=-1)
            )
            probs = torch.sigmoid(logits)
            scores_list.append(logits)
            probs_list.append(probs)

        all_scores = torch.cat(scores_list, dim=0)  # [32768]
        all_probs = torch.cat(probs_list, dim=0)  # [32768]

        topk = min(int(topk), all_codes.size(0))
        top_vals, top_idx = torch.topk(all_scores, k=topk, largest=True, sorted=True)

        top_codes = all_codes[top_idx]  # [topk, 3]
        top_probs = all_probs[top_idx]  # [topk]

        out = {
            "codes": top_codes.detach().cpu(),  # LongTensor [topk, 3]
            "probs": top_probs.detach().cpu(),  # FloatTensor [topk]
            "center": center.detach().cpu(),  # FloatTensor [3]
            "radius": float(radius),
        }
        if return_scores:
            out["scores"] = top_vals.detach().cpu()  # logits
        return out

    @torch.no_grad()
    def infer(
        self,
        receptor_data: Data,
        quantize: bool = True,
        num_patches: Optional[int] = None,
        center: tuple[float, float, float] = None,
        r: float = 7.0,
    ) -> dict[str, torch.Tensor]:
        rec_batch = receptor_data.batch
        rec_embed = self.encode_surface(receptor_data)
        d_emb = None
        if quantize:
            drug = receptor_data.codes  # [B * 3]
            drug = torch.stack(drug.chunk(receptor_data.batch_size, dim=0), dim=0)
            with torch.no_grad():
                code_emb = self.concise.d_encoder.embed(drug)
                code_emb = self.concise.d_project(code_emb)
                code_mixed, _ = self.concise.d_to_d_attention(
                    rearrange(code_emb, "b n k -> b k n")
                )
                code_mixed = rearrange(code_mixed, "n b k -> b n k")
                code_emb = code_emb + code_mixed
                code_emb = code_emb.squeeze(1)

            qd_emb = self.concise_projector(rearrange(code_emb, "b n d -> b (n d)"))
            qd_emb = F.normalize(qd_emb, p=2, dim=-1, eps=1e-8)
        else:
            drug = receptor_data.drug_ecfp
            drug = torch.stack(drug.chunk(receptor_data.batch_size, dim=0), dim=0)

            with torch.no_grad():
                qd_emb = self.concise.emb(drug).squeeze(1)

            qd_emb = self.concise_projector(rearrange(qd_emb, "b n d -> b (n d)"))
            qd_emb = F.normalize(qd_emb, p=2, dim=-1, eps=1e-8)
            d_emb = self.drug_projector(drug)
            d_emb = F.normalize(d_emb, p=2, dim=-1, eps=1e-8)

        (
            rec_patches,
            rec_original_patches,
            rec_centers,
            rec_patch_mask,
            rec_patch_indices,
        ) = create_patches_vectorized(
            rec_embed,
            receptor_data.z,
            receptor_data.pos,
            rec_batch,
            n_patches=num_patches if num_patches is not None else self.n_patches,
            radius=self.patch_radius,
        )  # [B, n_patches, D], [B, n_patches, 3], [B, n_patches], [B, n_patches]

        rec_logits_qd = einsum(rec_patches, qd_emb, "b n d, b d -> b n")
        if d_emb is not None:
            rec_logits_d = einsum(rec_patches, d_emb, "b n d, b d -> b n")
        else:
            rec_logits_d = None

        payload = {
            "logits_qd": rec_logits_qd,
            "logits_d": rec_logits_d,
            "mask": rec_patch_mask,
            "patch_centers": rec_centers,
            "patch_indices": rec_patch_indices,
        }

        return payload

    def encode_surface(
        self,
        receptor_data: Data,
    ) -> torch.Tensor:
        rec_h, rec_X = self.encoder(receptor_data)
        if self.mix:
            x_norm_node = rec_X.norm(dim=1)
            rec_h = rec_h + x_norm_node
        rec_embed = self.rec_node_embed_projector(rec_h)  # [N_rec*B, D]

        return rec_embed

    def _compute_loss_and_metrics(
        self,
        receptor_data,
        rec_batch: torch.Tensor,
        rec_logits_qd: torch.Tensor,
        rec_logits_d: torch.Tensor,
        rec_patch_indices: torch.Tensor,
        rec_patch_mask: torch.Tensor,
    ) -> dict:
        rec_patch_labels = aggregate_labels_to_patches(
            receptor_data.y, rec_batch, rec_patch_indices
        )

        has_positive = bool((rec_patch_labels > 0).any().item())

        if not has_positive:
            return {
                "loss": None,
                "auroc": 0.0,
                "auprc": 0.0,
                "auprc_gain": 0.0,
                "n_sampled": 0,
            }

        flat_rec_logits_qd = rearrange(rec_logits_qd, "b n -> (b n)")
        flat_rec_logits_d = rearrange(rec_logits_d, "b n -> (b n)")
        flat_rec_true_ = rearrange(rec_patch_labels, "b n -> (b n)")

        focal_loss_qd = self.focal_loss_fn(flat_rec_logits_qd, flat_rec_true_)
        focal_loss_d = self.focal_loss_fn(flat_rec_logits_d, flat_rec_true_)

        total_loss = focal_loss_qd + focal_loss_d

        with torch.no_grad():
            rec_metrics_qd = compute_metrics(
                flat_rec_logits_qd, flat_rec_true_, rec_patch_mask
            )
            rec_metrics_d = compute_metrics(
                flat_rec_logits_d, flat_rec_true_, rec_patch_mask
            )

        return {
            "focal_loss_qd": focal_loss_qd,
            "focal_loss_d": focal_loss_d,
            "loss": total_loss,
            "auroc_qd": rec_metrics_qd["auroc"],
            "auprc_qd": rec_metrics_qd["auprc"],
            "auprc_gain_qd": rec_metrics_qd["auprc_gain"],
            "auroc_d": rec_metrics_d["auroc"],
            "auprc_d": rec_metrics_d["auprc"],
            "auprc_gain_d": rec_metrics_d["auprc_gain"],
            "auroc": (rec_metrics_qd["auroc"] + rec_metrics_d["auroc"]) / 2,
            "auprc": (rec_metrics_qd["auprc"] + rec_metrics_d["auprc"]) / 2,
            "auprc_gain": (rec_metrics_qd["auprc_gain"] + rec_metrics_d["auprc_gain"])
            / 2,
            "n_sampled": rec_metrics_qd["n_sampled"],
        }
