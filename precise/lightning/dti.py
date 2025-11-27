import torch
import gc

import torch.nn as nn
import pytorch_lightning as pl
from transformers import get_cosine_schedule_with_warmup
from precise.models.precise_dti import PreciseDTI


class PreciseLightningModule(pl.LightningModule):
    def __init__(
        self,
        n_input_node_feats: int,
        drug_dim: int,
        n_atom_basis: int,
        learning_rate: float,
        concise: nn.Module,
    ):
        super().__init__()
        self.save_hyperparameters()

        self.model = PreciseDTI(
            n_input_node_feats=n_input_node_feats,
            drug_dim=drug_dim,
            n_atom_basis=n_atom_basis,
            concise=concise,
        )

        self.oom_count = 0
        self.oom_entries = set()

    def _handle_cuda_oom(self, batch_idx, batch, step_type="training"):
        torch.cuda.empty_cache()
        gc.collect()

        self.oom_count += 1

        if torch.cuda.is_available():
            print(
                f"GPU memory allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB"
            )
            print(f"GPU memory cached: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")

    def forward(self, batch):
        return self.model(batch)

    def infer(
        self,
        batch,
        quantize: bool = True,
        num_patches: int | None = None,
        center: tuple[float, float, float] = None,
        r: float = 7.0,
    ):
        return self.model.infer(
            receptor_data=batch,
            quantize=quantize,
            num_patches=num_patches,
            center=center,
            r=r,
        )

    def training_step(self, batch, batch_idx):
        try:
            output = self.forward(batch)
            loss = output["loss"]
            if loss is None:
                return loss
            auroc, auprc, auprc_gain = (
                output["auroc"],
                output["auprc"],
                output["auprc_gain"],
            )

            auroc_qd, auprc_qd, auprc_gain_qd = (
                output["auroc_qd"],
                output["auprc_qd"],
                output["auprc_gain_qd"],
            )
            auroc_d, auprc_d, auprc_gain_d = (
                output["auroc_d"],
                output["auprc_d"],
                output["auprc_gain_d"],
            )

            self.log(
                "train/loss",
                loss,
                on_step=True,
                on_epoch=True,
                prog_bar=True,
                logger=True,
                batch_size=batch.num_graphs,
            )

            self.log(
                "train/auroc",
                auroc,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auprc",
                auprc,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auprc_gain",
                auprc_gain,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )

            self.log(
                "train/auroc_qd",
                auroc_qd,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auprc_qd",
                auprc_qd,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auprc_gain_qd",
                auprc_gain_qd,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auroc_d",
                auroc_d,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auprc_d",
                auprc_d,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "train/auprc_gain_d",
                auprc_gain_d,
                on_step=True,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )

        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                self._handle_cuda_oom(batch_idx, batch, "training")
                return None
            else:
                raise e

    def validation_step(self, batch, batch_idx):
        try:
            output = self.forward(batch)
            loss = output["loss"]
            if loss is None:
                return loss

            auroc, auprc, auprc_gain = (
                output["auroc"],
                output["auprc"],
                output["auprc_gain"],
            )

            auroc_qd, auprc_qd, auprc_gain_qd = (
                output["auroc_qd"],
                output["auprc_qd"],
                output["auprc_gain_qd"],
            )
            auroc_d, auprc_d, auprc_gain_d = (
                output["auroc_d"],
                output["auprc_d"],
                output["auprc_gain_d"],
            )
            self.log(
                "val/loss",
                loss,
                on_step=False,
                on_epoch=True,
                logger=True,
                batch_size=batch.num_graphs,
            )

            self.log(
                "val/auroc",
                auroc,
                on_step=False,
                on_epoch=True,
                prog_bar=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "val/auprc",
                auprc,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "val/auprc_gain",
                auprc_gain,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )

            self.log(
                "val/auroc_qd",
                auroc_qd,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "val/auprc_qd",
                auprc_qd,
            )
            self.log(
                "val/auprc_gain_qd",
                auprc_gain_qd,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "val/auroc_d",
                auroc_d,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "val/auprc_d",
                auprc_d,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )
            self.log(
                "val/auprc_gain_d",
                auprc_gain_d,
                on_step=False,
                on_epoch=True,
                batch_size=batch.num_graphs,
            )

            return loss
        except RuntimeError as e:
            if "out of memory" in str(e).lower():
                self._handle_cuda_oom(batch_idx, batch, "validation")
                return None
            else:
                raise e
        except Exception as e:
            print(f"Error in validation step: {e}")
            return None

    def on_train_epoch_end(self):
        if self.oom_count > 0:
            self.log(
                "train/oom_count",
                self.oom_count,
                on_step=False,
                on_epoch=True,
                logger=True,
            )
            print(f"Training epoch completed with {self.oom_count} OOM events")
            print(f"Unique problematic entries: {len(self.oom_entries)}")

    def on_validation_epoch_end(self):
        if self.oom_count > 0:
            self.log(
                "val/oom_count",
                self.oom_count,
                on_step=False,
                on_epoch=True,
                logger=True,
            )
            print(f"Validation epoch completed with {self.oom_count} OOM events")

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=self.hparams.learning_rate)

        num_devices = self.trainer.num_devices
        num_training_steps = (
            len(self.trainer.datamodule.train_dataloader()) // num_devices
        ) * self.trainer.max_epochs

        scheduler = get_cosine_schedule_with_warmup(
            optimizer,
            num_warmup_steps=int(num_training_steps * 0.02),
            num_training_steps=num_training_steps,
        )

        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": scheduler,
                "interval": "step",
                "frequency": 1,
            },
        }
