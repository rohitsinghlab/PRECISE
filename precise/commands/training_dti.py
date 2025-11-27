from pathlib import Path

import typer
from rich.console import Console
from typing_extensions import Annotated


def train_model(
    data_dir: str,
    batch_size: int = 1,
    max_epochs: int = 50,
    learning_rate: float = 2e-4,
    num_workers: int = 4,
    devices: int = 1,
    seed: int = 42,
    checkpoint_path: str = None,
    log_dir: str = "lightning_logs",
    ecfp_n_bits: int = 2048,
    n_atom_basis: int = 384,
    wandb_project: str = "PRECISE",
    wandb_run_name: str = "DTI Precise Affinity",
):
    import torch
    import pytorch_lightning as pl
    from pytorch_lightning.callbacks import LearningRateMonitor
    from pytorch_lightning.loggers import WandbLogger
    from precise.lightning.dti import PreciseLightningModule
    from precise.dataset.dti_dataset import PreciseDataModule

    pl.seed_everything(seed, workers=True)

    torch.set_float32_matmul_precision("medium")

    concise = torch.hub.load(
        "rohitsinghlab/CoNCISE", "pretrained_concise_v2", pretrained=True
    )

    datamodule = PreciseDataModule(
        data_dir=data_dir,
        batch_size=batch_size,
        num_workers=num_workers,
        use_pyg_loader=True,
    )

    model = PreciseLightningModule(
        n_input_node_feats=4,
        drug_dim=ecfp_n_bits,
        n_atom_basis=n_atom_basis,
        learning_rate=learning_rate,
        concise=concise,
    )

    trainer = pl.Trainer(
        max_epochs=max_epochs,
        accelerator="auto",
        devices=devices,
        default_root_dir=log_dir,
        log_every_n_steps=10,
        logger=WandbLogger(name=wandb_run_name, project=wandb_project),
        val_check_interval=0.5,
        gradient_clip_val=1.0,
        callbacks=[
            pl.callbacks.ModelCheckpoint(
                monitor="val/loss",
                mode="min",
                save_top_k=20,
                filename="precisedti-{epoch:02d}-{val_loss:.2f}",
                save_on_train_epoch_end=False,
            ),
            LearningRateMonitor(logging_interval="step"),
        ],
    )

    if checkpoint_path:
        print(f"Resuming training from checkpoint: {checkpoint_path}")
        trainer.fit(model, datamodule, ckpt_path=checkpoint_path)
    else:
        print("Starting training from scratch")
        trainer.fit(model, datamodule)


def register_command(app: typer.Typer, console: Console):
    @app.command()
    def train_dti(
        data_dir: Annotated[
            Path,
            typer.Option(
                help="Root directory of the preprocessed dataset containing LMDB files",
                exists=True,
                file_okay=False,
                dir_okay=True,
                resolve_path=True,
            ),
        ],
        batch_size: Annotated[
            int, typer.Option(help="Batch size (forced to 1 due to model signature)")
        ] = 1,
        max_epochs: Annotated[
            int, typer.Option(help="Maximum number of training epochs")
        ] = 50,
        learning_rate: Annotated[
            float, typer.Option(help="Initial learning rate")
        ] = 2e-4,
        num_workers: Annotated[
            int, typer.Option(help="Number of data loading workers")
        ] = 4,
        devices: Annotated[
            int, typer.Option(help="Number of GPUs to use (if available)")
        ] = 1,
        seed: Annotated[int, typer.Option(help="Random seed for reproducibility")] = 42,
        checkpoint_path: Annotated[
            Path,
            typer.Option(
                help="Path to checkpoint file to resume training from",
                exists=True,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ] = None,
        log_dir: Annotated[
            Path, typer.Option(help="Directory for logs and checkpoints")
        ] = Path("lightning_logs"),
        ecfp_n_bits: Annotated[
            int, typer.Option(help="Number of bits for ECFP fingerprint")
        ] = 2048,
        n_atom_basis: Annotated[
            int, typer.Option(help="Dimension of atom basis in GotenNet")
        ] = 384,
        wandb_project: Annotated[
            str, typer.Option(help="Weights & Biases project name")
        ] = "PRECISE",
        wandb_run_name: Annotated[
            str, typer.Option(help="Weights & Biases run name")
        ] = "DTI Precise Affinity",
    ):
        """
        Train the PRECISE DTI.
        Example:
            precise train-dti --data-dir /path/to/data --max-epochs 100 --devices 2
        """
        if batch_size != 1:
            console.print(
                "[yellow]Warning: Batch size is forced to 1 due to the model's forward signature.[/yellow]"
            )
            batch_size = 1

        train_model(
            data_dir=str(data_dir),
            batch_size=batch_size,
            max_epochs=max_epochs,
            learning_rate=learning_rate,
            num_workers=num_workers,
            devices=devices,
            seed=seed,
            checkpoint_path=str(checkpoint_path) if checkpoint_path else None,
            log_dir=str(log_dir),
            ecfp_n_bits=ecfp_n_bits,
            n_atom_basis=n_atom_basis,
            wandb_project=wandb_project,
            wandb_run_name=wandb_run_name,
        )
