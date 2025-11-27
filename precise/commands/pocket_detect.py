import json
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from typing_extensions import Annotated


def _detect_pockets(
    pdb_path: str,
    smiles: str,
    checkpoint_path: str,
    output_path: Optional[str] = None,
    threshold: float = 0.1,
    bandwidth: float = 12.0,
    quantize: bool = False,
    device: str = "cpu",
) -> dict:
    import torch
    from precise.lightning.dti import PreciseLightningModule
    from precise.pocket_detection.inference import inference

    pdb_path_obj = Path(pdb_path)
    if not pdb_path_obj.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if output_path is None:
        output_path = str(pdb_path_obj.parent)
    else:
        output_path = str(Path(output_path))
        Path(output_path).mkdir(parents=True, exist_ok=True)

    concise = torch.hub.load(
        "rohitsinghlab/CoNCISE", "pretrained_concise_v2", pretrained=True
    )
    model = PreciseLightningModule.load_from_checkpoint(
        checkpoint_path, map_location=0, concise=concise
    )

    results = inference(
        model=model,
        ligand_smi=smiles,
        receptor_pdb_path=str(pdb_path),
        device=device,
        bandwidth=bandwidth,
        threshold=threshold,
    )

    return results


def register_command(app: typer.Typer, console: Console):
    @app.command()
    def detect_pockets(
        pdb_path: Annotated[
            Path,
            typer.Option(
                help="Path to the receptor PDB file",
                exists=True,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ],
        smiles: Annotated[str, typer.Option(help="SMILES string of the ligand")],
        checkpoint_path: Annotated[
            Path,
            typer.Option(
                help="Path to the Precise model checkpoint (downloads from HuggingFace Hub if not provided)",
                exists=True,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ] = None,
        output_path: Annotated[
            Path,
            typer.Option(
                help="Output directory for intermediate files (defaults to PDB directory)",
                exists=False,
                file_okay=False,
                dir_okay=True,
                resolve_path=True,
            ),
        ] = None,
        threshold: Annotated[
            float, typer.Option(help="Score threshold for hotspot detection")
        ] = 0.35,
        bandwidth: Annotated[
            float, typer.Option(help="Bandwidth for MeanShift clustering")
        ] = 12.0,
        quantize: Annotated[
            bool, typer.Option(help="Use quantized mode for inference")
        ] = False,
        device: Annotated[
            str, typer.Option(help="Device to run inference on (cpu/cuda)")
        ] = "cuda",
        output_json: Annotated[
            Path,
            typer.Option(
                help="Optional path to save results as JSON file",
                exists=False,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ] = None,
        visualize: Annotated[
            bool,
            typer.Option(help="Create colored surface visualization of binding scores"),
        ] = False,
    ):
        """
        Detect protein pockets using PRECISE.
        Example:
            precise detect-pockets --pdb-path protein.pdb --smiles "CCO" --checkpoint-path model.ckpt

        If --checkpoint-path is not provided, the default pretrained model will be
        downloaded from HuggingFace Hub automatically.
        """
        from precise.utils.model_hub import get_checkpoint_path

        final_checkpoint_path = get_checkpoint_path(
            str(checkpoint_path) if checkpoint_path else None
        )

        results = _detect_pockets(
            pdb_path=str(pdb_path),
            smiles=smiles,
            checkpoint_path=final_checkpoint_path,
            output_path=str(output_path) if output_path else None,
            threshold=threshold,
            bandwidth=bandwidth,
            quantize=quantize,
            device=device,
        )

        results_for_json = results.copy()
        if "vertex_scores_d" in results_for_json:
            import numpy as np

            if isinstance(results_for_json["vertex_scores_d"], np.ndarray):
                results_for_json["vertex_scores_d"] = results_for_json[
                    "vertex_scores_d"
                ].tolist()

        if output_json:
            output_json.parent.mkdir(parents=True, exist_ok=True)
            with open(output_json, "w") as f:
                json.dump(results_for_json, f, indent=2)
            console.print(f"\n[green]Results saved to: {output_json}[/green]")
        else:
            console.print("\n[bold]Results:[/bold]")
            display_results = {
                k: v for k, v in results_for_json.items() if k != "vertex_scores_d"
            }
            if "vertex_scores_d" in results_for_json:
                display_results["vertex_scores_d"] = (
                    f"<array of {len(results['vertex_scores_d'])} values>"
                )
            console.print(json.dumps(display_results, indent=2))

        if visualize and "vertex_scores_d" in results:
            from precise.utils.visualization import create_binding_score_surface

            console.print("\n[cyan]Creating binding score visualization...[/cyan]")

            surface_path = Path(pdb_path).parent / "surface.ply"

            if not surface_path.exists():
                console.print(
                    f"[yellow]Warning: Surface file not found at {surface_path}. Skipping visualization.[/yellow]"
                )
            else:
                output_vis_path = (
                    Path(output_path) / f"{Path(pdb_path).stem}_binding_scores.ply"
                    if output_path
                    else Path(pdb_path).parent
                    / f"{Path(pdb_path).stem}_binding_scores.ply"
                )

                viz_path = create_binding_score_surface(
                    surface_path,
                    results["vertex_scores_d"],
                    output_path=output_vis_path,
                    cmap_type="coolwarm",
                )
                console.print(
                    f"[green]Binding score visualization saved to: {viz_path}[/green]"
                )
