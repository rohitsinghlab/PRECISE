import logging
from functools import partial
from pathlib import Path
from typing import List

import typer
from rich.console import Console
from typing_extensions import Annotated


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def slack_func(depth: int, mode: str = "wide") -> float:
    len_ = 14

    # Slack values for thin search (aggressive filtering)
    thin_slack = [4, 3, 2, 2, 2, 1, 1, 0.75, 0.75, 0.5, 0.5, 0.25, 0.25, 0.125]

    # Slack values for wide search (more exploration)
    wide_slack = [
        4,
        3,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        1,
        1,
        0.75,
        0.75,
        0.75,
        0.75,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.25,
        0.25,
    ]

    depth = min(depth, len_ - 1)
    if mode == "wide":
        return wide_slack[depth]
    else:
        return thin_slack[depth]


def ligand_embedding(smiles: List[str], embedding_choice: str = "mf", device: int = 0):
    import numpy as np
    import torch
    from precise.models import (
        load_simple_coembedding_model,
        smiles_to_morgan_fingerprints,
    )

    # Generate Morgan fingerprints
    mf_embeddings = smiles_to_morgan_fingerprints(smiles)

    if embedding_choice == "mf":
        return torch.stack(mf_embeddings, dim=0).numpy()

    elif embedding_choice == "conplex":
        # Project through ConPLex model
        model = load_simple_coembedding_model().to(device)
        model.eval()

        batched_embeddings = []
        batch_size = 512

        with torch.no_grad():
            for start in range(0, len(mf_embeddings), batch_size):
                batch = mf_embeddings[start : start + batch_size]
                batch_tensor = torch.stack(batch).to(device)
                projected = model.drug_projector(batch_tensor)
                batched_embeddings.append(projected.cpu())

        embeddings = torch.cat(batched_embeddings, dim=0)
        return embeddings.numpy().astype(np.float32)

    else:
        raise ValueError(f"Unknown embedding choice: {embedding_choice}")


def compute_surface(rec_path: Path, rec_name: str, rec_chain: str, ply_path: Path):
    from precise.dataset.precise_dataset import get_data_from_ply
    from precise.preprocess.msms_pipeline import process_protein

    process_protein(rec_path, rec_name, rec_chain, ply_path)
    assert ply_path.exists(), f"Surface file not created: {ply_path}"
    rec_item = get_data_from_ply(ply_path)
    return rec_item


def score_pocket_find_best_smiles_codes(
    surface,
    x: float,
    y: float,
    z: float,
    db_path: str,
    precise_chpt: str,
    device: int,
    smiles_out_csv: str,
    no_codes_to_consider: int,
    no_smiles_per_codes: int,
):
    import pandas as pd
    from precise.screening.codes import codes_to_smiles, score_pocket

    pocket_info = score_pocket(
        surface, x, y, z, precise_chpt, device, no_codes_to_consider
    )
    codes = pocket_info["codes"].tolist()
    probs = pocket_info["probs"].tolist()

    codes_to_smiles(
        codes, db_path, smiles_out_csv, probs, no_smiles_per_code=no_smiles_per_codes
    )

    smiles_df = pd.read_csv(smiles_out_csv)
    smiles_df["smiles"] = smiles_df["smiles"].astype(str)
    return smiles_df


def get_smiles_clusters_and_dock(
    smiles_df,
    device: int,
    d_thres: float,
    pdb_path: Path,
    out_prefix: Path,
    box_center_sdf: Path = None,
    box_center_coords: tuple = None,
    embedding_choice: str = "conplex",
    slack_func_obj: callable = None,
    max_per_depth: int = -1,
):
    import numpy as np
    from scipy.cluster.hierarchy import linkage
    from precise.screening import docking
    from precise.screening.clustering import cluster_and_score

    smiles = smiles_df["smiles"].tolist()

    def dock_scoring(idx: int) -> float:
        if box_center_sdf is not None:
            return docking.dock_and_score(
                smiles[idx], pdb_path, out_prefix, box_center_sdf
            )
        else:
            x, y, z = box_center_coords
            scores = docking.dock_many_mols(
                [smiles[idx]], pdb_path, out_prefix, (x, y, z), (22, 22, 22)
            )
            return scores[0]

    def dock_filtering(nodes_scores_list: List[tuple], depth: int) -> List[int]:
        if len(nodes_scores_list) == 0:
            return []

        slack_allowed = slack_func_obj(depth)
        nodes_scores_min = min(nodes_scores_list, key=lambda x: x[1])[1]

        nodes_scores_filtered = [
            nsc
            for nsc in nodes_scores_list
            if nsc[1] <= nodes_scores_min + slack_allowed
        ]

        if max_per_depth > 0:
            nodes_scores_filtered = sorted(nodes_scores_filtered, key=lambda x: x[1])[
                :max_per_depth
            ]

        return [nx[0] for nx in nodes_scores_filtered]

    logger.info(f"Generating {embedding_choice} embeddings for {len(smiles)} molecules")
    embeddings = ligand_embedding(smiles, embedding_choice, device)

    logger.info("Building hierarchical clustering dendrogram")
    Z = linkage(embeddings, method="ward")
    np.save(str(out_prefix / "dendrogram.npy"), Z)

    logger.info("Starting tree-based virtual screening")
    clust_sc_dict = cluster_and_score(
        Z, dock_scoring, dock_filtering, outprefix=out_prefix, smiles_list=smiles
    )

    clust_scores = [
        clust_sc_dict[idx] if idx in clust_sc_dict else None
        for idx in range(len(smiles))
    ]
    smiles_df["vina_scores"] = clust_scores

    return smiles_df


def run_virtual_screen(
    pdb_path: str,
    output_dir: str,
    center_sdf: str = None,
    center_x: float = None,
    center_y: float = None,
    center_z: float = None,
    precise_chpt: str = None,
    db_path: str = None,
    dist_thres: float = 0.7,
    ligand_embedding: str = "conplex",
    search_mode: str = "wide",
    no_codes_to_consider: int = 500,
    no_smiles_per_codes: int = 100,
    max_per_depth: int = 20,
    device: int = 0,
) -> None:
    from precise.screening.docking import sdf_center

    rec_path = Path(pdb_path).resolve()
    out_path = Path(output_dir).resolve()

    if center_sdf is not None:
        csdf_path = Path(center_sdf).resolve()
        x, y, z = sdf_center(csdf_path)
    elif all(c is not None for c in [center_x, center_y, center_z]):
        x, y, z = center_x, center_y, center_z
        csdf_path = None
    else:
        raise ValueError(
            "Must provide either center_sdf OR all of center_x, center_y, center_z"
        )

    out_path.mkdir(parents=True, exist_ok=True)
    assert rec_path.exists(), f"Receptor file not found: {rec_path}"

    # Step 1: Generate surface
    logger.info(f"Generating surface for {rec_path}")
    surface = compute_surface(rec_path, "RECP", "A", out_path / "receptor.ply")

    # Step 2: Score pocket and retrieve SMILES
    smiles_csv = out_path / "selected_smiles.csv"
    logger.info(f"Scoring pocket and retrieving SMILES → {smiles_csv}")
    smiles_df = score_pocket_find_best_smiles_codes(
        surface,
        x,
        y,
        z,
        db_path,
        precise_chpt,
        device,
        str(smiles_csv),
        no_codes_to_consider,
        no_smiles_per_codes,
    )

    # Step 3: Cluster and dock
    slack_f = partial(slack_func, mode=search_mode)
    logger.info("Clustering ligands and performing tree-based docking")
    cluster_results_df = get_smiles_clusters_and_dock(
        smiles_df,
        device,
        dist_thres,
        rec_path,
        out_path,
        box_center_sdf=csdf_path,
        box_center_coords=(x, y, z) if csdf_path is None else None,
        embedding_choice=ligand_embedding,
        slack_func_obj=slack_f,
        max_per_depth=max_per_depth,
    )

    cluster_results_csv = out_path / "docking_results_cluster_representatives.csv"
    cluster_results_df.to_csv(cluster_results_csv, index=False)
    logger.info(f"Results saved to {cluster_results_csv}")
    logger.info("Virtual screening complete!")


def register_command(app: typer.Typer, console: Console):
    @app.command()
    def virtual_screen(
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
        output_dir: Annotated[
            Path,
            typer.Option(
                help="Output directory for results",
                exists=False,
                file_okay=False,
                dir_okay=True,
                resolve_path=True,
            ),
        ],
        center_sdf: Annotated[
            Path,
            typer.Option(
                help="SDF file defining the binding site center (mutually exclusive with center-x/y/z)",
                exists=True,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ] = None,
        center_x: Annotated[
            float,
            typer.Option(
                help="X coordinate of binding site center (must provide all x/y/z or use center-sdf)"
            ),
        ] = None,
        center_y: Annotated[
            float,
            typer.Option(
                help="Y coordinate of binding site center (must provide all x/y/z or use center-sdf)"
            ),
        ] = None,
        center_z: Annotated[
            float,
            typer.Option(
                help="Z coordinate of binding site center (must provide all x/y/z or use center-sdf)"
            ),
        ] = None,
        precise_chpt: Annotated[
            Path,
            typer.Option(
                help="Path to PRECISE model checkpoint (downloads from HuggingFace Hub if not provided)",
                exists=True,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ] = None,
        db_path: Annotated[
            Path,
            typer.Option(
                help="Path to ZINC/Enamine DuckDB database",
                exists=True,
                file_okay=True,
                dir_okay=False,
                resolve_path=True,
            ),
        ] = Path(
            "/hpc/group/singhlab/user/me196/projects/moleculerep/data/dbs/zincv2.db"
        ),
        dist_thres: Annotated[
            float, typer.Option(help="Clustering distance threshold")
        ] = 0.7,
        ligand_embedding: Annotated[
            str,
            typer.Option(
                help="Ligand embedding method: 'mf' (Morgan fingerprints) or 'conplex' (ConPLex embeddings)"
            ),
        ] = "conplex",
        search_mode: Annotated[
            str,
            typer.Option(
                help="Search strategy: 'wide' for broader exploration, 'thin' for focused search"
            ),
        ] = "wide",
        no_codes_to_consider: Annotated[
            int, typer.Option(help="Number of codes to retrieve from PRECISE model")
        ] = 500,
        no_smiles_per_codes: Annotated[
            int, typer.Option(help="Number of SMILES to retrieve per code")
        ] = 100,
        max_per_depth: Annotated[
            int, typer.Option(help="Maximum nodes to explore per tree depth")
        ] = 20,
        device: Annotated[int, typer.Option(help="GPU device ID")] = 0,
    ):
        """
        Run virtual screening using hierarchical clustering and docking.

        This command performs:
        1. Protein surface generation
        2. Pocket scoring with PRECISE to identify promising molecular codes
        3. SMILES retrieval from database
        4. Hierarchical clustering of ligands
        5. Tree-based docking with Uni-Dock

        Example with SDF:
            precise virtual-screen \\
                --pdb-path receptor.pdb \\
                --output-dir results/ \\
                --center-sdf ligand_center.sdf \\
                --db-path zinc.db \\
                --search-mode wide \\
                --no-codes-to-consider 500 \\
                --no-smiles-per-codes 100

        Example with coordinates:
            precise virtual-screen \\
                --pdb-path receptor.pdb \\
                --output-dir results/ \\
                --center-x 10.5 \\
                --center-y 20.3 \\
                --center-z 15.7 \\
                --db-path zinc.db

        If --precise-chpt is not provided, the default pretrained model will be
        downloaded from HuggingFace Hub automatically.
        """
        from precise.utils.model_hub import get_checkpoint_path

        # Validate center specification
        has_sdf = center_sdf is not None
        has_coords = all(c is not None for c in [center_x, center_y, center_z])

        if not has_sdf and not has_coords:
            console.print(
                "[red]Error: Must provide either --center-sdf OR all of --center-x/y/z[/red]"
            )
            raise typer.Exit(code=1)

        if has_sdf and has_coords:
            console.print(
                "[red]Error: Cannot provide both --center-sdf AND --center-x/y/z. Choose one method.[/red]"
            )
            raise typer.Exit(code=1)

        # Get checkpoint path (download from HuggingFace if not provided)
        final_checkpoint_path = get_checkpoint_path(
            str(precise_chpt) if precise_chpt else None
        )

        console.print("[cyan]Starting virtual screening pipeline...[/cyan]")
        console.print(f"  Receptor: {pdb_path}")
        console.print(f"  Output: {output_dir}")
        if has_sdf:
            console.print(f"  Center: {center_sdf}")
        else:
            console.print(f"  Center: ({center_x}, {center_y}, {center_z})")
        console.print(f"  Checkpoint: {final_checkpoint_path}")
        console.print(f"  Search mode: {search_mode}")
        console.print(f"  Ligand embedding: {ligand_embedding}")
        console.print(f"  Codes to consider: {no_codes_to_consider}")
        console.print(f"  SMILES per code: {no_smiles_per_codes}")

        try:
            run_virtual_screen(
                pdb_path=str(pdb_path),
                output_dir=str(output_dir),
                center_sdf=str(center_sdf) if center_sdf else None,
                center_x=center_x,
                center_y=center_y,
                center_z=center_z,
                precise_chpt=final_checkpoint_path,
                db_path=str(db_path),
                dist_thres=dist_thres,
                ligand_embedding=ligand_embedding,
                search_mode=search_mode,
                no_codes_to_consider=no_codes_to_consider,
                no_smiles_per_codes=no_smiles_per_codes,
                max_per_depth=max_per_depth,
                device=device,
            )
            console.print("\n[green]Virtual screening completed successfully![/green]")
        except Exception as e:
            console.print(f"\n[red]Error during virtual screening: {e}[/red]")
            raise typer.Exit(code=1)
