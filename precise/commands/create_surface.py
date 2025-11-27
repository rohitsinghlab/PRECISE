from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from typing_extensions import Annotated


def _create_surface(
    pdb_path: Path,
    output_path: Optional[Path] = None,
    chain_id: str = "A",
    color: bool = False,
    properties: str = "all",
) -> Path:
    """
    Create a molecular surface from a PDB file.

    Args:
        pdb_path: Path to the input PDB file
        output_path: Output path for the surface PLY file (defaults to PDB directory)
        chain_id: Chain identifier(s) to process (comma-separated for multimer)
        color: Create colored surface visualizations for physicochemical properties
        properties: Comma-separated properties to visualize (charge,hbond,hphob,shape_index,all)

    Returns:
        Path to the created surface file

    Raises:
        FileNotFoundError: If PDB file doesn't exist
        RuntimeError: If surface creation fails
    """
    from precise.preprocess.msms_pipeline import process_protein, process_multimer

    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    chain_ids = [c.strip() for c in chain_id.split(",")]
    is_multimer = len(chain_ids) > 1

    # Determine output path
    if output_path is None:
        if is_multimer:
            combined_chain_id = "_".join(chain_ids)
            output_path = (
                pdb_path.parent / f"{pdb_path.stem}_{combined_chain_id}_surface.ply"
            )
        else:
            output_path = pdb_path.parent / f"{pdb_path.stem}_surface.ply"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create surface
    if is_multimer:
        process_multimer(
            pdb_filename=pdb_path,
            pdb_id=str(pdb_path.stem),
            chain_ids=chain_ids,
            out_path=output_path,
        )
    else:
        process_protein(
            pdb_filename=pdb_path,
            pdb_id=str(pdb_path.stem),
            chain_id=chain_ids[0],
            out_path=output_path,
        )

    if not output_path.exists():
        raise RuntimeError(f"Surface file was not created: {output_path}")

    # Optionally create colored visualizations
    colored_files = {}
    if color:
        from precise.utils.visualization import create_colored_surfaces

        props = [p.strip() for p in properties.split(",")]
        colored_output_dir = output_path.parent / "colored_surfaces"
        colored_files = create_colored_surfaces(
            output_path, output_dir=colored_output_dir, properties=props
        )

    return output_path, colored_files


def register_command(app: typer.Typer, console: Console):
    @app.command()
    def create_surface(
        pdb_path: Annotated[
            Path,
            typer.Option(
                help="Path to the input PDB file",
                resolve_path=False,
            ),
        ],
        output_path: Annotated[
            Path,
            typer.Option(
                help="Output path for the surface PLY file",
                resolve_path=False,
            ),
        ] = None,
        chain_id: Annotated[
            str, typer.Option(help="Chain identifier to process")
        ] = "A",
        color: Annotated[
            bool,
            typer.Option(
                help="Create colored surface visualizations for physicochemical properties"
            ),
        ] = False,
        properties: Annotated[
            str,
            typer.Option(
                help="Comma-separated properties to visualize: charge,hbond,hphob,shape_index,all (default: all)"
            ),
        ] = "all",
    ):
        """
        Create a molecular surface from a PDB file.

        Example:
            precise create-surface --pdb-path protein.pdb --chain-id A --output-path surface.ply
        """
        pdb_path_resolved = pdb_path.resolve()

        console.print("[cyan]Creating molecular surface...[/cyan]")
        console.print(f"  Input PDB: {pdb_path_resolved}")
        console.print(f"  Chain ID(s): {chain_id}")
        if output_path:
            console.print(f"  Output: {output_path}")

        try:
            output_path_result, colored_files = _create_surface(
                pdb_path=pdb_path_resolved,
                output_path=output_path,
                chain_id=chain_id,
                color=color,
                properties=properties,
            )

            console.print("\n[green]Surface created successfully![/green]")
            console.print(f"  Output saved to: {output_path_result}")

            if colored_files:
                console.print(
                    f"\n[green]Created {len(colored_files)} colored surfaces:[/green]"
                )
                for prop_key, file_path in colored_files.items():
                    console.print(f"  {prop_key}: {file_path}")

        except FileNotFoundError as e:
            console.print(f"[red]Error: {e}[/red]")
            raise typer.Exit(code=1)
        except Exception as e:
            console.print(f"\n[red]Error creating surface: {e}[/red]")
            raise typer.Exit(code=1)
