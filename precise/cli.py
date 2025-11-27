import typer
from rich.console import Console


from precise.commands import training_dti, pocket_detect, screen
from precise.commands import create_surface, install_tools

app = typer.Typer(
    name="precise",
    help="PRECISE - Protein-ligand interaction prediction and virtual screening",
    add_completion=True,
    rich_markup_mode="rich",
)

console = Console()


@app.callback()
def callback():
    """
    PRECISE

    Use --help with any command to see detailed options.
    """
    pass


training_dti.register_command(app, console)
pocket_detect.register_command(app, console)
screen.register_command(app, console)
create_surface.register_command(app, console)
install_tools.register_command(app, console)


@app.command()
def version():
    console.print("[bold]Precise[/bold] version [cyan]0.1.0[/cyan]")


if __name__ == "__main__":
    app()
