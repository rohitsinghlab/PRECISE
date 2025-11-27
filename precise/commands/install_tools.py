from pathlib import Path

import typer
from rich.console import Console
from typing_extensions import Annotated


def install_external_tools(
    install_dir: Path,
    force: bool = False,
    skip_msms: bool = False,
    skip_apbs: bool = False,
) -> int:
    from precise.preprocess.install_tools import main as install_tools_main

    argv = ["--install-dir", str(install_dir)]

    if force:
        argv.append("--force")
    if skip_msms:
        argv.append("--skip-msms")
    if skip_apbs:
        argv.append("--skip-apbs")

    ret_code = install_tools_main(argv)

    return ret_code


def register_command(app: typer.Typer, console: Console):
    @app.command()
    def install_tools(
        install_dir: Annotated[
            Path,
            typer.Option(
                help="Directory to install tools (default: ./tools in current dir)",
                file_okay=False,
                dir_okay=True,
                resolve_path=True,
            ),
        ] = Path.home() / ".precise" / "tools",
        force: Annotated[
            bool, typer.Option(help="Force reinstall even if tools exist")
        ] = False,
        skip_msms: Annotated[bool, typer.Option(help="Skip MSMS installation")] = False,
        skip_apbs: Annotated[bool, typer.Option(help="Skip APBS installation")] = False,
    ):
        """
        Install external dependencies (MSMS, APBS) to a local tools directory.

        Example:
            precise install-tools
            precise install-tools --install-dir /opt/precise_tools --force
        """
        console.print(f"[cyan]Starting tools installation to: {install_dir}[/cyan]")

        try:
            ret_code = install_external_tools(
                install_dir=install_dir,
                force=force,
                skip_msms=skip_msms,
                skip_apbs=skip_apbs,
            )

            if ret_code == 0:
                console.print(
                    "\n[bold green]Installation completed successfully![/bold green]"
                )
                console.print(
                    f"Tools are located in: [underline]{install_dir}/bin[/underline]"
                )

                bin_path = install_dir / "bin"
                console.print("\n[yellow]Tip:[/yellow] You can add this to your PATH:")
                console.print(f'  export PATH="{bin_path}:$PATH"')
            else:
                console.print("[bold red]Installation failed with errors.[/bold red]")
                raise typer.Exit(code=ret_code)

        except Exception as e:
            console.print(
                f"[bold red]An error occurred during installation: {e}[/bold red]"
            )
            raise typer.Exit(code=1)
