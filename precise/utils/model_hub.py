from pathlib import Path
from typing import Optional
from rich.console import Console
from precise.utils.constants import HF_REPO_ID, MODEL_FILENAME

console = Console()


def get_cache_dir() -> Path:
    cache_dir = Path.home() / ".cache" / "precise" / "models"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def download_default_checkpoint(
    repo_id: str = HF_REPO_ID,
    filename: str = MODEL_FILENAME,
    force_download: bool = False,
) -> str:
    from huggingface_hub import hf_hub_download

    cache_dir = get_cache_dir()
    checkpoint_path = hf_hub_download(
        repo_id=repo_id,
        filename=filename,
        cache_dir=str(cache_dir),
        force_download=force_download,
    )
    return checkpoint_path


def get_checkpoint_path(checkpoint_path: Optional[str] = None) -> str:
    if checkpoint_path is not None:
        path = Path(checkpoint_path)
        if not path.exists():
            raise FileNotFoundError(f"Checkpoint file not found: {checkpoint_path}")
        return str(path)

    return download_default_checkpoint()
