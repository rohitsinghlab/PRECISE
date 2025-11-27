from pathlib import Path
from typing import List, Optional, Sequence, Union
import os
import urllib.request

import numpy as np
import torch
from torch import device as torch_device


from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from precise.models.conplex_dti.architectures import SimpleCoembedding
from precise.models.conplex_dti.featurizer.protein import ProtBertFeaturizer
from precise.utils.constants import CONPLEX_URL


def get_conplex_model_path() -> Path:
    model_filename = "BindingDB_ExperimentalValidModel.pt"

    # 1. Check environment variable
    if "CONPLEX_MODEL_PATH" in os.environ:
        env_path = Path(os.environ["CONPLEX_MODEL_PATH"])
        if env_path.exists():
            return env_path

    # 2. Check default cache directory
    cache_dir = Path.home() / ".precise" / "models"
    cache_path = cache_dir / model_filename

    if cache_path.exists():
        return cache_path

    # 3. Download from CONPLEX_URL
    print(f"ConPLex model not found in cache. Downloading from {CONPLEX_URL}...")
    cache_dir.mkdir(parents=True, exist_ok=True)

    try:
        urllib.request.urlretrieve(CONPLEX_URL, cache_path)
        print(f"ConPLex model downloaded to: {cache_path}")
        return cache_path
    except Exception as e:
        raise RuntimeError(f"Failed to download ConPLex model from {CONPLEX_URL}: {e}")


def load_simple_coembedding_model(
    model_path: Optional[Union[str, Path]] = None,
    map_location: Optional[Union[str, torch_device]] = "cpu",
    strict: bool = True,
) -> SimpleCoembedding:
    if model_path is None:
        checkpoint_path = get_conplex_model_path()
    else:
        checkpoint_path = Path(model_path).expanduser().resolve()
        if not checkpoint_path.exists():
            raise FileNotFoundError(
                f"SimpleCoembedding checkpoint not found at {checkpoint_path}"
            )

    checkpoint = torch.load(checkpoint_path, map_location=map_location)
    if isinstance(checkpoint, dict):
        # Handle checkpoints that wrap the state dict inside common keys.
        for key in ("state_dict", "model_state_dict", "model"):
            if key in checkpoint:
                checkpoint = checkpoint[key]
                break

    model = SimpleCoembedding()
    model.load_state_dict(checkpoint, strict=strict)
    model.eval()
    return model


feat = None


def seq_to_embedding(sequence: str) -> torch.Tensor:
    global feat
    if feat is None:
        feat = ProtBertFeaturizer()
    return feat.transform(sequence)


def smiles_to_morgan_fingerprints(
    smiles: Sequence[str],
    radius: int = 2,
    n_bits: int = 2048,
    use_chirality: bool = True,
) -> List[torch.Tensor]:
    fingerprints: List[torch.Tensor] = []
    for idx, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Invalid SMILES at position {idx}: {smi}")

        bit_vect = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius, nBits=n_bits, useChirality=use_chirality
        )
        np_fp = np.zeros((n_bits,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(bit_vect, np_fp)
        fingerprints.append(torch.from_numpy(np_fp.astype(np.float32)))
    return fingerprints


__all__ = [
    "load_simple_coembedding_model",
    "seq_to_embedding",
    "get_conplex_model_path",
    "smiles_to_morgan_fingerprints",
]
