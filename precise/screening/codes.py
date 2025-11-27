import pandas as pd
import numpy as np
import duckdb
from typing import List, Tuple
import torch
from precise.lightning.dti import PreciseLightningModule
from torch_geometric.data import Data


def pack3_32(tok1: int, tok2: int, tok3: int) -> np.uint16:
    alltoks = [tok1, tok2, tok3]
    assert max(alltoks) < 32 and min(alltoks) >= 0, "Tokens must be in range [0, 31]"
    return np.uint16(tok1 << 10 | tok2 << 5 | tok3)


def codes_to_smiles(
    codes: List[Tuple[int, int, int]],
    db_path: str,
    csv_outpath: str,
    probs: List[float],
    no_smiles_per_code: int = 20,
    max_number_smiles: int = 25000,
) -> pd.DataFrame:
    if len(codes) == 0:
        df = pd.DataFrame(columns=["smiles", "id", "code", "prob"])
        df.to_csv(csv_outpath, index=False)
        return df

    codes_to_search = [pack3_32(*code) for code in codes]
    codes_to_search = [int(code) for code in codes_to_search]
    codes_probs_maps = {c: p for c, p in zip(codes_to_search, probs)}

    query = f"""
    COPY (
        WITH code_table AS (
            SELECT unnest([{",".join(map(str, codes_to_search))}]) AS code
        ),
        limited_results AS (
            SELECT s.smiles, s.id, s.code
            FROM code_table c
            JOIN (
                SELECT * FROM sqlite_scan('{db_path}', 'smiles_codes')
            ) s ON c.code = s.code
            QUALIFY ROW_NUMBER() OVER (PARTITION BY s.code ORDER BY s.id) <= {no_smiles_per_code}
        )
        SELECT smiles, id, code
        FROM limited_results
        ORDER BY code, id
    ) TO '{csv_outpath}' (HEADER, DELIMITER ',')
    """

    con = duckdb.connect()
    _ = con.execute(query)
    con.close()

    df = pd.read_csv(csv_outpath)
    df["prob"] = df["code"].astype(int).apply(lambda x: codes_probs_maps[x])
    df = df.sort_values(by="prob", ascending=False).reset_index(drop=True)
    df = df.head(max_number_smiles)
    df.to_csv(csv_outpath, index=False)

    return df


def score_pocket(
    receptor: Data,
    x: float,
    y: float,
    z: float,
    precise_chpt: str,
    device: str = "cpu",
    no_codes: int = 500,
) -> dict:
    concise = torch.hub.load(
        "rohitsinghlab/CoNCISE", "pretrained_concise_v2", pretrained=True
    )
    precise_module = PreciseLightningModule.load_from_checkpoint(
        precise_chpt,
        map_location=device,
        concise=concise,
        strict=False,
    )
    precise_module.eval()

    with torch.no_grad():
        out = precise_module.model.find_best_codes(
            receptor.to(device), (x, y, z), 14.0, topk=no_codes
        )

    return out
