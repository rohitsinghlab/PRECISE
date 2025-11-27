import os
import time
import json
import asyncio
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import torch
from torch_geometric.data import Data, Dataset
from plyfile import PlyData
from tqdm import tqdm
import h5py as h5
from molfeat.trans.fp import FPVecTransformer
from torch.utils.data import Subset
import pyvista as pv
import trimesh

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def smiles_to_morgan_fingerprint(smiles: str, n_bits: int = 2048) -> torch.Tensor:
    """Convert SMILES string to Morgan fingerprint using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)

    fp_array = np.array(fp, dtype=np.float32)
    return torch.from_numpy(fp_array)


def build_single_example_batch(
    receptor_data: Data, drug_fp: torch.Tensor, device: torch.device
) -> Data:
    """Build a single example batch for DTI inference."""
    receptor_data = Data(
        pos=receptor_data.pos,
        z=receptor_data.z,
        faces=receptor_data.faces,
        drug_ecfp=drug_fp,
        batch=torch.zeros(receptor_data.pos.shape[0], dtype=torch.long),
        batch_size=1,
    )

    return receptor_data.to(device)


def data_transform(graph_data, device="cpu", use_curvature=True):
    pos_data = graph_data.pos
    vertdata = graph_data.z.clone()
    if hasattr(graph_data, "edge_index"):
        graph_data.faces = graph_data.edge_index

    pb_lower, pb_upper = -3.0, 3.0

    pb = vertdata[:, 0]
    hb = vertdata[:, 1]
    hp = vertdata[:, 2]

    pb = torch.clamp(pb, pb_lower, pb_upper)
    pb = (pb - pb_lower) / (pb_upper - pb_lower)
    pb = 2 * pb - 1
    hp = hp / 4.5

    if use_curvature:
        shape_index = calculate_shape_index_trimesh_pyvista(pos_data, graph_data.faces)
        shape_index_tensor = torch.from_numpy(shape_index).float()

        graph_data.z = torch.stack([pb, hb, hp, shape_index_tensor], dim=-1)
    else:
        graph_data.z = torch.stack([pb, hb, hp], dim=-1)

    return graph_data.to(device)


def calculate_shape_index_trimesh_pyvista(vertices, faces):
    """
    Calculate shape index using trimesh + pyvista approach.

    Args:
        vertices: torch.Tensor of shape (N, 3) - vertex positions
        faces: torch.Tensor of shape (3, M) or (M, 3) - face connectivity

    Returns:
        numpy array of shape indices for each vertex
    """

    if torch.is_tensor(vertices):
        vertices_np = vertices.detach().cpu().numpy()
    else:
        vertices_np = vertices

    if torch.is_tensor(faces):
        faces_np = faces.detach().cpu().numpy()
    else:
        faces_np = faces

    # Ensure faces are in the right format (M, 3)
    if faces_np.shape[0] == 3:
        faces_np = faces_np.T

    # Create trimesh mesh
    mesh = trimesh.Trimesh(vertices=vertices_np, faces=faces_np)

    # Convert to pyvista PolyData
    faces = np.hstack([3 * np.ones((len(mesh.faces), 1)), mesh.faces]).astype(np.int64)
    pv_mesh = pv.PolyData(mesh.vertices, faces)

    # Compute Gaussian & Mean curvature
    K = pv_mesh.curvature(curv_type="Gaussian")
    H = pv_mesh.curvature(curv_type="Mean")

    # Compute principal curvatures
    elem = np.square(H) - K
    elem[elem < 0] = 1e-8
    k1 = H + np.sqrt(elem)
    k2 = H - np.sqrt(elem)

    # Shape index
    si = (k1 + k2) / (k1 - k2)
    si = np.arctan(si) * (2.0 / np.pi)

    return si


def get_fingerprints(ligand_list):
    MORGAN_FINGERPRINT_DIMS = 2048
    transformer = FPVecTransformer(
        kind="ecfp:4", length=MORGAN_FINGERPRINT_DIMS, verbose=True
    )

    valid_features, valid_ids = transformer(ligand_list, ignore_errors=True)
    return valid_features, valid_ids


def sanitize_string(s):
    return s.replace("/", "|")


def construct_ligand_h5(ligand_list, output_h5_file, n_jobs=16):
    """
    Takes in the list of SMILES strings
    Saves the Morgan fingerprints produced in the `output_h5_file`
    """
    # Remove duplicates and sort to maintain consistent order
    ligand_list = sorted(set(ligand_list))
    valid_features, valid_ids = get_fingerprints(ligand_list)

    with h5.File(output_h5_file, "w") as hdf_file:
        for feat, idx in tqdm(zip(valid_features, valid_ids), desc="Saving h5"):
            dataset_name = ligand_list[idx]
            dataset_name = sanitize_string(dataset_name)
            if dataset_name in hdf_file:
                print(f"Dataset {dataset_name} already exists. Skipping.")
                continue
            try:
                hdf_file.create_dataset(dataset_name, data=feat)
            except Exception as e:
                print(f"Error saving {dataset_name}: {e}")

    return [ligand_list[idx] for idx in valid_ids]


def _read_smiles_from_info(path: str):
    try:
        with open(path, "r") as f:
            info = json.load(f)
        return info.get("smiles", None)
    except Exception:
        return None


async def _collect_smiles_async(
    info_paths,
    max_concurrency: int = 256,
    show_progress: bool = False,
    desc: str = "Reading SMILES",
):
    semaphore = asyncio.Semaphore(max_concurrency)

    async def read_one(p: str):
        async with semaphore:
            return await asyncio.to_thread(_read_smiles_from_info, p)

    tasks = [read_one(p) for p in info_paths]
    if not show_progress:
        return await asyncio.gather(*tasks)

    results = []
    for coro in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc=desc):
        res = await coro
        results.append(res)
    return results


def collect_smiles_asyncio(
    info_paths,
    max_concurrency: int = 256,
    show_progress: bool = False,
    desc: str = "Reading SMILES",
):
    try:
        results = asyncio.run(
            _collect_smiles_async(
                info_paths,
                max_concurrency=max_concurrency,
                show_progress=show_progress,
                desc=desc,
            )
        )
    except RuntimeError:
        # Fallback if an event loop is already running (e.g., notebooks)
        with ThreadPoolExecutor(max_workers=max_concurrency) as ex:
            iterator = ex.map(_read_smiles_from_info, info_paths)
            if show_progress:
                iterator = tqdm(iterator, total=len(info_paths), desc=desc)
            results = list(iterator)
    return [s for s in results if s is not None]


def get_data_from_ply(path, device="cpu"):
    plydata = PlyData.read(path)
    posdata = [torch.tensor(plydata["vertex"][axis]) for axis in ["x", "y", "z"]]
    posdata = torch.stack(posdata, dim=-1)

    vertdata = [
        torch.tensor(plydata["vertex"][axis])
        for axis in ["charge", "hbond", "hphob", "shape_index"]
    ]
    vertdata = torch.stack(vertdata, dim=-1)

    pb_lower, pb_upper = -3.0, 3.0
    pb = vertdata[:, 0]
    hb = vertdata[:, 1]
    hp = vertdata[:, 2]
    si = vertdata[:, 3]
    # iface = vertdata[:, 4]
    pb = torch.clamp(pb, pb_lower, pb_upper)
    pb = (pb - pb_lower) / (pb_upper - pb_lower)
    pb = 2 * pb - 1

    hp = hp / 4.5

    vertdata = torch.stack([pb, hb, hp, si], dim=-1)

    faceinfo = [
        torch.tensor(fa, dtype=torch.long) for fa in plydata["face"]["vertex_indices"]
    ]
    faceinfo = torch.stack(faceinfo, dim=0)
    data = Data(z=vertdata, faces=faceinfo, pos=posdata)  # , iface=iface)
    return data.to(device)


class PreciseOnDiskDataset(Dataset):
    def __init__(
        self,
        root_path: str,
        ecfp_n_bits: int = 2048,
        device: str = "cpu",
        transform=None,
        pre_transform=None,
        subsample_ratio=None,
        h5_path: str | None = None,
        precompute_h5: bool = True,
        smiles_io_concurrency: int = 256,
        show_smiles_progress: bool = True,
    ):
        super().__init__(root_path, transform, pre_transform)
        self.root_dir = root_path
        self.ecfp_n_bits = ecfp_n_bits
        self.device = device
        self._h5 = None  # Lazy-opened, persisted per-process handle
        self._precompute_h5 = precompute_h5
        self._smiles_io_concurrency = smiles_io_concurrency
        self._show_smiles_progress = show_smiles_progress
        # Enumerate valid protein entry directories
        entries = []
        split_map = {"train": [], "val": [], "test": [], "removed": []}

        for name in tqdm(sorted(os.listdir(self.root_dir)), desc="Processing entries"):
            p = os.path.join(self.root_dir, name)
            if not os.path.isdir(p):
                continue
            surface = os.path.join(p, "surface.ply")
            info = os.path.join(p, "info.json")
            if not (os.path.exists(surface) and os.path.exists(info)):
                continue

            try:
                with open(info, "r") as f:
                    meta = json.load(f)
            except Exception:
                meta = {}

            split = meta.get("split", "removed")
            if split not in split_map:
                split = "removed"

            idx = len(entries)
            entries.append({"dir": p, "surface": surface, "info": info, "split": split})
            split_map[split].append(idx)

        self.entries = entries
        self.split_map = split_map

        # Optional subsampling (deterministic; no shuffling)
        if subsample_ratio is not None:
            if not (0 < float(subsample_ratio) <= 1.0):
                raise ValueError("subsample_ratio must be in (0, 1].")
            n_keep = max(1, int(len(entries) * float(subsample_ratio)))
            entries = entries[:n_keep]
        self.entries = entries
        # Setup H5 path
        self.h5_path = (
            h5_path
            if h5_path is not None
            else os.path.join(self.root_dir, "ligands.h5")
        )

        # Precompute H5 only if requested and the file doesn't already exist
        if self._precompute_h5 and not os.path.exists(self.h5_path):
            # Only collect SMILES if we actually need to build the H5
            info_paths = [e["info"] for e in self.entries]
            smiles_list = collect_smiles_asyncio(
                info_paths,
                max_concurrency=self._smiles_io_concurrency,
                show_progress=self._show_smiles_progress,
                desc="Collecting SMILES",
            )
            if len(smiles_list) > 0:
                construct_ligand_h5(smiles_list, self.h5_path)

    def len(self):
        return len(self.entries)

    def _ensure_h5_open(self):
        if self._h5 is None:
            if not os.path.exists(self.h5_path):
                raise FileNotFoundError(
                    f"Fingerprint H5 file not found at {self.h5_path}. Please precompute."
                )
            # Open read-only; keep handle for reuse within this process
            self._h5 = h5.File(self.h5_path, "r")

    def __del__(self):
        try:
            if hasattr(self, "_h5") and self._h5 is not None:
                self._h5.close()
        except Exception:
            pass

    def get(self, idx):
        # Keep trying until we find a valid entry or exhaust all options
        original_idx = idx
        attempts = 0
        max_attempts = len(self.entries)

        while attempts < max_attempts:
            try:
                entry = self.entries[idx]

                # Load graph from surface.ply
                data = get_data_from_ply(
                    entry["surface"], device="cpu"
                )  # keep on CPU for DataLoader

                # Check for NaN in graph data
                if torch.isnan(data.z).any() or torch.isnan(data.pos).any():
                    print(
                        f"Warning: NaN detected in graph data for entry {idx} ({entry['dir']}), skipping..."
                    )
                    idx = (idx + 1) % len(self.entries)
                    attempts += 1
                    continue

                # Load metadata
                with open(entry["info"], "r") as f:
                    info = json.load(f)

                smiles = info.get("smiles", None)
                if smiles is None:
                    raise KeyError(f"'smiles' not found in {entry['info']}")

                # Read drug fingerprint from H5 (precomputed), using persistent handle
                self._ensure_h5_open()
                key = sanitize_string(smiles)

                if key not in self._h5:
                    raise KeyError(f"SMILES key {key} not found in {self.h5_path}.")
                vec = np.array(self._h5[key][...], dtype=np.float32)

                if vec.ndim != 1:
                    vec = vec.reshape(-1)

                # Check for NaN in fingerprint vector
                if np.isnan(vec).any():
                    print(
                        f"Warning: NaN detected in fingerprint for entry {idx} ({entry['dir']}), skipping..."
                    )
                    idx = (idx + 1) % len(self.entries)
                    attempts += 1
                    continue

                drug_vec = torch.from_numpy(vec)

                # Create node-wise labels from interface indices
                interface_indices = info.get("interface_vertex_indices", [])

                num_nodes = data.num_nodes
                y = torch.zeros((num_nodes,), dtype=torch.float32)
                if len(interface_indices) > 0:
                    idx_tensor = torch.tensor(interface_indices, dtype=torch.long)
                    # Clamp to valid range just in case
                    idx_tensor = idx_tensor[
                        (idx_tensor >= 0) & (idx_tensor < num_nodes)
                    ]
                    if idx_tensor.numel() > 0:
                        y[idx_tensor] = 1.0

                # Check for NaN in labels
                if torch.isnan(y).any():
                    print(
                        f"Warning: NaN detected in labels for entry {idx} ({entry['dir']}), skipping..."
                    )
                    idx = (idx + 1) % len(self.entries)
                    attempts += 1
                    continue

                # Attach attributes
                data.y = y
                data.drug = drug_vec  # shape [ecfp_n_bits]
                data.entry_dir = entry["dir"]

                data.pos = data.pos - data.pos.mean(dim=0)

                if self.transform is not None:
                    data = self.transform(data)

                    if torch.isnan(data.y).any() or torch.isnan(data.drug).any():
                        print(
                            f"Warning: NaN detected in transformed data for entry {idx} ({entry['dir']}), skipping..."
                        )
                        idx = (idx + 1) % len(self.entries)
                        attempts += 1
                        continue

                return data

            except Exception as e:
                print(f"Error loading entry {idx} ({entry['dir']}): {e}, skipping...")
                idx = (idx + 1) % len(self.entries)
                attempts += 1
                continue

        # If we've exhausted all attempts, raise an error
        raise RuntimeError(
            f"Could not find valid data after checking all {max_attempts} entries starting from index {original_idx}"
        )

    def get_split_subset(self, split: str):
        """
        Gets a torch.utils.data.Subset for a specific split, e.g., 'train', 'val', or 'test'.
        This method ensures that the indices for the requested split are valid.

        Args:
            split (str): The desired split name. Must be a key in self.split_map.

        Returns:
            torch.utils.data.Subset: A subset of the dataset for the specified split.
        """
        if split not in self.split_map:
            raise ValueError(
                f"Split '{split}' not found. Available splits are: {list(self.split_map.keys())}"
            )

        indices = self.split_map[split]
        if not indices:
            print(f"Warning: Split '{split}' contains no data.")
        # Verify that all indices correspond to the correct split
        for idx in indices:
            assert self.entries[idx]["split"] == split, (
                f"Index assignment error: Index {idx} is marked for split '{split}' but its data file says '{self.entries[idx]['split']}'"
            )

        return Subset(self, indices)
