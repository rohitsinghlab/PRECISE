from pathlib import Path
from typing import Optional

from torch.utils.data import Dataset, ConcatDataset
import pytorch_lightning as pl


class PrecisePickleLMDBDataset(Dataset):
    def __init__(self, root_path: str, transform: Optional[callable] = None):
        super().__init__()
        self.root_path = Path(root_path)
        self.transform = transform

        self.train_db_path = self.root_path / "train_lmdb_v2"
        self.val_db_path = self.root_path / "val_lmdb_v2"
        self.test_db_path = self.root_path / "test_lmdb"

        if not self.train_db_path.exists():
            raise FileNotFoundError(
                f"LMDB database not found at {self.train_db_path}. "
                "Please run the preprocessing script first."
            )
        if not self.val_db_path.exists():
            raise FileNotFoundError(
                f"LMDB database not found at {self.val_db_path}. "
                "Please run the preprocessing script first."
            )
        if not self.test_db_path.exists():
            raise FileNotFoundError(
                f"LMDB database not found at {self.test_db_path}. "
                "Please run the preprocessing script first."
            )

        self._load_metadata()

        self.train_env = None
        self.val_env = None
        self.test_env = None

    def _load_metadata(self):
        import lmdb

        with lmdb.open(str(self.train_db_path), readonly=True, lock=False) as env:
            with env.begin() as txn:
                self.train_len = int(txn.get(b"num_samples").decode())

        with lmdb.open(str(self.val_db_path), readonly=True, lock=False) as env:
            with env.begin() as txn:
                self.val_len = int(txn.get(b"num_samples").decode())

        with lmdb.open(str(self.test_db_path), readonly=True, lock=False) as env:
            with env.begin() as txn:
                self.test_len = int(txn.get(b"num_samples").decode())

        self.split_map = {
            "train": list(range(self.train_len)),
            "val": list(range(self.train_len, self.train_len + self.val_len)),
            "test": list(
                range(
                    self.train_len + self.val_len,
                    self.train_len + self.val_len + self.test_len,
                )
            ),
        }

        self._len = self.train_len + self.val_len + self.test_len

    def _init_db(self):
        import lmdb

        self.train_env = lmdb.open(
            str(self.train_db_path),
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        self.val_env = lmdb.open(
            str(self.val_db_path),
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        self.test_env = lmdb.open(
            str(self.test_db_path),
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )

    def __len__(self) -> int:
        return self._len

    def __getitem__(self, idx: int):
        import pickle

        if self.train_env is None:
            self._init_db()

        if idx < self.train_len:
            env = self.train_env
            local_idx = idx
        elif idx < self.train_len + self.val_len:
            env = self.val_env
            local_idx = idx - self.train_len
        else:
            env = self.test_env
            local_idx = idx - self.train_len - self.val_len

        with env.begin() as txn:
            serialized_data = txn.get(str(local_idx).encode())

        if serialized_data is None:
            raise KeyError(f"Key {local_idx} not found in database (global idx: {idx})")

        graph_data = pickle.loads(serialized_data)

        if self.transform:
            graph_data = self.transform(graph_data)

        return graph_data

    def get_split_subset(self, split: str):
        from torch.utils.data import Subset

        if split not in self.split_map:
            raise ValueError(
                f"Split '{split}' not found. Available splits: {list(self.split_map.keys())}"
            )

        indices = self.split_map[split]
        if not indices:
            print(f"Warning: Split '{split}' contains no data.")

        return Subset(self, indices)


class PreciseDataModule(pl.LightningDataModule):
    def __init__(
        self,
        data_dir: str,
        batch_size: int = 32,
        num_workers: int = 4,
        use_pyg_loader: bool = True,
    ):
        super().__init__()
        self.data_dir = data_dir
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.use_pyg_loader = use_pyg_loader

        self.full_dataset = None
        self.train_dataset = None
        self.val_dataset = None
        self.test_dataset = None

    def setup(self, stage: Optional[str] = None):
        if self.full_dataset is None:
            self.full_dataset = PrecisePickleLMDBDataset(
                root_path=self.data_dir,
                transform=None,
            )

        if stage == "fit" or stage is None:
            self.train_dataset = self.full_dataset.get_split_subset("train")
            val = self.full_dataset.get_split_subset("val")
            test = self.full_dataset.get_split_subset("test")

            self.val_dataset = ConcatDataset([val, test])

        if stage == "test" or stage is None:
            self.test_dataset = self.full_dataset.get_split_subset("test")

    def train_dataloader(self):
        from torch.utils.data import DataLoader
        from torch_geometric.loader import DataLoader as PyGDataLoader

        loader_class = PyGDataLoader if self.use_pyg_loader else DataLoader
        return loader_class(
            self.train_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            pin_memory=True,
            persistent_workers=True if self.num_workers > 0 else False,
        )

    def val_dataloader(self):
        from torch.utils.data import DataLoader
        from torch_geometric.loader import DataLoader as PyGDataLoader

        loader_class = PyGDataLoader if self.use_pyg_loader else DataLoader
        return loader_class(
            self.val_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=True,
            persistent_workers=True if self.num_workers > 0 else False,
        )

    def test_dataloader(self):
        from torch.utils.data import DataLoader
        from torch_geometric.loader import DataLoader as PyGDataLoader

        loader_class = PyGDataLoader if self.use_pyg_loader else DataLoader
        return loader_class(
            self.test_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=True,
            persistent_workers=True if self.num_workers > 0 else False,
        )
