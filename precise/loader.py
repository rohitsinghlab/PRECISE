from torch.utils.data import Dataset, DataLoader


class PreciseData(Dataset):
    def __init__(self, pdb_folder):
        