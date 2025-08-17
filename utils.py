from torch.utils.data import Dataset, DataLoader
import torch
from sksurv.util import Surv

class SupDataset(Dataset):
    def __init__(self, data, y):
        self.x = data
        self.y = y

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        sample = self.x[idx]
        labels = self.y[idx]
        return sample, labels

    def to_torch(self, device):
        self.x = torch.from_numpy(self.x.to_numpy()).float().to(device)
        self.y = torch.from_numpy(self.y.to_numpy()).float().to(device)

    def to_struct_array(self):
        self.y = Surv.from_dataframe("vital_status", "time", self.y)
        return

