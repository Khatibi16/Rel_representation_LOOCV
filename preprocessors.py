import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from pytorch_metric_learning.losses import SupConLoss
import torch

#normalizers
def load_normalizer(required_pp):
    """
    logs preprocessors that transform dat representations
    :param: string in { "log2", ...}
    :return:
    """

    if required_pp == "log2":
        return Log2Transform
    else:
        return IdentityProcessor

class Log2Transform():
    def __init__(self):
        super(Log2Transform, self).__init__()

    def fit(self, x):
        return #just for compatability

    def transform(self, x):
        return np.log2(x + 1)

    def fit_transform(self, x):
        return self.transform(x)

#data transformations
def load_data_transformer(required_dt):
    """
    Data reprensentation as input to training algorithm
    :param required_dt:  should be in {"identity", "rank", "bin"}
    :param params: dict with "n_bin"
    :return: data transformer object
    """
    if required_dt == "identity":
        return Identity
    elif required_dt == "rank":
        return Ranking
    elif required_dt == "bin":
        return Binning

class Identity():
    def __init__(self, params):
        super(Identity, self).__init__()

    def forward(self, x):
        return x

    def __call__(self, x):
        return self.forward(x)

class Ranking():
    def __init__(self, params):
        super(Ranking, self).__init__()

    def forward(self, x):
        return x.apply(lambda row: row.rank(), axis=1)

    def __call__(self, x):
        return self.forward(x)

class Binning():
    def __init__(self, params):
        super(Binning, self).__init__()
        self.n_bins = params["n_bins"]

    def forward(self, x):
        #replace zero values with Nan
        x.replace(x.min().min(), np.nan, inplace=True)

        #apply binning to each row zero values are in bin 0 then apply qcut on rest
        x = x.apply(lambda row: self.bin_row_qcut(row, self.n_bins), axis=1)
        x = x + 1
        x.fillna(0, inplace=True)

        return x

    def bin_row_qcut(self, row, bins):
        numeric_row = pd.to_numeric(row, errors="coerce")
        binned_row, bins_edges = pd.qcut(numeric_row.dropna(), q=bins, labels=False, retbins=True, duplicates="drop")
        return binned_row.reindex(row.index)

    def __call__(self, x):
        return self.forward(x)

#standardizer
def load_standardizer(required_pp):
    """
    logs preprocessors that transform dat representations
    :param: string in {"standard", "log2", ...}
    :return:
    """
    if required_pp == "standard":
        return Standardization
    else:
        return IdentityProcessor

class IdentityProcessor():
    def __init__(self):
        super(IdentityProcessor, self).__init__()

    def fit(self, x):
        return

    def transform(self, x):
        return x

    def fit_transform(self, x):
        return x

class Standardization():
    def __init__(self):
        self.scaler = StandardScaler()

    def fit(self, x):
        return self.scaler.fit(x)

    def transform(self,x):
        return pd.DataFrame(self.scaler.transform(x))

    def fit_transform(self, x):
        return pd.DataFrame(self.scaler.fit_transform(x))

#feature selections
def load_feature_selector(required_fs):
    """
    Feature Selectors (currently only highest variable gene selection)
    :param required_fs: string in {"hvg"}
    :return:
    """
    if required_fs == "hvg":
        return HvgSelector
    elif required_fs == "identity":
        return IdentitySelector

class HvgSelector():
    def __init__(self, n_features):
        self.n_features = n_features

    def fit(self, x):
        self.selected_features = x.std(axis=0).sort_values(ascending=False).index[:self.n_features]
        return

    def transform(self,x):
        return x.loc[:, self.selected_features]

    def fit_transform(self, x):
        self.fit(x)
        return self.transform(x)

class IdentitySelector():
    def __init__(self, n_features):
        return
    def fit(self, x):
        return
    def transform(self, x):
        return x
    def fit_transform(self,x):
        return x

#evaluation metrics
def load_evaluation_metric(required_metric):
    """
    Evaluation metrics currently only concordance index
    :param required_metric: string in {"c-index"}
    :return: class calculating c_index
    """
    if required_metric == "c-index":
        return HarrelC

class HarrelC():
    def __init__(self, from_torch=False, from_struct =False):
        self.from_torch = from_torch
        self.from_struct = from_struct

    def forward(self, labels, haz):
        if self.from_torch:
            return concordance_index_censored(labels.cpu().numpy()[:, 0].astype(bool), labels.cpu().numpy()[:, 1], haz.cpu().numpy()[:, 0])[0]
        elif self.from_struct:
            return concordance_index_censored(labels["vital_status"].astype(bool), labels["time"], haz)[0]

    def __call__(self, labels, haz):
        return self.forward(labels, haz)

