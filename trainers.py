from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.svm import FastKernelSurvivalSVM
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis
import torch
from torch.utils.data import DataLoader
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import  GaussianNB
from sklearn.ensemble import RandomForestClassifier
from models import DnnClf, DeepSurv, DnnSupCon, DeepAeSurv, DalSelfNet
from pytorch_metric_learning.losses import SupConLoss, NTXentLoss, SelfSupervisedLoss
from losses import nlplloss_vectorized
from utils import zero_mask_batch

def load_algo(required_clf):
    """
    Returns Wrapper class for different survival learning algorithms
    :param required_clf: string in {"CoxPH", "RSF", "SurvMLP"}
    :return: class of wrapper
    """
    #survival trainers
    if required_clf == "CoxPH":
        return CoxNetWrapper
    elif required_clf == "Rsf":
        return RsfWrapper
    elif required_clf == "DeepSurv":
        return DeepSurvWrapper
    elif required_clf == "SurvSVM":
        return SurvSVMWrapper
    elif required_clf == "GBSA":
        return GBSAWrapper


#survival model wrappers
class CoxNetWrapper():
    def __init__(self, hps):
        self.clf = CoxnetSurvivalAnalysis(alphas=hps["alphas"], l1_ratio=hps["l1_ratio"], max_iter=hps["max_iter"]) #n_alphas=hps["n_alphas"],

    def fit(self, ds):
        self.clf.fit(ds.x,ds.y)

    def evaluate(self, ds, evaluator):
        evaluator_obj = evaluator(from_struct=True)
        haz = self.clf.predict(ds.x)
        return evaluator_obj(ds.y, haz)

class SurvSVMWrapper():
    def __init__(self, hps):
        self.clf = FastKernelSurvivalSVM(alpha=hps["alpha"],
                                         rank_ratio=hps["rank_ratio"],
                                         kernel=hps["kernel"],
                                         degree=hps["degree"],
                                         gamma=hps["gamma"],
                                         max_iter=hps["max_iter"],
                                         random_state=hps["max_iter"])

    def fit(self, ds):
        self.clf.fit(ds.x,ds.y)

    def evaluate(self, ds, evaluator):
        evaluator_obj = evaluator(from_struct=True)
        haz = self.clf.predict(ds.x)
        return evaluator_obj(ds.y, haz)


class RsfWrapper():
    def __init__(self, hps):
        self.clf = RandomSurvivalForest(n_estimators=hps["n_estimators"],
                                        min_samples_split=hps["min_samples_split"],
                                        min_samples_leaf=hps["min_samples_leaf"],
                                        max_features=hps["max_features"],
                                        max_depth=hps["max_depth"],
                                        max_samples=hps["max_samples"],
                                        random_state=hps["random_state"],
                                        n_jobs=-1,
                                        )

    def fit(self, ds):
        self.clf.fit(ds.x, ds.y)

    def evaluate(self, ds, evaluator):
        evaluator = evaluator(from_struct=True)
        haz = self.clf.predict(ds.x)
        return evaluator(ds.y, haz)

class GBSAWrapper():
    def __init__(self, hps):
        self.clf = GradientBoostingSurvivalAnalysis(
                                        learning_rate=hps["learning_rate"],
                                        n_estimators=hps["n_estimators"],
                                        min_samples_split=hps["min_samples_split"],
                                        min_samples_leaf=hps["min_samples_leaf"],
                                        max_features=hps["max_features"],
                                        max_depth=hps["max_depth"],
                                        subsample =hps["subsample"],
                                        random_state=hps["random_state"]
                                        )

    def fit(self, ds):
        self.clf.fit(ds.x, ds.y)

    def evaluate(self, ds, evaluator):
        evaluator = evaluator(from_struct=True)
        haz = self.clf.predict(ds.x)
        return evaluator(ds.y, haz)

class DeepSurvWrapper():
    def __init__(self, hps):
        super(DeepSurvWrapper, self).__init__()
        self.hps = hps
        self.clf = DeepSurv(input_dim= hps["n_features"],
                            enc_dims=hps["enc_dims"],
                            surv_dims=hps["surv_dims"],
                            proj_dims=hps["proj_dims"],
                            dropout=hps["dropout"])
        self.clf.to(hps["device"])

        # init optimizer
        self.optimizer = torch.optim.Adam(self.clf.parameters(),
                                          lr=self.hps["learning_rate"],
                                          weight_decay=self.hps["weight_decay"],
                                          eps=1e-8)

    def fit(self, ds):
        train_loader = DataLoader(ds, batch_size=self.hps["batch_size"], shuffle =True)

        # train
        for epoch in range(self.hps["epochs"]):
            for batch, labels in train_loader:
                if batch.shape[0] > 1:
                    haz = self.clf(batch)
                    loss = nlplloss_vectorized(labels[:, 0], labels[:, 1], haz, self.hps["device"])
                    loss.backward()
                    self.optimizer.step()
                    self.optimizer.zero_grad()

    def evaluate(self, ds, evaluator):
        with torch.no_grad():
            eval_obj = evaluator(from_torch=True)
            self.clf.eval()
            haz = self.clf(ds.x)
            self.clf.train()

        return eval_obj(ds.y, haz)

