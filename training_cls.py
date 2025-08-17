from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from utils import SupDataset
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split, ParameterSampler
import numpy as np
import pandas as pd
import copy
import seaborn as sns
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap.plot
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import torch
import matplotlib.colors as mcolors
from preprocessors import load_data_transformer

class Pipeline():
    def __init__(self,hps, normalizer, data_transformer, feature_selector, standardizer):
        self.normalizer = normalizer()
        self.data_transformer = data_transformer(hps)
        self.feature_selector = feature_selector(hps["n_features"])
        self.standardizer = standardizer()

    def fit_transform(self, x):
        x = self.normalizer.fit_transform(x)
        x = self.data_transformer(x)
        x = self.feature_selector.fit_transform(x)
        x = self.standardizer.fit_transform(x)

        if isinstance(x, pd.DataFrame):
            return x
        else:
            return pd.DataFrame(x)

    def transform(self, x):
        x = self.normalizer.transform(x)
        x = self.data_transformer(x)
        x = self.feature_selector.transform(x)
        x = self.standardizer.transform(x)

        if isinstance(x, pd.DataFrame):
            return x
        else:
            return pd.DataFrame(x)

class CrossValidation():
    # todo: update code according to PretrainOnce
    def __init__(self,
                 n_splits,
                 n_repeats,
                 normalization,
                 data_transformer,
                 standardizer,
                 feature_selector,
                 algorithm,
                 evaluator,
                 n_iter,
                 hyper_params,
                 to_torch=False,
                 to_struct_array=False,
                 device = None,
                 random_seed=123):

        super(CrossValidation, self).__init__()
        self.n_splits = n_splits
        self.n_repeats =n_repeats
        self.normalization_cls = normalization
        self.data_transformer_cls = data_transformer
        self.standardizer_cls = standardizer
        self.feature_selector_cls = feature_selector
        self.algorithm = algorithm #algorithm class
        self.evaluator = evaluator
        self.n_iter = n_iter
        self.hyper_params = hyper_params
        self.to_torch = to_torch
        self.to_struct_array =to_struct_array
        self.device = device
        self.random_seed = random_seed

        if n_repeats > 1:
            self.data_splitting = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats)
        elif n_repeats == 1:
            self.data_splitting = StratifiedKFold(n_splits=n_splits)
        else:
            print("n_repeats smaller than 1 or not numerical")
            exit() #todo: proper way to raise error

    def make_pipeline(self, hps):
        pipeline = Pipeline(hps, self.normalization_cls, self.data_transformer_cls, self.feature_selector_cls, self.standardizer_cls)

        return pipeline

    def random_search(self, ds):
        #split train/val
        train_x, val_x, train_y, val_y = train_test_split(ds.x, ds.y)

        #random hyperparameters
        best_perform = 0
        best_hps = None
        print(self.random_seed)
        for hp_setup in list(ParameterSampler(self.hyper_params, n_iter=self.n_iter, random_state=self.random_seed)):
            #dataloader
            train_ds = SupDataset(train_x, train_y)
            val_ds = SupDataset(val_x, val_y)
            pipeline = self.make_pipeline(hp_setup)
            train_ds.x = pipeline.fit_transform(train_ds.x)
            val_ds.x = pipeline.transform(val_ds.x)

            #refactor
            if self.to_torch:
                train_ds.to_torch(self.device)
                val_ds.to_torch(self.device)
            elif self.to_struct_array:
                train_ds.to_struct_array()
                val_ds.to_struct_array()

            model = self.train_model(train_ds, hp_setup)
            val_perform = model.evaluate(val_ds, self.evaluator)

            if val_perform > best_perform:
                best_perform = val_perform
                best_hps = hp_setup

        return best_hps

    def train_model(self, ds, hps):
        model = self.algorithm(hps)
        model.fit(ds)
        return model

    def refit(self, ds, hps):
        # data transformation Note: inefficient for any except "bin"
        pipeline = self.make_pipeline(hps)
        ds.x = pipeline.fit_transform(ds.x)

        if self.to_torch:
            ds.to_torch(self.device)
        elif self.to_struct_array:
            ds.to_struct_array()

        model = self.train_model(ds, hps)
        return model, pipeline

    def run_experiment(self, x, y):

        #iterate over folds
        self.test_perform = []
        for train_ind, test_ind in self.data_splitting.split(x, y.iloc[:, 0]): #todo make agnostic to survival analysis

            #split data into train/val/test datasets
            train_data = SupDataset(x.iloc[train_ind,:], y.iloc[train_ind,:])
            test_data = SupDataset(x.iloc[test_ind,:], y.iloc[test_ind,:])

            #hp choice
            hps = self.random_search(train_data)

            # refit
            model, pipeline = self.refit(train_data, hps)
            test_data.x = pipeline.transform(test_data.x)

            # refactor
            if self.to_torch:
                test_data.to_torch(self.device)
            elif self.to_struct_array:
                test_data.to_struct_array()

            #evaluate
            self.test_perform.append(model.evaluate(test_data, self.evaluator))

        self.mean_test_perform = np.mean(self.test_perform)
        self.std_test_perform = np.std(self.test_perform)
        return

class DatasetCrossValidation():
    def __init__(self,
                 normalization,
                 data_transformer,
                 standardizer,
                 feature_selector,
                 algorithm,
                 evaluator,
                 n_iter,
                 hyper_params,
                 to_torch=False,
                 to_struct_array=False,
                 device = None,
                 random_seed=123):

        self.normalization_cls = normalization
        self.data_transformer_cls = data_transformer
        self.standardizer_cls = standardizer
        self.feature_selector_cls = feature_selector
        self.algorithm = algorithm #algorithm class
        self.evaluator = evaluator
        self.n_iter = n_iter
        self.hyper_params = hyper_params
        self.to_torch = to_torch
        self.to_struct_array =to_struct_array
        self.device = device
        self.random_seed = random_seed

    def make_pipeline(self, hps):
        pipeline = Pipeline(hps, self.normalization_cls, self.data_transformer_cls, self.feature_selector_cls, self.standardizer_cls)
        return pipeline

    def random_search(self, ds):
        #split train/val
        train_x, val_x, train_y, val_y = train_test_split(ds.x, ds.y)

        #random hyperparameters
        best_perform = 0
        best_hps = None

        for hp_setup in list(ParameterSampler(self.hyper_params, n_iter=self.n_iter, random_state=self.random_seed)):
            try:
                #dataset
                train_ds = SupDataset(train_x, train_y)
                val_ds = SupDataset(val_x, val_y)

                pipeline = self.make_pipeline(hp_setup)
                train_ds.x = pipeline.fit_transform(train_ds.x)
                val_ds.x = pipeline.transform(val_ds.x)

                #refactor
                if self.to_torch:
                    train_ds.to_torch(self.device)
                    val_ds.to_torch(self.device)
                elif self.to_struct_array:
                    train_ds.to_struct_array()
                    val_ds.to_struct_array()


                model = self.train_model(train_ds, hp_setup)

                val_perform = model.evaluate(val_ds, self.evaluator)

                if val_perform > best_perform:
                    best_perform = val_perform
                    best_hps = hp_setup
                    print("new best validation", best_perform)
            except:
                pass
        return best_hps

    def train_model(self, ds, hps):
        model = self.algorithm(hps)
        model.fit(ds)
        return model

    def refit(self, ds, hps):
        pipeline = self.make_pipeline(hps)
        ds.x = pipeline.fit_transform(ds.x)

        if self.to_torch:
            ds.to_torch(self.device)
        elif self.to_struct_array:
            ds.to_struct_array()

        model = self.train_model(ds, hps)
        return model, pipeline

    def run_experiment(self, datasets, target_cohorts):
        """
        :param x: dictionary with datasets x assumes pandas df(could be vectorized)
        :param y: dictionary with dataset y assumes pandas df (could be vectorized)
        :return:
        """
        self.test_performance={}

        for test_name in target_cohorts:
            print("Test Cohort", test_name)
            # make dataset oo
            test_ds = SupDataset(data=datasets[test_name]["x"], y=datasets[test_name]["y_tar"])

            # build train dataset
            for ds_name in target_cohorts:
                if ds_name != test_name:
                    if "tar_train_ds" in locals():
                        tar_train_ds.x = pd.concat((tar_train_ds.x, datasets[ds_name]["x"].copy()), axis=0)
                        tar_train_ds.y = pd.concat((tar_train_ds.y, datasets[ds_name]["y_tar"].copy()), axis=0)

                    else:
                        tar_train_ds = SupDataset(datasets[ds_name]["x"].copy(), datasets[ds_name]["y_tar"].copy())

            # hyperparameter search
            print("SHAPE train:", tar_train_ds.x.shape)
            best_hps = self.random_search(tar_train_ds)

            # refit
            best_model, pipeline = self.refit(tar_train_ds, best_hps)

            # evaluate
            test_ds.x = pipeline.transform(test_ds.x)
            if self.to_torch:
                test_ds.to_torch(self.device)

            self.test_performance[test_name] = best_model.evaluate(test_ds, self.evaluator)
            # clean up
            del tar_train_ds, best_model, best_hps, test_ds

        return
