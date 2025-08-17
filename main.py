import argparse
import yaml
import pandas as pd
from preprocessors import load_feature_selector, load_data_transformer, load_evaluation_metric, load_normalizer, load_standardizer
from trainers import load_algo
from training_cls import PretrainFineTuneDsCV, DatasetCrossValidation
import torch
import numpy as np
import warnings


def main(exp_params, train_params):
    #reproducability
    torch.cuda.manual_seed_all(exp_params["seed"])
    np.random.seed(exp_params["seed"])

    #device
    if exp_params["to_torch"]:
        if torch.cuda.is_available():
            train_params["device"]=[torch.device("cuda")]
        else:
            train_params["device"]=[torch.device("cpu")]
    else:
        train_params["device"]=[None]

    datasets={}
    mapping = {
        "HR+HER2-": 0,
        "HR+HER2+": 1,
        "HR-HER2+": 2,
        "TNBC": 3
    }

    #building the datasets
    for project, project_files in zip(exp_params["projects_tar"], exp_params["project_files"]):
        y_tar = pd.read_csv(exp_params["path_target"]
                            + project + "_os_labels.csv",
                            header=0, index_col=0)
        y_tar.columns = ["vital_status", "time"]
        y_tar.loc[y_tar['time'] <= 0, 'time'] = 1
        x = pd.read_csv(exp_params["path_data"]
                        + project_files,
                        header=0, index_col=0)
        print(f"{project} is NaN y", y_tar.isna().sum().sum())
        print(f"{project} is NaN x", x.isna().sum().sum())
        x.fillna(0, inplace=True)

        #align samples with x
        y_tar = y_tar.loc[x.index.intersection(y_tar.index),:]
        x = x.loc[y_tar.index,:]
        print(f"{project} Shape y", y_tar.shape)
        print(f"{project} Sahpe x", x.shape)

        #in case log2 was already performed
        if x.min().min() < 0:
            x = np.power(2, x) - 1

        datasets[project] = {}
        datasets[project]["y_tar"] = y_tar.copy()
        datasets[project]["x"] = x.copy()

        if "common_features" in locals():
            common_features = common_features.intersection(x.columns)
        else:
            common_features = x.columns

    # align feature spaces
    for project in exp_params["projects_tar"]:
        datasets[project]["x"] = datasets[project]["x"].loc[:, common_features]

    #setup preprocessing
    normalization = load_normalizer(exp_params["normalization"])
    data_transformer = load_data_transformer(exp_params["data_transformer"])
    standardizer = load_standardizer(exp_params["standardizer"])
    feature_selector = load_feature_selector(exp_params["feature_selector"])

    #setup algorithm
    algorithm = load_algo(exp_params["algorithm"])

    #load evaluation metric
    evaluator = load_evaluation_metric(exp_params["evaluation_metric"])

    #train
    cv_experiment = DatasetCrossValidation(
                        normalization=normalization,
                        data_transformer=data_transformer,
                        standardizer=standardizer,
                        feature_selector=feature_selector,
                        algorithm=algorithm,
                        evaluator=evaluator,
                        n_iter=exp_params["n_iter"],
                        hyper_params=train_params,
                        device=train_params["device"][0],
                        random_seed=exp_params["seed"],
                        to_torch=exp_params["to_torch"],
                        to_struct_array=exp_params["to_struct_array"]
    )

    cv_experiment.run_experiment(datasets, exp_params["projects_tar"])

    print(cv_experiment.test_performance)

if __name__ == "__main__":
    warnings.filterwarnings("ignore")

    parser = argparse.ArgumentParser(description="path to config file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--path_config", help="path to config file")
    parser.add_argument("-po", "--project", help="tcga project")
    parser.add_argument("-s", "--seed", help="seed for hp and random state")
    config = vars(parser.parse_args())

    with open(config["path_config"], "r") as fp:
        args = yaml.safe_load(fp)

    #assumes {exp_params: {}, hyper_params: {}}
    exp_params = args["exp_params"]
    train_params = args["train_params"]
    train_params["random_state"] = [int(config["seed"])]
    exp_params["seed"] = int(config["seed"])
    exp_params["to_torch"] = bool(exp_params["to_torch"])
    exp_params["to_struct_array"] = bool(exp_params["to_struct_array"])

    #print experimental conditions
    print("Project:", exp_params["project"])
    print("normalization", exp_params["normalization"])
    print("data_transformer", exp_params["data_transformer"])
    print("Feature Selection", exp_params["feature_selector"])
    print("standardizer", exp_params["standardizer"])
    print("Algorithm", exp_params["algorithm"])
    print("n_iter", exp_params["n_iter"])
    print("Projects", exp_params["project"])
    print("Seed", exp_params["seed"])
    print("Random State", train_params["random_state"])

    main(exp_params, train_params)