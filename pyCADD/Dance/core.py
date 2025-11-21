import json
import os
from typing import Any, Callable

import pandas as pd
from pandas import DataFrame, Series

from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV, RepeatedStratifiedKFold


def hyperparam_tuning(
    model: Any,
    param_gird: dict,
    X: DataFrame,
    y: Series,
    scoring: str = "roc_auc",
    cv: int = 5,
    n_jobs: int = -1,
    method: str = "grid",
    save_dir: str = None,
    model_name: str = None,
) -> dict:
    """Hyperparameter optimization for machine learning models.

    Args:
        model (Any): Model instance to optimize.
        param_gird (dict): Hyperparameter grid/distribution dictionary.
        X (DataFrame): Training feature data.
        y (Series): Training labels.
        scoring (str): Evaluation metric. Defaults to 'roc_auc'.
        cv (int): Number of cross-validation splits. Defaults to 5.
        n_jobs (int): Number of parallel jobs. Defaults to -1 (use all processors).
        method (str): Optimization method. Options:
            - 'grid': Grid search
            - 'random': Random search
        save_dir (str, optional): Directory to save parameter file. If None, no file is saved.
        model_name (str, optional): Name of the model for file naming.

    Returns:
        Dictionary containing optimized model parameters.

    Raises:
        ValueError: If method is not 'grid' or 'random'.
    """

    # Optimization method
    if method == "grid":
        cv_search = GridSearchCV(
            estimator=model, param_grid=param_gird, scoring=scoring, cv=cv, n_jobs=n_jobs
        )
    elif method == "random":
        cv_search = RandomizedSearchCV(
            estimator=model, param_distributions=param_gird, scoring=scoring, cv=cv, n_jobs=n_jobs
        )
    else:
        raise ValueError("Invalid method: %s" % method)

    cv_search.fit(X, y)
    best_params = cv_search.best_params_
    best_score = cv_search.best_score_
    print("Best parameters: %s" % best_params)
    print("AUC score: %s" % best_score)

    if save_dir is not None:
        params_file = "best_params_%s.json" % (
            model_name if model_name else model.__class__.__name__
        )
        with open(os.path.join(save_dir, params_file), "w") as f:
            json.dump(best_params, f)
    return best_params


def _score(
    model: Any,
    X_train: DataFrame,
    X_test: DataFrame,
    y_train: Series,
    y_test: Series,
    score_func: Callable,
) -> tuple:
    """Evaluate model performance.

    Args:
        model (Any): Model instance to evaluate.
        X_train (DataFrame): Training feature data.
        X_test (DataFrame): Testing feature data.
        y_train (Series): Training labels.
        y_test (Series): Testing labels.
        score_func (Callable): Evaluation function from sklearn.metrics.

    Returns:
        Tuple containing (train_score, test_score).
    """

    # Compatible with sklearn models and Consensus models
    # Consecutive fit operations don't record historical model weights - they overwrite directly
    # Therefore, no need to reinitialize the model
    model.fit(X_train, y_train)

    # Probability prediction
    y_train_predicted = model.predict_proba(X_train)[:, 1]
    y_test_predicted = model.predict_proba(X_test)[:, 1]

    # Sort by probability and calculate scores using the evaluation method
    train_score = score_func(y_train, y_train_predicted)
    test_score = score_func(y_test, y_test_predicted)

    return train_score, test_score


def calc_scp_score(
    X: DataFrame, y_true: Series, lower_is_better: bool = True, score_func: Callable = roc_auc_score
) -> dict:
    """Calculate Single Conformation Performance (SCP) scores.

    Args:
        X (DataFrame): Feature data with conformations as columns.
        y_true (Series): True labels.
        lower_is_better (bool): Whether the evaluation metric is a descending indicator.
        score_func (Callable): Evaluation function from sklearn.metrics. Defaults to roc_auc_score.

    Returns:
        Dictionary mapping conformation names to their SCP scores.
    """
    results_dict = {}

    # Performance of all conformations
    comformations = X.columns

    # Evaluate each conformation
    for comformation in comformations:
        if lower_is_better:
            scp = X[comformation] * -1
        else:
            scp = X[comformation]
        results_dict[comformation] = score_func(y_true, scp)

    return results_dict


def _get_cv_scp_score(
    splits: list, X: DataFrame, y: Series, score_func: Callable = roc_auc_score
) -> DataFrame:
    """Get SCP scores across cross-validation splits.

    Reports evaluation scores for the test set portion of each CV split.

    Args:
        splits (list): List of (train_indices, test_indices) tuples.
        X (DataFrame): Complete dataset used to split train and test sets based on indices.
        y (Series): Target labels.
        score_func (Callable): Evaluation function from sklearn.metrics. Defaults to roc_auc_score.

    Returns:
        DataFrame containing SCP strategy test scores across CV splits.
    """

    scp_aucs = []

    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        scp_score = calc_scp_score(X_test, y_test, True, score_func=score_func)
        scp_aucs.append(scp_score)

    cv_scp_auc_results = pd.DataFrame(scp_aucs)
    return cv_scp_auc_results


def _get_splits(
    X: DataFrame, y: Series, n_repeats: int = 30, k_folds: int = 4, random_state: int = 42
) -> list:
    """Create data indices for repeated cross-validation.

    Generates train/test indices for multiple rounds of cross-validation.

    Args:
        X (DataFrame): Feature data.
        y (Series): Target labels.
        n_repeats (int): Number of repetitions.
        k_folds (int): Number of folds per repetition.
        random_state (int): Random seed. Same seed produces identical splits.

    Returns:
        List of (train_indices, test_indices) tuples.
    """
    cv = RepeatedStratifiedKFold(n_repeats=n_repeats, n_splits=k_folds, random_state=random_state)
    return [*cv.split(X, y)]


def _repeat_cross_validation(
    model: Any, splits: list, X: DataFrame, y: Series, score_func: Callable
) -> list:
    """Perform repeated cross-validation on a model.

    Args:
        model (Any): Model instance to evaluate.
        splits (list): List of (train_indices, test_indices) tuples from RepeatedStratifiedKFold.
        X (DataFrame): Complete dataset used to split train and test sets based on indices.
        y (Series): Target labels.
        score_func (Callable): Evaluation function.

    Returns:
        List of test scores from n × k cross-validation runs.
    """
    validation_results = []

    for train_index, test_index in splits:
        # Get training and testing data for single cross-validation fold
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        # Fit model and perform evaluation scoring
        train_score, test_score = _score(model, X_train, X_test, y_train, y_test, score_func)
        validation_results.append(test_score)

    return validation_results
