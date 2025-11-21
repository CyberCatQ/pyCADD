import json
import logging
import os
import pickle
from typing import Any, Callable, Literal

import numpy as np
import pandas as pd
from pandas import DataFrame, Series
from pyCADD.Dance import core
from pyCADD.Dance.algorithm import default_params, consensus
from pyCADD.Dance.metrics import nef_score
from pyCADD.utils.tool import makedirs_from_list
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split

logger = logging.getLogger("pyCADD.Dance.base")

DATA_MINING_DIR = os.path.join(os.getcwd(), "dm_result")
DATA_PICKLE_DIR = os.path.join(DATA_MINING_DIR, "data_pickle")
DATA_CSV_DIR = os.path.join(DATA_MINING_DIR, "data_csv")
FIGURES_DIR = os.path.join(DATA_MINING_DIR, "figures")
MODELS_DIR = os.path.join(DATA_MINING_DIR, "models_pickle")
PARAMS_DIR = os.path.join(DATA_MINING_DIR, "models_params")


class Dancer:
    """Data analyzer for CADD (Computer-Aided Drug Design).

    A comprehensive data preprocessing and analysis tool for molecular datasets
    in computer-aided drug design workflows. Handles dataset merging, preprocessing,
    and preparation for machine learning models.
    """

    def __init__(self) -> None:
        """Initialize the Dancer instance.

        Creates necessary directories and initializes data storage attributes.
        """
        # Create necessary directories
        makedirs_from_list(self._required_dirs)
        self.current_datasets = []
        self.merged_data = None

    @property
    def _required_dirs(self):
        return [DATA_MINING_DIR, DATA_PICKLE_DIR, DATA_CSV_DIR, FIGURES_DIR, MODELS_DIR, PARAMS_DIR]

    def _add_dataset(
        self, csv_path: str, activity: Literal["positive", "negative"] = None, *args, **kwargs
    ) -> None:
        """Add a sample dataset.

        Args:
            csv_path (str): Path to the sample dataset CSV file.
            activity (str, optional): Activity type of the sample dataset. None for undefined.
            *args: Additional arguments passed to pd.read_csv().
            **kwargs: Additional keyword arguments passed to pd.read_csv().
        """
        dataset = pd.read_csv(csv_path, *args, **kwargs)
        dataset = dataset.sort_index(axis=1)
        # Drop rows with all features missing
        dataset.dropna(axis=0, how="all", inplace=True)
        self.current_datasets.append({"dataset": dataset, "activity": activity})

    def add_pos_dataset(self, csv_path: str) -> None:
        """Add a positive sample dataset.

        Args:
            csv_path (str): Path to the positive sample dataset CSV file.
        """
        self._add_dataset(csv_path, "positive", index_col=0)

    def add_neg_dataset(self, csv_path: str) -> None:
        """Add a negative sample dataset.

        Args:
            csv_path (str): Path to the negative sample dataset CSV file.
        """
        self._add_dataset(csv_path, "negative", index_col=0)

    def add_dataset(self, csv_path: str, *args, **kwargs) -> None:
        """Add a dataset.

        Args:
            csv_path (str): Path to the dataset CSV file.
            *args: Additional arguments passed to _add_dataset.
            **kwargs: Additional keyword arguments passed to _add_dataset.
        """
        self._add_dataset(csv_path, index_col=0, *args, **kwargs)

    def _add_label_col(self) -> None:
        """Add activity-corresponding label column.

        Converts activity labels to numerical values:
        - 'positive' -> 1
        - 'negative' -> 0
        - None -> 'Undefined'
        """
        for dataset in self.current_datasets:
            if dataset["activity"] is None:
                dataset["dataset"]["activity"] = "Undefined"
            elif dataset["activity"] == "positive":
                dataset["dataset"]["activity"] = 1
            elif dataset["activity"] == "negative":
                dataset["dataset"]["activity"] = 0

    def _concat_datasets(self) -> DataFrame:
        """Merge datasets.

        Concatenates all added datasets and stores the result in the
        merged_data attribute.

        Returns:
            The merged dataset.

        Raises:
            ValueError: If no datasets have been added.
        """
        if len(self.current_datasets) == 0:
            raise ValueError("No datasets added.")
        self._add_label_col()
        self.merged_data = pd.concat([dataset["dataset"] for dataset in self.current_datasets])
        return self.merged_data

    def get_merged_data(self) -> DataFrame:
        """Return the merged dataset.

        Returns:
            The merged dataset DataFrame.
        """
        return self.merged_data

    def _fill_nan(
        self, fill_na_value: Any = 0, dataset: DataFrame = None, inplace: bool = False
    ) -> None:
        """Fill missing values in the dataset.

        Args:
            fill_na_value (Any): Value to fill missing values with. Defaults to 0.
            dataset (DataFrame, optional): Dataset to fill missing values in. Defaults to self.merged_data.
            inplace (bool): Whether to modify the dataset in place.

        Returns:
            None if inplace is True, otherwise the modified dataset.
        """
        if dataset is None:
            dataset = self.merged_data
        if inplace:
            dataset.fillna(fill_na_value, inplace=True)
            return None
        else:
            return dataset.fillna(fill_na_value)

    def prepare_data(self, fill_nan: bool = True, *args, **kwargs) -> None:
        """Prepare the dataset for analysis.

        Args:
            fill_nan (bool): Whether to fill missing values. Defaults to True.
            *args: Additional arguments passed to self._fill_nan().
            **kwargs: Additional keyword arguments passed to self._fill_nan().
                The 'value' parameter can specify the fill value (default 0).
        """
        self._concat_datasets()
        if fill_nan:
            self._fill_nan(inplace=True, *args, **kwargs)

    def save_pickle(self, file_name: str, dataset: DataFrame = None) -> None:
        """Save dataset as pickle file.

        Args:
            file_name (str): Name of the pickle file to save.
            dataset (DataFrame, optional): Dataset to save. Defaults to self.merged_data.
        """
        dataset = dataset if dataset is not None else self.merged_data
        dataset.to_pickle(os.path.join(DATA_PICKLE_DIR, file_name))

    def save_csv(self, file_name: str, dataset: DataFrame = None) -> None:
        """Save dataset as CSV file.

        Args:
            file_name (str): Name of the CSV file to save.
            dataset (DataFrame, optional): Dataset to save. Defaults to self.merged_data.
        """
        dataset = dataset if dataset is not None else self.merged_data
        dataset.to_csv(os.path.join(DATA_CSV_DIR, file_name))

    def save(self, file_name: str, dataset: DataFrame = None) -> None:
        """Save dataset in the appropriate format based on file extension.

        Args:
            file_name (str): Name of the file to save. Must end with .csv, .pickle, or .pkl.
            dataset (DataFrame, optional): Dataset to save. Defaults to self.merged_data.

        Raises:
            ValueError: If file extension is not supported.
        """
        dataset = dataset if dataset is not None else self.merged_data
        if file_name.endswith(".csv"):
            self.save_csv(file_name, dataset)
        elif file_name.endswith(".pickle") or file_name.endswith(".pkl"):
            self.save_pickle(file_name, dataset)
        else:
            raise ValueError("file_name must end with .csv or .pickle")


class Matrix:
    """Result matrix for molecular data analysis.

    A data structure for handling molecular datasets with train/test splitting
    and data preprocessing capabilities for machine learning workflows.
    """

    def __init__(
        self, dataframe: DataFrame, test_size: float = 0.25, random_seed: int = 42
    ) -> None:
        """Initialize Matrix with molecular data.

        Args:
            dataframe (DataFrame): Input DataFrame containing molecular data.
            test_size (float): Proportion of data to use for testing. Defaults to 0.25.
            random_seed (int): Random seed for reproducible splits. Defaults to 42.
        """
        self.raw_data = dataframe
        self.data = self.raw_data.copy()
        self.test_size = test_size
        self.random_seed = random_seed

        self.train_data = None
        self.test_data = None
        self.X = None
        self.y = None

    @classmethod
    def from_pickle(cls, path: str, *args, **kwargs) -> "Matrix":
        """Create Matrix instance from pickle file.

        Args:
            path (str): Path to the pickle file.
            *args: Additional arguments passed to Matrix constructor.
            **kwargs: Additional keyword arguments passed to Matrix constructor.

        Returns:
            Matrix instance loaded from pickle file.
        """
        data = pd.read_pickle(path)
        return cls(data, *args, **kwargs)

    @classmethod
    def from_csv(cls, path: str, *args, **kwargs) -> "Matrix":
        """Create Matrix instance from CSV file.

        Args:
            path (str): Path to the CSV file.
            *args: Additional arguments passed to Matrix constructor.
            **kwargs: Additional keyword arguments passed to Matrix constructor.

        Returns:
            Matrix instance loaded from CSV file.
        """
        data = pd.read_csv(path)
        return cls(data, *args, **kwargs)

    @classmethod
    def from_splited_data(cls, train_data: DataFrame, test_data: DataFrame) -> "Matrix":
        """Create Matrix instance from pre-split data.

        Args:
            train_data (DataFrame): Training dataset.
            test_data (DataFrame): Testing dataset.

        Returns:
            Matrix instance with pre-split train and test data.
        """
        data = pd.concat([train_data, test_data], axis=0)
        cls_instance = cls(data)
        cls_instance.train_data = train_data
        cls_instance.test_data = test_data
        return cls_instance

    def split_train_test_data(
        self, test_size: float = None, random_seed: int = None, label_col: str = "activity"
    ) -> tuple:
        """Split data into training and testing sets.

        Args:
            test_size (float, optional): Proportion of the dataset to include in the test split.
                If None, uses the instance's test_size.
            random_seed (int, optional): Random seed for reproducible splits.
                If None, uses the instance's random_seed.
            label_col (str): Name of the label column. Defaults to 'activity'.

        Returns:
            Tuple containing (train_data, test_data).
        """
        if test_size is None:
            test_size = self.test_size
        if random_seed is None:
            random_seed = self.random_seed

        X, y = self.data.drop([label_col], axis=1), self.data[label_col]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_seed, stratify=y
        )
        self.X = X
        self.y = y
        self.train_data = pd.concat([X_train, y_train], axis=1)
        self.test_data = pd.concat([X_test, y_test], axis=1)

        return self.train_data, self.test_data

    def get_train_data(self, label_col: str = "activity") -> DataFrame:
        """Get training data.

        Args:
            label_col (str): Name of the label column. Defaults to 'activity'.

        Returns:
            Training dataset DataFrame.
        """
        if self.train_data is None:
            self.split_train_test_data(label_col=label_col)
        return self.train_data

    def get_test_data(self, label_col: str = "activity") -> DataFrame:
        """Get testing data.

        Args:
            label_col (str): Name of the label column. Defaults to 'activity'.

        Returns:
            Testing dataset DataFrame.
        """
        if self.test_data is None:
            self.split_train_test_data(label_col=label_col)
        return self.test_data


class Evaluator:
    """Model performance evaluator.

    A comprehensive evaluation toolkit for machine learning models in CADD workflows.
    Suitable for non-neural network models with cross-validation and test set evaluation.
    """

    def __init__(self, matrix: Matrix, label_col: str = "activity") -> None:
        """Initialize the Evaluator with a data matrix.

        Args:
            matrix (Matrix): Result matrix containing the molecular data.
            label_col (str): Name of the label column. Defaults to 'activity'.
        """
        self.matrix = matrix
        self.train_data = matrix.get_train_data()
        self.test_data = matrix.get_test_data()
        self.X = matrix.X
        self.y = matrix.y
        self.X_train = self.train_data.drop([label_col], axis=1)
        self.y_train = self.train_data[label_col]
        self.X_test = self.test_data.drop([label_col], axis=1)
        self.y_test = self.test_data[label_col]

        self.cv_results = {}
        self.clfs = {}
        self.scp_auc_df = None
        self.score_func_name = None

    @property
    def gbt_default_params(self) -> dict:
        """Get default GBT parameter space.

        Returns:
            Dictionary containing default Gradient Boosting Tree parameters.
        """
        return default_params.GBT_DEFAULT_PARAMS

    @property
    def lr_default_params(self) -> dict:
        """Get default Logistic Regression parameter space.

        Returns:
            Dictionary containing default Logistic Regression parameters.
        """
        return default_params.LR_DEFAULT_PARAMS

    @property
    def rf_default_params(self) -> dict:
        """Get default Random Forest parameter space.

        Returns:
            Dictionary containing default Random Forest parameters.
        """
        return default_params.RF_DEFAULT_PARAMS

    @staticmethod
    def get_weights(y: Series) -> list:
        """Calculate weights for each class label in imbalanced datasets.

        Args:
            y (Series): Series containing class labels.

        Returns:
            List of weights for each sample, inversely proportional to class frequency.
        """
        # Use weighted sampling for imbalanced datasets
        label_to_count = y.value_counts()
        _weights = []
        for label in y:
            _weights.append(1 / label_to_count[label])
        return _weights

    def get_lr_default_params(self) -> dict:
        """Get default Logistic Regression parameter space.

        Returns:
            Dictionary containing default LR parameters.
        """
        return self.lr_default_params

    def get_rf_default_params(self) -> dict:
        """Get default Random Forest parameter space.

        Returns:
            Dictionary containing default RF parameters.
        """
        return self.rf_default_params

    def get_gbt_default_params(self) -> dict:
        """Get default Gradient Boosting Tree parameter space.

        Returns:
            Dictionary containing default GBT parameters.
        """
        return self.gbt_default_params

    def load_params(self, path: str) -> dict:
        """Load parameters from file.

        Args:
            path (str): Path to the parameter file.

        Returns:
            Dictionary containing loaded parameters.
        """
        with open(path, "r") as f:
            params = json.load(f)
        return params

    def save_params(self, file_name: str, params: dict) -> None:
        """Save parameters to file.

        Args:
            file_name (str): Name of the parameter file.
            params (dict): Dictionary containing parameters to save.
        """
        with open(os.path.join(PARAMS_DIR, file_name), "w") as f:
            json.dump(params, f)

    def search_params(
        self, clf: Any, params_grid: dict, method: str = "grid", *args, **kwargs
    ) -> dict:
        """Perform hyperparameter search for the model.

        Args:
            clf (Any): Classifier instance.
            params_grid (dict): Parameter space dictionary.
            method (str): Optimization method. Options:
                - 'grid': Grid search
                - 'random': Random search
            *args: Additional arguments passed to hyperparam_tuning function.
            **kwargs: Additional keyword arguments passed to hyperparam_tuning function.

        Returns:
            Dictionary containing the best parameters found.
        """
        params = core.hyperparam_tuning(
            clf,
            params_grid,
            self.X_train,
            self.y_train,
            save_dir=PARAMS_DIR,
            method=method,
            *args,
            **kwargs,
        )
        return params

    def add_clf(self, clf: Any, clf_name: str = None) -> None:
        """Add classifier instance to the evaluation dictionary.

        Classifiers with the same name will be overwritten.

        Args:
            clf (Any): Classifier instance to add.
            clf_name (str, optional): Name for the classifier. If None, uses clf.__class__.__name__.
        """
        if clf_name is None:
            clf_name = clf.__class__.__name__

        self.clfs[clf_name] = clf

    def del_clf(self, clf_name: str) -> None:
        """Delete classifier instance.

        Args:
            clf_name (str): Name of the classifier to delete.
        """
        del self.clfs[clf_name]

    def get_clf(self, clf_name: str) -> Any:
        """Get a single classifier instance by name.

        Args:
            clf_name (str): Name of the classifier to retrieve.

        Returns:
            The classifier instance.
        """
        return self.clfs[clf_name]

    def get_clfs_dict(self) -> dict:
        """Get the dictionary of all classifier instances.

        Returns:
            Dictionary mapping classifier names to instances.
        """
        return self.clfs

    def print_classifier_info(self) -> None:
        """Print parameter information for all added classifiers.

        Displays a formatted table showing classifier names, parameters, and their values.
        """
        formatter = "{0:<25}{1:<30}{2:<40}"
        print("=" * 100)
        print(formatter.format("Classifier", "Parameters", "Values"))
        print("-" * 100)
        for clf_name, clf in self.clfs.items():
            parameters = [(k, v) for k, v in clf.get_params().items()]
            print(formatter.format(clf_name, str(parameters[0][0]), str(parameters[0][1])))
            for k, v in parameters[1:]:
                print(formatter.format("", str(k), str(v)))
            print("-" * 100)
        print("=" * 100)

    def repeat_cv(
        self,
        n_repeats: int = 30,
        k_folds: int = 4,
        random_seed: int = 42,
        score_func: Callable = roc_auc_score,
        use_train_set_only: bool = False,
    ) -> dict:
        """Perform repeated cross-validation on all added classifiers.

        Args:
            n_repeats (int): Number of repetitions for cross-validation.
            k_folds (int): Number of folds for cross-validation.
            random_seed (int): Random seed for reproducible results.
            score_func (Callable): Evaluation function. Defaults to roc_auc_score.
            use_train_set_only (bool): Whether to use only training set for cross-validation.
                If False, uses the complete dataset for cross-validation.

        Returns:
            Dictionary containing evaluation results:
                - SCP results for single conformation performance
                - clf_cv_results for cross-validation results of classifiers

        Raises:
            ValueError: If no classifiers have been added.
        """
        if len(self.clfs) == 0:
            raise ValueError("No classifier is added.")

        self.score_func_name = score_func.__name__

        print("\nRepeated Cross-Validation")
        print(f"Reapeat: {n_repeats} times")
        print(f"K-Folds: {k_folds} folds")
        print(f"Score Function: {score_func.__name__}")
        print("=" * 50)

        if use_train_set_only:
            cv_dataset = self.X_train
            cv_datalabel = self.y_train
        else:
            cv_dataset = self.X
            cv_datalabel = self.y

        splits = core._get_splits(cv_dataset, cv_datalabel, n_repeats, k_folds, random_seed)
        self.scp_score_df = core._get_cv_scp_score(
            splits, cv_dataset, cv_datalabel, score_func=score_func
        )
        best_scp_score = self.scp_score_df.max().max()
        mean_scp_score = self.scp_score_df.mean().mean()
        worst_scp_score = self.scp_score_df.min().min()
        std_scp_score = self.scp_score_df.values.flatten().std()

        self.cv_results["best_scp"] = best_scp_score
        self.cv_results["mean_scp"] = mean_scp_score
        self.cv_results["worst_scp"] = worst_scp_score
        self.cv_results["std_scp"] = std_scp_score

        print("SCP Best Score: %.4f" % best_scp_score)
        print("SCP Mean Score: %.4f" % mean_scp_score)
        print("SCP Worst Score: %.4f" % worst_scp_score)
        print("SCP Score std: %.4f" % std_scp_score)

        self.cv_results["clf_cv_results"] = {}
        for clf_name, clf in self.clfs.items():
            print("-" * 50)
            print(f"Validating {clf_name} ...")
            repeat_cv_result = core._repeat_cross_validation(
                clf, splits, cv_dataset, cv_datalabel, score_func
            )
            self.cv_results["clf_cv_results"][clf_name] = repeat_cv_result
            print("%s mean Score: %.4f" % (clf_name, np.mean(repeat_cv_result)))
            print("-" * 50)

        print("=" * 50)
        print(f"{n_repeats} x {k_folds} cross-validation finished.")
        print("=" * 50)

        with open(f"cv_results_{self.score_func_name}.json", "w") as f:
            json.dump(self.cv_results, f)

        return self.cv_results

    def print_cv_results(self) -> None:
        """Print cross-validation results.

        Displays formatted results including SCP scores and classifier performance.

        Raises:
            ValueError: If no cross-validation has been performed.
        """
        if len(self.cv_results) == 0:
            raise ValueError("No cross-validation is performed.")

        formatter = "{0:<25}{1:<30}"

        print(f"\nScore: {self.score_func_name}")
        print("=" * 50)
        print(formatter.format("Best SCP", "%.4f" % self.cv_results["best_scp"]))
        print(formatter.format("Mean SCP", "%.4f" % self.cv_results["mean_scp"]))
        print(formatter.format("Worst SCP", "%.4f" % self.cv_results["worst_scp"]))
        print(formatter.format("Std SCP", "%.4f" % self.cv_results["std_scp"]))
        print("-" * 50)
        print(formatter.format("Classifier", "Mean Score"))
        print("-" * 50)
        for clf_name, clf_cv_results in self.cv_results["clf_cv_results"].items():
            print(
                formatter.format(clf_name, "%.4f" % np.mean(clf_cv_results))
            )  # Print average AUC for each classifier
        print("=" * 50)

    @staticmethod
    def _eval_single_clf(y_test: Series, y_pred: Series, y_proba: Series) -> dict:
        """Evaluate the performance of a single classifier.

        Args:
            y_test (Series): True labels from the test set.
            y_pred (Series): Predicted labels.
            y_proba (Series): Predicted probabilities.

        Returns:
            Dictionary containing classifier evaluation metrics including:
            - auc: Area Under the ROC Curve
            - nef: Normalized Enrichment Factor
            - f1: F1 score
            - accuracy: Classification accuracy
            - precision: Precision score
            - recall: Recall score
            - fpr: False Positive Rate values
            - tpr: True Positive Rate values
        """
        auc = roc_auc_score(y_test, y_proba)
        nef = nef_score(y_test, y_proba)
        f1 = f1_score(y_test, y_pred)
        acc = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        fpr, tpr, thresholds = roc_curve(y_test, y_proba)

        return {
            "auc": auc,
            "nef": nef,
            "f1": f1,
            "accuracy": acc,
            "precision": precision,
            "recall": recall,
            "fpr": list(fpr),
            "tpr": list(tpr),
        }

    def testset_eval(self) -> dict:
        """Evaluate performance using test set data.

        Returns:
            Dictionary containing classifier evaluation results on the test set.

        Raises:
            ValueError: If no classifiers have been added.
        """
        if len(self.clfs) == 0:
            raise ValueError("No classifier is added.")

        testset_eval_results = {}
        print("\nTestset Evaluation")
        print("=" * 50)

        for clf_name, clf in self.clfs.items():
            clf.fit(self.X_train, self.y_train)  # Train model
            testset_eval_results.update(self._testset_eval(clf_name, clf))  # Evaluate model
            self._save_model(clf_name, clf)  # Save model

        print("=" * 50)
        print("Testset Evaluation Finished.")
        print("=" * 50)

        with open("testset_eval_results.json", "w") as f:
            json.dump(testset_eval_results, f)

        return testset_eval_results

    def _testset_eval(self, clf_name: str, clf: Any) -> dict:
        """Evaluate single classifier performance using test set data.

        Args:
            clf_name (str): Name of the classifier.
            clf (Any): Classifier instance.

        Returns:
            Dictionary containing the classifier's evaluation results.
        """

        _test_results = {}

        print("-" * 50)
        print(f"Evaluating {clf_name} ...")

        if isinstance(clf, consensus._Consensus):
            # Consensus method predicts the same number of samples as positive as true positive samples
            y_pred = clf.predict(self.X_test, limit_num=self.y_test.sum())
            y_proba = clf.predict_proba(self.X_test)[:, 1]
        else:
            y_pred = clf.predict(self.X_test)
            y_proba = clf.predict_proba(self.X_test)[:, 1]

        clf_eval_result = self._eval_single_clf(self.y_test, y_pred, y_proba)

        _test_results[clf_name] = clf_eval_result

        print("-" * 50)
        print(f"{clf_name} Evaluation Result:")
        print(f'AUC: {clf_eval_result["auc"]}')
        print(f'NEF: {clf_eval_result["nef"]}')
        print(f'F1: {clf_eval_result["f1"]}')
        print(f'Accuracy: {clf_eval_result["accuracy"]}')
        print(f'Precision: {clf_eval_result["precision"]}')
        print(f'Recall: {clf_eval_result["recall"]}')
        print("Confusion Matrix:")
        print(confusion_matrix(self.y_test, y_pred))
        print("-" * 50)

        return _test_results

    def _save_model(self, clf_name: str, clf: Any) -> None:
        """Save model as pickle file.

        Args:
            clf_name (str): Name of the classifier.
            clf (Any): Classifier instance to save.
        """
        file_path = os.path.join(MODELS_DIR, f"{clf_name}.pkl")
        with open(file_path, "wb") as f:
            pickle.dump(clf, f)
