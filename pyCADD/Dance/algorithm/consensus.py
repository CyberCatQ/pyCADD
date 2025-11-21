from typing import Literal

import numpy as np
import pandas as pd
from pandas import DataFrame, Series


class _Consensus:
    """Base class for consensus algorithms.

    New consensus methods (formula-based methods) should inherit from this class.
    Provides a common interface for consensus-based molecular screening approaches.
    """

    def __init__(self, lower_is_better: bool = False) -> None:
        """Initialize consensus algorithm.

        Args:
            lower_is_better (bool): Whether data features represent descending indicators
                (e.g., lower binding energy is better).
        """
        self.result = None
        self.lower_is_better = lower_is_better

        self.params = {
            "class": "Consensus Strategy",
            "lower_is_better": self.lower_is_better,
        }

    def get_params(self) -> dict:
        """Get parameters for sklearn compatibility.

        Returns:
            Dictionary containing algorithm parameters.
        """
        return self.params

    def fit(self) -> None:
        """Fit the consensus algorithm.

        To be overridden by subclasses with specific implementation.

        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError

    def predict(self, X: DataFrame, y: Series = None, limit_num: int = None) -> Series:
        """Make predictions using the consensus algorithm.

        Args:
            X (DataFrame): Input feature data.
            y (Series, optional): Input data labels (used to determine number of positives if limit_num is None).
            limit_num (int, optional): Number of top-ranked results to predict as positive.
                If None, uses the number of positive samples in y.

        Returns:
            Series containing binary predictions (0 or 1).

        Raises:
            RuntimeError: If result is None when trying to predict.
            ValueError: If y is required but not provided when limit_num is None.
        """
        # Consensus methods don't actually perform fitting
        # but calculate predictions based on formulas

        if X is not None:
            self.fit(X, y)
        if self.result is None:
            raise RuntimeError(
                f"Result is None while using {self.__class__.__name__} method to predict. X is required."
            )

        _predict_df = self.result.reset_index()
        _predict_df.columns = ["index", "probability"]
        _predict_df.sort_values(by="probability", ascending=False, inplace=True)
        _predict_df["prediction"] = np.nan

        if limit_num is None:
            if y is None:
                raise ValueError("y is required when limit_num is None.")
            print(
                f"No limit_num specified, TOP {y.sum()} results (same as the number of true positive samples) will be predicted as positive."
            )
            limit_num = y.sum()

        _predict_df.iloc[:limit_num, -1] = 1
        _predict_df.iloc[limit_num:, -1] = 0

        _predict_df.sort_index(inplace=True)
        _predict_df.set_index("index", inplace=True)

        return _predict_df.loc[:, "prediction"]

    def predict_proba(self, X: DataFrame = None, y: Series = None) -> np.ndarray:
        """Make probability predictions using the consensus algorithm.

        Args:
            X (DataFrame, optional): Input feature data.
            y (Series, optional): Input data labels (for compatibility).

        Returns:
            Numpy array with shape (n_samples, 2) containing probabilities.
            Column 0: probability of negative class
            Column 1: probability of positive class (consensus scores)

        Raises:
            RuntimeError: If result is None when trying to predict.
        """
        # Consensus methods don't actually perform fitting
        # but calculate results based on formulas - no actual probabilities
        # Directly return calculated results as Numpy array

        if X is not None:
            self.fit(X, y)
        if self.result is None:
            raise RuntimeError(
                f"Result is None while using {self.__class__.__name__} method to predict. X is required."
            )

        # Array format: [[name, value], [name, value], ...]
        # For sklearn predict_proba compatibility: dimension 0 = sample name (prediction in sklearn), dimension 1 = value (probability in sklearn)
        # Dimension 1 is the primary sorting dimension and used for generating AUC etc.
        return pd.DataFrame(self.result).reset_index().to_numpy()


class Average(_Consensus):
    """Arithmetic mean consensus model.

    Calculates arithmetic mean of molecular descriptor values across conformations.
    Arithmetic mean is defined as: mean = sum(x_i) / count
    """

    def __init__(self, lower_is_better: bool = False) -> None:
        """Initialize arithmetic mean consensus algorithm.

        Args:
            lower_is_better (bool): Whether lower values indicate better performance.
        """
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({"method": "Arithmetic mean"})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True) -> None:
        """Calculate arithmetic mean across conformations.

        Args:
            X (DataFrame): Feature data containing molecular descriptors across conformations.
            y (Series, optional): Data labels (not used in calculation, exists for compatibility).
            ignore_nan (bool): Whether to ignore missing values (all 0 values are treated as missing).
        """
        # In preprocessing stage, all NaN in X are filled with 0
        # In calculation stage, all 0s in X are refilled as NaN and ignored during mean calculation
        X = X.replace(0, np.nan) if ignore_nan else X

        # mean calculation automatically ignores NaN
        self.result = X.mean(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result


class Mean(Average):
    """Arithmetic mean consensus model.

    Alias for the Average class providing the same arithmetic mean functionality.
    """


class Geo_Average(_Consensus):
    """Geometric mean consensus model.

    Calculates geometric mean of molecular descriptor values across conformations.
    Geometric mean is defined as: geo_mean = (product(x_i))^(1/n)
    """

    def __init__(self, lower_is_better: bool = False) -> None:
        """Initialize geometric mean consensus algorithm.

        Args:
            lower_is_better (bool): Whether lower values indicate better performance.
        """
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({"method": "Geometric mean"})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True) -> None:
        """Calculate geometric mean across conformations.

        Takes absolute values for calculation, then applies sign.

        Args:
            X (DataFrame): Feature data containing molecular descriptors across conformations.
            y (Series, optional): Data labels (not used in calculation, exists for compatibility).
            ignore_nan (bool): Whether to ignore missing values (all 0 values are treated as missing).
        """
        # In preprocessing stage, all NaN in X are filled with 0
        # In calculation stage, all 0s in X are refilled as NaN and ignored during calculation
        X = X.replace(0, np.nan) if ignore_nan else X

        # Take absolute values for calculation, then apply sign
        self.result = Series(
            pow(
                # Absolute value continuous product
                X.prod(axis=1).abs(),
                # Take nth root where n is the number of valid values
                1 / X.notna().sum(axis=1),
            ),
            name="GEO",
        )

        self.result = self.result if self.lower_is_better else -1 * self.result


class GeoMean(Geo_Average):
    """Geometric mean consensus model.

    Alias for the Geo_Average class providing the same geometric mean functionality.
    """


class Minimum(_Consensus):
    """Minimum value consensus model.

    Selects the minimum value across conformations for each molecule.
    """

    def __init__(self, lower_is_better: bool = False) -> None:
        """Initialize minimum value consensus algorithm.

        Args:
            lower_is_better (bool): Whether lower values indicate better performance.
        """
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({"method": "Minimum"})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True) -> None:
        """Calculate minimum values across conformations.

        Args:
            X (DataFrame): Feature data containing molecular descriptors across conformations.
            y (Series, optional): Data labels (not used in calculation, exists for compatibility).
            ignore_nan (bool): Whether to ignore missing values (all 0 values are treated as missing).
        """
        X = X.replace(0, np.nan) if ignore_nan else X
        self.result = X.min(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result


class Maximum(_Consensus):
    """Maximum value consensus model.

    Selects the maximum value across conformations for each molecule.
    """

    def __init__(self, lower_is_better: bool = False) -> None:
        """Initialize maximum value consensus algorithm.

        Args:
            lower_is_better (bool): Whether lower values indicate better performance.
        """
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({"method": "Maximum"})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True) -> None:
        """Calculate maximum values across conformations.

        Args:
            X (DataFrame): Feature data containing molecular descriptors across conformations.
            y (Series, optional): Data labels (not used in calculation, exists for compatibility).
            ignore_nan (bool): Whether to ignore missing values (all 0 values are treated as missing).
        """
        X = X.replace(0, np.nan) if ignore_nan else X
        self.result = X.max(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result


# Functional interfaces


def average(data: DataFrame, method: Literal["ave", "geo"] = "ave") -> Series:
    """Calculate average values across molecular conformations.

    Provides functional interface for average calculations:
    - 'ave': Arithmetic mean
    - 'geo': Geometric mean

    Args:
        data (DataFrame): DataFrame containing molecular data to calculate averages from.
        method (str): Average calculation method ('ave' for arithmetic mean | 'geo' for geometric mean).

    Returns:
        Series containing average values.

    Raises:
        RuntimeError: If method is not 'ave' or 'geo'.
    """
    if method.upper() == "AVE":
        return Series(data.mean(axis=1), name="AVE")
    elif method.upper() == "GEO":
        # First eliminate signs then calculate geometric mean
        # Finally add sign
        return Series(-pow(data.prod(axis=1).abs(), 1 / data.notna().sum(axis=1)), name="GEO")
    else:
        raise RuntimeError("Invalid average value method.")


def minimum(data: DataFrame) -> Series:
    """Extract minimum values (best scores) from DataFrame columns.

    Useful for extracting optimal scores such as best docking scores.

    Args:
        data (DataFrame): DataFrame containing molecular data to extract minimums from.

    Returns:
        Series containing minimum values for each row.
    """
    return Series(data.min(axis=1), name="MIN")


def maximum(data: DataFrame) -> Series:
    """Extract maximum values from DataFrame columns.

    Args:
        data (DataFrame): DataFrame containing molecular data to extract maximums from.

    Returns:
        Series containing maximum values for each row.
    """
    return Series(data.max(axis=1), name="MAX")


def std(data: DataFrame, axis: int = 1) -> Series:
    """Calculate standard deviation.

    Args:
        data (DataFrame): DataFrame containing molecular data to calculate standard deviation from.
        axis (int): Calculation axis (0 for rows, 1 for columns).

    Returns:
        Series containing standard deviation values.
    """
    return Series(data.std(axis=axis), name="std")
