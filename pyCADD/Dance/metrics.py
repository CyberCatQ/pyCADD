from math import ceil

import pandas as pd


def nef_score(y_true, y_score, percent: int = None) -> float:
    """Calculate Normalized Enrichment Factor (NEF) score.

    NEF measures the enrichment of active compounds in the top-ranked subset
    compared to random selection.

    Args:
        y_true (array-like or pd.Series): True binary labels (0 or 1) for samples.
        y_score (array-like or pd.Series): Predicted scores/probabilities for samples.
        percent (int, optional): Early enrichment percentage. Calculate NEF for top `percent%` of samples.
            If None, uses the ratio of actives to total samples (Ra = actives / total).

    Returns:
        NEF score as a float value.
    """

    if isinstance(y_true, pd.Series):
        y_true = y_true.values
    if isinstance(y_score, pd.Series):
        y_score = y_score.values

    _df = pd.DataFrame(columns=["y_true", "y_score"])
    _df["y_true"] = y_true
    n_true = _df["y_true"].sum()

    _df["y_score"] = y_score
    _df.sort_values(by="y_score", ascending=False, inplace=True)

    Ra = n_true / _df["y_true"].shape[0]
    percent = percent / 100 if percent is not None else Ra
    top_n = ceil(percent * _df["y_true"].shape[0])

    true_in_top_n = _df["y_true"].iloc[:top_n].sum()

    nef = true_in_top_n / min(top_n, n_true)

    return nef
