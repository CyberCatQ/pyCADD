from math import ceil

import pandas as pd


def nef_score(y_true, y_score, percent: int = None):
    '''
    计算NEF值

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.
    y_score : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.
    percent : float
        EF percentage, i.e. EF_1% is enrichment factor for first percent of
        given samples. This function assumes that results are already sorted and
    '''

    if isinstance(y_true, pd.Series):
        y_true = y_true.values
    if isinstance(y_score, pd.Series):
        y_score = y_score.values

    _df = pd.DataFrame(columns=['y_true', 'y_score'])
    _df['y_true'] = y_true
    n_true = _df['y_true'].sum()

    _df['y_score'] = y_score
    _df.sort_values(by='y_score', ascending=False, inplace=True)

    Ra = n_true / _df['y_true'].shape[0]
    percent = percent / 100 if percent is not None else Ra
    top_n = ceil(percent * _df['y_true'].shape[0])

    true_in_top_n = _df['y_true'].iloc[:top_n].sum()

    nef = true_in_top_n / min(top_n, n_true)

    return nef
