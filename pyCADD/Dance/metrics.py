from math import ceil

import pandas as pd


def nef_score(y_true, y_score, percent: int = None):
    '''
    计算NEF值

    Parameters
    ----------
    y_true : array-like or pandas.Series
        样本的真实标签
    y_score : array-like or pandas.Series
        样本的预测得分(概率)
    percent : int
        早期富集百分率 取前 `percent%` 的样本计算NEF值
        默认为None 将使用阳性样本在总样本中的比率 Ra = actives / total
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
