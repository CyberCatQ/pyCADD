import json
import os
from typing import Any, Callable

import pandas as pd
from pandas import DataFrame, Series

from sklearn.metrics import roc_auc_score
from sklearn.model_selection import (GridSearchCV, RandomizedSearchCV,
                                     RepeatedStratifiedKFold)


def hyperparam_tuning(model, param_gird: dict, X: DataFrame, y: Series, scoring: str = 'roc_auc', cv: int = 5, n_jobs: int = -1, method: str = 'grid', save_dir: str = None, model_name: str = None):
    '''
    超参数调优

    Parameters
    ----------
    model : object
        需要调优的模型
    param_grid : dict
        超参数网格
    X : DataFrame
        训练集
    y : Series
        训练集标签
    scoring : str
        评估标准
    cv : int
        训练集分割数
    n_jobs : int
        训练进程数
    method : str
        调优方法
        * grid: 网格搜索
        * random: 随机搜索
    save_dir : str | None
        参数文件保存路径 为None则不保存
    model_name : str
        模型名称

    Return
    ----------
    dict
        调优后的模型参数
    '''

    # 调优方法
    if method == 'grid':
        cv_search = GridSearchCV(
            estimator=model, param_grid=param_gird, scoring=scoring, cv=cv, n_jobs=n_jobs)
    elif method == 'random':
        cv_search = RandomizedSearchCV(
            estimator=model, param_distributions=param_gird, scoring=scoring, cv=cv, n_jobs=n_jobs)
    else:
        raise ValueError('Invalid method: %s' % method)

    cv_search.fit(X, y)
    best_params = cv_search.best_params_
    best_score = cv_search.best_score_
    print('Best parameters: %s' % best_params)
    print('AUC score: %s' % best_score)

    if save_dir is not None:
        params_file = 'best_params_%s.json' % (
            model_name if model_name else model.__class__.__name__)
        with open(os.path.join(save_dir, params_file), 'w') as f:
            json.dump(best_params, f)
    return best_params


def _score(model, X_train: DataFrame, X_test: DataFrame, y_train: Series, y_test: Series, score_func: Callable):
    '''
    模型评估

    Parameters
    ----------
    model : object
        需要评估的模型
    X_train : DataFrame
        训练集
    X_test : DataFrame
        测试集
    y_train : Series
        训练集标签
    y_test : Series
        测试集标签
    score_func : callable
        评估方法
        即 sklearn.metrics 函数

    Return
    ----------
    tuple(float, float)
        训练集和测试集的ROC-AUC值
    '''

    # 适用于sklearn模型与Consensus模型
    model.fit(X_train, y_train)

    # 概率预测
    y_train_predicted = model.predict_proba(X_train)[:, 1]
    y_test_predicted = model.predict_proba(X_test)[:, 1]

    # 按照概率排序并以评分方法计算评分
    train_score = score_func(y_train, y_train_predicted)
    test_score = score_func(y_test, y_test_predicted)

    return train_score, test_score


def calc_scp_score(X: DataFrame, y_true: Series, lower_is_better: bool = True, score_func: Callable = roc_auc_score):
    '''
    获取单构象Performance计算ROC-AUC值

    Parameters
    ----------
    X : DataFrame
        数据特征
    y_true : Series
        标签
    lower_is_better : bool
        评估指标是否为下降性指标
    score_func : callable
        评估方法
        即 sklearn.metrics 函数 默认为ROC-AUC

    Return
    ----------
    dict
        {'PDB': SCP AUC}
    '''
    results_dict = {}

    # 所有构象的Performance
    comformations = X.columns

    # 对每个构象进行评估
    for comformation in comformations:
        if lower_is_better:
            scp = X[comformation] * -1
        else:
            scp = X[comformation]
        results_dict[comformation] = score_func(y_true, scp)

    return results_dict


def _get_cv_scp_score(splits: list, X: DataFrame, y: Series, score_func: Callable = roc_auc_score):
    '''
    在交叉测试cv中报告每一次cv的测试集部分的评分结果集

    Parameters
    ----------
    splits : list   
        索引列表(训练集, 测试集)
    X : DataFrame
        总数据集 依据索引从中划分训练集和测试集
    y : Series
        标签
    score_func : callable
        评估方法
        即 sklearn.metrics 函数 默认为ROC-AUC

    Return
    ----------
    DataFrame
        SCP策略测试AUC矩阵
    '''

    scp_aucs = []

    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        scp_score = calc_scp_score(X_test, y_test, True, score_func=score_func)
        scp_aucs.append(scp_score)

    cv_scp_auc_results = pd.DataFrame(scp_aucs)
    return cv_scp_auc_results


def _get_splits(X: DataFrame, y: Series, n_repeats: int = 30, k_folds: int = 4, random_state: int = 42):
    '''
    创建数据索引 
    即训练集与测试/验证集的索引号 用于多重交叉验证

    Parameters
    ----------
    X : DataFrame
        数据特征
    y : Series
        数据标签
    n_repeats : int
        重复划分次数 n
    k_folds : int
        数据划分折数 k
    random_state : int
        随机种子 相同的种子会产生相同的划分

    Returns
    -------
    list
        索引列表(训练集, 测试/验证集)
    '''
    cv = RepeatedStratifiedKFold(
        n_repeats=n_repeats, n_splits=k_folds, random_state=random_state)
    return [*cv.split(X, y)]


def _repeat_cross_validation(model: Any, splits: list, X: DataFrame, y: Series, score_func: Callable):
    '''
    对模型实施多重交叉验证

    Parameters
    ----------
    model : object
        需要评估的模型
    splits : list
        索引列表(训练集, 测试集)
        由RepeatedStratifiedKFold创建的索引
    X : DataFrame
        总数据集 依据索引从中划分训练集和测试集
    y : Series  
        标签
    score_func : callable
        评估方法

    Return
    ----------
    list
        n x m 次模型交叉验证结果
    '''
    validation_results = []

    for train_index, test_index in splits:
        # 单次交叉验证中取得训练集和测试集数据
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        # 拟合模型并进行评估打分
        train_score, test_score = _score(
            model, X_train, X_test, y_train, y_test, score_func)
        validation_results.append(test_score)

    return validation_results
