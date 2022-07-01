import json
from typing import Any, Callable

import numpy as np
import pandas as pd
from pandas import DataFrame, Series
from pyCADD.Dance.algorithm.default_params import (GBT_DEFAULT_PARAMS,
                                                   LR_DEFAULT_PARAMS,
                                                   RF_DEFAULT_PARAMS)
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score,
                             precision_score, recall_score, roc_auc_score,
                             roc_curve)
from sklearn.model_selection import (GridSearchCV, RandomizedSearchCV,
                                     RepeatedStratifiedKFold, train_test_split)


def hyperparam_tuning(model, param_gird: dict, X: DataFrame, y: Series, scoring: str = 'roc_auc', cv: int = 5, n_jobs: int = -1, method: str = 'grid', save: bool = True, model_name: str = None):
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
    save : bool
        是否保存调优结果
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

    if save:
        params_file = 'best_params_%s.json' % (
            model_name if model_name else model.__class__.__name__)
        with open(params_file, 'w') as f:
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


class Matrix:
    '''
    结果矩阵
    '''

    def __init__(self, dataframe: DataFrame, test_size=0.25, random_seed=42) -> None:
        '''
        Parameters
        ----------
        dataframe : DataFrame
            数据DataFrame
        test_size : float
            测试集比例
        random_seed : int
            随机种子
        '''
        self.raw_data = dataframe
        self.data = self.raw_data.copy()
        self.test_size = test_size
        self.random_seed = random_seed

        self.train_data = None
        self.test_data = None

    @classmethod
    def from_pickle(cls, path, *args, **kwargs) -> 'Matrix':
        '''
        从pickle文件中读取结果矩阵

        Parameters
        ----------
        path : str
            pickle文件路径

        '''
        data = pd.read_pickle(path)
        return cls(data, *args, **kwargs)

    @classmethod
    def from_csv(cls, path, *args, **kwargs) -> 'Matrix':
        '''
        从csv文件中读取结果矩阵

        Parameters
        ----------
        path : str
            csv文件路径
        '''
        data = pd.read_csv(path)
        return cls(data, *args, **kwargs)

    @classmethod
    def from_splited_data(cls, train_data: DataFrame, test_data: DataFrame):
        '''
        从划分好的数据中读取结果矩阵

        Parameters
        ----------
        train_data : DataFrame
            训练集数据
        test_data : DataFrame
            测试集数据
        '''
        data = pd.concat([train_data, test_data], axis=0)
        cls_instance = cls(data)
        cls_instance.train_data = train_data
        cls_instance.test_data = test_data
        return cls_instance

    def split_train_test_data(self, test_size: float = None, random_seed: int = None, label_col: str = 'activity') -> tuple[DataFrame, DataFrame]:
        '''
        划分训练集和测试集数据
        '''
        if test_size is None:
            test_size = self.test_size
        if random_seed is None:
            random_seed = self.random_seed

        X, y = self.data.drop([label_col], axis=1), self.data[label_col]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_seed, stratify=y)
        self.train_data = pd.concat([X_train, y_train], axis=1)
        self.test_data = pd.concat([X_test, y_test], axis=1)

        return self.train_data, self.test_data

    def get_train_data(self, label_col: str = 'activity') -> DataFrame:
        '''
        获取训练集数据
        '''
        if self.train_data is None:
            self.split_train_test_data()
        return self.train_data

    def get_test_data(self, label_col: str = 'activity') -> DataFrame:
        '''
        获取测试集数据
        '''
        if self.test_data is None:
            self.split_train_test_data()
        return self.test_data


class Evaluator:
    '''
    模型性能评估器
    仅适用于非神经网络模型
    '''

    def __init__(self, matrix: Matrix, label_col: str = 'activity') -> None:
        '''
        Parameters
        ----------
        matrix : Matrix
            结果矩阵
        '''
        self.matrix = matrix
        self.train_data = matrix.get_train_data()
        self.test_data = matrix.get_test_data()
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
        '''
        默认GBT参数空间
        '''
        return GBT_DEFAULT_PARAMS

    @property
    def lr_default_params(self) -> dict:
        '''
        默认LR参数空间
        '''
        return LR_DEFAULT_PARAMS

    @property
    def rf_default_params(self) -> dict:
        '''
        默认RF参数空间
        '''
        return RF_DEFAULT_PARAMS

    @staticmethod
    def get_weights(y: Series) -> Series:
        '''
        计算不平衡数据集中各类标签的权重
        '''
        # 样本不平衡时，采用权重采样
        label_to_count = y.value_counts()
        _weights = []
        for label in y:
            _weights.append(1 / label_to_count[label])
        return _weights

    def load_params(self, path: str) -> dict:
        '''
        加载参数文件
        '''
        with open(path, 'r') as f:
            params = json.load(f)
        return params

    def save_params(self, path: str, params: dict) -> None:
        '''
        保存参数文件
        '''
        with open(path, 'w') as f:
            json.dump(params, f)

    def search_params(self, clf: Any, params_grid: dict, *args, **kwargs) -> dict:
        '''
        对模型实施最佳超参数搜索

        Parameters
        ----------
        clf : Any
            分类器
        params_grid : dict
            参数空间
        args : tuple
            其他参数 传递给hyperparam_tuning函数
        kwargs : dict
            其他参数 传递给hyperparam_tuning函数
        '''
        params = hyperparam_tuning(
            clf, params_grid, self.X_train, self.y_train, *args, **kwargs
        )
        return params

    def add_clf(self, clf: Any, clf_name: str = None) -> None:
        '''
        添加分类器实例到待评估的self.clfs字典属性中
        同名分类器将被覆盖

        Parameters
        ----------
        clf : Any
            分类器实例
        clf_name : str
            分类器名称 如为None则使用clf的__class__.__name__
        '''
        if clf_name is None:
            clf_name = clf.__class__.__name__

        self.clfs[clf_name] = clf

    def del_clf(self, clf_name: str) -> None:
        '''
        删除分类器实例
        '''
        del self.clfs[clf_name]

    def get_clf(self, clf_name: str) -> Any:
        '''
        使用名称索引 获取已添加的单个分类器实例
        '''
        return self.clfs[clf_name]

    def get_clfs_list(self) -> dict:
        '''
        获取分类器实例列表
        '''
        return self.clfs

    def print_classifier_info(self) -> None:
        '''
        打印已添加的分类器参数信息
        '''
        formatter = '{0:<25}{1:<30}{2:<40}'
        print('='*100)
        print(formatter.format('Classifier', 'Parameters', 'Values'))
        print('-'*100)
        for clf_name, clf in self.clfs.items():
            parameters = [(k, v) for k, v in clf.get_params().items()]
            print(formatter.format(clf_name, str(
                parameters[0][0]), str(parameters[0][1])))
            for k, v in parameters[1:]:
                print(formatter.format('', str(k), str(v)))
            print('-'*100)
        print('='*100)

    def repeat_cv(self, n_repeats: int = 30, k_folds: int = 4, random_seed: int = 42, score_func: Callable = roc_auc_score) -> dict:
        '''
        对所有添加的分类器执行多重交叉验证

        Parameters
        ----------
        n_repeats : int
            多重交叉验证的次数
        k_folds : int
            交叉验证的折数
        random_seed : int
            随机种子
        score_func : Callable
            评估函数
            默认为roc_auc_score

        Returns
        -------
        dict
            评估结果字典
            SCP为单构象结果
            clf_cv_results为交叉验证结果
        '''
        if len(self.clfs) == 0:
            raise ValueError('No classifier is added.')

        self.score_func_name = score_func.__name__

        print('\nRepeated Cross-Validation')
        print(f'Reapeat: {n_repeats} times')
        print(f'K-Folds: {k_folds} folds')
        print(f'Score Function: {score_func.__name__}')
        print('=' * 50)

        splits = _get_splits(self.X_train, self.y_train,
                             n_repeats, k_folds, random_seed)
        self.scp_score_df = _get_cv_scp_score(
            splits, self.X_train, self.y_train, score_func=score_func)
        best_scp_score = self.scp_score_df.max().max()
        mean_scp_score = self.scp_score_df.mean().mean()
        worst_scp_score = self.scp_score_df.min().min()
        std_scp_score = self.scp_score_df.values.flatten().std()

        self.cv_results['best_scp'] = best_scp_score
        self.cv_results['mean_scp'] = mean_scp_score
        self.cv_results['worst_scp'] = worst_scp_score
        self.cv_results['std_scp'] = std_scp_score

        print('SCP Best Score: %.4f' % best_scp_score)
        print('SCP Mean Score: %.4f' % mean_scp_score)
        print('SCP Worst Score: %.4f' % worst_scp_score)
        print('SCP Score std: %.4f' % std_scp_score)

        self.cv_results['clf_cv_results'] = {}
        for clf_name, clf in self.clfs.items():
            print('-' * 50)
            print(f'Validating {clf_name} ...')
            repeat_cv_result = _repeat_cross_validation(
                clf, splits, self.X_train, self.y_train, score_func)
            self.cv_results['clf_cv_results'][clf_name] = repeat_cv_result
            print('%s mean Score: %.4f' % (
                clf_name,
                np.mean(repeat_cv_result)
            ))
            print('-' * 50)

        print('=' * 50)
        print(f'{n_repeats} x {k_folds} cross-validation finished.')
        print('=' * 50)

        with open('cv_results.json', 'w') as f:
            json.dump(self.cv_results, f)

        return self.cv_results

    def print_cv_results(self) -> None:
        '''
        打印交叉验证结果
        '''
        if len(self.cv_results) == 0:
            raise ValueError('No cross-validation is performed.')

        formatter = '{0:<25}{1:<30}'

        print(f'\nScore: {self.score_func_name}')
        print('='*50)
        print(formatter.format('Best SCP', '%.4f' %
              self.cv_results['best_scp']))
        print(formatter.format('Mean SCP', '%.4f' %
              self.cv_results['mean_scp']))
        print(formatter.format('Worst SCP', '%.4f' %
              self.cv_results['worst_scp']))
        print(formatter.format('Std SCP', '%.4f' % self.cv_results['std_scp']))
        print('-' * 50)
        print(formatter.format('Classifier', 'Mean Score'))
        print(formatter.format('-'*25, '-'*30))
        for clf_name, clf_cv_results in self.cv_results['clf_cv_results'].items():
            print(formatter.format(clf_name, '%.4f' %
                  np.mean(clf_cv_results)))  # 打印每个分类器的平均AUC
        print('=' * 50)

    @staticmethod
    def _eval_single_clf(y_test: Series, y_pred: Series, y_proba: Series):
        '''
        评估单个分类器的性能

        Parameters
        ----------
        y_test : Series
            测试集的标签
        y_pred : Series
            预测集的标签
        y_proba : Series
            预测集的概率

        Returns
        -------
        score_dict : dict
            分类器的评估结果
        '''
        auc = roc_auc_score(y_test, y_proba)
        f1 = f1_score(y_test, y_pred)
        acc = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        fpr, tpr, thresholds = roc_curve(y_test, y_proba)

        return {
            'auc': auc,
            'f1': f1,
            'accuracy': acc,
            'precision': precision,
            'recall': recall,
            'fpr': list(fpr),
            'tpr': list(tpr),
        }

    def testset_eval(self):
        '''
        使用测试集数据评估性能

        Returns
        -------
        testset_results : dict
            分类器的评估结果
        '''
        if len(self.clfs) == 0:
            raise ValueError('No classifier is added.')

        testset_eval_results = {}
        print('\nTestset Evaluation')
        print('='*50)
        for clf_name, clf in self.clfs.items():
            print('-' * 50)
            print(f'Evaluating {clf_name} ...')

            y_pred = clf.predict(self.X_test)
            y_proba = clf.predict_proba(self.X_test)[:, 1]
            clf_eval_result = self._eval_single_clf(
                self.y_test, y_pred, y_proba)
            testset_eval_results[clf_name] = clf_eval_result

            print('-' * 50)
            print(f'{clf_name} Evaluation Result:')
            print(f'AUC: {clf_eval_result["auc"]}')
            print(f'F1: {clf_eval_result["f1"]}')
            print(f'Accuracy: {clf_eval_result["accuracy"]}')
            print(f'Precision: {clf_eval_result["precision"]}')
            print(f'Recall: {clf_eval_result["recall"]}')
            print('Confusion Matrix:')
            print(confusion_matrix(self.y_test, y_pred))
            print('-' * 50)
        print('='*50)
        print('Testset Evaluation Finished.')
        print('='*50)

        with open('testset_eval_results.json', 'w') as f:
            json.dump(testset_eval_results, f)

        return testset_eval_results
