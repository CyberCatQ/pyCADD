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
from sklearn.metrics import (accuracy_score, confusion_matrix, f1_score,
                             precision_score, recall_score, roc_auc_score,
                             roc_curve)
from sklearn.model_selection import train_test_split

logger = logging.getLogger('pyCADD.Dance.base')

DATA_MINING_DIR = os.path.join(os.getcwd(), 'dm_result')
DATA_PICKLE_DIR = os.path.join(DATA_MINING_DIR, 'data_pickle')
DATA_CSV_DIR = os.path.join(DATA_MINING_DIR, 'data_csv')
FIGURES_DIR = os.path.join(DATA_MINING_DIR, 'figures')
MODELS_DIR = os.path.join(DATA_MINING_DIR, 'models_pickle')
PARAMS_DIR = os.path.join(DATA_MINING_DIR, 'models_params')


class Dancer:
    '''
    Data analyzer for CADD
    数据预处理及分析器
    '''

    def __init__(self) -> None:
        # 创建必要的目录
        makedirs_from_list(self._required_dirs)
        self.current_datasets = []
        self.merged_data = None

    @property
    def _required_dirs(self):
        return [
            DATA_MINING_DIR,
            DATA_PICKLE_DIR,
            DATA_CSV_DIR,
            FIGURES_DIR,
            MODELS_DIR,
            PARAMS_DIR
        ]

    def _add_dataset(self, csv_path: str, activity: Literal['positive', 'negative'] = None, *args, **kwargs) -> None:
        '''
        增加样本数据集

        Parameters
        ----------
        csv_path : str
            样本数据集的路径
        activity : Literal['positive', 'negative']
            样本数据集的活性类型 None则为空
        *args, **kwargs
            其他参数 传入pd.read_csv()
        '''
        dataset = pd.read_csv(csv_path, *args, **kwargs)
        dataset = dataset.sort_index(axis=1)
        # 丢弃所有特征缺失的行
        dataset.dropna(axis=0, how='all', inplace=True)
        self.current_datasets.append(
            {'dataset': dataset, 'activity': activity})

    def add_pos_dataset(self, csv_path: str) -> None:
        '''
        增加正样本数据集

        Parameters
        ----------
        csv_path : str
            正样本数据集的csv文件路径
        '''
        self._add_dataset(csv_path, 'positive', index_col=0)

    def add_neg_dataset(self, csv_path: str) -> None:
        '''
        增加负样本数据集

        Parameters
        ----------
        csv_path : str
            负样本数据集的csv文件路径
        '''
        self._add_dataset(csv_path, 'negative', index_col=0)

    def add_dataset(self, csv_path: str, *args, **kwargs) -> None:
        '''
        增加数据集

        Parameters
        ----------
        csv_path : str
            数据集的csv文件路径
        '''
        self._add_dataset(csv_path, index_col=0, *args, **kwargs)

    def _add_label_col(self) -> None:
        '''
        增加活性对应的标签列
        '''
        for dataset in self.current_datasets:
            if dataset['activity'] is None:
                dataset['dataset']['activity'] = 'Undefined'
            elif dataset['activity'] == 'positive':
                dataset['dataset']['activity'] = 1
            elif dataset['activity'] == 'negative':
                dataset['dataset']['activity'] = 0

    def _concat_datasets(self) -> DataFrame:
        '''
        合并数据集 保存在self.merged_data属性中 并返回合并数据
        '''
        if len(self.current_datasets) == 0:
            raise ValueError('No datasets added.')
        self._add_label_col()
        self.merged_data = pd.concat(
            [dataset['dataset'] for dataset in self.current_datasets])
        return self.merged_data

    def get_merged_data(self) -> DataFrame:
        '''
        返回合并数据集
        '''
        return self.merged_data

    def _fill_nan(self, fill_na_value: Any = 0, dataset: DataFrame = None, inplace: bool = False) -> None:
        '''
        返回填充缺失值的数据集

        Parameters
        ----------
        fill_na_value : Any
            填充缺失值的值 默认为0
        dataset : DataFrame
            填充缺失值的数据集 默认为self.merged_data
        '''
        if dataset is None:
            dataset = self.merged_data
        if inplace:
            dataset.fillna(fill_na_value, inplace=True)
            return None
        else:
            return dataset.fillna(fill_na_value)

    def prepare_data(self, fill_nan: bool = True, *args, **kwargs) -> None:
        '''
        准备数据集

        Parameters
        ----------
        fill_nan : bool
            是否填充缺失值 默认为True
        *args, **kwargs
            其他参数 传入self._fill_nan()
            value可指定填充值 默认0
        '''
        self._concat_datasets()
        if fill_nan:
            self._fill_nan(inplace=True, *args, **kwargs)

    def save_pickle(self, file_name: str, dataset: DataFrame = None) -> None:
        '''
        保存数据集
        '''
        dataset = dataset if dataset is not None else self.merged_data
        dataset.to_pickle(os.path.join(DATA_PICKLE_DIR, file_name))

    def save_csv(self, file_name: str, dataset: DataFrame = None) -> None:
        '''
        保存数据集
        '''
        dataset = dataset if dataset is not None else self.merged_data
        dataset.to_csv(os.path.join(DATA_CSV_DIR, file_name))

    def save(self, file_name: str, dataset: DataFrame = None) -> None:
        '''
        保存数据集
        '''
        dataset = dataset if dataset is not None else self.merged_data
        if file_name.endswith('.csv'):
            self.save_csv(file_name, dataset)
        elif file_name.endswith('.pickle') or file_name.endswith('.pkl'):
            self.save_pickle(file_name, dataset)
        else:
            raise ValueError('file_name must end with .csv or .pickle')


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

    def split_train_test_data(self, test_size: float = None, random_seed: int = None, label_col: str = 'activity') -> tuple:
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
        return default_params.GBT_DEFAULT_PARAMS

    @property
    def lr_default_params(self) -> dict:
        '''
        默认LR参数空间
        '''
        return default_params.LR_DEFAULT_PARAMS

    @property
    def rf_default_params(self) -> dict:
        '''
        默认RF参数空间
        '''
        return default_params.RF_DEFAULT_PARAMS

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

    def get_lr_default_params(self) -> dict:
        '''
        获取默认LR参数空间
        '''
        return self.lr_default_params

    def get_rf_default_params(self) -> dict:
        '''
        获取默认RF参数空间
        '''
        return self.rf_default_params

    def get_gbt_default_params(self) -> dict:
        '''
        获取默认GBT参数空间
        '''
        return self.gbt_default_params

    def load_params(self, path: str) -> dict:
        '''
        加载参数文件

        Parameters
        ----------
        path : str
            参数文件路径

        Returns
        -------
        params : dict
            参数字典
        '''
        with open(path, 'r') as f:
            params = json.load(f)
        return params

    def save_params(self, file_name: str, params: dict) -> None:
        '''
        保存参数文件

        Parameters
        ----------
        file_name : str
            参数文件名
        params : dict
            参数字典
        '''
        with open(os.path.join(PARAMS_DIR, file_name), 'w') as f:
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
        params = core.hyperparam_tuning(
            clf, params_grid, self.X_train, self.y_train, save_dir=PARAMS_DIR, *args, **kwargs
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

    def get_clfs_dict(self) -> dict:
        '''
        获取分类器实例字典
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

        splits = core._get_splits(self.X_train, self.y_train,
                                  n_repeats, k_folds, random_seed)
        self.scp_score_df = core._get_cv_scp_score(
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
            repeat_cv_result = core._repeat_cross_validation(
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

        with open(f'cv_results_{self.score_func_name}.json', 'w') as f:
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
        print(formatter.format('-' * 50))
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
        nef = nef_score(y_test, y_proba)
        f1 = f1_score(y_test, y_pred)
        acc = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        fpr, tpr, thresholds = roc_curve(y_test, y_proba)

        return {
            'auc': auc,
            'nef': nef,
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
            clf.fit(self.X_train, self.y_train)                             # 训练模型
            testset_eval_results.update(self._testset_eval(clf_name, clf))  # 评估模型
            self._save_model(clf_name, clf)                                 # 保存模型
            
        print('='*50)
        print('Testset Evaluation Finished.')
        print('='*50)

        with open('testset_eval_results.json', 'w') as f:
            json.dump(testset_eval_results, f)

        return testset_eval_results

    def _testset_eval(self, clf_name, clf) -> dict:
        '''
        利用测试集数据测试单个分类器的性能

        Parameters
        ----------
        clf_name : str  
            分类器名称
        clf : Classifier
            分类器

        Returns
        -------
        dict
            分类器的评估结果

        '''

        _test_results = {}

        print('-' * 50)
        print(f'Evaluating {clf_name} ...')

        if isinstance(clf, consensus._Consensus):
            # 共识性方法预测与真阳性样本相同数量的样本为阳性
            y_pred = clf.predict(self.X_test, limit_num=self.y_test.sum())
            y_proba = clf.predict_proba(self.X_test)[:, 1]
        else:
            y_pred = clf.predict(self.X_test)
            y_proba = clf.predict_proba(self.X_test)[:, 1]

        clf_eval_result = self._eval_single_clf(
                self.y_test, y_pred, y_proba)

        _test_results[clf_name] = clf_eval_result

        print('-' * 50)
        print(f'{clf_name} Evaluation Result:')
        print(f'AUC: {clf_eval_result["auc"]}')
        print(f'NEF: {clf_eval_result["nef"]}')
        print(f'F1: {clf_eval_result["f1"]}')
        print(f'Accuracy: {clf_eval_result["accuracy"]}')
        print(f'Precision: {clf_eval_result["precision"]}')
        print(f'Recall: {clf_eval_result["recall"]}')
        print('Confusion Matrix:')
        print(confusion_matrix(self.y_test, y_pred))
        print('-' * 50)

        return _test_results
    
    def _save_model(self, clf_name, clf):
        '''
        保存模型为pkl文件

        Parameters
        ----------
        clf_name : str
            分类器名称
        clf : Classifier
            分类器
        '''
        file_path = os.path.join(MODELS_DIR, f'{clf_name}.pkl')
        with open(file_path, 'wb') as f:
            pickle.dump(clf, f)