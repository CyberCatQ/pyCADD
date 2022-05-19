import json
import logging
import os
import pickle
from copy import deepcopy
from time import sleep

import hiddenlayer as hl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from pandas import DataFrame, Series
from pyCADD.Dance.algorithm import (Average, Consensus, Geo_Average, Minimum,
                                    MyMLP)
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, classification_report,
                             confusion_matrix, f1_score, precision_score,
                             recall_score, roc_auc_score, roc_curve)
from sklearn.model_selection import (GridSearchCV, RandomizedSearchCV,
                                     RepeatedStratifiedKFold, train_test_split)
from sklearn.naive_bayes import GaussianNB
from torch import nn
from torch.utils.data import DataLoader, TensorDataset
from torch.utils.data.sampler import WeightedRandomSampler
from xgboost import XGBClassifier

logger = logging.getLogger(__name__)


def read_matrix(file_path: str):
    '''
    读取数据矩阵

    Parameters
    ----------
    file_path : str
        数据矩阵文件路径
    Return
    ----------
    DataFrame
        pandas DF数据对象
    '''
    raw_data = pd.read_csv(file_path, index_col=0)
    return raw_data


'''
def read_docking_data(raw_data: DataFrame, label_col: str):
    
    提取对接分数部分
    Parameters
    ---------
    raw_data : DataFrame
        原始数据对象
    label_col : str
        阳性标签数据列名

    Return
    ---------
    DataFrame
        纯对接分数数据对象
    

    return raw_data.drop(label_col, axis=1)
'''


def _standrad_label(data: DataFrame, label_col: str, positive: ... = 'origin'):
    '''
    标准化标签为二进制

    Parameters
    ----------
    data : DataFrame
        待标准化的数据
    label_col : str
        标签列名
    positive : str | int | list
        正样本标签

    Return
    ----------
    DataFrame
        标准化后的数据
    '''
    if isinstance(positive, str) or isinstance(positive, int):
        positive = [positive]

    # 标准化阳性标签为1
    for _label in positive:
        data[label_col].replace(_label, value=1, inplace=True)
    # object转int
    data[label_col] = data[label_col].astype(int)

    return data


def split_data(data: DataFrame, label_col: str = None, preprocess: bool = True, positive: str = 'origin'):
    '''
    拆分数据与标签

    Parameters
    ----------
    data : DataFrame
        待拆分数据
    label_col : str
        标签列名
    preprocess : bool
        是否预处理数据(默认为True)
        填充Nan值为0, 并且将标签二进制化
    positive : str | int | list
        正样本标签(preprosses为True时有效)

    Return  
    ----------
    DataFrame, Series
        拆分后的数据和标签
    '''
    if preprocess:
        data.fillna(0, inplace=True)
        data = _standrad_label(data, label_col, positive)
    return data.drop(label_col, axis=1), data[label_col]


def merge(data_list: list):
    '''
    合并Series

    Parameters
    ----------
    data_list : list
        包含需要合并的Series的列表
    '''
    return pd.concat(data_list, axis=1)


def _format_data(data: DataFrame, label_col: str, pos_label, score_name: str = 'Docking_Score'):
    '''
    为生成统计图预处理矩阵数据
    '''

    # 标签列二进制化
    data = _standrad_label(data, label_col, pos_label)

    total_data = pd.DataFrame(columns=['PDB', score_name, label_col])

    for index in data.columns[:-1]:
        _data = data[[index, label_col]]
        _data = _data.rename(columns={index: score_name})
        _data.loc[:, 'PDB'] = index
        total_data = pd.concat([total_data, _data], ignore_index=True)

    return total_data


def get_roc(X: DataFrame, y_ture: Series, save: bool = False, lower_is_better: bool = True):
    '''
    ROC曲线下面积
    依据label列作为标签
    自动计算DataFrame的所有列的ROC曲线下面积AUC值

    Parameters
    ----------
    X : DataFrame  
        待计算的数据
    y_ture : Series
        标签列
    save : bool
        是否存储ROC曲线图片
    lower_is_better : bool
        是否为负样本的ROC曲线(默认为True)

    Return
    ---------
    Series
        曲线下面积AUC数据
    '''

    auc_dict = {}

    plt.figure(figsize=(10, 10), dpi=300.0)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('%s ROC curve' % os.path.basename(os.getcwd()))

    for index in X.columns:
        if lower_is_better:
            X[index] = X[index] * -1
        fpr, tpr, thersholds = roc_curve(y_ture, X[index])
        auc = roc_auc_score(y_ture, X[index])
        auc_dict[index] = auc
        plt.plot(fpr, tpr, label='%s (area = %.2f)' % (index, auc))
    plt.plot([0, 1], [0, 1], linestyle='--')
    plt.legend(loc='lower right')

    cwd_name = os.path.basename(os.getcwd())
    if save:
        plt.savefig('%s-ROC.jpg' % cwd_name)

    plt.show()

    return Series(auc_dict, name='AUC')


def get_scatter(data: DataFrame, label_col: str, pos_label, score_name: str = 'Docking_Score', save: bool = False):
    '''
    生成分布散点图

    Parameters
    ----------
    data : DataFrame
        待计算数据
    label : str
        阳性标签列名
    pos_lael : str | int | list
        显式指定阳性标签样式 如为列表则可指定多个标签
    score_name : str
        数据值名称
    save : bool
        是否保存散点图文件

    Return
    ----------
    str
        生成的散点图文件路径(save == True时)
    '''

    processed_data = _format_data(data, label_col, pos_label, score_name)
    processed_data[label_col].replace(1, 'Positive', inplace=True)
    processed_data[label_col].replace(0, 'Negative', inplace=True)
    processed_data.to_csv('scatter.csv')

    plt.figure(figsize=(10, 10), dpi=300.0)
    sns.scatterplot(data=processed_data, x='PDB',
                    y=score_name, hue='activity', hue_order=['Negative', 'Positive'], alpha=0.75, s=50)
    plt.xlabel('PDB Crystals')
    plt.ylabel(score_name)
    plt.legend(loc='upper right')

    cwd_name = os.path.basename(os.getcwd())
    plt.title('%s %s Scatter' % (cwd_name, score_name))

    if save:
        plt.savefig('%s-Scatter.jpg' % cwd_name)
        return os.path.abspath('%s-Scatter.jpg' % cwd_name)

    plt.show()


def correlate(data: DataFrame, method: str = 'pearson'):
    '''
    计算相关系数矩阵

    Parameters
    ----------
    data: DataFrame
        需要计算的数据
    method : str
        相关系数计算方法
            kendall是定类变量的统计
            pearson是对定距变量的统计
            spearman是对定序变量的统计

    Return
    ----------
    DataFrame
        相关系数矩阵
    '''
    return data.corr(method=method)


def heatmap(data: DataFrame, vmin: float = None, vmax: float = None, save: bool = True):
    '''
    绘制热力图

    Parameters
    ----------
    data : DataFrame
        数据源
    vmin : float
        数据值下限
    vmax : float
        数据值上限
    save : bool
        是否保存热力图文件

    Return
    ----------
    str
        生成的散点图文件路径(save == True时)
    '''
    plt.figure(figsize=(10, 10), dpi=300.0)
    sns.heatmap(data=data, square=True, cmap='RdBu_r', annot=True,
                fmt='.2f', linewidths=0.1, vmin=vmin, vmax=vmax)
    plt.xlabel('PDB Crystals')
    plt.ylabel('PDB Crystals')
    cwd_name = os.path.basename(os.getcwd())
    plt.title('%s Correlation Heatmap' % cwd_name)

    if save:
        plt.savefig('%s-Correlation.jpg' % cwd_name)
        return os.path.abspath('%s-Correlation.jpg' % cwd_name)

    plt.show()


def hyperparam_tuning(model, param_gird: dict, X: DataFrame, y: Series, scoring: str = 'roc_auc', cv: int = 5, n_jobs: int = -1, method: str = 'grid', save: bool = True, model_name: str = None):
    '''
    超参数调优

    Parameters
    ----------
    model : object
        需要调优的模型
    param_grid : dict

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
    logger.info('Best parameters: %s' % best_params)
    logger.info('Best CV score: %s' % best_score)

    if save:
        params_file = 'best_params_%s.json' % (
            model_name if model_name else model.__class__.__name__)
        with open(params_file, 'w') as f:
            json.dump(best_params, f)
        logger.info('Params file %s saved.' % params_file)
    return best_params


def get_splits(X: DataFrame, y: Series, n_repeats: int = 30, n_splits: int = 4, random_state: int = 42):
    '''
    为交叉验证创建训练集和测试集索引

    Parameters
    ----------
    X : DataFrame
        总数据集
    y : Series
        总标签集
    n_repeats : int
        重复次数 n
    n_splits : int
        拆分次数 m
    random_state : int
        随机种子

    Returns
    -------
    list
        索引列表(训练集, 测试集)
    '''
    cv = RepeatedStratifiedKFold(
        n_repeats=n_repeats, n_splits=n_splits, random_state=random_state)
    return [*cv.split(X, y)]


def get_score(model, X_train: DataFrame, X_test: DataFrame, y_train: Series, y_test: Series, score):
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
    score : callable
        评估方法

    Return
    ----------
    tuple(float, float)
        训练集和测试集的ROC-AUC值
    '''

    model.fit(X_train, y_train)
    y_train_predicted = model.predict_proba(X_train)[:, 1]
    y_test_predicted = model.predict_proba(X_test)[:, 1]
    train_score = score(y_train, y_train_predicted)
    test_score = score(y_test, y_test_predicted)
    return train_score, test_score


def get_best_SCP(X: DataFrame, y_true: Series, lower_is_better: bool = True):
    '''
    获取单构象最佳Performance及其ROC-AUC值

    Parameters
    ----------
    X : DataFrame
        训练集
    y_true : Series
        标签
    lower_is_better : bool
        是否为下降性指标

    Return
    ----------
    tuple(str, float)
        best_SCP, best_SCP_score
    '''
    results_dict = {}
    comformations = X.columns
    for comformation in comformations:
        if lower_is_better:
            scp = X[comformation] * -1
        else:
            scp = X[comformation]
        results_dict[comformation] = roc_auc_score(y_true, scp)

    return max(results_dict.items(), key=lambda x: x[1])


def get_SCP_report(splits, X, y):
    '''
    SCP策略交叉测试报告

    Parameters
    ----------
    splits : list   
        索引列表(训练集, 测试集)
    X : DataFrame
        总数据集
    y : Series
        标签

    Return
    ----------
    list
        SCP策略测试AUC值
    '''

    scp_performance = []

    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        best_scp_score = get_best_SCP(X_test, y_test)[1]
        scp_performance.append(best_scp_score)

    _mean = np.mean(scp_performance)
    _max = np.max(scp_performance)
    _min = np.min(scp_performance)
    _std = np.std(scp_performance)
    logger.info('SCP score: %.4f' % (_mean))
    logger.info('SCP Max score: %.4f' % (_max))
    logger.info('SCP Min score: %.4f' % (_min))
    logger.info('SCP std: %.4f' % (_std))

    return scp_performance


def cross_validation(model, splits, X, y, score):
    '''
    交叉验证

    Parameters
    ----------
    model : object
        需要评估的模型
    splits : list
        分割集合
    X : DataFrame
        总数据集
    y : Series  
        标签
    score : callable
        评估方法

    Return
    ----------
    list
        n x m 次模型交叉验证AUC值
    '''
    validation_results = []
    if isinstance(model, MyMLP):

        # 深度神经网络通常不使用交叉验证 但是可以使用
        for train_index, val_index in splits:
            X_train, X_val = X.iloc[train_index], X.iloc[val_index]
            y_train, y_val = y.iloc[train_index], y.iloc[val_index]
            train_dataset = TensorDataset(torch.tensor(
                X_train.values, dtype=torch.float), torch.tensor(y_train.values, dtype=torch.long))
            val_dataset = TensorDataset(torch.tensor(
                X_val.values, dtype=torch.float), torch.tensor(y_val.values, dtype=torch.long))
            train_dataloader = DataLoader(train_dataset, batch_size=128, sampler=WeightedRandomSampler(
                _Evaluator.get_weights(y_train), len(y_train)))
            val_dataloader = DataLoader(val_dataset, batch_size=128, sampler=WeightedRandomSampler(
                _Evaluator.get_weights(y_val), len(y_val)))
            # 新的MLP实例
            params = {
                'input_dim': X_train.shape[1],
                'output_dim': len(set(y_train)),
                'hidden_dim': [32, 16],
                'dropout1': 0.2,
                'dropout2': 0.5,
                'lr': 0.01,
                'weight_decay': 1e-4,
                'batch_size': 128,
                'epochs': 100,
                'device': torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            }
            _model = _Evaluator._get_mlp_classifier(params)
            _model.initialize()
            optimizer = torch.optim.Adam(
                _model.parameters(), lr=params['lr'], weight_decay=params['weight_decay'])
            loss_fn = nn.CrossEntropyLoss()
            _Evaluator.train_model(
                _model, optimizer, loss_fn, train_dataloader, val_dataloader, 100)

            with torch.no_grad():
                _model = _model.to('cpu')
                y_val_pred_proba = _model(torch.tensor(
                    X_val.values, dtype=torch.float)).softmax(dim=1).numpy()[:, 1]
            test_auc = roc_auc_score(y_val.values, y_val_pred_proba)
            validation_results.append(test_auc)
    else:
        for train_index, test_index in splits:
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            train_score, test_score = get_score(
                model, X_train, X_test, y_train, y_test, score)
            validation_results.append(test_score)

    return validation_results


def CV_model_evaluation(models: dict, X: DataFrame, y: Series, n_repeats=30, n_splits=4, random_state=42, plot: bool = False, score_name: str = 'AUC'):
    '''
    (30)x(4)模型评估

    Parameters
    ----------
    models : dict
        所有需要评估的模型 {'model_name': model}
    X : DataFrame
        训练数据集
    y : Series
        训练集标签
    n_repeats : int
        训练次数 n
    n_splits : int
        数据集分割数 k
    random_state : int
        随机种子
    plot : bool
        是否绘制ROC图

    Return
    ----------
    dict
        模型交叉评估AUC结果列
    '''

    splits = get_splits(X, y, n_repeats, n_splits, random_state)
    final_results = {}
    if score_name == 'AUC':
        score = roc_auc_score
    elif score_name == 'F1':
        score = f1_score
    elif score_name == 'Accuracy':
        score = accuracy_score
    elif score_name == 'Precision':
        score = precision_score
    elif score_name == 'Recall':
        score = recall_score
    else:
        raise ValueError(
            'score_name must be one of AUC, F1, Accuracy, Precision, Recall')

    for model_name, model in models.items():

        logger.info('Evaluating model: %s' % model_name)
        _current_result = cross_validation(model, splits, X, y, score)
        final_results[model_name] = _current_result

        _mean = np.mean(_current_result)
        logger.info('%s CV %s mean score: %s' %
                    (model_name, score_name, _mean))

    scp_performance = get_SCP_report(splits, X, y)
    final_results['SCP'] = scp_performance

    if plot:
        plt.figure(figsize=(20, 20))
        result_df = pd.DataFrame(final_results)
        plt.axhline(np.max(scp_performance), color='r', linestyle='--', lw=2)
        plt.axhline(np.mean(scp_performance), color='y', linestyle='--', lw=2)
        plt.axhline(np.min(scp_performance), color='b', linestyle='--', lw=2)
        plt.axhline(0.5, color='g', linestyle='--', lw=2)
        plt.legend(['Max SCP', 'Mean SCP', 'Worst SCP',
                   'Random'], loc='lower right')
        sns.violinplot(data=result_df, palette='Set2',
                       inner='point', split=True)
        plt.show()

    return final_results


class _Evaluator:

    def __init__(self, train_data: DataFrame, test_data: DataFrame):
        '''
        Parameters
        ----------
        train_data : DataFrame
            训练数据集
        test_data : DataFrame
            测试数据集
        '''
        self.train_data = train_data
        self.test_data = test_data
        self.X_train = train_data.drop(['activity'], axis=1)
        self.y_train = train_data['activity']
        self.X_test = test_data.drop(['activity'], axis=1)
        self.y_test = test_data['activity']

        self.log_step_interval = 10

        self.cv_results = None
        self.best_model = None

        self.classifiers = {}
        self.gbt_estimator_hyperparameters = {
            'gpu_id': [0],
            'tree_method': ['gpu_hist'],
            'n_estimators': [200, 250, 300],
            'max_depth': [3, 5, 10, 20],
            'gamma': [0.01, 0.1, 0.5, 1],
            'learning_rate': [0.01, 0.05, 0.1],
            'subsample': [0.3, 0.5, 0.6],
            'alpha': [0.01, 0.1, 0.5, 1],
            'colsample_bytree': [0.3, 0.5, 1],
            'objective': ['binary:logistic'],
            'eval_metric': ['auc'],
            'use_label_encoder': [False]
        }

        self.lr_estimator_hyperparameters = {
            'C': [0.01, 0.1, 1],
            'penalty': ['l1', 'l2']
        }

        self.rf_estimator_hyperparameters = {
            'n_estimators': [200, 250, 300],
            'max_depth': [3, 5, 10, 20],
            'max_features': ['auto', 'sqrt', 'log2'],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 5],
        }

    @staticmethod
    def _preprocess_data(data: DataFrame, test_size=0.25, random_state=42) -> None:
        '''
        划分数据集为特征和标签 并以0填充缺失值
        '''
        X, y = data.drop(['activity'], axis=1), data['activity']
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state)
        train_set = pd.concat([X_train, y_train], axis=1)
        test_set = pd.concat([X_test, y_test], axis=1)
        return train_set, test_set

    @classmethod
    def from_data(cls, data: DataFrame, test_size=0.25, random_state=42) -> '_Evaluator':
        train_set, test_set = cls._preprocess_data(
            data, test_size, random_state)
        return cls(train_set, test_set)

    @property
    def cs_classifiers(self):
        '''
        通识算法分类器
        '''
        return {
            'cs_mean': Average(lower_is_better=True),
            'cs_GeometricMean': Geo_Average(lower_is_better=True),
            'cs_Min': Minimum(lower_is_better=True)
        }

    @property
    def mlp_parameters(self):
        '''
        多层感知机参数 所有参数均在此调节
        '''
        return {
            'batch_size': 128,
            'epochs': 1000,
            'device': torch.device('cuda:0' if torch.cuda.is_available() else 'cpu'),
            'input_dim': self.X_train.shape[1],
            'output_dim': len(set(self.y_train)),
            'hidden_dim': (32, 16),
            'lr': 0.01,
            'weight_decay': 1e-4,
            'dropout1': 0.2,
            'dropout2': 0.5,
            'min_loss': 0.1,
            'min_epochs': 100
        }

    @staticmethod
    def get_weights(y):
        '''
        计算不平衡数据集中各类标签的权重
        '''
        # 样本不平衡时，采用权重采样
        label_to_count = y.value_counts()
        _weights = []
        for label in y:
            _weights.append(1 / label_to_count[label])
        return _weights

    @property
    def train_dataset(self):
        '''
        训练集
        '''
        return TensorDataset(torch.tensor(self.X_train.values, dtype=torch.float), torch.tensor(self.y_train.values, dtype=torch.long))

    @property
    def train_dataloader(self, weighted: bool = True):
        '''
        训练集数据加载器 按照权重采样
        '''
        weights = self.get_weights(self.y_train)
        sampler = WeightedRandomSampler(
            weights, len(weights)) if weighted else None
        batch_size = self.mlp_parameters['batch_size']
        return DataLoader(self.train_dataset, batch_size=batch_size, sampler=sampler)

    @property
    def test_dataset(self):
        '''
        测试集
        '''
        return TensorDataset(torch.tensor(self.X_test.values, dtype=torch.float), torch.tensor(self.y_test.values, dtype=torch.long))

    @property
    def test_dataloader(self, weighted: bool = True):
        '''
        测试集数据加载器 按照权重采样
        '''
        weights = self.get_weights(self.y_test)
        sampler = WeightedRandomSampler(
            weights, len(weights)) if weighted else None
        batch_size = self.mlp_parameters['batch_size']
        return DataLoader(self.test_dataset, batch_size=batch_size, sampler=sampler)

    @staticmethod
    def _get_mlp_classifier(mlp_parameters: dict = None):
        '''
        多层感知机MLP分类器 参数于mlp_parameters调节
        
        Parameters
        ----------
        mlp_parameters : dict
            多层感知机参数 
            {
            'input_dim': int,
            'output_dim': int,
            'hidden_dim': tuple[int, int],
            'lr': float,
            'weight_decay': float,
            'epochs': int,
            'dropout1': float,
            'dropout2': float,
            'batch_size': int,
            'device': torch.device
        }
        '''
        parameters = mlp_parameters
        input_dim = parameters['input_dim']
        output_dim = parameters['output_dim']
        hidden_dim = parameters['hidden_dim']
        dropout1 = parameters['dropout1']
        dropout2 = parameters['dropout2']
        lr = parameters['lr']
        weight_decay = parameters['weight_decay']
        batch_size = parameters['batch_size']
        epochs = parameters['epochs']
        device = parameters['device']

        return MyMLP(input_dim, hidden_dim, output_dim, dropout1, dropout2, lr, weight_decay, batch_size, epochs, device=device)

    def get_mlp_classifier(self):
        '''
        多层感知机MLP分类器 参数于mlp_parameters调节
        '''
        return self._get_mlp_classifier(self.mlp_parameters)

    def get_parameters(self, classifier_name, estimator=None, hyperparameters=None, method='random', force_search=False) -> dict:
        '''
        加载或搜索分类器参数  

        存在best_params_[classifier_name].json时, 加载参数 否则搜索参数

        Parameters
        ----------
        classifier_name : str
            分类器名称
            lr : LogisticRegression
            rf : RandomForestClassifier
            gbt : XGBClassifier

        estimator : sklearn.base.BaseEstimator
            分类器实例
        hyperparameters : dict
            分类器参数搜索范围
        method : str
            搜索方法['random', 'grid']
        force_search : bool
            是否强制重新进行网格搜索

        Returns
        -------
        dict
            分类器参数
        '''
        if estimator is None:
            if classifier_name == 'lr':
                clf_name = 'LogisticRegression'
                estimator = LogisticRegression(max_iter=1000)
                hyperparameters = self.lr_estimator_hyperparameters
            elif classifier_name == 'rf':
                clf_name = 'RandomForestClassifier'
                estimator = RandomForestClassifier()
                hyperparameters = self.rf_estimator_hyperparameters
            elif classifier_name == 'gbt':
                clf_name = 'XGBClassifier'
                estimator = XGBClassifier()
                hyperparameters = self.gbt_estimator_hyperparameters
            else:
                raise ValueError('Estimator is not specified.')
        else:
            if hyperparameters is None:
                raise ValueError('Hyperparameter grid is not specified.')
            clf_name = estimator.__class__.__name__

        # 已存在参数且不强制搜索 则直接读取
        if os.path.exists(f'best_params_{clf_name}.json') and not force_search:
            with open(f'best_params_{clf_name}.json', 'r') as f:
                parameters = json.load(f)
        # 否则进行超参搜索
        else:
            print(f'{classifier_name} parameters are not found, start searching...')
            parameters = hyperparam_tuning(
                estimator, hyperparameters, self.X_train, self.y_train, save=True, model_name=clf_name, method=method)
        return parameters

    def add_classifier(self, classifier_name, parameters=None, classifier_=None) -> None:
        '''
        添加分类器到分类器列表

        Parameters
        ----------
        classifier_name : str
            预设分类器名称
            lr : LogisticRegression
            rf : RandomForestClassifier
            gbt : XGBClassifier
            dummy : DummyClassifier
            nbc : GaussianNB
        parameters : dict | None
            分类器参数 将会被加载到分类器实例中
        classifier_ : sklearn.base.BaseEstimator
            其他非预设的分类器实例
        '''
        if classifier_name == 'lr':
            classifier = LogisticRegression(max_iter=1000)
            clf_name = 'LogisticRegression'
        elif classifier_name == 'rf':
            classifier = RandomForestClassifier()
            clf_name = 'RandomForestClassifier'
        elif classifier_name == 'gbt':
            classifier = XGBClassifier(
                eval_metric='auc', use_label_encoder=False)
            clf_name = 'XGBClassifier'
        elif classifier_name == 'dummy':
            classifier = DummyClassifier(
                strategy='stratified', random_state=42)
            clf_name = 'DummyClassifier'
        elif classifier_name == 'nbc':
            classifier = GaussianNB()
            clf_name = 'NaiveBayesClassifier'
        elif classifier_name == 'cs':
            self.classifiers.update(self.cs_classifiers)
            return None
        elif classifier_name == 'mlp':
            classifier = self.get_mlp_classifier()
            clf_name = 'MLPClassifier'
        else:
            if classifier_ is None:
                raise ValueError('Classifier is required')
            else:
                classifier = classifier_
                clf_name = classifier_name.upper()

        if parameters is not None:
            try:
                classifier.set_params(**parameters)
            except AttributeError:
                pass
        self.classifiers[clf_name] = classifier

    def cross_validation(self, plot=True, **kwargs):
        '''
        将分类器列表中的分类器执行30x4交叉验证

        Parameters
        ----------
        plot : bool
            是否绘制交叉验证结果图表
        kwargs : dict
            其他传递给CV_model_evaluation的参数
        '''
        results = CV_model_evaluation(
            self.classifiers, self.X_train, self.y_train, plot=plot, **kwargs)
        report_results = {}
        for clf_name, clf in self.classifiers.items():
            report_results[clf_name] = np.mean(results[clf_name])

        report_results['Best_SCP'] = np.max(results['SCP'])
        report_results['Mean_SCP'] = np.mean(results['SCP'])
        report_results['Worst_SCP'] = np.min(results['SCP'])
        report_results['SCP_std'] = np.std(results['SCP'])

        return report_results

    @staticmethod
    def init_weights(model):
        '''
        初始化模型权重
        '''
        for param in model.parameters():
            nn.init.normal_(param, mean=0, std=0.01)

    @staticmethod
    def train_model(model, optimizer, loss_function, train_dataloader, test_dataloader, epochs=500, device=None, min_epochs=100, min_loss=0.1, plot=False):
        '''
        训练模型

        保存最佳模型(Test loss低于min_loss, 训练epoch数大于min_epochs)参数

        Parameters
        ----------
        model : torch.nn.Module
            模型
        optimizer : torch.optim.Optimizer
            优化器
        loss_function : torch.nn.modules.loss
            损失函数
        train_dataloader : torch.utils.data.DataLoader
            训练集数据加载器
        test_dataloader : torch.utils.data.DataLoader
            测试集数据加载器
        epochs : int
            训练次数
        device : torch.device
            训练设备 CPU/GPU
        min_epochs : int
            最小训练次数
        min_loss : float
            最小损失
        plot : bool
            是否在训练过程中绘制状态(Loss/Acc/Auc)曲线

        Returns
        -------
        tuple[history, best_checkpoint]

        history : Hiddenlayer.History
            训练状态记录器 包含keys: train_loss, train_acc, train_auc, test_loss, test_acc, test_auc
            通过history[key].data获取数据
        best_checkpoint : torch.nn.Module
            最佳模型检查点参数字典

            {
                'model': model.state_dict(), 
                'optimizer': optimizer.state_dict(),
                'epoch': epoch + 1
            }
        '''
        if device is None:
            device = torch.device(
                'cuda:0' if torch.cuda.is_available() else 'cpu')
        history = hl.History()
        best_checkpoint = None
        model.to(device)
        for epoch in range(epochs):

            train_loss = []
            train_acc = []
            train_auc = []
            model.train()
            for i, (train_data_batch, train_label_batch) in enumerate(train_dataloader):
                train_data_batch_gpu = train_data_batch.float().to(device)
                train_label_batch_gpu = train_label_batch.to(device)
                train_outputs = model(train_data_batch_gpu)

                loss = loss_function(train_outputs, train_label_batch_gpu)
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                train_loss.append(loss.detach().cpu())
                train_outputs = train_outputs.detach().cpu()
                train_pred = train_outputs.argmax(dim=1)
                train_proba = train_outputs.softmax(dim=1)[:, 1]

                train_acc.append(accuracy_score(
                    train_label_batch.numpy(), train_pred.numpy()))
                train_auc.append(roc_auc_score(
                    train_label_batch.numpy(), train_proba.numpy()))

            train_loss = np.mean(train_loss)
            train_acc = np.mean(train_acc)
            train_auc = np.mean(train_auc)

            model.eval()
            with torch.no_grad():
                test_acc = []
                test_auc = []
                test_loss = []

                for i, (test_data_batch, test_label_batch) in enumerate(test_dataloader):
                    test_data_batch_gpu = test_data_batch.float().to(device)
                    test_outputs = model(test_data_batch_gpu)
                    test_outputs = test_outputs.detach().cpu()

                    batch_test_loss = loss_function(
                        test_outputs, test_label_batch)
                    test_pred = test_outputs.argmax(dim=1)
                    test_proba = test_outputs.softmax(dim=1)[:, 1]

                    test_acc.append(accuracy_score(
                        test_label_batch.numpy(), test_pred.numpy()))
                    test_auc.append(roc_auc_score(
                        test_label_batch.numpy(), test_proba.numpy()))
                    test_loss.append(batch_test_loss.item())

                test_acc = np.mean(test_acc)
                test_auc = np.mean(test_auc)
                test_loss = np.mean(test_loss)

                if test_loss < min_loss and epoch >= min_epochs:
                    min_loss = test_loss
                    best_checkpoint = {'model': deepcopy(model.state_dict()), 'optimizer': deepcopy(
                        optimizer.state_dict()), 'epoch': epoch + 1, 'test_loss': test_loss}

            history.log(
                epoch + 1,
                train_loss=loss.item(),
                test_loss=test_loss.item(),
                train_acc=train_acc,
                train_auc=train_auc,
                test_acc=test_acc,
                test_auc=test_auc
            )

            if plot:
                canvas = hl.Canvas()
                with canvas:
                    canvas.draw_plot(
                        [history['train_loss'], history['test_loss']], ylabel='Loss')
                    canvas.draw_plot(
                        [history['train_acc'], history['test_acc']], ylabel='Accuracy')
                    canvas.draw_plot(
                        [history['train_auc'], history['test_auc']], ylabel='AUC')
            else:
                print(f'Epoch {epoch + 1}/{epochs} | Train Loss: {train_loss:.4f} | Test Loss: {test_loss:.4f} | Train Acc: {train_acc:.4f} | Test Acc: {test_acc:.4f} | Train AUC: {train_auc:.4f} | Test AUC: {test_auc:.4f}')

        return history, best_checkpoint

    def mlp_train(self, plot=True):
        '''
        基于MLP模型训练

        训练完成后 history 和 best_checkpoint 属性会被赋值于evaluator实例history和best_checkpoint属性
        最终模型(完成所有epoch训练)检查点保存到evaluator实例final_model属性

        Parameters
        ----------
        plot : bool
            是否在训练过程中绘制状态(Loss/Acc/Auc)曲线

        Return
        ------
        tuple(nn.Module, nn.Module)
            最佳模型, 最终模型
        '''
        X_train, X_test, y_train, y_test = self.X_train, self.X_test, self.y_train, self.y_test
        device = self.mlp_parameters['device']
        epochs = self.mlp_parameters['epochs']
        weight_decay = self.mlp_parameters['weight_decay']
        lr = self.mlp_parameters['lr']
        min_epochs = self.mlp_parameters['min_epochs']
        min_loss = self.mlp_parameters['min_loss']

        # 实例化一个模型 该模型已在GPU上 所有向量需要转换到GPU
        model = self.get_mlp_classifier()
        # model.to(device)
        model.initialize()

        train_dataloader = self.test_dataloader
        test_dataloader = self.test_dataloader

        loss_function = nn.CrossEntropyLoss()
        # NOTE: 务必在初始化模型完成后再定义优化器
        optimizer = torch.optim.Adam(
            model.parameters(), lr=lr, weight_decay=weight_decay)

        with torch.no_grad():
            current_auc = roc_auc_score(y_train, np.argmax(model(torch.tensor(
                X_train.values, dtype=torch.float).to(device)).detach().cpu().numpy(), axis=1))
        print('Current AUC:', current_auc)
        print('Start Training...')
        sleep(1)

        history, best_checkpoint = self.train_model(
            model, optimizer, loss_function, train_dataloader, test_dataloader, epochs, device, min_epochs, min_loss, plot)
        self.history = history
        self.best_checkpoint = best_checkpoint
        self.final_checkpoint = {'model': model.state_dict(
        ), 'optimizer': optimizer.state_dict(), 'epoch': epochs}

        self.best_model = self.get_mlp_classifier()
        self.best_model.load_state_dict(best_checkpoint['model'])
        self.final_model = model

        return self.best_model, self.final_model

    def print_classifier_info(self):
        '''
        打印分类器参数信息
        '''
        formatter = '{0:<25}{1:<30}{2:<40}'
        print('='*100)
        print(formatter.format('Classifier', 'Parameters', 'Values'))
        print('-'*100)
        for clf_name, clf in self.classifiers.items():
            parameters = [(k, v) for k, v in clf.get_params().items()]
            print(formatter.format(clf_name, str(
                parameters[0][0]), str(parameters[0][1])))
            for k, v in parameters[1:]:
                print(formatter.format('', str(k), str(v)))
            print('-'*100)
        print('='*100)

    def print_cv_results(self):
        '''
        打印交叉验证结果
        '''
        formatter = '{0:<25}{1:<30}'
        print(formatter.format('Classifier', 'AUC'))
        print(formatter.format('-'*25, '-'*30))
        for k, v in self.cv_results.items():
            print(formatter.format(k, v))

    def draw_after_training(self):
        '''
        训练完成后 绘制训练过程状态曲线
        '''
        history = self.history
        canvas = hl.Canvas()
        with canvas:
            canvas.draw_plot(
                [history['train_loss'], history['test_loss']], ylabel='Loss')
            canvas.draw_plot(
                [history['train_acc'], history['test_acc']], ylabel='Accuracy')
            canvas.draw_plot(
                [history['train_auc'], history['test_auc']], ylabel='AUC')

    @staticmethod
    def evaluate(y_test, y_pred, y_proba, print_=True):
        '''
        评估分类器
        Parameters
        ----------
        y_test : numpy.ndarray
            测试集标签
        y_pred : numpy.ndarray
            预测集标签
        y_proba : numpy.ndarray
            预测集概率
        print_ : bool
            是否打印评估结果
        '''
        acc = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        auc = roc_auc_score(y_test, y_proba)
        fpr, tpr, _ = roc_curve(y_test, y_proba)

        if print_:
            print(classification_report(y_test, y_pred))
            print('Confusion Matrix\n', confusion_matrix(y_test, y_pred))
            print('Accuracy:', acc)
            print('Precision:', precision)
            print('Recall:', recall)
            print('F1 score:', f1)
            print('ROC AUC:', auc)

        return {
            'accuracy': acc,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'auc': auc,
            'roc': [fpr, tpr]
        }

    def report(self, plot=True):
        '''
        打印所有模型的评估结果
        Parameters
        ----------
        plot : bool
            是否绘制ROC结果图
        '''
        roc_curves_ = {}
        auc_ = {}
        for clf_name, clf in self.classifiers.items():
            print('-'*50)
            print('Model:', clf_name, '\n')
            if not isinstance(clf, MyMLP):
                clf.fit(self.X_train, self.y_train)
            y_pred = clf.predict(self.X_test)
            y_proba = clf.predict_proba(self.X_test)[:, 1]
            if not isinstance(clf, Consensus):
                _test_results = self.evaluate(
                    self.y_test, y_pred, y_proba, print_=True)
            else:
                _test_results['auc'] = roc_auc_score(self.y_test, y_proba)
                _test_results['roc'] = roc_curve(self.y_test, y_proba)
            roc_curves_[clf_name] = _test_results['roc']
            auc_[clf_name] = _test_results['auc']
            print('-'*50)

        if plot:
            plt.figure(figsize=(20, 20), dpi=300)
            plt.title('ROC Curve')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.0])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            for clf_name, roc_curve_ in roc_curves_.items():
                plt.plot(roc_curve_[0], roc_curve_[
                         1], label=f'{clf_name} (AUC={auc_[clf_name]:.3f})')
            plt.plot([0, 1], [0, 1], linestyle='--',
                     color='r', label='Random Guess (AUC=0.5)')
            plt.legend(loc='lower right')

    def _save_model(self, model_state_dict, save_path):
        '''
        保存模型状态参数

        Parameters
        ----------
        model_state_dict : dict
            模型检查点状态字典
        save_path : str
            保存文件路径
        '''
        model_state_dict.update(
            {
                'structure':
                {
                    'input_dim': self.mlp_parameters['input_dim'],
                    'hidden_dim': self.mlp_parameters['hidden_dim'],
                    'output_dim': self.mlp_parameters['output_dim'],
                    'dropout1': self.mlp_parameters['dropout1'],
                    'dropout2': self.mlp_parameters['dropout2'],
                }
            }
        )
        torch.save(model_state_dict, save_path)

    def save_model(self, dir_path='./models'):
        '''
        保存最佳模型和最终模型 并保存训练记录数据

        Parameters
        ----------
        dir_path : str
            保存目录路径
        '''
        self._save_model(self.best_checkpoint,
                         os.path.join(dir_path, 'best_model.pth'))
        self._save_model(self.final_checkpoint,
                         os.path.join(dir_path, 'final_model.pth'))
        with open(os.path.join(dir_path, 'history.pkl'), 'wb') as f:
            pickle.dump(self.history, f)

    def load_model(self, path):
        '''
        加载模型参数

        Parameters
        ----------
        path : str
            模型参数文件路径

        Returns
        -------
        tuple[nn.Module, torch.optim.Optimizer]
            加载状态后的MLP模型与优化器
        '''
        base_model = self.get_mlp_classifier()
        checkpoint = torch.load(path)
        base_model.load_state_dict(checkpoint['model'])
        base_optimizer = torch.optim.Adam(params=base_model.parameters())
        base_optimizer.load_state_dict(checkpoint['optimizer'])

        return base_model, base_optimizer

    def load_history(self, path):
        '''
        加载训练记录数据 赋值于evaluator实例history属性
        '''
        with open(path, 'rb') as f:
            self.history = pickle.load(f)

    def evaluate_workflow(self):
        '''
        评估数据集的标准工作流
            具体执行以下步骤:

                1. 添加预设的分类器(LR、XGBoost、RF、NBC、Dummy)以及通识算法(Mean、GeometricMean、Minimum)到待评估列表   
                2. 训练MLP分类器   
                3. 对待评估列表中的分类器执行30x4交叉验证   
                4. 报告评估结果 并保存于 cv_results.json 文件中   

        '''
        # lr_params = self.get_parameters('lr', method='grid', force_search=force_search)
        # gbt_params = self.get_parameters('gbt', method='random', force_search=force_search)
        # rf_params = self.get_parameters('rf', method='random', force_search=force_search)
        # self.add_classifier('lr', parameters=lr_params)
        # self.add_classifier('gbt', parameters=gbt_params)
        # self.add_classifier('rf', parameters=rf_params)
        self.add_classifier('lr')
        self.add_classifier('gbt')
        self.add_classifier('rf')
        self.add_classifier('nbc')
        self.add_classifier('dummy')
        self.add_classifier('cs')
        self.add_classifier('mlp')
        self.print_classifier_info()

        # 不需要对MLP进行CV时 可以取消注释下面的代码
        # self.classifiers.pop('MLPClassifier')

        print('='*100)
        print('Start Training MLP...')
        sleep(3)
        self.mlp_train(plot=False)
        self.save_model()
        print('Training finished. Use draw_after_training() to draw the results.')
        print('='*100)

        print('='*100)
        print('Start Cross-Validation...')
        sleep(1)
        self.cv_results = self.cross_validation()
        self.print_cv_results()
        with open('cv_results.json', 'w') as f:
            json.dump(self.cv_results, f)
        print('='*100)

        self.add_classifier('BEST MLP', classifier_=self.best_model)
        self.add_classifier('FINAL MLP', classifier_=self.final_model)
        self.report()
