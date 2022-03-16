import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import logging
import json

from pandas import DataFrame, Series
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.model_selection import RepeatedStratifiedKFold

logger = logging.getLogger(__name__)

def read_matrix(file_path: str):
    '''
    读取数据矩阵

    Parameter
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
    Parameter
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

def _standrad_label(data: DataFrame, label_col: str, positive: ...='origin'):
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

def split_data(data: DataFrame, label_col:str=None, preprocess:bool=True, positive:str='origin'):
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
    Parameter
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


def get_roc(X: DataFrame, y_ture: Series , save: bool = False, lower_is_better: bool = True):
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
                    y=score_name, hue='activity',hue_order=['Negative', 'Positive'], alpha=0.75, s=50)
    plt.xlabel('PDB Crystals')
    plt.ylabel(score_name)
    plt.legend(loc='upper right')

    cwd_name = os.path.basename(os.getcwd())
    plt.title('%s %s Scatter' %(cwd_name ,score_name))

    if save:
        plt.savefig('%s-Scatter.jpg' % cwd_name)
        return os.path.abspath('%s-Scatter.jpg' % cwd_name)

    plt.show()

def correlate(data: DataFrame, method:str='pearson'):
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

def heatmap(data:DataFrame, vmin:float=None, vmax:float=None, save:bool=True):
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
    plt.figure(figsize=(10, 10),dpi=300.0)
    sns.heatmap(data=data, square=True, cmap='RdBu_r', annot=True, fmt='.2f', linewidths=0.1, vmin=vmin, vmax=vmax)
    plt.xlabel('PDB Crystals')
    plt.ylabel('PDB Crystals')
    cwd_name = os.path.basename(os.getcwd())
    plt.title('%s Correlation Heatmap' % cwd_name)

    if save:
        plt.savefig('%s-Correlation.jpg' % cwd_name)
        return os.path.abspath('%s-Correlation.jpg' % cwd_name)

    plt.show()

def hyperparam_tuning(model, param_gird:dict, X:DataFrame, y:Series, scoring:str='roc_auc', cv:int=5, n_jobs:int=-1, method:str='grid', save:bool=True, model_name:str=None):
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
            grid: 网格搜索
            random: 随机搜索
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
        cv_search = GridSearchCV(estimator=model, param_grid=param_gird, scoring=scoring, cv=cv, n_jobs=n_jobs)
    elif method == 'random':
        cv_search = RandomizedSearchCV(estimator=model, param_distributions=param_gird, scoring=scoring, cv=cv, n_jobs=n_jobs)
    else:
        raise ValueError('Invalid method: %s' % method)
    
    cv_search.fit(X, y)
    best_params = cv_search.best_params_
    best_score = cv_search.best_score_
    logger.info('Best parameters: %s' % best_params)
    logger.info('Best CV score: %s' % best_score)

    if save:
        params_file = 'best_params_%s.json' % (model_name if model_name else model.__class__.__name__ )
        with open(params_file, 'w') as f:
            json.dump(best_params, f)
        logger.info('Params file %s saved.' % params_file)
    return best_params

def get_splits(X:DataFrame, y:Series, n_repeats:int=30, n_splits:int=4, random_state:int=42):
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
    cv = RepeatedStratifiedKFold(n_repeats=n_repeats, n_splits=n_splits, random_state=random_state)
    return [*cv.split(X, y)]

def cross_validation(model, splits, X, y):
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
    
    Return
    ----------
    list
        n x m 次模型交叉验证AUC值
    '''
    validation_results = []

    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        train_score, test_score = get_auc_score(model, X_train, X_test, y_train, y_test)
        validation_results.append(test_score)
    
    return validation_results

def get_auc_score(model, X_train:DataFrame, X_test:DataFrame, y_train:Series, y_test:Series):
    '''
    模型ROC-AUC评估

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
    
    Return
    ----------
    tuple(float, float)
        训练集和测试集的ROC-AUC值
    '''

    model.fit(X_train, y_train)
    y_train_predicted = model.predict_proba(X_train)[:,1]
    y_test_predicted = model.predict_proba(X_test)[:,1]
    train_score = roc_auc_score(y_train, y_train_predicted)
    test_score = roc_auc_score(y_test, y_test_predicted)
    return train_score, test_score

def get_best_SCP(X:DataFrame, y_true:Series, lower_is_better:bool=True):
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
    _std = np.std(scp_performance)
    logger.info('SCP score: %.4f' % (_mean))
    logger.info('SCP Max score: %.4f' % (_max))
    logger.info('SCP std: %.4f' % (_std))

    return scp_performance
    
def CV_model_evaluation(models:dict, X:DataFrame, y:Series, n_repeats=30, n_splits=4, random_state=42, plot:bool=False):
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

    for model_name, model in models.items():

        logger.info('Evaluating model: %s' % model_name)
        _current_result = []

        if model_name.startswith('ml_'):
            _current_result = cross_validation(model, splits, X, y)
        elif model_name.startswith('cs_'):
            for train_index, test_index in splits:
                X_train, X_test = X.iloc[train_index], X.iloc[test_index]
                y_train, y_test = y.iloc[train_index], y.iloc[test_index]
                y_predicted = model(X_test.replace(0, np.nan)) * -1
                test_score = roc_auc_score(y_test, y_predicted)
                _current_result.append(test_score)

        final_results[model_name] = _current_result
        _mean = np.mean(_current_result)
        logger.info('%s CV score: %s' % (model_name, _mean))

    return final_results
