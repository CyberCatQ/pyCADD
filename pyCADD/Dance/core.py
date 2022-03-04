import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import logging
import json

from pandas import DataFrame, Series
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV
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


def read_docking_data(raw_data: DataFrame, label_col: str):
    '''
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
    '''

    return raw_data.drop(label_col, axis=1)


def merge(data_list: list):
    '''
    合并Series
    Parameter
    ----------
    data_list : list
        包含需要合并的Series的列表
    '''
    return pd.concat(data_list, axis=1)


def _standrad_label(data: DataFrame, label_col: str, pos_label: str):
    '''
    标准化标签为二进制
    '''
    if isinstance(pos_label, str) or isinstance(pos_label, int):
        pos_label = [pos_label]

    for _label in pos_label:
        data[label_col].replace(_label, value=1, inplace=True)

    data[label_col] = pd.to_numeric(
        data[label_col], errors='coerce').fillna(0).astype('int32')

    # 标签列移至末尾
    _label_col = data.pop(label_col)
    data.insert(data.shape[1], column=label_col, value=_label_col)

    return data


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


def get_auc(data: DataFrame, label_col: str, pos_label, save: bool = False, ascending: bool = False):
    '''
    ROC曲线下面积
    依据label列作为标签
    自动计算DataFrame的所有列的ROC曲线下面积AUC值

    Parameters
    ----------
    data : DataFrame
        待计算数据
    label : str
        阳性标签列名
    pos_lael : str | int | list
        显式指定阳性标签样式 如为列表则可指定多个标签
    save : bool
        是否存储ROC曲线图片
    ascending : bool
        是否以升序方式排序(如score为负数 则应为True)

    Return
    ---------
    Series
        曲线下面积AUC数据
    '''

    # 标签列二进制化
    data = _standrad_label(data, label_col, pos_label)

    auc_dict = {}

    plt.figure(figsize=(10, 10), dpi=300.0)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('%s ROC curve' % os.path.basename(os.getcwd()))

    for index in data.columns[:-1]:
        if ascending:
            data[index] = - data[index]
        _column = data[[label_col, index]].dropna(how='any')
        fpr, tpr, thersholds = roc_curve(_column[label_col], _column[index])
        auc = roc_auc_score(_column[label_col], _column[index])
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

def split_data(data, label_col:str=None):
    '''
    拆分数据与标签
    
    Parameters
    ----------
    data : DataFrame
        待拆分数据

    Return  
    ----------
    DataFrame, Series
        拆分后的数据和标签
    '''
    if label_col is None:
        return data.iloc[:,:-1], data.iloc[:,-1]
    else:
        return data.drop(label_col, axis=1), data[label_col]

def hyperparam_tuning(model, params, X, y, scoring='roc_auc', cv=5, n_jobs=-1, method='grid', save=False, model_name=None):
    '''
    超参数调优

    Parameters
    ----------
    model : object
        需要调优的模型
    params : dict
        模型参数
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
        cv_search = GridSearchCV(estimator=model, param_grid=params, scoring=scoring, cv=cv, n_jobs=n_jobs)
    elif method == 'random':
        cv_search = RandomizedSearchCV(estimator=model, param_distributions=params, scoring=scoring, cv=cv, n_jobs=n_jobs)
    else:
        raise ValueError('Invalid method: %s' % method)
    
    cv_search.fit(X, y)
    best_params = cv_search.best_params_
    best_score = cv_search.best_score_
    logger.info('Best parameters: %s' % best_params)
    logger.info('Best CV score: %s' % best_score)

    if save:
        with open('best_params_%s.json' % (model_name if model_name else model.__class__.__name__), 'w') as f:
            json.dump(best_params, f)
    return best_params


def _GS_result_report(model, params, X_train, X_test, y_train, y_test):
    '''
    网格搜索结果报告(ROC-AUC评估)

    Parameters
    ----------
    model : object
        需要评估的模型
    params : dict
        调优后的模型参数
    X_train : DataFrame
        训练集
    X_test : DataFrame
        测试集
    y_train : Series
        训练集标签
    y_test : Series
        测试集标签

    '''
    if params:
        model.set_params(**params)
    model.fit(X_train, y_train)
    y_train_predicted = model.predict_proba(X_train)[:,1]
    y_test_predicted = model.predict_proba(X_test)[:,1]
    train_score = roc_auc_score(y_train, y_train_predicted)
    test_score = roc_auc_score(y_test, y_test_predicted)
    return train_score, test_score

def get_best_SCP(X, y_true):
    '''
    获取单构象最佳Performance ROC-AUC

    Parameters
    ----------
    X : DataFrame
        训练集
    y_true : Series
        标签

    Return
    ----------
    float
        best_SCP
    '''
    results_dict = {}
    comformations = X.columns
    for comformation in comformations:
        scp = X[comformation] * -1
        results_dict[comformation] = roc_auc_score(y_true, scp)
    
    best_scp = max(results_dict.values())
    return best_scp   

def get_SCP_report(X, y, n_repeats:int=30, n_splits:int=4, random_state=42):
    '''
    最佳SCP策略报告

    '''
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    splits = [*cv.split(X, y)]
    scp_performance = []

    for train_index, test_index in splits:
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]
        best_scp = get_best_SCP(X_test, y_test)
        scp_performance.append(best_scp)
        
    logger.info('SCP score: %.4f' % (np.mean(scp_performance)))
    logger.info('SCP Max score: %.4f' % (np.max(scp_performance)))
    logger.info('SCP std: %.4f' % (np.std(scp_performance)))
    

def CV_model_evaluation(models:dict, X:pd.DataFrame, y:pd.Series, n_repeats=30, n_splits=4, random_state=42, scoring='roc_auc'):
    '''
    (30)x(4)模型评估

    Parameters
    ----------
    models : dict
        所有需要评估的模型 {'model_name': model}
    X : DataFrame
        总数据集
    y : Series
        总标签
    n_repeats : int
        训练次数
    n_splits : int
        数据集分割数
    random_state : int
        随机种子
    scoring : str
        评估标准
    
    Return
    ----------
    dict
        模型评估结果
    '''
    cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    splits = [*cv.split(X, y)]
    final_results = {}

    for model_name, model in models.items():

        logger.info('Evaluating model: %s' % model_name)
        _current_result = []

        # n_repeats x n_splits 交叉验证 
        for train_index, test_index in splits:
            X_train, X_test = X.iloc[train_index], X.iloc[test_index]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            if model_name.startswith('ml_'):
                train_score, test_score = _GS_result_report(model, None, X_train, X_test, y_train, y_test)
            elif model_name.startswith('cs_'):
                y_predicted = model(X_test.replace(0, np.nan)) * -1
                test_score = roc_auc_score(y_test, y_predicted)
            _current_result.append(test_score)

        final_results[model_name] = np.mean(_current_result)
        logger.info('%s CV score: %s' % (model_name, final_results[model_name]))
    logger.info('Final results: %s' % final_results)
    return final_results

def _get_model_params(model_name, params_file):
    '''
    获取模型参数

    Parameters
    ----------
    model_name : str
        模型名称
    params_file : str
        模型参数文件

    Return
    ----------
    dict
        模型参数
    '''
    if model_name.startswith('ml_'):
        params = pd.read_csv(params_file).to_dict()
    elif model_name.startswith('cs_'):
        params = None
    else:
        raise ValueError('Invalid model name: %s' % model_name)
    return params

