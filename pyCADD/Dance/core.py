import os

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pandas import DataFrame, Series
from sklearn.metrics import roc_auc_score, roc_curve


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
