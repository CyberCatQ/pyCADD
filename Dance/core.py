import os
import pandas as pd
from pandas import DataFrame, Series
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt

def read_matrix(file_path:str):
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
    
def read_docking_data(raw_data:DataFrame, label_col:str):
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

def merge(data_list:list):
    '''
    合并Series
    Parameter
    ----------
    data_list : list
        包含需要合并的Series的列表
    '''
    return pd.concat(data_list, axis=1)

def get_auc(data:DataFrame, label_col:str, pos_label, save:bool=False, ascending:bool=False):
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
    # 标准化标签列为二进制
    if isinstance(pos_label, str) or isinstance(pos_label, int):
        pos_label = [pos_label]

    for _label in pos_label:
        data[label_col].replace(_label, value=1, inplace=True)

    data[label_col] = pd.to_numeric(data[label_col], errors='coerce').fillna(0).astype('int32')
    
    # 标签列移至末尾
    _label_col = data.pop(label_col)
    data.insert(data.shape[1], column=label_col, value=_label_col)
    auc_dict = {}

    plt.figure(figsize=(10,10), dpi=300.0)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC curve')

    for index in data.columns[:-1]:
        if ascending:
            data[index] = - data[index]
        fpr, tpr, thersholds = roc_curve(data[label_col], data[index])
        auc = roc_auc_score(data[label_col], data[index])
        auc_dict[index] = auc
        plt.plot(fpr, tpr, label='%s (area = %.2f)' % (index, auc))
    plt.plot([0 ,1], [0, 1], linestyle='--')
    plt.legend(loc='lower right')
    
    cwd_name = os.path.basename(os.getcwd())
    if save:
        plt.savefig('%s-ROC.jpg' % cwd_name)

    plt.show()

    return Series(auc_dict, name='AUC')
        




