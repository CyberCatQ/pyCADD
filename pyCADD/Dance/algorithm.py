# 数据融合规则算法
import pandas as pd
import numpy as np
from pandas import DataFrame, Series

# 包装为类
class Consensus:
    '''
    共识算法基类
    '''
    def __init__(self, X:DataFrame=None, y:Series=None, lower_is_better=False):
        '''
        Parameters
        ----------
        data : DataFrame
            待计算数据
        lower_is_better : bool
            是否为下降性指标
        '''
        self.X = X
        self.y = y
        self.result = None
        self.lower_is_better = lower_is_better
    
    def fit(self):
        '''
        共识算法计算
        '''
        raise NotImplementedError
    
    def predict(self, X:DataFrame=None, y:Series=None):
        if X is not None:
            self.fit(X, y)
        return self.result if self.result is not None else None

    def predict_proba(self, X:DataFrame=None, y:Series=None):
        '''
        按照计算结果进行预测
        '''
        if X is not None:
            self.fit(X, y)
        return pd.DataFrame(self.result).reset_index().to_numpy() if self.result is not None else None

    def score(self, method, labels=None):
        '''
        评分
        Parameters
        ----------
        method : callable
            评分函数
        labels : list | ndarray | Series
            实际标签
        '''
        if labels is not None:
            self.y = labels
        if self.y is None:
            raise ValueError('labels is None')
        return method(self.y, self.result) if self.result is not None else None

class Average(Consensus):
    '''
    算数平均值模型
    '''
    def __init__(self, X=None, y=None, lower_is_better=False):
        super().__init__(X, y, lower_is_better)
        self.result = None
        
    def fit(self, X:DataFrame=None, y:Series=None):
        '''
        计算算数平均
        '''
        self.X = X.replace(0, np.nan) if X is not None else self.X
        self.y = y if y is not None else self.y
        self.result = self.X.mean(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result

class Geo_Average(Consensus):
    '''
    几何平均值模型
    '''
    def __init__(self, X=None, y=None, lower_is_better=False):
        super().__init__(X, y, lower_is_better)
        self.result = None
        
    def fit(self, X:DataFrame=None, y:Series=None):
        '''
        计算几何平均
        '''
        self.X = X.replace(0, np.nan) if X is not None else self.X
        self.y = y if y is not None else self.y
        self.result = Series(-1 * pow(self.X.prod(axis=1).abs(), 1/self.X.notna().sum(axis=1)), name='GEO')
        if self.lower_is_better:
            self.result = -1 * self.result

class Minimum(Consensus):
    '''
    最小值模型
    '''
    def __init__(self, X=None, y=None, lower_is_better=False):
        super().__init__(X, y, lower_is_better)
        self.result = None
        
    def fit(self, X:DataFrame=None, y:Series=None):
        '''
        计算几何平均
        '''
        self.X = X.replace(0, np.nan) if X is not None else self.X
        self.y = y if y is not None else self.y
        self.result = self.X.min(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result

class Maximum(Consensus):
    '''
    最大值模型
    '''
    def __init__(self, X=None, y=None, lower_is_better=False):
        super().__init__(X, y, lower_is_better)
        self.result = None
        
    def fit(self, X:DataFrame=None, y:Series=None):
        '''
        计算几何平均
        '''
        self.X = X.replace(0, np.nan) if X is not None else self.X
        self.y = y if y is not None else self.y
        self.result = self.X.max(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result

# 函数
def average(data: DataFrame, method: str = 'ave'):
    '''
    平均值算法 : 计算并生成DataFrame对接数据的平均值列  
    GEO 几何平均
    AVE 算数平均
    Parameters
    ----------
    data : DataFrame
        待计算数据
    method : str
        平均值计算方法 ( 算术平均值ave | 几何平均值geo )

    Return
    ----------
    Series
        平均值结果数据列
    '''
    if method.upper() == 'AVE':
        return Series(data.mean(axis=1), name='AVE')
    elif method.upper() == 'GEO':
        # 首先消除符号再计算几何平均
        # 最后添加符号
        return Series(- pow(data.prod(axis=1).abs(), 1/data.notna().sum(axis=1)), name='GEO')
    else:
        raise RuntimeError('Invalid average value method.')

def geo_average(data: DataFrame):
    '''
    几何平均 : 计算DataFrame数据的几何平均值
    Parameters
    ----------
    data : DataFrame
        待计算数据

    Return
    ----------
    Series
        几何平均结果数据列
    '''
    return Series(- pow(data.prod(axis=1).abs(), 1/data.notna().sum(axis=1)), name='GEO')

def minimum(data: DataFrame):
    '''
    最小值(最优评分) : 提取DataFrame数据列中的最小值(如最佳对接分数)

    Parameters
    ----------
    data : DataFrame
        待计算数据

    Return
    ----------
    Series
        最小值结果数据列
    '''
    return Series(data.min(axis=1), name='MIN')


def maximum(data: DataFrame):
    '''
    最大值 : 提取DataFrame数据列中的最大值

    Parameters
    ----------
    data : DataFrame
        待计算数据

    Return
    ----------
    Series
        最大值结果数据列
    '''
    return Series(data.max(axis=1), name='MAX')


def std(data: DataFrame, axis: int = 1):
    '''
    计算标准偏差

    Parameters
    ----------
    data : DataFrame
        待计算数据
    axis : int
        计算坐标轴 row: 0, column: 1

    Return
    ----------
    Series
        标准差结果数据列
    '''
    return Series(data.std(axis=axis), name='std')


def z_score(data: DataFrame, ratio: tuple = (0.7, 0.3)):
    '''
    标准分数(Z-score) : 对DataFrame数据进行标准分数变换
    Z-score-combined = 0.7 * Z-score-receptor + 0.3 * Z-score-ligand

    Parameters
    ----------
    data : DataFrame
        待计算数据
    ratio : tuple
        z_score组合权重比例(A : B)
            A: 受体Z_score分数权重 default: 0.7
            B: 配体Z_score分数权重 default: 0.3

    Return
    ----------
    tuple(DataFrame, DataFrame, DataFrame)
        z_score_receptor, z_score_ligand, z_score_combined
    '''
    receptor_ratio, ligand_ratio = ratio
    z_score_receptor = Series(
        ((data - data.mean()) / data.std()).mean(axis=1), name='z_score_receptor')
    z_score_ligand = Series(
        (((data.T - data.T.mean()) / data.T.std()).T).mean(axis=1), name='z_score_ligand')
    z_score_combined = Series(
        receptor_ratio * z_score_receptor + ligand_ratio * z_score_ligand.fillna(0), name='Z_score_combined')

    return z_score_receptor, z_score_ligand, z_score_combined

'''
def relative(data: DataFrame):
    
    相对分数 : 与(单个)参考值的相对比值
    Parameters
    ----------
    data : DataFrame
        待计算数据

    Return
    ----------
    DataFrame
        相对分数结果数据
    
    _processed = []
    for column in data.columns:
        try:
            #ref_value = float(reference_data.loc['Reference', column])
            result = data[column] #/ ref_value
            _processed.append(result)
        except KeyError:
            _processed.append(data[column])
        
    return pd.DataFrame(_processed).T
'''
# Machine Learning

def gbt_classifier(**params):
    '''
    梯度提升树
    Parameters
    ----------
    params : dict
        超参数
    '''
    try:
        from xgboost import XGBClassifier
    except ImportError:
        raise RuntimeError('XGBoost is not installed.')
    
    return XGBClassifier(**params)

def logistic_regression(**params):
    '''
    逻辑回归
    Parameters
    ----------
    params : dict
        超参数
    '''
    try:
        from sklearn.linear_model import LogisticRegression
    except ImportError:
        raise RuntimeError('Sklearn is not installed.')
    
    return LogisticRegression(**params)

def dummy_classifier(**params):
    '''
    简单规则预测分类器

    Parameters
    ----------
    params : dict
        超参数
    '''
    try:
        from sklearn.dummy import DummyClassifier
    except ImportError:
        raise RuntimeError('Sklearn is not installed.')
    
    return DummyClassifier(**params)

def randomforest_classifier(**params):
    '''
    随机森林分类器

    Parameters
    ----------
    params : dict
        超参数
    '''
    try:
        from sklearn.ensemble import RandomForestClassifier
    except ImportError:
        raise RuntimeError('Sklearn is not installed.')
    
    return RandomForestClassifier(**params)

def naive_bayes_classifier(**params):
    '''
    朴素贝叶斯分类器

    Parameters
    ----------
    params : dict
        超参数
    '''
    try:
        from sklearn.naive_bayes import GaussianNB
    except ImportError:
        raise RuntimeError('Sklearn is not installed.')
    
    return GaussianNB(**params)