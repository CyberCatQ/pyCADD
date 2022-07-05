from typing import Literal

import numpy as np
import pandas as pd
from pandas import DataFrame, Series


class _Consensus:
    '''
    共识算法基类
    新的共识性方法(基于公式的方法)应该继承此类
    '''

    def __init__(self, lower_is_better=False):
        '''
        Parameters
        ----------
        lower_is_better : bool
            数据特征是否为下降性指标
        '''
        self.result = None
        self.lower_is_better = lower_is_better

        self.params = {
            'class' : 'Consensus Strategy',
            'lower_is_better': self.lower_is_better,
            }

    def get_params(self):
        '''
        获取参数 为兼容sklearn的参数设计
        '''
        return self.params

    def fit(self):
        '''
        共识算法拟合 待重载实现
        '''
        raise NotImplementedError

    def predict(self, X: DataFrame, y: Series = None, limit_num: int = None):
        '''
        按照算法进行预测

        Parameters
        ----------
        X : DataFrame
            输入数据
        y : Series
            输入数据的标签
        limit_num : int
            预测结果的数量n 排名前n的结果将被断言为阳性

        Returns
        -------
        Series
            预测结果
        '''
        # 共识性方法实际并不进行拟合
        # 而是根据公式计算出预测结果

        if X is not None:
            self.fit(X, y)
        if self.result is None:
            raise RuntimeError(
                f'Result is None while using {self.__class__.__name__} method to predict. X is required.')

        _predict_df = self.result.reset_index()
        _predict_df.columns = ['index', 'probability']
        _predict_df.sort_values(by='probability', ascending=False, inplace=True)
        _predict_df['prediction'] = np.nan

        if limit_num is None:
            if y is None:
                raise ValueError('y is required when limit_num is None.')
            print(
                f'No limit_num specified, TOP {y.sum()} results (same as the number of true positive samples) will be predicted as positive.')
            limit_num = y.sum()

        _predict_df.iloc[:limit_num, -1] = 1
        _predict_df.iloc[limit_num:, -1] = 0

        _predict_df.sort_index(inplace=True)
        _predict_df.set_index('index', inplace=True)

        return _predict_df.loc[:, 'prediction']

    def predict_proba(self, X: DataFrame = None, y: Series = None):
        '''
        按照算法进行概率预测
        '''
        # 共识性方法实际并不进行拟合
        # 而是根据公式计算出预测结果 没有概率
        # 直接返回计算结果的Numpy array

        if X is not None:
            self.fit(X, y)
        if self.result is None:
            raise RuntimeError(
                f'Result is None while using {self.__class__.__name__} method to predict. X is required.')

        # array [[name, value], [name, value], ...]
        # 为兼容在sklearn中的predict_proba 第0维度为样本名(sklearn中为预测值) 第1维度为值(sklearn中为概率)
        # 第1维度为主要排序的维度 并用于生成AUC等
        return pd.DataFrame(self.result).reset_index().to_numpy()


class Average(_Consensus):
    '''
    算数平均值模型
    算术平均数定义为：
        平均数 = 求和(x_i) / 总数
    '''

    def __init__(self, lower_is_better=False):
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({'method': 'Arithmetic mean'})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True):
        '''
        计算算数平均
        Parameters
        ----------
        X : DataFrame
            数据特征 即用于计算平均值的数据
        y : Series
            数据标签 实际不使用在当前计算中 为了兼容性而存在
        ignore_nan : bool
            是否忽视缺失值(所有0值被视为缺失值)
        '''
        # 在预处理阶段 X的所有NaN被填充为0
        # 在计算阶段 X的所有0被重新填充为NaN 在计算Mean时将被忽略
        X = X.replace(0, np.nan) if ignore_nan else X

        # mean计算自动忽略NaN
        self.result = X.mean(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result


class Mean(Average):
    '''
    算术平均值模型
    Mean作为Average的别称
    '''


class Geo_Average(_Consensus):
    '''
    几何平均值模型
    几何平均值定义为：
        几何平均数 = 连续乘积(x_i)^(1/n)
    '''

    def __init__(self, lower_is_better=False):
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({'method': 'Geometric mean'})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True):
        '''
        计算几何平均
        取绝对值进行计算后 再赋予符号
        Parameters
        ----------
        X : DataFrame
            数据特征 即用于计算平均值的数据
        y : Series
            数据标签 实际不使用在当前计算中 为了兼容性而存在
        ignore_nan : bool
            是否忽视缺失值(所有0值被视为缺失值)
        '''
        # 在预处理阶段 X的所有NaN被填充为0
        # 在计算阶段 X的所有0被重新填充为NaN 在计算Mean时将被忽略
        X = X.replace(0, np.nan) if ignore_nan else X

        # 取绝对值进行计算后 再赋予符号
        self.result = Series(
            pow(
                # 绝对值连续乘积
                X.prod(axis=1).abs(),
                # 开项数次根
                1/X.notna().sum(axis=1)
            ),
            name='GEO'
        )

        self.result = self.result if self.lower_is_better else -1 * self.result


class GeoMean(Geo_Average):
    '''
    几何平均值模型
    GeoMean作为Geo_Average的别称
    '''


class Minimum(_Consensus):
    '''
    最小值模型
    '''

    def __init__(self, lower_is_better=False):
        super().__init__(lower_is_better)
        self.result = None
        self.params.update({'method': 'Minimum'})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True):
        '''
        取得最小值
        '''
        X = X.replace(0, np.nan) if ignore_nan else X
        self.result = X.min(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result


class Maximum(_Consensus):
    '''
    最大值模型
    '''

    def __init__(self, X=None, y=None, lower_is_better=False):
        super().__init__(X, y, lower_is_better)
        self.result = None
        self.params.update({'method': 'Maximum'})

    def fit(self, X: DataFrame, y: Series = None, ignore_nan: bool = True):
        '''
        取得最大值
        '''
        X = X.replace(0, np.nan) if ignore_nan else X
        self.result = X.max(axis=1)
        if self.lower_is_better:
            self.result = -1 * self.result

# 函数接口


def average(data: DataFrame, method: Literal['ave', 'geo'] = 'ave'):
    '''
    平均值算法 : 计算并生成DataFrame对接数据的平均值列  
    ave 算数平均
    geo 几何平均

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
