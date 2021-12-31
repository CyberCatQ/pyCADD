# 数据融合规则算法
import pandas as pd
from pandas import DataFrame, Series


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


def relative(data: DataFrame):
    '''
    相对分数 : 与(单个)参考值的相对比值
    Parameters
    ----------
    data : DataFrame
        待计算数据

    Return
    ----------
    DataFrame
        相对分数结果数据
    '''
    _processed = []
    for column in data.columns:
        try:
            #ref_value = float(reference_data.loc['Reference', column])
            result = data[column] #/ ref_value
            _processed.append(result)
        except KeyError:
            _processed.append(data[column])
        
    return pd.DataFrame(_processed).T