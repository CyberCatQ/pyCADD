# 数据融合规则算法

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
    '''
    return Series(data.min(axis=1), name='Min')


def maximum(data: DataFrame):
    '''
    最大值 : 提取DataFrame数据列中的最大值

    Parameters
    ----------
    data : DataFrame
        待计算数据
    '''
    return data.max(axis=1)

def std(data: DataFrame, axis:int=1):
    '''
    计算标准偏差

    Parameters
    ----------
    data : DataFrame
        待计算数据
    axis : int
        计算坐标轴 row: 0, column: 1
    '''
    return data.std(axis=axis)

def z_score(data: DataFrame):
    '''
    标准分数(Z-score) : 对DataFrame数据进行标准分数变换
    Z-score-combined = 0.7 * Z-score-receptor + 0.3 * Z-score-ligand

    Parameters
    ----------
    data : DataFrame
        待计算数据
    
    Return
    tuple(DataFrame, DataFrame, DataFrame)
        z_score_receptor, z_score_ligand, z_score_combined
    '''

    z_score_receptor = (data - data.mean()) / data.std()
    z_score_ligand = ((data.T - data.T.mean()) / data.T.std()).T
    z_score_combined = 0.7 * z_score_receptor + 0.3 * z_score_ligand.fillna(0)

    # 对于NaN值 填充10作为惩罚值
    z_score_combined['mean'] = z_score_combined.fillna(10).mean(axis=1)
    return z_score_receptor, z_score_ligand, z_score_combined


def P_value(data: DataFrame):
    '''
    pearson相关系数矩阵

    Kendall是定类变量的统计
    pearson是对定距变量的统计
    spearman是对定序变量的统计
    '''
    return data.corr(method='pearson')


def relative(data: DataFrame):
    '''
    相对分数 : 与(单个)参考值的相对比值
    '''


def relative_ave(data: DataFrame):
    '''
    相对分数(平均) : 与多个参考值的平均的相对比值
    '''
