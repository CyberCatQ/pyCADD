import numpy as np
import logging
import json
import os
from pandas import Series

from pyCADD.Dance import algorithm, core
from rich.prompt import Prompt, Confirm
from sklearn.model_selection import train_test_split

logger = logging.getLogger('pyCADD.Dance.base')


class Dancer:
    '''
    Data analyzer for CADD
    '''

    def __init__(self, file_path=None, label_col=None) -> None:
        self.raw_data = None                            # 矩阵原始数据
        self.docking_data = None                        # 对接分数数据(X)
        self.merge_data = None                          # 合并后的数据
        self.current_data = []                          # 计算完成的数据
        self.current_data_dic = {}                      # 计算完成的数据字典{name: data}
        self.label_col = label_col                      # 标签列名
        self.label = None                               # 标签列(y)
        if not file_path:
            raise ValueError('Require matrix file path.')
        self.file_path = file_path                      # 矩阵文件路径
        self.read_data(self.file_path)

    @property
    def activity_data(self):
        return self.label

    def read_data(self, file_path: str):
        '''
        读取对接结果数据

        Parameters
        ---------
        file_path : str
            数据文件路径(matrix文件)
        '''

        self.raw_data = core.read_matrix(file_path)
        if not self.label_col:
            self.label_col = Prompt.ask('Enter the column name of label', choices=list(
                self.raw_data.columns), default=self.raw_data.columns[-1])
        self.docking_data, self.label = core.split_data(
            self.raw_data, self.label_col)

    def mean(self, method: str = 'ave', ignore_nan: bool = True):
        '''
        计算对接分数均值项

        Parameters
        ---------
        method : str
            均值计算方法
                * ave: 算术平均值
                * geo: 几何平均值
        ignore_nan : bool
            是否忽略NaN值(默认True)
        '''
        _data_to_mean = self.docking_data.replace(
            0, np.nan) if ignore_nan else self.docking_data
        self.mean_data = algorithm.average(_data_to_mean, method)
        self.current_data.append(self.mean_data.name)
        self.current_data_dic[self.mean_data.name] = self.mean_data
        logger.debug(
            'Mean value(method: %s) has been appended to current data.' % method)

    def min(self, ignore_nan: bool = True):
        '''
        计算最小值项

        Parameters
        ---------
        ignore_nan : bool
            是否忽略NaN值(默认True)
        '''
        _data_to_min = self.docking_data.replace(
            0, np.nan) if ignore_nan else self.docking_data
        self.min_data = algorithm.minimum(_data_to_min)
        self.current_data.append(self.min_data.name)
        self.current_data_dic[self.min_data.name] = self.min_data
        logger.debug('Minimum value has been appended to current data.')

    def max(self, ignore_nan: bool = True):
        '''
        计算最大值项

        Parameters
        ---------
        ignore_nan : bool
            是否忽略NaN值(默认True)
        '''
        _data_to_max = self.docking_data.replace(
            0, np.nan) if ignore_nan else self.docking_data
        self.max_data = algorithm.maximum(_data_to_max)
        self.current_data.append(self.max_data.name)
        self.current_data_dic[self.max_data.name] = self.max_data
        logger.debug('Maximum value has been appended to current data.')

    def z_score(self, ratio: tuple = (0.7, 0.3)):
        '''
        计算z_score

        Parameters 
        ---------
        ratio : tuple
            组合Z值中 receptor平均Z值 与 ligand平均Z值的比例
        '''

        self.z_score_receptor, self.z_score_ligand, self.z_score_combined = algorithm.z_score(
            self.docking_data, ratio)

        for dataset in [self.z_score_receptor, self.z_score_ligand, self.z_score_combined]:
            self.current_data.append(dataset.name)
            self.current_data_dic[dataset.name] = dataset.fillna(10)

        '''            
        self.current_data.append(self.z_score_combined.name) 
        self.current_data_dic[self.z_score_combined.name] = self.z_score_combined
        '''
        logger.debug('Z-score dataset has been appended to current data.')

    def best_scp(self):
        '''
        原始最佳单晶体对接得分(best single-conformation performance) AUC数据
        全数据集验证
        '''

        self.best_SCP, self.best_SCP_score = core.get_best_SCP(
            self.docking_data, self.label)
        self.current_data.append('Best_SCP-%s' % self.best_SCP)
        self.current_data_dic['Best_SCP-%s' %
                              self.best_SCP] = self.docking_data[self.best_SCP]
        logger.debug(
            'Best SCP-%s data has been appended to current data.' % self.best_SCP)

    def merge(self, data_list):
        '''
        合并指定的数据

        Parameters
        ---------
        data_list : list
            需要合并的数据组成的列表
        '''
        logger.debug('Datasets required merge:\n%s' % '\n'.join('Name: ' + str(dataset.name) +
                     ', Length: ' + str(dataset.size) + ' ,dtype: ' + str(dataset.dtype) for dataset in data_list))
        self.merge_data = core.merge(data_list)

        return self.merge_data

    def auc(self, save: bool = False, lower_is_better: bool = True):
        '''
        生成ROC曲线并计算AUC

        Parameters
        ----------
        pos_label : str | int | list
            支持的表示阳性标签的类型
        save : bool
            是否保存ROC曲线图文件
        lower_is_better : bool
            是否为负样本的ROC曲线(默认为True)
        '''

        self.auc_data = core.get_roc(
            self.merge_data, self.activity_data, save, lower_is_better)

    def scatter(self, pos_label, score_name: str = 'Docking_Score', save: bool = False):
        '''
        生成散点分布图

        Parameters
        ----------
        pos_label : str | int | list
            支持的表示阳性标签的类型
        score_name : str
            数据值的名称
        save : bool
            是否保存散点图文件

        '''
        scatter_file = core.get_scatter(self.raw_data, self.label_col,
                                        pos_label, score_name, save)
        if save:
            logger.debug('%s saved.' % scatter_file)

    def correlate(self, method: str = 'pearson', save: bool = False):
        '''
        生成相关性系数热力图

        Parameters
        ----------
        method : str
            计算相关系数的方法
                pearson, kendall, spearman
        save : bool
            是否保存热力图文件

        '''
        logger.debug('correlation method: %s' % method)
        self.corr_data = core.correlate(self.docking_data)
        # self.corr_data.to_csv('results/correlation.csv')
        heatmap_file = core.heatmap(self.corr_data, 0, 1, save)
        if save:
            logger.debug('%s saved.' % heatmap_file)


class Dancer_ML(Dancer):
    '''
    Dancer for Machine Learning
    '''

    def __init__(self, file_path=None, label_col=None) -> None:
        super().__init__(file_path, label_col)
        self.support_methods = ['Dummy', 'GBT', 'LR', 'RF', 'NBC']
        self.method = None
        self.model = None
        self.params_file = None
        self.params = None
        self.splits = None
        self.RANDOM_STATE = 42
        self.test_size = 0.2
        self._train_test_split(self.test_size)

        self.eva_models = {}
        self.evaluation_results = {}

    @property
    def X(self):
        return self.docking_data

    @property
    def y(self):
        return self.activity_data

    def _train_test_split(self, test_size: float = 0.2):
        '''
        训练集 测试集拆分
        '''
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            self.docking_data, self.activity_data, test_size=test_size, random_state=self.RANDOM_STATE, stratify=self.activity_data)

    def get_splits(self, repeats: int = 30, cv: int = 4):
        '''
        得到多重交叉验证的训练集和测试集的索引
        '''
        self.splits = core.get_splits(
            self.docking_data, self.activity_data, repeats, cv, self.RANDOM_STATE)

    def set_method(self, method: str):
        '''
        设定ML方法
        '''
        self.model = None
        self.params_file = None
        self.parms = None

        if method not in self.support_methods:
            raise ValueError('Method %s is not supported.' % method)

        if method == 'GBT':
            self.method = 'ml_XGBClassifier'
            self.params_file = 'best_params_XGBClassifier.json'
            self.model = algorithm.gbt_classifier()
        elif method == 'LR':
            self.method = 'ml_LogisticRegression'
            self.params_file = 'best_params_LogisticRegression.json'
            self.model = algorithm.logistic_regression()
        elif method == 'Dummy':
            self.method = 'ml_DummyClassifier'
            self.params_file = 'best_params_DummyClassifier.json'
            self.model = algorithm.dummy_classifier(strategy='stratified')
        elif method == 'RF':
            self.method = 'ml_RandomForestClassifier'
            self.params_file = 'best_params_RandomForestClassifier.json'
            self.model = algorithm.randomforest_classifier()
        elif method == 'NBC':
            self.method = 'ml_GaussianNB'
            self.params_file = 'best_params_GaussianNB.json'
            self.model = algorithm.naive_bayes_classifier()

        self.eva_models[self.method] = self.model
        logger.debug('Method %s has been set.' % self.method)

    def check_params_file(self):
        '''
        尝试检查参数文件是否存在
        '''
        if os.path.exists(self.params_file):
            if Confirm.ask('Params file %s exists, Use it?' % self.params_file, default=True):
                self.params = json.load(open(self.params_file))
                self.set_params(self.params)
                return True
            else:
                logger.info('Params file %s has been ignored.' %
                            self.params_file)
                return False
        else:
            logger.info('Params file %s does not exists.' % self.params_file)
            return False

    def set_params(self, params: dict = None):
        '''
        设定超参数

        Parameters
        ---------
        params : dict 
            超参数字典
        '''
        if params is not None:
            self.params = params
        else:
            params = self.params
        self.model.set_params(**params)
        logger.info('Params have been set:\n %s' % self.model.get_params())

    def hyperparam_tuning(self, param_grid: dict, method: str = 'grid', *args, **kwargs):
        '''
        超参数调优

        Parameters
        ----------
        param_grid : dict
            超参数取值范围字典
        method : str
            调优方法
                * grid
                * random
        '''
        logger.debug('Current Model: %s' % self.method)
        logger.debug('Params Grid:\n%s' % param_grid)
        self.params = core.hyperparam_tuning(
            self.model, param_grid, self.X_train, self.y_train, method=method, save=True, **kwargs)
        logger.debug('Best Params:\n%s' % self.params)
        self.set_params(self.params)

    def train(self):
        '''
        训练模型
        '''
        if not self.params:
            logger.info('Params have not been set, use default params.')

        self.model.fit(self.X_train, self.y_train)
        logger.info('Model %s has been trained.' % self.method)

    def add_consensus_scoring(self):
        '''
        结果中加入共识评分方法:
            * 平均值
            * 几何平均值
            * 最小值
        '''
        self.eva_models['cs_mean'] = algorithm.Average(lower_is_better=True)
        self.eva_models['cs_geo'] = algorithm.Geo_Average(lower_is_better=True)
        self.eva_models['cs_min'] = algorithm.Minimum(lower_is_better=True)

    def scp_report(self):
        '''
        SCP策略交叉测试报告
        '''
        self.scp_performance = core.get_SCP_report(self.splits, self.X, self.y)

    def models_evaluation_report(self):
        '''
        30x4模型交叉验证评估报告
        '''
        if not self.eva_models:
            logger.error('No model has been set.')
            return
        self.evaluation_results = core.CV_model_evaluation(
            self.eva_models, self.X, self.y, random_state=self.RANDOM_STATE)

    def roc_auc(self, save: bool = False):
        '''
        ROC曲线
        '''
        model_predicts = []
        # 测试集验证
        X = self.X_test
        y = self.y_test

        for model_name, model in self.eva_models.items():
            if model_name.startswith('cs_'):
                model_predicts.append(model(X.replace(0, np.nan)))
            elif model_name.startswith('ml_'):
                model_predicts.append(Series(model.predict_proba(
                    X)[:, 1] * -1, name=model_name, index=X.index))

        core.get_roc(self.merge(model_predicts), y, save=save)
