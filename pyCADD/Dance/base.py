
import logging

from pyCADD.Dance import algorithm, core
from rich.prompt import Prompt

logger = logging.getLogger('pyCADD.Dance.base')


class Dancer:
    '''
    Data analyzer for CADD
    '''

    def __init__(self, file_path=None, label_col=None) -> None:
        self.raw_data = None                            # 矩阵原始数据
        self.docking_data = None                        # 对接分数数据
        self.merge_data = None
        self.current_data = []                          # 计算完成的数据
        self.current_data_dic = {}
        self.label_col = label_col
        if not file_path:
            raise ValueError('Require matrix file path.')
        self.file_path = file_path                      # 矩阵文件路径
        self.read_data(self.file_path)

    @property
    def activity_data(self):
        return self.raw_data[self.label_col]

    def read_data(self, file_path: str):
        '''
        读取对接结果数据

        Parameter
        ---------
        file_path : str
            数据文件路径(matrix文件)
        '''

        self.raw_data = core.read_matrix(file_path)
        if not self.label_col:
            self.label_col = Prompt.ask('Enter the column name of label', choices=list(
            self.raw_data.columns), default=self.raw_data.columns[-1])
        self.docking_data = core.read_docking_data(
            self.raw_data, self.label_col)

    def mean(self, method: str = 'ave'):
        '''
        计算对接分数均值项
        Parameter
        ---------
        method : str
            均值计算方法
                ave: 算术平均值
                geo: 几何平均值
        '''
        self.mean_data = algorithm.average(self.docking_data, method)
        self.current_data.append(self.mean_data.name)
        self.current_data_dic[self.mean_data.name] = self.mean_data
        logger.debug(
            'Mean value(method: %s) has been appended to current data.' % method)

    def min(self):
        '''
        计算最小值项
        '''
        self.min_data = algorithm.minimum(self.docking_data)
        self.current_data.append(self.min_data.name)
        self.current_data_dic[self.min_data.name] = self.min_data
        logger.debug('Minimum value has been appended to current data.')

    def max(self):
        '''
        计算最大值项
        '''
        self.max_data = algorithm.maximum(self.docking_data)
        self.current_data.append(self.max_data.name)
        self.current_data_dic[self.max_data.name] = self.max_data
        logger.debug('Maximum value has been appended to current data.')

    def z_score(self, ratio: tuple = (0.7, 0.3)):
        '''
        计算z_score
        Parameter 
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
    
    def scp(self):
        '''
        原始单晶体对接得分(single-conformation performance)排序数据

        Parameter
        ---------
        reference_data : str
            参考数据文件路径(csv)
        '''
        
        self.relative_data = algorithm.relative(self.docking_data)
        for column in self.relative_data.columns:
            self.current_data.append(column)
            self.current_data_dic[column] = self.relative_data[column]
        logger.debug('SCP data has been appended to current data.')

    def merge(self, data_list):
        '''
        合并指定的数据
        Parameter
        ---------
        data_list : list
            需要合并的数据组成的列表
        '''
        logger.debug('Datasets required merge:\n%s' % '\n'.join('Name: ' + str(dataset.name) +
                     ', Length: ' + str(dataset.size) + ' ,dtype: ' + str(dataset.dtype) for dataset in data_list))
        self.merge_data = core.merge(data_list)

    def auc(self, pos_label, save: bool = False, ascending: bool = False):
        '''
        生成ROC曲线并计算AUC
        Parameters
        ----------
        pos_label : str | int | list
            支持的表示阳性标签的类型
        save : bool
            是否保存ROC曲线图文件
        ascending : bool
            是否以升序方式排序数据(针对数据为负值 越负者越优的情况)

        '''

        self.auc_data = core.get_auc(
            self.merge_data, self.label_col, pos_label, save, ascending)

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
        self.corr_data.to_csv('results/correlation.csv')
        heatmap_file = core.heatmap(self.corr_data, 0, 1, save)
        if save:
            logger.debug('%s saved.' % heatmap_file)
