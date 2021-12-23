
import logging

from pyCADD.Dance import algorithm, core
from rich.prompt import Prompt

logger = logging.getLogger('pyCADD.Dance.base')


class Dancer:
    '''
    Data analyzer for CADD
    '''

    def __init__(self, file_path=None) -> None:
        self.raw_data = None                            # 矩阵原始数据
        self.docking_data = None                        # 对接分数数据
        self.merge_data = None
        self.current_data = []                          # 计算完成的数据
        self.current_data_dic = {}
        self.file_path = file_path                      # 矩阵文件路径
        self.read_data(self.file_path)

    @property
    def activity_data(self):
        return self.raw_data[self.label_col]

    def read_data(self, file_path: str):
        '''
        读取对接结果数据
        '''

        self.raw_data = core.read_matrix(file_path)
        self.label_col = Prompt.ask('Enter the column name of label', choices=list(
            self.raw_data.columns), default=self.raw_data.columns[-1])
        self.docking_data = core.read_docking_data(
            self.raw_data, self.label_col)

    def mean(self, method: str = 'ave'):
        '''
        计算对接分数均值项
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

    def z_score(self, ratio:tuple=(0.7, 0.3)):
        '''
        计算z_score
        '''
        
        self.z_score_receptor, self.z_score_ligand, self.z_score_combined = algorithm.z_score(self.docking_data, ratio)
        
          
        for dataset in [self.z_score_receptor, self.z_score_ligand, self.z_score_combined]:
            self.current_data.append(dataset.name)
            self.current_data_dic[dataset.name] = dataset.fillna(10)

        '''            
        self.current_data.append(self.z_score_combined.name) 
        self.current_data_dic[self.z_score_combined.name] = self.z_score_combined
        '''
        logger.debug('Z-score dataset has been appended to current data.')
        
    def merge(self, data_list):
        '''
        合并指定的数据
        '''
        logger.debug('Datasets required merge:\n%s' % '\n'.join('Name: ' + str(dataset.name) +
                     ', Length: ' + str(dataset.size) + ' ,dtype: ' + str(dataset.dtype) for dataset in data_list))
        self.merge_data = core.merge(data_list)

    def auc(self, pos_label, save: bool = False, ascending: bool = False):
        '''
        生成ROC曲线并计算AUC
        '''

        self.auc_data = core.get_auc(
            self.merge_data, self.label_col, pos_label, save, ascending)
