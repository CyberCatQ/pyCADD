import logging
import os

from pyCADD.Dance.base import Dancer
from pyCADD.utils.ui import UI

logger = logging.getLogger('pyCADD.Dance')
enter_text = '[bold]Enter the Code of Options'

class UI_Dance(UI):
    '''
    数据处理程序UI
    '''

    def __init__(self, menu_name: str = 'Data Analyzer') -> None:
        super().__init__(menu_name=menu_name)
        self.merged = False
        self.dancer = None
        self.check_label = False
        self.main_options = [
            '1. Read docking data matrix',
            '2. Calculate',
            '3. Draw Plots',
            '0. Exit'
        ]

        self.calc_options = [
            '1. Calculate the arithmetic mean',
            '2. Calculate the geometric mean',
            '3. Calculate the minimum value',
            '4. Calculate the maximum value',
            '5. Calculate Z-score',
            '6. Calculate single-conformation performance',
            '0. Back'
        ]

        self.draw_options = [
            '1. Draw ROC curve and calculate the AUC',
            '2. Draw Scatter',
            '3. Draw Correlation Heatmap',
            '0. Back'
        ]

    def _print_current_data(self):
        self.clear_info()
        self.create_panel(additional_info='Current datasets: %s' %
                          self.dancer.current_data)

    def _merge_data(self):
        '''
        合并已计算数据集
        '''
        logger.info('Current data: %s' % self.dancer.current_data)
        if not self.get_confirm('Merge all current data?'):
            _need_merge = input(
                'Enter the datasets need to be merge(separated by commas): ').split(',')
        else:
            _need_merge = self.dancer.current_data
        logger.info('Merging datasets: %s' % _need_merge)

        data_list = [self.dancer.current_data_dic[index]
                     for index in _need_merge]
        data_list.append(self.dancer.activity_data)
        self.dancer.merge(data_list)
        self.merged = True

        self.create_panel(additional_info='Merged datasets: %s' %
                          _need_merge, show_panel=False)

    def _get_label_info(self):
        '''
        获取标签相关基本信息
        '''
        logger.info('Label column: %s' % self.dancer.label_col)
        labels = list(self.dancer.activity_data.value_counts().index)
        logger.info('Labels in column: %s' % labels)
        self.pos_label = input(
            'Enter the positive labels(separated by commas): ').split(',')

        # 简易校验
        for label in self.pos_label:
            if label not in labels:
                self.create_panel()
                logger.error('%s is not in the label column.' % label)
                return
        logger.info('Positive labels: %s' % self.pos_label)
        self.check_label = True

    def run(self, flag):
        '''
        主菜单
        '''
        if flag == '1':
            default_matrix_file = os.path.abspath('results/matrix.csv')
            if os.path.exist(default_matrix_file):
                if self.get_confirm('Dectected results file: %s \nUse it?' % default_matrix_file):
                    file_path = default_matrix_file
                else:
                    file_path = input('Enter the matrix file path').strip()
            else:
                file_path = input('Enter the matrix file path: ').strip()

            if not os.path.exist(file_path):
                self.create_panel()
                logger.error('File %s not found.' % file_path)
                return

            self.dancer = Dancer(file_path)
            self.create_panel(
                additional_info='Read the matrix file: %s' % file_path)
        else:
            if self.dancer is None:
                self.create_panel()
                logger.error('No matrix dataset read.')
                return

        if flag == '2':
            self.create_panel(self.calc_options)
            while True:
                flag = self.get_input(enter_text, choices=[str(
                    i) for i in range(len(self.calc_options))], default='0')
                if flag == '0':
                    self.create_panel(self.main_options)
                    return
                self.calculate(flag)

        elif flag == '3':
            self.create_panel(self.draw_options)
            while True:
                flag = self.get_input(enter_text, choices=[str(
                    i) for i in range(len(self.draw_options))], default='0')
                if flag == '0':
                    self.create_panel(self.main_options)
                    return
                self.plot(flag)

    def calculate(self, flag):
        '''
        计算菜单
        '''
        if flag == '1':
            logger.info('Calculating the arithmetic mean of matrix')
            self.dancer.mean(method='ave')
            logger.info('Arithmetic mean calculate done.')

            self._print_current_data()

        elif flag == '2':
            logger.info('Calculating the geometric mean of matrix')
            self.dancer.mean(method='geo')
            logger.info('Geometric mean calculate done.')

            self._print_current_data()

        elif flag == '3':
            logger.info('Extracting the minimum value')
            self.dancer.min()
            logger.info('Minimum value Extracted.')
            self._print_current_data()

        elif flag == '4':
            logger.info('Extracting the maximum value')
            self.dancer.max()
            logger.info('Maximum value Extracted.')
            self._print_current_data()

        elif flag == '5':
            receptor_ratio = float(self.get_input(
                'Enter the ratio of receptor z_score', default='7'))
            ligand_ratio = float(self.get_input(
                'Enter the ratio of ligand z_score', default='3'))
            total = receptor_ratio + ligand_ratio

            receptor_ratio /= total
            ligand_ratio /= total
            ratio = (receptor_ratio, ligand_ratio)
            self.dancer.z_score(ratio)

            self._print_current_data()
            logger.info('Calculating Z-score with %s : %s' % ratio)
            logger.info('Z-score calculate done.')
        
        elif flag == '6':
            logger.info('Calculating single-conformation performance')
            ''' 
            default_reference = os.path.abspath('results/reference.csv')
            if os.path.exist(default_reference):
                if self.get_confirm('Detected reference file: %s \nUse it?' % default_reference):
                    reference_file = default_reference
                else:
                    reference_file = self.get_input('Enter the referenced data file path')
            else:
                reference_file = self.get_input('Enter the referenced data file path')
            '''
            self.dancer.scp()
            self._print_current_data()
            logger.info('SCP calculate done.')

    def plot(self, flag):
        '''
        绘图菜单
        '''
        while not self.check_label:
            self._get_label_info()

        if flag == '1':
            self._merge_data()
            if not self.merged:
                self.create_panel()
                logger.error('No merged dataset.')
                return

            ascending = self.get_confirm(
                'Ascending sort(for negative number)?')
            logger.info('Sort by ascending order: %s' % ascending)

            self.dancer.auc(self.pos_label, True, ascending)
            self.create_panel()
            logger.info('ROC curve image file %s-ROC.jpg saved.' %
                        os.path.basename(os.getcwd()))

        elif flag == '2':

            score_name = self.get_input(
                'Enter the name of score', default='Docking_Score')
            logger.info('Score Name: %s' % score_name)

            self.dancer.scatter(self.pos_label, score_name, True)
            self.create_panel()
            logger.info('Scatter plot %s-Scatter.jpg saved.' %
                        os.path.basename(os.getcwd()))

        elif flag == '3':

            method = self.get_input('Enter the method of correlation', choices=[
                                    'pearson', 'kendall', 'spearman'], default='pearson', show_choices=True)
            logger.info('Correlation method: %s' % method)
            self.dancer.correlate(method, True)
            self.create_panel()
            logger.info('Heatmap %s-Correlation.jpg saved.' %
                        os.path.basename(os.getcwd()))
