import logging
from pyCADD.ui import UI
from pyCADD.Dance.base import Dancer
from pyCADD.utils.check import check_file

logger = logging.getLogger('pyCADD.Dance.UI')

class UI_Dance(UI):
    '''
    数据处理程序UI
    '''

    def __init__(self, menu_name: str = 'Data Analyzer') -> None:
        super().__init__(menu_name=menu_name)
        self.dancer = Dancer()
        self.main_option = [
            '1. Read docking data matrix',
            '2. Calculate the arithmetic mean',
            '3. Calculate the geometric mean',
            '4. Calculate the minimum value',
            '5. Calculate the maximum value',
            '6. Merge datasets',
            '7. Draw ROC curve and calculate the AUC',
            '0. Exit'
        ]

    def run(self, flag):

        if flag == '1':
            file_path = input('Enter the matrix file path: ').strip()
            if not check_file(file_path):
                logger.error('File %s not found.' % file_path)
                return

            self.dancer.read_data(file_path)
            self.create_panel(additional_info='Read the matrix file: %s' % file_path)
        
        elif flag == '2':
            logger.info('Calculating the arithmetic mean of matrix')
            self.dancer.mean(method='ave')
            logger.info('Arithmetic mean calculate done.')

            self.clear_info()
            self.create_panel(additional_info='Current datasets: %s' % self.dancer.current_data)
        
        elif flag == '3':
            logger.info('Calculating the geometric mean of matrix')
            self.dancer.mean(method='geo')
            logger.info('Geometric mean calculate done.')

            self.clear_info()
            self.create_panel(additional_info='Current datasets: %s' % self.dancer.current_data)

        
        elif flag == '4':
            logger.info('Extracting the minimum value')
            self.dancer.min()
            logger.info('Minimum value Extracted.')

            self.clear_info()
            self.create_panel(additional_info='Current datasets: %s' % self.dancer.current_data)

        elif flag == '5':
            logger.info('Extracting the maximum value')
            self.dancer.max()
            logger.info('Maximum value Extracted.')

            self.clear_info()
            self.create_panel(additional_info='Current datasets: %s' % self.dancer.current_data)

        elif flag == '6':
            logger.info('Current data: %s' % self.dancer.current_data)
            if not self.get_confirm('Merge all current data?'):
                _need_merge = input('Enter the datasets need to be merge(separated by commas): ').split(',')
            else:
                _need_merge = self.dancer.current_data
            logger.info('Merging datasets: %s' % _need_merge)

            data_list = [self.dancer.current_data_dic[index] for index in _need_merge]
            data_list.append(self.dancer.activity_data)
            self.dancer.merge(data_list)

            self.create_panel(additional_info='Merged datasets: %s' % _need_merge)
        
        elif flag == '7':

            logger.info('Label column: %s' % self.dancer.label_col)
            labels = list(self.dancer.activity_data.value_counts().index)
            logger.info('Labels in column: %s' % labels)

            pos_label = input('Enter the positive labels(separated by commas): ').split(',')
            self.dancer.auc(pos_label, True, self.get_confirm('Ascending sort(for negative number)?'))
            logger.info('ROC curve image file %s-ROC.jpg saved.' % self.dancer.label_col)

            self.create_panel()
        
if __name__ == '__main__':
    enter_text = '[bold]Enter the Code of Options'
    ui_dance = UI_Dance()
    ui_dance.create_panel(ui_dance.main_option)

    while True:
        flag = ui_dance.get_input(enter_text, [str(i) for i in range(0, 8)])
        if flag == '0':
            break
        
        ui_dance.run(flag)






        




    