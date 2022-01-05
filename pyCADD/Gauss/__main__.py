import os
import logging
from pyCADD.utils.check import check_file
from pyCADD.ui import UI
from pyCADD.Gauss.base import Gauss

logger = logging.getLogger('pyCADD.Gauss')
access_DFT = ['B3LYP', 'PBE0', 'CAM-B3LYP']
access_basis_set = ['6-31g*', '6-311g*', 'def2SVP', 'def2TZVP', 'def2TZVPP']
access_solvent = ['water', 'methanol', 'ethanol']


class UI_Gauss(UI):
    '''
    Gauss 计算脚本UI
    '''

    def __init__(self, menu_name: str = 'Gauss Calculate') -> None:
        super().__init__(menu_name=menu_name)

        self.main_options = [
            '1. Set Charge & spin multiplicity',
            '2. Select DFT, basis set and PCM solvent',
            '3. Gaussian structure optimization',
            '4. Single point energy calculation',
            '5. Excited state calculation',
            '0. Exit'
        ]

        self._init_setting = False
        self._system_init()
        self._get_original_st()
        self.gauss = Gauss(self.origin_st, self.cpu_count, self.mem)

    def _system_init(self):
        '''
        设定计算核心数量与内存大小
        '''
        self.cpu_count = self.get_input('Enter the number of CPU to be used', 
                                        choices=[str(i + 1) for i in range(os.cpu_count())], 
                                        default=str(os.cpu_count()))
        self.mem = self.get_input(
            'Enter the memory size to be used(GB)', default='16') + 'GB'
        self.create_panel(additional_info='[bright_cyan]CPU: %s   Memory: %s' % (
            self.cpu_count, self.mem), show_panel=False)

    def _get_original_st(self):
        '''
        获取并检查原始结构文件路径
        '''

        origin_st = self.get_input(
            'Need an original structure for calculation.\nEnter the file path of molecule structure')
        if check_file(origin_st):
            self.origin_st = origin_st
            logger.debug('Origin file: %s' % self.origin_st)
            self.create_panel(
                additional_info='Original structure file: %s' % self.origin_st, show_panel=False)
        else:
            self.create_panel()
            logger.error('File %s not found.' % origin_st)
            return

    def _basic_init(self):
        '''
        获取与设定电荷量与自旋多重度
        '''

        charge = self.get_input(
            'Enter the charge of molecule(default to 0)', default='0')
        spin_multi = self.get_input(
            'Enter the spin multiplicity(default to 1)', default='1')

        self.gauss.set_charge(int(charge))
        self.gauss.set_multispin(int(spin_multi))

        self._init_setting = True

    def _method_init(self):
        '''
        选择泛函、基组与溶剂模型
        '''
        self.dft = self.get_input(
            'Select DFT', choices=access_DFT, default='B3LYP', show_choices=True)
        self.basis_set = self.get_input(
            'Select basis set', choices=access_basis_set, default='6-31g*', show_choices=True)
        self.solvent = self.get_input(
            'Select solvent for PCM', choices=access_solvent, default='water', show_choices=True)

    def run(self, flag):
        '''
        业务逻辑
        '''

        if flag == '1':
            self._basic_init()
            self.create_panel(
                additional_info='[bright_cyan]Charge: %s' % self.gauss.charge, show_panel=False)
            self.create_panel(
                additional_info='[bright_cyan]Spin multiplicity: %s' % self.gauss.multispin)
        elif flag == '2':
            self._method_init()
            self.create_panel(additional_info='[bright_cyan]DFT: %s  Basis set: %s  Solvent: %s' % (
                self.dft, self.basis_set, self.solvent))

        if not self._init_setting:
            logger.info('Charge and spin multiplicity are not defined.')
            if not self.get_confirm('Use default charge(0)/spin multiplicity(1)?'):
                self.create_panel()
                logger.error('Please set charge/spin multiplicity firstly.')
                return
            else:
                self._basic_init(0, 1)

        if flag == '3':

            logger.info('Gauss input file %s is created.')

if __name__ == '__main__':
    enter_text = '[bold]Enter the Code of Options'
    ui_gauss = UI_Gauss()
    ui_gauss.create_panel(ui_gauss.main_options)

    while True:
        flag = ui_gauss.get_input(enter_text, choices=[str(i) for i in range(len(ui_gauss.main_options))], default='0')
        if flag == '0':
            break
        
        ui_gauss.run(flag)