from pyCADD.utils.ui import UI
from pyCADD.Gauss.base import Gauss
import os
import logging
logger = logging.getLogger('pyCADD.Gauss')

access_DFT = ['B3LYP', 'PBEPBE', 'CAM-B3LYP', 'PBE1PBE']
access_basis_set = ['6-31g*', '6-311g*', 'def2SVP', 'def2TZVP', 'def2TZVPP']
access_solvent = ['water', 'methanol', 'ethanol', 'None']

class UI_Gauss(UI):
    '''
    Gauss 计算脚本UI
    '''

    def __init__(self, menu_name: str = 'Gaussian Calculate', original_st:str = None) -> None:
        super().__init__(menu_name=menu_name)

        self.main_options = [
            '1. Set Charge & spin multiplicity',
            '2. Select DFT, basis set and PCM solvent',
            '3. Calculation system setting',
            '4. Gaussian structure optimization',
            '5. Single point energy calculation',
            '6. Absorption energy of excited states calculation',
            '7. Emission energy of excited states calculation',
            '8. Extract MO cube files',
            '0. Exit'
        ]

        self._init_setting = False
        self._init_method = False
        self._current_system_setting()
        if not original_st:
            self._get_original_st()
        else:
            self.origin_st = original_st
        logger.debug('Origin file: %s' % self.origin_st)
        self.create_panel(additional_info={ "origin" : 'Original structure file: [bright_cyan]%s[/]' % self.origin_st}, show_panel=False)
        
        self.gauss = Gauss(self.origin_st)
        self.create_panel(self.main_options, additional_info = {"file_type": 'File type: [bright_cyan]%s[/]' % self.gauss.file_type})

    @property
    def _system_loading(self):
        '''
        当前系统计算资源信息
        Return
        ----------
        cpu(s), memory
        '''
        return Gauss.system_info
    
    @property
    def cpu_count(self):
        '''
        CPU核心数量
        '''
        return self._system_loading[0]
    
    @property
    def mem(self):
        '''
        内存大小
        '''
        return self._system_loading[1]

    def _current_system_setting(self):
        '''
        当前实际系统计算资源信息
        '''
        self.create_panel(additional_info={'system': 'CPU: [bright_cyan]%s[/]   Memory: [bright_cyan]%s[/]' % (self.cpu_count, self.mem)}, show_panel=False)

    def set_system(self):
        '''
        设定计算核心数量与内存大小
        '''
        cpu_count = self.get_input('Enter the number of CPU to be used', 
                                        choices=[str(i + 1) for i in range(os.cpu_count())], 
                                        default=str(os.cpu_count()))
        mem = self.get_input(
            'Enter the memory size to be used(GB)', default='16') + 'GB'
        self.gauss.set_system(cpu_count, mem)
        self.create_panel(additional_info={'system' : 'CPU: [bright_cyan]%s[/]   Memory: [bright_cyan]%s[/]' % (cpu_count, mem)}, show_panel=False)

    def _get_original_st(self):
        '''
        获取并检查原始结构文件路径
        '''

        origin_st = self.get_input(
            'Need an original structure for calculation.\nEnter the file path of molecule structure')
        if os.path.exists(origin_st):
            self.origin_st = origin_st
        else:
            raise FileNotFoundError('File %s not found.' % origin_st)

    def _basic_init(self, charge:'int | str'=None, spin_multi:'int | str'=None):
        '''
        获取与设定电荷量与自旋多重度
        '''
        if charge is None:
            charge = self.get_input('Enter the charge of molecule(default to 0)', default='0')
        if spin_multi is None:
            spin_multi = self.get_input('Enter the spin multiplicity(default to 1)', default='1')

        self.gauss.set_charge(int(charge))
        self.gauss.set_multiplicity(int(spin_multi))
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
        
        self.gauss.set_DFT(self.dft)
        self.gauss.set_basis_set(self.basis_set)
        self.gauss.set_solvent(self.solvent)

        self._init_method = True

    def run(self, flag):
        '''
        业务逻辑
        '''

        if flag == '1':
            self._basic_init()
            self.create_panel(
                additional_info={"charge" : 'Charge: [bright_cyan]%s[/]' % self.gauss.charge}, show_panel=False)
            self.create_panel(
                additional_info={"spin_multi" : 'Spin multiplicity: [bright_cyan]%s[/]' % self.gauss.spin_multi})
        elif flag == '2':
            self._method_init()
            self.create_panel(additional_info={"method":'DFT: [bright_cyan]%s[/]  Basis set: [bright_cyan]%s[/]  Solvent: [bright_cyan]%s[/]' % (self.dft, self.basis_set, self.solvent)})
        elif flag == '3':
            self.set_system()
            self.create_panel()
            logger.info('System setting changed.')
        elif flag == '8':
            pass
        else:
            if not self._init_setting:
                logger.info('Charge and spin multiplicity are not defined.')
                if not self.get_confirm('Use default charge(0)/spin multiplicity(1)?'):
                    self.create_panel()
                    logger.error('Please set charge/spin multiplicity firstly.')
                    return
                else:
                    self._basic_init(0, 1)

            if not self._init_method:
                logger.info('DFT and basis set are not defined.')
                self._method_init()
                
        if flag == '4':
            loose = self.get_confirm('Use loose mode?')
            self.input_file = self.gauss.create_inputfile('opt', loose)
            logger.info('Gauss input file %s is created.' % self.input_file)
            if self.get_confirm('Running optimization job?'):
                self.gauss.run()
            #self.create_panel()
        
        elif flag == '5':
            self.input_file = self.gauss.create_inputfile('energy')
            logger.info('Gauss input file %s is created.' % self.input_file)
            if self.get_confirm('Running energy calculation job?'):
                self.gauss.run()
            #self.create_panel()
        
        elif flag == '6':
            self.input_file = self.gauss.create_inputfile('absorb')
            logger.info('Gauss input file %s is created.' % self.input_file)
            if self.get_confirm('Running absorption energy calculation job?'):
                self.gauss.run()
            #self.create_panel()

        elif flag == '7':
            self.input_file = self.gauss.create_inputfile('emission')
            logger.info('Gauss input file %s is created.' % self.input_file)
            if self.get_confirm('Running emission energy calculation job?'):
                self.gauss.run()
            #self.create_panel()
        
        elif flag == '8':
            if self.gauss.file_type != 'format check point(.fchk)':
                logger.error('Require a fchk file.')
                return

            self.gauss.get_mo_info()
            logger.info('Current HOMO: %s, energy: %s eV' % (self.gauss.homo_index, self.gauss.homo_energy))
            logger.info('Current LUMO: %s, energy: %s eV' % (self.gauss.lumo_index, self.gauss.lumo_energy))
            logger.info('HOMO-LUMO Gap: %s eV' % self.gauss.gap)
            if self.get_confirm('Extract HOMO/LUMO cube file?', default=True):
                mos = [self.gauss.homo_index, self.gauss.lumo_index]
            else:
                mos = self.get_input('Enter the MO(s) need to be extracted(separated by comma)').split(',')
            for mo in mos:
                self.gauss.extract_cube(mo)
            #self.create_panel()