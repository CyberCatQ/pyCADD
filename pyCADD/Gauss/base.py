import os
import logging
logger = logging.getLogger(__name__)

from pyCADD.Gauss import core

class Gauss:
    '''
    Gaussian 计算调用模块
    '''

    def __init__(self, st_path: str) -> None:
        self.charge = None                  # 电荷数
        self.spin_multi = None              # 自旋多重度
        self.dft = None                     # 泛函数
        self.basis_set = None               # 基组
        self.solvent = None                 # PCM模型溶剂
        self.job = None                     # 任务名

        self.file_type = None               # 文件类型
        self.st_path = st_path              # 原始结构文件路径
        self.read_origin_st()

    @classmethod
    @property
    def gauss(cls):
        return core.get_gaussian()    # 高斯可执行文件路径

    @classmethod
    @property
    def system_info(cls):
        '''
        当前系统计算资源信息
        
        Return
        ----------
        cpu(s), memory
        '''
        cls.cpu_count, cls.mem = core._get_system_info(cls.gauss)
        return cls.cpu_count, cls.mem
    
    def read_origin_st(self):
        '''
        读取原始结构文件名与格式信息
        '''

        self.base_name = os.path.basename(self.st_path).split('.')[0]
        if self.st_path.endswith('.out'):
            self.output_file = self.st_path
            self.file_type = 'Gaussian out(.out)'
        elif self.st_path.endswith('.chk'):
            self.chk_file = self.st_path
            self.file_type = 'check point(.chk)'
        elif self.st_path.endswith('.fchk'):
            self.fchk_file = self.st_path
            self.file_type = 'format check point(.fchk)'
        elif self.st_path.endswith('.gjf'):
            self.input_file = self.st_path
            self.file_type = 'Gaussian input(.gjf)'
        else:
            self.file_type = 'Other'

    @classmethod
    def set_system(cls, cpu_count, mem):
        '''
        计算资源设定

        Parameters
        ----------
        cpu_count : int
            计算使用核心数量
        mem : str
            计算使用内存大小
        '''
        core.system_default(cls.gauss, cpu_count, mem)
        logger.debug('CPU & Memory usage has been changed.')

    def set_charge(self, charge:int=0):
        '''
        设定电荷量
        '''
        self.charge = charge
        logger.debug('Charge has been set to %s' % self.charge)

    def set_multiplicity(self, spin_multi:int=1):
        '''
        设定自旋多重度
        '''
        self.spin_multi = spin_multi
        logger.debug('Spin multiplicity has been set to %s' % self.spin_multi)

    def set_DFT(self, dft: str):
        '''
        设定泛函数
        '''
        self.dft = dft
        logger.debug('DFT has been set to %s' % self.dft)

    def set_basis_set(self, basis_set: str):
        '''
        设定基组
        '''
        self.basis_set = basis_set
        logger.debug('Basis set has been set to %s' % self.basis_set)
    
    def set_solvent(self, solvent: str):
        '''
        设定PCM模型溶剂
        '''
        self.solvent = solvent
        logger.debug('PCM solvent has been set to %s' % self.solvent)

    def _print_current_info(self):
        '''
        显示当前设定状态
        '''

        logger.info('Current structure file: %s' % self.st_path)
        logger.info('Current Setting: Charge = %s  Multiplicity = %s  DFT = %s  Basis_set = %s Solvent = %s' % (self.charge, self.spin_multi, self.dft, self.basis_set, self.solvent))
        
    def create_inputfile(self, job_name:str, loose:bool=True):
        '''
        创建当前设定状态下的输入文件
        Parameter
        ---------
        job_name : str
            任务名
                * opt 结构优化
                * energy 单点能量计算
                * absorb 激发态激发能(吸收)
                * emission 激发态发射能
        loose : bool
            是否提高优化任务中的收敛限 更快收敛
        
        Return
        ----------
        str
            创建的高斯输入文件名称
        '''
        self._print_current_info()
        if job_name == 'opt':
            self.job = 'Optimize'
            self.input_file, self.chk_file = core.generate_opt(self.st_path, self.charge, self.spin_multi, self.dft, self.basis_set, self.solvent, loose)
        elif job_name == 'energy':
            self.job = 'Single Point Energy'
            self.input_file, self.chk_file = core.generate_energy(self.st_path, self.charge, self.spin_multi, self.dft, self.basis_set, self.solvent)
        elif job_name == 'absorb':
            self.job = 'Absorption Energy of Excited States'
            # 计算吸收(激发)能量 KEYWORD即计算TDDFT下的单点能
            self.input_file, self.chk_file = core.generate_energy(self.st_path, self.charge, self.spin_multi, self.dft, self.basis_set, self.solvent, correct=False, td=True)
        elif job_name == 'emission':
            self.job = 'Emission Energy of Excited States'
            # 计算发射能量 Keyword即计算TDDFT下的结构优化
            self.input_file, self.chk_file = core.generate_opt(self.st_path, self.charge, self.spin_multi, self.dft, self.basis_set, self.solvent, loose=False, correct=True, td=True)

        else:
            logger.error('Invaild job name. No file is created.')
            return

        logger.info('Current job name: %s' % self.job)
        logger.info('Input file %s saved.' % self.input_file)

        return self.input_file
    
    def get_mo_info(self):
        '''
        获取HOMO/LUMO分子轨道信息
        '''
        info_dict = core.get_mo(self.fchk_file)
        homo_dict = info_dict['homo']
        lumo_dict = info_dict['lumo']
        gap_value = info_dict['gap']

        self.homo_index = homo_dict['index']
        self.homo_energy = homo_dict['energy']
        self.lumo_index = lumo_dict['index']
        self.lumo_energy = lumo_dict['energy']
        self.gap = gap_value

    def extract_cube(self, mo:int):
        '''
        提取分子轨道cube格点文件

        Parameters
        ----------
        mo : int
            分子轨道(MO)编号
        '''
        logger.debug('Current fchk file: %s' % self.fchk_file)
        logger.info('Extracting MO: %s' % mo)
        cube_file = core.cube_file_generate(self.fchk_file, mo)
        logger.info('MO %s cube file %s created.' % (mo, cube_file))

    def run(self):
        '''
        启动任务
        '''
        logger.info('Curren job: %s' % self.job)
        logger.info('Current input file: %s' % self.input_file)
        logger.info('Current system usage: CPU = %s  Mem = %s' % (self.cpu_count, self.mem) )
        logger.info('Prepare to running calculation')
        self.output_file = self.input_file.split('.')[0] + '.out'

        
        daemon = core.Daemon(cmd = "%s < %s > %s && formchk %s > /dev/null" %(self.gauss, self.input_file, self.output_file, self.chk_file))
        daemon.start()

        logger.info('Job has been submitted. You can safely exit the shell at any time.')

        #logger.info('Start tracing output file ...')
        #core.tail_gauss_job(self.output_file)
        #logger.info('Calculation done. %s is saved.' % self.output_file)

        #self.fchk_file = core.generate_fchk(self.chk_file)
        #fchk_file = self.chk_file.split('.')[0] + '.fchk'
        #logger.info('Formchk file created %s.' % fchk_file)
        

        


