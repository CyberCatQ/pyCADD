import os
import logging
import re
logger = logging.getLogger(__name__)


class Gauss:
    '''
    Gauss 计算调用模块
    '''

    def __init__(self, st_path: str, cpu_count: int = None, mem: str = '4GB') -> None:
        self.charge = None                  # 电荷数
        self.spin_multi = None              # 自旋多重度
        self.dft = None                     # 泛函数
        self.basis_set = None               # 基组
        self.st_path = st_path              # 原始结构文件路径
        self.read_origin_st()

        if not cpu_count:
            cpu_count = os.cpu_count()
        self.cpu_count = cpu_count          # 计算将使用的CPU核心数量
        self.mem = mem                      # 计算将使用的内存大小

    def read_origin_st(self):
        '''
        读取原始结构文件名与基础信息
        '''
        # 基于Multiwfn
        st_info = os.popen('''
        Multiwfn %s << EOF
        q
        EOF''' % self.st_path).read().strip()
        self.base_name = os.path.basename(self.st_path).split('.')[0]
        self.molecule_wt = re.search(
            r'(?<=Molecule weight:)[0-9. ]+', st_info).group().strip()
        self.formula = re.search(
            r'(?<=Formula:)[A-Za-z0-9 ]+', st_info).group().strip()
        self.atom_count = re.search(
            r'(?<=Totally)[ 0-9]+', st_info).group().strip()

    def set_charge(self, charge=0):
        '''
        设定电荷量
        '''
        self.charge = str(charge)
        logger.debug('Charge has been set to %s' % self.charge)

    def set_multispin(self, spin_multi=1):
        '''
        设定自旋多重度
        '''
        self.multispin = str(spin_multi)
        logger.debug('Spin multiplicity has been set to %s' % self.multispin)

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

